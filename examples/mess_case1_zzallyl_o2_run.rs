use std::collections::HashMap;

use MarXus::masterequation::api::{
    run_multiwell_chemical_activation_with_diagnostics, MultiWellChemicalActivationInput,
};
use MarXus::masterequation::mess_input::{parse_mess_input, MessBuildOptions, MessDeck};
use MarXus::masterequation::microcanonical_builder::{
    ArrayMicrocanonicalProvider, MicrocanonicalNetworkData, TransitionStateModel,
};
use MarXus::masterequation::network_builder::build_microcanonical_network_data_from_named_models;
use MarXus::masterequation::reaction_network::{
    ChemicalActivationDefinition, MicrocanonicalProvider, WellDefinition,
};
use MarXus::masterequation::report::{
    build_species_species_rate_matrix_from_solution, print_species_species_rate_matrix,
};
use MarXus::molecule::{MolType, MoleculeBuilder};
use MarXus::utils::atomic_masses::mass_vector_from_symbols_amu;

fn main() -> Result<(), String> {
    // Parse the original MESS deck (placed in `examples/`).
    let mess = include_str!("Case1_ZZ_allyl+O2_2025_12_02.inp");
    let deck = parse_mess_input(mess)?;

    // Build topology + micro models. Keep max energy modest so the current dense ME engine fits in RAM.
    let built = deck.build_multiwell_network(MessBuildOptions {
        energy_grain_width_cm1: Some(20.0),
        max_energy_kcal_mol: 50.0,
        collision_band_half_width: 20,
        nonreactive_grain_count: 10,
    })?;

    println!("MESS import + multiwell ME: ZZAllylOH + O2 (Case1)");
    println!(
        "Temperature = {:.1} K    Pressure = {:.1} torr",
        deck.global.temperature_kelvin, deck.global.pressure_torr
    );
    println!("Wells: {}", built.wells.len());
    for (i, w) in built.wells.iter().enumerate() {
        println!(
            "  well[{}] = {:<4}  grains [{}..{})  channels={}",
            i,
            w.well_name,
            w.lowest_included_grain_index,
            w.one_past_highest_included_grain_index,
            w.channels.len()
        );
    }

    // Microcanonical arrays ρ(E) and k(E).
    let micro_data: MicrocanonicalNetworkData =
        build_microcanonical_network_data_from_named_models(
            &built.wells,
            &built.well_models_by_name,
            &built.channel_models_by_name,
        )?;
    let micro = ArrayMicrocanonicalProvider::new(&micro_data);

    // Chemical activation: MESS reactant `R` captures into `G2` via barrier `B12`.
    let activated_well_index = built
        .wells
        .iter()
        .position(|w| w.well_name == "G2")
        .ok_or("Expected well 'G2'")?;
    let recombination_channel_index = built.wells[activated_well_index]
        .channels
        .iter()
        .position(|ch| ch.name == "G2_to_R")
        .ok_or("Expected channel 'G2_to_R' (from MESS barrier B12)")?;

    let activation = ChemicalActivationDefinition {
        activated_well_index,
        recombination_channel_index,
    };

    let (solution, _diag) = run_multiwell_chemical_activation_with_diagnostics(
        MultiWellChemicalActivationInput {
            wells: built.wells.clone(),
            settings: built.settings.clone(),
            activation,
        },
        &micro,
    )?;

    // Compute barrierless thermal capture coefficient (bimolecular): R -> G2.
    let (k_capture_cm3_molecule_s, capture_diag) =
        compute_capture_r_to_g2_from_mess_deck(&deck, &built.channel_models_by_name)?;

    // -------------------------
    // High-pressure table
    // -------------------------
    let (mut hp_species, mut hp_matrix) = build_high_pressure_rate_matrix(
        deck.global.temperature_kelvin,
        0.695_034_76,
        &built.wells,
        &micro,
    )?;
    insert_entry_named(
        &mut hp_species,
        &mut hp_matrix,
        "R",
        "G2",
        k_capture_cm3_molecule_s,
    );

    print_species_species_rate_matrix(
        deck.global.temperature_kelvin,
        None,
        &hp_species,
        &hp_matrix,
    )?;

    // -------------------------
    // Pressure-dependent table (ME macroscopic rates)
    // -------------------------
    let (mut species_names, mut me_matrix) =
        build_species_species_rate_matrix_from_solution(&built.wells, &solution)?;
    insert_entry_named(
        &mut species_names,
        &mut me_matrix,
        "R",
        "G2",
        k_capture_cm3_molecule_s,
    );

    print_species_species_rate_matrix(
        deck.global.temperature_kelvin,
        Some(&format!("{:.0} torr", deck.global.pressure_torr)),
        &species_names,
        &me_matrix,
    )?;

    // Dedicated capture section at end.
    println!();
    println!("================================================================================");
    println!("Barrierless capture (PST) — dedicated summary");
    println!("Reaction: R -> G2");
    println!("================================================================================");
    println!(
        "T={:>6.1} K  P={:>7.1} torr   R -> G2: k_cap={:>12.6e} cm^3 molecule^-1 s^-1   k_cap[M]={:>12.6e} s^-1   k_cap*[R_excess]={:>12.6e} s^-1",
        capture_diag.temperature_kelvin,
        capture_diag.pressure_torr,
        capture_diag.k_capture_cm3_molecule_s,
        capture_diag.pseudo_first_order_at_pressure,
        capture_diag.pseudo_first_order_at_excess
    );

    Ok(())
}

#[derive(Clone, Debug)]
struct CaptureDiagnostics {
    temperature_kelvin: f64,
    pressure_torr: f64,
    k_capture_cm3_molecule_s: f64,
    pseudo_first_order_at_pressure: f64,
    pseudo_first_order_at_excess: f64,
}

fn compute_capture_r_to_g2_from_mess_deck(
    deck: &MessDeck,
    channel_models_by_name: &HashMap<
        String,
        MarXus::masterequation::microcanonical_builder::ChannelMicroModel,
    >,
) -> Result<(f64, CaptureDiagnostics), String> {
    let temperature_kelvin = deck.global.temperature_kelvin;
    let pressure_torr = deck.global.pressure_torr;
    let pressure_pa = pressure_torr_to_pa(pressure_torr);
    let excess_reactant_number_density_cm3 =
        deck.global.excess_reactant_concentration_cm3.unwrap_or(0.0);

    let reactant_name = deck
        .global
        .reactant_name
        .clone()
        .unwrap_or_else(|| "R".to_string());
    let bimol = deck
        .bimolecular
        .get(&reactant_name)
        .ok_or_else(|| format!("Missing bimolecular reactant '{}'", reactant_name))?;

    let ch_model = channel_models_by_name
        .get("G2_to_R")
        .ok_or("Missing channel model 'G2_to_R'")?;

    let (pst_core, ts_vib_cm1, ts_electronic_degeneracy) = match &ch_model.transition_state {
        TransitionStateModel::PhaseSpaceTheoryRRHO {
            pst_core,
            vibrational_frequencies_cm1,
            electronic_degeneracy,
        } => (
            pst_core,
            vibrational_frequencies_cm1,
            *electronic_degeneracy,
        ),
        _ => return Err("Channel 'G2_to_R' is not a PhaseSpaceTheoryRRHO TS model.".into()),
    };

    // Internal weights for the two reactant fragments.
    let (mass_a_amu, brot_a_cm1) = rotational_constants_from_geometry_cm1(
        &bimol.fragment_a.geometry_symbols,
        &bimol.fragment_a.geometry_angstrom,
    )?;
    let (mass_b_amu, brot_b_cm1) = rotational_constants_from_geometry_cm1(
        &bimol.fragment_b.geometry_symbols,
        &bimol.fragment_b.geometry_angstrom,
    )?;

    let is_linear_a = bimol.fragment_a.geometry_symbols.len() == 2;
    let is_linear_b = bimol.fragment_b.geometry_symbols.len() == 2;

    let freq_cutoff_cm1 = 100.0;

    let mut frag_a = MoleculeBuilder::new(bimol.fragment_a.name.clone(), MolType::mol)
        .nlin(is_linear_a)
        .freq(bimol.fragment_a.vibrational_frequencies_cm1.clone())
        .brot(brot_a_cm1)
        .mass(mass_a_amu)
        .dh0(0.0)
        .multi(bimol.fragment_a.electronic_degeneracy_ground)
        .chiral(1.0)
        .symnum(bimol.fragment_a.symmetry_factor)
        .build();

    let mut frag_b = MoleculeBuilder::new(bimol.fragment_b.name.clone(), MolType::mol)
        .nlin(is_linear_b)
        .freq(bimol.fragment_b.vibrational_frequencies_cm1.clone())
        .brot(brot_b_cm1)
        .mass(mass_b_amu)
        .dh0(0.0)
        .multi(bimol.fragment_b.electronic_degeneracy_ground)
        .chiral(1.0)
        .symnum(bimol.fragment_b.symmetry_factor)
        .build();

    frag_a.eval_all_therm_func(temperature_kelvin, pressure_pa, freq_cutoff_cm1);
    frag_b.eval_all_therm_func(temperature_kelvin, pressure_pa, freq_cutoff_cm1);

    let q_int_a = frag_a.thermo.pfelec * frag_a.thermo.pfrot * frag_a.thermo.pfvib;
    let q_int_b = frag_b.thermo.pfelec * frag_b.thermo.pfrot * frag_b.thermo.pfvib;

    let temperature_au = k_b_t_hartree(temperature_kelvin);
    let q_pst = pst_core.canonical_weight_from_kbt_hartree(temperature_au)?;
    let q_vib_ts = vibrational_partition_function_quantum(ts_vib_cm1, temperature_kelvin);
    let w_ts_total = q_pst * q_vib_ts * ts_electronic_degeneracy;

    let w_ab = bimolecular_weight_per_volume_atomic_units(
        mass_a_amu,
        mass_b_amu,
        q_int_a,
        q_int_b,
        temperature_au,
    );
    let k_capture_au = (temperature_au / (2.0 * std::f64::consts::PI)) * (w_ts_total / w_ab);
    let k_capture_cm3_molecule_s = k_capture_au * AU_VOLUME_PER_TIME_CM3_S;

    let m_number_density_cm3 = number_density_cm3(pressure_pa, temperature_kelvin);
    let pseudo_first_order = k_capture_cm3_molecule_s * m_number_density_cm3;
    let pseudo_first_order_excess = k_capture_cm3_molecule_s * excess_reactant_number_density_cm3;

    Ok((
        k_capture_cm3_molecule_s,
        CaptureDiagnostics {
            temperature_kelvin,
            pressure_torr,
            k_capture_cm3_molecule_s,
            pseudo_first_order_at_pressure: pseudo_first_order,
            pseudo_first_order_at_excess: pseudo_first_order_excess,
        },
    ))
}

fn build_high_pressure_rate_matrix(
    temperature_kelvin: f64,
    boltzmann_cm1_per_k: f64,
    wells: &[WellDefinition],
    micro: &dyn MicrocanonicalProvider,
) -> Result<(Vec<String>, Vec<Vec<Option<f64>>>), String> {
    // Collect species = all wells + all sink endpoints found in channel names.
    let mut names: Vec<String> = Vec::new();
    for w in wells {
        if !names.contains(&w.well_name) {
            names.push(w.well_name.clone());
        }
    }
    for w in wells {
        for ch in &w.channels {
            if let Some((_from, to)) = ch.name.split_once("_to_") {
                let to = to.to_string();
                if !names.contains(&to) {
                    names.push(to);
                }
            }
        }
    }

    let mut index_of: HashMap<String, usize> = HashMap::new();
    for (i, s) in names.iter().enumerate() {
        index_of.insert(s.clone(), i);
    }

    let n = names.len();
    let mut matrix: Vec<Vec<Option<f64>>> = vec![vec![None; n]; n];

    for (w_idx, w) in wells.iter().enumerate() {
        let from_i = *index_of
            .get(&w.well_name)
            .ok_or_else(|| "Missing well in species list".to_string())?;

        let denom = canonical_partition_weight_for_well(
            temperature_kelvin,
            boltzmann_cm1_per_k,
            w_idx,
            w,
            micro,
        )?;
        if denom <= 0.0 {
            continue;
        }

        for (ch_idx, ch) in w.channels.iter().enumerate() {
            let to_name = ch
                .name
                .split_once("_to_")
                .map(|(_, to)| to)
                .unwrap_or(ch.name.as_str())
                .to_string();
            let to_j = match index_of.get(&to_name) {
                Some(v) => *v,
                None => continue,
            };

            let num = canonical_weighted_channel_rate_sum(
                temperature_kelvin,
                boltzmann_cm1_per_k,
                w_idx,
                w,
                micro,
                ch_idx,
            )?;
            let k_inf = num / denom;
            if k_inf.is_finite() && k_inf > 0.0 {
                matrix[from_i][to_j] = Some(k_inf);
            }
        }
    }

    Ok((names, matrix))
}

fn canonical_partition_weight_for_well(
    temperature_kelvin: f64,
    boltzmann_cm1_per_k: f64,
    well_index: usize,
    well: &WellDefinition,
    micro: &dyn MicrocanonicalProvider,
) -> Result<f64, String> {
    let mut sum = 0.0;
    let t = temperature_kelvin;
    let kb = boltzmann_cm1_per_k;
    let d_e = well.energy_grain_width_cm1;

    for local_grain in well.lowest_included_grain_index..well.one_past_highest_included_grain_index
    {
        let rho = micro.density_of_states(well_index, local_grain);
        if rho <= 0.0 {
            continue;
        }
        let abs_grain = (local_grain as isize) + well.alignment_offset_in_grains;
        let e_abs = (abs_grain as f64) * d_e;
        let w = (-e_abs / (kb * t)).exp();
        sum += rho * w;
    }

    Ok(sum)
}

fn canonical_weighted_channel_rate_sum(
    temperature_kelvin: f64,
    boltzmann_cm1_per_k: f64,
    well_index: usize,
    well: &WellDefinition,
    micro: &dyn MicrocanonicalProvider,
    channel_index: usize,
) -> Result<f64, String> {
    let mut sum = 0.0;
    let t = temperature_kelvin;
    let kb = boltzmann_cm1_per_k;
    let d_e = well.energy_grain_width_cm1;

    for local_grain in well.lowest_included_grain_index..well.one_past_highest_included_grain_index
    {
        if local_grain < well.nonreactive_grain_count {
            continue;
        }
        let rho = micro.density_of_states(well_index, local_grain);
        if rho <= 0.0 {
            continue;
        }
        let k_e = micro
            .microcanonical_rate(well_index, channel_index, local_grain)
            .max(0.0);
        if k_e <= 0.0 {
            continue;
        }
        let abs_grain = (local_grain as isize) + well.alignment_offset_in_grains;
        let e_abs = (abs_grain as f64) * d_e;
        let w = (-e_abs / (kb * t)).exp();
        sum += rho * w * k_e;
    }

    Ok(sum)
}

fn insert_entry_named(
    species: &mut Vec<String>,
    matrix: &mut Vec<Vec<Option<f64>>>,
    from: &str,
    to: &str,
    value: f64,
) {
    let from_i = species.iter().position(|s| s == from);
    let to_j = species.iter().position(|s| s == to);
    if let (Some(i), Some(j)) = (from_i, to_j) {
        matrix[i][j] = Some(value);
    }
}

fn k_b_t_hartree(temperature_kelvin: f64) -> f64 {
    const BOLTZMANN_CM1_PER_K: f64 = 0.695_034_76;
    const CM1_TO_HARTREE: f64 = 4.55635e-6;
    temperature_kelvin * BOLTZMANN_CM1_PER_K * CM1_TO_HARTREE
}

const AMU_TO_ELECTRON_MASS: f64 = 1822.888_486_209;
const BOHR_TO_CM: f64 = 0.529_177_210_903e-8;
const ATOMIC_TIME_TO_S: f64 = 2.418_884_326_585_7e-17;
const AU_VOLUME_PER_TIME_CM3_S: f64 = (BOHR_TO_CM * BOHR_TO_CM * BOHR_TO_CM) / ATOMIC_TIME_TO_S;

fn bimolecular_weight_per_volume_atomic_units(
    mass_a_amu: f64,
    mass_b_amu: f64,
    internal_weight_a: f64,
    internal_weight_b: f64,
    temperature_au: f64,
) -> f64 {
    let mass_a_au = mass_a_amu * AMU_TO_ELECTRON_MASS;
    let mass_b_au = mass_b_amu * AMU_TO_ELECTRON_MASS;
    let reduced_mass_au = (mass_a_au * mass_b_au) / (mass_a_au + mass_b_au);

    let translational_prefactor_per_volume =
        (reduced_mass_au * temperature_au / (2.0 * std::f64::consts::PI)).powf(1.5);
    translational_prefactor_per_volume * internal_weight_a * internal_weight_b
}

fn number_density_cm3(pressure_pa: f64, temperature_kelvin: f64) -> f64 {
    const BOLTZMANN_SI: f64 = 1.380_649e-23; // J/K
    (pressure_pa / (BOLTZMANN_SI * temperature_kelvin)) / 1.0e6
}

fn pressure_torr_to_pa(pressure_torr: f64) -> f64 {
    // exact by definition: 1 atm = 760 Torr = 101325 Pa
    pressure_torr * (101_325.0 / 760.0)
}

fn rotational_constants_from_geometry_cm1(
    symbols: &[String],
    coordinates_angstrom: &[[f64; 3]],
) -> Result<(f64, Vec<f64>), String> {
    if symbols.is_empty()
        || coordinates_angstrom.is_empty()
        || symbols.len() != coordinates_angstrom.len()
    {
        return Err("Geometry must have matching non-empty symbols and coordinates.".into());
    }

    let masses_amu = mass_vector_from_symbols_amu(symbols)?;
    let mass_amu_total = masses_amu.iter().sum::<f64>();

    if symbols.len() == 1 {
        return Ok((mass_amu_total, Vec::new()));
    }

    let coords: Vec<[f64; 3]> = coordinates_angstrom.to_vec();
    let brot = MarXus::inertia::inertia::get_brot(&coords, &masses_amu);

    let mut finite: Vec<f64> = brot
        .into_iter()
        .filter(|b| b.is_finite() && *b > 0.0 && *b < 1.0e6)
        .collect();

    if finite.len() == 2 {
        let b = 0.5 * (finite[0] + finite[1]);
        return Ok((mass_amu_total, vec![b]));
    }

    if finite.len() >= 3 {
        finite.truncate(3);
        return Ok((mass_amu_total, finite));
    }

    Err("Failed to compute usable rotational constants from geometry.".into())
}

fn vibrational_partition_function_quantum(frequencies_cm1: &[f64], temperature_kelvin: f64) -> f64 {
    const KB_CM1_PER_K: f64 = 0.695_034_76;
    let t = temperature_kelvin.max(1e-12);
    let beta = 1.0 / (KB_CM1_PER_K * t);

    let mut q = 1.0;
    for &w in frequencies_cm1 {
        if !(w > 0.0) || !w.is_finite() {
            continue;
        }
        let x = w * beta;
        let e = (-x).exp();
        let denom = 1.0 - e;
        if denom > 0.0 {
            q *= 1.0 / denom;
        }
    }
    q
}
