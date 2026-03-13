use std::collections::{BTreeSet, HashMap};

use MarXus::masterequation::mess_input::{parse_mess_input, MessBuildOptions};
use MarXus::masterequation::microcanonical_builder::{
    ArrayMicrocanonicalProvider, MicrocanonicalNetworkData,
};
use MarXus::masterequation::network_builder::build_microcanonical_network_data_from_named_models;
use MarXus::masterequation::report::print_species_species_rate_matrix;

fn main() -> Result<(), String> {
    // End-to-end demo:
    //   MESS deck -> MarXus network topology -> microcanonical arrays (rho(E), k(E))
    //   -> canonical high-pressure averaged rate table (all species in the topology)
    //
    // This example uses a *small embedded* excerpt deck for compilation + smoke testing.
    // Replace the include_str! with your real MESS deck content when ready.
    let mess = include_str!("mess_zzallyl_o2_case1_excerpt.inp");
    let deck = parse_mess_input(mess)?;

    let built = deck.build_multiwell_network(MessBuildOptions {
        energy_grain_width_cm1: Some(20.0),
        max_energy_kcal_mol: 50.0,
        collision_band_half_width: 20,
        nonreactive_grain_count: 10,
    })?;

    // Build microcanonical arrays for every well and every channel.
    let micro_data: MicrocanonicalNetworkData =
        build_microcanonical_network_data_from_named_models(
            &built.wells,
            &built.well_models_by_name,
            &built.channel_models_by_name,
        )?;
    let micro = ArrayMicrocanonicalProvider::new(&micro_data);

    // Collect "species" list for reporting:
    // - all wells
    // - all sink endpoints that appear in channel names as "..._to_<sink>"
    let mut species: BTreeSet<String> = BTreeSet::new();
    for w in &built.wells {
        species.insert(w.well_name.clone());
        for ch in &w.channels {
            if let Some((_from, to)) = ch.name.split_once("_to_") {
                species.insert(to.to_string());
            }
        }
    }
    let species: Vec<String> = species.into_iter().collect();
    let index_of: HashMap<String, usize> = species
        .iter()
        .enumerate()
        .map(|(i, s)| (s.clone(), i))
        .collect();

    // High-pressure (canonical) rate coefficients:
    // For each well w and channel c, compute:
    //   k_inf(w->c) = < k(E) >_eq,w
    // where the equilibrium weight is proportional to rho(E)*exp(-E_abs/kBT).
    //
    // NOTE: This gives 1/s for unimolecular channels (well->well or well->sink).
    // Bimolecular capture rates (e.g. R->G2) are not produced by this averaging;
    // they are handled separately via PST or TST capture models.
    let hp_matrix = build_high_pressure_rate_matrix(
        deck.global.temperature_kelvin,
        deck.global.boltzmann_constant_wavenumber_per_kelvin(),
        &built.wells,
        &micro,
        &species,
        &index_of,
    )?;

    print_species_species_rate_matrix(deck.global.temperature_kelvin, None, &species, &hp_matrix)?;

    Ok(())
}

fn build_high_pressure_rate_matrix(
    temperature_kelvin: f64,
    boltzmann_cm1_per_k: f64,
    wells: &[MarXus::masterequation::reaction_network::WellDefinition],
    micro: &dyn MarXus::masterequation::reaction_network::MicrocanonicalProvider,
    species: &[String],
    index_of: &HashMap<String, usize>,
) -> Result<Vec<Vec<Option<f64>>>, String> {
    let n = species.len();
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

    Ok(matrix)
}

fn canonical_partition_weight_for_well(
    temperature_kelvin: f64,
    boltzmann_cm1_per_k: f64,
    well_index: usize,
    well: &MarXus::masterequation::reaction_network::WellDefinition,
    micro: &dyn MarXus::masterequation::reaction_network::MicrocanonicalProvider,
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
    well: &MarXus::masterequation::reaction_network::WellDefinition,
    micro: &dyn MarXus::masterequation::reaction_network::MicrocanonicalProvider,
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

trait MessGlobalExt {
    fn boltzmann_constant_wavenumber_per_kelvin(&self) -> f64;
}

impl MessGlobalExt for MarXus::masterequation::mess_input::MessGlobal {
    fn boltzmann_constant_wavenumber_per_kelvin(&self) -> f64 {
        0.695_034_76
    }
}
