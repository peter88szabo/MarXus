use MarXus::barrierless::phasespace::phase_space_theory::PhaseSpaceTheoryModel;
use MarXus::barrierless::phasespace::types::{
    CaptureFragment, CaptureFragmentRotorModel, PhaseSpaceTheoryInput,
};
use MarXus::masterequation::input_deck::parse_master_equation_input;
use MarXus::masterequation::report::print_species_species_rate_matrix;
use MarXus::molecule::{MolType, MoleculeBuilder};
use MarXus::utils::atomic_masses::mass_vector_from_symbols_amu;

fn main() -> Result<(), String> {
    // -------------------------------------------------------------------------
    // 1) Mirror the provided MESS deck into a MarXus network deck (parsing demo)
    // -------------------------------------------------------------------------
    let deck = include_str!("zzallyloh_o2_multiwell_input.txt");
    let parsed = parse_master_equation_input(deck)?;
    let (_settings, wells) = parsed.build_settings_and_wells()?;

    println!("MESS mirror: ZZAllylOH + O2 (selected subset)");
    println!("Parsed MarXus wells: {}", wells.len());
    for (w, well) in wells.iter().enumerate() {
        println!(
            "  well[{}] = {:<4}  grains [{}..{})  channels={}",
            w,
            well.well_name,
            well.lowest_included_grain_index,
            well.one_past_highest_included_grain_index,
            well.channels.len()
        );
    }

    // -------------------------------------------------------------------------
    // 2) High-pressure capture through the barrierless PST entrance channel
    //    (MESS: Barrier B12, R -> G2, Core PhaseSpaceTheory)
    // -------------------------------------------------------------------------
    //
    // MESS (from the user-provided deck):
    //   SymmetryFactor = 2.0
    //   PotentialPrefactor[au] = 2.4
    //   PotentialPowerExponent = 6
    //
    // We evaluate at the deck's single temperature (300 K) and at 760 Torr for
    // pseudo-first-order collision-frequency reporting.
    // MESS deck lists:
    //   TemperatureList[K] 300.
    //   PressureList[torr] 760.
    let temperature_list_kelvin: Vec<f64> = vec![300.0];
    let pressure_list_torr: Vec<f64> = vec![760.0];

    // deck: ExcessReactantConcentration[molecule/cm^3]
    let excess_reactant_number_density_cm3 = 1.0e6;

    // Fragment geometries (Å) taken verbatim from the provided MESS input.
    // Fragment A: C5H9O3 (17 atoms)
    let fragment_a_symbols: Vec<String> = vec![
        "C", "C", "C", "O", "O", "C", "C", "O", "H", "H", "H", "H", "H", "H", "H", "H", "H",
    ]
    .into_iter()
    .map(|s| s.to_string())
    .collect();
    let fragment_a_coords_angstrom: Vec<[f64; 3]> = vec![
        [0.027010, -0.088846, -0.007340],
        [-0.023222, 0.182884, 1.463884],
        [1.301630, 0.432703, 2.108625],
        [2.094484, -0.746154, 2.300944],
        [1.457234, -1.503513, 3.332338],
        [-1.237541, 0.249803, 2.132462],
        [-1.506516, 0.286525, 3.485481],
        [-0.628454, 0.069974, 4.480515],
        [-2.128881, 0.267327, 1.514056],
        [1.935622, 1.015359, 1.435028],
        [1.221768, 0.969770, 3.054489],
        [-0.968336, -0.222369, -0.426690],
        [0.522985, 0.724696, -0.545477],
        [0.603651, -0.995312, -0.211746],
        [-2.504049, 0.491423, 3.845718],
        [0.131836, -0.432106, 4.144932],
        [2.211926, -1.712875, 3.895438],
    ];

    // Fragment B: O2 (2 atoms)
    let fragment_b_symbols: Vec<String> =
        vec!["O", "O"].into_iter().map(|s| s.to_string()).collect();
    let fragment_b_coords_angstrom: Vec<[f64; 3]> =
        vec![[0.0, 0.0, -0.6060625836], [0.0, 0.0, 0.6060625836]];

    let pst_model = PhaseSpaceTheoryModel::new(PhaseSpaceTheoryInput {
        fragment_a: CaptureFragment {
            mass_amu: None,
            rotor: CaptureFragmentRotorModel::GeometryAngstrom {
                symbols: fragment_a_symbols.clone(),
                coordinates_angstrom: fragment_a_coords_angstrom.clone(),
            },
        },
        fragment_b: CaptureFragment {
            mass_amu: None,
            rotor: CaptureFragmentRotorModel::GeometryAngstrom {
                symbols: fragment_b_symbols.clone(),
                coordinates_angstrom: fragment_b_coords_angstrom.clone(),
            },
        },
        symmetry_operations: 2.0,
        potential_prefactor_au: 2.4,
        potential_power_exponent: 6.0,
    })?;

    // Internal partition-function-like weights (rot * vib * elec) for the fragments.
    let (frag_a_mass_amu, frag_a_brot_cm1) =
        rotational_constants_from_geometry_cm1(&fragment_a_symbols, &fragment_a_coords_angstrom)?;
    let (frag_b_mass_amu, frag_b_brot_cm1) =
        rotational_constants_from_geometry_cm1(&fragment_b_symbols, &fragment_b_coords_angstrom)?;

    // Vibrational frequencies (cm^-1) from the MESS deck.
    // Fragment A has 45; O2 has 1.
    let fragment_a_vib_cm1: Vec<f64> = vec![
        100.00, 100.00, 114.51, 162.81, 201.14, 217.00, 297.35, 314.49, 350.11, 413.89, 463.02,
        579.46, 674.17, 704.44, 750.08, 811.43, 866.03, 879.96, 969.51, 988.41, 1001.53, 1027.73,
        1054.33, 1188.13, 1233.83, 1264.40, 1309.72, 1340.17, 1353.30, 1368.94, 1375.79, 1432.35,
        1440.67, 1454.79, 1470.78, 1498.97, 2953.23, 2985.55, 2998.72, 3037.68, 3063.59, 3089.73,
        3149.13, 3566.30, 3729.52,
    ];
    let o2_vib_cm1: Vec<f64> = vec![1585.09];
    let mut ts_vib_cm1 = fragment_a_vib_cm1.clone();
    ts_vib_cm1.extend_from_slice(&o2_vib_cm1);
    let ts_electronic_degeneracy = 2.0; // MESS: B12 ElectronicLevels: 0 2

    let freq_cutoff_cm1 = 100.0;

    let mut fragment_a = MoleculeBuilder::new("C5H9O3".to_string(), MolType::mol)
        .nlin(false)
        .freq(fragment_a_vib_cm1)
        .brot(frag_a_brot_cm1)
        .mass(frag_a_mass_amu)
        .dh0(0.0)
        .multi(2.0) // MESS: ElectronicLevels 0 2
        .chiral(1.0)
        .symnum(1.0) // MESS: SymmetryFactor 1.0 in the fragment Core
        .build();

    let mut o2 = MoleculeBuilder::new("O2".to_string(), MolType::mol)
        .nlin(true)
        .freq(o2_vib_cm1)
        .brot(frag_b_brot_cm1)
        .mass(frag_b_mass_amu)
        .dh0(0.0)
        .multi(3.0) // MESS: ElectronicLevels 0 3
        .chiral(1.0)
        .symnum(2.0) // MESS: Core SymmetryFactor 2 for O2
        .build();

    // Well G2 (RRHO) from the provided MESS deck excerpt.
    // We build this only so we can also report the reverse dissociation rate k_inf(G2 -> R)
    // in the exact same "Species-Species Rate Tables" format as `steadystate_me_c2h3.rs`.
    let g2_symbols: Vec<String> = vec![
        "O", "C", "C", "C", "C", "O", "C", "O", "O", "O", "H", "H", "H", "H", "H", "H", "H", "H",
        "H",
    ]
    .into_iter()
    .map(|s| s.to_string())
    .collect();
    let g2_coords_angstrom: Vec<[f64; 3]> = vec![
        [1.786799, 2.065626, -0.512394],
        [0.527477, 1.927276, 0.104262],
        [-2.798845, 0.590322, -0.935288],
        [-1.898020, 1.312758, -0.264514],
        [-0.565576, 1.599453, -0.890692],
        [-2.688356, 0.118225, -2.191651],
        [-2.217405, 1.878927, 1.090543],
        [2.149287, 0.778821, -1.004991],
        [-0.645530, 2.770579, -1.775338],
        [-1.216617, 2.476535, -2.897703],
        [0.346902, 2.896765, 0.571240],
        [0.563033, 1.153925, 0.874079],
        [-3.740156, 0.309857, -0.477509],
        [-0.239388, 0.779279, -1.530500],
        [-2.065091, 0.666407, -2.691314],
        [-3.220931, 1.588486, 1.395992],
        [-2.180109, 2.970777, 1.081837],
        [-1.523628, 1.534747, 1.859399],
        [2.961554, 0.612926, -0.512358],
    ];

    let (g2_mass_amu, g2_brot_cm1) =
        rotational_constants_from_geometry_cm1(&g2_symbols, &g2_coords_angstrom)?;

    let g2_vib_cm1: Vec<f64> = vec![
        100.00, 100.00, 100.00, 157.45, 174.24, 197.98, 205.92, 226.46, 231.87, 264.20, 307.19,
        370.96, 420.48, 513.45, 538.75, 586.39, 594.91, 618.61, 788.23, 860.96, 927.60, 936.88,
        984.21, 1024.11, 1052.71, 1079.05, 1136.42, 1179.89, 1206.58, 1267.94, 1272.09, 1278.34,
        1311.33, 1339.36, 1364.46, 1375.36, 1388.24, 1397.94, 1419.49, 1447.25, 1462.03, 1697.32,
        2971.61, 2991.04, 3027.57, 3034.71, 3054.36, 3061.86, 3119.57, 3627.95, 3730.28,
    ];

    let mut g2 = MoleculeBuilder::new("G2".to_string(), MolType::mol)
        .nlin(false)
        .freq(g2_vib_cm1)
        .brot(g2_brot_cm1)
        .mass(g2_mass_amu)
        .dh0(0.0)
        .multi(2.0) // MESS: ElectronicLevels 0 2
        .chiral(1.0)
        .symnum(0.5) // MESS: Core SymmetryFactor 0.5 (mirrored literally)
        .build();

    println!();
    println!("Barrierless capture (MESS Barrier B12: R -> G2, PhaseSpaceTheory core)");
    println!("PST params: symmetry_ops = 2.0, V0 = 2.4 au, n = 6");
    println!();
    println!("Diagnostics (capture only):");
    println!("T/K    P/torr     Q_PST(T)        q_int(A)        q_int(O2)        R -> G2: k_cap (cm^3 molecule^-1 s^-1)     k_cap[M] (s^-1)        k_cap*[R_excess] (s^-1)");

    let mut capture_rows: Vec<CaptureRow> = Vec::new();

    for &temperature_kelvin in &temperature_list_kelvin {
        for &pressure_torr in &pressure_list_torr {
            let pressure_pa = pressure_torr_to_pa(pressure_torr);

            fragment_a.eval_all_therm_func(temperature_kelvin, pressure_pa, freq_cutoff_cm1);
            o2.eval_all_therm_func(temperature_kelvin, pressure_pa, freq_cutoff_cm1);
            g2.eval_all_therm_func(temperature_kelvin, pressure_pa, freq_cutoff_cm1);

            let q_int_a =
                fragment_a.thermo.pfelec * fragment_a.thermo.pfrot * fragment_a.thermo.pfvib;
            let q_int_b = o2.thermo.pfelec * o2.thermo.pfrot * o2.thermo.pfvib;
            let _q_int_g2 = g2.thermo.pfelec * g2.thermo.pfrot * g2.thermo.pfvib;

            let temperature_au = k_b_t_hartree(temperature_kelvin);
            let q_pst = pst_model.canonical_weight_from_kbt_hartree(temperature_au)?;
            let q_vib_ts = vibrational_partition_function_quantum(&ts_vib_cm1, temperature_kelvin);

            // IMPORTANT:
            // In the MESS input, barrier B12 is an `RRHO` model with:
            //   Core PhaseSpaceTheory  (this contributes `Q_PST(T)`)
            //   Frequencies[1/cm]      (this contributes `Q_vib(T)`)
            //   ElectronicLevels       (this contributes `Q_elec(T)`)
            //
            // Therefore the TS canonical weight entering the thermal high-pressure capture rate is:
            //   W_TS(T) = Q_PST(T) * Q_vib(T) * Q_elec(T)
            //
            // Previously we used only `Q_PST(T)`, which underestimates capture by ~2 orders of magnitude
            // for this system at 300 K because the many low-frequency modes give Q_vib(T) >> 1.
            let w_ts_total = q_pst * q_vib_ts * ts_electronic_degeneracy;

            // Thermal high-pressure capture coefficient:
            // k(T) = (T/2π) * W_TS(T) / W_AB(T)   (MESS convention, atomic units)
            let w_ab = bimolecular_weight_per_volume_atomic_units(
                frag_a_mass_amu,
                frag_b_mass_amu,
                q_int_a,
                q_int_b,
                temperature_au,
            );
            let k_capture_au =
                (temperature_au / (2.0 * std::f64::consts::PI)) * (w_ts_total / w_ab);
            let k_capture_cm3_molecule_s = k_capture_au * AU_VOLUME_PER_TIME_CM3_S;

            let m_number_density_cm3 = number_density_cm3(pressure_pa, temperature_kelvin);
            let pseudo_first_order = k_capture_cm3_molecule_s * m_number_density_cm3;
            let pseudo_first_order_excess =
                k_capture_cm3_molecule_s * excess_reactant_number_density_cm3;

            capture_rows.push(CaptureRow {
                temperature_kelvin,
                pressure_torr,
                k_capture_cm3_molecule_s,
                pseudo_first_order_at_pressure: pseudo_first_order,
                pseudo_first_order_at_excess: pseudo_first_order_excess,
            });

            println!(
                "{:<5.0} {:>7.1} {:>14.6e} {:>14.6e} {:>14.6e} {:>22.6e} {:>18.6e} {:>22.6e}",
                temperature_kelvin,
                pressure_torr,
                q_pst,
                q_int_a,
                q_int_b,
                k_capture_cm3_molecule_s,
                pseudo_first_order,
                pseudo_first_order_excess
            );

            // -----------------------------------------------------------------
            // Full "Species-Species Rate Tables" matrix (MESS-style)
            // -----------------------------------------------------------------
            //
            // This mirror example currently computes only the entrance capture `R -> G2`
            // from the PST model. All other entries are left as "***".
            // Keep the same ordering as in your requested MESS-like table.
            let species: Vec<String> = vec![
                "G2".to_string(),
                "G3".to_string(),
                "G4".to_string(),
                "G6".to_string(),
                "G13".to_string(),
                "R".to_string(),
                "P1".to_string(),
                "P5".to_string(),
                "P7".to_string(),
            ];

            let n = species.len();
            let mut hp: Vec<Vec<Option<f64>>> = vec![vec![None; n]; n];

            let idx = |name: &str| -> usize {
                species
                    .iter()
                    .position(|s| s == name)
                    .expect("species list missing name")
            };

            // High-pressure capture: R -> G2 (bimolecular, cm^3 molecule^-1 s^-1)
            hp[idx("R")][idx("G2")] = Some(k_capture_cm3_molecule_s);

            // NOTE: we intentionally do not fill `G2 -> R` here; that needs a proper
            // exit-channel microcanonical model + detailed balance.

            print_species_species_rate_matrix(temperature_kelvin, None, &species, &hp)?;
            print_species_species_rate_matrix(temperature_kelvin, Some("760 torr"), &species, &hp)?;
        }
    }

    // Dedicated capture section (explicit label, clean at the end)
    println!();
    println!("================================================================================");
    println!("Barrierless capture (PST) — dedicated summary");
    println!("Reaction: R -> G2");
    println!("================================================================================");
    for row in &capture_rows {
        println!(
            "T={:>6.1} K  P={:>7.1} torr   R -> G2: k_cap={:>12.6e} cm^3 molecule^-1 s^-1   k_cap[M]={:>12.6e} s^-1   k_cap*[R_excess]={:>12.6e} s^-1",
            row.temperature_kelvin,
            row.pressure_torr,
            row.k_capture_cm3_molecule_s,
            row.pseudo_first_order_at_pressure,
            row.pseudo_first_order_at_excess
        );
    }

    Ok(())
}

#[derive(Clone, Debug)]
struct CaptureRow {
    temperature_kelvin: f64,
    pressure_torr: f64,
    k_capture_cm3_molecule_s: f64,
    pseudo_first_order_at_pressure: f64,
    pseudo_first_order_at_excess: f64,
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
    // Quantum harmonic-oscillator partition function per mode:
    //   q = 1 / (1 - exp(-h c ω / (k_B T))) = 1 / (1 - exp(-ω / (kB_cm1_per_K * T)))
    //
    // This matches MESS' RRHO weight() vibrational factor.
    const KB_CM1_PER_K: f64 = 0.695_034_76;
    let t = temperature_kelvin.max(1e-12);
    let beta = 1.0 / (KB_CM1_PER_K * t);

    let mut q = 1.0;
    for &w in frequencies_cm1 {
        if !(w > 0.0) || !w.is_finite() {
            continue;
        }
        let x = w * beta;
        // Guard for extreme x: exp(-x) underflows to 0 -> q -> 1, which is fine.
        let e = (-x).exp();
        let denom = 1.0 - e;
        if denom > 0.0 {
            q *= 1.0 / denom;
        }
    }
    q
}
