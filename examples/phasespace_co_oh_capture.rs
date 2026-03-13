use MarXus::barrierless::phasespace::phase_space_theory::PhaseSpaceTheoryModel;
use MarXus::barrierless::phasespace::types::{
    CaptureFragment, CaptureFragmentRotorModel, PhaseSpaceTheoryInput,
};
use MarXus::molecule::{MolType, MoleculeBuilder};

fn main() -> Result<(), String> {
    // Reuse the same basic CO + OH constants as `examples/co_oh.rs`.
    let temperatures_kelvin = vec![300.0, 500.0, 1000.0, 1500.0];

    // Diatomic vibrational frequencies and equilibrium bond lengths (from the existing example).
    let oh_frequency_cm1 = 3737.76;
    let oh_bond_length_angstrom = 0.969;
    let co_frequency_cm1 = 2169.814;
    let co_bond_length_angstrom = 1.128;

    // Derived rigid-rotor rotational constants (cm^-1).
    let oh_brot_cm1 = diatomic_brot_cm1(15.999, 1.008, oh_bond_length_angstrom);
    let co_brot_cm1 = diatomic_brot_cm1(12.0, 15.999, co_bond_length_angstrom);

    // --------------------------
    // Barrierless capture inputs
    // --------------------------
    //
    // The PST model requires long-range potential parameters:
    //   V(R) = V0 / R^n
    // with V0 in atomic units and n > 1.
    //
    // NOTE:
    // The existing CO+OH examples in this repo provide fragment properties (masses, B constants,
    // vibrational frequencies, etc.) but do not define a long-range potential.
    //
    // For now we keep V0 and n as explicit example constants (so the example runs out of the box).
    // If you have ab initio / fitted long-range parameters for CO+OH, replace them here.
    let potential_power_exponent: f64 = 6.0;
    let potential_prefactor_au: f64 = 10.0;
    let symmetry_operations: f64 = 1.0;

    // Used only to convert bimolecular k(T) to a pseudo-first-order collision frequency k(T)[M].
    let pressure_pa: f64 = 101_325.0;

    // Used only for an order-of-magnitude collision-limit check.
    let collision_diameter_angstrom: f64 = 3.5;

    // Build the PST (phase-space) TS model.
    let pst_model = PhaseSpaceTheoryModel::new(PhaseSpaceTheoryInput {
        fragment_a: CaptureFragment {
            mass_amu: Some(28.0),
            rotor: CaptureFragmentRotorModel::LinearRigidRotor {
                rotational_constant_cm1: co_brot_cm1,
            },
        },
        fragment_b: CaptureFragment {
            mass_amu: Some(17.0),
            rotor: CaptureFragmentRotorModel::LinearRigidRotor {
                rotational_constant_cm1: oh_brot_cm1,
            },
        },
        symmetry_operations,
        potential_prefactor_au,
        potential_power_exponent,
    })?;

    // Build reactant molecules to get internal partition functions (rot+vib+elec).
    //
    // We intentionally use the internal PF:
    //   q_int = q_elec * q_rot * q_vib
    // and then use the standard bimolecular Eyring factor (kBT/P) explicitly, so that
    // translation is handled consistently.
    let freq_cutoff_cm1 = 100.0;

    let mut co = MoleculeBuilder::new("CO".to_string(), MolType::mol)
        .nlin(true)
        .freq(vec![co_frequency_cm1])
        .brot(vec![co_brot_cm1])
        .mass(28.0)
        .dh0(0.0)
        .multi(1.0)
        .chiral(1.0)
        .symnum(1.0)
        .build();

    let mut oh = MoleculeBuilder::new("OH".to_string(), MolType::mol)
        .nlin(true)
        .freq(vec![oh_frequency_cm1])
        .brot(vec![oh_brot_cm1])
        .mass(17.0)
        .dh0(0.0)
        .multi(1.0)
        .chiral(1.0)
        .symnum(1.0)
        .build();

    println!("CO + OH: PST-based high-pressure capture (illustrative)");
    println!(
        "PST params: n = {:.3}, V0 = {:.6} au, symmetry_ops = {:.1}",
        potential_power_exponent, potential_prefactor_au, symmetry_operations
    );
    println!(
        "Pressure: {:.3} atm (used only for pseudo-first-order collision-frequency reporting)",
        pressure_pa / 101_325.0
    );
    println!();
    println!(
        "T/K     q_int(CO)    q_int(OH)    Q_PST(T)     k_cap (cm^3 molecule^-1 s^-1)   k_cap*[M] (s^-1)     k_coll_limit (cm^3 molecule^-1 s^-1)"
    );

    for &t in &temperatures_kelvin {
        co.eval_all_therm_func(t, pressure_pa, freq_cutoff_cm1);
        oh.eval_all_therm_func(t, pressure_pa, freq_cutoff_cm1);

        let q_int_co = co.thermo.pfelec * co.thermo.pfrot * co.thermo.pfvib;
        let q_int_oh = oh.thermo.pfelec * oh.thermo.pfrot * oh.thermo.pfvib;

        let temperature_au = k_b_t_hartree(t);
        let q_pst = pst_model.canonical_weight_from_kbt_hartree(temperature_au)?;

        // IMPORTANT (units/physics):
        // This example follows the same *atomic-units* canonical TST scaling used in MESS:
        //
        //   k(T) = (T / 2π) * W^‡(T) / W_AB(T)    [a0^3 / t0]
        //
        // where:
        // - T is k_B*T in Hartree (atomic units),
        // - W^‡ is the barrier/TS weight (here: PST weight),
        // - W_AB is the bimolecular reactant weight including relative translation
        //   (per unit volume) and the two internal weights.
        //
        // This avoids mixing SI Eyring prefactors with a.u. PST weights, which can
        // easily produce nonsense magnitudes.
        let reactant_weight_bimolecular_per_volume = bimolecular_weight_per_volume_atomic_units(
            28.0,
            17.0,
            q_int_co,
            q_int_oh,
            temperature_au,
        );
        let k_capture_au = (temperature_au / (2.0 * std::f64::consts::PI))
            * (q_pst / reactant_weight_bimolecular_per_volume);
        let k_capture = k_capture_au * AU_VOLUME_PER_TIME_CM3_S;

        let number_density_cm3 = number_density_cm3(pressure_pa, t);
        let pseudo_first_order = k_capture * number_density_cm3;

        let k_coll =
            hard_sphere_collision_limit_cm3_molecule_s(collision_diameter_angstrom, 28.0, 17.0, t);

        println!(
            "{:<6.0} {:>10.3e} {:>10.3e} {:>10.3e} {:>22.3e} {:>16.3e} {:>22.3e}",
            t, q_int_co, q_int_oh, q_pst, k_capture, pseudo_first_order, k_coll
        );

        if k_capture.is_finite() && k_coll.is_finite() && k_capture > k_coll {
            eprintln!(
                "WARNING: at T={:.0} K, k_cap ({:.3e}) exceeds hard-sphere collision limit ({:.3e}). \
                 This usually indicates inconsistent long-range parameters (V0,n) and/or an effective collision diameter that is too small.",
                t, k_capture, k_coll
            );
        }
    }

    Ok(())
}

fn reduced_mass_amu(m1: f64, m2: f64) -> f64 {
    m1 * m2 / (m1 + m2)
}

fn diatomic_brot_cm1(m1_amu: f64, m2_amu: f64, bond_angstrom: f64) -> f64 {
    let mu = reduced_mass_amu(m1_amu, m2_amu);
    // B (cm^-1) = h / (8 pi^2 c I), with I in amu*Å^2 and the constant folded below.
    let inertia_amu_ang2 = mu * bond_angstrom * bond_angstrom;
    16.857629 / inertia_amu_ang2
}

fn k_b_t_hartree(temperature_kelvin: f64) -> f64 {
    // k_B in cm^-1/K (same convention used in the ME modules), converted to Hartree.
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
    // MESS-like bimolecular reactant weight per unit volume (atomic units):
    //
    //   W_AB(T) = (μ T / 2π)^(3/2) * W_A,int(T) * W_B,int(T)
    //
    // where T is k_B*T in Hartree and μ is the reduced mass in electron-mass units.
    let mass_a_au = mass_a_amu * AMU_TO_ELECTRON_MASS;
    let mass_b_au = mass_b_amu * AMU_TO_ELECTRON_MASS;
    let reduced_mass_au = (mass_a_au * mass_b_au) / (mass_a_au + mass_b_au);

    let translational_prefactor_per_volume =
        (reduced_mass_au * temperature_au / (2.0 * std::f64::consts::PI)).powf(1.5);
    translational_prefactor_per_volume * internal_weight_a * internal_weight_b
}

fn number_density_cm3(pressure_pa: f64, temperature_kelvin: f64) -> f64 {
    // ideal gas: n = P/(k_B T) [m^-3] -> [cm^-3]
    const BOLTZMANN_SI: f64 = 1.380_649e-23; // J/K
    (pressure_pa / (BOLTZMANN_SI * temperature_kelvin)) / 1.0e6
}

fn hard_sphere_collision_limit_cm3_molecule_s(
    collision_diameter_angstrom: f64,
    mass_a_amu: f64,
    mass_b_amu: f64,
    temperature_kelvin: f64,
) -> f64 {
    // Hard-sphere collision-limit estimate:
    //   k_coll = sigma * <v_rel>
    // where sigma = pi d^2 and <v_rel> = sqrt(8 k_B T / (pi mu)).
    //
    // Units: returns cm^3 molecule^-1 s^-1.
    const BOLTZMANN_SI: f64 = 1.380_649e-23; // J/K
    const AMU_TO_KG: f64 = 1.660_539_066_60e-27;

    let d_m = collision_diameter_angstrom * 1.0e-10;
    let sigma_m2 = std::f64::consts::PI * d_m * d_m;

    let m1 = mass_a_amu * AMU_TO_KG;
    let m2 = mass_b_amu * AMU_TO_KG;
    let mu = (m1 * m2) / (m1 + m2);

    let v_rel = (8.0 * BOLTZMANN_SI * temperature_kelvin / (std::f64::consts::PI * mu)).sqrt();
    let k_m3_s = sigma_m2 * v_rel;
    k_m3_s * 1.0e6
}
