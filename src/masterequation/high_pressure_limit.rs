use crate::molecule::MoleculeStruct;
use crate::tunneling::tunneling::{eckart, wigner};

const KB_OVER_H_PER_KELVIN_PER_SECOND: f64 = 2.083_661_912e10;
const BOLTZMANN_SI: f64 = 1.380_649e-23; // J/K
const PLANCK_SI: f64 = 6.626_070_15e-34; // J*s
const KCAL_PER_HARTREE: f64 = 627.51;
const R_KCAL_PER_MOL_PER_KELVIN: f64 = 1.987_204_258_640_83e-3;
const ATM_TO_PASCAL: f64 = 101_325.0;
const CM1_TO_HARTREE: f64 = 4.55635e-6;
const RGAS_AU: f64 = 8.31446261815324 / 1000.0 / 2625.5; // Hartree/mol/K

#[derive(Clone, Copy, Debug)]
pub enum ReactionMolecularity {
    Unimolecular,
    Bimolecular,
}

/// Input for simple high-pressure unimolecular TST/Eyring:
///
/// k_inf(T) = kappa * (k_B T / h) * exp( -ΔG‡ / (R T) )
#[derive(Clone, Copy, Debug)]
pub struct HighPressureTstInput {
    pub temperature_kelvin: f64,
    pub delta_g_dagger_kcal_mol: f64,
    pub tunneling_kappa: f64,
    /// Needed for bimolecular TST; ignored for unimolecular.
    pub pressure_atm: f64,
    pub molecularity: ReactionMolecularity,
}

/// Output bundle for high-pressure TST from free energies.
#[derive(Clone, Copy, Debug)]
pub struct HighPressureTstResult {
    pub rate_constant_s_inv: f64,
    pub prefactor_s_inv: f64,
    pub exponential_factor: f64,
}

/// Compute high-pressure unimolecular TST rate using a free-energy barrier in kcal/mol.
pub fn high_pressure_unimolecular_tst_rate(
    input: HighPressureTstInput,
) -> Result<HighPressureTstResult, String> {
    if input.temperature_kelvin <= 0.0 {
        return Err("temperature_kelvin must be positive".into());
    }
    if input.tunneling_kappa <= 0.0 || !input.tunneling_kappa.is_finite() {
        return Err("tunneling_kappa must be positive and finite".into());
    }

    let prefactor_base =
        input.tunneling_kappa * KB_OVER_H_PER_KELVIN_PER_SECOND * input.temperature_kelvin;
    let vmolar_factor = match input.molecularity {
        ReactionMolecularity::Unimolecular => 1.0,
        ReactionMolecularity::Bimolecular => {
            if input.pressure_atm <= 0.0 {
                return Err("pressure_atm must be positive for bimolecular TST".into());
            }
            // Vmolar factor for bimolecular steps:
            // Vmolar = (k_B T / P) converted from m^3/molecule to cm^3/molecule.
            let pressure_pa = input.pressure_atm * ATM_TO_PASCAL;
            (BOLTZMANN_SI * input.temperature_kelvin / pressure_pa) * 1.0e6
        }
    };
    let prefactor = prefactor_base * vmolar_factor;
    let exponent =
        -input.delta_g_dagger_kcal_mol / (R_KCAL_PER_MOL_PER_KELVIN * input.temperature_kelvin);
    let exponential_factor = exponent.exp();
    let k = prefactor * exponential_factor;

    if !k.is_finite() || k < 0.0 {
        return Err("Computed high-pressure rate is invalid".into());
    }

    Ok(HighPressureTstResult {
        rate_constant_s_inv: k,
        prefactor_s_inv: prefactor,
        exponential_factor,
    })
}

/// Output from computing high-pressure TST using `MoleculeStruct` thermochemistry.
#[derive(Clone, Copy, Debug)]
pub struct HighPressureFromMoleculesResult {
    pub g_reactant_kcal_mol: f64,
    pub g_ts_kcal_mol: f64,
    pub delta_g_dagger_kcal_mol: f64,
    pub rate_constant_s_inv: f64,
}

/// Reaction-level thermodynamic output:
/// dG = G_TS - sum(G_reactants), dZPE, dH0, and rate prefactor pieces.
#[derive(Clone, Copy, Debug)]
pub struct HighPressureThermoResult {
    pub d_g_kcal_mol: f64,
    pub d_zpe_kcal_mol: f64,
    pub d_h0_kcal_mol: f64,
    pub prefactor: f64,
    pub rate_constant: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct EckartTunnelingInput {
    pub temperature_kelvin: f64,
    pub imaginary_frequency_cm1: f64,
    pub forward_barrier_kcal_mol: f64,
    pub reverse_barrier_kcal_mol: f64,
    pub integration_step_kcal_mol: f64,
    pub integration_max_kcal_mol: f64,
}

/// Evaluate thermochemistry on reactant/TS and compute k_inf via TST.
///
/// Notes:
/// - Pressure is provided in atm and converted to Pa for `eval_all_therm_func`.
/// - Both species are evaluated at the same T and P.
pub fn high_pressure_unimolecular_tst_from_molecules(
    reactant: &mut MoleculeStruct,
    transition_state: &mut MoleculeStruct,
    temperature_kelvin: f64,
    pressure_atm: f64,
    freq_cutoff_cm1: f64,
    tunneling_kappa: f64,
) -> Result<HighPressureFromMoleculesResult, String> {
    if pressure_atm <= 0.0 {
        return Err("pressure_atm must be positive".into());
    }

    let pressure_pa = pressure_atm * ATM_TO_PASCAL;

    reactant.eval_all_therm_func(temperature_kelvin, pressure_pa, freq_cutoff_cm1);
    transition_state.eval_all_therm_func(temperature_kelvin, pressure_pa, freq_cutoff_cm1);

    let g_reactant_kcal_mol = reactant.thermo.gtot * KCAL_PER_HARTREE;
    let g_ts_kcal_mol = transition_state.thermo.gtot * KCAL_PER_HARTREE;
    let delta_g_dagger_kcal_mol = g_ts_kcal_mol - g_reactant_kcal_mol;

    let rate = high_pressure_unimolecular_tst_rate(HighPressureTstInput {
        temperature_kelvin,
        delta_g_dagger_kcal_mol,
        tunneling_kappa,
        pressure_atm,
        molecularity: ReactionMolecularity::Unimolecular,
    })?;

    Ok(HighPressureFromMoleculesResult {
        g_reactant_kcal_mol,
        g_ts_kcal_mol,
        delta_g_dagger_kcal_mol,
        rate_constant_s_inv: rate.rate_constant_s_inv,
    })
}

/// TST rate for A (+B) -> TS:
///   dG  = G_TS - (G_A + G_B)
///   dZPE = ZPE_TS - (ZPE_A + ZPE_B)
///   dH0 = Eelec_TS - (Eelec_A + Eelec_B) + dZPE
///   rate = kappa * prefac * exp(-dG/(R*T))
///
/// `reactant_b` is `None` for unimolecular reactions.
pub fn high_pressure_tst_with_thermo(
    reactant_a: &mut MoleculeStruct,
    mut reactant_b: Option<&mut MoleculeStruct>,
    transition_state: &mut MoleculeStruct,
    temperature_kelvin: f64,
    pressure_atm: f64,
    freq_cutoff_cm1: f64,
    tunneling_kappa: f64,
    molecularity: ReactionMolecularity,
) -> Result<HighPressureThermoResult, String> {
    if temperature_kelvin <= 0.0 {
        return Err("temperature_kelvin must be positive".into());
    }
    if pressure_atm <= 0.0 {
        return Err("pressure_atm must be positive".into());
    }

    let pressure_pa = pressure_atm * ATM_TO_PASCAL;
    reactant_a.eval_all_therm_func(temperature_kelvin, pressure_pa, freq_cutoff_cm1);
    transition_state.eval_all_therm_func(temperature_kelvin, pressure_pa, freq_cutoff_cm1);
    if let Some(rb) = reactant_b.as_deref_mut() {
        rb.eval_all_therm_func(temperature_kelvin, pressure_pa, freq_cutoff_cm1);
    }

    let g_a = reactant_a.thermo.gtot;
    let zpe_a = reactant_a.zpe * CM1_TO_HARTREE;
    let eelec_a = reactant_a.ene * CM1_TO_HARTREE;

    let (g_b, zpe_b, eelec_b) = match reactant_b {
        Some(rb) => (
            rb.thermo.gtot,
            rb.zpe * CM1_TO_HARTREE,
            rb.ene * CM1_TO_HARTREE,
        ),
        None => (0.0, 0.0, 0.0),
    };

    let g_ts = transition_state.thermo.gtot;
    let zpe_ts = transition_state.zpe * CM1_TO_HARTREE;
    let eelec_ts = transition_state.ene * CM1_TO_HARTREE;

    let d_g = g_ts - (g_a + g_b);
    let d_zpe = zpe_ts - (zpe_a + zpe_b);
    let d_h0 = eelec_ts - (eelec_a + eelec_b) + d_zpe;

    let vmolar = match molecularity {
        ReactionMolecularity::Unimolecular => 1.0,
        ReactionMolecularity::Bimolecular => {
            (BOLTZMANN_SI * temperature_kelvin / pressure_pa) * 1.0e6
        }
    };
    let prefactor = tunneling_kappa * (BOLTZMANN_SI / PLANCK_SI) * temperature_kelvin * vmolar;
    let d_g_kcal = d_g * KCAL_PER_HARTREE;
    let rate = prefactor * (-(d_g_kcal) / (R_KCAL_PER_MOL_PER_KELVIN * temperature_kelvin)).exp();

    Ok(HighPressureThermoResult {
        d_g_kcal_mol: d_g_kcal,
        d_zpe_kcal_mol: d_zpe * KCAL_PER_HARTREE,
        d_h0_kcal_mol: d_h0 * KCAL_PER_HARTREE,
        prefactor,
        rate_constant: rate,
    })
}

/// Compute Eckart tunneling correction kappa using the existing tunneling module.
pub fn eckart_tunneling_kappa(input: EckartTunnelingInput) -> Result<f64, String> {
    if input.temperature_kelvin <= 0.0 {
        return Err("temperature_kelvin must be positive".into());
    }
    if input.imaginary_frequency_cm1 <= 0.0 {
        return Err("imaginary_frequency_cm1 must be positive".into());
    }
    if input.forward_barrier_kcal_mol <= 0.0 || input.reverse_barrier_kcal_mol <= 0.0 {
        return Err("Eckart barriers must be positive".into());
    }
    if input.integration_step_kcal_mol <= 0.0 || input.integration_max_kcal_mol <= 0.0 {
        return Err("Eckart integration settings must be positive".into());
    }

    let beta = 1.0 / (RGAS_AU * input.temperature_kelvin);
    let omega = input.imaginary_frequency_cm1 * CM1_TO_HARTREE;
    let vf = input.forward_barrier_kcal_mol / KCAL_PER_HARTREE;
    let vb = input.reverse_barrier_kcal_mol / KCAL_PER_HARTREE;
    let de = input.integration_step_kcal_mol / KCAL_PER_HARTREE;
    let emax = input.integration_max_kcal_mol / KCAL_PER_HARTREE;

    let (kappa1, _kappa2) = eckart(beta, omega, vf, vb, de, emax);
    if kappa1.is_finite() && kappa1 > 0.0 && kappa1 < 1.0e4 {
        return Ok(kappa1);
    }

    // Guardrail: if Eckart is numerically unstable in this parameterization,
    // use Wigner as a controlled near-barrier approximation.
    let kappa_wigner = wigner(beta, omega);
    if !kappa_wigner.is_finite() || kappa_wigner <= 0.0 {
        return Err("Both Eckart and Wigner tunneling corrections are invalid.".into());
    }
    Ok(kappa_wigner)
}
