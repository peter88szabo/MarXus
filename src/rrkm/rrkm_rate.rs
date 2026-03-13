#![allow(non_snake_case)]
use crate::rrkm::sum_and_density::get_rovib_WE_or_rhoE;

//compute the RRKM formula: k(E) = sigma * W_ts(E)/rho(E) / hplanc
// where
// W_ts: the sum of states at the TS
// rho:  density of states for reactants (the complex to dissociate)
// sigma is the symmetry number

pub fn get_kE(
    nebin: usize,
    dE: f64,
    nvib_ts: usize,
    nvib_cpx: usize,
    omega_ts: &[f64],
    omega_cpx: &[f64],
    nrot_ts: usize,
    nrot_cpx: usize,
    Brot_ts: &[f64],
    Brot_cpx: &[f64],
    sigma_ts: f64,
    sigma_cpx: f64,
    dH0: f64,
) -> Vec<f64> {
    let mut freq_bin_cpx = vec![0; nvib_cpx];
    for i in 0..nvib_cpx {
        freq_bin_cpx[i] = (omega_cpx[i] / dE + 0.5) as usize;
    }

    let mut freq_bin_ts = vec![0; nvib_ts];

    for i in 0..nvib_ts {
        freq_bin_ts[i] = (omega_ts[i] / dE + 0.5) as usize;
    }

    let WE_ts = get_rovib_WE_or_rhoE(
        "sum".to_string(),
        nvib_ts,
        nebin,
        dE,
        nrot_ts,
        &freq_bin_ts,
        &Brot_ts,
    );

    let rhoE_cpx = get_rovib_WE_or_rhoE(
        "den".to_string(),
        nvib_cpx,
        nebin,
        dE,
        nrot_cpx,
        &freq_bin_cpx,
        &Brot_cpx,
    );

    const H_PLANCK: f64 = 3.3356E-11;
    let sigma = sigma_ts / sigma_cpx;

    //Minimum energy (including ZPE) of the reaction as integer energy bin
    let nbin_dH0 = (dH0 / dE + 0.5) as usize; //rate is compuated from the top of the barrier

    // The RRKM formula: k(E) = sigma * W_ts(E)/rho(E) / hplanck
    let mut kE = vec![0.0; nebin + 1];

    for i in nbin_dH0..=nebin {
        kE[i] = sigma * WE_ts[i - nbin_dH0] / rhoE_cpx[i] / H_PLANCK;
    }

    return kE;
}

// RRKM microcanonical rate from precomputed TS sum-of-states and well density-of-states.
//
// This is a small helper to avoid re-implementing the same RRKM bin-shift logic in other modules.
// Conventions match `get_kE`:
// - Energies on a uniform grid with bin width `dE` (cm^-1)
// - `ts_sum_states` is W‡(E) on bins (dimensionless cumulative states)
// - `well_density_of_states` is ρ(E) on bins (states per cm^-1)
// - `threshold_energy_cm1` is E0 (cm^-1); k(E)=0 for E<E0
//
// k(E) [1/s] = W‡(E - E0) / ρ(E) / H_PLANCK,
// where H_PLANCK is in s/cm so 1/H_PLANCK = c (cm/s).
pub fn rrkm_rate_from_sum_states_and_density(
    dE: f64,
    ts_sum_states: &[f64],
    well_density_of_states: &[f64],
    threshold_energy_cm1: f64,
) -> Result<Vec<f64>, String> {
    if dE <= 0.0 || !dE.is_finite() {
        return Err("dE must be positive and finite.".into());
    }
    if threshold_energy_cm1 < 0.0 || !threshold_energy_cm1.is_finite() {
        return Err("threshold_energy_cm1 must be non-negative and finite.".into());
    }
    if ts_sum_states.len() != well_density_of_states.len() {
        return Err("ts_sum_states and well_density_of_states must have the same length.".into());
    }

    const H_PLANCK: f64 = 3.3356E-11; // s/cm

    let shift = (threshold_energy_cm1 / dE + 0.5) as usize;
    let mut kE = vec![0.0; ts_sum_states.len()];

    for i in shift..ts_sum_states.len() {
        let rho = well_density_of_states[i];
        if rho <= 0.0 || !rho.is_finite() {
            continue;
        }
        let w = ts_sum_states[i - shift];
        if w <= 0.0 || !w.is_finite() {
            continue;
        }
        kE[i] = w / rho / H_PLANCK;
    }

    Ok(kE)
}
