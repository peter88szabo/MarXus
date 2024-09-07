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
    dH0: f64) -> Vec<f64> {

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
