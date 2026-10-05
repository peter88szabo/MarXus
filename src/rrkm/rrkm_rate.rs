#![allow(non_snake_case)]
use crate::constants::H_PLANCK_CM;
use crate::rrkm::sum_and_density::get_rovib_WE_or_rhoE;

//compute the RRKM formula: k(E) = sigma * W_ts(E)/rho(E) / hplanc
// where
// W_ts: the sum of states at the TS
// rho:  density of states for reactants (the complex to dissociate), states per cm-1
// sigma = sigma_cpx / sigma_ts is the reaction path degeneracy

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

    // Reaction path degeneracy alpha = sigma_cpx / sigma_ts (Forst, Theory of Unimolecular
    // Reactions, 1973, Ch. 4, Sec. 5). This ratio is only a special case: when the transition
    // state has several paths back to the reactant, or only one of the two has a symmetry element
    // other than a rotation (Schlag), pass sigma_cpx = alpha and sigma_ts = 1 instead.
    let sigma = sigma_cpx / sigma_ts;

    //Minimum energy (including ZPE) of the reaction as integer energy bin
    let nbin_dH0 = (dH0 / dE + 0.5) as usize; //rate is compuated from the top of the barrier

    // The RRKM formula: k(E) = sigma * W_ts(E)/rho(E) / hplanck
    let mut kE = vec![0.0; nebin + 1];

    for i in nbin_dH0..=nebin {
        kE[i] = sigma * WE_ts[i - nbin_dH0] / rhoE_cpx[i] / H_PLANCK_CM;
    }

    return kE;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rate_scales_with_reaction_path_degeneracy_sigma_over_sigma_ts() {
        // Reaction path degeneracy alpha = sigma / sigma_ts (reactant over transition state),
        // Forst, Theory of Unimolecular Reactions (1973), Ch. 4, Sec. 5.
        let (nebin, d_e, dh0) = (1000, 10.0, 2000.0);
        let omega_ts = [500.0, 1200.0];
        let omega_cpx = [300.0, 900.0, 1500.0];
        let brot = [1.0, 0.5, 0.25];
        let k_11 = get_kE(nebin, d_e, 2, 3, &omega_ts, &omega_cpx, 3, 3, &brot, &brot, 1.0, 1.0, dh0);
        let k_12 = get_kE(nebin, d_e, 2, 3, &omega_ts, &omega_cpx, 3, 3, &brot, &brot, 1.0, 2.0, dh0);
        let ratio = k_12[nebin] / k_11[nebin];
        assert!((ratio - 2.0).abs() < 1e-12, "k(sigma_cpx=2)/k(sigma_cpx=1) = {ratio}");
    }
}
