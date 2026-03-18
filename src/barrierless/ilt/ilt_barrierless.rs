use crate::constants::{H_PLANCK_CM, KB_CM};
use crate::numeric::lanczos_gamma::gamma_func;
use crate::rrkm::bimolecular_sum_and_density::bimol_get_rovib_WE_or_rhoE;
use crate::rrkm::sum_and_density::get_rovib_WE_or_rhoE;

// All energies are in cm-1.
// Output W(E) is the RRKM-compatible sum of states.

fn freq_to_bins(freqs: &[f64], dE: f64) -> Vec<usize> {
    let mut bins = vec![0; freqs.len()];
    for (i, omega) in freqs.iter().enumerate() {
        bins[i] = (omega / dE + 0.5) as usize;
    }
    bins
}

fn density_value_at_energy(densitys: &[f64], energy: f64, dE: f64) -> f64 {
    if energy < 0.0 {
        return 0.0;
    }

    let idx = (energy / dE + 0.5) as usize;
    densitys.get(idx).copied().unwrap_or(0.0)
}

fn unimolecular_rho(
    nebin: usize,
    dE: f64,
    nvib: usize,
    nrot: usize,
    omega: &[f64],
    Brot: &[f64],
) -> Vec<f64> {
    let freq_bin = freq_to_bins(omega, dE);
    get_rovib_WE_or_rhoE("den".to_string(), nvib, nebin, dE, nrot, &freq_bin, Brot)
}

fn w_ilt_unimol_at_energy(
    e: f64,
    ea: f64,
    dE: f64,
    n: f64,
    prefactor_krho: f64,
    rho: &[f64],
) -> f64 {
    if e <= ea {
        return 0.0;
    }

    let n_steps = 4000usize;
    let h = (e - ea) / n_steps as f64;
    let mut sum = 0.0;

    for i in 0..=n_steps {
        let x = i as f64 * h;
        let tau = ea + x;
        let weight = if i == 0 || i == n_steps { 0.5 } else { 1.0 };

        let kernel = if x == 0.0 && n > 1.0 {
            0.0
        } else {
            x.powf(n - 1.0)
        };

        let val = density_value_at_energy(rho, e - tau, dE) * kernel;
        sum += weight * val;
    }

    // ILT integral gives k(E)*rho(E); convert to RRKM sum of states W(E)
    H_PLANCK_CM * prefactor_krho * sum * h
}

fn w_ilt_bimol_at_energy(
    e: f64,
    eth: f64,
    dE: f64,
    n: f64,
    prefactor_krho: f64,
    rho_r: &[f64],
) -> f64 {
    if e <= eth {
        return 0.0;
    }

    let n_steps = 4000usize;
    let h = (e - eth) / n_steps as f64;
    let mut sum = 0.0;

    for i in 0..=n_steps {
        let x = i as f64 * h;
        let tau = eth + x;
        let weight = if i == 0 || i == n_steps { 0.5 } else { 1.0 };
        let kernel = x.powf(n + 0.5);
        let val = density_value_at_energy(rho_r, e - tau, dE) * kernel;

        sum += weight * val;
    }

    // ILT integral gives k(E)*rho(E); convert to RRKM sum of states W(E)
    H_PLANCK_CM * prefactor_krho * sum * h
}

/// Unimolecular barrierless ILT expression
///
/// First forms k(E) * rho(E) by ILT, then converts it to the RRKM-compatible
/// sum of states W(E) using W(E) = H_PLANCK * k(E) * rho(E).
///
/// Formula:
/// k(E) * rho(E) = (A * beta0^n / Gamma(n))
///                * ∫_{Ea}^{E} rho(E - tau) * (tau - Ea)^(n-1) d tau
pub fn W_ilt_unimol(
    a: f64,
    ea: f64,
    t0: f64,
    n: f64,
    nebin: usize,
    dE: f64,
    nvib: usize,
    nrot: usize,
    omega: &[f64],
    Brot: &[f64],
) -> Vec<f64> {
    let rho = unimolecular_rho(nebin, dE, nvib, nrot, omega, Brot);
    let beta0 = 1.0 / (KB_CM * t0);
    let prefactor_krho = a * beta0.powf(n) / gamma_func(n);
    let mut W_e = vec![0.0; nebin + 1];

    for (i, w) in W_e.iter_mut().enumerate() {
        let e = (i as f64) * dE;
        *w = w_ilt_unimol_at_energy(e, ea, dE, n, prefactor_krho, &rho);
    }

    W_e
}

/// Bimolecular barrierless ILT expression
///
/// First forms k(E) * rho_r(E) by ILT, then converts it to the RRKM-compatible
/// sum of states W(E) using W(E) = H_PLANCK * k(E) * rho_r(E).
///
/// Formula:
/// k(E) * rho_r(E) = bimol_prefactor_cm
///                 * (A * beta0^n / Gamma(n + 1.5))
///                 * ∫_{Ea + delta_h0}^{E}
///                     rho_r(E - tau) * (tau - Ea - delta_h0)^(n + 0.5) d tau
pub fn W_ilt_bimol(
    a: f64,
    ea: f64,
    delta_h0: f64,
    t0: f64,
    n: f64,
    bimol_prefactor_cm: f64,
    nebin: usize,
    dE: f64,
    nvib_frag1: usize,
    nrot_frag1: usize,
    omega_frag1: &[f64],
    Brot_frag1: &[f64],
    nvib_frag2: usize,
    nrot_frag2: usize,
    omega_frag2: &[f64],
    Brot_frag2: &[f64],
) -> Vec<f64> {
    let eth = ea + delta_h0;

    let rho_r = bimol_get_rovib_WE_or_rhoE(
        "den".to_string(),
        nebin,
        dE,
        nvib_frag1,
        nrot_frag1,
        omega_frag1,
        Brot_frag1,
        nvib_frag2,
        nrot_frag2,
        omega_frag2,
        Brot_frag2,
    );

    let beta0 = 1.0 / (KB_CM * t0);
    let prefactor_krho = bimol_prefactor_cm * a * beta0.powf(n) / gamma_func(n + 1.5);
    let mut W_e = vec![0.0; nebin + 1];

    for (i, w) in W_e.iter_mut().enumerate() {
        let e = (i as f64) * dE;
        *w = w_ilt_bimol_at_energy(e, eth, dE, n, prefactor_krho, &rho_r);
    }

    W_e
}

#[cfg(test)]
mod tests {
    use super::{W_ilt_bimol, W_ilt_unimol};

    #[test]
    fn smoke_test_ilt_calls() {
        let w_ilt_unimol_vec = W_ilt_unimol(
            1.0e12,
            100.0,
            300.0,
            2.0,
            200,
            10.0,
            2,
            3,
            &[500.0, 750.0],
            &[0.2, 0.1, 0.08],
        );

        println!("unimolecular ILT w(E)[100] = {}", w_ilt_unimol_vec[100]);
        println!("unimolecular ILT len = {}", w_ilt_unimol_vec.len());

        let w_ilt_bimol_vec = W_ilt_bimol(
            1.0e12,
            50.0,
            25.0,
            300.0,
            2.0,
            1.0,
            200,
            10.0,
            1,
            3,
            &[400.0],
            &[0.25, 0.12, 0.09],
            1,
            2,
            &[350.0],
            &[0.18, 0.11],
        );

        println!("bimolecular ILT w(E)[120] = {}", w_ilt_bimol_vec[120]);
        println!("bimolecular ILT len = {}", w_ilt_bimol_vec.len());
    }
}
