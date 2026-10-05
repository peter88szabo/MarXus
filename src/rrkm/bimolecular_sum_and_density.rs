use crate::rrkm::sum_and_density::get_rovib_WE_or_rhoE;

fn freq_to_bins(freqs: &[f64], dE: f64) -> Vec<usize> {
    let mut bins = vec![0; freqs.len()];
    for (i, omega) in freqs.iter().enumerate() {
        bins[i] = (omega / dE + 0.5) as usize;
    }
    bins
}

fn discrete_convolution(lhs: &[f64], rhs: &[f64], nebin: usize) -> Vec<f64> {
    let mut out = vec![0.0; nebin + 1];

    for i in 0..=nebin {
        let mut sum = 0.0;
        for j in 0..=i {
            sum += lhs[j] * rhs[i - j];
        }
        out[i] = sum;
    }

    out
}

pub fn bimol_get_rovib_WE_or_rhoE(
    what: String,
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
    // Internal (rovibrational) states of two independent fragments, combined by convolution
    // (Forst, Chem. Rev. 71, 339 (1971), eqs. 30 and 33; relative translation is not included):
    //   rho_12(E) = int_0^E rho_1(x) rho_2(E - x) dx,   W_12(E) = int_0^E rho_1(x) W_2(E - x) dx
    // Fragment 1 therefore always enters as a density (states per cm-1).
    let freq_bin_frag1 = freq_to_bins(omega_frag1, dE);
    let freq_bin_frag2 = freq_to_bins(omega_frag2, dE);

    let rho_frag1 = get_rovib_WE_or_rhoE(
        "den".to_string(),
        nvib_frag1,
        nebin,
        dE,
        nrot_frag1,
        &freq_bin_frag1,
        Brot_frag1,
    );

    let states_frag2 = get_rovib_WE_or_rhoE(
        what,
        nvib_frag2,
        nebin,
        dE,
        nrot_frag2,
        &freq_bin_frag2,
        Brot_frag2,
    );

    let mut out = discrete_convolution(&rho_frag1, &states_frag2, nebin);
    for v in out.iter_mut() {
        *v *= dE;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_oscillators_count_exactly() {
        // Two independent 1000 cm-1 oscillators: states with v1 + v2 <= 4 at E = 4000 cm-1 -> 15.
        for d_e in [10.0, 5.0] {
            let nebin = (5000.0 / d_e) as usize;
            let i4000 = (4000.0 / d_e) as usize;
            let w = bimol_get_rovib_WE_or_rhoE(
                "sum".to_string(), nebin, d_e, 1, 0, &[1000.0], &[], 1, 0, &[1000.0], &[],
            );
            assert!((w[i4000] - 15.0).abs() < 1e-9, "dE={d_e}: W(4000) = {}", w[i4000]);
        }
    }

    #[test]
    fn two_linear_rotors_follow_forst_eq43() {
        // Two 2D rotors: r = 4, Q'_4 = 1/(B1 B2), G(E) = Q'_4 E^2 / Gamma(3) = E^2 / (2 B1 B2)
        // (Forst, Chem. Rev. 71, 339 (1971), eqs. 33, 40, 43); grain error of order dE/E.
        let (b1, b2) = (1.5, 0.3);
        let d_e = 2.0; // dE != 1 so that a missing dE factor cannot go unnoticed
        let nebin = 2000;
        let e = nebin as f64 * d_e;
        let w = bimol_get_rovib_WE_or_rhoE(
            "sum".to_string(), nebin, d_e, 0, 2, &[], &[b1, b1], 0, 2, &[], &[b2, b2],
        );
        let rho = bimol_get_rovib_WE_or_rhoE(
            "den".to_string(), nebin, d_e, 0, 2, &[], &[b1, b1], 0, 2, &[], &[b2, b2],
        );
        let w_ref = e * e / (2.0 * b1 * b2);
        let rho_ref = e / (b1 * b2);
        assert!(((w[nebin] - w_ref) / w_ref).abs() < 2e-3, "W = {} vs {}", w[nebin], w_ref);
        assert!(((rho[nebin] - rho_ref) / rho_ref).abs() < 2e-3, "rho = {} vs {}", rho[nebin], rho_ref);
    }
}
