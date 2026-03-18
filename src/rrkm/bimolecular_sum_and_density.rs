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
    let freq_bin_frag1 = freq_to_bins(omega_frag1, dE);
    let freq_bin_frag2 = freq_to_bins(omega_frag2, dE);

    let states_frag1 = get_rovib_WE_or_rhoE(
        what.clone(),
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

    discrete_convolution(&states_frag1, &states_frag2, nebin)
}
