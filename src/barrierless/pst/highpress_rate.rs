use crate::constants::{HPLANCK_AU, KB_AU_PER_K, TWO_PI};

/// High-pressure PST capture rate
/// Ref: Troe & Ushakov J. Phys. Chem. A, 110, 2006, 6732-6741
/// Eq 2.9 in Ref.

pub fn kcapture_pst_highpress(temp: f64, mu: f64, vmax: &Vec<f64>) -> f64 {
    let kT = KB_AU_PER_K * temp;

    // Translational partition function
    let qtr = (mu * kT / TWO_PI).powf(1.5);

    let mut sum = 0.0;

    for j in 0..vmax.len() {
        sum += (2.0 * ((j + 1) as f64)) * f64::exp((-vmax[j] / kT));
    }

    let kcap = (kT / (HPLANCK_AU * qtr)) * sum;

    kcap
}
