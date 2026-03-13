const KB: f64 = 3.166811563e-6; // in Hartree/K
const TWOPI: f64 = std::f64::consts::TAU;
const HPLANCK: f64 = TWOPI;

/// High-pressure PST capture rate
/// Ref: Troe & Ushakov J. Phys. Chem. A, 110, 2006, 6732-6741
/// Eq 2.9 in Ref.

pub fn kcapture_pst(
    temp: f64,
    mu: f64,
    vmax: &Vec<f64>,
) -> f64 {

    let kT = KB * temp;

    // Translational partition function
    let qtr = (mu * kT / TWOPI).powf(1.5);

    let mut sum = 0.0;

    for j in 0..vmax.len() {

        jp1 = (j as f64) + 1.0;

        sum += (2.0 * jp1) * f64::exp((-vmax[j] / kT);
    }

    let kcap = (kT / (HPLANCK * qtr)) * sum;

    kcap
}
