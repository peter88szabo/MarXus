/// Special functions used across numerical physics/chemistry modules.
///
/// This module intentionally has no external dependencies to keep MarXus easy to build.
///
/// Currently provided:
/// - `ln_gamma(z)`: natural logarithm of the Gamma function Γ(z)
/// - `gamma(z)`: Gamma function Γ(z)
///
/// Implementation notes:
/// - Uses a standard Lanczos approximation for ln(Γ) with reflection formula for z < 0.5.
/// - Accuracy is sufficient for smooth physics prefactors and partition-function style work.
/// - This is *not* intended for extreme-precision special-function applications.
use std::f64::consts::PI;

const LANCZOS_G: f64 = 7.0;

// Coefficients for g=7, n=9 (common choice).
const LANCZOS_COEFFS: [f64; 9] = [
    0.999_999_999_999_809_93,
    676.520_368_121_885_1,
    -1259.139_216_722_402_8,
    771.323_428_777_653_13,
    -176.615_029_162_140_59,
    12.507_343_278_686_905,
    -0.138_571_095_265_720_12,
    9.984_369_578_019_572e-6,
    1.505_632_735_149_311_6e-7,
];

/// Natural logarithm of the Gamma function, ln(Γ(z)).
pub fn ln_gamma(z: f64) -> Result<f64, String> {
    if !z.is_finite() {
        return Err("ln_gamma: input must be finite.".into());
    }
    if z <= 0.0 && (z.fract() == 0.0) {
        return Err("ln_gamma: pole at non-positive integer.".into());
    }

    // Reflection formula to handle small z:
    //   Γ(z) Γ(1-z) = π / sin(π z)
    // => ln Γ(z) = ln π - ln sin(π z) - ln Γ(1-z)
    if z < 0.5 {
        let sin_term = (PI * z).sin();
        if sin_term == 0.0 {
            return Err("ln_gamma: sin(pi*z)=0 (pole).".into());
        }
        let lg = ln_gamma(1.0 - z)?;
        return Ok(PI.ln() - sin_term.abs().ln() - lg);
    }

    // Lanczos approximation for z >= 0.5
    let z_minus_1 = z - 1.0;
    let mut x = LANCZOS_COEFFS[0];
    for (i, c) in LANCZOS_COEFFS.iter().enumerate().skip(1) {
        x += c / (z_minus_1 + (i as f64));
    }

    let t = z_minus_1 + LANCZOS_G + 0.5;
    let half_ln_2pi = 0.5 * (2.0 * PI).ln();
    Ok(half_ln_2pi + (z_minus_1 + 0.5) * t.ln() - t + x.ln())
}

/// Gamma function Γ(z).
pub fn gamma(z: f64) -> Result<f64, String> {
    Ok(ln_gamma(z)?.exp())
}
