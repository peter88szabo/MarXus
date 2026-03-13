// Morse potential function (minimum = 0, asymptote = D_e)

// src/barrierless/pst/morsepot.rs
// Morse potential (minimum = 0, asymptote = D_e)

/// Returns (V(r), dV/dr) for the Morse potential:
/// V(r)  = D_e * (1 - exp(-beta*(r - r_e)))^2
/// dV/dr = 2 * D_e * beta * exp(-beta*(r - r_e)) * (1 - exp(-beta*(r - r_e)))
pub fn morse_value_and_derivative(r: f64, de: f64, beta: f64, re: f64) -> (f64, f64) {
    let x = (-beta * (r - re)).exp();      // x = exp(-beta*(r-re))
    let one_minus_x = 1.0 - x;

    let v = de * one_minus_x * one_minus_x;
    let dv = 2.0 * de * beta * x * one_minus_x;

    (v, dv)
}
