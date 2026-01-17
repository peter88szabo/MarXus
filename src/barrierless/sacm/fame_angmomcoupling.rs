use crate::numeric::lanczos_gamma::gamma_func;

// This function calculates angular momentum coupling factors Fame(E,J).
// Ref.: JCP 79, 6017 (1983), eqs. A17, C1, C11, prolate symmetric top
// expression (C10) modified after Ber. Bunsenges. Phys. Chem. 98, 1563
// (1994).
// 
// KZR1: parameter deciding the case (linear or non-linear molecule)
// ene: energy
// jtot: rotational quantum number
// sstar: spin angular momentum
// Brot_A: rotational constant for fragment A
// Brot_B: rotational constant for fragment B
// Returns the calculated Fame value.
pub fn factor_angmom(kzr1: i32, ene: f64, jot: f64, mut sstar: f64, rot_a: f64, rot_b: f64) -> f64 {
    let qgam = gamma_func(sstar) * gamma_func(1.5) / gamma_func(sstar + 0.5);
    let arg = rot_a / ene;
    let argp = rot_a / (rot_a - rot_b);

    match kzr1 {
        3 => {
            // Linear molecule case
            sstar += 1.0;  // SSTAR adjustment
            
            let f1 = gamma_func(sstar + 1.0) / gamma_func(sstar - 1.0)
                * (1.0 + jot) * (1.0 + jot / 2.0) * arg * arg;
            
            f1 / (f1 + 1.0)
        }
        4 => {
            // Fame calculation for KZR1 = 4
            (2.0 * jot + 1.0) * arg.sqrt() / qgam
        }
        5 => {
            // Fame calculation for KZR1 = 5
            if jot == 0.0 {
                // Special case where JOT == 0
                return arg.sqrt() / qgam;
            }

            let x = (rot_a - rot_b) * (jot + 1.0) * (jot + 1.0) / ene;
            let arg_tanh = (1.0 / (2.0 * qgam)) * (1.0 + (2.0 + sstar) / 10.0 * x) * x.sqrt() * 2.0;
            argp.sqrt() * arg_tanh.tanh()
        }
        _ => panic!("Invalid value for KZR1"),  // Handle unexpected KZR1 values
    }
}
