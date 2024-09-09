use std::f64;
mod lanczos_gamma;
use lanczos_gamma::gamma_func;

// This function calculates the correction factor between the partition
// function of a hindered and a free internal rotor.
// 
// KT: Boltzmann constant times temperature (in energy units)
// V0: barrier height for the hindered rotor
// QFREE: partition function of the free internal rotor
// N: number of hindered rotors
// 
// Returns the correction factor FQHIND.
pub fn fqhind(kt: f64, v0: f64, qfree: f64, n: f64) -> f64 {
    if v0 != 0.0 {
        let ktv0 = kt / v0;
        let hnykt = n * (std::f64::consts::PI / ktv0).sqrt() / qfree;
        let term1 = (-1.2 * ktv0).exp() / (qfree * (1.0 - (-hnykt).exp()));
        let term2 = (1.0 - (-ktv0).exp()).powf(1.2);
        term1 + term2
    } else {
        1.0
    }
}


// This function calculates the correction factor Whind(E)/Wfree(E)
// for the sums of states of an ensemble of S oscillators and one
// hindered rotor with a hindrance potential V0.
//
// S: number of oscillators
// V0: hindrance potential
// E: energy
//
// Returns the correction factor RWHIND.
pub fn rwhind(s: i32, v0: f64, e: f64) -> f64 {
    const PI: f64 = std::f64::consts::PI;  // Using the value of pi
    let mut ic: i32 = 0;
    let mut gs: f64;
    let mut a: f64;
    let mut b = vec![0.0; s];  // Array to store B coefficients

    if v0 > 0.0 && e > 0.0 {
        ic += 1;
        if ic == 1 {
            gs = gamma_func(s as f64 + 1.5);
            a = gs / gamma_func(PI * (s as f64 + 2.0));

            // Calculate B coefficients
            for ny in 0..s {
                b[ny as usize] = (-1.0f64).powi(ny)
                    * gs
                    * (ny as f64 + 2.5)
                    / (gamma_func(ny as f64 + 1.0) * gamma_func(s as f64 - ny as f64) * PI * (ny as f64 + 2.0) * (ny as f64 + 1.5));
            }
        }

        let ve = v0 / e;
        if ve < 1.0 {
            let mut whind = 1.0;
            for i in 0..s {
                whind -= b[i as usize] * ve.powf(i as f64 + 1.5);
            }
            whind
        } else {
            let whind = a / ve.sqrt();
            whind
        }
    } else {
        1.0
    }
}

// This function calculates the correction factor rhohind(E)/rhofree(E)
// for the densities of states of an ensemble of S oscillators and one
// hindered rotor with a hindrance potential V0.
//
// S: number of oscillators
// V0: hindrance potential
// E: energy
//
// Returns the correction factor RRHIND.
pub fn rrhind(s: i32, v0: f64, e: f64) -> f64 {
    const PI: f64 = std::f64::consts::PI;  // Using the value of pi
    let mut ic: i32 = 0;
    let mut gs: f64;
    let mut a: f64;
    let mut b = vec![0.0; s];  // Array to store B coefficients

    if v0 > 0.0 && e > 0.0 {
        ic += 1;
        if ic == 1 {
            gs = gamma_func(s as f64 + 0.5);
            a = gs / gamma_func(PI * (s as f64 + 1.0));

            // Calculate B coefficients
            for ny in 0..s - 1 {
                b[ny as usize] = (-1.0f64).powi(ny)
                    * gs
                    * (ny as f64 + 2.5)
                    / (gamma_func(ny as f64 + 1.0) * gamma_func(s as f64 - 1.0 - ny as f64) * PI * (ny as f64 + 2.0) * (ny as f64 + 1.5));
            }
        }

        let ve = v0 / e;
        if ve < 1.0 {
            let mut rhind = 1.0;
            for i in 0..s - 1 {
                rhind -= b[i as usize] * ve.powf(i as f64 + 1.5);
            }
            rhind
        } else {
            let rhind = a / ve.sqrt();
            rhind
        }
    } else {
        1.0
    }
}






