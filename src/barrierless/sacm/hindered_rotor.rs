use crate::numeric::lanczos_gamma::gamma_func;

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
pub fn hindered_rotor_sum_factor(s: i32, v0: f64, e: f64) -> f64 {
    if v0 <= 0.0 || e <= 0.0 || s <= 0 {
        return 1.0;
    }

    let s_usize = s as usize;
    let sqrt_pi = std::f64::consts::PI.sqrt();
    let gs = gamma_func(s as f64 + 1.5);
    let a = gs / (sqrt_pi * gamma_func(s as f64 + 2.0));

    let mut b = vec![0.0; s_usize];
    for ny in 0..s_usize {
        let ny_f = ny as f64;
        b[ny] = (-1.0_f64).powi(ny as i32)
            * gs
            * (ny_f + 2.5)
            / (gamma_func(ny_f + 1.0)
                * gamma_func(s as f64 - ny_f)
                * sqrt_pi
                * (ny_f + 2.0)
                * (ny_f + 1.5));
    }

    let ve = v0 / e;
    if ve >= 1.0 {
        return a / ve.sqrt();
    }

    let mut whind = 1.0;
    for i in 0..s_usize {
        whind -= b[i] * ve.powf(i as f64 + 1.5);
    }
    whind
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
pub fn hindered_rotor_density_factor(s: i32, v0: f64, e: f64) -> f64 {
    if v0 <= 0.0 || e <= 0.0 || s <= 1 {
        return 1.0;
    }

    let s_usize = s as usize;
    let sqrt_pi = std::f64::consts::PI.sqrt();
    let gs = gamma_func(s as f64 + 0.5);
    let a = gs / (sqrt_pi * gamma_func(s as f64 + 1.0));

    let mut b = vec![0.0; s_usize - 1];
    for ny in 0..s_usize - 1 {
        let ny_f = ny as f64;
        b[ny] = (-1.0_f64).powi(ny as i32)
            * gs
            * (ny_f + 2.5)
            / (gamma_func(ny_f + 1.0)
                * gamma_func(s as f64 - 1.0 - ny_f)
                * sqrt_pi
                * (ny_f + 2.0)
                * (ny_f + 1.5));
    }

    let ve = v0 / e;
    if ve >= 1.0 {
        return a / ve.sqrt();
    }

    let mut rhind = 1.0;
    for i in 0..s_usize - 1 {
        rhind -= b[i] * ve.powf(i as f64 + 1.5);
    }
    rhind
}




