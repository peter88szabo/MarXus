#![allow(non_snake_case)]

//=============================================================================================
// Beyer-Swinehart direct counting of rho(E) or W(E)
//=============================================================================================
// nvib      -->  number of vibrational modes
// nebin     -->  number of energy bins
// freq_bin  -->  which energy bin of the i'th oscillator
// res       -->  both input and output, the results of counting
// --------------------------------------------------------------------------------------------
// res can be either density or number of states depending on its starting (input) value
//
// when initial (input) value is
// res = [1, 0, 0, 0,...., 0] --> res = pure vibrational density of states
// res = [1, 1, 1, 1,...., 1] --> res = pure vibrational number of states
//
// or alternatively it can be initialized as pure rotational density or number of states
// to obtain the ro-vibrational energy-dependent rho(E) and W(E)
//=============================================================================================
pub fn beyer_swinehart_counting(nvib: usize, nebin: usize, freq_bin: &[usize], res: &Vec<f64>) -> Vec<f64> {
    let mut results = res.clone();

    for i in 0..nvib {
        let iosc = freq_bin[i];
        for j in iosc..=nebin {
            results[j] += results[j - iosc];
        }
    }
    return results;
}
//=============================================================================================



//=============================================================================================
// Calculating the ro-vibrational W(E) or rho(E)
//=============================================================================================
pub fn get_rovib_WE_or_rhoE(
    what: String,
    nvib: usize,
    nebin: usize,
    dE: f64,
    nrot: usize,
    freq_bin: &[usize],
    Brot: &[f64],
) -> Vec<f64> {
    // Initialize sum or density of states

    let mut res = vec![0.0; nebin + 1];

    if nrot == 0 {
        // only vibrations, no rotations

        match what.as_str() {
            "sum" => res = vec![1.0; nebin + 1],
            "den" => res[0] = 1.0,
            _ => println!("No Match Found in get_rovib_WE_or_rhoE()"),
        };
    } else {
        // rotations initialize the W(E) or rho(E) to be convoluted with vibrations later on

        res = get_pure_rotational_WE_or_rhoE(what, nebin, dE, nrot, &Brot);
    }

    // Convoluting the vibrations with rotations
    let result = beyer_swinehart_counting(nvib, nebin, &freq_bin, &res);

    return result;
}
//==========================================================================================

use crate::numeric::lanczos_gamma::gamma_func;
//==========================================================================================
pub fn get_pure_rotational_WE_or_rhoE(
    what: String,
    nebin: usize,
    dE: f64,
    nrot: usize,
    Brot: &[f64],
) -> Vec<f64> {
    // Product of rotational constants
    let mut prod_Brot = 1.0;
    for i in 0..nrot {
        prod_Brot *= Brot[i];
    }
    prod_Brot = prod_Brot.sqrt();

    let rdim = (nrot as f64) / 2.0;

    let sqrt_pi = f64::sqrt(std::f64::consts::PI);
    let crt = (f64::powf(sqrt_pi, nrot as f64)) / prod_Brot;

    let const_W = crt / gamma_func(1.0 + rdim);
    let const_rho = crt / gamma_func(rdim);

    let mut res = vec![0.0; nebin + 1];
    res[0] = 1.0;

    // Create W(E) or Rho(E) for nrot-dimensional classical rotor
    for i in 1..=nebin {
        let Ei = (i as f64) * dE;
        let Eicenter = ((i as f64) - 0.5) * dE;

        match what.as_str() {
            "sum" => res[i] = const_W * (f64::powf(Ei, rdim)),
            "den" => res[i] = const_rho * (f64::powf(Eicenter, rdim - 1.0)),
            _ => println!("Wrong mode used in get_pure_rotational_WE_or_rhoE()"),
        }
    }
    return res;
}
//=============================================================================================
