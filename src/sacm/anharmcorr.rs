// This function calculates anharmonicity factors for the densities of
// states of the reactant molecule.
// 
// kanh: parameter controlling the behavior of the function
// ene: energy
// nvib: number of vibrational states
// nyr: array of vibrational constants
// dr: array of anharmonic constants
// 
// Returns the anharmonicity factor for the reactant molecule.
pub fn anharmonic_reactant_rhoE(kanh: i32, ene: f64, nvib: usize, nyr: &[f64], dr: &[f64]) -> f64 {
    match kanh + 1 {
        3 => {
            let nenner = 2 * nvib - 3;
            let mut anh_rho = 1.0;
            for i in 0..nvibr {
                anh_rho *= 1.0 + (e + nyr[i] / 2.0) / (dr[i] * nenner as f64);
            }
            anh_rho
        }
        2 | _ => 1.0,  // Handle default or KANH+1 = 2 case
    }
}



// This function calculates anharmonicity factors for the sums of
// states of the product molecules.
//
// KANH: parameter controlling the behavior of the function
// EI: energy of the interaction
// NVIB1: number of vibrational states of the first product
// NVIB2: number of vibrational states of the second product
// NY1, NY2: arrays of vibrational constants for the two products
// D1, D2: arrays of anharmonic constants for the two products
//
// Returns the anharmonicity factor for the product molecules.
pub fn anharmonic_product_wE(
    kanh: i32,
    ei: f64,
    nvib1: usize,
    nvib2: usize,
    ny1: &[f64],
    ny2: &[f64],
    d1: &[f64],
    d2: &[f64]
) -> f64 {
    match kanh + 1 {
        3 => {
            let nenner = 2 * (nvib1 + nvib2) - 1;
            let mut anh_we = 1.0;
            for i in 0..nvib1 {
                anh_we *= 1.0 + (ei + ny1[i] / 2.0) / (d1[i] * nenner as f64);
            }
            for i in 0..nvib2 {
                anh_we *= 1.0 + (ei + ny2[i] / 2.0) / (d2[i] * nenner as f64);
            }
            anh_we
        }
        2 | _ => 1.0,  // Handle default or KANH+1 = 2 case
    }
}

