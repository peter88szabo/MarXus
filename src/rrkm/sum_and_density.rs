#![allow(non_snake_case)]

use crate::numeric::lanczos_gamma::gamma_func;

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
//

#[derive(Debug, Clone, Copy)]
pub enum RotorSymmetry {
    SphericalTop,
    OblateSymmetricTop,
    ProlateSymmetricTop,
}

#[derive(Debug, Clone)]
pub struct JResolvedStates {
    pub rho_ej: Vec<f64>,
    pub wej: Vec<f64>,
}

fn round_to_bin(energy: f64, dE: f64) -> f64 {
    ((energy / dE) + 0.5).floor() * dE
}

fn bin_index(energy: f64, dE: f64) -> usize {
    ((energy / dE) + 0.5).floor().max(0.0) as usize
}

//==========================================================================================
pub fn get_Jres_rovib_WEJ_or_rhoEJ(
//==========================================================================================
    rotor: RotorSymmetry,
    jtot: usize,
    dE: f64,
    n_ebin: usize,
    brot: &[f64],
    rho_e: &[f64],
    we: &[f64],
    b_effective: Option<f64>,
) -> JResolvedStates {
    // ======================================================================================
    // Computation of the E,J-dependent sum and density of states.
    // In the case of symmetric tops, the K-rotors are averaged.
    // This follows the SACM Fortran logic but uses clearer Rust names.
    // ======================================================================================
    assert!(rho_e.len() > n_ebin, "rho_e length must be n_ebin+1");
    assert!(we.len() > n_ebin, "we length must be n_ebin+1");

    let make_j_resolved = |base: &[f64]| -> Vec<f64> {
        let mut out = vec![0.0; n_ebin + 1];
        let j = jtot as f64;

        match rotor {
            RotorSymmetry::SphericalTop => {
                let b = b_effective.unwrap_or_else(|| brot.get(1).copied().unwrap_or(brot[0]));
                let rot_energy = round_to_bin(b * j * (j + 1.0), dE);
                let min_bin = bin_index(rot_energy, dE);
                for i in min_bin..=n_ebin {
                    let ecorr = i as f64 * dE - rot_energy;
                    let idx = bin_index(ecorr, dE);
                    out[i] = (2 * jtot + 1) as f64 * base[idx];
                }
            }
            RotorSymmetry::OblateSymmetricTop => {
                let b = b_effective.unwrap_or_else(|| brot.get(1).copied().unwrap_or(brot[0]));
                let c = brot.get(2).copied().unwrap_or(b);
                let delta = b - c;
                let rot_energy = round_to_bin(b * j * (j + 1.0), dE);
                let min_energy = b * j * (j + 1.0) - delta * j * j;
                let min_bin = bin_index(min_energy, dE);

                for i in 0..=n_ebin {
                    if jtot == 0 {
                        out[i] = base[i];
                        continue;
                    }
                    if i < min_bin {
                        continue;
                    }

                    let ei = i as f64 * dE;
                    let base_energy = ei - b * j * (j + 1.0) + delta * j * j;
                    let base_idx = bin_index(base_energy, dE);
                    out[i] = 2.0 * base[base_idx];

                    let (kmin, include_k0) = if ei < rot_energy {
                        let kmin = (((b * j * (j + 1.0) - ei) / delta).sqrt() + 1.0).floor() as usize;
                        (kmin, false)
                    } else {
                        (1, true)
                    };

                    if jtot > 1 && kmin <= jtot - 1 {
                        for k in (kmin..=jtot - 1).rev() {
                            let kf = k as f64;
                            let k_energy = ei - b * j * (j + 1.0) + delta * kf * kf;
                            let k_idx = bin_index(k_energy, dE);
                            out[i] += 2.0 * base[k_idx];
                        }
                    }

                    if include_k0 {
                        let k0_idx = bin_index(ei - rot_energy, dE);
                        out[i] += base[k0_idx];
                    }
                }
            }
            RotorSymmetry::ProlateSymmetricTop => {
                let b = b_effective.unwrap_or_else(|| brot.get(1).copied().unwrap_or(brot[0]));
                let a = brot.get(0).copied().unwrap_or(b);
                let delta = a - b;
                let rot_energy = round_to_bin(b * j * (j + 1.0), dE);
                let min_bin = bin_index(rot_energy, dE);

                for i in min_bin..=n_ebin {
                    let ecorr = i as f64 * dE - rot_energy;
                    let idx = bin_index(ecorr, dE);
                    out[i] = base[idx];

                    if delta > 0.0 {
                        let kx = (ecorr / delta).sqrt().floor() as usize;
                        let kmax = jtot.min(kx);
                        for k in 1..=kmax {
                            let kf = k as f64;
                            let k_energy = ecorr - round_to_bin(delta * kf * kf, dE);
                            let k_idx = bin_index(k_energy, dE);
                            out[i] += 2.0 * base[k_idx];
                        }
                    }
                }
            }
        }

        out
    };

    JResolvedStates {
        rho_ej: make_j_resolved(rho_e),
        wej: make_j_resolved(we),
    }
}
//==========================================================================================
