#![allow(non_snake_case)]

use crate::constants::PI;
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
pub fn beyer_swinehart_counting(
    nvib: usize,
    nebin: usize,
    freq_bin: &[usize],
    res: &Vec<f64>,
) -> Vec<f64> {
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
// Separable classical rigid free rotors, rotational constants in cm-1.
// Each rotor contributes one factor to Q'_r of Forst, Chem. Rev. 71, 339 (1971), eqs. 37-40:
//   OneD  (internal or K-rotor):  Q' = sqrt(pi/B) / sigma           [Baer & Hase eq. 6.25]
//   TwoD  (linear molecule):      Q' = 1 / (sigma B)                [Baer & Hase eq. 6.28]
//   Top3D (overall rotation):     Q' = sqrt(pi) / (sigma sqrt(ABC)) [Baer & Hase eqs. 6.37-6.39]
// Top3D is the exact classical asymmetric top, i.e. a 2D rotor with sqrt(BC) times a 1D rotor
// with A (Forst eq. 35), NOT three independent 1D rotors (Forst eq. 36), which is pi too large.
//
// Symmetry numbers: use sigma = 1 when the symmetry is already contained in the reaction path
// degeneracy, otherwise it is counted twice (Forst, Theory of Unimolecular Reactions, 1973, p. 91).
//=============================================================================================
#[derive(Debug, Clone, Copy)]
pub enum Rotor {
    OneD { b: f64, sigma: f64 },
    TwoD { b: f64, sigma: f64 },
    Top3D { a: f64, b: f64, c: f64, sigma: f64 },
}

impl Rotor {
    pub fn dimension(&self) -> usize {
        match self {
            Rotor::OneD { .. } => 1,
            Rotor::TwoD { .. } => 2,
            Rotor::Top3D { .. } => 3,
        }
    }

    pub fn q_prime(&self) -> f64 {
        match *self {
            Rotor::OneD { b, sigma } => (PI / b).sqrt() / sigma,
            Rotor::TwoD { b, sigma } => 1.0 / (sigma * b),
            Rotor::Top3D { a, b, c, sigma } => PI.sqrt() / (sigma * (a * b * c).sqrt()),
        }
    }
}

// The (nrot, Brot) convention used throughout MarXus: nrot = 0 atom, 1 one-dimensional rotor,
// 2 linear molecule (one 2D rotor), 3 nonlinear molecule (3D top). Symmetry numbers are handled
// by the callers (reaction path degeneracy), so sigma = 1 here.
pub fn rotors_from_brot(nrot: usize, Brot: &[f64]) -> Vec<Rotor> {
    match nrot {
        0 => vec![],
        1 => vec![Rotor::OneD { b: Brot[0], sigma: 1.0 }],
        2 => vec![Rotor::TwoD { b: (Brot[0] * Brot[1]).sqrt(), sigma: 1.0 }],
        3 => vec![Rotor::Top3D { a: Brot[0], b: Brot[1], c: Brot[2], sigma: 1.0 }],
        _ => panic!(
            "nrot = {nrot} does not identify a rotor set; build a Vec<Rotor> and call \
             get_rovib_rotors_WE_or_rhoE() instead"
        ),
    }
}

//==========================================================================================
// Rotational W(E) or rho(E) for any set of separable rotors of total dimension r
// (Forst 1971, eq. 43):   G_r(E) = Q'_r E^{r/2} / Gamma(1 + r/2)
//
// "sum": W[i]   = G_r(i dE)
// "den": rho[i] = (G_r(i dE) - G_r((i-1) dE)) / dE   states per cm-1, grain i = ((i-1)dE, i dE]
//        rho[0] = G_r(0) / dE                        (non-zero only for r = 0: one state at E=0)
// so that sum_{j<=i} rho[j] dE = W[i] holds exactly, also after convolution with vibrations.
//==========================================================================================
pub fn get_rotors_WE_or_rhoE(what: String, nebin: usize, dE: f64, rotors: &[Rotor]) -> Vec<f64> {
    let rdim: usize = rotors.iter().map(|r| r.dimension()).sum();
    let q_prime: f64 = rotors.iter().map(|r| r.q_prime()).product();

    let half_r = (rdim as f64) / 2.0;
    let const_G = q_prime / gamma_func(1.0 + half_r);
    let G = |e: f64| const_G * f64::powf(e, half_r);

    let mut res = vec![0.0; nebin + 1];

    match what.as_str() {
        "sum" => {
            for i in 0..=nebin {
                res[i] = G((i as f64) * dE);
            }
        }
        "den" => {
            res[0] = G(0.0) / dE;
            for i in 1..=nebin {
                res[i] = (G((i as f64) * dE) - G(((i - 1) as f64) * dE)) / dE;
            }
        }
        _ => panic!("Wrong mode '{what}' in get_rotors_WE_or_rhoE(), use \"sum\" or \"den\""),
    }
    return res;
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
    get_rotors_WE_or_rhoE(what, nebin, dE, &rotors_from_brot(nrot, Brot))
}
//=============================================================================================

//=============================================================================================
// Calculating the ro-vibrational W(E) or rho(E) (rho in states per cm-1)
//=============================================================================================
pub fn get_rovib_rotors_WE_or_rhoE(
    what: String,
    nvib: usize,
    nebin: usize,
    dE: f64,
    freq_bin: &[usize],
    rotors: &[Rotor],
) -> Vec<f64> {
    // The rotational W(E) or rho(E) is the starting array convoluted with the vibrations
    let res = get_rotors_WE_or_rhoE(what, nebin, dE, rotors);

    beyer_swinehart_counting(nvib, nebin, &freq_bin, &res)
}

pub fn get_rovib_WE_or_rhoE(
    what: String,
    nvib: usize,
    nebin: usize,
    dE: f64,
    nrot: usize,
    freq_bin: &[usize],
    Brot: &[f64],
) -> Vec<f64> {
    get_rovib_rotors_WE_or_rhoE(what, nvib, nebin, dE, freq_bin, &rotors_from_brot(nrot, Brot))
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
                        let kmin =
                            (((b * j * (j + 1.0) - ei) / delta).sqrt() + 1.0).floor() as usize;
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

#[cfg(test)]
mod tests {
    use super::*;

    fn rel_err(x: f64, reference: f64) -> f64 {
        ((x - reference) / reference).abs()
    }

    #[test]
    fn classical_rotor_sums_match_baer_hase() {
        // Baer & Hase, Unimolecular Reaction Dynamics (1996), B in cm-1:
        //   1D rotor        N(E) = 2 sqrt(E/B)                     eq. 6.25
        //   2D (linear)     N(E) = E/B                             eq. 6.28
        //   3D top          N(E) = sqrt(pi)/(3/2)! E^{3/2}/sqrt(ABC) = (4/3) E^{3/2}/sqrt(ABC)
        //                   (eq. 6.37; exact classical result for the asymmetric top)
        let d_e = 1.0;
        let nebin = 1000;
        let e = nebin as f64 * d_e;

        let w1 = get_pure_rotational_WE_or_rhoE("sum".to_string(), nebin, d_e, 1, &[5.0]);
        assert!(rel_err(w1[nebin], 2.0 * (e / 5.0).sqrt()) < 1e-12, "1D: {}", w1[nebin]);

        let w2 = get_pure_rotational_WE_or_rhoE("sum".to_string(), nebin, d_e, 2, &[0.5, 0.5]);
        assert!(rel_err(w2[nebin], e / 0.5) < 1e-12, "2D: {}", w2[nebin]);

        let (a, b, c) = (1.0, 0.5, 0.25);
        let w3 = get_pure_rotational_WE_or_rhoE("sum".to_string(), nebin, d_e, 3, &[a, b, c]);
        let w3_ref = 4.0 / 3.0 * e.powf(1.5) / (a * b * c).sqrt();
        assert!(rel_err(w3[nebin], w3_ref) < 1e-12, "3D: {} vs {}", w3[nebin], w3_ref);
    }

    #[test]
    fn classical_top_density_matches_baer_hase() {
        // Baer & Hase eq. 6.39: rho(E) = 2 sqrt(E) / sqrt(ABC)  (states per cm-1)
        let d_e = 1.0;
        let nebin = 2000;
        let (a, b, c) = (1.0, 0.5, 0.25);
        let rho = get_pure_rotational_WE_or_rhoE("den".to_string(), nebin, d_e, 3, &[a, b, c]);
        let e_mid = (nebin as f64 - 0.5) * d_e;
        let rho_ref = 2.0 * e_mid.sqrt() / (a * b * c).sqrt();
        assert!(rel_err(rho[nebin], rho_ref) < 1e-6, "{} vs {}", rho[nebin], rho_ref);
    }

    #[test]
    fn vibrational_density_is_per_wavenumber() {
        // One oscillator of 1000 cm-1 on a 10 cm-1 grain: one state in the grain at 1000 cm-1,
        // i.e. a density of 1/dE = 0.1 states per cm-1 there, independent of the rotor count.
        let d_e = 10.0;
        let nebin = 200;
        let rho = get_rovib_WE_or_rhoE("den".to_string(), 1, nebin, d_e, 0, &[100], &[]);
        assert!(rel_err(rho[100], 1.0 / d_e) < 1e-12, "rho(1000) = {}", rho[100]);
        assert_eq!(rho[50], 0.0);
    }

    #[test]
    fn general_rotor_set_follows_forst_eq43() {
        // Forst, Chem. Rev. 71, 339 (1971), eqs. 39-43: separable rotors of total dimension r,
        //   G_r(E) = Q'_r E^{r/2} / Gamma(1 + r/2),   Q'_r = product of the per-rotor factors.
        // 3D top (A,B,C; sigma=2) + one internal 1D rotor (B_int; sigma=3): r = 4,
        //   Q'_4 = [sqrt(pi)/(2 sqrt(ABC))] * [sqrt(pi/B_int)/3],  Gamma(3) = 2.
        let d_e = 1.0;
        let nebin = 3000;
        let e = nebin as f64 * d_e;
        let (a, b, c, b_int) = (1.0, 0.5, 0.25, 5.0);
        let rotors = [
            Rotor::Top3D { a, b, c, sigma: 2.0 },
            Rotor::OneD { b: b_int, sigma: 3.0 },
        ];
        let q_prime = (PI.sqrt() / (2.0 * (a * b * c).sqrt())) * ((PI / b_int).sqrt() / 3.0);

        let w = get_rotors_WE_or_rhoE("sum".to_string(), nebin, d_e, &rotors);
        let w_ref = q_prime * e * e / 2.0;
        assert!(rel_err(w[nebin], w_ref) < 1e-12, "{} vs {}", w[nebin], w_ref);

        // rho_r(E) = Q'_r E^{r/2-1} / Gamma(r/2) = Q'_4 E (grain average over the last grain)
        let rho = get_rotors_WE_or_rhoE("den".to_string(), nebin, d_e, &rotors);
        let rho_ref = q_prime * (e - 0.5 * d_e);
        assert!(rel_err(rho[nebin], rho_ref) < 1e-9, "{} vs {}", rho[nebin], rho_ref);
    }

    #[test]
    fn density_integrates_to_sum_of_states() {
        // With rho in states per cm-1, sum_{j<=i} rho[j]*dE must reproduce W[i] (rotors + vibrations).
        let d_e = 10.0;
        let nebin = 600;
        let freq_bin = [50, 120];
        let brot = [1.2, 0.4, 0.3];
        let rho = get_rovib_WE_or_rhoE("den".to_string(), 2, nebin, d_e, 3, &freq_bin, &brot);
        let w = get_rovib_WE_or_rhoE("sum".to_string(), 2, nebin, d_e, 3, &freq_bin, &brot);
        let mut cumulative = 0.0;
        for i in 0..=nebin {
            cumulative += rho[i] * d_e;
            assert!(
                (cumulative - w[i]).abs() <= 1e-9 * w[i].max(1.0),
                "i={i}: sum rho*dE = {cumulative}, W = {}",
                w[i]
            );
        }
    }
}
