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

use std::f64::consts::PI;
use other_module::gamma; // Assuming the gamma function is defined in another module

fn whitten_rabinovitch_we_and_rhoe(
    nvib: usize,
    n_ebin: usize,
    de: f64,
    omega: &[f64],
    nrot: usize,
    brot: &[f64],
    rho_e: &mut [f64],
    we: &mut [f64],
) {
    let mut prod_b = 1.0;
    let mut evib_prod = 1.0;
    let rdim = nrot as f64 / 2.0;

    // Compute prodB = sqrt(Brot(1) * Brot(2) * ... * Brot(nrot))
    for &b in brot.iter() {
        prod_b *= b.sqrt();
    }

    // Compute Evib_prod = omega(1) * omega(2) * ... * omega(nvib)
    for &w in omega.iter() {
        evib_prod *= w;
    }

    let ezpe = omega.iter().sum::<f64>() / 2.0;
    let esum_sq = ezpe * ezpe * 4.0;
    let esq_sum = omega.iter().map(|&w| w * w).sum::<f64>();

    let gama_half = gamma(0.5).powi(nrot as i32);
    let const_rhoe = gamma(nvib as f64 + rdim);
    let const_we = gamma(nvib as f64 + 1.0 + rdim);
    let factor = (nvib as f64 - 1.0) * (nvib as f64 + rdim) / nvib as f64;
    let beta = factor * esq_sum / esum_sq;
    let dc = gama_half / (prod_b * evib_prod);

    we[0] = 1.0;
    rho_e[0] = 1.0;

    // Create W(E) and Rho(E) for nrot-dimensional classical rotor + nvib harmonic oscillator
    for i in 1..=n_ebin {
        let ene = (i as f64 - 0.5) * de;
        let ered = ene / ezpe;
        let const_a = 1.0 - beta * wfunc_wr(ered);
        let dum = ene + const_a * ezpe;

        we[i] = dc * dum.powf(nvib as f64 + rdim) / const_we;
        rho_e[i] = dc * dum.powf(nvib as f64 + rdim - 1.0) / const_rhoe;
    }
}

// Whitten-Rabinovitch W(Ered) function
fn wfunc_wr(ered: f64) -> f64 {
    if ered >= 1.0 {
        (-2.4191 * ered.powf(0.25)).exp()
    } else {
        1.0 / (5.0 * ered + 2.73 * ered.sqrt() + 3.51)
    }
}


//==========================================================================================
fn get_Jres_rovib_WEJ_or_rhoEJ -> Vec<f64>(
//==========================================================================================
    what: String,
    jtot: usize,
    de: f64,
    n_ebin: usize,
    nrot: usize,
    brot: &[f64],
    rho_e: &[f64],
    we: &[f64],
) {
    // ======================================================================================
    // Computation of the E,J-dependent Sum and Density of states
    // In the case of symmetric top (prolate or oblate), the K-rotors are averaged
    //
    // This routine requires the ro-vibrational W(E) and rho(E),
    // then the J-dependence is introduced as
    // W(E) --> W(E - Erot(J))   and   rho(E) --> rho(E - Erot(J))
    // ======================================================================================

    let mut rot_j_bin: f64;
    let mut min_ebin: usize;
    let mut ecorr: f64;
    let mut ecorr_a: f64;
    let mut a: f64;
    let mut b: f64;
    let mut kx: usize;
    let mut kmax: usize;

    // -----------------------------------------------------------------------
    // molecule is treated as spherical top
    // B = Brot(3)
    // rotJ_bin = int(B*J*(J+1)/dE+0.5)*dE
    // minEbin = int(rotJ_bin/dE)
    //
    // do i=minEbin,nEbin
    //     Ecorr = i*dE - rotJ_bin
    //     RhoEJ(i) = float(2*J+1)*RhoE(int(Ecorr/dE))
    // enddo
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // molecule is prolate [A >> B = C] or oblate [A = B >> C] symmetric top
    // The energy of the rotor E(J) = B*J*(J+1) + (A-B)*K2
    // B calculated as the geometric mean of the two small rotational constants
    //
    // Ref: W. Forst, Unimolecular Reactions (2003)
    //      Page 92-94, Chapter 4.5.9
    //
    // Although I used here Matthias Olzmanns's SACM routines as reference
    // not the formulas from the book of W. Forst
    // -----------------------------------------------------------------------


    // Prolate case: A >> B = C symmetric top
    a = brot[0];
    b = (brot[1] * brot[2]).sqrt();

    // Oblate:
    // A = Brot(3)
    // B = sqrt(Brot(1)*Brot(2))

    // Spherical top approximation
    //
    rot_j_bin = (b * (jtot * (jtot + 1)) as f64 / de + 0.5).floor() * de;

    min_ebin = (rot_j_bin / de).floor() as usize;

    for i in min_ebin..=n_ebin {
        ecorr = i as f64 * de - rot_j_bin;

        rho_ej[i] = rho_e[(ecorr / de).floor() as usize];
        wej[i] = we[(ecorr / de).floor() as usize];

        // Prolate: see Eq. 4.106 in reference
        kx = ((ecorr / (a - b)).sqrt().floor()) as usize;
        kmax = jtot.min(kx);

        for k in 1..=kmax {
            ecorr_a = ecorr - ((a - b) * (k * k) as f64 / de + 0.5).floor() * de;

            rho_ej[i] += 2.0 * rho_e[(ecorr_a / de).floor() as usize];
            wej[i] += 2.0 * we[(ecorr_a / de).floor() as usize];
        }
    }
}
//==========================================================================================


