//ALL Wl(E,J) values are meant for rotating contribution of collision partners


/// This subroutine calculates the overall sum of states by convolution
/// of Wl(E,J) with the vibrational degrees of freedom of the fragments
/// via the Beyer-Swinehart algorithm followed by the convolution with
/// eventually existing free internal rotors of the fragments.
pub fn wEJ_vib_rot_pst(
    kzp: usize,
    nce: usize,
    emax: f64,
    length: f64,
    j: usize,
    b1: f64,
    b2: f64,
    nvib1: usize,
    nvib2: usize,
    ny1: &[f64],
    ny2: &[f64],
    nroti1: usize,
    nroti2: usize,
    mw: &mut [f64],
    mt: &mut [f64],
    rhoir: &mut [f64],
    air1: &[f64],
    air2: &[f64],
) {
    match kzp {
        1 => wlla(emax, length, j, b1, mt),
        2 => wlsa(emax, length, j, b1, mt),
        3 => W0_linear_and_linear(nce, emax, length, j, b1, b2, mt),
        4 => W0_linear_and_sphericaltop(nce, emax, length, j, b1, b2, mt),
        5 => W0_sphericaltop_and_sphericaltop(nce, emax, length, j, b1, b2, mt),
        _ => panic!("invalid KZP value: {}", kzp),
    }

    let nvibg = nvib1 + nvib2;

    let mut ny = vec![0.0_f64; nvibg];
    ny[..nvib1].copy_from_slice(&ny1[..nvib1]);
    ny[nvib1..nvibg].copy_from_slice(&ny2[..nvib2]);

    let mut mr = vec![0usize; nvibg];
    for i in 0..nvibg {
        mr[i] = (ny[i] / length + 0.5).floor() as usize;
    }

    let n = (emax / length + 0.5).floor() as usize;

    count(nvibg, n, &mr, mt);

    for i in 1..=n {
        mt[i] += mt[i - 1];
    }

    if nroti1 == 0 && nroti2 == 0 {
        mw[..(n + 1)].copy_from_slice(&mt[..(n + 1)]);
        return;
    }

    rhoiro(n, length, nroti1, nroti2, air1, air2, rhoir);

    mw[0] = 1.0;
    for i in 1..=n {
        mw[i] = 0.0;
        for l in 0..=i {
            mw[i] += mt[l] * rhoir[i - l];
        }
        mw[i] *= length;
    }
}

//--------------------------------------------------------------------------------
pub fn W0_atom_and_linear(aemax: f64, length: f64, j: usize, b: f64, mt: &mut [f64]) {

    fn wlliat(e: f64, j: usize, b: f64) -> f64 {
        let jmax = ((0.25 + e / b).sqrt() - 0.5).floor() as usize;

        if j > jmax {
            ((jmax + 1) * (jmax + 1)) as f64
        } else {
            ((2 * j + 1) * (jmax + 1) - j * (j + 1)) as f64
        }
    }

    let n = (aemax / length + 0.5).floor() as usize;

    mt[0] = 1.0;

    for i in 1..=n {
        let e = (i as f64) * length;
        let em = e - length;
        mt[i] = wlliat(e, j, b) - wlliat(em, j, b);
    }
}
//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------
pub fn W0_atom_and_sphericaltop(aemax: f64, length: f64, j: usize, b: f64, mt: &mut [f64]) {

    /// This function calculates exact values of Wl(E,J) for systems
    /// spherical top/atom.
    fn wlstat(e: f64, j: usize, b: f64) -> f64 {
        let jmax = ((0.25 + e / b).sqrt() - 0.5).floor() as usize;

        if j > jmax {
            ((jmax + 1) * (2 * jmax + 3) * (2 * jmax + 1)) as f64 / 3.0
        } else {
            ((2 * j + 1) * (jmax + 1) * (jmax + 1)) as f64
                - ((2 * j + 1) * (j + 1) * j) as f64 / 3.0
        }
    }

    let n = (aemax / length + 0.5).floor() as usize;

    mt[0] = 1.0;

    for i in 1..=n {
        let e = (i as f64) * length;
        let em = e - length;
        mt[i] = wlstat(e, j, b) - wlstat(em, j, b);
    }
}
//--------------------------------------------------------------------------------


/// This subroutine estimates Wl(E,J) for systems linear/linear and
/// creates the starting vector for the Beyer-Swinehart count.
pub fn W0_linear_and_linear(
    nce: i32,
    aemax: f64,
    length: f64,
    jg: usize,
    b1: f64,
    b2: f64,
    mt: &mut [f64],
) {
    let n = (aemax / length + 0.5).floor() as usize;

    mt[0] = 1.0;

    if nce <= 0 {
        let mut nsch = 0;

        for i in 1..=n {
            let e = (i as f64) * length;
            let em = e - length;

            if nsch == 0 {
                let exact = wllili(e, jg, b1, b2);

                let orm = wlhill(e, b1, b2);
                let arg = ((2 * jg + 1) as f64) * wlloll(e, b1, b2) / orm;
                let class = orm * polat(arg);

                let abw = ((exact - class).abs() / exact) * 100.0;

                if abw <= 2.0 {
                    nsch = 1;
                }

                mt[i] = exact - wllili(em, jg, b1, b2);
            } else {
                let orm = wlhill(e, b1, b2);
                let arg = ((2 * jg + 1) as f64) * wlloll(e, b1, b2) / orm;
                let class = orm * polat(arg);

                let orm = wlhill(em, b1, b2);
                let arg = ((2 * jg + 1) as f64) * wlloll(em, b1, b2) / orm;
                let classn = orm * polat(arg);

                mt[i] = class - classn;
            }
        }
    } else {
        for i in 1..=n {
            let e = (i as f64) * length;
            let em = e - length;
            mt[i] = wllili(e, jg, b1, b2) - wllili(em, jg, b1, b2);
        }
    }
}



pub fn W0_linear_and_sphericaltop(
    nce: i32,
    aemax: f64,
    length: f64,
    jg: usize,
    bst: f64,
    bli: f64,
    mt: &mut [f64],
) {
    let n = (aemax / length + 0.5).floor() as usize;

    mt[0] = 1.0;

    if nce <= 0 {
        let mut nsch = 0;

        for i in 1..=n {
            let e = (i as f64) * length;
            let em = e - length;

            if nsch == 0 {
                let exact = wlstli(e, jg, bst, bli);

                let orm = wlhisl(e, bst, bli);
                let arg = ((2 * jg + 1) as f64) * wllosl(e, bst, bli) / orm;
                let class = orm * polat(arg);

                let abw = ((exact - class).abs() / exact) * 100.0;

                if abw <= 2.0 {
                    nsch = 1;
                }

                mt[i] = exact - wlstli(em, jg, bst, bli);
            } else {
                let orm = wlhisl(e, bst, bli);
                let arg = ((2 * jg + 1) as f64) * wllosl(e, bst, bli) / orm;
                let class = orm * polat(arg);

                let orm = wlhisl(em, bst, bli);
                let arg = ((2 * jg + 1) as f64) * wllosl(em, bst, bli) / orm;
                let classn = orm * polat(arg);

                mt[i] = class - classn;
            }
        }
    } else {
        for i in 1..=n {
            let e = (i as f64) * length;
            let em = e - length;
            mt[i] = wlstli(e, jg, bst, bli) - wlstli(em, jg, bst, bli);
        }
    }
}




pub fn W0_sphericaltop_and_sphericaltop(
    nce: i32,
    aemax: f64,
    length: f64,
    jg: usize,
    b1: f64,
    b2: f64,
    mt: &mut [f64],
) {
    let n = (aemax / length + 0.5).floor() as usize;

    mt[0] = 1.0;

    if nce <= 0 {
        let mut nsch = 0;

        for i in 1..=n {
            let e = (i as f64) * length;
            let em = e - length;

            if nsch == 0 {
                let exact = wlstst(e, jg, b1, b2);

                let orm = wlhiss(e, b1, b2);
                let arg = ((2 * jg + 1) as f64) * wlloss(e, b1, b2) / orm;
                let class = orm * polat(arg);

                let abw = ((exact - class).abs() / exact) * 100.0;

                if abw <= 2.0 {
                    nsch = 1;
                }

                mt[i] = exact - wlstst(em, jg, b1, b2);
            } else {
                let orm = wlhiss(e, b1, b2);
                let arg = ((2 * jg + 1) as f64) * wlloss(e, b1, b2) / orm;
                let class = orm * polat(arg);

                let orm = wlhiss(em, b1, b2);
                let arg = ((2 * jg + 1) as f64) * wlloss(em, b1, b2) / orm;
                let classn = orm * polat(arg);

                mt[i] = class - classn;
            }
        }
    } else {
        for i in 1..=n {
            let e = (i as f64) * length;
            let em = e - length;
            mt[i] = wlstst(e, jg, b1, b2) - wlstst(em, jg, b1, b2);
        }
    }
}



/// This function counts exact values of Wl(E,J) for systems linear/
/// linear.
pub fn wllili(e: f64, jg: usize, b1: f64, b2: f64) -> f64 {
    let mut wllili = 0.0;

    let j1max = ((0.25 + e / b1).sqrt() - 0.5).floor() as usize;

    for j1 in 0..=j1max {
        let e_rem = e - b1 * (j1 * (j1 + 1)) as f64;

        if e_rem < 0.0 {
            continue;
        }

        let j2max = ((0.25 + e_rem / b2).sqrt() - 0.5).floor() as usize;

        for j2 in 0..=j2max {
            let j_min = j1.abs_diff(j2);
            let j_max = j1 + j2;

            for j in j_min..=j_max {
                let l_min = jg.abs_diff(j);
                let l_max = jg + j;

                for _l in l_min..=l_max {
                    wllili += 1.0;
                }
            }
        }
    }

    wllili
}

/// This function calculates classically Wl(E,J=0) for systems linear/
/// linear.
pub fn wlloll(e: f64, b1: f64, b2: f64) -> f64 {
    2.0 * e * e.sqrt() * (b1.sqrt() + b2.sqrt() - (b1 + b2).sqrt()) / (3.0 * b1 * b2)
}

/// This function calculates classically Wl(E,highJ) for systems linear/
/// linear.
pub fn wlhill(e: f64, b1: f64, b2: f64) -> f64 {
    e * e / (2.0 * b1 * b2)
}

/// This function counts exact values of Wl(E,J) for systems spherical
/// top/linear.
pub fn wlstli(e: f64, jg: usize, bst: f64, bli: f64) -> f64 {
    let mut wlstli = 0.0;

    let j1max = ((0.25 + e / bst).sqrt() - 0.5).floor() as usize;

    for j1 in 0..=j1max {
        let e_rem = e - bst * (j1 * (j1 + 1)) as f64;

        if e_rem < 0.0 {
            continue;
        }

        let j2max = ((0.25 + e_rem / bli).sqrt() - 0.5).floor() as usize;
        let swl = (2 * j1 + 1) as f64;

        for j2 in 0..=j2max {
            let j_min = j1.abs_diff(j2);
            let j_max = j1 + j2;

            for j in j_min..=j_max {
                let l_min = jg.abs_diff(j);
                let l_max = jg + j;

                for _l in l_min..=l_max {
                    wlstli += swl;
                }
            }
        }
    }

    wlstli
}


/// This function calculates classically Wl(E,J=0) for systems
/// spherical top/linear.
pub fn wllosl(e: f64, bst: f64, bli: f64) -> f64 {
    e * e * (bst / (bst + bli)).sqrt().asin()
        / (2.0 * bst * bst.sqrt() * bli.sqrt())
}

/// This function calculates classically Wl(E,highJ) for systems
/// spherical top/linear.
pub fn wlhisl(e: f64, bst: f64, bli: f64) -> f64 {
    8.0 * e * e * e.sqrt() / (15.0 * bst * bst.sqrt() * bli)
}

/// This function counts exact values of Wl(E,J) for systems spherical
/// top/spherical top.
pub fn wlstst(e: f64, jg: usize, bst1: f64, bst2: f64) -> f64 {
    let mut wlstst = 0.0;

    let j1max = ((0.25 + e / bst1).sqrt() - 0.5).floor() as usize;

    for j1 in 0..=j1max {
        let e_rem = e - bst1 * (j1 * (j1 + 1)) as f64;

        if e_rem < 0.0 {
            continue;
        }

        let j2max = ((0.25 + e_rem / bst2).sqrt() - 0.5).floor() as usize;
        let swl1 = (2 * j1 + 1) as f64;

        for j2 in 0..=j2max {
            let swl = (2 * j2 + 1) as f64 * swl1;

            let j_min = j1.abs_diff(j2);
            let j_max = j1 + j2;

            for j in j_min..=j_max {
                let l_min = jg.abs_diff(j);
                let l_max = jg + j;

                for _l in l_min..=l_max {
                    wlstst += swl;
                }
            }
        }
    }

    wlstst
}

/// This function calculates classically Wl(E,J=0) for systems spherical
/// top/spherical top.
pub fn wlloss(e: f64, bst1: f64, bst2: f64) -> f64 {
    8.0 * e * e * e.sqrt() / (15.0 * bst1 * bst2 * (bst1 + bst2).sqrt())
}

/// This function calculates classically Wl(E,highJ) for systems
/// spherical top/spherical top.
pub fn wlhiss(e: f64, bst1: f64, bst2: f64) -> f64 {
    3.1416 * e * e * e / (6.0 * bst1 * bst2 * (bst1 * bst2).sqrt())
}




/// This function interpolates Wl(E,J) between Wl(E,J=0) and Wl(E,highJ)
/// in a doubly reduced representation.
pub fn polat(x: f64) -> f64 {
    x.tanh() * (0.08 * x.sqrt() * (-1.1 * (x - 1.0) * (x - 1.0)).exp() + 1.0)
}

/// This function calculates the classical density of states of the
/// internal rotations of the fragments.
///
pub fn rhoiro(
    n: usize,
    length: f64,
    nroti1: usize,
    nroti2: usize,
    air1: &[f64],
    air2: &[f64],
    rhos: &mut [f64],
) {
    let r = (nroti1 + nroti2) as f64;

    rhos[0] = 0.0;

    let mut sqrb = 1.0;

    for i in 0..nroti1 {
        sqrb *= air1[i].sqrt();
    }

    for i in 0..nroti2 {
        sqrb *= air2[i].sqrt();
    }

    let gd = r / 2.0 - 1.0;

    let gamma = gamma_func(r / 2.0);

    for i in 1..=n {
        let e = (i as f64) * length;

        rhos[i] = 1.7725_f64.powf(r) / (gamma * sqrb) * e.powf(gd);
    }
}
use crate::numeric::lanczos_gamma::gamma_func;
