use std::f64;

/// This function estimates by an iteration procedure the J-dependent
/// threshold energy E0(J) with a Morse-type radial part. If the return value is -1,
/// then the molecule is in a rotationally repulsive state for the J
/// under consideration.
///
/// dmorse: dissociation energy of the Morse potential
/// ezp: zero-point energy of the reaction products
/// dez: ezn - ezp - erc / 2
/// alpha_beta: alpha/beta
/// be: mean of the two smallest rotational constants of the reactant
/// a1, a2: centrifugal parameters (see ref.)
/// j: rotational quantum number
/// 
/// Returns the calculated E0(J) or -1 if the molecule is in a repulsive state.
/// 
/// [Ref.: J. Troe, J.Chem.Phys. 75, 226 (1981)]
pub fn morse_threshold_energy(
    dmorse: f64,
    ezp: f64,
    dez: f64,
    alpha_beta: f64,
    be: f64,
    a1: f64,
    a2: f64,
    j: i32
) -> f64 {
    let mut znew: f64 = 85.0;
    let b = dez / dmorse;

    if j == 0 {
        // Handle the J = 0 case separately
        handle_j_zero(b, alpha_beta, dmorse, ezp, &mut znew)
    } else {
        // Handle the general J > 0 case
        handle_j_nonzero(j, b, alpha_beta, dmorse, ezp, be, a1, a2, &mut znew)
    }
}

// Helper function for the case when J == 0
fn handle_j_zero(b: f64, alpha_beta: f64, dmorse: f64, ezp: f64, znew: &mut f64) -> f64 {
    let mut iter = 0;
    loop {
        let zold = *znew;
        *znew = (2.0 * (1.0 - (-zold).exp()) / (b * alpha_beta * (-alpha_beta * zold).exp())).ln();
        if !znew.is_finite() {
            return dmorse + ezp;
        }
        if *znew > zold {
            return dmorse + ezp;
        }

        let tol = (*znew - zold).abs();
        if tol <= 0.001 {
            break;
        }
        iter += 1;
        if iter > 1000 {
            return dmorse + ezp;
        }
    }

    // Default return if the loop finishes normally
    dmorse + ezp
}

// Helper function for the case when J > 0
fn handle_j_nonzero(
    j: i32,
    b: f64,
    alpha_beta: f64,
    dmorse: f64,
    ezp: f64,
    be: f64,
    a1: f64,
    a2: f64,
    znew: &mut f64 ) -> f64 {

    let a = be * (j as f64) * (j as f64 + 1.0) / dmorse;

    let mut iter = 0;
    loop {
        let zold = *znew;
        *znew = calculate_znew(a, b, alpha_beta, zold, a1, a2);

        if !znew.is_finite() || *znew < 0.0 {
            return -1.0; // Return -1 for the repulsive state
        }

        let tol = (*znew - zold).abs();
        if tol <= 0.001 {
            break;
        }
        iter += 1;
        if iter > 1000 {
            return -1.0;
        }
    }

    let enullj = calculate_enullj(*znew, be, j, dmorse, a1, a2, b, alpha_beta);
    enullj * dmorse + ezp
}

// Calculates the new ZNEW value during the iteration
fn calculate_znew(a: f64, b: f64, alpha_beta: f64, zold: f64, a1: f64, a2: f64) -> f64 {
    -(a * (a1 + 2.0 * a2 * zold) / (2.0 * (1.0 - (-zold).exp()) * (1.0 + a1 * zold + a2 * zold.powi(2)) * (1.0 + a1 * zold + a2 * zold.powi(2))) +
      b * alpha_beta * (-alpha_beta * zold).exp() / (2.0 * (1.0 - (-zold).exp()))).ln()
}

// Calculates the final value of ENULLJ
fn calculate_enullj(znew: f64, be: f64, j: i32, dmorse: f64, a1: f64, a2: f64, b: f64, alpha_beta: f64) -> f64 {
    let rotational_part = be * (j as f64) * (j as f64 + 1.0) / (dmorse * (1.0 + a1 * znew + a2 * znew.powi(2)));
    let exp_part = (1.0 - (-znew).exp()).powi(2);
    let additional_part = b * (-alpha_beta * znew).exp();

    exp_part + rotational_part + additional_part
}
