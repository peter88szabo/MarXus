/// Solve A x = b where A is symmetric tridiagonal SPD.
///
/// A is represented by:
/// - diag[i] = A_{i,i}, length n
/// - offdiag[i] = A_{i+1,i} = A_{i,i+1}, length n-1
pub(crate) fn solve_spd_symmetric_tridiagonal(
    diag: &[f64],
    offdiag: &[f64],
    rhs: &[f64],
) -> Result<Vec<f64>, String> {
    let n = diag.len();
    if rhs.len() != n {
        return Err("RHS length mismatch in tridiagonal solve.".into());
    }
    if n == 0 {
        return Err("Empty system.".into());
    }
    if offdiag.len() + 1 != n && n > 1 {
        return Err("offdiag length mismatch in tridiagonal solve.".into());
    }

    let mut l_diag = vec![0.0; n];
    let mut l_sub = vec![0.0; n.saturating_sub(1)];

    l_diag[0] = diag[0].sqrt();
    if !l_diag[0].is_finite() || l_diag[0] <= 0.0 {
        return Err("Cholesky failed at i=0 (non-SPD).".into());
    }

    for i in 0..(n - 1) {
        let m = offdiag[i] / l_diag[i];
        l_sub[i] = m;

        let next = diag[i + 1] - m * m;
        if next <= 0.0 || !next.is_finite() {
            return Err(format!("Cholesky failed at i={} (non-SPD).", i + 1));
        }
        l_diag[i + 1] = next.sqrt();
    }

    let mut y = vec![0.0; n];
    y[0] = rhs[0] / l_diag[0];
    for i in 1..n {
        y[i] = (rhs[i] - l_sub[i - 1] * y[i - 1]) / l_diag[i];
    }

    let mut x = vec![0.0; n];
    x[n - 1] = y[n - 1] / l_diag[n - 1];
    for i_rev in 0..(n - 1) {
        let i = (n - 2) - i_rev;
        x[i] = (y[i] - l_sub[i] * x[i + 1]) / l_diag[i];
    }

    Ok(x)
}

/// Solve A x = b for a (not-necessarily-SPD) tridiagonal system via LU/Thomas.
///
/// A is represented by:
/// - subdiag[i] = A_{i+1,i}, length n-1
/// - diag[i] = A_{i,i}, length n
/// - superdiag[i] = A_{i,i+1}, length n-1
pub(crate) fn solve_tridiagonal_general(
    subdiag: &[f64],
    diag: &[f64],
    superdiag: &[f64],
    rhs: &[f64],
) -> Result<Vec<f64>, String> {
    let n = diag.len();
    if rhs.len() != n {
        return Err("RHS length mismatch in tridiagonal solve.".into());
    }
    if n == 0 {
        return Err("Empty system.".into());
    }
    if n > 1 && (subdiag.len() + 1 != n || superdiag.len() + 1 != n) {
        return Err("subdiag/superdiag length mismatch in tridiagonal solve.".into());
    }

    let mut c_prime = vec![0.0; n.saturating_sub(1)];
    let mut d_prime = vec![0.0; n];

    let mut denom = diag[0];
    if !denom.is_finite() || denom.abs() < 1e-30 {
        return Err("Tridiagonal LU failed at i=0 (zero/invalid pivot).".into());
    }

    if n > 1 {
        c_prime[0] = superdiag[0] / denom;
    }
    d_prime[0] = rhs[0] / denom;

    for i in 1..n {
        denom = diag[i] - subdiag[i - 1] * c_prime.get(i - 1).copied().unwrap_or(0.0);
        if !denom.is_finite() || denom.abs() < 1e-30 {
            return Err(format!(
                "Tridiagonal LU failed at i={} (zero/invalid pivot).",
                i
            ));
        }
        if i < n - 1 {
            c_prime[i] = superdiag[i] / denom;
        }
        d_prime[i] = (rhs[i] - subdiag[i - 1] * d_prime[i - 1]) / denom;
    }

    let mut x = vec![0.0; n];
    x[n - 1] = d_prime[n - 1];
    for i_rev in 0..(n - 1) {
        let i = (n - 2) - i_rev;
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    Ok(x)
}

