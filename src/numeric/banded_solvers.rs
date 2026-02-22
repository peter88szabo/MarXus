use super::linear_algebra::{cholesky_solve_spd, DenseMatrix};

fn banded_matvec_general(
    diag: &[f64],
    upper: &[Vec<f64>],
    lower: &[Vec<f64>],
    x: &[f64],
) -> Vec<f64> {
    let n = diag.len();
    let mut y = vec![0.0; n];
    for i in 0..n {
        y[i] = diag[i] * x[i];
    }
    let bw = upper.len().min(lower.len());
    for d in 1..=bw {
        for i in 0..n.saturating_sub(d) {
            let j = i + d;
            y[i] += upper[d - 1][i] * x[j];
            y[j] += lower[d - 1][i] * x[i];
        }
    }
    y
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn l2_norm(v: &[f64]) -> f64 {
    v.iter().copied().map(|x| x * x).sum::<f64>().sqrt()
}

pub(crate) fn solve_banded_bicgstab_left_jacobi(
    diag: &[f64],
    upper: &[Vec<f64>],
    lower: &[Vec<f64>],
    rhs: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<Vec<f64>, String> {
    let n = diag.len();
    if rhs.len() != n {
        return Err("RHS length mismatch in BiCGSTAB solve.".into());
    }
    if n == 0 {
        return Err("Empty system in BiCGSTAB solve.".into());
    }
    if tol <= 0.0 || !tol.is_finite() {
        return Err("Invalid tolerance in BiCGSTAB solve.".into());
    }
    if diag.iter().any(|d| !d.is_finite() || *d == 0.0) {
        return Err("Invalid diagonal in BiCGSTAB solve.".into());
    }

    // Left Jacobi preconditioner: solve A' x = b', A' = M^{-1}A, b' = M^{-1}b with M=diag(A).
    let mut b_prime = vec![0.0; n];
    for i in 0..n {
        b_prime[i] = rhs[i] / diag[i];
    }
    let b_norm = l2_norm(&b_prime).max(1e-300);

    let mut x = vec![0.0; n];
    let mut r = b_prime.clone(); // x=0 => r=b'
    let r0 = r.clone();

    let mut p = vec![0.0; n];
    let mut v = vec![0.0; n];
    let mut s = vec![0.0; n];
    let mut t = vec![0.0; n];

    let mut rho_old = 1.0_f64;
    let mut alpha = 1.0_f64;
    let mut omega = 1.0_f64;

    if l2_norm(&r) / b_norm <= tol {
        return Ok(x);
    }

    for _iter in 0..max_iter {
        let rho_new = dot(&r0, &r);
        if rho_new.abs() < 1e-300 {
            return Err("BiCGSTAB breakdown (rho ~ 0).".into());
        }

        let beta = (rho_new / rho_old) * (alpha / omega);
        for i in 0..n {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }

        // v = A' p = M^{-1}(A p)
        let ap = banded_matvec_general(diag, upper, lower, &p);
        for i in 0..n {
            v[i] = ap[i] / diag[i];
        }

        let denom = dot(&r0, &v);
        if denom.abs() < 1e-300 {
            return Err("BiCGSTAB breakdown (r0·v ~ 0).".into());
        }
        alpha = rho_new / denom;

        for i in 0..n {
            s[i] = r[i] - alpha * v[i];
        }
        if l2_norm(&s) / b_norm <= tol {
            for i in 0..n {
                x[i] += alpha * p[i];
            }
            return Ok(x);
        }

        // t = A' s
        let as_vec = banded_matvec_general(diag, upper, lower, &s);
        for i in 0..n {
            t[i] = as_vec[i] / diag[i];
        }

        let t_dot_t = dot(&t, &t);
        if t_dot_t.abs() < 1e-300 {
            return Err("BiCGSTAB breakdown (t·t ~ 0).".into());
        }
        omega = dot(&t, &s) / t_dot_t;
        if omega.abs() < 1e-300 {
            return Err("BiCGSTAB breakdown (omega ~ 0).".into());
        }

        for i in 0..n {
            x[i] += alpha * p[i] + omega * s[i];
        }
        for i in 0..n {
            r[i] = s[i] - omega * t[i];
        }

        if l2_norm(&r) / b_norm <= tol {
            return Ok(x);
        }

        rho_old = rho_new;
    }

    Err("BiCGSTAB did not converge within max_iter.".into())
}

pub(crate) fn solve_banded_gauss_seidel(
    diag: &[f64],
    upper: &[Vec<f64>],
    lower: &[Vec<f64>],
    rhs: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<Vec<f64>, String> {
    let n = diag.len();
    if rhs.len() != n {
        return Err("RHS length mismatch in Gauss-Seidel solve.".into());
    }
    if n == 0 {
        return Err("Empty system in Gauss-Seidel solve.".into());
    }
    if diag.iter().any(|d| !d.is_finite() || *d == 0.0) {
        return Err("Invalid diagonal in Gauss-Seidel solve.".into());
    }

    let bw = upper.len().min(lower.len());
    let rhs_norm = l2_norm(rhs).max(1e-300);
    let mut x = vec![0.0; n];

    for _iter in 0..max_iter {
        let x_old = x.clone();
        for i in 0..n {
            let mut sum = rhs[i];

            // columns < i: (row i, col i-d) stored in lower[d-1][i-d]
            let dmax = bw.min(i);
            for d in 1..=dmax {
                let col = i - d;
                sum -= lower[d - 1][col] * x[col];
            }

            // columns > i: (row i, col i+d) stored in upper[d-1][i]
            let dmax = bw.min(n - 1 - i);
            for d in 1..=dmax {
                let col = i + d;
                sum -= upper[d - 1][i] * x_old[col];
            }

            x[i] = sum / diag[i];
        }

        let ax = banded_matvec_general(diag, upper, lower, &x);
        let mut r = vec![0.0; n];
        for i in 0..n {
            r[i] = ax[i] - rhs[i];
        }
        if l2_norm(&r) <= tol * rhs_norm {
            return Ok(x);
        }
    }

    Err("Gauss-Seidel did not converge within max_iter.".into())
}

pub(crate) fn solve_spd_symmetric_banded_cholesky(
    diag: &[f64],
    bands: &[Vec<f64>],
    rhs: &[f64],
) -> Result<Vec<f64>, String> {
    let n = diag.len();
    if rhs.len() != n {
        return Err("RHS length mismatch in banded Cholesky solve.".into());
    }
    if n == 0 {
        return Err("Empty system in banded Cholesky solve.".into());
    }
    if diag.iter().any(|d| !d.is_finite() || *d <= 0.0) {
        return Err("Invalid diagonal in banded Cholesky solve.".into());
    }
    if bands.iter().any(|b| b.iter().any(|x| !x.is_finite())) {
        return Err("Invalid band entries in banded Cholesky solve.".into());
    }

    // A is represented by:
    // - diag[i] = A_{i,i}
    // - bands[d-1][i] = A_{i, i+d} = A_{i+d, i} for d=1..bw
    let bw = bands.len();
    for (d, b) in bands.iter().enumerate() {
        let expected = n.saturating_sub(d + 1);
        if b.len() != expected {
            return Err("Band length mismatch in banded Cholesky solve.".into());
        }
    }

    // Cholesky factor L (lower-triangular banded), stored similarly:
    // - l_diag[i] = L_{i,i}
    // - l_sub[d-1][i] = L_{i+d, i}
    let mut l_diag = vec![0.0; n];
    let mut l_sub = (0..bw)
        .map(|d| vec![0.0; n.saturating_sub(d + 1)])
        .collect::<Vec<Vec<f64>>>();

    for i in 0..n {
        let k_min = i.saturating_sub(bw);

        // Compute L_{i,i}
        let mut sum = diag[i];
        for k in k_min..i {
            let d = i - k;
            if d == 0 || d > bw {
                continue;
            }
            let lik = l_sub[d - 1][k];
            sum -= lik * lik;
        }
        if !(sum > 0.0) || !sum.is_finite() {
            return Err(format!(
                "Banded Cholesky failed at i={} (non-SPD pivot).",
                i
            ));
        }
        l_diag[i] = sum.sqrt();

        // Compute column i below diagonal within bandwidth
        let j_max = (i + bw).min(n.saturating_sub(1));
        for j in (i + 1)..=j_max {
            let dj = j - i;
            let mut a_ji = bands[dj - 1][i];

            // subtract Σ_k L_{j,k} L_{i,k}, with k in intersection of bands
            let k_start = k_min.max(j.saturating_sub(bw));
            for k in k_start..i {
                let dik = i - k;
                let djk = j - k;
                if dik == 0 || djk == 0 || dik > bw || djk > bw {
                    continue;
                }
                let lik = l_sub[dik - 1][k];
                let ljk = l_sub[djk - 1][k];
                a_ji -= ljk * lik;
            }

            l_sub[dj - 1][i] = a_ji / l_diag[i];
        }
    }

    // Forward solve: L y = rhs
    let mut y = vec![0.0; n];
    for i in 0..n {
        let k_min = i.saturating_sub(bw);
        let mut sum = rhs[i];
        for k in k_min..i {
            let d = i - k;
            if d == 0 || d > bw {
                continue;
            }
            sum -= l_sub[d - 1][k] * y[k];
        }
        y[i] = sum / l_diag[i];
    }

    // Back solve: L^T x = y
    let mut x = vec![0.0; n];
    for i_rev in 0..n {
        let i = n - 1 - i_rev;
        let j_max = (i + bw).min(n.saturating_sub(1));
        let mut sum = y[i];
        for j in (i + 1)..=j_max {
            let d = j - i;
            sum -= l_sub[d - 1][i] * x[j];
        }
        x[i] = sum / l_diag[i];
    }

    Ok(x)
}

pub(crate) fn solve_spd_symmetric_banded_dense_fallback(
    diag: &[f64],
    bands: &[Vec<f64>],
    rhs: &[f64],
) -> Result<Vec<f64>, String> {
    let n = diag.len();
    if rhs.len() != n {
        return Err("RHS length mismatch in dense fallback solve.".into());
    }
    let mut a = DenseMatrix::zeros(n);
    for i in 0..n {
        a.set(i, i, diag[i]);
    }
    for d in 1..=bands.len() {
        let b = &bands[d - 1];
        for i in 0..b.len() {
            let j = i + d;
            let v = b[i];
            a.set(i, j, v);
            a.set(j, i, v);
        }
    }
    cholesky_solve_spd(&a, rhs).map_err(|e| format!("Dense Cholesky fallback failed: {e}"))
}
