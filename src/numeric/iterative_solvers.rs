use super::linear_algebra::DenseMatrix;

pub struct BiCgStabDiagnostics {
    pub iterations: usize,
    pub final_residual_relative_l2: f64,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn l2_norm(v: &[f64]) -> f64 {
    v.iter().copied().map(|x| x * x).sum::<f64>().sqrt()
}

fn axpy_in_place(y: &mut [f64], a: f64, x: &[f64]) {
    for (yy, xx) in y.iter_mut().zip(x.iter()) {
        *yy += a * (*xx);
    }
}

fn sub_scaled(out: &mut [f64], a: &[f64], scale: f64, b: &[f64]) {
    for ((o, aa), bb) in out.iter_mut().zip(a.iter()).zip(b.iter()) {
        *o = *aa - scale * (*bb);
    }
}

fn jacobi_precondition(diag: &[f64], v: &[f64], out: &mut [f64]) {
    for ((o, vv), d) in out.iter_mut().zip(v.iter()).zip(diag.iter()) {
        *o = *vv / *d;
    }
}

pub(crate) fn solve_bicgstab_left_jacobi_dense(
    a: &DenseMatrix,
    rhs: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<(Vec<f64>, BiCgStabDiagnostics), String> {
    let n = a.size();
    if rhs.len() != n {
        return Err("RHS length mismatch in dense BiCGSTAB solve.".into());
    }
    if n == 0 {
        return Err("Empty system in dense BiCGSTAB solve.".into());
    }
    if tol <= 0.0 || !tol.is_finite() {
        return Err("Invalid tolerance in dense BiCGSTAB solve.".into());
    }

    let diag = a.diagonal()?;
    if diag.iter().any(|d| !d.is_finite() || *d == 0.0) {
        return Err("Invalid diagonal for Jacobi preconditioner in BiCGSTAB.".into());
    }

    // Solve A' x = b' with A' = M^{-1}A and b' = M^{-1}b (left Jacobi).
    let mut b_prime = vec![0.0; n];
    jacobi_precondition(&diag, rhs, &mut b_prime);
    let b_norm = l2_norm(&b_prime).max(1e-300);

    let mut x = vec![0.0; n];
    let mut r = b_prime.clone();
    let r0 = r.clone();

    let mut p = vec![0.0; n];
    let mut v = vec![0.0; n];
    let mut s = vec![0.0; n];
    let mut t = vec![0.0; n];

    let mut rho_old = 1.0_f64;
    let mut alpha = 1.0_f64;
    let mut omega = 1.0_f64;

    let mut rel = l2_norm(&r) / b_norm;
    if rel <= tol {
        return Ok((
            x,
            BiCgStabDiagnostics {
                iterations: 0,
                final_residual_relative_l2: rel,
            },
        ));
    }

    for iter in 0..max_iter {
        let rho_new = dot(&r0, &r);
        if rho_new.abs() < 1e-300 {
            return Err("BiCGSTAB breakdown (rho ~ 0).".into());
        }

        let beta = (rho_new / rho_old) * (alpha / omega);
        for i in 0..n {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }

        // v = A' p = M^{-1}(A p)
        let ap = a.matvec(&p)?;
        jacobi_precondition(&diag, &ap, &mut v);

        let denom = dot(&r0, &v);
        if denom.abs() < 1e-300 {
            return Err("BiCGSTAB breakdown (r0·v ~ 0).".into());
        }
        alpha = rho_new / denom;

        sub_scaled(&mut s, &r, alpha, &v);
        rel = l2_norm(&s) / b_norm;
        if rel <= tol {
            axpy_in_place(&mut x, alpha, &p);
            return Ok((
                x,
                BiCgStabDiagnostics {
                    iterations: iter + 1,
                    final_residual_relative_l2: rel,
                },
            ));
        }

        // t = A' s
        let as_vec = a.matvec(&s)?;
        jacobi_precondition(&diag, &as_vec, &mut t);

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
        sub_scaled(&mut r, &s, omega, &t);

        rel = l2_norm(&r) / b_norm;
        if rel <= tol {
            return Ok((
                x,
                BiCgStabDiagnostics {
                    iterations: iter + 1,
                    final_residual_relative_l2: rel,
                },
            ));
        }

        rho_old = rho_new;
    }

    Err(format!(
        "BiCGSTAB did not converge within max_iter (final rel residual {}).",
        rel
    ))
}
