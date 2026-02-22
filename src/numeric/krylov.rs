pub trait LinearOperator {
    fn dim(&self) -> usize;
    fn matvec(&mut self, x: &[f64], y: &mut [f64]) -> Result<(), String>;
}

pub trait Preconditioner {
    fn dim(&self) -> usize;
    fn apply(&self, rhs: &[f64], out: &mut [f64]) -> Result<(), String>;
}

pub struct IdentityPreconditioner {
    n: usize,
}

impl IdentityPreconditioner {
    pub fn new(n: usize) -> Self {
        Self { n }
    }
}

impl Preconditioner for IdentityPreconditioner {
    fn dim(&self) -> usize {
        self.n
    }

    fn apply(&self, rhs: &[f64], out: &mut [f64]) -> Result<(), String> {
        if rhs.len() != self.n || out.len() != self.n {
            return Err("Dimension mismatch in IdentityPreconditioner.".into());
        }
        out.copy_from_slice(rhs);
        Ok(())
    }
}

pub struct JacobiPreconditioner {
    inv_diag: Vec<f64>,
}

impl JacobiPreconditioner {
    pub fn from_diagonal(diag: &[f64]) -> Result<Self, String> {
        if diag.is_empty() {
            return Err("Empty diagonal for Jacobi preconditioner.".into());
        }
        let mut inv = vec![0.0; diag.len()];
        for (i, (o, d)) in inv.iter_mut().zip(diag.iter()).enumerate() {
            if !d.is_finite() || d.abs() < 1e-300 {
                return Err(format!(
                    "Invalid Jacobi diagonal entry at i={} (value={}).",
                    i, d
                ));
            }
            *o = 1.0 / *d;
        }
        Ok(Self { inv_diag: inv })
    }
}

impl Preconditioner for JacobiPreconditioner {
    fn dim(&self) -> usize {
        self.inv_diag.len()
    }

    fn apply(&self, rhs: &[f64], out: &mut [f64]) -> Result<(), String> {
        let n = self.inv_diag.len();
        if rhs.len() != n || out.len() != n {
            return Err("Dimension mismatch in Jacobi preconditioner.".into());
        }
        for i in 0..n {
            out[i] = rhs[i] * self.inv_diag[i];
        }
        Ok(())
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn l2_norm(v: &[f64]) -> f64 {
    v.iter().copied().map(|x| x * x).sum::<f64>().sqrt()
}

fn axpy(y: &mut [f64], a: f64, x: &[f64]) {
    for (yy, xx) in y.iter_mut().zip(x.iter()) {
        *yy += a * (*xx);
    }
}

fn sub_scaled(out: &mut [f64], a: &[f64], scale: f64, b: &[f64]) {
    for ((o, aa), bb) in out.iter_mut().zip(a.iter()).zip(b.iter()) {
        *o = *aa - scale * (*bb);
    }
}

pub struct KrylovDiagnostics {
    pub iterations: usize,
    pub final_residual_relative_l2: f64,
}

pub fn solve_bicgstab_left_preconditioned(
    op: &mut dyn LinearOperator,
    m_inv: &dyn Preconditioner,
    rhs: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<(Vec<f64>, KrylovDiagnostics), String> {
    let n = op.dim();
    if rhs.len() != n || m_inv.dim() != n {
        return Err("Dimension mismatch in BiCGSTAB.".into());
    }
    if n == 0 {
        return Err("Empty system in BiCGSTAB.".into());
    }
    if tol <= 0.0 || !tol.is_finite() {
        return Err("Invalid tolerance in BiCGSTAB.".into());
    }

    // Left preconditioning: solve A' x = b', with A' = M^{-1}A and b' = M^{-1}b.
    let mut b_prime = vec![0.0; n];
    m_inv.apply(rhs, &mut b_prime)?;
    let b_norm = l2_norm(&b_prime).max(1e-300);

    let mut x = vec![0.0; n];

    let mut r = b_prime.clone();
    let r0 = r.clone();

    let mut p = vec![0.0; n];
    let mut v = vec![0.0; n];
    let mut s = vec![0.0; n];
    let mut t = vec![0.0; n];

    let mut ap = vec![0.0; n];
    let mut as_vec = vec![0.0; n];

    let mut rho_old = 1.0_f64;
    let mut alpha = 1.0_f64;
    let mut omega = 1.0_f64;

    let mut rel = l2_norm(&r) / b_norm;
    if rel <= tol {
        return Ok((
            x,
            KrylovDiagnostics {
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
        op.matvec(&p, &mut ap)?;
        m_inv.apply(&ap, &mut v)?;

        let denom = dot(&r0, &v);
        if denom.abs() < 1e-300 {
            return Err("BiCGSTAB breakdown (r0·v ~ 0).".into());
        }
        alpha = rho_new / denom;

        sub_scaled(&mut s, &r, alpha, &v);
        rel = l2_norm(&s) / b_norm;
        if rel <= tol {
            axpy(&mut x, alpha, &p);
            return Ok((
                x,
                KrylovDiagnostics {
                    iterations: iter + 1,
                    final_residual_relative_l2: rel,
                },
            ));
        }

        // t = A' s
        op.matvec(&s, &mut as_vec)?;
        m_inv.apply(&as_vec, &mut t)?;

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
                KrylovDiagnostics {
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

pub fn solve_gmres_left_preconditioned_restarted(
    op: &mut dyn LinearOperator,
    m_inv: &dyn Preconditioner,
    rhs: &[f64],
    tol: f64,
    max_iter: usize,
    restart: usize,
) -> Result<(Vec<f64>, KrylovDiagnostics), String> {
    let n = op.dim();
    if rhs.len() != n || m_inv.dim() != n {
        return Err("Dimension mismatch in GMRES.".into());
    }
    if n == 0 {
        return Err("Empty system in GMRES.".into());
    }
    if tol <= 0.0 || !tol.is_finite() {
        return Err("Invalid tolerance in GMRES.".into());
    }
    if restart == 0 {
        return Err("GMRES restart must be > 0.".into());
    }

    let mut x = vec![0.0; n];

    // Precondition RHS: b' = M^{-1} b
    let mut b_prime = vec![0.0; n];
    m_inv.apply(rhs, &mut b_prime)?;
    let b_norm = l2_norm(&b_prime).max(1e-300);

    let mut total_iter = 0usize;

    // r = b' - A' x, but x=0 => r=b'
    let mut r = b_prime.clone();
    let mut rel = l2_norm(&r) / b_norm;
    if rel <= tol {
        return Ok((
            x,
            KrylovDiagnostics {
                iterations: 0,
                final_residual_relative_l2: rel,
            },
        ));
    }

    // Workspace
    let m = restart.min(max_iter.max(1));
    let mut v: Vec<Vec<f64>> = (0..(m + 1)).map(|_| vec![0.0; n]).collect();
    let mut h: Vec<Vec<f64>> = (0..(m + 1)).map(|_| vec![0.0; m]).collect();
    let mut cs = vec![0.0; m];
    let mut sn = vec![0.0; m];
    let mut g = vec![0.0; m + 1];
    let mut w = vec![0.0; n];
    let mut tmp = vec![0.0; n];
    let mut y = vec![0.0; m];

    while total_iter < max_iter {
        // v0 = r / beta
        let beta = l2_norm(&r);
        if beta <= 0.0 || !beta.is_finite() {
            return Err("GMRES residual norm invalid.".into());
        }
        for i in 0..n {
            v[0][i] = r[i] / beta;
        }
        for x in &mut g {
            *x = 0.0;
        }
        g[0] = beta;

        let mut inner = 0usize;
        for j in 0..m {
            // w = A' vj = M^{-1}(A vj)
            op.matvec(&v[j], &mut tmp)?;
            m_inv.apply(&tmp, &mut w)?;

            // Modified Gram-Schmidt
            for i in 0..=j {
                h[i][j] = dot(&w, &v[i]);
                axpy(&mut w, -h[i][j], &v[i]);
            }
            h[j + 1][j] = l2_norm(&w);
            if h[j + 1][j] > 0.0 {
                for i in 0..n {
                    v[j + 1][i] = w[i] / h[j + 1][j];
                }
            } else {
                for i in 0..n {
                    v[j + 1][i] = 0.0;
                }
            }

            // Apply existing Givens rotations to the new column
            for i in 0..j {
                let temp = cs[i] * h[i][j] + sn[i] * h[i + 1][j];
                h[i + 1][j] = -sn[i] * h[i][j] + cs[i] * h[i + 1][j];
                h[i][j] = temp;
            }

            // Create and apply new Givens rotation
            let (c, s) = givens(h[j][j], h[j + 1][j]);
            cs[j] = c;
            sn[j] = s;
            let temp = c * h[j][j] + s * h[j + 1][j];
            h[j + 1][j] = 0.0;
            h[j][j] = temp;

            let temp_g = c * g[j] + s * g[j + 1];
            g[j + 1] = -s * g[j] + c * g[j + 1];
            g[j] = temp_g;

            inner = j + 1;
            total_iter += 1;

            rel = g[inner].abs() / b_norm;
            if rel <= tol || total_iter >= max_iter {
                break;
            }
        }

        // Solve upper-triangular system for y
        for i in 0..inner {
            y[i] = g[i];
        }
        for i_rev in 0..inner {
            let i = inner - 1 - i_rev;
            let mut sum = y[i];
            for k in (i + 1)..inner {
                sum -= h[i][k] * y[k];
            }
            if h[i][i].abs() < 1e-300 || !h[i][i].is_finite() {
                return Err("GMRES encountered near-singular Hessenberg diagonal.".into());
            }
            y[i] = sum / h[i][i];
        }

        // Update x = x + V(:,0..inner-1) y
        for i in 0..n {
            let mut acc = 0.0;
            for j in 0..inner {
                acc += v[j][i] * y[j];
            }
            x[i] += acc;
        }

        // Recompute residual r = b' - A' x
        op.matvec(&x, &mut tmp)?;
        m_inv.apply(&tmp, &mut w)?;
        for i in 0..n {
            r[i] = b_prime[i] - w[i];
        }

        rel = l2_norm(&r) / b_norm;
        if rel <= tol {
            return Ok((
                x,
                KrylovDiagnostics {
                    iterations: total_iter,
                    final_residual_relative_l2: rel,
                },
            ));
        }
    }

    Err(format!(
        "GMRES did not converge within max_iter (final rel residual {}).",
        rel
    ))
}

fn givens(a: f64, b: f64) -> (f64, f64) {
    if b == 0.0 {
        return (1.0, 0.0);
    }
    if a == 0.0 {
        return (0.0, 1.0);
    }
    let r = (a * a + b * b).sqrt();
    (a / r, b / r)
}
