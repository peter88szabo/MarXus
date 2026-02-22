use super::linear_algebra::DenseMatrix;

pub struct LdltDiagnostics {
    pub min_pivot_abs: f64,
    pub two_by_two_pivot_count: usize,
}

fn solve_2x2(a11: f64, a21: f64, a22: f64, b1: f64, b2: f64) -> Result<(f64, f64), String> {
    let det = a11 * a22 - a21 * a21;
    if !det.is_finite() || det.abs() < 1e-300 {
        return Err("Singular/ill-conditioned 2x2 pivot in LDLT.".into());
    }
    let x1 = (a22 * b1 - a21 * b2) / det;
    let x2 = (-a21 * b1 + a11 * b2) / det;
    Ok((x1, x2))
}

pub(crate) fn solve_symmetric_indefinite_ldlt_bunch_kaufman(
    a_in: &DenseMatrix,
    rhs_in: &[f64],
) -> Result<(Vec<f64>, LdltDiagnostics), String> {
    let n = a_in.size();
    if rhs_in.len() != n {
        return Err("RHS length mismatch in LDLT solve.".into());
    }
    if n == 0 {
        return Err("Empty system in LDLT solve.".into());
    }

    // Work on a mutable copy.
    let mut a = a_in.clone();
    let mut rhs = rhs_in.to_vec();
    let mut perm: Vec<usize> = (0..n).collect();

    // Track pivot blocks: block_size[k] = 1 or 2 (2 means block starts at k).
    let mut block_size = vec![1usize; n];

    // Bunch–Kaufman threshold.
    let alpha = (1.0 + 17.0_f64.sqrt()) / 8.0;
    let mut min_pivot_abs = f64::INFINITY;
    let mut two_by_two = 0usize;

    let mut k = 0usize;
    while k < n {
        // Choose pivot (1x1 or 2x2) with possible symmetric permutation.
        let mut imax = k;
        let mut colmax = 0.0_f64;
        for i in (k + 1)..n {
            let v = a.get(i, k).abs();
            if v > colmax {
                colmax = v;
                imax = i;
            }
        }

        let akk = a.get(k, k).abs();
        if colmax == 0.0 || akk >= alpha * colmax {
            // 1x1 pivot at k
        } else {
            // Consider pivot at imax or 2x2
            let mut rowmax = 0.0_f64;
            for j in k..n {
                if j == imax {
                    continue;
                }
                rowmax = rowmax.max(a.get(imax, j).abs());
            }

            let aii = a.get(imax, imax).abs();
            if aii >= alpha * rowmax {
                // 1x1 pivot at imax: swap k <-> imax
                a.swap_rows_cols_symmetric(k, imax);
                rhs.swap(k, imax);
                perm.swap(k, imax);
            } else {
                // 2x2 pivot using (k, imax): bring imax to k+1 if needed.
                if k + 1 >= n {
                    return Err("LDLT needs a 2x2 pivot at last index.".into());
                }
                if imax != k + 1 {
                    a.swap_rows_cols_symmetric(k + 1, imax);
                    rhs.swap(k + 1, imax);
                    perm.swap(k + 1, imax);
                }
                block_size[k] = 2;
            }
        }

        if block_size[k] == 1 {
            let d = a.get(k, k);
            if !d.is_finite() || d.abs() < 1e-300 {
                return Err("Zero/invalid pivot encountered in LDLT.".into());
            }
            min_pivot_abs = min_pivot_abs.min(d.abs());

            // Compute L column k and update trailing submatrix.
            for i in (k + 1)..n {
                let lik = a.get(i, k) / d;
                a.set(i, k, lik);
            }

            for i in (k + 1)..n {
                let lik = a.get(i, k);
                for j in i..n {
                    let ljk = a.get(j, k);
                    let new = a.get(i, j) - lik * d * ljk;
                    a.set(i, j, new);
                    a.set(j, i, new);
                }
            }

            k += 1;
        } else {
            // 2x2 pivot block at k,k+1
            let a11 = a.get(k, k);
            let a21 = a.get(k + 1, k);
            let a22 = a.get(k + 1, k + 1);
            let det = a11 * a22 - a21 * a21;
            if !det.is_finite() || det.abs() < 1e-300 {
                return Err("Singular/invalid 2x2 pivot encountered in LDLT.".into());
            }
            min_pivot_abs = min_pivot_abs.min(det.abs().sqrt());
            two_by_two += 1;

            // Compute L columns k and k+1 below the pivot block.
            for i in (k + 2)..n {
                let b1 = a.get(i, k);
                let b2 = a.get(i, k + 1);
                let (l1, l2) = solve_2x2(a11, a21, a22, b1, b2)?;
                a.set(i, k, l1);
                a.set(i, k + 1, l2);
            }

            // Update trailing submatrix A22 := A22 - L21 * D * L21^T.
            for i in (k + 2)..n {
                let li1 = a.get(i, k);
                let li2 = a.get(i, k + 1);
                for j in i..n {
                    let lj1 = a.get(j, k);
                    let lj2 = a.get(j, k + 1);
                    let corr =
                        a11 * li1 * lj1 + a21 * (li1 * lj2 + li2 * lj1) + a22 * li2 * lj2;
                    let new = a.get(i, j) - corr;
                    a.set(i, j, new);
                    a.set(j, i, new);
                }
            }

            k += 2;
        }
    }

    // Forward solve: L y = rhs (L unit lower, stored in strict lower part of a).
    let mut y = rhs;
    for i in 0..n {
        let mut sum = y[i];
        for j in 0..i {
            sum -= a.get(i, j) * y[j];
        }
        y[i] = sum;
    }

    // Solve D z = y (D is block diagonal, entries remain on diagonal/subdiagonal).
    let mut z = vec![0.0; n];
    let mut i = 0usize;
    while i < n {
        if block_size[i] == 1 {
            let d = a.get(i, i);
            if !d.is_finite() || d.abs() < 1e-300 {
                return Err("Invalid diagonal in LDLT backsolve.".into());
            }
            z[i] = y[i] / d;
            i += 1;
        } else {
            let a11 = a.get(i, i);
            let a21 = a.get(i + 1, i);
            let a22 = a.get(i + 1, i + 1);
            let (z1, z2) = solve_2x2(a11, a21, a22, y[i], y[i + 1])?;
            z[i] = z1;
            z[i + 1] = z2;
            i += 2;
        }
    }

    // Back solve: L^T x = z
    let mut x = z;
    for i_rev in 0..n {
        let i = n - 1 - i_rev;
        let mut sum = x[i];
        for j in (i + 1)..n {
            sum -= a.get(j, i) * x[j];
        }
        x[i] = sum;
    }

    // Unpermute to original ordering.
    let mut x_out = vec![0.0; n];
    for i in 0..n {
        x_out[perm[i]] = x[i];
    }

    Ok((
        x_out,
        LdltDiagnostics {
            min_pivot_abs: min_pivot_abs.max(0.0),
            two_by_two_pivot_count: two_by_two,
        },
    ))
}

