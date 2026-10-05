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
            if akk * rowmax >= alpha * colmax * colmax {
                // 1x1 pivot at k without interchange
                // (Bunch & Kaufman, Math. Comp. 31, 163 (1977); LAPACK dsytf2)
            } else if aii >= alpha * rowmax {
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
                    let corr = a11 * li1 * lj1 + a21 * (li1 * lj2 + li2 * lj1) + a22 * li2 * lj2;
                    let new = a.get(i, j) - corr;
                    a.set(i, j, new);
                    a.set(j, i, new);
                }
            }

            k += 2;
        }
    }

    // Inside a 2x2 pivot block (rows k, k+1) L is the identity: a(k+1, k) is the off-diagonal
    // element of D, not of L, and must be skipped in both triangular solves.
    let in_same_block = |row: usize, col: usize| col + 1 == row && block_size[col] == 2;

    // Forward solve: L y = rhs (L unit lower, stored in strict lower part of a).
    let mut y = rhs;
    for i in 0..n {
        let mut sum = y[i];
        for j in 0..i {
            if in_same_block(i, j) {
                continue;
            }
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
            if in_same_block(j, i) {
                continue;
            }
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

#[cfg(test)]
mod tests {
    use super::*;

    fn dense(rows: &[&[f64]]) -> DenseMatrix {
        let n = rows.len();
        let mut a = DenseMatrix::zeros(n);
        for i in 0..n {
            for j in 0..n {
                a.set(i, j, rows[i][j]);
            }
        }
        a
    }

    fn assert_solves(a: &DenseMatrix, x_ref: &[f64]) {
        let n = x_ref.len();
        let b: Vec<f64> = (0..n).map(|i| (0..n).map(|j| a.get(i, j) * x_ref[j]).sum()).collect();
        let (x, _) = solve_symmetric_indefinite_ldlt_bunch_kaufman(a, &b).unwrap();
        for i in 0..n {
            assert!((x[i] - x_ref[i]).abs() < 1e-10, "x[{i}] = {}, expected {}", x[i], x_ref[i]);
        }
    }

    #[test]
    fn two_by_two_pivot_block_is_solved_correctly() {
        // Zero diagonal forces a 2x2 pivot; inside the block L is the identity.
        assert_solves(&dense(&[&[0.0, 1.0], &[1.0, 0.0]]), &[2.0, 1.0]);
    }

    #[test]
    fn bunch_kaufman_takes_1x1_pivot_when_akk_times_rowmax_is_large() {
        // |a00| * rowmax = 0.5 * 4 >= alpha * colmax^2 = 0.64 * 1: 1x1 pivot without interchange
        // (Bunch & Kaufman, Math. Comp. 31, 163 (1977)). The 2x2 block of rows/cols 0,1 is singular.
        let a = dense(&[&[0.5, 1.0, 0.0], &[1.0, 2.0, 4.0], &[0.0, 4.0, 1.0]]);
        assert_solves(&a, &[1.0, 2.0, 3.0]);
    }

    #[test]
    fn random_symmetric_indefinite_system() {
        let n = 8;
        let mut seed: u64 = 12345;
        let mut next = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((seed >> 11) as f64 / (1u64 << 53) as f64) - 0.5
        };
        let mut a = DenseMatrix::zeros(n);
        for i in 0..n {
            for j in 0..=i {
                let v = next();
                a.set(i, j, v);
                a.set(j, i, v);
            }
        }
        let x_ref: Vec<f64> = (0..n).map(|i| 1.0 + i as f64).collect();
        assert_solves(&a, &x_ref);
    }
}
