//! Eigenvalues and eigenvectors of real symmetric matrices.
//!
//! The master equation is solved through the eigenvalues and eigenvectors of the symmetrized relaxation
//! matrix as in Olzmann's work, where they "were calculated by using the EISPACK routine tql2 after
//! symmetrization" (Pfeifle, Olzmann, Int. J. Chem. Kinet. 46, 231 (2014), p. 235; Olzmann, Phys. Chem.
//! Chem. Phys. 4, 3614 (2002): routine tql1 "applied after symmetrization of J"). The algorithms are the
//! classical ones of these references:
//!   - Householder reduction of a symmetric matrix to tridiagonal form with accumulation of the
//!     transformations (tred2);
//!   - QL algorithm with implicit shifts for the symmetric tridiagonal matrix, transforming the
//!     accumulated eigenvectors (tql2 / tqli);
//! B. T. Smith et al., Matrix Eigensystem Routines - EISPACK Guide (Springer, 1976);
//! W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery, Numerical Recipes in Fortran,
//! 2nd ed. (Cambridge University Press, 1992), Sec. 11.2-11.3.
//!
//! For the lowest eigenpairs of a large symmetric positive definite band matrix, inverse iteration with
//! its banded Cholesky factor is used (Numerical Recipes Sec. 11.7): x <- S^-1 x converges to the
//! eigenvector of the smallest eigenvalue with the ratio lambda_1/lambda_2 per step, fastest when
//! lambda_1 is well separated, i.e. when S is nearly singular. The factor is that of S + sigma I with a
//! small shift sigma (inverse iteration with a shift, NR92 Sec. 11.7), which keeps the factorization
//! computable when lambda_1 lies below the double-precision resolution (`cholesky_safety_shift`).
//!
//! References:
//!   NR92: W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery, Numerical Recipes in Fortran,
//!         2nd ed. (Cambridge University Press, 1992), Sec. 11.2-11.3 (tred2, tqli), Sec. 11.7 (inverse
//!         iteration, with a shift).
//!   B. T. Smith et al., Matrix Eigensystem Routines - EISPACK Guide, 2nd ed. (Springer, 1976).
//!   J. H. Wilkinson, The Algebraic Eigenvalue Problem (Clarendon Press, Oxford, 1965): backward stability of
//!         the Householder and QL reductions.
//!   N. J. Higham, Accuracy and Stability of Numerical Algorithms, 2nd ed. (SIAM, Philadelphia, 2002),
//!         ch. 10 (Cholesky factorization: backward error, condition for the factorization to run to
//!         completion).
//!   J. W. Demmel, On floating point errors in Cholesky, LAPACK Working Note 14 (1989).
//!   H. Weyl, Math. Ann. 71, 441 (1912): eigenvalue perturbation bound |lambda(S + dS) - lambda(S)| <= ||dS||.
//!   C. Davis, W. M. Kahan, SIAM J. Numer. Anal. 7, 1 (1970): eigenvector perturbation, of order
//!         ||dS|| / (gap to the other eigenvalues).

use super::banded_solvers::SymmetricBandMatrix;

/// Eigenvalues (ascending) and orthonormal eigenvectors (vectors[k] belongs to values[k]) of the
/// symmetric matrix given by its rows (only the lower triangle is used).
pub fn symmetric_eigen(matrix: &[Vec<f64>]) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
    let n = matrix.len();
    if matrix.iter().any(|row| row.len() != n) {
        return Err("symmetric_eigen: the matrix is not square.".into());
    }
    if n == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    // Work on a full copy built from the lower triangle.
    let mut a = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..=i {
            a[i][j] = matrix[i][j];
            a[j][i] = matrix[i][j];
        }
    }
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n];
    householder_tridiagonalize(&mut a, &mut d, &mut e);
    tridiagonal_ql_implicit(&mut d, &mut e, &mut a)?;
    Ok(sorted_eigenpairs(d, &a))
}

/// Householder reduction to tridiagonal form (tred2). On return `d` holds the diagonal, `e[i]` the
/// element coupling i - 1 and i (e[0] = 0), and `a` the orthogonal transformation (columns).
fn householder_tridiagonalize(a: &mut [Vec<f64>], d: &mut [f64], e: &mut [f64]) {
    let n = d.len();
    for i in (1..n).rev() {
        let l = i - 1;
        let mut h = 0.0;
        if l > 0 {
            let scale: f64 = (0..=l).map(|k| a[i][k].abs()).sum();
            if scale == 0.0 {
                e[i] = a[i][l];
            } else {
                for k in 0..=l {
                    a[i][k] /= scale;
                    h += a[i][k] * a[i][k];
                }
                let f = a[i][l];
                let g = if f >= 0.0 { -h.sqrt() } else { h.sqrt() };
                e[i] = scale * g;
                h -= f * g;
                a[i][l] = f - g;
                let mut f = 0.0;
                for j in 0..=l {
                    a[j][i] = a[i][j] / h;
                    let mut g = 0.0;
                    for k in 0..=j {
                        g += a[j][k] * a[i][k];
                    }
                    for k in (j + 1)..=l {
                        g += a[k][j] * a[i][k];
                    }
                    e[j] = g / h;
                    f += e[j] * a[i][j];
                }
                let hh = f / (h + h);
                for j in 0..=l {
                    let f = a[i][j];
                    let g = e[j] - hh * f;
                    e[j] = g;
                    for k in 0..=j {
                        a[j][k] -= f * e[k] + g * a[i][k];
                    }
                }
            }
        } else {
            e[i] = a[i][l];
        }
        d[i] = h;
    }
    d[0] = 0.0;
    e[0] = 0.0;
    // Accumulation of the transformations.
    for i in 0..n {
        if d[i] != 0.0 {
            for j in 0..i {
                let g: f64 = (0..i).map(|k| a[i][k] * a[k][j]).sum();
                for k in 0..i {
                    a[k][j] -= g * a[k][i];
                }
            }
        }
        d[i] = a[i][i];
        a[i][i] = 1.0;
        for j in 0..i {
            a[j][i] = 0.0;
            a[i][j] = 0.0;
        }
    }
}

/// QL algorithm with implicit shifts (tql2/tqli) for the symmetric tridiagonal matrix with diagonal `d`
/// and off-diagonal `e` (e[i] couples i - 1 and i, e[0] unused); the eigenvectors are accumulated into
/// the columns of `z` (identity for the tridiagonal matrix itself, the Householder transformation for
/// a full matrix). On return `d` holds the eigenvalues (unsorted).
fn tridiagonal_ql_implicit(d: &mut [f64], e: &mut [f64], z: &mut [Vec<f64>]) -> Result<(), String> {
    let n = d.len();
    for i in 1..n {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;
    for l in 0..n {
        let mut iterations = 0;
        loop {
            // Look for a negligible off-diagonal element to split the matrix.
            let mut m = l;
            while m + 1 < n {
                let dd = d[m].abs() + d[m + 1].abs();
                if e[m].abs() <= f64::EPSILON * dd {
                    break;
                }
                m += 1;
            }
            if m == l {
                break;
            }
            iterations += 1;
            if iterations > 60 {
                return Err("QL algorithm: no convergence after 60 iterations.".into());
            }
            let mut g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let mut r = g.hypot(1.0);
            g = d[m] - d[l] + e[l] / (g + if g >= 0.0 { r.abs() } else { -r.abs() });
            let (mut s, mut c, mut p) = (1.0, 1.0, 0.0);
            let mut deflated_early = false;
            let mut i = m as isize - 1;
            while i >= l as isize {
                let iu = i as usize;
                let f = s * e[iu];
                let b = c * e[iu];
                r = f.hypot(g);
                e[iu + 1] = r;
                if r == 0.0 {
                    d[iu + 1] -= p;
                    e[m] = 0.0;
                    deflated_early = true;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[iu + 1] - p;
                r = (d[iu] - g) * s + 2.0 * c * b;
                p = s * r;
                d[iu + 1] = g + p;
                g = c * r - b;
                for row in z.iter_mut() {
                    let f = row[iu + 1];
                    row[iu + 1] = s * row[iu] + c * f;
                    row[iu] = c * row[iu] - s * f;
                }
                i -= 1;
            }
            if deflated_early {
                continue;
            }
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        }
    }
    Ok(())
}

/// Eigenvalues ascending with the eigenvectors (columns of z) as rows.
fn sorted_eigenpairs(d: Vec<f64>, z: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = d.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| d[a].partial_cmp(&d[b]).unwrap_or(std::cmp::Ordering::Equal));
    let values = order.iter().map(|&k| d[k]).collect();
    let vectors = order.iter().map(|&k| (0..n).map(|i| z[i][k]).collect()).collect();
    (values, vectors)
}

/// Eigenvalues (ascending) and orthonormal eigenvectors of the symmetric tridiagonal matrix with the
/// diagonal `diagonal` and the off-diagonal `off_diagonal` (off_diagonal[i] couples i and i + 1).
pub fn symmetric_tridiagonal_eigen(diagonal: &[f64], off_diagonal: &[f64]) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
    let n = diagonal.len();
    if n == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    if off_diagonal.len() + 1 != n {
        return Err("symmetric_tridiagonal_eigen: the off-diagonal must have n - 1 elements.".into());
    }
    let mut d = diagonal.to_vec();
    let mut e = vec![0.0; n];
    e[1..].copy_from_slice(off_diagonal);
    let mut z: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect()).collect();
    tridiagonal_ql_implicit(&mut d, &mut e, &mut z)?;
    Ok(sorted_eigenpairs(d, &z))
}

/// Shift sigma = n eps max|S_ii| that keeps the Cholesky factor of S + sigma I computable when S is
/// positive (semi)definite but its lowest eigenvalue lies below the double-precision resolution.
///
/// The Cholesky factorization runs to completion when the smallest eigenvalue of the matrix, relative to
/// its diagonal, exceeds a multiple of n eps (Higham 2002, ch. 10; Demmel 1989; see the module
/// references); without the shift, a lowest eigenvalue below eps S_ii turns pivots into rounding noise.
/// sigma is far below the second eigenvalue in a master equation, so the inverse iteration converges as
/// fast as without it (ratio (lambda_1 + sigma)/(lambda_2 + sigma)), and the eigenvectors are those of S.
/// Not suitable when several lowest eigenvalues lie within sigma of each other (e.g. a reducible operator):
/// the iteration cannot separate them.
pub fn cholesky_safety_shift(matrix: &SymmetricBandMatrix) -> f64 {
    matrix.dimension() as f64 * f64::EPSILON * matrix.max_abs_diagonal()
}

/// The `count` lowest eigenpairs of a symmetric positive definite band matrix by inverse iteration
/// with its Cholesky factor; each further eigenvector is kept orthogonal to the previous ones
/// (deflation). `start` is the initial vector of the first iteration. Converged when the relative change
/// of the eigenvalue and the change of the unit eigenvector are both below `relative_tolerance`.
pub fn lowest_eigenpairs_banded(
    matrix: &SymmetricBandMatrix,
    count: usize,
    start: &[f64],
    shift: f64,
    relative_tolerance: f64,
    max_iterations: usize,
) -> Result<Vec<(f64, Vec<f64>)>, String> {
    let n = matrix.dimension();
    if start.len() != n || count == 0 || count > n {
        return Err(format!("Inverse iteration: start vector of length {} for dimension {n}, {count} eigenpairs.", start.len()));
    }
    // Inverse iteration with a shift (NR92 Sec. 11.7): (S + shift I) has the eigenvectors of S and the
    // eigenvalues lambda + shift; the iteration converges with the ratio (lambda_1 + shift)/(lambda_2 + shift).
    let factor = matrix.shifted(shift).cholesky()?;
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    let orthogonalize = |x: &mut Vec<f64>, found: &[(f64, Vec<f64>)]| {
        for (_, v) in found {
            let c = dot(x, v);
            x.iter_mut().zip(v).for_each(|(a, b)| *a -= c * b);
        }
    };
    let normalize = |x: &mut Vec<f64>| -> f64 {
        let norm = dot(x, x).sqrt();
        if norm > 0.0 {
            x.iter_mut().for_each(|a| *a /= norm);
        }
        norm
    };

    let mut found: Vec<(f64, Vec<f64>)> = Vec::with_capacity(count);
    for k in 0..count {
        let mut x = start.to_vec();
        orthogonalize(&mut x, &found);
        if normalize(&mut x) < 1e-8 {
            // The start vector lies in the span of the eigenvectors already found: use an alternating one.
            x = (0..n).map(|i| if (i + k) % 2 == 0 { 1.0 } else { -0.5 }).collect();
            orthogonalize(&mut x, &found);
            normalize(&mut x);
        }
        let mut lambda = f64::NAN;
        let mut converged = false;
        for _ in 0..max_iterations {
            // y = (S + shift I)^-1 x; for an eigenvector x, x.y = 1/(lambda + shift).
            let mut y = factor.solve(&x)?;
            orthogonalize(&mut y, &found);
            let new_lambda = 1.0 / dot(&x, &y) - shift;
            normalize(&mut y);
            // Vector change measured directly (1 - |x.y| ~ theta^2/2 would hide an error theta of
            // sqrt(tolerance)); the sign of y is aligned with x first.
            if dot(&x, &y) < 0.0 {
                y.iter_mut().for_each(|a| *a = -*a);
            }
            let change = x.iter().zip(&y).map(|(a, b)| (a - b) * (a - b)).sum::<f64>().sqrt();
            x = y;
            // The eigenvalue change is measured relative to the shifted eigenvalue, the quantity that is
            // actually iterated (a lambda far below the shift is resolved only to about the shift times the
            // tolerance).
            if lambda.is_finite()
                && (new_lambda - lambda).abs() <= relative_tolerance * (new_lambda + shift).abs()
                && change <= relative_tolerance
            {
                lambda = new_lambda;
                converged = true;
                break;
            }
            lambda = new_lambda;
        }
        if !converged {
            return Err(format!("Inverse iteration: eigenpair {k} did not converge in {max_iterations} iterations."));
        }
        // Sign convention: the component of largest magnitude is positive.
        let largest = x.iter().cloned().fold(0.0_f64, |m, v| if v.abs() > m.abs() { v } else { m });
        if largest < 0.0 {
            x.iter_mut().for_each(|a| *a = -*a);
        }
        found.push((lambda, x));
    }
    Ok(found)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pseudo_random(seed: &mut u64) -> f64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64 / (1u64 << 53) as f64) - 0.5
    }

    fn random_symmetric(n: usize, seed: u64) -> Vec<Vec<f64>> {
        let mut s = seed;
        let mut a = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..=i {
                let v = pseudo_random(&mut s);
                a[i][j] = v;
                a[j][i] = v;
            }
        }
        a
    }

    fn check_decomposition(a: &[Vec<f64>], values: &[f64], vectors: &[Vec<f64>], tol: f64) {
        let n = a.len();
        let norm = a.iter().flatten().map(|x| x * x).sum::<f64>().sqrt();
        assert!(values.windows(2).all(|w| w[0] <= w[1]), "eigenvalues not ascending");
        for k in 0..n {
            // ||A v - lambda v|| and orthonormality.
            for i in 0..n {
                let av: f64 = (0..n).map(|j| a[i][j] * vectors[k][j]).sum();
                assert!((av - values[k] * vectors[k][i]).abs() < tol * norm, "residual of eigenpair {k}");
            }
            for l in 0..n {
                let dot: f64 = (0..n).map(|i| vectors[k][i] * vectors[l][i]).sum();
                let expected = if k == l { 1.0 } else { 0.0 };
                assert!((dot - expected).abs() < tol, "orthonormality {k},{l}: {dot}");
            }
        }
    }

    #[test]
    fn eigenvalues_of_a_small_matrix_are_exact() {
        // [[2,1,0],[1,2,1],[0,1,2]]: eigenvalues 2 - sqrt(2), 2, 2 + sqrt(2).
        let a = vec![vec![2.0, 1.0, 0.0], vec![1.0, 2.0, 1.0], vec![0.0, 1.0, 2.0]];
        let (values, vectors) = symmetric_eigen(&a).unwrap();
        let expected = [2.0 - 2f64.sqrt(), 2.0, 2.0 + 2f64.sqrt()];
        for k in 0..3 {
            assert!((values[k] - expected[k]).abs() < 1e-14, "{:?}", values);
        }
        check_decomposition(&a, &values, &vectors, 1e-13);
    }

    #[test]
    fn random_symmetric_matrices_are_diagonalized() {
        for (n, seed) in [(1, 1), (2, 2), (7, 3), (40, 4), (120, 5)] {
            let a = random_symmetric(n, seed);
            let (values, vectors) = symmetric_eigen(&a).unwrap();
            check_decomposition(&a, &values, &vectors, 1e-11);
            // The trace is the sum of the eigenvalues.
            let trace: f64 = (0..n).map(|i| a[i][i]).sum();
            assert!((trace - values.iter().sum::<f64>()).abs() < 1e-10);
        }
    }

    #[test]
    fn tridiagonal_matrices_are_diagonalized_directly() {
        // Second-difference matrix, eigenvalues 2 - 2 cos(k pi/(n+1)).
        let n = 30;
        let (values, vectors) = symmetric_tridiagonal_eigen(&vec![2.0; n], &vec![-1.0; n - 1]).unwrap();
        for k in 0..n {
            let exact = 2.0 - 2.0 * ((k + 1) as f64 * std::f64::consts::PI / (n + 1) as f64).cos();
            assert!((values[k] - exact).abs() < 1e-13);
        }
        let mut dense = vec![vec![0.0; n]; n];
        for i in 0..n {
            dense[i][i] = 2.0;
            if i + 1 < n {
                dense[i][i + 1] = -1.0;
                dense[i + 1][i] = -1.0;
            }
        }
        check_decomposition(&dense, &values, &vectors, 1e-12);
    }

    #[test]
    fn inverse_iteration_finds_the_lowest_eigenpairs_of_a_band_matrix() {
        // Diagonally dominant symmetric band matrix; compare with the dense decomposition.
        let (n, bw) = (80, 4);
        let mut seed = 17u64;
        let mut band = SymmetricBandMatrix::zeros(n, bw);
        let mut dense = vec![vec![0.0; n]; n];
        for i in 0..n {
            for k in i.saturating_sub(bw)..i {
                let v = pseudo_random(&mut seed);
                band.set_lower(i, k, v);
                dense[i][k] = v;
                dense[k][i] = v;
            }
        }
        for i in 0..n {
            let d = 2.5 + 0.1 * i as f64;
            band.set_diagonal(i, d);
            dense[i][i] = d;
        }
        let (values, vectors) = symmetric_eigen(&dense).unwrap();
        let pairs = lowest_eigenpairs_banded(&band, 2, &vec![1.0; n], 0.0, 1e-13, 500).unwrap();
        for k in 0..2 {
            assert!(((pairs[k].0 - values[k]) / values[k]).abs() < 1e-10, "lambda_{k}: {} vs {}", pairs[k].0, values[k]);
            let overlap: f64 = (0..n).map(|i| pairs[k].1[i] * vectors[k][i]).sum();
            assert!((overlap.abs() - 1.0).abs() < 1e-8, "eigenvector {k}: overlap {overlap}");
        }
    }

    #[test]
    fn inverse_iteration_resolves_an_eigenvalue_many_orders_below_the_others() {
        // diag(1e-12, 1, 2, ...) with weak coupling: a nearly singular matrix, the case of a deep well.
        let n = 20;
        let mut band = SymmetricBandMatrix::zeros(n, 1);
        band.set_diagonal(0, 1.0e-12 + 1.0e-8);
        for i in 1..n {
            band.set_diagonal(i, i as f64 + 1.0e-8);
        }
        for i in 1..n {
            band.set_lower(i, i - 1, -1.0e-8);
        }
        let mut dense = vec![vec![0.0; n]; n];
        dense[0][0] = 1.0e-12 + 1.0e-8;
        for i in 1..n {
            dense[i][i] = i as f64 + 1.0e-8;
            dense[i][i - 1] = -1.0e-8;
            dense[i - 1][i] = -1.0e-8;
        }
        let pairs = lowest_eigenpairs_banded(&band, 1, &vec![1.0; n], 0.0, 1e-14, 100).unwrap();
        // lambda_1 = 1e-12 + 1e-8 - 1e-16/(1 - 1e-12) + ..., to second order in the coupling.
        let expected = 1.0e-12 + 1.0e-8 - 1.0e-16 / (1.0 - 1.0e-12);
        assert!(((pairs[0].0 - expected) / expected).abs() < 1e-9, "{} vs {expected}", pairs[0].0);
    }

    #[test]
    fn a_shift_leaves_the_eigenpairs_unchanged() {
        // (S + sigma I) has the eigenvectors of S and the eigenvalues lambda + sigma (NR92 Sec. 11.7).
        let n = 60;
        let bw = 3;
        let mut seed = 77u64;
        let mut band = SymmetricBandMatrix::zeros(n, bw);
        for i in 0..n {
            for k in i.saturating_sub(bw)..i {
                band.set_lower(i, k, pseudo_random(&mut seed));
            }
            band.set_diagonal(i, 2.5 + 0.1 * i as f64);
        }
        let plain = lowest_eigenpairs_banded(&band, 2, &vec![1.0; n], 0.0, 1e-13, 500).unwrap();
        let shifted = lowest_eigenpairs_banded(&band, 2, &vec![1.0; n], 0.3, 1e-13, 500).unwrap();
        for k in 0..2 {
            assert!((shifted[k].0 / plain[k].0 - 1.0).abs() < 1e-12, "{} vs {}", shifted[k].0, plain[k].0);
            let overlap: f64 = (0..n).map(|i| shifted[k].1[i] * plain[k].1[i]).sum();
            assert!((overlap - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn the_safety_shift_scales_with_dimension_precision_and_the_largest_diagonal_element() {
        let mut band = SymmetricBandMatrix::zeros(4, 1);
        for (i, d) in [1.0, -7.0, 3.0, 2.0].iter().enumerate() {
            band.set_diagonal(i, *d);
        }
        assert_eq!(cholesky_safety_shift(&band), 4.0 * f64::EPSILON * 7.0);
    }

    #[test]
    fn a_shift_lets_inverse_iteration_find_the_null_vector_of_a_singular_matrix() {
        // Weighted path Laplacian, sum_i w_i (e_i - e_{i+1})(e_i - e_{i+1})^T with weights spanning 6 orders
        // of magnitude: positive semidefinite with the null vector (1, ..., 1)/sqrt(n), the analogue of a
        // conservative relaxation operator. Its Cholesky factor does not exist in exact arithmetic (zero last
        // pivot); with the safety shift the inverse iteration converges to the null vector. The eigenvector is
        // determined to about eps ||S|| / lambda_2 (Davis, Kahan, SIAM J. Numer. Anal. 7, 1 (1970)), which the
        // 6 orders keep far below the tolerance of the test.
        let n = 40;
        let mut band = SymmetricBandMatrix::zeros(n, 1);
        let mut diagonal = vec![0.0; n];
        for i in 0..n - 1 {
            let w = 10f64.powf(-6.0 * i as f64 / (n - 2) as f64);
            diagonal[i] += w;
            diagonal[i + 1] += w;
            band.set_lower(i + 1, i, -w);
        }
        for (i, d) in diagonal.iter().enumerate() {
            band.set_diagonal(i, *d);
        }
        let shift = cholesky_safety_shift(&band);
        let pairs = lowest_eigenpairs_banded(&band, 1, &vec![1.0; n], shift, 1e-12, 1000).unwrap();
        let expected = 1.0 / (n as f64).sqrt();
        for v in &pairs[0].1 {
            assert!((v - expected).abs() < 1e-6 * expected, "{v} vs {expected}");
        }
        assert!(pairs[0].0.abs() < 10.0 * shift, "{}", pairs[0].0);
    }
}
