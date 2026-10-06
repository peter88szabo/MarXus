//! Interface to LAPACK for dense symmetric eigenproblems.
//!
//! The Fortran LAPACK routines are called directly through their C symbols (`dsyevd_`), with
//! column-major storage and Fortran integer arguments; the library is the system OpenBLAS (libopenblas),
//! linked when the crate is built with the `openblas` feature (default). Without the feature every
//! routine returns an error; nothing is replaced silently by another algorithm.
//!
//! DSYEVD: all eigenvalues and eigenvectors of a real symmetric matrix by Householder reduction to
//! tridiagonal form and the divide-and-conquer method (Cuppen, Numer. Math. 36, 177 (1981); Gu, Eisenstat,
//! SIAM J. Matrix Anal. Appl. 16, 172 (1995); E. Anderson et al., LAPACK Users' Guide, 3rd ed. (SIAM,
//! 1999), Sec. 2.4.4). Same result as `symmetric_eigen::symmetric_eigen` (Householder + QL), with blocked
//! (level-3 BLAS) reduction and a faster tridiagonal stage for large matrices.

#[cfg(feature = "openblas")]
use std::os::raw::{c_char, c_int};

#[cfg(feature = "openblas")]
#[link(name = "openblas")]
extern "C" {
    fn dsyevd_(
        jobz: *const c_char,
        uplo: *const c_char,
        n: *const c_int,
        a: *mut f64,
        lda: *const c_int,
        w: *mut f64,
        work: *mut f64,
        lwork: *const c_int,
        iwork: *mut c_int,
        liwork: *const c_int,
        info: *mut c_int,
    );
}

#[cfg(feature = "openblas")]
#[link(name = "openblas")]
extern "C" {
    fn openblas_set_num_threads(num_threads: c_int);
    fn openblas_get_num_threads() -> c_int;
}

/// Number of threads that OpenBLAS uses inside one call (process-wide), e.g. 1 when the conditions of a
/// run are already computed in parallel. Returns the number now in use; 1 without the `openblas` feature.
pub fn set_blas_threads(threads: usize) -> usize {
    #[cfg(feature = "openblas")]
    {
        // SAFETY: OpenBLAS functions without pointer arguments; the thread count is process-wide state of
        // the library.
        unsafe {
            openblas_set_num_threads(threads.max(1) as c_int);
            openblas_get_num_threads().max(1) as usize
        }
    }
    #[cfg(not(feature = "openblas"))]
    {
        let _ = threads;
        1
    }
}

/// Eigenvalues (ascending) and orthonormal eigenvectors (vectors[k] belongs to values[k]) of the
/// symmetric matrix given by its rows, by LAPACK DSYEVD. Only the lower triangle is used, as in
/// `symmetric_eigen::symmetric_eigen`.
pub fn symmetric_eigen_lapack(matrix: &[Vec<f64>]) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
    let n = matrix.len();
    if matrix.iter().any(|row| row.len() != n) {
        return Err("symmetric_eigen_lapack: the matrix is not square.".into());
    }
    dsyevd_all_eigenpairs(matrix)
}

#[cfg(not(feature = "openblas"))]
fn dsyevd_all_eigenpairs(_matrix: &[Vec<f64>]) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
    Err("symmetric_eigen_lapack: MarXus was built without LAPACK; build with the `openblas` feature.".into())
}

#[cfg(feature = "openblas")]
fn dsyevd_all_eigenpairs(matrix: &[Vec<f64>]) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
    let n = matrix.len();
    if n == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    if n > c_int::MAX as usize {
        return Err("symmetric_eigen_lapack: the dimension exceeds the LAPACK integer range.".into());
    }
    // Row-major rows flattened are the column-major transpose; with UPLO = 'U' LAPACK reads the upper
    // triangle of that transpose, i.e. the lower triangle of the rows.
    let mut a: Vec<f64> = matrix.iter().flat_map(|row| row.iter().copied()).collect();
    let mut eigenvalues = vec![0.0; n];
    let order = n as c_int;
    let jobz = b'V' as c_char;
    let uplo = b'U' as c_char;
    let mut info: c_int = 0;
    // Workspace query (LWORK = LIWORK = -1); the optimal sizes enable the blocked reduction.
    let mut work_query = [0.0];
    let mut integer_query: [c_int; 1] = [0];
    let query: c_int = -1;
    unsafe {
        dsyevd_(
            &jobz,
            &uplo,
            &order,
            a.as_mut_ptr(),
            &order,
            eigenvalues.as_mut_ptr(),
            work_query.as_mut_ptr(),
            &query,
            integer_query.as_mut_ptr(),
            &query,
            &mut info,
        );
    }
    if info != 0 || !(work_query[0] >= 1.0 && work_query[0] <= c_int::MAX as f64) || integer_query[0] < 1 {
        return Err(format!("symmetric_eigen_lapack: invalid DSYEVD workspace query (info = {info})."));
    }
    let lwork = work_query[0].ceil() as c_int;
    let liwork = integer_query[0];
    let mut work = vec![0.0; lwork as usize];
    let mut integer_work: Vec<c_int> = vec![0; liwork as usize];
    unsafe {
        dsyevd_(
            &jobz,
            &uplo,
            &order,
            a.as_mut_ptr(),
            &order,
            eigenvalues.as_mut_ptr(),
            work.as_mut_ptr(),
            &lwork,
            integer_work.as_mut_ptr(),
            &liwork,
            &mut info,
        );
    }
    if info < 0 {
        return Err(format!("symmetric_eigen_lapack: DSYEVD rejected argument {}.", -info));
    }
    if info > 0 {
        return Err(format!("symmetric_eigen_lapack: DSYEVD did not converge (info = {info})."));
    }
    // Column-major eigenvectors: column k, contiguous, belongs to eigenvalue k (ascending).
    let vectors = a.chunks(n).map(|column| column.to_vec()).collect();
    Ok((eigenvalues, vectors))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numeric::symmetric_eigen::symmetric_eigen;

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

    #[cfg(feature = "openblas")]
    #[test]
    fn the_blas_thread_count_can_be_set_and_lapack_still_agrees() {
        assert_eq!(set_blas_threads(1), 1);
        let a = random_symmetric(40, 7);
        let (lapack, _) = symmetric_eigen_lapack(&a).unwrap();
        let (reference, _) = symmetric_eigen(&a).unwrap();
        for (x, y) in lapack.iter().zip(&reference) {
            assert!((x - y).abs() < 1e-10, "{x} vs {y}");
        }
    }

    #[cfg(feature = "openblas")]
    #[test]
    fn divide_and_conquer_agrees_with_the_householder_ql_decomposition() {
        for (n, seed) in [(1, 11), (2, 12), (7, 13), (60, 14), (150, 15)] {
            let a = random_symmetric(n, seed);
            let (values, vectors) = symmetric_eigen_lapack(&a).unwrap();
            let (reference_values, reference_vectors) = symmetric_eigen(&a).unwrap();
            assert_eq!(values.len(), n);
            assert!(values.windows(2).all(|w| w[0] <= w[1]), "eigenvalues not ascending");
            for k in 0..n {
                assert!((values[k] - reference_values[k]).abs() < 1e-12, "{} vs {}", values[k], reference_values[k]);
                // A v_k = lambda_k v_k, |v_k| = 1, and v_k = +-(reference vector) (simple eigenvalues).
                let v = &vectors[k];
                for i in 0..n {
                    let av: f64 = (0..n).map(|j| a[i][j] * v[j]).sum();
                    assert!((av - values[k] * v[i]).abs() < 1e-12);
                }
                assert!((v.iter().map(|x| x * x).sum::<f64>() - 1.0).abs() < 1e-12);
                let overlap: f64 = v.iter().zip(&reference_vectors[k]).map(|(a, b)| a * b).sum();
                assert!((overlap.abs() - 1.0).abs() < 1e-9, "overlap {overlap}");
            }
        }
    }

    #[cfg(feature = "openblas")]
    #[test]
    fn only_the_lower_triangle_is_read() {
        // The same convention as `symmetric_eigen`: the strict upper triangle is ignored.
        let mut a = random_symmetric(9, 21);
        let (values, _) = symmetric_eigen_lapack(&a).unwrap();
        for i in 0..9 {
            for j in i + 1..9 {
                a[i][j] = 1.0e3;
            }
        }
        let (garbled, _) = symmetric_eigen_lapack(&a).unwrap();
        assert_eq!(values, garbled);
    }

    #[test]
    fn a_non_square_matrix_is_an_error() {
        assert!(symmetric_eigen_lapack(&[vec![1.0, 0.0], vec![0.0]]).is_err());
    }

    #[cfg(not(feature = "openblas"))]
    #[test]
    fn without_the_openblas_feature_the_lapack_route_is_an_error() {
        let error = symmetric_eigen_lapack(&random_symmetric(3, 1)).unwrap_err();
        assert!(error.contains("openblas"), "{error}");
    }
}
