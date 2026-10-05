//! Inverse of a small dense matrix by Gauss-Jordan elimination with partial pivoting (W. H. Press,
//! S. A. Teukolsky, W. T. Vetterling, B. P. Flannery, Numerical Recipes in Fortran, 2nd ed. (1992), Sec. 2.1).

/// Inverse of the square matrix given by its rows. An error if the matrix is not square or is singular
/// (a pivot of relative size below 1e-14 of the largest element).
pub fn invert_dense(matrix: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n = matrix.len();
    if matrix.iter().any(|row| row.len() != n) {
        return Err("invert_dense: the matrix is not square.".into());
    }
    let scale = matrix.iter().flatten().fold(0.0_f64, |m, x| m.max(x.abs()));
    // Augmented matrix [A | I].
    let mut a: Vec<Vec<f64>> = matrix
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut r = row.clone();
            r.extend((0..n).map(|j| if i == j { 1.0 } else { 0.0 }));
            r
        })
        .collect();
    for col in 0..n {
        // Partial pivoting: the largest element of the column at or below the diagonal.
        let pivot_row = (col..n)
            .max_by(|&x, &y| a[x][col].abs().partial_cmp(&a[y][col].abs()).unwrap())
            .unwrap();
        if !(a[pivot_row][col].abs() > 1e-14 * scale) {
            return Err(format!("invert_dense: the matrix is singular (column {col})."));
        }
        a.swap(col, pivot_row);
        let pivot = a[col][col];
        a[col].iter_mut().for_each(|x| *x /= pivot);
        let pivot_row_values = a[col].clone();
        for (r, row) in a.iter_mut().enumerate() {
            if r != col {
                let factor = row[col];
                if factor != 0.0 {
                    row.iter_mut().zip(&pivot_row_values).for_each(|(x, p)| *x -= factor * p);
                }
            }
        }
    }
    Ok(a.into_iter().map(|row| row[n..].to_vec()).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pseudo_random(seed: &mut u64) -> f64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64 / (1u64 << 53) as f64) - 0.5
    }

    #[test]
    fn the_inverse_times_the_matrix_is_the_identity() {
        let mut seed = 3u64;
        for n in [1usize, 2, 4, 8] {
            let a: Vec<Vec<f64>> = (0..n).map(|_| (0..n).map(|_| pseudo_random(&mut seed)).collect()).collect();
            let inv = invert_dense(&a).unwrap();
            for i in 0..n {
                for j in 0..n {
                    let product: f64 = (0..n).map(|k| inv[i][k] * a[k][j]).sum();
                    assert!((product - if i == j { 1.0 } else { 0.0 }).abs() < 1e-12, "n {n}: ({i},{j}) {product}");
                }
            }
        }
    }

    #[test]
    fn a_zero_leading_element_is_handled_by_pivoting() {
        let inv = invert_dense(&[vec![0.0, 2.0], vec![3.0, 0.0]]).unwrap();
        assert_eq!(inv, vec![vec![0.0, 1.0 / 3.0], vec![0.5, 0.0]]);
    }

    #[test]
    fn singular_and_non_square_matrices_are_errors() {
        assert!(invert_dense(&[vec![1.0, 2.0], vec![2.0, 4.0]]).is_err());
        assert!(invert_dense(&[vec![1.0, 2.0], vec![2.0]]).is_err());
    }
}
