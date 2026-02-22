/// Dense matrix (row-major) with minimal operations we need.
#[derive(Clone)]
pub struct DenseMatrix {
    n: usize,
    data: Vec<f64>,
}

impl DenseMatrix {
    pub(crate) fn zeros(n: usize) -> Self {
        Self {
            n,
            data: vec![0.0; n * n],
        }
    }
    pub(crate) fn size(&self) -> usize {
        self.n
    }

    pub(crate) fn get(&self, row: usize, col: usize) -> f64 {
        self.data[row * self.n + col]
    }
    pub(crate) fn set(&mut self, row: usize, col: usize, value: f64) {
        self.data[row * self.n + col] = value;
    }
    pub(crate) fn add(&mut self, row: usize, col: usize, delta: f64) {
        let idx = row * self.n + col;
        self.data[idx] += delta;
    }

    pub(crate) fn diagonal(&self) -> Result<Vec<f64>, String> {
        let n = self.n;
        let mut d = vec![0.0; n];
        for i in 0..n {
            d[i] = self.get(i, i);
        }
        Ok(d)
    }

    pub(crate) fn matvec(&self, x: &[f64]) -> Result<Vec<f64>, String> {
        let n = self.n;
        if x.len() != n {
            return Err("Dimension mismatch in DenseMatrix::matvec".into());
        }
        let mut y = vec![0.0; n];
        for row in 0..n {
            let mut sum = 0.0;
            let base = row * n;
            for col in 0..n {
                sum += self.data[base + col] * x[col];
            }
            y[row] = sum;
        }
        Ok(y)
    }

    pub(crate) fn frobenius_norm(&self) -> f64 {
        self.data
            .iter()
            .copied()
            .map(|x| x * x)
            .sum::<f64>()
            .sqrt()
    }

    pub(crate) fn symmetry_relative_frobenius(&self) -> f64 {
        let n = self.n;
        let norm_a = self.frobenius_norm().max(1e-300);
        let mut diff2 = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let d = self.get(i, j) - self.get(j, i);
                diff2 += 2.0 * d * d;
            }
        }
        diff2.sqrt() / norm_a
    }

    pub(crate) fn swap_rows_cols_symmetric(&mut self, i: usize, j: usize) {
        if i == j {
            return;
        }
        let n = self.n;

        // Swap columns.
        for row in 0..n {
            let idx_i = row * n + i;
            let idx_j = row * n + j;
            self.data.swap(idx_i, idx_j);
        }

        // Swap rows.
        let row_i = i * n;
        let row_j = j * n;
        for col in 0..n {
            self.data.swap(row_i + col, row_j + col);
        }
    }

    pub(crate) fn scaled(&self, factor: f64) -> Self {
        let mut out = self.clone();
        for v in &mut out.data {
            *v *= factor;
        }
        out
    }

    /// Similarity transform A = W * M * W^{-1} for diagonal W.
    pub(crate) fn similarity_transform(&self, scale: &DiagonalScale) -> Self {
        let n = self.n;
        let mut out = DenseMatrix::zeros(n);
        for row in 0..n {
            for col in 0..n {
                let m = self.get(row, col);
                let a = m * (scale.diagonal[row] / scale.diagonal[col]);
                out.set(row, col, a);
            }
        }
        out
    }
}

/// Diagonal scale W.
pub struct DiagonalScale {
    pub(crate) diagonal: Vec<f64>,
}

impl DiagonalScale {
    pub(crate) fn apply_to_vector(&self, vec_in: &[f64]) -> Vec<f64> {
        vec_in
            .iter()
            .zip(self.diagonal.iter())
            .map(|(x, w)| x * w)
            .collect()
    }
    pub(crate) fn inverse_apply_to_vector(&self, vec_in: &[f64]) -> Vec<f64> {
        vec_in
            .iter()
            .zip(self.diagonal.iter())
            .map(|(x, w)| x / w)
            .collect()
    }
}

pub struct CholeskyDiagnostics {
    pub min_pivot: f64,
}

/// Solve Ax = b for SPD A using Cholesky (dense).
pub(crate) fn cholesky_solve_spd(a: &DenseMatrix, b: &[f64]) -> Result<Vec<f64>, String> {
    cholesky_solve_spd_with_diagnostics(a, b).map(|(x, _)| x)
}

pub(crate) fn cholesky_solve_spd_with_diagnostics(
    a: &DenseMatrix,
    b: &[f64],
) -> Result<(Vec<f64>, CholeskyDiagnostics), String> {
    let n = a.size();
    if b.len() != n {
        return Err("Dimension mismatch in cholesky_solve_spd".into());
    }

    let mut l = DenseMatrix::zeros(n);
    let mut min_pivot = f64::INFINITY;

    for i in 0..n {
        for j in 0..=i {
            let mut sum = a.get(i, j);
            for k in 0..j {
                sum -= l.get(i, k) * l.get(j, k);
            }
            if i == j {
                if sum <= 0.0 {
                    return Err("Matrix not SPD (negative pivot)".into());
                }
                min_pivot = min_pivot.min(sum);
                l.set(i, j, sum.sqrt());
            } else {
                let diag = l.get(j, j);
                if diag == 0.0 {
                    return Err("Zero diagonal in Cholesky factor".into());
                }
                l.set(i, j, sum / diag);
            }
        }
    }

    // Forward solve: L y = b
    let mut y = vec![0.0; n];
    for i in 0..n {
        let mut sum = b[i];
        for k in 0..i {
            sum -= l.get(i, k) * y[k];
        }
        y[i] = sum / l.get(i, i);
    }

    // Back solve: L^T x = y
    let mut x = vec![0.0; n];
    for i_rev in 0..n {
        let i = n - 1 - i_rev;
        let mut sum = y[i];
        for k in (i + 1)..n {
            sum -= l.get(k, i) * x[k];
        }
        x[i] = sum / l.get(i, i);
    }

    Ok((
        x,
        CholeskyDiagnostics {
            min_pivot: min_pivot.max(0.0),
        },
    ))
}
