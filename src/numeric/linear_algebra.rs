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

/// Solve Ax = b for SPD A using Cholesky (dense).
pub(crate) fn cholesky_solve_spd(a: &DenseMatrix, b: &[f64]) -> Result<Vec<f64>, String> {
    let n = a.size();
    if b.len() != n {
        return Err("Dimension mismatch in cholesky_solve_spd".into());
    }

    let mut l = DenseMatrix::zeros(n);

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

    Ok(x)
}
