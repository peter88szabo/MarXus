pub fn scalar_x_vector(scalar: f64, vec: &[f64]) -> Vec<f64> {
    vec.iter_mut().for_each(|element| *element *= scalar  )
}

pub fn scalar_x_matrix(mat: &[Vec<f64>], scalar: f64) -> Vec<Vec<f64>> {
    mat.iter().map(|row| row.iter().map(|&x| x * scalar).collect()).collect()
}

pub fn vector_x_vector_elementwise(vec1: &[f64], vec2: &[f64]) -> Vec<f64> {
    assert_eq!(vec1.len(), vec2.len());
    vec1.iter().zip(vec2.iter()).map(|(&x, &y)| x * y).collect()
}

pub fn dot_product(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn diad(a: &[f64], b: &[f64], c: &mut Vec<Vec<f64>>) {
    // Calculate the outer/diadic/tensor product of vectors a0 and b0 and store it in matrix c
    let alen = a.len();
    let blen = b.len();

    if alen != blen {
        panic!("Two vectors in diad() function must have the same size");
    } 

    for i in 0..alen {
        for j in 0..alen {
            c[i][j] = a[i] * b[j];
        }
    }
}



//A*B*C
//A^T*B*C
//A*B*C^T
//<a|B|C>
//exponential of matrix
//transpose of a matrix
//inverse of a matrix
//A*B^-1
//A^-1*B

fn matrix_x_vector(vector: &[f64], matrix: &[Vec<f64>]) -> Option<Vec<f64>> {
    let vector_len = vector.len();
    let matrix_rows = matrix.len();
    let matrix_cols = matrix[0].len();

    if vector_len != matrix_rows {
        return None;
    }

    let mut result = vec![0.0; matrix_cols];

    for i in 0..matrix_cols {
        for j in 0..vector_len {
            result[i] += vector[j] * matrix[j][i];
        }
    }

    Some(result)
}



pub fn matrix_add(sgn: i8, mat_a: &[Vec<f64>], mat_b: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let rows1 = mat_a.len();
    let cols1 = mat_a[0].len();
    let rows2 = mat_b.len();
    let cols2 = mat_b[0].len();

    if sgn != 1 || sgn != -1 {
        panic!("Wrong mode for sign in matrix_add function");
    }

    if rows1 != rows2 || cols1 != cols2 {
        return None;
    }

    let mut result = vec![vec![0.0; cols1]; rows1];

    for i in 0..rows1 {
        for j in 0..cols1 {
            result[i][j] = mat_a[i][j] + sgn * mat_b[i][j];
        }
    }

    Some(result)
}



pub fn linearcomb_two_matrix(sc1: f64, sc2: f64, mat_a: &[Vec<f64>], mat_b: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let rows1 = mat_a.len();
    let cols1 = mat_a[0].len();
    let rows2 = mat_b.len();
    let cols2 = mat_b[0].len();

    if sgn != 1 || sgn != -1 {
        panic!("Wrong mode for sign in matrix_add function");
    }

    if rows1 != rows2 || cols1 != cols2 {
        return None;
    }

    let mut result = vec![vec![0.0; cols1]; rows1];

    for i in 0..rows1 {
        for j in 0..cols1 {
            result[i][j] = (sc1 * mat_a[i][j]) + (sc2 * mat_b[i][j]);
        }
    }

    Some(result)
}


pub fn matmul_elementwise(scalar: f64, mat_a: &[Vec<f64>], mat_b: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let rows1 = mat_a.len();
    let cols1 = mat_a[0].len();
    let rows2 = mat_b.len();
    let cols2 = mat_b[0].len();

    if rows1 != rows2 || cols1 != cols2 {
        return None;
    }

    let mut result = vec![vec![0.0; cols1]; rows1];

    for i in 0..rows1 {
        for j in 0..cols1 {
            result[i][j] = scalar * mat_a[i][j] * mat_b[i][j];
        }
    }

    Some(result)
}

pub fn matdiv_elementwise(scalar: f64, mat_a: &[Vec<f64>], mat_b: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let rows1 = mat_a.len();
    let cols1 = mat_a[0].len();
    let rows2 = mat_b.len();
    let cols2 = mat_b[0].len();

    if rows1 != rows2 || cols1 != cols2 {
        return None;
    }

    let mut result = vec![vec![0.0; cols1]; rows1];

    for i in 0..rows1 {
        for j in 0..cols1 {
            result[i][j] = scalar * mat_a[i][j] / mat_b[i][j];
        }
    }

    Some(result)
}


pub fn matmul(mat_a: &[Vec<f64>], mat_b: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let rows1 = mat_a.len();
    let cols1 = mat_a[0].len();
    let rows2 = mat_b.len();
    let cols2 = mat_b[0].len();

    // Check if the matrices can be multiplied
    if cols1 != rows2 {
        return None;
    }

    let mut res = vec![vec![0.0; cols2]; rows1];

    for i in 0..rows1 {
        for j in 0..cols2 {
            for k in 0..cols1 {
                res[i][j] += mat_a[i][k] * mat_b[k][j];
            }
        }
    }

    Some(res)
}

