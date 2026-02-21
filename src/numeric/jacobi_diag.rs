pub fn jacobi(
    a: &mut Vec<Vec<f64>>,
    itmax: usize,
    th: f64,
    iord: i32,
) -> (Vec<Vec<f64>>, Vec<f64>) {
    //================================================================================================
    // a(:,:)--> input matrix to be diagonalized. (It's restored at the end of the procedure)
    // th    --> treshold for the iteration
    // itmax --> maximum number of the jacobian rotations
    // iord  --> sorting variable, +1: increasing order, -1: decrecasing order
    // eigvec(:,:) --> matrix of the eigenvectors
    // eigval(:)   --> array of the eigenvalues
    //================================================================================================
    let n = a.len();
    let mut eigvec: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
    let mut eigval: Vec<f64> = vec![0.0; n];

    let mut xj: f64;
    let mut amax: f64;
    let mut aii: f64;
    let mut ajj: f64;
    let mut aij: f64;
    let mut alpha: f64;
    let mut t: f64;
    let mut c: f64;
    let mut s: f64;

    let mut iter = 0;

    // Initialize x to the identity matrix
    for i in 0..n {
        for j in 0..n {
            eigvec[i][j] = if i == j { 1.0 } else { 0.0 };
        }
    }

    // Initialize d with diagonal elements of a
    for i in 0..n {
        eigval[i] = a[i][i];
    }

    // Jacobi iteration
    while iter < itmax {
        amax = 0.0;

        for i in 1..n {
            for j in 0..i {
                aii = eigval[i];
                ajj = eigval[j];
                aij = a[i][j];

                if aij.abs() > amax {
                    amax = aij.abs();
                }

                if aij.abs() <= th {
                    continue;
                }

                alpha = 0.5 * (aii - ajj) / aij;

                t = if alpha < 0.0 {
                    -alpha - (1.0 + alpha * alpha).sqrt()
                } else {
                    -alpha + (1.0 + alpha * alpha).sqrt()
                };

                c = 1.0 / (1.0 + t * t).sqrt();
                s = c * t;

                for k in 0..n {
                    xj = c * eigvec[k][j] - s * eigvec[k][i];
                    eigvec[k][i] = s * eigvec[k][j] + c * eigvec[k][i];
                    eigvec[k][j] = xj;

                    if k == j {
                        continue;
                    }

                    if k < j {
                        xj = c * a[j][k] - s * a[i][k];
                        a[i][k] = s * a[j][k] + c * a[i][k];
                        a[j][k] = xj;
                        continue;
                    }

                    if k == i {
                        continue;
                    }

                    if k < i {
                        xj = c * a[k][j] - s * a[i][k];
                        a[i][k] = s * a[k][j] + c * a[i][k];
                        a[k][j] = xj;
                        continue;
                    }

                    xj = c * a[k][j] - s * a[k][i];
                    a[k][i] = s * a[k][j] + c * a[k][i];
                    a[k][j] = xj;
                }

                eigval[i] = c * c * aii + s * s * ajj + 2.0 * s * c * aij;
                eigval[j] = c * c * ajj + s * s * aii - 2.0 * s * c * aij;
                a[i][j] = 0.0;
            }
        }

        if amax <= th {
            break;
        }
        iter += 1;
    }

    if iter == itmax {
        panic!("Maximum iteration reached in Jacobi diagonalizetion. Increase itmax!");
    }

    if iord == 1 {
        // Arrange eigenvalues in increasing order
        for k in 0..n - 1 {
            let mut dmn = eigval[k];
            let mut kmin = k;

            for j in k + 1..n {
                if dmn > eigval[j] {
                    kmin = j;
                    dmn = eigval[j];
                }
            }

            if k != kmin {
                for j in 0..n {
                    let temp = eigvec[j][kmin];
                    eigvec[j][kmin] = eigvec[j][k];
                    eigvec[j][k] = temp;
                }
                let temp = eigval[kmin];
                eigval[kmin] = eigval[k];
                eigval[k] = temp;
            }
        }
    } else if iord == -1 {
        // Arrange eigenvalues in decreasing order
        for k in 0..n - 1 {
            let mut dmx = eigval[k];
            let mut kmax = k;

            for j in k + 1..n {
                if dmx < eigval[j] {
                    kmax = j;
                    dmx = eigval[j];
                }
            }

            if k != kmax {
                for j in 0..n {
                    let temp = eigvec[j][kmax];
                    eigvec[j][kmax] = eigvec[j][k];
                    eigvec[j][k] = temp;
                }
                let temp = eigval[kmax];
                eigval[kmax] = eigval[k];
                eigval[k] = temp;
            }
        }
    }

    // Restore the original a(i,j) input matrix
    for i in 0..n {
        for j in 0..i {
            a[i][j] = a[j][i];
        }
    }

    println!("Iterative diagonalization is done in {} steps.", iter);

    return (eigvec, eigval);
}
//==================================================================================

pub fn print_eigval_and_eigvec(eigval: &[f64], eigvec: &[Vec<f64>]) {
    println!("Eigenvalues and Eigenvectors");
    println!();
    let (minval_1d, maxval_1d) = min_max_abs_value_1d_array(&eigval);
    let (minval_2d, maxval_2d) = min_max_abs_value_2d_array(&eigvec);

    // Print eigenvalues in a row
    for val in eigval {
        if minval_1d < 1e-4 || maxval_1d > 1e5 {
            print!("{:11.5e}    ", val);
        } else {
            print!("{:>11.5}    ", val);
        }
    }
    println!(); // Move to the next line

    // Empty line
    println!();

    // Print eigenvectors
    for row in eigvec {
        for &element in row {
            if minval_2d < 1e-4 || maxval_2d > 1e5 {
                print!("{:11.5e}    ", element);
            } else {
                print!("{:>11.5}    ", element);
            }
        }
        println!(); // Move to the next row
    }
}

pub fn print_matrix(matrix: &[Vec<f64>]) {
    for row in matrix {
        for &element in row {
            // Adjust the formatting as needed, e.g., {:.2} for two decimal places
            print!("{:11.5e}  ", element);
        }
        println!(); // Move to the next row
    }
}

fn min_max_abs_value_2d_array(matrix: &[Vec<f64>]) -> (f64, f64) {
    // Initialize the minimum and maximum absolute values
    let mut min_abs_val = f64::INFINITY;
    let mut max_abs_val = f64::NEG_INFINITY;

    // Iterate through the matrix to find the minimum and maximum absolute values
    for row in matrix.iter() {
        for &value in row.iter() {
            let abs_value = value.abs();
            if abs_value < min_abs_val {
                min_abs_val = abs_value;
            }
            if abs_value > max_abs_val {
                max_abs_val = abs_value;
            }
        }
    }
    return (min_abs_val, max_abs_val);
}

fn min_max_abs_value_1d_array(vector: &[f64]) -> (f64, f64) {
    let (max_abs, min_abs) =
        vector
            .iter()
            .fold((f64::NEG_INFINITY, f64::INFINITY), |(max, min), &x| {
                let abs_x = x.abs();
                (max.max(abs_x), min.min(abs_x))
            });

    return (min_abs, max_abs);
}
