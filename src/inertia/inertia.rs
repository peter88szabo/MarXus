use crate::constants::{INERTIA_AV, INERTIA_PH, INERTIA_SL, PI};
use crate::numeric::jacobi_diag::jacobi;

const TOCM1: f64 = 1.0e2 * (INERTIA_PH * INERTIA_AV) / (8.0 * PI * PI * INERTIA_SL);
const TOMHZ: f64 = 1.0e6 * (INERTIA_PH * INERTIA_AV) / (8.0 * PI * PI);

pub fn get_brot(xyz: &Vec<[f64; 3]>, mass: &Vec<f64>) -> [f64; 3] {
    let mut ixx = 0.0;
    let mut iyy = 0.0;
    let mut izz = 0.0;
    let mut ixy = 0.0;
    let mut ixz = 0.0;
    let mut iyz = 0.0;

    // Calculate elements of the inertia tensor
    for (k, &mass) in mass.iter().enumerate() {
        let x = xyz[k][0];
        let y = xyz[k][1];
        let z = xyz[k][2];

        let x2 = x * x;
        let y2 = y * y;
        let z2 = z * z;

        ixx += mass * (y2 + z2);
        iyy += mass * (x2 + z2);
        izz += mass * (x2 + y2);
        ixy += mass * x * y;
        ixz += mass * x * z;
        iyz += mass * y * z;
    }

    let mut inertia_tensor: Vec<Vec<f64>> = vec![
        vec![ixx, ixy, ixz],
        vec![ixy, iyy, iyz],
        vec![ixz, iyz, izz],
    ];

    let (_eigvec, eigval) = jacobi(&mut inertia_tensor, 100, 1.0e-10, 1);

    let mut inertia_amuang2 = [0.0; 3]; // Initialize with a size of 3
    let mut brot_cm1 = [0.0; 3]; // Initialize with a size of 3
    let mut brot_mhz = [0.0; 3]; // Initialize with a size of 3

    for i in 0..3 {
        brot_cm1[i] = TOCM1 / eigval[i];
        brot_mhz[i] = TOMHZ / eigval[i];
        inertia_amuang2[i] = eigval[i] / INERTIA_AV * 10.0;
    }

    //if you need brot or intertia in other units, then just print them
    return brot_cm1;
}
