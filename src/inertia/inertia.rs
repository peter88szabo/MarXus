use crate::constants::{INERTIA_AV, INERTIA_PH, INERTIA_SL, PI};
use crate::numeric::jacobi_diag::jacobi;

const TOCM1: f64 = 1.0e2 * (INERTIA_PH * INERTIA_AV) / (8.0 * PI * PI * INERTIA_SL);
const TOMHZ: f64 = 1.0e6 * (INERTIA_PH * INERTIA_AV) / (8.0 * PI * PI);

pub fn get_brot(xyz: &Vec<[f64; 3]>, mass: &Vec<f64>) -> [f64; 3] {
    // Centre of mass: the inertia tensor must be taken about the centre of mass,
    // otherwise the parallel-axis term M*|r_cm|^2 inflates the moments.
    let total_mass: f64 = mass.iter().sum();
    let mut cm = [0.0; 3];
    for (k, &mass) in mass.iter().enumerate() {
        for i in 0..3 {
            cm[i] += mass * xyz[k][i];
        }
    }
    for i in 0..3 {
        cm[i] /= total_mass;
    }

    let mut ixx = 0.0;
    let mut iyy = 0.0;
    let mut izz = 0.0;
    let mut ixy = 0.0;
    let mut ixz = 0.0;
    let mut iyz = 0.0;

    // Inertia tensor I = sum_k m_k (|r_k|^2 * 1 - r_k r_k^T), r_k relative to the centre of mass;
    // the off-diagonal elements carry a minus sign.
    for (k, &mass) in mass.iter().enumerate() {
        let x = xyz[k][0] - cm[0];
        let y = xyz[k][1] - cm[1];
        let z = xyz[k][2] - cm[2];

        let x2 = x * x;
        let y2 = y * y;
        let z2 = z * z;

        ixx += mass * (y2 + z2);
        iyy += mass * (x2 + z2);
        izz += mass * (x2 + y2);
        ixy -= mass * x * y;
        ixz -= mass * x * z;
        iyz -= mass * y * z;
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

#[cfg(test)]
mod tests {
    use super::*;

    // Rotation matrix R = Rz(gamma) * Ry(beta) * Rz(alpha) for generic Euler angles,
    // so that all three products of inertia become non-zero after rotation.
    fn euler_rotation(alpha: f64, beta: f64, gamma: f64) -> [[f64; 3]; 3] {
        let (sa, ca) = alpha.sin_cos();
        let (sb, cb) = beta.sin_cos();
        let (sg, cg) = gamma.sin_cos();
        let rz_a = [[ca, -sa, 0.0], [sa, ca, 0.0], [0.0, 0.0, 1.0]];
        let ry_b = [[cb, 0.0, sb], [0.0, 1.0, 0.0], [-sb, 0.0, cb]];
        let rz_g = [[cg, -sg, 0.0], [sg, cg, 0.0], [0.0, 0.0, 1.0]];
        let mul = |p: [[f64; 3]; 3], q: [[f64; 3]; 3]| {
            let mut r = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    for k in 0..3 {
                        r[i][j] += p[i][k] * q[k][j];
                    }
                }
            }
            r
        };
        mul(rz_g, mul(ry_b, rz_a))
    }

    #[test]
    fn brot_is_invariant_under_rotation_and_translation() {
        // Pairs of equal masses at +-a x, +-b y, +-c z: the centre of mass is the origin and the
        // inertia tensor is diagonal, so the principal moments are known exactly:
        //   I_x = 2(m2 b^2 + m3 c^2),  I_y = 2(m1 a^2 + m3 c^2),  I_z = 2(m1 a^2 + m2 b^2)
        let (m1, m2, m3) = (1.0, 12.0, 16.0);
        let (a, b, c) = (1.1, 0.7, 1.9);
        let mass = vec![m1, m1, m2, m2, m3, m3];
        let body = [
            [a, 0.0, 0.0],
            [-a, 0.0, 0.0],
            [0.0, b, 0.0],
            [0.0, -b, 0.0],
            [0.0, 0.0, c],
            [0.0, 0.0, -c],
        ];

        let mut moments = [
            2.0 * (m2 * b * b + m3 * c * c),
            2.0 * (m1 * a * a + m3 * c * c),
            2.0 * (m1 * a * a + m2 * b * b),
        ];
        moments.sort_by(|p, q| p.partial_cmp(q).unwrap());
        let expected = [TOCM1 / moments[0], TOCM1 / moments[1], TOCM1 / moments[2]];

        // Rotate to a generic orientation and move the centre of mass away from the origin.
        let rot = euler_rotation(0.3, 0.7, 1.1);
        let shift = [2.5, -1.3, 0.8];
        let xyz: Vec<[f64; 3]> = body
            .iter()
            .map(|r| {
                let mut out = [0.0; 3];
                for i in 0..3 {
                    out[i] = rot[i][0] * r[0] + rot[i][1] * r[1] + rot[i][2] * r[2] + shift[i];
                }
                out
            })
            .collect();

        let brot = get_brot(&xyz, &mass);

        for i in 0..3 {
            let rel = (brot[i] - expected[i]).abs() / expected[i];
            assert!(
                rel < 1.0e-8,
                "B[{i}] = {} cm-1, expected {} cm-1 (rel. err. {rel:e})",
                brot[i],
                expected[i]
            );
        }
    }
}
