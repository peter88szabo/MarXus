pub fn cenmass(m: &[usize], q: &mut [f64], p: &mut [f64], mass: &[f64]) {
    // Initialize variables
    let mut vcm = [0.0; 3];
    let mut qcm = [0.0; 3];
    let mut tot_mass = 0.0;
    let n = m.len();

    // Calculate total mass of fragment
    for i in 0..n {
        let j = m[i];
        tot_mass += mass[j];
    }

    // Calculate center of mass velocity and coordinates
    for i in 0..n {
        let j = m[i];
        let jx = 3 * j;
        let jy = jx + 1;
        let jz = jx + 2;

        vcm[0] += p[jx];
        vcm[1] += p[jy];
        vcm[2] += p[jz];

        qcm[0] += mass[j] * q[jx];
        qcm[1] += mass[j] * q[jy];
        qcm[2] += mass[j] * q[jz];
    }

    for i in 0..3 {
        vcm[i] /= tot_mass;
        qcm[i] /= tot_mass;
    }

    // Update momenta and coordinates
    for i in 0..n {
        let j = m[i];
        let jx = 3 * j;
        let jy = jx + 1;
        let jz = jx + 2;

        let pcm_x = mass[j] * vcm[0];
        let pcm_y = mass[j] * vcm[1];
        let pcm_z = mass[j] * vcm[2];

        p[jx] -= pcm_x;
        p[jy] -= pcm_y;
        p[jz] -= pcm_z;

        q[jx] -= qcm[0];
        q[jy] -= qcm[1];
        q[jz] -= qcm[2];
    }
}

