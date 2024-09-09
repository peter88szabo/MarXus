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

    let inertia_tensor: [[f64; 3]; 3] = [
        [ixx, ixy, ixz],
        [ixy, iyy, iyz],
        [ixz, iyz, izz],
    ];

    let eigvals = [0.0, 0.0, 0.0];

    //TODO: it must be converted to cm-1 before return

    return eigvals;
}

