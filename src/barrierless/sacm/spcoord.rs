#[derive(Debug, Clone, Copy)]
pub struct AtomCoord {
    pub mass: f64,
    pub position: [f64; 3],
}

#[derive(Debug, Clone)]
pub struct FragmentSpec {
    pub atom_indices: Vec<usize>,
}

#[derive(Debug, Clone)]
pub struct SpCoordInput {
    pub atoms: Vec<AtomCoord>,
    pub fragment_x: FragmentSpec,
    pub fragment_y: FragmentSpec,
    pub fragment_z: Option<FragmentSpec>,
}

#[derive(Debug, Clone)]
pub struct SpCoordResult {
    pub center_x: [f64; 3],
    pub center_y: [f64; 3],
    pub center_z: Option<[f64; 3]>,
    pub distance_xy: f64,
    pub distance_zy: Option<f64>,
    pub angle_deg: Option<f64>,
    pub mass_x: f64,
    pub mass_y: f64,
    pub mass_z: Option<f64>,
}

/// Compute pseudo-fragment coordinates for diatomic or triatomic capture.
pub fn compute_spcoord(input: &SpCoordInput) -> SpCoordResult {
    let (center_x, mass_x) = center_of_mass(&input.atoms, &input.fragment_x.atom_indices);
    let (center_y, mass_y) = center_of_mass(&input.atoms, &input.fragment_y.atom_indices);

    let vec_xy = [
        center_x[0] - center_y[0],
        center_x[1] - center_y[1],
        center_x[2] - center_y[2],
    ];
    let distance_xy = norm(vec_xy);

    if let Some(ref frag_z) = input.fragment_z {
        let (center_z, mass_z) = center_of_mass(&input.atoms, &frag_z.atom_indices);
        let vec_zy = [
            center_z[0] - center_y[0],
            center_z[1] - center_y[1],
            center_z[2] - center_y[2],
        ];
        let distance_zy = norm(vec_zy);
        let cos_angle = if distance_xy > 0.0 && distance_zy > 0.0 {
            dot(vec_xy, vec_zy) / (distance_xy * distance_zy)
        } else {
            1.0
        };
        let angle_deg = cos_angle.clamp(-1.0, 1.0).acos() * 180.0 / std::f64::consts::PI;

        SpCoordResult {
            center_x,
            center_y,
            center_z: Some(center_z),
            distance_xy,
            distance_zy: Some(distance_zy),
            angle_deg: Some(angle_deg),
            mass_x,
            mass_y,
            mass_z: Some(mass_z),
        }
    } else {
        SpCoordResult {
            center_x,
            center_y,
            center_z: None,
            distance_xy,
            distance_zy: None,
            angle_deg: None,
            mass_x,
            mass_y,
            mass_z: None,
        }
    }
}

/// Center of mass for a fragment defined by atom indices.
fn center_of_mass(atoms: &[AtomCoord], indices: &[usize]) -> ([f64; 3], f64) {
    let mut total_mass = 0.0;
    let mut acc = [0.0; 3];

    for &idx in indices {
        let atom = atoms.get(idx).expect("atom index out of range");
        total_mass += atom.mass;
        acc[0] += atom.mass * atom.position[0];
        acc[1] += atom.mass * atom.position[1];
        acc[2] += atom.mass * atom.position[2];
    }

    if total_mass == 0.0 {
        return (acc, total_mass);
    }

    (
        [acc[0] / total_mass, acc[1] / total_mass, acc[2] / total_mass],
        total_mass,
    )
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(v: [f64; 3]) -> f64 {
    dot(v, v).sqrt()
}
