//! One-dimensional internal rotors (hindered and free) for the state counts of a species.
//!
//! The torsion is separated from the other degrees of freedom; its only coupling to the overall rotation is the
//! reduced moment of inertia (`reduced_moment_pitzer`: the Kilpatrick-Pitzer reduced moment I(1,3), K. S. Pitzer,
//! J. Chem. Phys. 14, 239 (1946); J. E. Kilpatrick, K. S. Pitzer, J. Chem. Phys. 17, 1064 (1949); or
//! `reduced_moment_bond_axis`: I_A I_B/(I_A + I_B) of the two groups about the bond axis, East, Radom, J. Chem.
//! Phys. 106, 6655 (1997), I(2,1)). The rotor levels are the eigenvalues of
//!   H = -B d^2/dphi^2 + V(phi),   B = h/(8 pi^2 c I_red),
//! in a Fourier basis of the functions of period 2 pi/sigma (sigma the rotor symmetry number): only the
//! sigma-periodic states are counted, which accounts for the symmetry without a division by sigma. The potential
//! V(phi) = c_0 + sum_k [a_k cos(k sigma phi) + b_k sin(k sigma phi)] is given by its Fourier coefficients or
//! interpolated through equidistant points on [0, 2 pi/sigma) (exact trigonometric interpolation). The levels are
//! counted from the rotor ground level E_0, and convolved with the cell counts of the other degrees of freedom:
//!   rho(E) = sum_n rho_rest(E - (E_n - E_0)),
//! the same for the sum of states. This is the uncoupled one-dimensional rotor model of the MESS deck format
//! (`Rotor Hindered`, `Rotor Free`).
//! All energies in cm-1, moments of inertia in amu Angstrom^2.

/// Fourier series of a torsional potential in the angle x = sigma phi:
/// V(x) = constant + sum_k [cosine[k-1] cos(k x) + sine[k-1] sin(k x)] (cm-1).
#[derive(Debug, Clone, PartialEq)]
pub struct TorsionalPotential {
    pub constant: f64,
    pub cosine: Vec<f64>,
    pub sine: Vec<f64>,
}

/// Default size of the Fourier basis (odd: the constant and K cosine-sine pairs, K = 499).
pub const DEFAULT_BASIS_SIZE: usize = 999;
/// Largest basis used when the energy range needs more levels.
pub const MAX_BASIS_SIZE: usize = 1999;

impl TorsionalPotential {
    /// V at x = sigma phi (cm-1).
    pub fn evaluate(&self, x: f64) -> f64 {
        self.constant
            + self.cosine.iter().enumerate().map(|(k, a)| a * ((k + 1) as f64 * x).cos()).sum::<f64>()
            + self.sine.iter().enumerate().map(|(k, b)| b * ((k + 1) as f64 * x).sin()).sum::<f64>()
    }

    /// Minimum and maximum of V on a fine grid of x (cm-1).
    pub fn extrema(&self) -> (f64, f64) {
        let n = 20_000;
        (0..n)
            .map(|i| self.evaluate(2.0 * std::f64::consts::PI * i as f64 / n as f64))
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| (lo.min(v), hi.max(v)))
    }
}

/// Trigonometric interpolation through N equidistant points V_j at x_j = 2 pi j/N (cm-1), the first at x = 0:
/// c_0 = mean, a_k = (2/N) sum_j V_j cos(k x_j), b_k = (2/N) sum_j V_j sin(k x_j) for k < N/2, and for even N the
/// last cosine coefficient a_(N/2) = (1/N) sum_j (-1)^j V_j (exact through every point).
pub fn potential_from_equidistant_points(values_cm1: &[f64]) -> Result<TorsionalPotential, String> {
    let n = values_cm1.len();
    if n == 0 || values_cm1.iter().any(|v| !v.is_finite()) {
        return Err("Torsional potential: no points, or a value that is not finite.".into());
    }
    let x = |j: usize| 2.0 * std::f64::consts::PI * j as f64 / n as f64;
    let constant = values_cm1.iter().sum::<f64>() / n as f64;
    let (mut cosine, mut sine) = (Vec::new(), Vec::new());
    for k in 1..=n / 2 {
        let a: f64 = values_cm1.iter().enumerate().map(|(j, v)| v * (k as f64 * x(j)).cos()).sum::<f64>();
        let b: f64 = values_cm1.iter().enumerate().map(|(j, v)| v * (k as f64 * x(j)).sin()).sum::<f64>();
        if 2 * k == n {
            cosine.push(a / n as f64);
        } else {
            cosine.push(2.0 * a / n as f64);
            sine.push(2.0 * b / n as f64);
        }
    }
    Ok(TorsionalPotential { constant, cosine, sine })
}

/// Fourier coefficients given in the order c_0, a_1, b_1, a_2, b_2, ... (cm-1).
pub fn potential_from_fourier_expansion(values_cm1: &[f64]) -> Result<TorsionalPotential, String> {
    let (&constant, rest) = values_cm1.split_first().ok_or("Torsional potential: no Fourier coefficients.")?;
    Ok(TorsionalPotential {
        constant,
        cosine: rest.iter().step_by(2).copied().collect(),
        sine: rest.iter().skip(1).step_by(2).copied().collect(),
    })
}

/// Eigenvalues (cm-1, ascending) of H = -B d^2/dphi^2 + V in the real Fourier basis of the functions of period
/// 2 pi/sigma: 1/sqrt(2 pi), cos(k x)/sqrt(pi), sin(k x)/sqrt(pi), k = 1..K, x = sigma phi, size 2K + 1.
/// Matrix elements: kinetic B sigma^2 k^2 on the diagonal; with A(0) = 2 c_0, A(m) = a_|m|, S(m) = sign(m) b_|m|,
///   <1|V|1> = c_0, <1|V|cos_k> = a_k/sqrt 2, <1|V|sin_k> = b_k/sqrt 2,
///   <cos_k|V|cos_l> = [A(k-l) + A(k+l)]/2, <sin_k|V|sin_l> = [A(k-l) - A(k+l)]/2, <cos_k|V|sin_l> = [S(l+k) + S(l-k)]/2.
pub fn rotor_levels(b_cm1: f64, sigma: u32, potential: &TorsionalPotential, basis_size: usize) -> Result<Vec<f64>, String> {
    if !(b_cm1 > 0.0) || sigma == 0 || basis_size % 2 == 0 {
        return Err(format!(
            "Internal rotor: B = {b_cm1} cm-1 must be positive, the symmetry number at least 1 and the basis size odd \
             ({basis_size})."
        ));
    }
    let k_max = (basis_size - 1) / 2;
    let a = |m: usize| if m == 0 { 2.0 * potential.constant } else { potential.cosine.get(m - 1).copied().unwrap_or(0.0) };
    let s = |m: isize| {
        let value = if m == 0 { 0.0 } else { potential.sine.get(m.unsigned_abs() - 1).copied().unwrap_or(0.0) };
        if m < 0 { -value } else { value }
    };
    // Index: 0 the constant; 2k - 1 cos(k x); 2k sin(k x).
    let n = basis_size;
    let mut h = vec![vec![0.0; n]; n];
    let kinetic = b_cm1 * (sigma * sigma) as f64;
    h[0][0] = potential.constant;
    for k in 1..=k_max {
        let (c, sn) = (2 * k - 1, 2 * k);
        h[0][c] = a(k) / std::f64::consts::SQRT_2;
        h[0][sn] = s(k as isize) / std::f64::consts::SQRT_2;
        for l in 1..=k_max {
            let (cl, sl) = (2 * l - 1, 2 * l);
            let diff = k.abs_diff(l);
            h[c][cl] = 0.5 * (a(diff) + a(k + l));
            h[sn][sl] = 0.5 * (a(diff) - a(k + l));
            h[c][sl] = 0.5 * (s((l + k) as isize) + s(l as isize - k as isize));
        }
        h[c][c] += kinetic * (k * k) as f64;
        h[sn][sn] += kinetic * (k * k) as f64;
    }
    for i in 0..n {
        for j in 0..i {
            let value = if h[j][i] != 0.0 { h[j][i] } else { h[i][j] };
            h[i][j] = value;
            h[j][i] = value;
        }
    }
    symmetric_eigenvalues(&h)
}

fn symmetric_eigenvalues(matrix: &[Vec<f64>]) -> Result<Vec<f64>, String> {
    #[cfg(feature = "openblas")]
    let (mut values, _) = crate::numeric::lapack_interface::symmetric_eigen_lapack(matrix)?;
    #[cfg(not(feature = "openblas"))]
    let (mut values, _) = crate::numeric::symmetric_eigen::symmetric_eigen(matrix)?;
    values.sort_by(|x, y| x.total_cmp(y));
    Ok(values)
}

/// A one-dimensional internal rotor: its ground level and its levels above the ground (cm-1).
#[derive(Debug, Clone, PartialEq)]
pub struct HinderedRotor {
    pub rotational_constant_cm1: f64,
    pub symmetry_number: u32,
    /// Lowest level, on the energy scale of the potential (cm-1).
    pub ground_energy_cm1: f64,
    /// E_n - E_0 of the levels kept: those below the top of the basis, B (K sigma)^2 + V_min (cm-1).
    pub levels_above_ground_cm1: Vec<f64>,
}

impl HinderedRotor {
    pub fn new(b_cm1: f64, sigma: u32, potential: TorsionalPotential, basis_size: usize) -> Result<Self, String> {
        let levels = rotor_levels(b_cm1, sigma, &potential, basis_size)?;
        let (v_min, _) = potential.extrema();
        let k_max = ((basis_size - 1) / 2) as f64;
        let top = b_cm1 * (k_max * sigma as f64).powi(2) + v_min;
        let ground = levels[0];
        Ok(Self {
            rotational_constant_cm1: b_cm1,
            symmetry_number: sigma,
            ground_energy_cm1: ground,
            levels_above_ground_cm1: levels.iter().filter(|&&e| e <= top).map(|e| e - ground).collect(),
        })
    }

    /// The rotor with levels up to `energy_above_ground_cm1`: the basis of `size_min` functions, enlarged to
    /// 2K + 1 with K = ceil(sqrt((E + E_0 - V_min)/B)/sigma) when that is larger, so that the top of the basis
    /// B (K sigma)^2 + V_min lies at or above E + E_0 (the rule of the MESS deck format, `HamiltonSizeMin`,
    /// `HamiltonSizeMax`). Above the barrier the levels lie near B (k sigma)^2 + c_0, so the highest level kept can
    /// still end below E; K then grows by a tenth until the levels reach E. A size above `size_max` is an error, not
    /// a truncation of the levels.
    pub fn for_energy_range(
        b_cm1: f64,
        sigma: u32,
        potential: TorsionalPotential,
        size_min: usize,
        size_max: usize,
        energy_above_ground_cm1: f64,
    ) -> Result<Self, String> {
        let mut rotor = Self::new(b_cm1, sigma, potential.clone(), size_min)?;
        let (v_min, _) = potential.extrema();
        let span = energy_above_ground_cm1 + rotor.ground_energy_cm1 - v_min;
        let mut size = if span > 0.0 { 2 * ((span / b_cm1).sqrt() / sigma as f64).ceil() as usize + 1 } else { 0 };
        size = size.max(size_min);
        loop {
            if size > size_max {
                return Err(format!(
                    "Internal rotor (B = {b_cm1} cm-1, symmetry {sigma}): the levels up to {energy_above_ground_cm1:.0} \
                     cm-1 above the ground need a Fourier basis of {size} functions or more, above HamiltonSizeMax = \
                     {size_max}. Raise HamiltonSizeMax in the Rotor block."
                ));
            }
            if size > size_min {
                rotor = Self::new(b_cm1, sigma, potential.clone(), size)?;
            }
            if rotor.highest_level_above_ground_cm1() >= energy_above_ground_cm1 {
                return Ok(rotor);
            }
            let k = (size - 1) / 2;
            size = 2 * (k + (k / 10).max(1)) + 1;
        }
    }

    /// Energy above the ground up to which the levels are complete (cm-1).
    pub fn highest_level_above_ground_cm1(&self) -> f64 {
        self.levels_above_ground_cm1.last().copied().unwrap_or(0.0)
    }
}

/// B = h/(8 pi^2 c I) in cm-1 for I in amu Angstrom^2 (the conversion of `inertia::get_brot`).
pub fn rotational_constant_cm1(reduced_moment_amu_angstrom2: f64) -> f64 {
    crate::inertia::inertia::ROTATIONAL_CONSTANT_CM1_AMU_ANGSTROM2 / reduced_moment_amu_angstrom2
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn check_rotor_atoms(n: usize, group: &[usize], axis: (usize, usize)) -> Result<(), String> {
    if group.is_empty() || axis.0 == axis.1 || axis.0 >= n || axis.1 >= n || group.iter().any(|&a| a >= n) {
        return Err(format!("Internal rotor: group {group:?} or axis {axis:?} does not fit the {n} atoms."));
    }
    if group.contains(&axis.0) || group.contains(&axis.1) {
        return Err(format!("Internal rotor: the group {group:?} must not contain an axis atom {axis:?}."));
    }
    Ok(())
}

/// Kilpatrick-Pitzer reduced moment I(1,3) (amu Angstrom^2) of the group of atoms rotating about the axis through
/// atoms `axis.0` -> `axis.1` (0-based): the motion v0_a = e x (r_a - r_A1) of the group atoms, corrected for the
/// overall rotation and translation it induces (zero total linear and angular momentum),
///   I_red = sum_(a in G) m_a |v0_a|^2 - L^T I^-1 L - P^2/M,   P = sum m v0,   L = sum m (r - R_cm) x v0,
/// with I the inertia tensor of the whole molecule about its centre of mass.
pub fn reduced_moment_pitzer(masses_amu: &[f64], coords_angstrom: &[[f64; 3]], group: &[usize], axis: (usize, usize)) -> Result<f64, String> {
    let n = masses_amu.len();
    check_rotor_atoms(n, group, axis)?;
    let total: f64 = masses_amu.iter().sum();
    let mut com = [0.0; 3];
    for (m, r) in masses_amu.iter().zip(coords_angstrom) {
        (0..3).for_each(|i| com[i] += m * r[i] / total);
    }
    let rel = |r: [f64; 3]| [r[0] - com[0], r[1] - com[1], r[2] - com[2]];
    let (r1, r2) = (coords_angstrom[axis.0], coords_angstrom[axis.1]);
    let d = [r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]];
    let norm = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    let e = [d[0] / norm, d[1] / norm, d[2] / norm];
    let (mut kinetic, mut p, mut l) = (0.0, [0.0; 3], [0.0; 3]);
    for &a in group {
        let m = masses_amu[a];
        let r = coords_angstrom[a];
        let v0 = cross(e, [r[0] - r1[0], r[1] - r1[1], r[2] - r1[2]]);
        kinetic += m * (v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2]);
        let c = cross(rel(r), v0);
        (0..3).for_each(|i| {
            p[i] += m * v0[i];
            l[i] += m * c[i];
        });
    }
    let mut inertia = [[0.0; 3]; 3];
    for (m, r) in masses_amu.iter().zip(coords_angstrom) {
        let x = rel(*r);
        let r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
        for i in 0..3 {
            for j in 0..3 {
                inertia[i][j] += m * (if i == j { r2 } else { 0.0 } - x[i] * x[j]);
            }
        }
    }
    let rows: Vec<Vec<f64>> = inertia.iter().map(|row| row.to_vec()).collect();
    let inverse = crate::numeric::dense_inverse::invert_dense(&rows).map_err(|e| format!("Internal rotor: inertia tensor: {e}"))?;
    let l_i_l: f64 = (0..3).map(|i| (0..3).map(|j| l[i] * inverse[i][j] * l[j]).sum::<f64>()).sum();
    let p2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
    let reduced = kinetic - l_i_l - p2 / total;
    if !(reduced > 0.0) {
        return Err(format!("Internal rotor: non-positive reduced moment {reduced} for group {group:?}."));
    }
    Ok(reduced)
}

/// Reduced moment I_A I_B/(I_A + I_B) (amu Angstrom^2) of the rotating group (A) and the rest of the molecule (B),
/// both about the axis line through atoms `axis.0` and `axis.1` (0-based).
pub fn reduced_moment_bond_axis(masses_amu: &[f64], coords_angstrom: &[[f64; 3]], group: &[usize], axis: (usize, usize)) -> Result<f64, String> {
    check_rotor_atoms(masses_amu.len(), group, axis)?;
    let (r1, r2) = (coords_angstrom[axis.0], coords_angstrom[axis.1]);
    let d = [r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]];
    let norm = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    let e = [d[0] / norm, d[1] / norm, d[2] / norm];
    let moment = |a: usize| {
        let r = coords_angstrom[a];
        let c = cross(e, [r[0] - r1[0], r[1] - r1[1], r[2] - r1[2]]);
        masses_amu[a] * (c[0] * c[0] + c[1] * c[1] + c[2] * c[2])
    };
    let i_a: f64 = group.iter().map(|&a| moment(a)).sum();
    let i_b: f64 = (0..masses_amu.len()).filter(|a| !group.contains(a)).map(moment).sum();
    if !(i_a > 0.0 && i_b > 0.0) {
        return Err(format!("Internal rotor: a group without moment about the axis (I_A = {i_a}, I_B = {i_b})."));
    }
    Ok(i_a * i_b / (i_a + i_b))
}

/// Definition of the reduced moment of inertia of an internal rotor computed from the geometry.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum ReducedMomentModel {
    /// Kilpatrick-Pitzer I(1,3) (`reduced_moment_pitzer`), the reduced moment of the MESS deck format.
    #[default]
    Pitzer,
    /// I_A I_B/(I_A + I_B) about the bond axis (`reduced_moment_bond_axis`).
    BondAxis,
}

impl ReducedMomentModel {
    /// Reduced moment (amu Angstrom^2) of the group rotating about the axis (0-based atoms).
    pub fn reduced_moment(self, masses_amu: &[f64], coords_angstrom: &[[f64; 3]], group: &[usize], axis: (usize, usize)) -> Result<f64, String> {
        match self {
            Self::Pitzer => reduced_moment_pitzer(masses_amu, coords_angstrom, group, axis),
            Self::BondAxis => reduced_moment_bond_axis(masses_amu, coords_angstrom, group, axis),
        }
    }
}

/// rho(E) = sum_n rho_rest(E - E_n) on cells of width `cell_cm1`: every level E_n above the ground (E_0 = 0
/// included) shifts the counts by ceil(E_n/cell) cells. The same for sums of states.
pub fn convolve_rotor_levels(counts: &[f64], levels_above_ground_cm1: &[f64], cell_cm1: f64) -> Vec<f64> {
    let mut out = vec![0.0; counts.len()];
    for &e in levels_above_ground_cm1 {
        let shift = (e / cell_cm1 - 1e-9).ceil().max(0.0) as usize;
        if shift >= counts.len() {
            continue;
        }
        for i in shift..counts.len() {
            out[i] += counts[i - shift];
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::CM1_TO_KCAL;

    fn kcal(x: f64) -> f64 {
        x / CM1_TO_KCAL
    }

    fn partition_function(levels_above_ground: &[f64], temperature: f64) -> f64 {
        let kt = crate::constants::KB_CM * temperature;
        levels_above_ground.iter().map(|e| (-e / kt).exp()).sum()
    }

    #[test]
    fn a_free_rotor_has_the_levels_b_sigma_squared_k_squared() {
        let b = 5.0;
        let sigma = 3;
        let free = TorsionalPotential { constant: 0.0, cosine: Vec::new(), sine: Vec::new() };
        let levels = rotor_levels(b, sigma, &free, 21).unwrap();
        let mut expected: Vec<f64> = (0..=10).flat_map(|k| if k == 0 { vec![0.0] } else { vec![(k * k) as f64; 2] }).collect();
        expected.iter_mut().for_each(|x| *x *= b * (sigma * sigma) as f64);
        for (e, x) in levels.iter().zip(&expected) {
            assert!((e - x).abs() < 1e-9, "{e} vs {x}");
        }
    }

    #[test]
    fn equidistant_points_are_interpolated_exactly() {
        // Trigonometric interpolation through N equidistant points on [0, 2 pi) (in x = sigma phi): the series
        // passes through every point, for even and odd N.
        for values in [vec![0.0, 1.2, 3.4, 2.2, 0.7, 0.1], vec![0.0, 0.5, 2.0, 1.0, 0.3]] {
            let v = potential_from_equidistant_points(&values).unwrap();
            let n = values.len();
            for (j, &value) in values.iter().enumerate() {
                let x = 2.0 * std::f64::consts::PI * j as f64 / n as f64;
                assert!((v.evaluate(x) - value).abs() < 1e-12, "N = {n}, j = {j}: {} vs {value}", v.evaluate(x));
            }
        }
    }

    #[test]
    fn the_fourier_expansion_lines_are_constant_then_cosine_and_sine_pairs() {
        // Deck order of `FourierExpansion`: c_0, a_1, b_1, a_2, b_2, ...
        let v = potential_from_fourier_expansion(&[1.0, 2.0, 3.0, 4.0]).unwrap();
        assert_eq!(v, TorsionalPotential { constant: 1.0, cosine: vec![2.0, 4.0], sine: vec![3.0] });
    }

    /// Rotors of CH2CH2OOH (deck and log of the MESS example test_c2h5o2): 12 equidistant points in kcal/mol.
    const CH2_POINTS: [f64; 12] = [0., 0.1161, 0.4338, 0.7406, 0.8283, 0.7185, 0.5845, 0.4924, 0.3602, 0.2390, 0.1397, 0.0540];
    const OOH_POINTS: [f64; 12] = [0., 2.1290, 2.6963, 0.9531, 0.3697, 1.8389, 2.7351, 0.8833, 1.5379, 4.8807, 5.7630, 2.4307];

    #[test]
    fn hindered_rotor_ground_level_and_partition_function_reproduce_the_reference_log() {
        // Reference values printed by MESS for these rotors (B from the geometry, ground energy, quantum statistical
        // weight relative to the ground level): CH2 (sigma 2) B = 9.93852 cm-1, E_0 = 0.155905 kcal/mol,
        // Q = 1.43332 / 3.0102 / 6.64165 at 100 / 300 / 1000 K; OOH (sigma 1) B = 1.94297, E_0 = 0.236392,
        // Q = 1.56215 / 4.43451 / 16.4305.
        for (points, sigma, b, e0, q) in [
            (&CH2_POINTS, 2u32, 9.93852, 0.155905, [1.43332, 3.0102, 6.64165]),
            (&OOH_POINTS, 1u32, 1.94297, 0.236392, [1.56215, 4.43451, 16.4305]),
        ] {
            let v = potential_from_equidistant_points(&points.map(kcal)).unwrap();
            let rotor = HinderedRotor::new(b, sigma, v, DEFAULT_BASIS_SIZE).unwrap();
            assert!((rotor.ground_energy_cm1 * CM1_TO_KCAL - e0).abs() < 2e-6, "E_0 = {}", rotor.ground_energy_cm1 * CM1_TO_KCAL);
            for (t, q_ref) in [100.0, 300.0, 1000.0].into_iter().zip(q) {
                let q_ours = partition_function(&rotor.levels_above_ground_cm1, t);
                assert!((q_ours / q_ref - 1.0).abs() < 2e-5, "T = {t}: {q_ours} vs {q_ref}");
            }
        }
    }

    #[test]
    fn threefold_rotor_levels_agree_with_an_independent_calculation() {
        // Ethane-like rotor, B = 10.70410 cm-1, V = 512 - 512 cos(3 phi) cm-1: the levels of a full 2 pi basis are
        // 150.754 (A), 150.760 (E), 438.302 (E), 438.528 (A), 692.466 (A), 696.193 (E), 900.342 (E), 931.234 (A)
        // cm-1 (reference values of the threefold-rotor unit test of the MESMER distribution). The 2 pi/3-periodic
        // basis keeps the A levels: 150.754, 438.528, 692.466, 931.234.
        let v = TorsionalPotential { constant: 512.0, cosine: vec![-512.0], sine: Vec::new() };
        let rotor = HinderedRotor::new(10.70410, 3, v, DEFAULT_BASIS_SIZE).unwrap();
        assert!((rotor.ground_energy_cm1 - 150.754).abs() < 1.5e-3, "{}", rotor.ground_energy_cm1);
        for (e, reference) in rotor.levels_above_ground_cm1[1..4].iter().zip([438.528, 692.466, 931.234]) {
            assert!((e + rotor.ground_energy_cm1 - reference).abs() < 1.5e-3, "{} vs {reference}", e + rotor.ground_energy_cm1);
        }
    }

    /// Geometry of CH2CH2OOH (MESS example test_c2h5o2, Angstrom) and its masses.
    fn ch2ch2ooh() -> (Vec<f64>, Vec<[f64; 3]>) {
        let symbols = ["C", "C", "O", "O", "H", "H", "H", "H", "H"];
        let coords = vec![
            [-0.5568842010, -0.1350255050, -1.5432632352],
            [0.5679781841, 0.3168021972, -0.6781811609],
            [0.5949497985, -0.3351970823, 0.5931806658],
            [-0.5410965158, 0.1809160868, 1.3404121762],
            [-0.7405637729, -1.1949606730, -1.6648586005],
            [-1.1557232231, 0.5675872950, -2.1045331845],
            [1.5384329412, 0.0351418462, -1.1059503285],
            [0.5542726244, 1.3994065074, -0.5349405010],
            [-1.1834515699, -0.5243371984, 1.1892410097],
        ];
        let symbols: Vec<String> = symbols.iter().map(|s| s.to_string()).collect();
        (crate::utils::atomic_masses::mass_vector_from_symbols_amu(&symbols).unwrap(), coords)
    }

    #[test]
    fn pitzer_reduced_moments_reproduce_the_reference_rotational_constants() {
        // Kilpatrick-Pitzer reduced moment I(1,3) (rotor motion corrected for overall rotation and translation);
        // reference effective rotational constants of the MESS log: CH2 (group 5 6, axis 2 1) 9.93852, OOH (group 4 9,
        // axis 2 3) 1.94297, OH (group 9, axis 3 4) 19.5404 cm-1. Atoms 0-based here.
        // The reference converts with rounded atomic-unit constants: B = 1/(2 I), amu = 1822.8885 m_e, Angstrom =
        // 1.88971616463 bohr, cm-1 = 4.55633e-6 hartree, i.e. B I = 16.8578262 cm-1 amu Angstrom^2; here
        // 16.8576304 (CODATA 2018: 16.8576292). The comparison removes that factor.
        let reference_conversion = 1.0 / (2.0 * 1822.8885 * 1.88971616463_f64.powi(2) * 4.55633e-6);
        let (masses, coords) = ch2ch2ooh();
        // Agreement within half a unit of the last printed digit.
        for (group, axis, b_ref, half_digit) in
            [(vec![4, 5], (1, 0), 9.93852, 5e-6), (vec![3, 8], (1, 2), 1.94297, 5e-6), (vec![8], (2, 3), 19.5404, 5e-5)]
        {
            let i = reduced_moment_pitzer(&masses, &coords, &group, axis).unwrap();
            let b = rotational_constant_cm1(i) * reference_conversion / rotational_constant_cm1(1.0);
            assert!((b - b_ref).abs() <= half_digit, "group {group:?}: B = {b} vs {b_ref}");
        }
    }

    #[test]
    fn bond_axis_reduced_moment_is_the_product_over_the_sum_of_the_group_moments() {
        // I_A I_B/(I_A + I_B) with I_A, I_B the moments of the two groups about the axis line: for ethane (two equal
        // methyl groups) it is I_CH3/2, and both definitions agree for a symmetric top with the axis through the
        // centre of mass.
        let symbols: Vec<String> = ["C", "C", "H", "H", "H", "H", "H", "H"].iter().map(|s| s.to_string()).collect();
        let masses = crate::utils::atomic_masses::mass_vector_from_symbols_amu(&symbols).unwrap();
        let (r, h, z) = (1.02, 1.09_f64, 0.765);
        let mut coords = vec![[0.0, 0.0, -z], [0.0, 0.0, z]];
        for (k, sign) in [(0, -1.0), (1, 1.0)] {
            for j in 0..3 {
                let angle = 2.0 * std::f64::consts::PI * j as f64 / 3.0 + k as f64 * std::f64::consts::PI / 3.0;
                coords.push([r * angle.cos(), r * angle.sin(), sign * (z + (h * h - r * r).sqrt())]);
            }
        }
        let i_ch3 = 3.0 * masses[2] * r * r;
        let bond_axis = reduced_moment_bond_axis(&masses, &coords, &[2, 3, 4], (0, 1)).unwrap();
        assert!((bond_axis / (0.5 * i_ch3) - 1.0).abs() < 1e-12, "{bond_axis} vs {}", 0.5 * i_ch3);
        let pitzer = reduced_moment_pitzer(&masses, &coords, &[2, 3, 4], (0, 1)).unwrap();
        assert!((pitzer / (0.5 * i_ch3) - 1.0).abs() < 1e-12, "{pitzer} vs {}", 0.5 * i_ch3);
    }

    #[test]
    fn rotor_levels_are_convolved_into_the_cell_counts() {
        // rho(E) = sum_n rho_rest(E - E_n), the levels E_n above the ground placed at cell ceil(E_n / dE).
        let rest = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let out = convolve_rotor_levels(&rest, &[0.0, 1.5, 3.0], 1.0);
        // shifts 0, 2, 3
        assert_eq!(out, vec![1.0, 2.0, 4.0, 7.0, 10.0, 13.0]);
    }

    #[test]
    fn the_basis_is_enlarged_until_the_levels_cover_the_energy_range() {
        // heavy twofold rotor: the default basis reaches B (499 sigma)^2 = 0.05 * 998^2 = 49800 cm-1
        let (b, sigma) = (0.05, 2);
        let potential = TorsionalPotential { constant: 400.0, cosine: vec![-400.0], sine: Vec::new() };
        let small = HinderedRotor::for_energy_range(b, sigma, potential.clone(), DEFAULT_BASIS_SIZE, MAX_BASIS_SIZE, 30_000.0).unwrap();
        assert_eq!(small, HinderedRotor::new(b, sigma, potential.clone(), DEFAULT_BASIS_SIZE).unwrap());

        // 60000 cm-1 needs K = ceil(sqrt((60000 + E0 - V_min)/B)/sigma) = 548, size 2K + 1 = 1097
        let large = HinderedRotor::for_energy_range(b, sigma, potential.clone(), DEFAULT_BASIS_SIZE, MAX_BASIS_SIZE, 60_000.0).unwrap();
        assert_eq!(large, HinderedRotor::new(b, sigma, potential.clone(), 1097).unwrap());
        assert!(large.highest_level_above_ground_cm1() >= 60_000.0);

        let err = HinderedRotor::for_energy_range(b, sigma, potential, DEFAULT_BASIS_SIZE, 1001, 60_000.0).unwrap_err();
        assert!(err.contains("HamiltonSizeMax"), "{err}");

        // Above the barrier the levels lie near B k^2 + c_0: with K = ceil(sqrt((E + E0 - V_min)/B)) = 71 the
        // highest level kept (k = 63) ends below E = 5000 cm-1; the basis grows until the levels reach E.
        let potential = TorsionalPotential { constant: 1000.0, cosine: vec![-1000.0], sine: Vec::new() };
        let rotor = HinderedRotor::for_energy_range(1.0, 1, potential.clone(), 101, MAX_BASIS_SIZE, 5000.0).unwrap();
        assert!(HinderedRotor::new(1.0, 1, potential, 143).unwrap().highest_level_above_ground_cm1() < 5000.0);
        assert!(rotor.highest_level_above_ground_cm1() >= 5000.0, "{}", rotor.highest_level_above_ground_cm1());
    }

    #[test]
    fn the_reduced_moment_model_selects_the_definition() {
        // CH2CH2OOH: the OOH group off the principal axes, where the two definitions differ
        let (masses, coords) = ch2ch2ooh();
        let (group, axis) = (vec![3, 8], (1, 2));
        assert_eq!(ReducedMomentModel::default(), ReducedMomentModel::Pitzer);
        let pitzer = ReducedMomentModel::Pitzer.reduced_moment(&masses, &coords, &group, axis).unwrap();
        let bond_axis = ReducedMomentModel::BondAxis.reduced_moment(&masses, &coords, &group, axis).unwrap();
        assert_eq!(pitzer, reduced_moment_pitzer(&masses, &coords, &group, axis).unwrap());
        assert_eq!(bond_axis, reduced_moment_bond_axis(&masses, &coords, &group, axis).unwrap());
        assert!((pitzer / bond_axis - 1.0).abs() > 1e-3, "{pitzer} vs {bond_axis}");
    }
}
