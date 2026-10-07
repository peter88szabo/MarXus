//! Statistical partitioning of the excess energy of a dissociation between the fragment ion, the neutral fragment
//! and their relative translation (SBB10 eqs. 3-5): ABC+(E) -> AB+(E_i) + C(E_n), E - E0 = E_i + E_n + E_trans.
//! Only the vibrational frequencies and the rotational degrees of freedom of the fragments enter; there are no
//! adjustable parameters besides the number of translational degrees of freedom (SBB10 p. 1237).
//!
//! Reference: B. Sztáray, A. Bodi, T. Baer, J. Mass Spectrom. 45, 1233 (2010) (SBB10).

/// Classical density of `dof` translational degrees of freedom (1, 2 or 3; phase space theory gives 2, SBB10
/// p. 1237) on cells: rho ~ E^(dof/2 - 1), cell i holding the states in ((i-1) dE, i dE], cell 0 empty. The
/// constant factor cancels in the normalized partitioning.
pub fn translational_density(cells: usize, cell_cm1: f64, dof: u32) -> Result<Vec<f64>, String> {
    if !(1..=3).contains(&dof) {
        return Err(format!("Translational density: {dof} degrees of freedom; 1, 2 or 3 are possible."));
    }
    let cumulative = |i: usize| (i as f64 * cell_cm1).powf(dof as f64 / 2.0);
    Ok((0..cells).map(|i| if i == 0 { 0.0 } else { cumulative(i) - cumulative(i - 1) }).collect())
}

/// Discrete convolution on cells, c[m] = sum_{j=0}^{m} a[j] b[m - j], for m < len.
pub fn convolve(a: &[f64], b: &[f64], len: usize) -> Vec<f64> {
    let mut c = vec![0.0; len];
    for (j, x) in a.iter().enumerate().take(len).filter(|(_, x)| **x != 0.0) {
        for (k, y) in b.iter().enumerate().take(len - j) {
            c[j + k] += x * y;
        }
    }
    c
}

/// SBB10 eq. 5: the normalized internal energy distribution P(E_i | E - E0) of the fragment ion for an excess energy
/// of `excess_cells` cells, E_i = 0 .. E - E0:
///   P(E_i, E - E0) = rho_F(E_i) [rho_N (x) rho_tr](E - E0 - E_i) / sum over E_i,
/// with rho_F, rho_N the rovibrational densities of the fragment ion and the neutral fragment and rho_tr the
/// translational density. Where the discrete count has no product states (zero excess energy, and the lowest cells,
/// where classical rotational and translational densities vanish at E = 0) the ion is formed in its lowest cell.
pub fn fragment_energy_distribution(rho_f: &[f64], rho_n: &[f64], rho_tr: &[f64], excess_cells: usize) -> Result<Vec<f64>, String> {
    let n = excess_cells;
    check_lengths(rho_f, rho_n, rho_tr, n + 1)?;
    let rest = convolve(rho_n, rho_tr, n + 1);
    let mut p: Vec<f64> = (0..=n).map(|i| rho_f[i] * rest[n - i]).collect();
    let total: f64 = p.iter().sum();
    if total == 0.0 {
        let mut ground = vec![0.0; n + 1];
        ground[0] = 1.0;
        return Ok(ground);
    }
    if !total.is_finite() || total < 0.0 {
        return Err(format!("Product energy distribution: invalid densities of states at {n} cells of excess energy."));
    }
    p.iter_mut().for_each(|x| *x /= total);
    Ok(p)
}

/// The internal energy distribution of the fragment ion formed from a parent ion distribution: SBB10 eq. 5 "calculated
/// for each E in the energy distribution of the parent ion and then summed over the parent P(E) distribution" (SBB10
/// p. 1236):
///   D(E_i) = sum_{E >= E0} w(E) P(E_i, E - E0),
/// with w(E) the parent probability per cell times the fraction at E that forms this fragment, and E0 the
/// dissociation limit in cells on the parent's energy scale. D is not normalized: its sum is sum_{E >= E0} w(E).
/// Cell i of the result is the fragment ion energy i dE above its ground state.
pub fn daughter_distribution(weights: &[f64], e0_cells: usize, rho_f: &[f64], rho_n: &[f64], rho_tr: &[f64]) -> Result<Vec<f64>, String> {
    let len = weights.len().saturating_sub(e0_cells);
    ProductPartitioning::new(rho_f, rho_n, rho_tr, len)?.daughter(weights, e0_cells)
}

/// The convolutions of SBB10 eq. 5 for excess energies up to `len` cells, computed once for many parent
/// distributions: R = rho_N (x) rho_tr and the normalization rho_F (x) R.
#[derive(Debug, Clone)]
pub struct ProductPartitioning {
    rho_f: Vec<f64>,
    rest: Vec<f64>,
    all: Vec<f64>,
}

impl ProductPartitioning {
    pub fn new(rho_f: &[f64], rho_n: &[f64], rho_tr: &[f64], len: usize) -> Result<Self, String> {
        check_lengths(rho_f, rho_n, rho_tr, len)?;
        let rest = convolve(rho_n, rho_tr, len);
        let all = convolve(rho_f, &rest, len);
        Ok(Self { rho_f: rho_f[..len].to_vec(), rest, all })
    }

    /// `daughter_distribution` with the precomputed convolutions: D(i) = rho_F(i) sum_{n >= i} q(n) R(n - i), with
    /// q(n) = w(E0 + n) / (rho_F (x) R)(n); excess energies without discrete product states go to the lowest cell.
    pub fn daughter(&self, weights: &[f64], e0_cells: usize) -> Result<Vec<f64>, String> {
        if e0_cells >= weights.len() {
            return Ok(Vec::new());
        }
        let len = weights.len() - e0_cells;
        if len > self.all.len() {
            return Err(format!(
                "Product energy distribution: excess energies up to {len} cells, the partitioning covers {}.",
                self.all.len()
            ));
        }
        let mut daughter = vec![0.0; len];
        let mut q = vec![0.0; len];
        for n in 0..len {
            let w = weights[e0_cells + n];
            if w == 0.0 {
                continue;
            }
            if self.all[n] > 0.0 {
                q[n] = w / self.all[n];
            } else {
                // no discrete product states at this excess energy: lowest cell (`fragment_energy_distribution`)
                daughter[0] += w;
            }
        }
        let Some(high) = q.iter().rposition(|x| *x != 0.0) else {
            return Ok(daughter);
        };
        let low = q.iter().position(|x| *x != 0.0).unwrap_or(0);
        for i in 0..=high {
            if self.rho_f[i] != 0.0 {
                daughter[i] += self.rho_f[i] * (i.max(low)..=high).map(|n| q[n] * self.rest[n - i]).sum::<f64>();
            }
        }
        Ok(daughter)
    }
}

fn check_lengths(rho_f: &[f64], rho_n: &[f64], rho_tr: &[f64], len: usize) -> Result<(), String> {
    if rho_f.len() < len || rho_n.len() < len || rho_tr.len() < len {
        return Err(format!(
            "Product energy distribution: the densities of states must cover {len} cells (fragment ion {}, neutral {}, \
             translation {}).",
            rho_f.len(),
            rho_n.len(),
            rho_tr.len()
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn brute_force(rho_f: &[f64], rho_n: &[f64], rho_tr: &[f64], n: usize) -> Vec<f64> {
        let mut p: Vec<f64> = (0..=n)
            .map(|i| rho_f[i] * (0..=n - i).map(|x| rho_n[x] * rho_tr[n - i - x]).sum::<f64>())
            .collect();
        let total: f64 = p.iter().sum();
        p.iter_mut().for_each(|v| *v /= total);
        p
    }

    fn rough(len: usize, seed: u64) -> Vec<f64> {
        (0..len).map(|i| ((i as u64 * 2654435761 + seed) % 97) as f64 + 0.5).collect()
    }

    #[test]
    fn translational_densities_are_classical_power_laws_on_cells() {
        // d translational degrees of freedom: rho ~ E^(d/2 - 1); cell i holds the states in ((i-1) dE, i dE]
        let two = translational_density(6, 1.0, 2).unwrap();
        assert_eq!(two, vec![0.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
        let three = translational_density(6, 1.0, 3).unwrap();
        for i in 1..6 {
            let expected = (i as f64).powf(1.5) - ((i - 1) as f64).powf(1.5);
            assert!((three[i] / three[1] - expected).abs() < 1e-12);
        }
        let one = translational_density(6, 1.0, 1).unwrap();
        assert!((one[4] / one[1] - (2.0 - 3.0_f64.sqrt())).abs() < 1e-12);
        assert!(translational_density(6, 1.0, 0).is_err() && translational_density(6, 1.0, 4).is_err());
    }

    #[test]
    fn the_fragment_ion_distribution_is_the_normalized_product_with_the_convolved_rest() {
        // SBB10 eq. 5: P(E_i, E - E0) = rho_F(E_i) [rho_N (x) rho_tr](E - E0 - E_i) / normalization
        let (rho_f, rho_n, rho_tr) = (rough(300, 1), rough(300, 7), translational_density(300, 1.0, 2).unwrap());
        for n in [0, 1, 37, 250] {
            let p = fragment_energy_distribution(&rho_f, &rho_n, &rho_tr, n).unwrap();
            let expected = if n == 0 { vec![1.0] } else { brute_force(&rho_f, &rho_n, &rho_tr, n) };
            assert_eq!(p.len(), n + 1);
            for (a, b) in p.iter().zip(&expected) {
                assert!((a - b).abs() < 1e-12, "n = {n}");
            }
        }
    }

    #[test]
    fn atom_loss_with_two_translational_degrees_of_freedom_keeps_the_shape_of_the_fragment_density() {
        // Neutral atom (one state at zero) and 2D translation: P(E_i | E) ~ rho_F(E_i) for E_i < E.
        let n = 200;
        let rho_f: Vec<f64> = (0..=n).map(|i| (i as f64).powi(3)).collect();
        let mut atom = vec![0.0; n + 1];
        atom[0] = 1.0;
        let p = fragment_energy_distribution(&rho_f, &atom, &translational_density(n + 1, 1.0, 2).unwrap(), n).unwrap();
        let total: f64 = (0..n).map(|i| (i as f64).powi(3)).sum();
        for i in [0, 5, 100, 199] {
            assert!((p[i] - (i as f64).powi(3) / total).abs() < 1e-15);
        }
        assert_eq!(p[n], 0.0);
    }

    #[test]
    fn excess_energies_without_discrete_product_states_form_the_fragment_in_its_ground_cell() {
        // Classical rotors leave the ground cell of the fragment empty and the translational density is zero at E = 0:
        // at 1 cell of excess energy the discrete count has no product states; the fragment ion is formed in its lowest
        // cell, as at zero excess energy.
        let rho_f = vec![0.0, 0.3, 0.6, 0.8, 1.0];
        let mut atom = vec![0.0; 5];
        atom[0] = 4.0;
        let tr = translational_density(5, 1.0, 2).unwrap();
        assert_eq!(fragment_energy_distribution(&rho_f, &atom, &tr, 1).unwrap(), vec![1.0, 0.0]);
        let mut parent = vec![0.0; 5];
        parent[3] = 0.7;
        let daughter = daughter_distribution(&parent, 2, &rho_f, &atom, &tr).unwrap();
        assert_eq!(daughter, vec![0.7, 0.0, 0.0]);
    }

    #[test]
    fn the_daughter_distribution_sums_the_partitioning_over_the_parent_distribution() {
        // "These distributions are calculated for each E in the energy distribution of the parent ion and then
        // summed over the parent P(E) distribution" (SBB10 p. 1236).
        let (rho_f, rho_n, rho_tr) = (rough(400, 3), rough(400, 11), translational_density(400, 1.0, 3).unwrap());
        let parent = rough(400, 5);
        let e0 = 120;
        let daughter = daughter_distribution(&parent, e0, &rho_f, &rho_n, &rho_tr).unwrap();
        let mut expected = vec![0.0; 400 - e0];
        for e in e0..400 {
            let partition = if e == e0 { vec![1.0] } else { brute_force(&rho_f, &rho_n, &rho_tr, e - e0) };
            for (i, p) in partition.iter().enumerate() {
                expected[i] += parent[e] * p;
            }
        }
        assert_eq!(daughter.len(), expected.len());
        for (a, b) in daughter.iter().zip(&expected) {
            assert!((a - b).abs() < 1e-10 * b.abs().max(1.0));
        }
        let weight: f64 = parent[e0..].iter().sum();
        assert!((daughter.iter().sum::<f64>() - weight).abs() < 1e-9 * weight);
    }

    #[test]
    fn a_thermal_parent_with_a_statistical_rate_gives_a_thermal_daughter() {
        // (derived) With k(E) rho(E) proportional to [rho_F (x) rho_N (x) rho_tr](E - E0) and a Boltzmann parent, the
        // weights k rho exp(-E/kT) give a daughter distribution proportional to rho_F(E_i) exp(-E_i/kT).
        // the Boltzmann tail cut by the grid top is below 1e-20 of the retained part
        let (len, e0, kt) = (12_000, 1500, 150.0);
        let rho_f: Vec<f64> = (0..len).map(|i| ((i + 1) as f64).powf(2.5)).collect();
        let rho_n: Vec<f64> = (0..len).map(|i| ((i + 1) as f64).powf(0.5)).collect();
        let rho_tr = translational_density(len, 1.0, 2).unwrap();
        let rest = convolve(&rho_n, &rho_tr, len);
        let all = convolve(&rho_f, &rest, len);
        let weights: Vec<f64> =
            (0..len).map(|e| if e < e0 { 0.0 } else { all[e - e0] * (-(e as f64) / kt).exp() }).collect();
        let daughter = daughter_distribution(&weights, e0, &rho_f, &rho_n, &rho_tr).unwrap();
        let thermal: Vec<f64> = (0..daughter.len()).map(|i| rho_f[i] * (-(i as f64) / kt).exp()).collect();
        let (sd, st): (f64, f64) = (daughter[..2000].iter().sum(), thermal[..2000].iter().sum());
        for i in [10, 400, 1000, 1999] {
            assert!((daughter[i] / sd / (thermal[i] / st) - 1.0).abs() < 1e-9, "cell {i}");
        }
    }
}
