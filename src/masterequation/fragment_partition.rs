//! Dissociation of a well C into a fragment B that is itself a well of the network and a partner A in excess
//! (reports/nonthermal_sources_design.md, Section 15.1). Green, Robertson, Chem. Phys. Lett. 605-606, 44 (2014) (GR14):
//! the fragment receives the energy e_j <= X_i of the C grain i above the pair asymptote with the probability P(e_j | X_i)
//! (their partition matrix Q, columns summing to 1, eq. 15); the reverse association follows from detailed balance, which
//! holds for any Q (eqs. 16-17), with the partner at its concentration [A] in the weight of the fragment states:
//!   f_B(j) = rho_B(e_j) exp(-(E_B0 + e_j)/kT) phi_A,   phi_A = Q_A,int exp(-E_A0/kT) C'(mu) (kT)^(3/2) / [A],
//! so that [C]/[B] = K_c [A] in equilibrium. Kernels:
//!   - `Prior`: P ~ rho_B(e) [rho_A (x) rho_t](X - e), rho_t ~ E^(1/2) the relative translation (GR14 eq. 18; the same
//!     form as Sztaray, Bodi, Baer, J. Mass Spectrom. 45, 1233 (2010), eq. 5, `photoion::product_energy`);
//!   - `ModifiedPrior`: rho_B^n with n = n0 (T/T_ref)^m (Shannon, Blitz, Seakins, J. Phys. Chem. A 128, 1501 (2024), and
//!     their ref. 5);
//!   - `TwoPieceGaussian`: A exp(-(e - mu)^2/(2 sigma_L^2)) for e <= mu and A exp(-(e - mu)^2/(2 sigma_R^2)) for e > mu,
//!     mu and the widths linear in X (Shannon et al. 2024, eq. 4 and SI; the sides named explicitly, because the paper's
//!     sigma_1/sigma_2 labels conflict with its deck).
//! Grain j of the fragment covers [(j - 1/2) dE, (j + 1/2) dE) above its bottom, clipped to [0, X]; P is normalized over
//! the grains of the fragment grid (energy beyond the grid top is removed before the normalization).

use crate::barrierless::ilt::ilt_barrierless::translational_partition_constant;
use crate::constants::KB_CM;
use crate::numeric::special_functions::normal_interval_probability;
use crate::photoion::product_energy::{convolve, translational_density};

/// The partner in excess: its concentration and what its partition function needs.
#[derive(Debug, Clone, PartialEq)]
pub struct ExcessPartner {
    pub name: String,
    /// [A] in molecule cm-3.
    pub concentration_cm3: f64,
    /// Ground energy of the partner on the energy scale of the deck (cm-1).
    pub ground_energy_cm1: f64,
    /// Internal density of states of the partner (per cm-1) on cells from its ground state.
    pub density_cells: Vec<f64>,
    pub cell_width_cm1: f64,
    /// Reduced mass of fragment and partner (amu), for the relative translation.
    pub reduced_mass_amu: f64,
}

impl ExcessPartner {
    /// ln phi_A(T) = ln Q_A,int - E_A0/kT + ln(C'(mu) (kT)^(3/2)) - ln [A].
    pub fn log_weight_factor(&self, temperature_kelvin: f64) -> f64 {
        let kt = KB_CM * temperature_kelvin;
        let d = self.cell_width_cm1;
        let internal: f64 = self.density_cells.iter().enumerate().map(|(c, r)| r * (-(c as f64 * d) / kt).exp()).sum::<f64>() * d;
        internal.ln() - self.ground_energy_cm1 / kt + (translational_partition_constant(self.reduced_mass_amu) * kt.powf(1.5)).ln()
            - self.concentration_cm3.ln()
    }

    pub fn validate(&self) -> Result<(), String> {
        if !(self.concentration_cm3 > 0.0 && self.concentration_cm3.is_finite()) {
            return Err(format!("Partner '{}': the concentration must be positive.", self.name));
        }
        if !(self.cell_width_cm1 > 0.0 && self.reduced_mass_amu > 0.0) || self.density_cells.iter().all(|r| *r == 0.0) {
            return Err(format!("Partner '{}': needs a cell width, a reduced mass and a density of states.", self.name));
        }
        Ok(())
    }
}

/// The energy partitioning P(e | X) of the fragment.
#[derive(Debug, Clone, PartialEq)]
pub enum FragmentKernel {
    /// GR14 eq. 18: `fragment_cells` rho_B and `remainder_cells` rho_A (x) rho_t (`prior_remainder`) on cells from 0.
    Prior { fragment_cells: Vec<f64>, remainder_cells: Vec<f64>, cell_width_cm1: f64 },
    /// rho_B^n rho_A (x) rho_t with n = order (T/T_ref)^temperature_exponent.
    ModifiedPrior {
        fragment_cells: Vec<f64>,
        remainder_cells: Vec<f64>,
        cell_width_cm1: f64,
        order: f64,
        temperature_exponent: f64,
        reference_temperature_kelvin: f64,
    },
    /// [intercept, gradient] of mu, sigma_L and sigma_R in cm-1 against X in cm-1.
    TwoPieceGaussian { mu: [f64; 2], sigma_left: [f64; 2], sigma_right: [f64; 2] },
    /// Everything into grain 0: the kernel of a lumped thermal reactant state, a fragment "well" of one grain
    /// (reports/nonthermal_sources_design.md, Section 15.2).
    LowestGrain,
}

/// rho_A (x) rho_t on `cells` cells: the partner density convolved with the 3D relative translational density (both on
/// cells from 0; constant factors cancel in the normalized kernel).
pub fn prior_remainder(partner_cells: &[f64], cells: usize, cell_width_cm1: f64) -> Result<Vec<f64>, String> {
    let rho_t = translational_density(cells, cell_width_cm1, 3)?;
    Ok(convolve(partner_cells, &rho_t, cells))
}

/// Grain of the fragment that contains the energy `e` above its bottom.
fn grain_of(e: f64, grain_width_cm1: f64) -> usize {
    ((e + 0.5 * grain_width_cm1) / grain_width_cm1).floor().max(0.0) as usize
}

/// Normalized (grain, probability) pairs from unnormalized grain masses; all in the lowest grain if they vanish.
fn normalized(masses: Vec<f64>) -> Vec<(usize, f64)> {
    let total: f64 = masses.iter().sum();
    if !(total > 0.0) {
        return vec![(0, 1.0)];
    }
    masses.into_iter().enumerate().filter(|(_, m)| *m > 0.0).map(|(j, m)| (j, m / total)).collect()
}

impl FragmentKernel {
    /// P over the fragment grains 0 .. `fragment_grains` for the energy `excess_cm1` above the pair asymptote; empty below
    /// the asymptote.
    pub fn distribution(&self, excess_cm1: f64, grain_width_cm1: f64, fragment_grains: usize, temperature_kelvin: f64) -> Result<Vec<(usize, f64)>, String> {
        if excess_cm1 < 0.0 {
            return Ok(Vec::new());
        }
        let top = fragment_grains.min(grain_of(excess_cm1, grain_width_cm1) + 1);
        match self {
            FragmentKernel::Prior { fragment_cells, remainder_cells, cell_width_cm1 } => {
                prior(fragment_cells, remainder_cells, *cell_width_cm1, 1.0, excess_cm1, grain_width_cm1, top)
            }
            FragmentKernel::ModifiedPrior { fragment_cells, remainder_cells, cell_width_cm1, order, temperature_exponent, reference_temperature_kelvin } => {
                let n = order * (temperature_kelvin / reference_temperature_kelvin).powf(*temperature_exponent);
                prior(fragment_cells, remainder_cells, *cell_width_cm1, n, excess_cm1, grain_width_cm1, top)
            }
            FragmentKernel::LowestGrain => Ok(vec![(0, 1.0)]),
            FragmentKernel::TwoPieceGaussian { mu, sigma_left, sigma_right } => {
                let centre = mu[0] + mu[1] * excess_cm1;
                let (left, right) = (sigma_left[0] + sigma_left[1] * excess_cm1, sigma_right[0] + sigma_right[1] * excess_cm1);
                if !(left > 0.0 && right > 0.0) {
                    return Err(format!(
                        "Two-piece Gaussian fragment energy: the widths at X = {excess_cm1:.1} cm-1 are {left:.1} and {right:.1} cm-1; \
                         both must be positive."
                    ));
                }
                // Mass of [a, b] relative to A sqrt(2 pi): sigma_L P_N(mu, sigma_L)([a, min(b, mu)]) + sigma_R P_N(mu, sigma_R)([max(a, mu), b]).
                let mass = |a: f64, b: f64| -> f64 {
                    let l = if a < centre { left * normal_interval_probability(a, b.min(centre), centre, left) } else { 0.0 };
                    let r = if b > centre { right * normal_interval_probability(a.max(centre), b, centre, right) } else { 0.0 };
                    l + r
                };
                let masses = (0..top)
                    .map(|j| {
                        let a = ((j as f64 - 0.5) * grain_width_cm1).max(0.0);
                        let b = ((j as f64 + 0.5) * grain_width_cm1).min(excess_cm1);
                        if b > a { mass(a, b) } else if j == 0 { 1.0 } else { 0.0 }
                    })
                    .collect();
                Ok(normalized(masses))
            }
        }
    }
}

/// rho_B(c)^n R(N - c) on the cells c = 0 .. N (N the excess energy in cells), summed into the grains below `top`.
fn prior(fragment: &[f64], remainder: &[f64], cell: f64, n: f64, excess_cm1: f64, grain_width_cm1: f64, top: usize) -> Result<Vec<(usize, f64)>, String> {
    let cells = (excess_cm1 / cell).round() as usize;
    if cells >= fragment.len() || cells >= remainder.len() {
        return Err(format!(
            "Prior fragment energy: X = {excess_cm1:.1} cm-1 lies beyond the {} cells of the densities of states.",
            fragment.len().min(remainder.len())
        ));
    }
    let mut masses = vec![0.0; top];
    for c in 0..=cells {
        let j = grain_of(c as f64 * cell, grain_width_cm1);
        if j < top {
            let rho = if n == 1.0 { fragment[c] } else { fragment[c].powf(n) };
            masses[j] += rho * remainder[cells - c];
        }
    }
    Ok(normalized(masses))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::barrierless::ilt::ilt_barrierless::translational_partition_constant;
    use crate::constants::KB_CM;
    use crate::photoion::product_energy::{fragment_energy_distribution, translational_density};

    /// 1 kJ/mol in cm-1 (1 kcal = 4.184 kJ).
    const KJ: f64 = 1.0 / (2.85914e-3 * 4.184);

    fn sum(p: &[(usize, f64)]) -> f64 {
        p.iter().map(|(_, x)| x).sum()
    }

    #[test]
    fn the_two_piece_gaussian_reproduces_the_prompt_fractions_of_the_oh_glyoxal_deck() {
        // Shannon, Blitz, Seakins, J. Phys. Chem. A 128, 1501 (2024), eq. 4 with the coefficients of their SI deck (cm-1,
        // linear in the energy X above the HC(O)CO + H2O asymptote; the wider width at high energy on the right, as their
        // code applies it). X = 130.2 kJ/mol + E above TS1; the fraction of HC(O)CO above TS2 (37.3 kJ/mol) is the prompt
        // fraction: 0.3319, 0.7053, 0.8997, 0.9960 at 5, 15, 25, 50 kJ/mol (closed form, both widths).
        let kernel = FragmentKernel::TwoPieceGaussian {
            mu: [-2685.4603641118783, 0.50158003],
            sigma_left: [-132.00846000618344, 0.053485],
            sigma_right: [-1122.110663126146, 0.13401113],
        };
        let de = 10.0;
        for (e, expected) in [(5.0, 0.3319), (15.0, 0.7053), (25.0, 0.8997), (50.0, 0.9960)] {
            let x = (130.2 + e) * KJ;
            let p = kernel.distribution(x, de, 2000, 298.0).unwrap();
            assert!((sum(&p) - 1.0).abs() < 1e-12);
            assert!(p.iter().all(|&(j, _)| (j as f64 - 0.5) * de <= x), "no energy above X");
            let prompt: f64 = p.iter().filter(|&&(j, _)| j as f64 * de > 37.3 * KJ).map(|(_, x)| x).sum();
            assert!((prompt - expected).abs() < 3e-3, "E = {e} kJ/mol: {prompt} vs {expected}");
        }
    }

    #[test]
    fn the_prior_is_the_normalized_product_of_the_fragment_density_and_the_remainder() {
        // Green, Robertson, Chem. Phys. Lett. 605-606, 44 (2014), eq. 18: P(e|X) ~ rho_B(e) [rho_A (x) rho_t](X - e), the same
        // form as the statistical partitioning of the photoion module (Sztaray, Bodi, Baer 2010, eq. 5), on the cells and
        // summed into the grains of the fragment.
        let cells = 6000;
        let rho_b: Vec<f64> = (0..cells).map(|c| (c as f64).powi(3)).collect();
        let rho_a: Vec<f64> = (0..cells).map(|c| (c as f64).powi(1)).collect();
        let remainder = prior_remainder(&rho_a, cells, 1.0).unwrap();
        let kernel = FragmentKernel::Prior { fragment_cells: rho_b.clone(), remainder_cells: remainder, cell_width_cm1: 1.0 };
        let (x, de) = (5000.0, 20.0);
        let p = kernel.distribution(x, de, 300, 300.0).unwrap();
        assert!((sum(&p) - 1.0).abs() < 1e-12);
        let cellwise = fragment_energy_distribution(&rho_b, &rho_a, &translational_density(cells, 1.0, 3).unwrap(), 5000).unwrap();
        for &(j, pj) in &p {
            let from_cells: f64 = cellwise.iter().enumerate().filter(|(c, _)| ((*c as f64 + 0.5 * de) / de).floor() as usize == j).map(|(_, v)| v).sum();
            assert!((pj - from_cells).abs() < 1e-12, "grain {j}: {pj} vs {from_cells}");
        }
        // Classical limit: rho_B ~ e^3, rho_A (x) rho_t ~ (X - e)^(2 + 3/2): a Beta distribution with the mean fraction
        // 4/(4 + 3.5) of X in the fragment.
        let mean: f64 = p.iter().map(|&(j, pj)| j as f64 * de * pj).sum::<f64>() / x;
        assert!((mean - 4.0 / 7.5).abs() < 2e-3, "{mean}");
    }

    #[test]
    fn a_modified_prior_of_order_one_is_the_prior_and_a_lower_order_puts_less_energy_into_the_fragment() {
        // Shannon et al. 2024 (their ref. 5): rho_B^n, n = n0 (T/T_ref)^m.
        let cells = 4000;
        let rho_b: Vec<f64> = (0..cells).map(|c| 1.0 + (c as f64).powi(4)).collect();
        let rho_a: Vec<f64> = (0..cells).map(|c| 1.0 + c as f64).collect();
        let remainder = prior_remainder(&rho_a, cells, 1.0).unwrap();
        let prior = FragmentKernel::Prior { fragment_cells: rho_b.clone(), remainder_cells: remainder.clone(), cell_width_cm1: 1.0 };
        let modified = |n0: f64, m: f64| FragmentKernel::ModifiedPrior {
            fragment_cells: rho_b.clone(),
            remainder_cells: remainder.clone(),
            cell_width_cm1: 1.0,
            order: n0,
            temperature_exponent: m,
            reference_temperature_kelvin: 298.0,
        };
        let (x, de) = (3000.0, 10.0);
        let a = prior.distribution(x, de, 400, 500.0).unwrap();
        let b = modified(1.0, 0.0).distribution(x, de, 400, 500.0).unwrap();
        assert_eq!(a.len(), b.len());
        assert!(a.iter().zip(&b).all(|(p, q)| p.0 == q.0 && (p.1 - q.1).abs() < 1e-14));
        let mean = |p: &[(usize, f64)]| p.iter().map(|&(j, v)| j as f64 * v).sum::<f64>();
        let low = modified(0.27, 0.0).distribution(x, de, 400, 500.0).unwrap();
        assert!(mean(&low) < mean(&a), "{} vs {}", mean(&low), mean(&a));
        // n = n0 (T/T_ref)^m: order 0.5 with m = 1 at 2 T_ref is order 1.
        let scaled = modified(0.5, 1.0).distribution(x, de, 400, 596.0).unwrap();
        assert!(a.iter().zip(&scaled).all(|(p, q)| (p.1 - q.1).abs() < 1e-14));
    }

    #[test]
    fn energies_at_and_below_the_asymptote_go_into_the_lowest_grain_or_nowhere() {
        let kernel = FragmentKernel::TwoPieceGaussian { mu: [500.0, 0.0], sigma_left: [100.0, 0.0], sigma_right: [100.0, 0.0] };
        assert!(kernel.distribution(-1.0, 10.0, 100, 300.0).unwrap().is_empty());
        assert_eq!(kernel.distribution(3.0, 10.0, 100, 300.0).unwrap(), vec![(0, 1.0)]);
        // All grains below X, also when X lies beyond the fragment grid: the part above the grid is removed and the rest
        // renormalized.
        let p = kernel.distribution(5000.0, 10.0, 40, 300.0).unwrap();
        assert!((sum(&p) - 1.0).abs() < 1e-12 && p.iter().all(|&(j, _)| j < 40));
    }

    #[test]
    fn the_partner_weight_is_its_partition_function_per_volume_over_its_concentration() {
        // phi_A = Q_A,int exp(-E_A0/kT) C'(mu) (kT)^(3/2) / [A] (Green, Robertson 2014, detailed balance of eq. 15 with the
        // partner in excess); one oscillator of 1000 cm-1: Q_int = 1/(1 - exp(-1000/kT)).
        let mut density_cells = vec![0.0; 20001];
        for v in 0..=20 {
            density_cells[v * 1000] = 1.0;
        }
        let partner = ExcessPartner {
            name: "A".into(),
            concentration_cm3: 1e14,
            ground_energy_cm1: 150.0,
            density_cells,
            cell_width_cm1: 1.0,
            reduced_mass_amu: 12.0,
        };
        let t = 400.0;
        let kt = KB_CM * t;
        let expected = -(1.0 - (-1000.0 / kt).exp()).ln() - 150.0 / kt + (translational_partition_constant(12.0) * kt.powf(1.5)).ln() - 1e14_f64.ln();
        assert!((partner.log_weight_factor(t) - expected).abs() < 1e-10, "{} vs {expected}", partner.log_weight_factor(t));
    }

    #[test]
    fn the_lowest_grain_kernel_puts_everything_into_grain_zero() {
        // The kernel of a lumped thermal reactant state (one grain): reports/nonthermal_sources_design.md, Section 15.2.
        let k = FragmentKernel::LowestGrain;
        assert_eq!(k.distribution(1234.0, 10.0, 1, 300.0).unwrap(), vec![(0, 1.0)]);
        assert!(k.distribution(-1.0, 10.0, 1, 300.0).unwrap().is_empty());
    }
}
