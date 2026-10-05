//! Steady-state solution J N^s = R F of the chemical-activation master equation.
//!
//! Steady state, dN/dt = 0 in PO14 eq. 2: R F = J N^s, N^s = R J^-1 F (PO14 eq. 5; O02 eq. 7; GO10
//! eq. 7). Olzmann and co-workers obtain N^s "by directly solving R F = J N^s" with a band solver
//! (PO14 p. 235; GO10 p. 12293); here R = 1 and F is normalized, so N^s are the populations per unit
//! formation rate and sum_r (K_r N^s) are the yields (O02 eq. 10).
//!
//! Symmetrization: with D = diag(sqrt(f)), f = rho exp(-E/kT) on the absolute energy scale,
//!   S = D^-1 J D,   S_rc = J_rc sqrt(f_c/f_r),
//! is symmetric when J obeys detailed balance, J_rc f_c = J_cr f_r (see `chemical_activation_operator.rs`;
//! R19 eqs. 5.74-5.77: S = F^-1 M F with F = diag(b^1/2) is a similarity transform, so S and J have the
//! same eigenvalues; Olzmann also works with the symmetrized matrix, PO14 p. 235). S is then positive
//! definite (the eigenvalues of J are positive, GO10 text before eq. 12) and banded in
//! the energy ordering of the states, so S y = D^-1 F is solved by a banded Cholesky factorization and
//! N = D y. A few steps of iterative refinement on the original equation J N = F remove the effect
//! of rounding in the symmetrized matrix. Without exact detailed balance (e.g. isomerization rates that
//! are not microscopically reversible) the symmetric solver is refused and BiCGSTAB on D^-1 J D can be
//! used instead.

use crate::numeric::banded_solvers::SymmetricBandMatrix;
use crate::numeric::krylov::{solve_bicgstab_left_preconditioned, JacobiPreconditioner, LinearOperator};

use super::chemical_activation_operator::ChemicalActivationOperator;

/// Largest relative asymmetry |S_rc - S_cr| / max(|S_rc|, |S_cr|) accepted by the Cholesky solver.
pub const SYMMETRY_TOLERANCE: f64 = 1.0e-8;
/// Maximum number of iterative-refinement steps after the Cholesky solve.
const MAX_REFINEMENT_STEPS: usize = 3;

/// Linear solver for J N = F.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LinearSolver {
    /// Banded Cholesky of the symmetrized operator (requires detailed balance), with iterative refinement.
    BandedCholesky,
    /// BiCGSTAB with Jacobi preconditioning on D^-1 J D (no symmetry required).
    BiCgStab { relative_tolerance: f64, max_iterations: usize },
}

/// Source distribution mapped onto the retained states of the operator.
#[derive(Debug, Clone)]
pub struct ProjectedSource {
    /// F on the retained states (in the order of `ChemicalActivationOperator::states`).
    pub on_states: Vec<f64>,
    /// Source fraction that is formed directly in the absorbed grains of each well (directly
    /// stabilized; intermediate steady state only).
    pub absorbed_per_well: Vec<f64>,
}

/// Normalize a source given on the grain grid of every well (source[w][i]) to unit total and map it
/// onto the retained states.
pub fn project_source(op: &ChemicalActivationOperator, source: &[Vec<f64>]) -> Result<ProjectedSource, String> {
    if source.len() != op.index_of.len() {
        return Err(format!("Source given for {} wells, the network has {}.", source.len(), op.index_of.len()));
    }
    let mut total = 0.0;
    for (w, dist) in source.iter().enumerate() {
        if dist.len() != op.index_of[w].len() {
            return Err(format!(
                "Source of well {w} has {} grains, the well has {}.",
                dist.len(),
                op.index_of[w].len()
            ));
        }
        if dist.iter().any(|x| !(*x >= 0.0) || !x.is_finite()) {
            return Err(format!("Source of well {w}: entries must be finite and >= 0."));
        }
        total += dist.iter().sum::<f64>();
    }
    if !(total > 0.0) {
        return Err("The source distribution is zero everywhere.".into());
    }

    // F is the normalized nascent distribution (PO14 eq. 1 and text: "f(E) denotes the normalized
    // nascent distribution"), so that N^s are populations per unit formation rate.
    let mut on_states = vec![0.0; op.dimension()];
    let mut absorbed_per_well = vec![0.0; source.len()];
    for (w, dist) in source.iter().enumerate() {
        for (i, &x) in dist.iter().enumerate() {
            match op.index_of[w][i] {
                Some(s) => on_states[s] = x / total,
                None => absorbed_per_well[w] += x / total,
            }
        }
    }
    Ok(ProjectedSource { on_states, absorbed_per_well })
}

/// Steady-state populations and numerical diagnostics.
#[derive(Debug, Clone)]
pub struct SteadyStateSolution {
    /// N^s on the retained states for R = 1 and the given F.
    pub population: Vec<f64>,
    /// ||F - J N||_2 / ||F||_2, evaluated with the unsymmetrized J.
    pub relative_residual: f64,
    /// max |S_rc - S_cr| / max(|S_rc|, |S_cr|) of the symmetrized operator.
    pub max_relative_asymmetry: f64,
}

/// Solve J N = F on the retained states.
pub fn solve_steady_state(
    op: &ChemicalActivationOperator,
    source_on_states: &[f64],
    solver: &LinearSolver,
) -> Result<SteadyStateSolution, String> {
    let n = op.dimension();
    if source_on_states.len() != n {
        return Err(format!("Source of length {} for {n} states.", source_on_states.len()));
    }
    let f_norm = l2_norm(source_on_states);
    if n == 0 || f_norm == 0.0 {
        // Everything formed below the absorbing barrier: no population above it.
        return Ok(SteadyStateSolution { population: vec![0.0; n], relative_residual: 0.0, max_relative_asymmetry: 0.0 });
    }

    let symmetrized = symmetrize(op);
    let (s_rows, half_log_d, d) = (&symmetrized.rows, &symmetrized.half_log_d, &symmetrized.d);
    let element = |r: usize, c: usize| symmetrized.element(r, c);
    let max_relative_asymmetry = symmetrized.max_relative_asymmetry;

    // Right-hand side of S y = D^-1 F, from logarithms (F/D can exceed the floating-point range only
    // through the factor 1/D).
    let scaled_rhs = |rhs: &[f64]| -> Vec<f64> {
        rhs.iter()
            .zip(half_log_d.iter())
            .map(|(&x, &h)| if x == 0.0 { 0.0 } else { x.signum() * (x.abs().ln() - h).exp() })
            .collect()
    };

    let population = match *solver {
        LinearSolver::BandedCholesky => {
            if max_relative_asymmetry > SYMMETRY_TOLERANCE {
                return Err(format!(
                    "The master-equation operator is not symmetrizable (relative asymmetry {max_relative_asymmetry:e} > \
                     {SYMMETRY_TOLERANCE:e}): the isomerization rates violate detailed balance, \
                     rho_a k_ab = rho_b k_ba. Fix the rates or use the BiCGSTAB solver."
                ));
            }
            let factor = symmetrized.band_matrix().cholesky().map_err(|e| {
                format!(
                    "{e} The symmetrized operator must be positive definite; in the final steady state this \
                     fails when the thermal rate coefficient is negligible compared with the collision frequency \
                     (near-singular J); use the intermediate steady state then."
                )
            })?;

            // N = D y, then iterative refinement on J N = F: N <- N + D S^-1 D^-1 (F - J N).
            let to_population = |y: Vec<f64>| -> Vec<f64> { y.iter().zip(d).map(|(a, b)| a * b).collect() };
            let mut population = to_population(factor.solve(&scaled_rhs(source_on_states))?);
            let mut residual_norm = l2_norm(&residual(op, &population, source_on_states));
            for _ in 0..MAX_REFINEMENT_STEPS {
                let r = residual(op, &population, source_on_states);
                let correction = to_population(factor.solve(&scaled_rhs(&r))?);
                let candidate: Vec<f64> = population.iter().zip(&correction).map(|(a, b)| a + b).collect();
                let candidate_norm = l2_norm(&residual(op, &candidate, source_on_states));
                if !(candidate_norm < residual_norm) {
                    break;
                }
                population = candidate;
                residual_norm = candidate_norm;
            }
            population
        }
        LinearSolver::BiCgStab { relative_tolerance, max_iterations } => {
            struct Scaled<'a> {
                rows: &'a [Vec<(usize, f64)>],
            }
            impl LinearOperator for Scaled<'_> {
                fn dim(&self) -> usize {
                    self.rows.len()
                }
                fn matvec(&mut self, x: &[f64], y: &mut [f64]) -> Result<(), String> {
                    for (yi, row) in y.iter_mut().zip(self.rows) {
                        *yi = row.iter().map(|&(c, v)| v * x[c]).sum();
                    }
                    Ok(())
                }
            }
            let diagonal: Vec<f64> = (0..n).map(|r| element(r, r)).collect();
            let preconditioner = JacobiPreconditioner::from_diagonal(&diagonal)?;
            let (y, _) = solve_bicgstab_left_preconditioned(
                &mut Scaled { rows: &s_rows },
                &preconditioner,
                &scaled_rhs(source_on_states),
                relative_tolerance,
                max_iterations,
            )?;
            y.iter().zip(d.iter()).map(|(a, b)| a * b).collect()
        }
    };

    let relative_residual = l2_norm(&residual(op, &population, source_on_states)) / f_norm;
    Ok(SteadyStateSolution { population, relative_residual, max_relative_asymmetry })
}

/// The symmetrized operator S = D^-1 J D with D = diag(sqrt(f)), S_rc = J_rc sqrt(f_c/f_r)
/// (R19 eq. 5.75), shared by the steady-state solution and the eigenvalue analysis.
pub(crate) struct SymmetrizedOperator {
    /// Rows of S, sorted by column like J.
    pub rows: Vec<Vec<(usize, f64)>>,
    /// ln D_r = (ln f_r - max ln f)/2: D relative to the largest weight (the overall scale of D cancels
    /// in S); logarithms avoid under- and overflow of rho exp(-E/kT).
    pub half_log_d: Vec<f64>,
    /// D_r.
    pub d: Vec<f64>,
    /// max |S_rc - S_cr| / max(|S_rc|, |S_cr|).
    pub max_relative_asymmetry: f64,
}

impl SymmetrizedOperator {
    pub fn element(&self, r: usize, c: usize) -> f64 {
        self.rows[r].binary_search_by_key(&c, |&(col, _)| col).map(|k| self.rows[r][k].1).unwrap_or(0.0)
    }

    /// S as a symmetric band matrix, (S_rc + S_cr)/2 off the diagonal.
    pub fn band_matrix(&self) -> SymmetricBandMatrix {
        let n = self.rows.len();
        let bandwidth = self
            .rows
            .iter()
            .enumerate()
            .flat_map(|(r, row)| row.iter().map(move |&(c, _)| r.abs_diff(c)))
            .max()
            .unwrap_or(0);
        let mut band = SymmetricBandMatrix::zeros(n, bandwidth);
        for (r, row) in self.rows.iter().enumerate() {
            for &(c, v) in row {
                if c == r {
                    band.set_diagonal(r, v);
                } else if c < r {
                    band.set_lower(r, c, 0.5 * (v + self.element(c, r)));
                }
            }
        }
        band
    }

    /// S as a dense symmetric matrix, (S_rc + S_cr)/2 off the diagonal.
    pub fn dense(&self) -> Vec<Vec<f64>> {
        let n = self.rows.len();
        let mut a = vec![vec![0.0; n]; n];
        for (r, row) in self.rows.iter().enumerate() {
            for &(c, v) in row {
                a[r][c] += 0.5 * v;
                a[c][r] += 0.5 * v;
            }
        }
        a
    }
}

pub(crate) fn symmetrize(op: &ChemicalActivationOperator) -> SymmetrizedOperator {
    let log_max = op.log_boltzmann_weight.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let half_log_d: Vec<f64> = op.log_boltzmann_weight.iter().map(|l| 0.5 * (l - log_max)).collect();
    let d = half_log_d.iter().map(|h| h.exp()).collect();
    let rows: Vec<Vec<(usize, f64)>> = op
        .rows
        .iter()
        .enumerate()
        .map(|(r, row)| row.iter().map(|&(c, v)| (c, v * (half_log_d[c] - half_log_d[r]).exp())).collect())
        .collect();
    let mut symmetrized = SymmetrizedOperator { rows, half_log_d, d, max_relative_asymmetry: 0.0 };
    let mut asymmetry: f64 = 0.0;
    for (r, row) in symmetrized.rows.iter().enumerate() {
        for &(c, v) in row {
            if c != r {
                let w = symmetrized.element(c, r);
                asymmetry = asymmetry.max((v - w).abs() / v.abs().max(w.abs()));
            }
        }
    }
    symmetrized.max_relative_asymmetry = asymmetry;
    symmetrized
}

/// F - J N.
fn residual(op: &ChemicalActivationOperator, population: &[f64], source: &[f64]) -> Vec<f64> {
    op.apply(population).iter().zip(source).map(|(jn, f)| f - jn).collect()
}

fn l2_norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::masterequation::chemical_activation_network::{AbsorbingBarrier, ChemicalActivationOptions, CollisionModel, SteadyState};
    use crate::masterequation::chemical_activation_operator::assemble_operator;
    use crate::masterequation::chemical_activation_operator::tests::{all_option_combinations, conditions, two_well_network};

    const BICGSTAB: LinearSolver = LinearSolver::BiCgStab { relative_tolerance: 1e-13, max_iterations: 20_000 };

    /// Nascent distribution in well A between grains 320 and 399 (a hot, chemically activated source).
    fn hot_source(op: &ChemicalActivationOperator) -> ProjectedSource {
        let source = vec![
            (0..400).map(|i| if i >= 320 { (-((i - 320) as f64) / 15.0).exp() } else { 0.0 }).collect(),
            vec![0.0; 460],
        ];
        project_source(op, &source).unwrap()
    }

    fn residual(op: &ChemicalActivationOperator, n: &[f64], f: &[f64]) -> f64 {
        let jn = op.apply(n);
        let num: f64 = jn.iter().zip(f).map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt();
        num / f.iter().map(|x| x * x).sum::<f64>().sqrt()
    }

    #[test]
    fn steady_state_satisfies_j_n_equals_f_for_both_solvers_and_all_options() {
        let network = two_well_network();
        for options in all_option_combinations() {
            let op = assemble_operator(&network, &conditions(), &options).unwrap();
            let f = hot_source(&op);
            for solver in [LinearSolver::BandedCholesky, BICGSTAB] {
                let solution = solve_steady_state(&op, &f.on_states, &solver).unwrap();
                let r = residual(&op, &solution.population, &f.on_states);
                assert!(r < 1e-10, "{options:?} {solver:?}: residual {r:e}");
                assert!((solution.relative_residual - r).abs() < 1e-12);
                let largest = solution.population.iter().cloned().fold(0.0, f64::max);
                assert!(solution.population.iter().all(|&x| x > -1e-10 * largest), "{options:?} {solver:?}: negative population");
            }
        }
    }

    #[test]
    fn banded_cholesky_and_bicgstab_give_the_same_populations() {
        let network = two_well_network();
        let options = ChemicalActivationOptions {
            collision_model: CollisionModel::ExponentialDown { cutoff_in_mean_down: 10.0 },
            steady_state: SteadyState::Intermediate { barrier: AbsorbingBarrier::default() },
        };
        let op = assemble_operator(&network, &conditions(), &options).unwrap();
        let f = hot_source(&op);
        let a = solve_steady_state(&op, &f.on_states, &LinearSolver::BandedCholesky).unwrap();
        let b = solve_steady_state(&op, &f.on_states, &BICGSTAB).unwrap();
        let scale = a.population.iter().cloned().fold(0.0, f64::max);
        for (x, y) in a.population.iter().zip(&b.population) {
            assert!((x - y).abs() < 1e-8 * scale, "{x:e} vs {y:e}");
        }
    }

    #[test]
    fn cholesky_refuses_an_operator_without_detailed_balance() {
        let mut network = two_well_network();
        for k in network.wells[1].channels[1].rate_constant_s_inv.iter_mut() {
            *k *= 2.0; // B -> A no longer the microscopic reverse of A -> B
        }
        let options = all_option_combinations().remove(1);
        let op = assemble_operator(&network, &conditions(), &options).unwrap();
        let f = hot_source(&op);
        assert!(solve_steady_state(&op, &f.on_states, &LinearSolver::BandedCholesky).is_err());
        let solution = solve_steady_state(&op, &f.on_states, &BICGSTAB).unwrap();
        assert!(solution.max_relative_asymmetry > 0.1);
        assert!(residual(&op, &solution.population, &f.on_states) < 1e-10);
    }

    #[test]
    fn source_in_absorbed_grains_counts_as_directly_stabilized() {
        let network = two_well_network();
        let options = ChemicalActivationOptions {
            collision_model: CollisionModel::Stepladder,
            steady_state: SteadyState::Intermediate { barrier: AbsorbingBarrier::AtGrains(vec![100, 100]) },
        };
        let op = assemble_operator(&network, &conditions(), &options).unwrap();
        let mut a = vec![0.0; 400];
        a[50] = 1.0; // absorbed
        a[350] = 3.0; // retained
        let projected = project_source(&op, &[a, vec![0.0; 460]]).unwrap();
        assert!((projected.absorbed_per_well[0] - 0.25).abs() < 1e-15);
        assert_eq!(projected.absorbed_per_well[1], 0.0);
        assert!((projected.on_states.iter().sum::<f64>() - 0.75).abs() < 1e-15);
        assert_eq!(projected.on_states[op.index_of[0][350].unwrap()], 0.75);
    }

    #[test]
    fn project_source_rejects_negative_or_empty_distributions() {
        let network = two_well_network();
        let op = assemble_operator(&network, &conditions(), &all_option_combinations().remove(0)).unwrap();
        let mut a = vec![0.0; 400];
        assert!(project_source(&op, &[a.clone(), vec![0.0; 460]]).is_err());
        a[350] = -1.0;
        assert!(project_source(&op, &[a, vec![0.0; 460]]).is_err());
        assert!(project_source(&op, &[vec![1.0; 399], vec![0.0; 460]]).is_err());
    }
}
