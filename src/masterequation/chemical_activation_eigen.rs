//! Eigenvalue analysis of the master equation (Olzmann's solution method).
//!
//! The relaxation matrix J of PO14 eq. 2 is similar to the symmetric matrix S = D^-1 J D, D = diag(sqrt(f))
//! (R19 eqs. 5.74-5.77; "after symmetrization", PO14 p. 235), whose eigenvalues lambda_i > 0 (GO10 text
//! before eq. 12) and orthonormal eigenvectors u_i give the eigenvectors E_i = D u_i of J. With these:
//!
//!   thermal (unimolecular) rate coefficient    k^th = lambda_1, the lowest eigenvalue of J      (GO10 eq. 12)
//!   channel rate coefficients    k_j^th = sum_i k_j(E_i) Ñ_i with Ñ the normalized eigenvector of lambda_1,
//!                                "analogous to eqn (9) but with Ñs = Ñs_th" (GO10 after eq. 12); with the
//!                                bimolecular sink, sum_j k_j^th + k_c[D] = lambda_1 (column sums of J)
//!   time-dependent populations for a source R F switched on at t = 0, N(0) = 0:
//!                                N(t) = R sum_i e_i (1 - exp(-lambda_i t))/lambda_i E_i,  F = sum_i e_i E_i
//!                                                                                    (PO14 eqs. 3-4)
//!   validity of the steady-state (D/S) picture: the absorbing-barrier result is reproduced by the final
//!   steady state with a bimolecular sink for 0.01 omega > k_c[D] > 10 k_uni, and the intermediate steady
//!   state exists for (0.1 lambda_F)^-1 < t < (10 k_uni)^-1, lambda_F being the eigenvalue "that
//!   corresponds to the eigenvector most closely resembling the initial distribution f(E)" (O02, discussion
//!   of Fig. 2, after Schranz, Nordholm, Chem. Phys. 87, 163 (1984)).
//!
//! Three solvers (`EigenSolver`):
//!   - inverse iteration with the banded Cholesky factor of S + sigma I (default; sigma = n eps max S_ii,
//!     `cholesky_safety_shift`, keeps the factor computable when lambda_1 is below the double-precision
//!     resolution; same eigenvectors): the lowest eigenpairs only, O(n bw^2), the most accurate thermal
//!     eigenvector;
//!   - full decomposition (Householder reduction + QL with implicit shifts, the EISPACK tred2/tql2 route
//!     used by Olzmann, PO14 ref. 39; O02 ref. 34): all eigenpairs, O(n^3), needed for N(t);
//!   - the same full decomposition by LAPACK DSYEVD (Householder reduction + divide and conquer,
//!     `numeric/lapack_interface.rs`): all eigenpairs, much faster for large n; requires the `openblas`
//!     build feature (an error otherwise).
//! Both require exact detailed balance (a symmetrizable J).
//!
//! Sum rule. The column sums of J are the loss rates of the grains (product channels and sink; collisions
//! and isomerization conserve population), so for the normalized thermal eigenvector
//!   lambda_1 = 1^T J Ñ = sum_j k_j^th + k_c[D]                                       (GO10 eq. 12)
//! where the right-hand side is GO10's "averaging procedure analogous to eqn (9) but with Ñs = Ñs^th".
//! The right-hand side is a sum of positive terms, free of cancellation. lambda_1 from the eigen-solver is
//! an eigenvalue of S + dS with ||dS|| of order eps ||S|| (backward error of the Cholesky factorization,
//! Higham, Accuracy and Stability of Numerical Algorithms, 2nd ed. (2002), ch. 10, and of the Householder/QL
//! reduction, Wilkinson, The Algebraic Eigenvalue Problem (1965)), hence carries an absolute error of order
//! eps ||S|| (Weyl, Math. Ann. 71, 441 (1912)); the eigenvector is perturbed only by about eps ||S|| /
//! (lambda_2 - lambda_1) (Davis, Kahan, SIAM J. Numer. Anal. 7, 1 (1970)). The reported unimolecular rate
//! coefficient is therefore the right-hand side, k_uni = sum_j k_j^th + k_c[D] (GO10's averaging
//! procedure), with lambda_1 reported beside it. Their relative difference |lambda_1 - k_uni| / k_uni
//! measures whether the eigenpair is resolved in double precision; above the tolerance it produces a
//! warning (the result is kept).
//!
//! The analysis uses the operator without absorbing barrier (the final-steady-state operator), as GO10
//! eq. 12; no absorbing barrier is needed for thermal or pressure-dependent unimolecular rate
//! coefficients.
//!
//! References: see `chemical_activation_network.rs`.

use crate::numeric::lapack_interface::symmetric_eigen_lapack;
use crate::numeric::symmetric_eigen::{cholesky_safety_shift, lowest_eigenpairs_banded, symmetric_eigen};

use super::chemical_activation_network::{ChannelDestination, ChemicalActivationNetwork};
use super::chemical_activation_operator::ChemicalActivationOperator;
use super::chemical_activation_steady_state::{symmetrize, SYMMETRY_TOLERANCE};

/// Eigen-solver of the analysis.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum EigenSolver {
    /// Lowest eigenpairs by inverse iteration with the banded Cholesky factor (default).
    #[default]
    InverseIteration,
    /// All eigenpairs by Householder reduction and the QL algorithm.
    FullDecomposition,
    /// All eigenpairs by LAPACK DSYEVD (Householder reduction and divide and conquer).
    FullDecompositionLapack,
}

/// All eigenpairs of the symmetrized operator by one of the full decompositions.
fn full_decomposition(dense: &[Vec<f64>], solver: EigenSolver) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
    match solver {
        EigenSolver::FullDecomposition => symmetric_eigen(dense),
        EigenSolver::FullDecompositionLapack => symmetric_eigen_lapack(dense),
        EigenSolver::InverseIteration => Err(
            "Eigenvalue analysis: inverse iteration gives only the lowest eigenpairs; the time-dependent \
             solution needs all of them (full decomposition)."
                .into(),
        ),
    }
}

/// Rate coefficient of one channel from the thermal (lowest) eigenvector.
#[derive(Debug, Clone)]
pub struct ChannelThermalRate {
    pub well: usize,
    pub channel: usize,
    pub name: String,
    pub destination: ChannelDestination,
    /// k_j^th = sum_i k_j(E_i) Ñ_i (s-1), Ñ normalized over all states of the network.
    pub thermal_rate_s_inv: f64,
    /// High-pressure limit k_j^inf = sum_i k_j(E_i) f_i / sum_i f_i over the well, f = rho exp(-E/kT).
    pub high_pressure_rate_s_inv: f64,
}

/// Thermal rate coefficients of a network at one temperature and pressure.
#[derive(Debug, Clone)]
pub struct ThermalRateCoefficients {
    /// k_uni = sum_j k_j^th + k_c[D] over the product channels and sinks (s-1): the specific rate
    /// coefficients averaged over the normalized thermal eigenvector (GO10, text after eq. 12). The
    /// reported unimolecular (thermal) rate coefficient.
    pub k_uni_s_inv: f64,
    /// Lowest eigenvalue of J (s-1), equal to k_uni in exact arithmetic (GO10 eq. 12).
    pub lambda_1_s_inv: f64,
    /// |lambda_1 - k_uni| / k_uni.
    pub sum_rule_relative_deviation: f64,
    /// Set when the sum-rule deviation exceeds the tolerance: explanation and the numbers.
    pub warning: Option<String>,
    /// Shift sigma (s-1) of the factorization S + sigma I in the inverse iteration (`cholesky_safety_shift`);
    /// 0 for the full decompositions.
    pub inverse_iteration_shift_s_inv: f64,
    /// Double-precision floor eps max|S_ii| (s-1): the order of the absolute rounding error of lambda_1.
    /// A lambda_1 at or below it is rounding noise (it may even come out negative).
    pub precision_floor_s_inv: f64,
    /// lambda_2 (s-1): lambda_2/k_uni measures the separation of the thermal decay from relaxation.
    /// NaN if inverse iteration cannot resolve it (separation beyond about 1e16).
    pub lambda_2_s_inv: f64,
    pub channels: Vec<ChannelThermalRate>,
    /// k_c[D] times the population fraction of each well.
    pub sink_rates_s_inv: Vec<f64>,
    /// Population fraction of each well in the thermal eigenvector.
    pub population_fractions: Vec<f64>,
    /// Normalized thermal distribution of each well on its grid.
    pub distributions: Vec<Vec<f64>>,
}

/// Default bound on the relative sum-rule deviation |lambda_1 - k_uni| / k_uni above which a warning is given.
pub const DEFAULT_SUM_RULE_TOLERANCE: f64 = 1.5e-2;

/// Thermal rate coefficients from the lowest eigenpair of J (GO10 eq. 12): k_uni as the eigenvector
/// average, lambda_1, and a warning if they differ by more than `sum_rule_tolerance` (relative). An error
/// only if no eigenpair can be computed or the thermal eigenvector gives no positive loss rate.
pub fn thermal_rate_coefficients(
    network: &ChemicalActivationNetwork,
    op: &ChemicalActivationOperator,
    solver: EigenSolver,
    sum_rule_tolerance: f64,
) -> Result<ThermalRateCoefficients, String> {
    let symmetrized = symmetrize(op);
    require_symmetry(symmetrized.max_relative_asymmetry)?;
    let n = op.dimension();
    if n == 0 {
        return Err("Eigenvalue analysis: the operator has no states.".into());
    }
    // Double-precision floor of the eigenvalues: eps max|S_ii| (Weyl's inequality with a backward error of
    // order eps ||S||, see the module documentation).
    let precision_floor = f64::EPSILON * symmetrized.band_matrix().max_abs_diagonal();
    // Lowest two eigenpairs of S.
    let mut shift = 0.0;
    let (lambda_1, u_1, lambda_2) = match solver {
        EigenSolver::InverseIteration => {
            // Start from sqrt(f), the symmetrized Boltzmann distribution, close to the thermal eigenvector.
            let band = symmetrized.band_matrix();
            // The factor of S + sigma I, sigma = n eps max S_ii: same eigenvectors, and the factor exists
            // also where lambda_1 of S lies below the double-precision resolution (`cholesky_safety_shift`).
            shift = cholesky_safety_shift(&band);
            let first = lowest_eigenpairs_banded(&band, 1, &symmetrized.d, shift, 1e-11, 5000).map_err(|e| {
                format!(
                    "{e} (Cholesky factor of S + sigma I with the safety shift sigma = {shift:e} s-1.) What to do: \
                     use a full decomposition, which does not need the factor (EigenSolver::FullDecompositionLapack, \
                     `--eigen-solver lapack`, or FullDecomposition); it gives k_uni from the thermal eigenvector."
                )
            })?;
            // lambda_2 is only a diagnostic of the separation lambda_2/lambda_1. The deflated iteration
            // cannot resolve it when lambda_2/lambda_1 exceeds the floating-point range of the
            // re-orthogonalization (about 1e16); it is then reported as NaN ("not resolved") and lambda_1
            // is kept.
            let lambda_2 = if n >= 2 {
                lowest_eigenpairs_banded(&band, 2, &symmetrized.d, shift, 1e-11, 5000).map_or(f64::NAN, |p| p[1].0)
            } else {
                f64::NAN
            };
            (first[0].0, first[0].1.clone(), lambda_2)
        }
        EigenSolver::FullDecomposition | EigenSolver::FullDecompositionLapack => {
            let (values, vectors) = full_decomposition(&symmetrized.dense(), solver)?;
            (values[0], vectors[0].clone(), values.get(1).copied().unwrap_or(f64::NAN))
        }
    };

    // Thermal eigenvector of J, E_1 = D u_1, normalized to unit population.
    let mut population: Vec<f64> = u_1.iter().zip(&symmetrized.d).map(|(u, d)| u * d).collect();
    let total: f64 = population.iter().sum();
    population.iter_mut().for_each(|x| *x /= total);

    let n_wells = network.wells.len();
    let mut per_well: Vec<Vec<f64>> = network.wells.iter().map(|w| vec![0.0; w.grain_count()]).collect();
    for (s, &(w, i)) in op.states.iter().enumerate() {
        per_well[w][i] = population[s];
    }
    let population_fractions: Vec<f64> = per_well.iter().map(|p| p.iter().sum()).collect();
    let kt = op.kt_cm1;
    let mut channels = Vec::new();
    for (w, well) in network.wells.iter().enumerate() {
        // Boltzmann weights of the well relative to their maximum.
        let log_f: Vec<f64> = (0..well.grain_count())
            .map(|i| well.density_of_states[i].ln() - network.absolute_energy_cm1(w, i) / kt)
            .collect();
        let max = log_f.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let f: Vec<f64> = log_f.iter().map(|l| (l - max).exp()).collect();
        let f_total: f64 = f.iter().sum();
        for (c, channel) in well.channels.iter().enumerate() {
            let k = &channel.rate_constant_s_inv;
            channels.push(ChannelThermalRate {
                well: w,
                channel: c,
                name: channel.name.clone(),
                destination: channel.destination.clone(),
                thermal_rate_s_inv: k.iter().zip(&per_well[w]).map(|(k, p)| k * p).sum(),
                high_pressure_rate_s_inv: k.iter().zip(&f).map(|(k, f)| k * f).sum::<f64>() / f_total,
            });
        }
    }
    let sink_rates_s_inv: Vec<f64> =
        (0..n_wells).map(|w| network.wells[w].bimolecular_sink_s_inv * population_fractions[w]).collect();
    let k_uni = channels
        .iter()
        .filter(|c| matches!(c.destination, ChannelDestination::Products { .. }))
        .map(|c| c.thermal_rate_s_inv)
        .sum::<f64>()
        + sink_rates_s_inv.iter().sum::<f64>();
    if !(k_uni > 0.0 && k_uni.is_finite()) {
        return Err(format!(
            "Eigenvalue analysis: the thermal eigenvector gives no positive loss rate (k_uni = {k_uni:e} s-1, \
             lambda_1 = {lambda_1:e} s-1): no product channel or sink is reached from the network, or the \
             eigenvector is not resolved."
        ));
    }
    let sum_rule_relative_deviation = (lambda_1 - k_uni).abs() / k_uni;
    let warning = (!(sum_rule_relative_deviation <= sum_rule_tolerance)).then(|| {
        let separation = lambda_2 / k_uni;
        let sign = if lambda_1 > 0.0 {
            String::new()
        } else {
            format!(
                " The lowest eigenvalue is not positive, which is impossible for J when population can leave the \
                 network (all eigenvalues are positive, GO10 text before eq. 12): lambda_1 = {lambda_1:e} s-1 is \
                 rounding noise below the double-precision floor eps*max(S_ii) = {precision_floor:e} s-1. This is \
                 not a merging of eigenvalues (thermal decay indistinguishable from relaxation): the separation \
                 lambda_2/k_uni = {separation:e} is large; merging would show as a small lambda_2/k_uni."
            )
        };
        format!(
            "sum rule lambda_1 = k_uni violated: lambda_1 = {lambda_1:e} s-1, k_uni = {k_uni:e} s-1, relative \
             deviation {sum_rule_relative_deviation:.3e} > {sum_rule_tolerance:e}.{sign} k_uni is the average of the \
             specific rate coefficients over the normalized thermal eigenvector (Gonzalez-Garcia, Olzmann, Phys. \
             Chem. Chem. Phys. 12, 12290 (2010), text after eq. 12), equal to the lowest eigenvalue lambda_1 (eq. \
             12) in exact arithmetic. Here lambda_1 is not resolved in double precision (it lies far below the \
             collision frequency and the largest k(E)); the eigenvector, and with it k_uni, is much less sensitive. \
             What to do: use k_uni, and confirm it with the inverse iteration (the most accurate solver for the \
             thermal eigenvector); if the solvers give different k_uni, the condition is beyond double precision."
        )
    });
    let distributions = per_well
        .iter()
        .zip(&population_fractions)
        .map(|(p, &total)| if total > 0.0 { p.iter().map(|x| x / total).collect() } else { p.clone() })
        .collect();
    Ok(ThermalRateCoefficients {
        k_uni_s_inv: k_uni,
        lambda_1_s_inv: lambda_1,
        sum_rule_relative_deviation,
        warning,
        inverse_iteration_shift_s_inv: shift,
        precision_floor_s_inv: precision_floor,
        lambda_2_s_inv: lambda_2,
        channels,
        sink_rates_s_inv,
        population_fractions,
        distributions,
    })
}

fn require_symmetry(asymmetry: f64) -> Result<(), String> {
    if asymmetry > SYMMETRY_TOLERANCE {
        return Err(format!(
            "Eigenvalue analysis: the operator is not symmetrizable (relative asymmetry {asymmetry:e}); the \
             isomerization rates violate detailed balance."
        ));
    }
    Ok(())
}

/// All eigenpairs of J for the time-dependent solution.
#[derive(Debug, Clone)]
pub struct EigenSystem {
    /// lambda_i ascending (s-1).
    pub eigenvalues: Vec<f64>,
    /// Orthonormal eigenvectors u_i of S; E_i = D u_i.
    symmetric_vectors: Vec<Vec<f64>>,
    /// D_r and ln D_r.
    d: Vec<f64>,
    half_log_d: Vec<f64>,
}

impl EigenSystem {
    /// Full decomposition of the symmetrized J (`FullDecomposition` or `FullDecompositionLapack`).
    pub fn new(op: &ChemicalActivationOperator, solver: EigenSolver) -> Result<Self, String> {
        let symmetrized = symmetrize(op);
        require_symmetry(symmetrized.max_relative_asymmetry)?;
        let (eigenvalues, symmetric_vectors) = full_decomposition(&symmetrized.dense(), solver)?;
        Ok(Self { eigenvalues, symmetric_vectors, d: symmetrized.d, half_log_d: symmetrized.half_log_d })
    }

    /// Expansion coefficients e_i of F = sum_i e_i E_i: e_i = u_i . D^-1 F.
    fn expansion(&self, source_on_states: &[f64]) -> Vec<f64> {
        let scaled: Vec<f64> = source_on_states
            .iter()
            .zip(&self.half_log_d)
            .map(|(&x, &h)| if x == 0.0 { 0.0 } else { x.signum() * (x.abs().ln() - h).exp() })
            .collect();
        self.symmetric_vectors.iter().map(|u| u.iter().zip(&scaled).map(|(a, b)| a * b).sum()).collect()
    }

    /// N(t) for the source R F (F on the states of the operator) switched on at t = 0 (PO14 eqs. 3-4).
    pub fn time_dependent_population(&self, source_on_states: &[f64], formation_rate: f64, time_s: f64) -> Vec<f64> {
        let e = self.expansion(source_on_states);
        let n = self.d.len();
        let mut y = vec![0.0; n];
        for (k, u) in self.symmetric_vectors.iter().enumerate() {
            // (1 - exp(-lambda t))/lambda, written with expm1 for small lambda t.
            let lambda = self.eigenvalues[k];
            let weight = formation_rate * e[k] * (-(-lambda * time_s).exp_m1()) / lambda;
            y.iter_mut().zip(u).for_each(|(a, b)| *a += weight * b);
        }
        y.iter().zip(&self.d).map(|(a, d)| a * d).collect()
    }

    /// lambda_F: the eigenvalue whose eigenvector has the largest weight in F (O02).
    pub fn lambda_f(&self, source_on_states: &[f64]) -> f64 {
        // The eigenvector with the largest expansion coefficient in the orthonormal (symmetrized) basis.
        let e = self.expansion(source_on_states);
        let k = (0..e.len()).max_by(|&a, &b| e[a].abs().partial_cmp(&e[b].abs()).unwrap()).unwrap_or(0);
        self.eigenvalues[k]
    }
}

/// Validity window of the steady-state (decomposition/stabilization) picture (O02, Fig. 2 discussion).
#[derive(Debug, Clone, Copy)]
pub struct SteadyStateWindow {
    /// 0.01 omega > k_c[D] > 10 k_uni: the final steady state with the sink reproduces the absorbing-barrier
    /// (intermediate steady-state) branching.
    pub sink_within_window: bool,
    /// (0.1 lambda_F)^-1 and (10 k_uni)^-1 (s): the time window of the intermediate steady state.
    pub intermediate_steady_state_from_s: f64,
    pub intermediate_steady_state_until_s: f64,
}

pub fn steady_state_window(collision_frequency_s_inv: f64, sink_s_inv: f64, k_uni_s_inv: f64, lambda_f_s_inv: f64) -> SteadyStateWindow {
    SteadyStateWindow {
        sink_within_window: 0.01 * collision_frequency_s_inv > sink_s_inv && sink_s_inv > 10.0 * k_uni_s_inv,
        intermediate_steady_state_from_s: 1.0 / (0.1 * lambda_f_s_inv),
        intermediate_steady_state_until_s: 1.0 / (10.0 * k_uni_s_inv),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::KB_CM;
    use crate::masterequation::chemical_activation_network::tests::test_well;
    use crate::masterequation::chemical_activation_network::{
        ChemicalActivationOptions, CollisionModel, Conditions, SteadyState,
    };
    use crate::masterequation::chemical_activation_operator::assemble_operator;
    use crate::masterequation::chemical_activation_operator::tests::{conditions, two_well_network};
    use crate::masterequation::chemical_activation_steady_state::{project_source, solve_steady_state, LinearSolver};

    fn final_options(collision_model: CollisionModel) -> ChemicalActivationOptions {
        ChemicalActivationOptions { collision_model, steady_state: SteadyState::Final }
    }

    const EXPONENTIAL: CollisionModel = CollisionModel::ExponentialDown { cutoff_in_mean_down: 10.0 };

    #[test]
    fn both_solvers_give_the_lowest_eigenpairs_of_j() {
        let network = two_well_network();
        let op = assemble_operator(&network, &conditions(), &final_options(EXPONENTIAL)).unwrap();
        let ii = thermal_rate_coefficients(&network, &op, EigenSolver::InverseIteration, DEFAULT_SUM_RULE_TOLERANCE).unwrap();
        let full = thermal_rate_coefficients(&network, &op, EigenSolver::FullDecomposition, DEFAULT_SUM_RULE_TOLERANCE).unwrap();
        assert!((ii.lambda_1_s_inv / full.lambda_1_s_inv - 1.0).abs() < 1e-8, "{} vs {}", ii.lambda_1_s_inv, full.lambda_1_s_inv);
        assert!((ii.lambda_2_s_inv / full.lambda_2_s_inv - 1.0).abs() < 1e-6, "{} vs {}", ii.lambda_2_s_inv, full.lambda_2_s_inv);
        assert!(ii.lambda_2_s_inv > ii.lambda_1_s_inv);
        // J E1 = lambda_1 E1 for the population vector of the thermal eigenvector.
        let e1: Vec<f64> = op.states.iter().map(|&(w, i)| ii.distributions[w][i] * ii.population_fractions[w]).collect();
        let je1 = op.apply(&e1);
        let scale = je1.iter().map(|x| x.abs()).fold(0.0, f64::max);
        for (a, b) in je1.iter().zip(&e1) {
            assert!((a - ii.lambda_1_s_inv * b).abs() < 1e-8 * scale);
        }
    }

    #[cfg(feature = "openblas")]
    #[test]
    fn the_lapack_decomposition_agrees_with_the_other_solvers() {
        let network = two_well_network();
        let op = assemble_operator(&network, &conditions(), &final_options(EXPONENTIAL)).unwrap();
        let ii = thermal_rate_coefficients(&network, &op, EigenSolver::InverseIteration, DEFAULT_SUM_RULE_TOLERANCE).unwrap();
        let lapack =
            thermal_rate_coefficients(&network, &op, EigenSolver::FullDecompositionLapack, DEFAULT_SUM_RULE_TOLERANCE).unwrap();
        assert!((lapack.lambda_1_s_inv / ii.lambda_1_s_inv - 1.0).abs() < 1e-8);
        assert!((lapack.lambda_2_s_inv / ii.lambda_2_s_inv - 1.0).abs() < 1e-6);
        // The full spectrum for N(t) from either decomposition.
        let native = EigenSystem::new(&op, EigenSolver::FullDecomposition).unwrap();
        let system = EigenSystem::new(&op, EigenSolver::FullDecompositionLapack).unwrap();
        let largest = native.eigenvalues.last().unwrap().abs();
        for (a, b) in system.eigenvalues.iter().zip(&native.eigenvalues) {
            assert!((a - b).abs() < 1e-10 * largest);
        }
    }

    #[test]
    fn the_time_dependent_solution_needs_a_full_decomposition() {
        let network = two_well_network();
        let op = assemble_operator(&network, &conditions(), &final_options(EXPONENTIAL)).unwrap();
        assert!(EigenSystem::new(&op, EigenSolver::InverseIteration).is_err());
    }

    #[test]
    fn channel_rates_and_sink_add_up_to_lambda_1() {
        // GO10 eq. 12: k^th = sum_j k_j^th (+ k_c[D] with a sink) = lambda_1.
        let network = two_well_network();
        let op = assemble_operator(&network, &conditions(), &final_options(CollisionModel::Stepladder)).unwrap();
        let th = thermal_rate_coefficients(&network, &op, EigenSolver::InverseIteration, DEFAULT_SUM_RULE_TOLERANCE).unwrap();
        let products: f64 = th
            .channels
            .iter()
            .filter(|c| matches!(c.destination, ChannelDestination::Products { .. }))
            .map(|c| c.thermal_rate_s_inv)
            .sum();
        let sinks: f64 = th.sink_rates_s_inv.iter().sum();
        assert!(((products + sinks) / th.lambda_1_s_inv - 1.0).abs() < 1e-8, "{} + {} vs {}", products, sinks, th.lambda_1_s_inv);
        assert!((th.population_fractions.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn k_uni_is_the_eigenvector_average_and_the_sum_rule_deviation_is_reported() {
        // GO10: k^th from the "averaging procedure analogous to eqn (9) but with Ñs = Ñs^th"; equal to lambda_1
        // (eq. 12). The reported k_uni is the average; lambda_1 and the relative difference are reported with it.
        let network = two_well_network();
        let op = assemble_operator(&network, &conditions(), &final_options(EXPONENTIAL)).unwrap();
        let th = thermal_rate_coefficients(&network, &op, EigenSolver::InverseIteration, DEFAULT_SUM_RULE_TOLERANCE).unwrap();
        let products: f64 = th
            .channels
            .iter()
            .filter(|c| matches!(c.destination, ChannelDestination::Products { .. }))
            .map(|c| c.thermal_rate_s_inv)
            .sum();
        let sinks: f64 = th.sink_rates_s_inv.iter().sum();
        assert_eq!(th.k_uni_s_inv, products + sinks);
        assert!((th.k_uni_s_inv / th.lambda_1_s_inv - 1.0).abs() < 1e-8);
        assert_eq!(th.sum_rule_relative_deviation, (th.lambda_1_s_inv - th.k_uni_s_inv).abs() / th.k_uni_s_inv);
        assert!(th.warning.is_none());
    }

    #[test]
    fn a_violated_sum_rule_gives_a_warning_and_no_error() {
        // A deep well (threshold 6000 cm-1, exponential down) at 200 K: lambda_1 differs from the eigenvector
        // average by about 6e-5, below the default tolerance and above 1e-5.
        let network = ChemicalActivationNetwork { grain_width_cm1: 10.0, wells: vec![test_well("A", 800, 0, 600)] };
        let cold = Conditions { temperature_kelvin: 200.0, pressure_torr: 10.0 };
        let op = assemble_operator(&network, &cold, &final_options(EXPONENTIAL)).unwrap();
        let quiet = thermal_rate_coefficients(&network, &op, EigenSolver::InverseIteration, DEFAULT_SUM_RULE_TOLERANCE).unwrap();
        assert!(quiet.sum_rule_relative_deviation > 1e-5 && quiet.sum_rule_relative_deviation < DEFAULT_SUM_RULE_TOLERANCE);
        assert!(quiet.warning.is_none());
        let warned = thermal_rate_coefficients(&network, &op, EigenSolver::InverseIteration, 1e-5).unwrap();
        let warning = warned.warning.as_ref().expect("a warning above the tolerance");
        assert!(warning.contains("sum rule") && warning.contains("Phys. Chem. Chem. Phys. 12, 12290"), "{warning}");
        assert_eq!(warned.k_uni_s_inv, quiet.k_uni_s_inv);
    }

    #[test]
    fn a_non_positive_lambda_1_is_a_warning_when_the_eigenvector_average_is_positive() {
        // At 150 K the lowest eigenvalue of the full decomposition is below the double-precision resolution
        // and comes out negative; the thermal eigenvector still gives a positive loss rate, the same as the
        // shifted inverse iteration.
        let network = ChemicalActivationNetwork { grain_width_cm1: 10.0, wells: vec![test_well("A", 800, 0, 600)] };
        let cold = Conditions { temperature_kelvin: 150.0, pressure_torr: 10.0 };
        let op = assemble_operator(&network, &cold, &final_options(EXPONENTIAL)).unwrap();
        let th = thermal_rate_coefficients(&network, &op, EigenSolver::FullDecomposition, DEFAULT_SUM_RULE_TOLERANCE).unwrap();
        assert!(th.lambda_1_s_inv <= 0.0);
        assert!(th.k_uni_s_inv > 0.0);
        let ii = thermal_rate_coefficients(&network, &op, EigenSolver::InverseIteration, DEFAULT_SUM_RULE_TOLERANCE).unwrap();
        assert!((ii.k_uni_s_inv / th.k_uni_s_inv - 1.0).abs() < 1e-5, "{} vs {}", ii.k_uni_s_inv, th.k_uni_s_inv);
        // The warning says that a non-positive lambda_1 is numerical (GO10: all eigenvalues are positive) and
        // what to do.
        let warning = th.warning.as_ref().expect("a warning");
        assert!(warning.contains("not positive") && warning.contains("inverse iteration"), "{warning}");
        // It names the cause: rounding noise below the double-precision floor (printed), not a merging of
        // eigenvalues (the separation lambda_2/k_uni is printed).
        assert!(th.precision_floor_s_inv > 0.0 && th.lambda_1_s_inv.abs() < 100.0 * th.precision_floor_s_inv);
        assert!(
            warning.contains("rounding noise below the double-precision floor")
                && warning.contains(&format!("{:e}", th.precision_floor_s_inv))
                && warning.contains("not a merging of eigenvalues"),
            "{warning}"
        );
    }

    #[test]
    fn the_shifted_factorization_gives_k_uni_where_the_plain_cholesky_factor_does_not_exist() {
        // At 125 K the Cholesky factor of S itself does not exist in double precision (non-positive pivot);
        // with the safety shift sigma = n eps max S_ii the inverse iteration finds the thermal eigenvector.
        let network = ChemicalActivationNetwork { grain_width_cm1: 10.0, wells: vec![test_well("A", 800, 0, 600)] };
        let options = final_options(EXPONENTIAL);
        let cold = Conditions { temperature_kelvin: 125.0, pressure_torr: 10.0 };
        let op = assemble_operator(&network, &cold, &options).unwrap();
        assert!(symmetrize(&op).band_matrix().cholesky().is_err(), "the unshifted factor exists after all");
        let th = thermal_rate_coefficients(&network, &op, EigenSolver::InverseIteration, DEFAULT_SUM_RULE_TOLERANCE).unwrap();
        let band = symmetrize(&op).band_matrix();
        assert_eq!(th.inverse_iteration_shift_s_inv, cholesky_safety_shift(&band));
        assert!(th.k_uni_s_inv > 0.0);
        // The same k_uni as the full decomposition, which needs no factor.
        let full = thermal_rate_coefficients(&network, &op, EigenSolver::FullDecomposition, DEFAULT_SUM_RULE_TOLERANCE).unwrap();
        assert!((th.k_uni_s_inv / full.k_uni_s_inv - 1.0).abs() < 1e-5, "{} vs {}", th.k_uni_s_inv, full.k_uni_s_inv);
        // Falls with temperature: below the value at 150 K, where the unshifted factor exists.
        let warmer = Conditions { temperature_kelvin: 150.0, pressure_torr: 10.0 };
        let op_150 = assemble_operator(&network, &warmer, &options).unwrap();
        let th_150 =
            thermal_rate_coefficients(&network, &op_150, EigenSolver::InverseIteration, DEFAULT_SUM_RULE_TOLERANCE).unwrap();
        assert!(th.k_uni_s_inv < th_150.k_uni_s_inv);
    }

    #[test]
    fn the_shift_does_not_change_a_resolved_result() {
        // Where lambda_1 is resolved, k_uni and lambda_1 agree with the full decomposition (no shift).
        let network = two_well_network();
        let op = assemble_operator(&network, &conditions(), &final_options(EXPONENTIAL)).unwrap();
        let ii = thermal_rate_coefficients(&network, &op, EigenSolver::InverseIteration, DEFAULT_SUM_RULE_TOLERANCE).unwrap();
        let full = thermal_rate_coefficients(&network, &op, EigenSolver::FullDecomposition, DEFAULT_SUM_RULE_TOLERANCE).unwrap();
        assert!(ii.inverse_iteration_shift_s_inv > 0.0);
        assert_eq!(full.inverse_iteration_shift_s_inv, 0.0);
        assert!((ii.k_uni_s_inv / full.k_uni_s_inv - 1.0).abs() < 1e-8);
        assert!((ii.lambda_1_s_inv / full.lambda_1_s_inv - 1.0).abs() < 1e-8);
    }

    #[test]
    fn a_network_without_losses_is_an_error() {
        // No product channel and no sink: nothing leaves the network, k_uni = 0.
        let mut well = test_well("A", 300, 0, 200);
        well.channels.clear();
        let network = ChemicalActivationNetwork { grain_width_cm1: 10.0, wells: vec![well] };
        let op = assemble_operator(&network, &conditions(), &final_options(EXPONENTIAL)).unwrap();
        assert!(thermal_rate_coefficients(&network, &op, EigenSolver::FullDecomposition, DEFAULT_SUM_RULE_TOLERANCE).is_err());
    }

    #[test]
    fn at_high_pressure_lambda_1_is_the_boltzmann_average_of_k() {
        // Single well, collisions fast compared with reaction: the thermal distribution is Boltzmann and
        // lambda_1 -> k_inf = sum k f / sum f.
        let network = ChemicalActivationNetwork { grain_width_cm1: 10.0, wells: vec![test_well("A", 400, 0, 200)] };
        let high = Conditions { temperature_kelvin: 600.0, pressure_torr: 1.0e9 };
        let op = assemble_operator(&network, &high, &final_options(EXPONENTIAL)).unwrap();
        let th = thermal_rate_coefficients(&network, &op, EigenSolver::InverseIteration, DEFAULT_SUM_RULE_TOLERANCE).unwrap();
        let kt = KB_CM * 600.0;
        let well = &network.wells[0];
        let f: Vec<f64> = (0..400).map(|i| well.density_of_states[i] * (-(i as f64) * 10.0 / kt).exp()).collect();
        let k_inf = f.iter().zip(&well.channels[0].rate_constant_s_inv).map(|(f, k)| f * k).sum::<f64>() / f.iter().sum::<f64>();
        assert!((th.lambda_1_s_inv / k_inf - 1.0).abs() < 1e-3, "{} vs {k_inf}", th.lambda_1_s_inv);
        assert!((th.channels[0].high_pressure_rate_s_inv / k_inf - 1.0).abs() < 1e-12);
        // Falloff: at low pressure lambda_1 lies below the high-pressure limit.
        let low = Conditions { temperature_kelvin: 600.0, pressure_torr: 1.0 };
        let op = assemble_operator(&network, &low, &final_options(EXPONENTIAL)).unwrap();
        let th_low = thermal_rate_coefficients(&network, &op, EigenSolver::InverseIteration, DEFAULT_SUM_RULE_TOLERANCE).unwrap();
        assert!(th_low.lambda_1_s_inv < 0.9 * k_inf);
    }

    #[test]
    fn time_dependent_population_approaches_the_steady_state() {
        // PO14 eq. 3: N(t -> inf) = R J^-1 F and N(t) ~ R F t for t << 1/lambda_max.
        let network = two_well_network();
        let op = assemble_operator(&network, &conditions(), &final_options(EXPONENTIAL)).unwrap();
        let source: Vec<Vec<f64>> = vec![
            (0..400).map(|i| if i >= 320 { (-((i - 320) as f64) / 15.0).exp() } else { 0.0 }).collect(),
            vec![0.0; 460],
        ];
        let f = project_source(&op, &source).unwrap();
        let steady = solve_steady_state(&op, &f.on_states, &LinearSolver::BandedCholesky).unwrap();
        let system = EigenSystem::new(&op, EigenSolver::FullDecomposition).unwrap();
        let rate = 2.0;
        let late = system.time_dependent_population(&f.on_states, rate, 1.0e3 / system.eigenvalues[0]);
        let scale = steady.population.iter().cloned().fold(0.0, f64::max);
        for (a, b) in late.iter().zip(&steady.population) {
            assert!((a - rate * b).abs() < 1e-6 * rate * scale, "{a} vs {}", rate * b);
        }
        let t = 1.0e-6 / system.eigenvalues.last().unwrap();
        let early = system.time_dependent_population(&f.on_states, rate, t);
        let fmax = f.on_states.iter().cloned().fold(0.0, f64::max);
        for (a, b) in early.iter().zip(&f.on_states) {
            assert!((a - rate * b * t).abs() < 1e-4 * rate * fmax * t);
        }
    }

    #[test]
    fn lambda_f_belongs_to_the_spectrum_and_defines_the_window() {
        let network = two_well_network();
        let op = assemble_operator(&network, &conditions(), &final_options(EXPONENTIAL)).unwrap();
        let system = EigenSystem::new(&op, EigenSolver::FullDecomposition).unwrap();
        let source: Vec<Vec<f64>> = vec![(0..400).map(|i| if i == 350 { 1.0 } else { 0.0 }).collect(), vec![0.0; 460]];
        let f = project_source(&op, &source).unwrap();
        let lambda_f = system.lambda_f(&f.on_states);
        assert!(system.eigenvalues.iter().any(|&l| l == lambda_f));
        assert!(lambda_f > system.eigenvalues[0]);
        let window = steady_state_window(1.0e9, 1.0e4, 1.0e-2, lambda_f);
        assert!(window.sink_within_window); // 1e7 > 1e4 > 0.1
        assert!((window.intermediate_steady_state_from_s - 10.0 / lambda_f).abs() < 1e-12 / lambda_f);
        assert!((window.intermediate_steady_state_until_s - 10.0).abs() < 1e-12);
        assert!(!steady_state_window(1.0e9, 1.0e8, 1.0e-2, lambda_f).sink_within_window);
        assert!(!steady_state_window(1.0e9, 1.0e-2, 1.0e-2, lambda_f).sink_within_window);
    }
}


