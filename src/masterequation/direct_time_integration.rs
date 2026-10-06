//! Direct time integration of the energy-grained master equation: the third solution method of MarXus.
//!
//! The populations N of the grains of all wells and the yield Y_x accumulated in every exit of the network
//! are integrated in time,
//!   dN/dt = R F - J N,   dY_x/dt = sum_E k_x(E) N(E),
//! with J the operator of the steady-state method (`chemical_activation_operator.rs`; PO14 eq. 2) and the
//! exits x: the product channels (k_x(E)), the bimolecular sinks (k_c[D]) and, with the absorbing-barrier
//! operator of the intermediate steady state, the stabilization flux into the absorbed grains of each well.
//! The column sums of J are the losses out of the network, so sum N + sum Y changes only by the formation
//! R: the total is conserved for a pulse. No steady state is assumed and no separation of chemical and
//! relaxation modes is needed: the early relaxation and the later chemistry follow from one integration.
//!
//! Integrator: the adaptive, L-stable Rosenbrock methods of `numeric::integrators` (adapted from KPP). The
//! Jacobian is constant, -J for N and K^T for Y; the matrix G = s I - Jac of a stage (s = 1/(h gamma)) is
//! block lower triangular, so a stage is solved in two parts:
//!   (s I + J) x_N = b_N, through the banded Cholesky factor of s I + S, S = D^-1 J D symmetric
//!   (detailed balance; `chemical_activation_steady_state::symmetrize`), x_N = D (s I + S)^-1 D^-1 b_N;
//!   x_Y = (b_Y + K^T x_N) / s.
//!
//! Identities used as tests: the yields of a pulse at t -> infinity equal k_x^T J^-1 F, the yields of the
//! steady state with the same F (final steady state; intermediate steady state with the absorbing-barrier
//! operator); continuous formation approaches N^s = J^-1 R F; N(t) equals the eigenvector expansion of PO14
//! eqs. 3-4; a thermalized population decays at late times with the lowest eigenvalue of J (GO10 eq. 12).

use super::chemical_activation_network::{
    ChannelDestination, ChemicalActivationNetwork, ChemicalActivationOptions, Conditions,
};
use super::chemical_activation_operator::{assemble_operator, ChemicalActivationOperator};
use super::chemical_activation_steady_state::{project_source, symmetrize, SYMMETRY_TOLERANCE};
use crate::numeric::banded_solvers::{BandedCholeskyFactor, SymmetricBandMatrix};
use crate::numeric::integrators::rosenbrock::{
    integrate, IntegrationStatistics, RosenbrockOptions, StiffSystem,
};
use crate::numeric::integrators::rosenbrock_methods::RosenbrockMethod;

/// How the population starts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum InitialState {
    /// A pulse: N(0) = F (the normalized source), no further formation.
    Pulse,
    /// Continuous formation R F with R = 1 (per second) from N(0) = 0.
    ContinuousFormation,
}

/// Settings of the integration.
#[derive(Debug, Clone, PartialEq)]
pub struct TimeIntegrationSettings {
    pub method: RosenbrockMethod,
    pub relative_tolerance: f64,
    pub absolute_tolerance: f64,
    /// Output times in s, increasing, > 0.
    pub times_s: Vec<f64>,
}

impl Default for TimeIntegrationSettings {
    fn default() -> Self {
        Self {
            method: RosenbrockMethod::Rodas4,
            relative_tolerance: 1e-6,
            absolute_tolerance: 1e-14,
            times_s: Vec::new(),
        }
    }
}

/// The state at one output time.
#[derive(Debug, Clone, PartialEq)]
pub struct TimePoint {
    pub time_s: f64,
    /// Total population of every well.
    pub well_populations: Vec<f64>,
    /// Yield accumulated in every exit (`TimeEvolution::exits`).
    pub exit_yields: Vec<f64>,
}

/// The time evolution at one condition.
#[derive(Debug, Clone)]
pub struct TimeEvolution {
    pub conditions: Conditions,
    /// Names of the exits: `W->X` (product channel), `escape(W)` (bimolecular sink), `stab(W)` (absorbing
    /// barrier).
    pub exits: Vec<String>,
    pub points: Vec<TimePoint>,
    /// Populations of all retained states at the last output time.
    pub final_state_population: Vec<f64>,
    pub statistics: IntegrationStatistics,
    /// Cholesky factorizations of s I + S actually computed (the others were reused from the cache).
    pub factorizations_computed: usize,
}

/// Logarithmically spaced times from `t_min` to `t_max`, `per_decade` points per decade.
pub fn log_spaced_times(t_min: f64, t_max: f64, per_decade: usize) -> Vec<f64> {
    if !(t_min > 0.0 && t_max >= t_min) || per_decade == 0 {
        return Vec::new();
    }
    let intervals = ((t_max / t_min).log10() * per_decade as f64 + 1e-9).floor() as usize;
    (0..=intervals)
        .map(|k| t_min * 10f64.powf(k as f64 / per_decade as f64))
        .collect()
}

/// The master equation with the accumulated exit yields, as a stiff system for the Rosenbrock integrator.
struct MasterEquationSystem<'a> {
    op: &'a ChemicalActivationOperator,
    /// S = D^-1 J D as a band matrix, ln D and D.
    band: SymmetricBandMatrix,
    half_log_d: Vec<f64>,
    d: Vec<f64>,
    /// Rates of every exit: (state, k) with k in s-1.
    exit_rates: Vec<Vec<(usize, f64)>>,
    /// Formation R F on the states, and directly into the exits (formation below an absorbing barrier).
    formation: Vec<f64>,
    exit_formation: Vec<f64>,
    /// Cholesky factors of s I + S for the last shifts s (at most `FACTOR_CACHE`), and the one in use.
    factors: Vec<(f64, BandedCholeskyFactor)>,
    current: usize,
    /// Factorizations computed.
    computed: usize,
}

/// Number of factorizations kept for reuse (the steps are powers of two, so step sizes repeat).
const FACTOR_CACHE: usize = 4;

impl StiffSystem for MasterEquationSystem<'_> {
    fn dimension(&self) -> usize {
        self.op.dimension() + self.exit_rates.len()
    }

    fn rhs(&self, _t: f64, y: &[f64], dydt: &mut [f64]) {
        let n = self.op.dimension();
        let jn = self.op.apply(&y[..n]);
        for s in 0..n {
            dydt[s] = self.formation[s] - jn[s];
        }
        for (x, rates) in self.exit_rates.iter().enumerate() {
            dydt[n + x] =
                self.exit_formation[x] + rates.iter().map(|&(s, k)| k * y[s]).sum::<f64>();
        }
    }

    fn prepare(&mut self, _t: f64, _y: &[f64], shift: f64) -> Result<(), String> {
        if let Some(k) = self.factors.iter().position(|(cached, _)| *cached == shift) {
            self.current = k;
            return Ok(());
        }
        let factor = self.band.shifted(shift).cholesky()?;
        self.computed += 1;
        if self.factors.len() == FACTOR_CACHE {
            // Replace the oldest entry.
            self.factors.remove(0);
        }
        self.factors.push((shift, factor));
        self.current = self.factors.len() - 1;
        Ok(())
    }

    fn solve(&self, b: &mut [f64]) -> Result<(), String> {
        let (shift, factor) = self
            .factors
            .get(self.current)
            .ok_or("Time integration: solve before prepare.")?;
        let n = self.op.dimension();
        // x_N = D (s I + S)^-1 D^-1 b_N, with D^-1 b_N from logarithms (D spans many orders of magnitude).
        let scaled: Vec<f64> = b[..n]
            .iter()
            .zip(&self.half_log_d)
            .map(|(&x, &h)| {
                if x == 0.0 {
                    0.0
                } else {
                    x.signum() * (x.abs().ln() - h).exp()
                }
            })
            .collect();
        let y = factor.solve(&scaled)?;
        for s in 0..n {
            b[s] = y[s] * self.d[s];
        }
        // x_Y = (b_Y + K^T x_N) / s.
        for (x, rates) in self.exit_rates.iter().enumerate() {
            b[n + x] = (b[n + x] + rates.iter().map(|&(s, k)| k * b[s]).sum::<f64>()) / shift;
        }
        Ok(())
    }
}

/// Integrates the master equation of `network` at `conditions` with the operator of `options` (final or
/// intermediate steady state) from `initial`, with the source `source` (grain distribution of every well,
/// normalized here) as the pulse or as the continuous formation.
pub fn integrate_master_equation(
    network: &ChemicalActivationNetwork,
    conditions: &Conditions,
    options: &ChemicalActivationOptions,
    source: &[Vec<f64>],
    initial: InitialState,
    settings: &TimeIntegrationSettings,
) -> Result<TimeEvolution, String> {
    if settings.times_s.is_empty()
        || settings.times_s.iter().any(|&t| !(t > 0.0))
        || settings.times_s.windows(2).any(|w| w[1] <= w[0])
    {
        return Err("Time integration: the output times must be positive and increasing.".into());
    }
    let op = assemble_operator(network, conditions, options)?;
    let n = op.dimension();
    let symmetrized = symmetrize(&op);
    if symmetrized.max_relative_asymmetry > SYMMETRY_TOLERANCE {
        return Err(format!(
            "Time integration: the operator is not symmetrizable (relative asymmetry {:e}); the isomerization \
             rates violate detailed balance.",
            symmetrized.max_relative_asymmetry
        ));
    }
    let projected = project_source(&op, source)?;

    // Exits: product channels, bimolecular sinks, stabilization into the absorbed grains of each well.
    let mut exits = Vec::new();
    let mut exit_rates: Vec<Vec<(usize, f64)>> = Vec::new();
    let mut exit_formation = Vec::new();
    for (w, well) in network.wells.iter().enumerate() {
        for channel in &well.channels {
            if let ChannelDestination::Products { name } = &channel.destination {
                exits.push(format!("{}->{name}", well.name));
                exit_rates.push(
                    op.states
                        .iter()
                        .enumerate()
                        .filter(|&(_, &(sw, i))| sw == w && channel.rate_constant_s_inv[i] > 0.0)
                        .map(|(s, &(_, i))| (s, channel.rate_constant_s_inv[i]))
                        .collect(),
                );
                exit_formation.push(0.0);
            }
        }
    }
    for (w, well) in network.wells.iter().enumerate() {
        if well.bimolecular_sink_s_inv > 0.0 {
            exits.push(format!("escape({})", well.name));
            exit_rates.push(
                op.states
                    .iter()
                    .enumerate()
                    .filter(|&(_, &(sw, _))| sw == w)
                    .map(|(s, _)| (s, well.bimolecular_sink_s_inv))
                    .collect(),
            );
            exit_formation.push(0.0);
        }
    }
    let absorbing = op.stabilization.iter().any(|targets| !targets.is_empty())
        || projected.absorbed_per_well.iter().any(|&a| a > 0.0);
    if absorbing {
        for (w, well) in network.wells.iter().enumerate() {
            exits.push(format!("stab({})", well.name));
            exit_rates.push(
                op.stabilization
                    .iter()
                    .enumerate()
                    .flat_map(|(s, targets)| {
                        targets
                            .iter()
                            .filter(|&&(target, _)| target == w)
                            .map(move |&(_, k)| (s, k))
                    })
                    .collect(),
            );
            exit_formation.push(projected.absorbed_per_well[w]);
        }
    }

    // Initial state and formation.
    let m = exits.len();
    let mut y = vec![0.0; n + m];
    let (formation, exit_formation) = match initial {
        InitialState::Pulse => {
            y[..n].copy_from_slice(&projected.on_states);
            for (x, &f) in exit_formation.iter().enumerate() {
                y[n + x] = f;
            }
            (vec![0.0; n], vec![0.0; m])
        }
        InitialState::ContinuousFormation => (projected.on_states.clone(), exit_formation),
    };

    let mut system = MasterEquationSystem {
        op: &op,
        band: symmetrized.band_matrix(),
        half_log_d: symmetrized.half_log_d.clone(),
        d: symmetrized.d.clone(),
        exit_rates,
        formation,
        exit_formation,
        factors: Vec::new(),
        current: 0,
        computed: 0,
    };
    let mut rosenbrock = RosenbrockOptions::new(
        settings.method,
        settings.relative_tolerance,
        settings.absolute_tolerance,
    );
    rosenbrock.autonomous = true;
    // Constant Jacobian: steps of powers of two let the factorizations of s I + S be reused.
    rosenbrock.power_of_two_steps = true;
    let mut statistics = IntegrationStatistics::default();
    let mut points = Vec::new();
    let mut t = 0.0;
    let mut h_next = 1e-3 * settings.times_s[0];
    for &t_out in &settings.times_s {
        rosenbrock.h_start = h_next;
        let stats = integrate(&mut system, &mut y, t, t_out, &rosenbrock).map_err(|e| {
            format!(
                "Time integration at T = {} K, p = {} Torr: {e}",
                conditions.temperature_kelvin, conditions.pressure_torr
            )
        })?;
        statistics.function_evaluations += stats.function_evaluations;
        statistics.factorizations += stats.factorizations;
        statistics.solves += stats.solves;
        statistics.steps += stats.steps;
        statistics.accepted_steps += stats.accepted_steps;
        statistics.rejected_steps += stats.rejected_steps;
        statistics.last_step = stats.last_step;
        statistics.next_step = stats.next_step;
        if stats.next_step > 0.0 {
            h_next = stats.next_step;
        }
        t = t_out;
        let mut well_populations = vec![0.0; network.wells.len()];
        for (s, &(w, _)) in op.states.iter().enumerate() {
            well_populations[w] += y[s];
        }
        points.push(TimePoint {
            time_s: t_out,
            well_populations,
            exit_yields: y[n..].to_vec(),
        });
    }
    Ok(TimeEvolution {
        conditions: conditions.clone(),
        exits,
        points,
        final_state_population: y[..n].to_vec(),
        statistics,
        factorizations_computed: system.computed,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::KB_CM;
    use crate::masterequation::chemical_activation_driver::{
        run_chemical_activation, ChemicalActivationRun, SourceSpecification,
    };
    use crate::masterequation::chemical_activation_eigen::{
        thermal_rate_coefficients, EigenSolver, EigenSystem, DEFAULT_SUM_RULE_TOLERANCE,
    };
    use crate::masterequation::chemical_activation_network::{
        AbsorbingBarrier, Channel, ChannelDestination, CollisionModel, SteadyState,
    };
    use crate::masterequation::chemical_activation_operator::assemble_operator;
    use crate::masterequation::chemical_activation_operator::tests::two_well_network;
    use crate::masterequation::chemical_activation_sources::{
        thermal_distribution, thermal_entrance_source,
    };
    use crate::masterequation::chemical_activation_steady_state::{project_source, LinearSolver};

    const MODEL: CollisionModel = CollisionModel::ExponentialDown {
        cutoff_in_mean_down: 10.0,
    };

    /// Two wells A, B (B with a sink), entrance A <- R opening at grain 320.
    fn network() -> ChemicalActivationNetwork {
        let mut network = two_well_network();
        network.wells[0].channels.push(Channel {
            name: "A->reactants".into(),
            destination: ChannelDestination::Products { name: "R".into() },
            threshold_grain: None,
            rate_constant_s_inv: (0..400)
                .map(|i| {
                    if i >= 320 {
                        2.0e7 * ((i - 320) as f64 + 1.0)
                    } else {
                        0.0
                    }
                })
                .collect(),
        });
        network
    }

    fn conditions(t: f64, p: f64) -> Conditions {
        let mut c = crate::masterequation::chemical_activation_operator::tests::conditions();
        c.temperature_kelvin = t;
        c.pressure_torr = p;
        c
    }

    fn options(steady_state: SteadyState) -> ChemicalActivationOptions {
        ChemicalActivationOptions {
            collision_model: MODEL,
            steady_state,
        }
    }

    fn source(network: &ChemicalActivationNetwork, t: f64) -> Vec<Vec<f64>> {
        thermal_entrance_source(network, &[(0, 2)], KB_CM * t).unwrap()
    }

    fn settings(times_s: Vec<f64>) -> TimeIntegrationSettings {
        TimeIntegrationSettings {
            method: RosenbrockMethod::Rodas4,
            relative_tolerance: 1e-8,
            absolute_tolerance: 1e-16,
            times_s,
        }
    }

    #[test]
    fn log_spaced_times_cover_the_decades() {
        let t = log_spaced_times(1e-9, 1e-6, 2);
        assert_eq!(t.len(), 7);
        assert!((t[0] - 1e-9).abs() < 1e-24 && (t[6] / 1e-6 - 1.0).abs() < 1e-12);
        assert!((t[1] / (1e-9 * 10f64.sqrt()) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn a_pulse_conserves_the_total_and_ends_in_the_final_steady_state_yields() {
        // Y_x(infinity) = k_x^T J^-1 F: the yields of the final steady state with the same F.
        let network = network();
        let (t, p) = (300.0, 760.0);
        let f = source(&network, t);
        let evolution = integrate_master_equation(
            &network,
            &conditions(t, p),
            &options(SteadyState::Final),
            &f,
            InitialState::Pulse,
            &settings(log_spaced_times(1e-12, 10.0, 4)),
        )
        .unwrap();
        for point in &evolution.points {
            let total: f64 =
                point.well_populations.iter().sum::<f64>() + point.exit_yields.iter().sum::<f64>();
            assert!(
                (total - 1.0).abs() < 1e-9,
                "t = {:e}: total {total}",
                point.time_s
            );
        }
        let last = evolution.points.last().unwrap();
        assert!(
            last.well_populations.iter().sum::<f64>() < 1e-9,
            "{:?}",
            last.well_populations
        );
        // Steps of powers of two: the factorizations of s I + S are reused.
        assert!(
            evolution.factorizations_computed >= 1
                && evolution.factorizations_computed * 2 < evolution.statistics.steps,
            "{} factorizations for {} steps",
            evolution.factorizations_computed,
            evolution.statistics.steps
        );
        let run = ChemicalActivationRun {
            temperatures_kelvin: vec![t],
            pressures_torr: vec![p],
            options: options(SteadyState::Final),
            solver: LinearSolver::BandedCholesky,
            source: SourceSpecification::Fixed(f.clone()),
            tolerance: 1e-8,
        };
        let steady = run_chemical_activation(&network, &run)
            .unwrap()
            .remove(0)
            .result;
        // Exits: the product channels of the wells in order, then the sinks.
        let mut expected: Vec<f64> = steady
            .channels
            .iter()
            .filter(|c| matches!(c.destination, ChannelDestination::Products { .. }))
            .map(|c| c.flux)
            .collect();
        expected.extend(
            network
                .wells
                .iter()
                .zip(&steady.wells)
                .filter(|(w, _)| w.bimolecular_sink_s_inv > 0.0)
                .map(|(_, r)| r.bimolecular_sink_yield),
        );
        assert_eq!(
            evolution.exits.len(),
            expected.len(),
            "{:?}",
            evolution.exits
        );
        for ((name, y), e) in evolution.exits.iter().zip(&last.exit_yields).zip(&expected) {
            assert!((y / e - 1.0).abs() < 1e-6, "{name}: {y} vs {e}");
        }
    }

    #[test]
    fn continuous_formation_approaches_the_steady_state_population() {
        let network = network();
        let (t, p) = (300.0, 760.0);
        let f = source(&network, t);
        let evolution = integrate_master_equation(
            &network,
            &conditions(t, p),
            &options(SteadyState::Final),
            &f,
            InitialState::ContinuousFormation,
            &settings(log_spaced_times(1e-12, 10.0, 2)),
        )
        .unwrap();
        let op =
            assemble_operator(&network, &conditions(t, p), &options(SteadyState::Final)).unwrap();
        let projected = project_source(&op, &f).unwrap();
        let steady = crate::masterequation::chemical_activation_steady_state::solve_steady_state(
            &op,
            &projected.on_states,
            &LinearSolver::BandedCholesky,
        )
        .unwrap()
        .population;
        let norm: f64 = steady.iter().sum();
        let deviation: f64 = evolution
            .final_state_population
            .iter()
            .zip(&steady)
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>()
            / norm;
        assert!(deviation < 1e-6, "{deviation}");
    }

    #[test]
    fn populations_follow_the_eigenvector_expansion() {
        // PO14 eqs. 3-4: N(t) = sum_i e_i (1 - exp(-lambda_i t))/lambda_i E_i for continuous formation.
        let network = network();
        let (t, p) = (300.0, 760.0);
        let f = source(&network, t);
        let times = vec![1e-11, 1e-9, 1e-7, 1e-5];
        let evolution = integrate_master_equation(
            &network,
            &conditions(t, p),
            &options(SteadyState::Final),
            &f,
            InitialState::ContinuousFormation,
            &settings(times.clone()),
        )
        .unwrap();
        let op =
            assemble_operator(&network, &conditions(t, p), &options(SteadyState::Final)).unwrap();
        let projected = project_source(&op, &f).unwrap();
        let eigen = EigenSystem::new(&op, EigenSolver::FullDecomposition).unwrap();
        for (point, &time) in evolution.points.iter().zip(&times) {
            let expansion = eigen.time_dependent_population(&projected.on_states, 1.0, time);
            let mut wells = vec![0.0; network.wells.len()];
            for (s, &(w, _)) in op.states.iter().enumerate() {
                wells[w] += expansion[s];
            }
            let scale: f64 = wells.iter().sum();
            for (a, b) in point.well_populations.iter().zip(&wells) {
                assert!((a - b).abs() < 1e-6 * scale, "t = {time:e}: {a} vs {b}");
            }
        }
    }

    #[test]
    fn a_thermalized_well_decays_with_the_lowest_eigenvalue() {
        // Late-time decay rate -d ln(sum N)/dt -> lambda_1 (GO10 eq. 12).
        let network = network();
        let (t, p) = (300.0, 760.0);
        let mut f: Vec<Vec<f64>> = network
            .wells
            .iter()
            .map(|w| vec![0.0; w.grain_count()])
            .collect();
        f[0] = thermal_distribution(
            &network.wells[0].density_of_states,
            network.grain_width_cm1,
            KB_CM * t,
        )
        .unwrap();
        let op =
            assemble_operator(&network, &conditions(t, p), &options(SteadyState::Final)).unwrap();
        let lambda_1 = thermal_rate_coefficients(
            &network,
            &op,
            EigenSolver::FullDecomposition,
            DEFAULT_SUM_RULE_TOLERANCE,
        )
        .unwrap()
        .k_uni_s_inv;
        let (t1, t2) = (5.0 / lambda_1, 6.0 / lambda_1);
        let evolution = integrate_master_equation(
            &network,
            &conditions(t, p),
            &options(SteadyState::Final),
            &f,
            InitialState::Pulse,
            &settings(vec![t1, t2]),
        )
        .unwrap();
        let n1: f64 = evolution.points[0].well_populations.iter().sum();
        let n2: f64 = evolution.points[1].well_populations.iter().sum();
        let rate = (n1 / n2).ln() / (t2 - t1);
        assert!((rate / lambda_1 - 1.0).abs() < 1e-4, "{rate} vs {lambda_1}");
    }

    #[test]
    fn with_the_absorbing_barrier_a_pulse_ends_in_the_intermediate_steady_state_yields() {
        let network = network();
        let (t, p) = (300.0, 760.0);
        let f = source(&network, t);
        let intermediate = SteadyState::Intermediate {
            barrier: AbsorbingBarrier::default(),
        };
        let evolution = integrate_master_equation(
            &network,
            &conditions(t, p),
            &options(intermediate.clone()),
            &f,
            InitialState::Pulse,
            &settings(log_spaced_times(1e-12, 1.0, 4)),
        )
        .unwrap();
        let run = ChemicalActivationRun {
            temperatures_kelvin: vec![t],
            pressures_torr: vec![p],
            options: options(intermediate),
            solver: LinearSolver::BandedCholesky,
            source: SourceSpecification::Fixed(f.clone()),
            tolerance: 1e-8,
        };
        let steady = run_chemical_activation(&network, &run)
            .unwrap()
            .remove(0)
            .result;
        let last = evolution.points.last().unwrap();
        for (w, well) in network.wells.iter().enumerate() {
            let k = evolution
                .exits
                .iter()
                .position(|x| *x == format!("stab({})", well.name))
                .expect("a stabilization exit");
            let e = steady.wells[w].stabilization_yield;
            assert!(
                (last.exit_yields[k] - e).abs() < 1e-6 * e.max(1e-12),
                "stab({}): {} vs {e}",
                well.name,
                last.exit_yields[k]
            );
        }
    }
}
