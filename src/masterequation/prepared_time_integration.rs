//! Time integration of a prepared experiment: an initial population, any number of source channels with their
//! own distributions and time profiles, and a piecewise-constant bath history (design note N,
//! papers/Reactant_flux_initiation/MarXus_Nonthermal_Sources.tex, eqs. main-normal, product, balance; Secs. 8, 9.1,
//! 10.3; reports/nonthermal_sources_design.md, Section 4.4):
//!   dn/dt = -J[T_b(t), p(t)] n + sum_a R_a(t) F_a,   n(0) = N_0 F_0,   dY_x/dt = k_x^T n.
//! Events, where the integrator stops and restarts its step control: the output times, the impulses (exact jumps
//! n(t_p+) = n(t_p-) + N_p F_a), the breakpoints of the rate profiles and the starts of the bath segments. Between
//! events the Rosenbrock integrator runs autonomous where every amplitude is constant and nonautonomous otherwise
//! (df/dt = dS/dt at fixed n, J constant inside a segment). At a bath boundary the operator is assembled for the new
//! segment on the same grain grid; populations are carried over through the grains (a reservoir state is spread
//! with the Boltzmann weights of the old bath and the grains are summed into the new partition), so the population
//! vector is continuous across the jump (N eq. shockhistory). Every segment has its own system and its own
//! factorization cache, so no factor of an earlier operator is reused (N Sec. 10.3).
//!
//! Source mass that falls below an absorbing barrier is counted at once in `stab(W)` (N Sec. 9.4); source mass in
//! the grains of a low-energy reservoir keeps its total but takes the reservoir's Boltzmann shape (N, keypoint of
//! Sec. 8), which `ProjectionNote` reports for every projection.

use super::chemical_activation_network::{
    ChemicalActivationNetwork, ChemicalActivationOptions, CollisionModel, Conditions, SteadyState,
};
use super::chemical_activation_operator::{assemble_operator, ChemicalActivationOperator};
use super::direct_time_integration::{exit_definitions, MasterEquationSystem, SourceTerm, TimeIntegrationSettings};
use super::prepared_distributions::GrainDistribution;
use super::source_profiles::TimeProfile;
use crate::numeric::integrators::rosenbrock::{integrate, IntegrationStatistics, RosenbrockOptions};

/// The population present at t = 0: an amount with its distribution.
#[derive(Debug, Clone, PartialEq)]
pub struct InitialPopulation {
    pub amount: f64,
    pub distribution: GrainDistribution,
}

/// A source channel: a normalized distribution and the time profile of its amplitude.
#[derive(Debug, Clone, PartialEq)]
pub struct SourceChannel {
    pub name: String,
    pub distribution: GrainDistribution,
    pub profile: TimeProfile,
}

/// A bath condition from `start_s` until the next segment.
#[derive(Debug, Clone)]
pub struct BathSegment {
    pub start_s: f64,
    pub conditions: Conditions,
}

/// A prepared experiment (N Sec. 1): initial population, source channels and bath history.
#[derive(Debug, Clone)]
pub struct Preparation {
    pub initial: Option<InitialPopulation>,
    pub channels: Vec<SourceChannel>,
    /// Starts increasing, the first at 0.
    pub bath: Vec<BathSegment>,
}

/// What a projection of a distribution onto the states of a segment kept: per well, the fraction of the
/// distribution in the grains of a low-energy reservoir (its shape there is replaced by the reservoir's Boltzmann
/// shape) and the fraction below an absorbing barrier (stabilized at once), and the mean energy above the well
/// bottom before and after the projection.
#[derive(Debug, Clone, PartialEq)]
pub struct ProjectionNote {
    /// "initial" or the channel name.
    pub what: String,
    pub segment: usize,
    pub reservoir_fraction: Vec<f64>,
    pub absorbed_fraction: Vec<f64>,
    pub mean_energy_before_cm1: Vec<f64>,
    pub mean_energy_after_cm1: Vec<f64>,
}

/// The state at one output time (N Sec. observables).
#[derive(Debug, Clone, PartialEq)]
pub struct TransientPoint {
    pub time_s: f64,
    /// C_w = sum_i n_wi.
    pub well_populations: Vec<f64>,
    /// Y_x, cumulative (`TransientResult::exits`).
    pub exit_yields: Vec<f64>,
    /// q_x = k_x^T n, instantaneous.
    pub exit_fluxes: Vec<f64>,
    /// N_a,in(0, t) of every channel (N eq. injected).
    pub injected: Vec<f64>,
    /// Mean energy above the bottom of every well of its surviving population (NaN for an empty well).
    pub mean_energy_above_bottom_cm1: Vec<f64>,
    /// Fraction of the population of every well at or above its lowest reaction threshold (N eq. shockobs).
    pub tail_fraction: Vec<f64>,
    /// k_inst = sum_x q_x / sum_w C_w (equals -d ln N/dt while no source is active; N eq. instant).
    pub loss_hazard_s_inv: f64,
    /// (sum n + sum Y - N_0 - sum_a N_a,in) / (N_0 + sum_a N_a,in) (N eq. balance).
    pub balance_deviation: f64,
}

/// Result of `integrate_preparation`.
#[derive(Debug, Clone)]
pub struct TransientResult {
    pub exits: Vec<String>,
    pub channels: Vec<String>,
    pub points: Vec<TransientPoint>,
    /// The bath segments that were integrated (start, conditions).
    pub segments: Vec<(f64, Conditions)>,
    pub projections: Vec<ProjectionNote>,
    pub statistics: IntegrationStatistics,
    pub factorizations_computed: usize,
}

/// Projection of `amount` x `distribution` onto the states of `op`: (on states, absorbed per well, note).
fn project(
    network: &ChemicalActivationNetwork,
    op: &ChemicalActivationOperator,
    distribution: &GrainDistribution,
    amount: f64,
    what: &str,
    segment: usize,
) -> Result<(Vec<f64>, Vec<f64>, ProjectionNote), String> {
    if distribution.mass.len() != network.wells.len()
        || distribution.mass.iter().zip(&network.wells).any(|(m, w)| m.len() != w.grain_count())
    {
        return Err(format!("Preparation '{what}': the distribution is not on the grid of the network."));
    }
    let mut on_states = vec![0.0; op.dimension()];
    let mut absorbed = vec![0.0; network.wells.len()];
    let mut in_reservoir = vec![0.0; network.wells.len()];
    for (w, masses) in distribution.mass.iter().enumerate() {
        let reservoir_state = op.reservoirs[w].as_ref().map(|r| r.state);
        for (i, &x) in masses.iter().enumerate() {
            match op.index_of[w][i] {
                Some(s) => {
                    on_states[s] += amount * x;
                    if Some(s) == reservoir_state {
                        in_reservoir[w] += x;
                    }
                }
                None => absorbed[w] += amount * x,
            }
        }
    }
    let before = distribution.mean_energy_above_bottom_cm1(network);
    let after_grains = op.grain_populations(&on_states);
    let after = GrainDistribution { mass: after_grains, lost_fraction: 0.0, description: String::new() }.mean_energy_above_bottom_cm1(network);
    let note = ProjectionNote {
        what: what.to_string(),
        segment,
        reservoir_fraction: in_reservoir,
        absorbed_fraction: absorbed.iter().map(|a| if amount > 0.0 { a / amount } else { 0.0 }).collect(),
        mean_energy_before_cm1: before,
        mean_energy_after_cm1: after,
    };
    Ok((on_states, absorbed, note))
}

/// One bath segment: its operator, exits and projected source shapes.
struct Segment {
    op: ChemicalActivationOperator,
    exits: Vec<String>,
    exit_rates: Vec<Vec<(usize, f64)>>,
    /// Per channel: F on the states and the absorbed part per well (amounts per unit of the channel).
    shapes: Vec<(Vec<f64>, Vec<f64>)>,
}

/// Integrates `preparation` on `network` with the collision model and the operator variant `steady_state` (final:
/// no absorbing barrier; intermediate: with it) to the output times of `settings`.
pub fn integrate_preparation(
    network: &ChemicalActivationNetwork,
    collision_model: CollisionModel,
    steady_state: SteadyState,
    preparation: &Preparation,
    settings: &TimeIntegrationSettings,
) -> Result<TransientResult, String> {
    let times = &settings.times_s;
    if times.is_empty() || times.iter().any(|&t| !(t > 0.0) || !t.is_finite()) || times.windows(2).any(|w| w[1] <= w[0]) {
        return Err("Time integration: the output times must be positive and increasing.".into());
    }
    let bath = &preparation.bath;
    if bath.is_empty() || bath[0].start_s != 0.0 || bath.windows(2).any(|w| w[1].start_s <= w[0].start_s) {
        return Err("Preparation: the bath segments must start at 0 s and follow in increasing order.".into());
    }
    for channel in &preparation.channels {
        channel.profile.validate().map_err(|e| format!("Source '{}': {e}", channel.name))?;
    }
    if let Some(initial) = &preparation.initial {
        if !(initial.amount >= 0.0) || !initial.amount.is_finite() {
            return Err("Preparation: the initial amount must be finite and >= 0.".into());
        }
    }
    let t_end = times[times.len() - 1];

    // Events: output times, impulses, rate breakpoints and segment starts inside (0, t_end].
    let mut events: Vec<f64> = times.clone();
    for channel in &preparation.channels {
        events.extend(channel.profile.impulses().iter().map(|&(t, _)| t));
        events.extend(channel.profile.breakpoints());
    }
    events.extend(bath.iter().map(|b| b.start_s));
    events.retain(|&t| t > 0.0 && t <= t_end);
    events.sort_by(|a, b| a.total_cmp(b));
    events.dedup();

    let options_of = |_: &Conditions| ChemicalActivationOptions { collision_model, steady_state: steady_state.clone() };
    let absorbing_barrier = matches!(steady_state, SteadyState::Intermediate { .. });
    let mut projections = Vec::new();
    let build_segment = |k: usize, projections: &mut Vec<ProjectionNote>| -> Result<Segment, String> {
        let conditions = &bath[k].conditions;
        let op = assemble_operator(network, conditions, &options_of(conditions))?;
        let (exits, exit_rates) = exit_definitions(network, &op, absorbing_barrier);
        let mut shapes = Vec::new();
        for channel in &preparation.channels {
            let (on_states, absorbed, note) = project(network, &op, &channel.distribution, 1.0, &channel.name, k)?;
            projections.push(note);
            shapes.push((on_states, absorbed));
        }
        Ok(Segment { op, exits, exit_rates, shapes })
    };

    // Global exit list (by name) and the population state carried between segments, on the grains.
    let mut segment_index = 0;
    let mut segment = build_segment(0, &mut projections)?;
    let exits = segment.exits.clone();
    let stab_index = |w: usize| exits.iter().position(|e| *e == format!("stab({})", network.wells[w].name));
    let mut y_states = vec![0.0; segment.op.dimension()];
    let mut yields = vec![0.0; exits.len()];
    let mut amount_in = 0.0;
    let add_absorbed = |yields: &mut Vec<f64>, absorbed: &[f64], scale: f64| -> Result<(), String> {
        for (w, &a) in absorbed.iter().enumerate() {
            if a > 0.0 {
                let x = stab_index(w).ok_or("Preparation: source mass below an absorbing barrier without a stab exit.")?;
                yields[x] += scale * a;
            }
        }
        Ok(())
    };
    if let Some(initial) = &preparation.initial {
        let (on_states, absorbed, note) = project(network, &segment.op, &initial.distribution, initial.amount, "initial", 0)?;
        projections.push(note);
        y_states.iter_mut().zip(&on_states).for_each(|(y, x)| *y += x);
        add_absorbed(&mut yields, &absorbed, 1.0)?;
        amount_in += initial.amount;
    }
    let apply_impulses = |t: f64, segment: &Segment, y_states: &mut Vec<f64>, yields: &mut Vec<f64>| -> Result<(), String> {
        for (channel, (on_states, absorbed)) in preparation.channels.iter().zip(&segment.shapes) {
            for (tp, amount) in channel.profile.impulses() {
                if tp == t {
                    y_states.iter_mut().zip(on_states).for_each(|(y, x)| *y += amount * x);
                    add_absorbed(yields, absorbed, amount)?;
                }
            }
        }
        Ok(())
    };
    apply_impulses(0.0, &segment, &mut y_states, &mut yields)?;

    let mut statistics = IntegrationStatistics::default();
    let mut factorizations = 0;
    let mut points = Vec::new();
    let mut segments_used = vec![(0.0, bath[0].conditions.clone())];
    let mut t = 0.0;
    let mut h_next = 1e-3 * times[0].min(events.first().copied().unwrap_or(times[0]));
    for &t_event in &events {
        // Integrate (t, t_event) with the current segment.
        {
            let exit_map: Vec<usize> = segment.exits.iter().map(|e| exits.iter().position(|g| g == e).unwrap()).collect();
            let mut system = MasterEquationSystem::new(&segment.op, segment.exit_rates.clone())?;
            for (channel, (on_states, absorbed)) in preparation.channels.iter().zip(&segment.shapes) {
                let mut into_exits = vec![0.0; segment.exits.len()];
                for (w, &a) in absorbed.iter().enumerate() {
                    if a > 0.0 {
                        let name = format!("stab({})", network.wells[w].name);
                        let x = segment.exits.iter().position(|e| *e == name).ok_or("Preparation: absorbed source without a stab exit.")?;
                        into_exits[x] = a;
                    }
                }
                system.sources.push(SourceTerm { on_states: on_states.clone(), into_exits, profile: channel.profile.clone() });
            }
            let n = segment.op.dimension();
            let mut y = vec![0.0; n + segment.exits.len()];
            y[..n].copy_from_slice(&y_states);
            for (x, &g) in exit_map.iter().enumerate() {
                y[n + x] = yields[g];
            }
            let mut rosenbrock = RosenbrockOptions::new(settings.method, settings.relative_tolerance, settings.absolute_tolerance);
            rosenbrock.autonomous = preparation.channels.iter().all(|c| c.profile.is_constant_on(t, t_event));
            rosenbrock.power_of_two_steps = true;
            rosenbrock.h_start = h_next.min(t_event - t);
            let conditions = &bath[segment_index].conditions;
            let stats = integrate(&mut system, &mut y, t, t_event, &rosenbrock).map_err(|e| {
                format!("Time integration at T = {} K, p = {} Torr: {e}", conditions.temperature_kelvin, conditions.pressure_torr)
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
            factorizations += system.factorizations_computed();
            y_states.copy_from_slice(&y[..n]);
            for (x, &g) in exit_map.iter().enumerate() {
                yields[g] = y[n + x];
            }
        }
        t = t_event;
        // A new bath segment starts here: rebuild the operator and carry the population over the grains.
        if let Some(k) = bath.iter().position(|b| b.start_s == t) {
            let grains = segment.op.grain_populations(&y_states);
            segment_index = k;
            segment = build_segment(k, &mut projections)?;
            let mut moved = 0.0;
            for (w, g) in grains.iter().enumerate() {
                for (i, &x) in g.iter().enumerate() {
                    if segment.op.index_of[w][i].is_none() && x != 0.0 {
                        let s = stab_index(w).ok_or("Preparation: population below an absorbing barrier without a stab exit.")?;
                        yields[s] += x;
                        moved += x;
                    }
                }
            }
            let _ = moved;
            y_states = segment.op.state_populations(&grains);
            segments_used.push((t, bath[k].conditions.clone()));
        }
        apply_impulses(t, &segment, &mut y_states, &mut yields)?;
        if times.contains(&t) {
            points.push(observe(network, &segment, &exits, &y_states, &yields, preparation, amount_in, t));
        }
    }
    Ok(TransientResult {
        exits,
        channels: preparation.channels.iter().map(|c| c.name.clone()).collect(),
        points,
        segments: segments_used,
        projections,
        statistics,
        factorizations_computed: factorizations,
    })
}

/// The single source shape of a preparation for the steady-state methods: with open-ended feeds, the shape of the
/// steady continuous source F_eff = sum_a R_a F_a / sum_a R_a over those feeds (N eq. effectivesource); otherwise the
/// amount-weighted shape of everything injected, N_0 F_0 + sum_a N_a F_a (normalized), whose steady-state yields equal
/// the yields of the whole preparation at infinite time in a constant bath (N eqs. pulseyield, ssyield). Returns the
/// shape and how it was formed.
pub fn steady_source_shape(preparation: &Preparation) -> Result<(GrainDistribution, String), String> {
    let feeds: Vec<(f64, &GrainDistribution, &str)> = preparation
        .channels
        .iter()
        .filter_map(|c| match &c.profile {
            TimeProfile::Feed { end_s: None, rate, .. } if *rate > 0.0 => Some((*rate, &c.distribution, c.name.as_str())),
            _ => None,
        })
        .collect();
    let (parts, how): (Vec<(f64, &GrainDistribution)>, String) = if !feeds.is_empty() {
        (
            feeds.iter().map(|&(r, d, _)| (r, d)).collect(),
            format!("open-ended feeds {}, weighted by their rates", feeds.iter().map(|f| format!("'{}'", f.2)).collect::<Vec<_>>().join(", ")),
        )
    } else {
        let mut parts = Vec::new();
        let mut names = Vec::new();
        if let Some(initial) = &preparation.initial {
            parts.push((initial.amount, &initial.distribution));
            names.push("initial population".to_string());
        }
        for c in &preparation.channels {
            let amount = c.profile.total_amount().ok_or_else(|| format!("Source '{}': no total amount.", c.name))?;
            parts.push((amount, &c.distribution));
            names.push(format!("'{}'", c.name));
        }
        (parts, format!("everything injected ({}), weighted by its amount", names.join(", ")))
    };
    let total: f64 = parts.iter().map(|(w, _)| w).sum();
    if !(total > 0.0) {
        return Err("Preparation: nothing is injected, so there is no source shape for a steady state.".into());
    }
    let shape = super::prepared_distributions::mixture(&parts.iter().map(|&(w, d)| (w / total, d.clone())).collect::<Vec<_>>())?;
    Ok((shape, how))
}

/// The observables of N Sec. observables at time `t`.
#[allow(clippy::too_many_arguments)]
fn observe(
    network: &ChemicalActivationNetwork,
    segment: &Segment,
    exits: &[String],
    y_states: &[f64],
    yields: &[f64],
    preparation: &Preparation,
    initial_amount: f64,
    t: f64,
) -> TransientPoint {
    let op = &segment.op;
    let mut well_populations = vec![0.0; network.wells.len()];
    for (s, &(w, _)) in op.states.iter().enumerate() {
        well_populations[w] += y_states[s];
    }
    let mut exit_fluxes = vec![0.0; exits.len()];
    for (x, rates) in segment.exit_rates.iter().enumerate() {
        let g = exits.iter().position(|e| *e == segment.exits[x]).unwrap();
        exit_fluxes[g] = rates.iter().map(|&(s, k)| k * y_states[s]).sum();
    }
    let grains = op.grain_populations(y_states);
    let de = network.grain_width_cm1;
    let mut mean_energy = Vec::new();
    let mut tail = Vec::new();
    for (w, g) in grains.iter().enumerate() {
        let total: f64 = g.iter().sum();
        if total > 0.0 {
            mean_energy.push(g.iter().enumerate().map(|(i, x)| i as f64 * de * x).sum::<f64>() / total);
            let threshold = network.wells[w].lowest_threshold_grain().unwrap_or(g.len());
            tail.push(g.iter().skip(threshold).sum::<f64>() / total);
        } else {
            mean_energy.push(f64::NAN);
            tail.push(f64::NAN);
        }
    }
    let injected: Vec<f64> = preparation.channels.iter().map(|c| c.profile.injected_until(t)).collect();
    let total_in = initial_amount + injected.iter().sum::<f64>();
    let present: f64 = y_states.iter().sum::<f64>() + yields.iter().sum::<f64>();
    let population: f64 = well_populations.iter().sum();
    TransientPoint {
        time_s: t,
        well_populations,
        exit_yields: yields.to_vec(),
        loss_hazard_s_inv: if population > 0.0 { exit_fluxes.iter().sum::<f64>() / population } else { f64::NAN },
        exit_fluxes,
        injected,
        mean_energy_above_bottom_cm1: mean_energy,
        tail_fraction: tail,
        balance_deviation: if total_in > 0.0 { (present - total_in) / total_in } else { present },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::masterequation::chemical_activation_network::{
        AbsorbingBarrier, ChemicalActivationNetwork, CollisionModel, Conditions, SteadyState,
    };
    use crate::masterequation::chemical_activation_operator::tests::{conditions, two_well_network};
    use crate::masterequation::direct_time_integration::{integrate_master_equation, log_spaced_times, InitialState, TimeIntegrationSettings};
    use crate::masterequation::chemical_activation_network::ChemicalActivationOptions;
    use crate::masterequation::prepared_distributions::{gaussian, thermal, EnergyReference, GrainDistribution, Representation};
    use crate::masterequation::source_profiles::TimeProfile;

    const MODEL: CollisionModel = CollisionModel::ExponentialDown { cutoff_in_mean_down: 10.0 };

    fn settings(times: Vec<f64>) -> TimeIntegrationSettings {
        TimeIntegrationSettings { relative_tolerance: 1e-8, absolute_tolerance: 1e-16, times_s: times, ..TimeIntegrationSettings::default() }
    }

    fn hot_source(network: &ChemicalActivationNetwork) -> GrainDistribution {
        gaussian(network, 0, 3300.0, 150.0, Representation::Density, EnergyReference::AboveWellGround).unwrap()
    }

    fn constant_bath() -> Vec<BathSegment> {
        vec![BathSegment { start_s: 0.0, conditions: conditions() }]
    }

    fn close(a: f64, b: f64, rel: f64, abs: f64) -> bool {
        (a - b).abs() <= rel * a.abs().max(b.abs()) + abs
    }

    fn assert_points_agree(a: &[TransientPoint], b: &[TransientPoint], rel: f64, abs: f64) {
        assert_eq!(a.len(), b.len());
        for (p, q) in a.iter().zip(b) {
            for (x, y) in p.well_populations.iter().zip(&q.well_populations).chain(p.exit_yields.iter().zip(&q.exit_yields)) {
                assert!(close(*x, *y, rel, abs), "t = {:e}: {x:e} vs {y:e}", p.time_s);
            }
        }
    }

    #[test]
    fn an_initial_population_equals_an_impulse_at_zero_and_the_pulse_of_the_existing_driver() {
        let network = two_well_network();
        let f = hot_source(&network);
        let times = log_spaced_times(1e-10, 1e-2, 2);
        let old = integrate_master_equation(
            &network,
            &conditions(),
            &ChemicalActivationOptions { collision_model: MODEL, steady_state: SteadyState::Final },
            &f.mass,
            InitialState::Pulse,
            &settings(times.clone()),
        )
        .unwrap();
        let initial = Preparation {
            initial: Some(InitialPopulation { amount: 1.0, distribution: f.clone() }),
            channels: Vec::new(),
            bath: constant_bath(),
        };
        let impulse = Preparation {
            initial: None,
            channels: vec![SourceChannel { name: "laser".into(), distribution: f.clone(), profile: TimeProfile::Impulse { time_s: 0.0, amount: 1.0 } }],
            bath: constant_bath(),
        };
        let a = integrate_preparation(&network, MODEL, SteadyState::Final, &initial, &settings(times.clone())).unwrap();
        let b = integrate_preparation(&network, MODEL, SteadyState::Final, &impulse, &settings(times.clone())).unwrap();
        assert_eq!(a.exits, old.exits);
        for (p, q) in a.points.iter().zip(&old.points) {
            for (x, y) in p.well_populations.iter().zip(&q.well_populations).chain(p.exit_yields.iter().zip(&q.exit_yields)) {
                assert!(close(*x, *y, 1e-12, 1e-18), "t = {:e}: {x:e} vs {y:e}", p.time_s);
            }
        }
        assert_points_agree(&a.points, &b.points, 1e-12, 1e-18);
        assert_eq!(b.points.last().unwrap().injected, vec![1.0]);
    }

    #[test]
    fn a_feed_from_zero_equals_the_continuous_formation_of_the_existing_driver() {
        let network = two_well_network();
        let f = hot_source(&network);
        let times = log_spaced_times(1e-9, 1e-3, 2);
        let old = integrate_master_equation(
            &network,
            &conditions(),
            &ChemicalActivationOptions { collision_model: MODEL, steady_state: SteadyState::Final },
            &f.mass,
            InitialState::ContinuousFormation,
            &settings(times.clone()),
        )
        .unwrap();
        let feed = Preparation {
            initial: None,
            channels: vec![SourceChannel { name: "feed".into(), distribution: f, profile: TimeProfile::Feed { start_s: 0.0, end_s: None, rate: 1.0 } }],
            bath: constant_bath(),
        };
        let a = integrate_preparation(&network, MODEL, SteadyState::Final, &feed, &settings(times)).unwrap();
        for (p, q) in a.points.iter().zip(&old.points) {
            for (x, y) in p.well_populations.iter().zip(&q.well_populations).chain(p.exit_yields.iter().zip(&q.exit_yields)) {
                assert!(close(*x, *y, 1e-12, 1e-24), "t = {:e}: {x:e} vs {y:e}", p.time_s);
            }
        }
    }

    #[test]
    fn independently_propagated_channels_add_up_to_their_simultaneous_propagation() {
        // N tab. validation, "Source superposition": the system is linear.
        let network = two_well_network();
        let hot = hot_source(&network);
        let warm = thermal(&network, 1, 900.0).unwrap();
        let a = SourceChannel { name: "a".into(), distribution: hot, profile: TimeProfile::Impulse { time_s: 2e-9, amount: 0.7 } };
        let b = SourceChannel {
            name: "b".into(),
            distribution: warm,
            profile: TimeProfile::Gaussian { centre_s: 5e-8, sigma_s: 2e-8, amount: 0.3 },
        };
        let times = log_spaced_times(1e-9, 1e-3, 3);
        let run = |channels: Vec<SourceChannel>| {
            integrate_preparation(&network, MODEL, SteadyState::Final, &Preparation { initial: None, channels, bath: constant_bath() }, &settings(times.clone()))
                .unwrap()
        };
        let both = run(vec![a.clone(), b.clone()]);
        let (only_a, only_b) = (run(vec![a]), run(vec![b]));
        for ((p, qa), qb) in both.points.iter().zip(&only_a.points).zip(&only_b.points) {
            for k in 0..p.well_populations.len() {
                assert!(close(p.well_populations[k], qa.well_populations[k] + qb.well_populations[k], 1e-6, 1e-14), "t = {:e}", p.time_s);
            }
            for k in 0..p.exit_yields.len() {
                assert!(close(p.exit_yields[k], qa.exit_yields[k] + qb.exit_yields[k], 1e-6, 1e-14), "t = {:e}", p.time_s);
            }
        }
    }

    #[test]
    fn a_narrowing_rectangular_pulse_approaches_the_impulse_with_the_same_amount() {
        // N tab. validation, "Finite-pulse limit".
        let network = two_well_network();
        let f = hot_source(&network);
        let times = vec![1e-7, 1e-6];
        let run = |profile: TimeProfile| {
            let prep = Preparation { initial: None, channels: vec![SourceChannel { name: "p".into(), distribution: f.clone(), profile }], bath: constant_bath() };
            integrate_preparation(&network, MODEL, SteadyState::Final, &prep, &settings(times.clone())).unwrap()
        };
        // To first order a rectangular pulse of width w is the impulse delayed by w/2: Y_rect(t) - Y_imp(t) =
        // -(w/2) dY_imp/dt, so the deviation falls in proportion to the width.
        let impulse = run(TimeProfile::Impulse { time_s: 0.0, amount: 1.0 });
        let mut deviations = Vec::new();
        for width in [1e-8, 1e-9, 1e-10] {
            let pulse = run(TimeProfile::Rectangular { start_s: 0.0, end_s: width, amount: 1.0 });
            assert!((pulse.points[1].injected[0] - 1.0).abs() < 1e-15);
            deviations.push(
                pulse.points[1].exit_yields.iter().zip(&impulse.points[1].exit_yields).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max),
            );
        }
        for pair in deviations.windows(2) {
            assert!((pair[0] / pair[1] / 10.0 - 1.0).abs() < 0.05, "{deviations:?}");
        }
    }

    #[test]
    fn the_population_balance_holds_with_impulses_profiles_and_an_initial_population() {
        // N eq. balance: sum n + sum Y = N_0 + sum_a N_a,in(0, t).
        let network = two_well_network();
        let prep = Preparation {
            initial: Some(InitialPopulation { amount: 0.5, distribution: thermal(&network, 0, 300.0).unwrap() }),
            channels: vec![
                SourceChannel {
                    name: "train".into(),
                    distribution: hot_source(&network),
                    profile: TimeProfile::Train(vec![
                        TimeProfile::Impulse { time_s: 1e-8, amount: 0.2 },
                        TimeProfile::Impulse { time_s: 3e-7, amount: 0.2 },
                        TimeProfile::Rectangular { start_s: 1e-6, end_s: 5e-6, amount: 0.1 },
                    ]),
                },
                SourceChannel {
                    name: "precursor".into(),
                    distribution: thermal(&network, 1, 1200.0).unwrap(),
                    profile: TimeProfile::PrecursorDecay { start_s: 0.0, formation_rate_s_inv: 1e5, total_loss_rate_s_inv: 2e5, precursor_amount: 1.0 },
                },
            ],
            bath: constant_bath(),
        };
        let result = integrate_preparation(&network, MODEL, SteadyState::Final, &prep, &settings(log_spaced_times(1e-9, 1e-2, 4))).unwrap();
        for p in &result.points {
            assert!(p.balance_deviation.abs() < 1e-7, "t = {:e}: {:e}", p.time_s, p.balance_deviation);
        }
        let last = result.points.last().unwrap();
        assert!((last.injected[0] - 0.5).abs() < 1e-15);
        assert!((last.injected[1] - 0.5 * (1.0 - (-2e5f64 * 1e-2).exp())).abs() < 1e-12);
    }

    #[test]
    fn an_absorbing_barrier_counts_source_mass_below_it_at_once() {
        // N Sec. 9.4: mass formed below the barrier is stabilized at once, as initial population and as an impulse.
        let network = two_well_network();
        // At 300 K the barrier of A lies at grain 41 (10 k_BT below its lowest threshold), above its reservoir.
        let cold = gaussian(&network, 0, 100.0, 30.0, Representation::Density, EnergyReference::AboveWellGround).unwrap();
        let barrier = SteadyState::Intermediate { barrier: AbsorbingBarrier::default() };
        let times = vec![1e-12];
        let initial = Preparation { initial: Some(InitialPopulation { amount: 1.0, distribution: cold.clone() }), channels: Vec::new(), bath: constant_bath() };
        let impulse = Preparation {
            initial: None,
            channels: vec![SourceChannel { name: "c".into(), distribution: cold, profile: TimeProfile::Impulse { time_s: 0.0, amount: 1.0 } }],
            bath: constant_bath(),
        };
        let a = integrate_preparation(&network, MODEL, barrier.clone(), &initial, &settings(times.clone())).unwrap();
        let b = integrate_preparation(&network, MODEL, barrier, &impulse, &settings(times)).unwrap();
        let stab = a.exits.iter().position(|e| e == "stab(A)").unwrap();
        let note = a.projections.iter().find(|n| n.what == "initial").unwrap();
        assert!(note.absorbed_fraction[0] > 0.99, "{note:?}");
        assert!((a.points[0].exit_yields[stab] - note.absorbed_fraction[0]).abs() < 1e-9, "{:?}", a.points[0].exit_yields);
        assert!((a.points[0].exit_yields[stab] - b.points[0].exit_yields[stab]).abs() < 1e-15);
    }

    #[test]
    fn a_closed_well_relaxes_to_its_bath_distribution_from_cold_and_hot_preparations() {
        // N tab. validation, "Collision-only relaxation": no losses; the population is conserved and relaxes to the
        // Boltzmann distribution of the bath, whatever the preparation temperature.
        let mut network = two_well_network();
        network.wells.truncate(1);
        network.wells[0].channels.clear();
        let bath = Conditions { temperature_kelvin: 1000.0, pressure_torr: 760.0 };
        let thermal_mean = thermal(&network, 0, 1000.0).unwrap().mean_energy_above_bottom_cm1(&network)[0];
        for t_prep in [300.0, 2000.0] {
            let prep = Preparation {
                initial: Some(InitialPopulation { amount: 1.0, distribution: thermal(&network, 0, t_prep).unwrap() }),
                channels: Vec::new(),
                bath: vec![BathSegment { start_s: 0.0, conditions: bath.clone() }],
            };
            let result = integrate_preparation(&network, MODEL, SteadyState::Final, &prep, &settings(log_spaced_times(1e-10, 1e-5, 2))).unwrap();
            let first = &result.points[0];
            let last = result.points.last().unwrap();
            assert!((last.well_populations[0] - 1.0).abs() < 1e-9);
            assert!(
                (last.mean_energy_above_bottom_cm1[0] / thermal_mean - 1.0).abs() < 1e-6,
                "T_prep {t_prep}: {} vs {thermal_mean}",
                last.mean_energy_above_bottom_cm1[0]
            );
            // It starts on the side of its preparation temperature.
            assert_eq!(first.mean_energy_above_bottom_cm1[0] < thermal_mean, t_prep < 1000.0);
        }
    }

    #[test]
    fn a_constant_history_split_into_segments_equals_one_segment_and_a_jump_conserves_population() {
        let network = two_well_network();
        let f = hot_source(&network);
        let times = log_spaced_times(1e-9, 1e-3, 2);
        let prep = |bath: Vec<BathSegment>| Preparation {
            initial: Some(InitialPopulation { amount: 1.0, distribution: f.clone() }),
            channels: Vec::new(),
            bath,
        };
        let one = integrate_preparation(&network, MODEL, SteadyState::Final, &prep(constant_bath()), &settings(times.clone())).unwrap();
        let split = integrate_preparation(
            &network,
            MODEL,
            SteadyState::Final,
            &prep(vec![BathSegment { start_s: 0.0, conditions: conditions() }, BathSegment { start_s: 3e-7, conditions: conditions() }]),
            &settings(times.clone()),
        )
        .unwrap();
        assert_points_agree(&one.points, &split.points, 1e-6, 1e-14);
        // A temperature jump 300 K -> 1500 K: the total is conserved through the remapping.
        let jump = integrate_preparation(
            &network,
            MODEL,
            SteadyState::Final,
            &prep(vec![
                BathSegment { start_s: 0.0, conditions: conditions() },
                BathSegment { start_s: 1e-6, conditions: Conditions { temperature_kelvin: 1500.0, pressure_torr: 760.0 } },
            ]),
            &settings(times),
        )
        .unwrap();
        for p in &jump.points {
            assert!(p.balance_deviation.abs() < 1e-7, "t = {:e}: {:e}", p.time_s, p.balance_deviation);
        }
        assert_eq!(jump.segments.len(), 2);
    }

    #[test]
    fn a_cold_preparation_in_a_hot_bath_reacts_late_and_a_hot_one_at_once() {
        // N Sec. 8.1 (shock heating) and the two-grain example (N Sec. toy): with the same operator, a cold
        // preparation must first be activated (no initial reactive flux), a hot one reacts at once; both end with the
        // same slow decay. One well with its threshold at 8000 cm-1, far above the low-energy reservoir of the 1500 K
        // bath (which reaches about 2350 cm-1 with this density of states).
        let mut network = two_well_network();
        network.wells.truncate(1);
        network.wells[0] = crate::masterequation::chemical_activation_network::tests::test_well("A", 1000, 0, 800);
        network.wells[0].density_of_states = (0..1000).map(|i| (1.0 + 0.02 * i as f64).powi(8)).collect();
        let hot_bath = vec![BathSegment { start_s: 0.0, conditions: Conditions { temperature_kelvin: 1500.0, pressure_torr: 760.0 } }];
        // The first output at 1e-13 s, about 1e-3 collisions (1e-11 s would already be 0.1 collisions at 760 Torr).
        let times = log_spaced_times(1e-13, 1e-5, 2);
        let run = |distribution: GrainDistribution| {
            let prep = Preparation { initial: Some(InitialPopulation { amount: 1.0, distribution }), channels: Vec::new(), bath: hot_bath.clone() };
            integrate_preparation(&network, MODEL, SteadyState::Final, &prep, &settings(times.clone())).unwrap()
        };
        let cold = run(thermal(&network, 0, 300.0).unwrap());
        let hot = run(gaussian(&network, 0, 9000.0, 200.0, Representation::Density, EnergyReference::AboveWellGround).unwrap());
        let products = |r: &TransientResult, k: usize| r.points[k].exit_yields.iter().sum::<f64>();
        assert!(products(&hot, 0) > 1e3 * products(&cold, 0), "{:e} vs {:e}", products(&hot, 0), products(&cold, 0));
        // The cold preparation has 4.7e-10 of its 300 K population above the threshold. After 1e-13 s (about 1e-3
        // collisions) its tail is 5e-6: the part of it in the reservoir state carries the bath's Boltzmann shape and is
        // activated at once at 1500 K rates (N, keypoint of Sec. 8). Still far below the hot preparation.
        let (tc, th) = (cold.points[0].tail_fraction[0], hot.points[0].tail_fraction[0]);
        assert!(tc < 1e-4 * th && th > 0.5, "{tc:e} {th:e}");
        // The same operator: once relaxed, both decay with the same hazard (the slowest eigenvalue).
        let k = times.iter().position(|&t| (t - 1e-6).abs() < 1e-12).unwrap();
        let (kc, kh) = (cold.points[k].loss_hazard_s_inv, hot.points[k].loss_hazard_s_inv);
        assert!((kc / kh - 1.0).abs() < 1e-3, "{kc:e} vs {kh:e}");
        // The cold preparation falls largely into the reservoir of the hot bath, which replaces its shape by the
        // bath's Boltzmann shape there (N, keypoint of Sec. 8): the projection reports it.
        let note = cold.projections.iter().find(|n| n.what == "initial").unwrap();
        assert!(note.reservoir_fraction[0] > 0.5, "{note:?}");
        assert!(note.mean_energy_after_cm1[0] > note.mean_energy_before_cm1[0], "{note:?}");
    }

    #[test]
    fn two_grains_follow_the_analytical_cold_and_hot_solutions() {
        // N Sec. toy, eqs. toy-toyhot: grains L and H, L -> H at rate a, H -> L at rate b, H -> products at rate k:
        //   d/dt (n_L, n_H) = -[[a, -b], [-a, b + k]] (n_L, n_H),  dY/dt = k n_H,
        //   lambda_pm = (a + b + k +- sqrt((a + b + k)^2 - 4 a k))/2.
        // Cold start (n_L = 1): Y = a k/(l+ - l-) [phi(l-) - phi(l+)], phi(l) = (1 - e^{-l t})/l, Y ~ a k t^2/2. Hot start
        // (n_H = 1): n_H = c+ e^{-l+ t} + c- e^{-l- t} with c+ = (b + k - l-)/(l+ - l-), c- = (l+ - b - k)/(l+ - l-),
        // Y = k [c+ phi(l+) + c- phi(l-)], Y ~ k t. a and b are those of the assembled collision operator.
        use crate::masterequation::chemical_activation_network::tests::test_well;
        use crate::masterequation::chemical_activation_operator::assemble_operator;
        let mut well = test_well("A", 2, 0, 1);
        well.density_of_states = vec![1.0, 1.5];
        let network = ChemicalActivationNetwork { grain_width_cm1: 100.0, wells: vec![well] };
        let bath = Conditions { temperature_kelvin: 300.0, pressure_torr: 0.1 };
        let op = assemble_operator(&network, &bath, &ChemicalActivationOptions { collision_model: MODEL, steady_state: SteadyState::Final }).unwrap();
        assert!(op.reservoirs[0].is_none() && op.dimension() == 2);
        let entry = |r: usize, c: usize| op.rows[r].iter().find(|&&(col, _)| col == c).map_or(0.0, |&(_, v)| v);
        let k = network.wells[0].channels[0].rate_constant_s_inv[1];
        let (a, b) = (entry(0, 0), -entry(0, 1));
        assert!((entry(1, 0) + a).abs() < 1e-9 * a && (entry(1, 1) - (b + k)).abs() < 1e-9 * (b + k));
        let sum = a + b + k;
        let root = (sum * sum - 4.0 * a * k).sqrt();
        let (lp, lm) = (0.5 * (sum + root), 0.5 * (sum - root));
        let phi = |l: f64, t: f64| -(-l * t).exp_m1() / l;
        let y_cold = |t: f64| a * k / (lp - lm) * (phi(lm, t) - phi(lp, t));
        let (cp, cm) = ((b + k - lm) / (lp - lm), (lp - b - k) / (lp - lm));
        let y_hot = |t: f64| k * (cp * phi(lp, t) + cm * phi(lm, t));
        let times = log_spaced_times(1e-3 / sum, 30.0 / lm, 4);
        let run = |mass: Vec<f64>| {
            let distribution = GrainDistribution { mass: vec![mass], lost_fraction: 0.0, description: "grain".into() };
            let prep = Preparation { initial: Some(InitialPopulation { amount: 1.0, distribution }), channels: Vec::new(), bath: vec![BathSegment { start_s: 0.0, conditions: bath.clone() }] };
            integrate_preparation(&network, MODEL, SteadyState::Final, &prep, &settings(times.clone())).unwrap()
        };
        let (cold, hot) = (run(vec![1.0, 0.0]), run(vec![0.0, 1.0]));
        for (p, q) in cold.points.iter().zip(&hot.points) {
            let t = p.time_s;
            assert!((p.exit_yields[0] / y_cold(t) - 1.0).abs() < 1e-5, "cold t = {t:e}: {:e} vs {:e}", p.exit_yields[0], y_cold(t));
            assert!((q.exit_yields[0] / y_hot(t) - 1.0).abs() < 1e-5, "hot t = {t:e}: {:e} vs {:e}", q.exit_yields[0], y_hot(t));
        }
        // Short times: the cold start forms products quadratically, the hot start linearly.
        let t0 = times[0];
        assert!((cold.points[0].exit_yields[0] / (0.5 * a * k * t0 * t0) - 1.0).abs() < 2e-3);
        assert!((hot.points[0].exit_yields[0] / (k * t0) - 1.0).abs() < 2e-3);
        // Both end with unit yield.
        assert!((cold.points.last().unwrap().exit_yields[0] - 1.0).abs() < 1e-6 && (hot.points.last().unwrap().exit_yields[0] - 1.0).abs() < 1e-6);
    }


    #[test]
    fn the_steady_state_shape_is_the_feed_shape_or_the_amount_weighted_injected_shape() {
        let network = two_well_network();
        let (a, b) = (thermal(&network, 0, 300.0).unwrap(), hot_source(&network));
        let impulse = SourceChannel { name: "i".into(), distribution: b.clone(), profile: TimeProfile::Impulse { time_s: 1e-9, amount: 3.0 } };
        let prep = Preparation { initial: Some(InitialPopulation { amount: 1.0, distribution: a.clone() }), channels: vec![impulse.clone()], bath: constant_bath() };
        let (shape, how) = steady_source_shape(&prep).unwrap();
        assert!(how.contains("amount"), "{how}");
        for (x, (y, z)) in shape.mass[0].iter().zip(a.mass[0].iter().zip(&b.mass[0])) {
            assert!((x - (0.25 * y + 0.75 * z)).abs() < 1e-14);
        }
        // An open-ended feed defines the steady source alone.
        let feed = SourceChannel { name: "f".into(), distribution: a.clone(), profile: TimeProfile::Feed { start_s: 0.0, end_s: None, rate: 2.0 } };
        let prep = Preparation { initial: None, channels: vec![impulse, feed], bath: constant_bath() };
        let (shape, how) = steady_source_shape(&prep).unwrap();
        assert!(how.contains("'f'"), "{how}");
        assert!(shape.mass[0].iter().zip(&a.mass[0]).all(|(x, y)| (x - y).abs() < 1e-14));
    }
}
