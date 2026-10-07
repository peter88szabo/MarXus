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
use super::chemical_activation_eigen::lowest_eigenvector_grain_populations;
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
    /// Flux coefficient r_wc = sum_i k_wc(E_i) n_wi / C_w of every channel of every well (in the order of
    /// `Well::channels`: products and isomerization), s-1 (BFG 2015 eq. A5; Miller et al. 2016 eq. 4). NaN for an empty
    /// well.
    pub flux_coefficients_s_inv: Vec<Vec<f64>>,
    /// R_w = sum_c r_wc + k_c[D]_w, s-1: all reactive loss out of the well per molecule in it (BFG 2015 eq. 1).
    pub loss_coefficients_s_inv: Vec<f64>,
    /// Standard deviation of the energy of every well (cm-1).
    pub energy_spread_cm1: Vec<f64>,
    /// d<E>/dt of every well (cm-1 s-1), exact from dn/dt = -J n + sum_a R_a(t) F_a (continuous sources included).
    pub mean_energy_rate_cm1_s: Vec<f64>,
    /// E_f: the mean energy above the bottom of every well in the lowest eigenvector of J of the bath segment, the final
    /// steady state (Barker, King 1995, eq. 11); the Boltzmann mean for a closed well. NaN in the intermediate variant.
    pub final_mean_energy_cm1: Vec<f64>,
    /// tau_vib = -(<E> - E_f)/(d<E>/dt) of every well, s (Barker, King 1995, eq. 11); negative while <E> moves away
    /// from E_f; NaN once |<E> - E_f| <= 10 rtol E_f (rtol the relative tolerance of the integration), where both are
    /// rounding and integration error.
    pub vibrational_relaxation_time_s: Vec<f64>,
    /// t_inc from the start s of the bath segment: (t - s) + ln(N(t)/N_ref)/k_inst(t), N_ref the population after the
    /// start plus the impulses since. Its late-time plateau is the back-extrapolation of the first-order decay to
    /// N/N_ref = 1 (Barker, King 1995, p. 4960 and eq. 9). NaN once a continuous source has injected in the segment.
    pub incubation_time_s: f64,
    /// Z_w = integral of omega_w dt from t = 0: the time in collisions of every well (Eng et al. 2001, Fig. 10).
    pub collision_numbers: Vec<f64>,
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
    /// Diagnostics that could not be computed (E_f of a segment), with the reason.
    pub warnings: Vec<String>,
}

/// Projection of `amount` x `distribution` onto the states of `op`: (on states, absorbed per well, note).
pub(crate) fn project(
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
    /// Start of the segment (s).
    start_s: f64,
    /// E_f of every well (`TransientPoint::final_mean_energy_cm1`).
    final_mean_energy: Vec<f64>,
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
    let mut warnings = Vec::new();
    let build_segment = |k: usize, projections: &mut Vec<ProjectionNote>, warnings: &mut Vec<String>| -> Result<Segment, String> {
        let conditions = &bath[k].conditions;
        let op = assemble_operator(network, conditions, &options_of(conditions))?;
        let (exits, exit_rates) = exit_definitions(network, &op, absorbing_barrier);
        let mut shapes = Vec::new();
        for channel in &preparation.channels {
            let (on_states, absorbed, note) = project(network, &op, &channel.distribution, 1.0, &channel.name, k)?;
            projections.push(note);
            shapes.push((on_states, absorbed));
        }
        let final_mean_energy = if absorbing_barrier {
            vec![f64::NAN; network.wells.len()]
        } else {
            match lowest_eigenvector_grain_populations(&op) {
                Ok(grains) => grains.iter().map(|g| mean_and_spread(g, network.grain_width_cm1).0).collect(),
                Err(e) => {
                    warnings.push(format!("E_f of bath segment {} not available: {e}", k + 1));
                    vec![f64::NAN; network.wells.len()]
                }
            }
        };
        Ok(Segment { op, exits, exit_rates, shapes, start_s: bath[k].start_s, final_mean_energy })
    };

    // Global exit list (by name) and the population state carried between segments, on the grains.
    let mut segment_index = 0;
    let mut segment = build_segment(0, &mut projections, &mut warnings)?;
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
    // Applies the impulses at t; returns the amount put on the states.
    let apply_impulses = |t: f64, segment: &Segment, y_states: &mut Vec<f64>, yields: &mut Vec<f64>| -> Result<f64, String> {
        let mut added = 0.0;
        for (channel, (on_states, absorbed)) in preparation.channels.iter().zip(&segment.shapes) {
            for (tp, amount) in channel.profile.impulses() {
                if tp == t {
                    y_states.iter_mut().zip(on_states).for_each(|(y, x)| *y += amount * x);
                    added += amount * on_states.iter().sum::<f64>();
                    add_absorbed(yields, absorbed, amount)?;
                }
            }
        }
        Ok(added)
    };
    apply_impulses(0.0, &segment, &mut y_states, &mut yields)?;
    // N_ref of the incubation time and the collision numbers at the start of the current segment.
    let mut reference_population: f64 = y_states.iter().sum();
    let mut collisions_at_start = vec![0.0; network.wells.len()];

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
        let switched = bath.iter().position(|b| b.start_s == t);
        if let Some(k) = switched {
            for (z, w) in collisions_at_start.iter_mut().zip(&segment.op.wells) {
                *z += w.collision_frequency_s_inv * (t - segment.start_s);
            }
            let grains = segment.op.grain_populations(&y_states);
            segment_index = k;
            segment = build_segment(k, &mut projections, &mut warnings)?;
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
        let added = apply_impulses(t, &segment, &mut y_states, &mut yields)?;
        if switched.is_some() {
            reference_population = y_states.iter().sum();
        } else {
            reference_population += added;
        }
        if times.contains(&t) {
            let history = History {
                reference_population,
                collisions_at_start: &collisions_at_start,
                energy_resolution: 10.0 * settings.relative_tolerance,
            };
            points.push(observe(network, &segment, &exits, &y_states, &yields, preparation, amount_in, t, &history));
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
        warnings,
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

/// What the observables need from the integration so far: N_ref of the incubation time and the collision numbers at the
/// start of the current segment.
struct History<'a> {
    reference_population: f64,
    collisions_at_start: &'a [f64],
    /// |<E> - E_f| below this times E_f is not resolved (10 times the relative tolerance of the integration).
    energy_resolution: f64,
}

/// Mean energy above the bottom and its standard deviation of populations on the grains of one well (NaN if empty).
fn mean_and_spread(grains: &[f64], grain_width_cm1: f64) -> (f64, f64) {
    let total: f64 = grains.iter().sum();
    if !(total > 0.0) {
        return (f64::NAN, f64::NAN);
    }
    let mean = grains.iter().enumerate().map(|(i, x)| i as f64 * grain_width_cm1 * x).sum::<f64>() / total;
    let variance = grains.iter().enumerate().map(|(i, x)| (i as f64 * grain_width_cm1 - mean).powi(2) * x).sum::<f64>() / total;
    (mean, variance.max(0.0).sqrt())
}

/// The observables of N Sec. observables and the diagnostics of reports/nonthermal_sources_design.md, Section 12, at
/// time `t`.
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
    history: &History,
) -> TransientPoint {
    let op = &segment.op;
    let n_wells = network.wells.len();
    let mut well_populations = vec![0.0; n_wells];
    for (s, &(w, _)) in op.states.iter().enumerate() {
        well_populations[w] += y_states[s];
    }
    let mut exit_fluxes = vec![0.0; exits.len()];
    for (x, rates) in segment.exit_rates.iter().enumerate() {
        let g = exits.iter().position(|e| *e == segment.exits[x]).unwrap();
        exit_fluxes[g] = rates.iter().map(|&(s, k)| k * y_states[s]).sum();
    }
    let grains = op.grain_populations(y_states);
    // dn/dt = -J n + sum_a R_a(t) F_a, on the grains.
    let mut rate_of_change: Vec<f64> = op.apply(y_states).iter().map(|x| -x).collect();
    for (channel, (on_states, _)) in preparation.channels.iter().zip(&segment.shapes) {
        let r = channel.profile.rate(t);
        if r != 0.0 {
            rate_of_change.iter_mut().zip(on_states).for_each(|(d, f)| *d += r * f);
        }
    }
    let grains_rate = op.grain_populations(&rate_of_change);
    let de = network.grain_width_cm1;
    let (mut mean_energy, mut spread, mut tail, mut energy_rate, mut tau_vib) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let (mut flux_coefficients, mut loss_coefficients) = (Vec::new(), Vec::new());
    for (w, g) in grains.iter().enumerate() {
        let well = &network.wells[w];
        let total: f64 = g.iter().sum();
        let (mean, sd) = mean_and_spread(g, de);
        mean_energy.push(mean);
        spread.push(sd);
        if total > 0.0 {
            let threshold = well.lowest_threshold_grain().unwrap_or(g.len());
            tail.push(g.iter().skip(threshold).sum::<f64>() / total);
            let r: Vec<f64> = well.channels.iter().map(|c| c.rate_constant_s_inv.iter().zip(g).map(|(k, x)| k * x).sum::<f64>() / total).collect();
            loss_coefficients.push(r.iter().sum::<f64>() + well.bimolecular_sink_s_inv);
            flux_coefficients.push(r);
            let d_total: f64 = grains_rate[w].iter().sum();
            let d_energy: f64 = grains_rate[w].iter().enumerate().map(|(i, x)| i as f64 * de * x).sum();
            let rate = (d_energy - mean * d_total) / total;
            energy_rate.push(rate);
            // Once <E> is at E_f within the resolution, both <E> - E_f and d<E>/dt are rounding and integration error.
            let excess = mean - segment.final_mean_energy[w];
            tau_vib.push(if excess.abs() > history.energy_resolution * segment.final_mean_energy[w].abs() { -excess / rate } else { f64::NAN });
        } else {
            tail.push(f64::NAN);
            loss_coefficients.push(f64::NAN);
            flux_coefficients.push(vec![f64::NAN; well.channels.len()]);
            energy_rate.push(f64::NAN);
            tau_vib.push(f64::NAN);
        }
    }
    let injected: Vec<f64> = preparation.channels.iter().map(|c| c.profile.injected_until(t)).collect();
    let total_in = initial_amount + injected.iter().sum::<f64>();
    let present: f64 = y_states.iter().sum::<f64>() + yields.iter().sum::<f64>();
    let population: f64 = well_populations.iter().sum();
    let loss_hazard = if population > 0.0 { exit_fluxes.iter().sum::<f64>() / population } else { f64::NAN };
    // Continuous injection since the segment start: the injected amount less the impulses in (start, t].
    let start = segment.start_s;
    let continuous: f64 = preparation
        .channels
        .iter()
        .map(|c| {
            let impulses: f64 = c.profile.impulses().iter().filter(|&&(tp, _)| tp > start && tp <= t).map(|&(_, a)| a).sum();
            c.profile.injected_until(t) - c.profile.injected_until(start) - impulses
        })
        .sum();
    let incubation_time = if continuous > 1e-12 * total_in.max(f64::MIN_POSITIVE) || !(loss_hazard > 0.0) || !(history.reference_population > 0.0) {
        f64::NAN
    } else {
        (t - start) + (population / history.reference_population).ln() / loss_hazard
    };
    let collision_numbers = op
        .wells
        .iter()
        .zip(history.collisions_at_start)
        .map(|(w, z)| z + w.collision_frequency_s_inv * (t - start))
        .collect();
    TransientPoint {
        time_s: t,
        well_populations,
        exit_yields: yields.to_vec(),
        loss_hazard_s_inv: loss_hazard,
        exit_fluxes,
        injected,
        mean_energy_above_bottom_cm1: mean_energy,
        tail_fraction: tail,
        balance_deviation: if total_in > 0.0 { (present - total_in) / total_in } else { present },
        flux_coefficients_s_inv: flux_coefficients,
        loss_coefficients_s_inv: loss_coefficients,
        energy_spread_cm1: spread,
        mean_energy_rate_cm1_s: energy_rate,
        final_mean_energy_cm1: segment.final_mean_energy.clone(),
        vibrational_relaxation_time_s: tau_vib,
        incubation_time_s: incubation_time,
        collision_numbers,
    }
}

/// The effective coefficient of BFG 2015 eq. A7, k^e(t_0, t_k) = (t_k - t_0)^-1 integral_{t_0}^{t_k} r dt, of a series r at
/// the output times, for every t_k after t_0 = the first output time (NaN at t_0, and from the first NaN of r on). Between
/// two output times r is taken linear in u = ln t, and integral r dt = integral r e^u du is integrated exactly:
/// integral_{t_k}^{t_k+1} r dt = r_k (t_k+1 - t_k) + (r_k+1 - r_k)/du ((du - 1) t_k+1 + t_k). Exact for a constant r and for
/// r linear in ln t, second order in the output spacing otherwise.
pub fn effective_coefficient(times: &[f64], r: &[f64]) -> Vec<f64> {
    let mut out = vec![f64::NAN; times.len()];
    let mut integral = 0.0;
    for k in 1..times.len() {
        let (t0, t1) = (times[k - 1], times[k]);
        let du = (t1 / t0).ln();
        integral += r[k - 1] * (t1 - t0) + (r[k] - r[k - 1]) / du * ((du - 1.0) * t1 + t0);
        out[k] = integral / (t1 - times[0]);
    }
    out
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
    use crate::masterequation::chemical_activation_network::ChannelDestination;
    use crate::masterequation::chemical_activation_eigen::{thermal_rate_coefficients, EigenSolver, EigenSystem};

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

    /// Well A of `two_well_network` alone, with its product channel (threshold at grain 300).
    fn one_reactive_well() -> ChemicalActivationNetwork {
        let mut network = two_well_network();
        network.wells.truncate(1);
        network.wells[0].channels.retain(|c| matches!(c.destination, ChannelDestination::Products { .. }));
        network
    }

    fn final_operator(network: &ChemicalActivationNetwork, conditions: &Conditions) -> ChemicalActivationOperator {
        assemble_operator(network, conditions, &ChemicalActivationOptions { collision_model: MODEL, steady_state: SteadyState::Final }).unwrap()
    }

    fn pulse(distribution: GrainDistribution, bath: Vec<BathSegment>) -> Preparation {
        Preparation { initial: Some(InitialPopulation { amount: 1.0, distribution }), channels: Vec::new(), bath }
    }

    #[test]
    fn flux_coefficients_add_up_to_the_exit_fluxes_and_end_at_the_thermal_rate_coefficients() {
        // BFG eqs. A5 and 1, Comment eq. 4: r_wc = sum_i k_wc(E_i) n_wi / C_w, R_w = sum_c r_wc + k_c[D]_w. Once only
        // the slowest mode is left the populations are the thermal eigenvector, whose average of k_wc is the thermal
        // rate coefficient (GO10, text after eq. 12; per unit of the whole population there, so divided by the well
        // fraction here). The lowest eigenvector also gives E_f (Barker, King 1995, eq. 11).
        // The second mode dies out as exp(-(lambda_2 - lambda_1) t) relative to the first; at t_late its share is
        // e^-15, while the population, exp(-lambda_1 t), stays far above the absolute tolerance. (This network has
        // little separation: lambda_2/lambda_1 = 2.1 at 1000 K.)
        let mut network = two_well_network();
        network.wells[1].bimolecular_sink_s_inv = 1e3;
        let bath = Conditions { temperature_kelvin: 1000.0, pressure_torr: 760.0 };
        let op = final_operator(&network, &bath);
        let thermal_eigen = thermal_rate_coefficients(&network, &op, EigenSolver::FullDecomposition, 1.0).unwrap();
        let (lambda_1, lambda_2) = (thermal_eigen.lambda_1_s_inv, thermal_eigen.lambda_2_s_inv);
        let t_late = 15.0 / (lambda_2 - lambda_1);
        assert!(lambda_1 * t_late < 20.0, "{lambda_1:e} {lambda_2:e}");
        let mut times = log_spaced_times(1e-14, 0.3 * t_late, 2);
        times.push(t_late);
        // Start: the Boltzmann distribution of the bath in both wells, with the well fractions of the thermal eigenvector
        // (so that the second mode starts small). r_wc then goes from the Boltzmann average of k_wc (the high-pressure
        // value) to the eigenvector average (the fall-off value).
        let fractions = &thermal_eigen.population_fractions;
        let start = crate::masterequation::prepared_distributions::mixture(&[
            (fractions[0], thermal(&network, 0, 1000.0).unwrap()),
            (fractions[1], thermal(&network, 1, 1000.0).unwrap()),
        ])
        .unwrap();
        let prep = pulse(start, vec![BathSegment { start_s: 0.0, conditions: bath }]);
        let result = integrate_preparation(&network, MODEL, SteadyState::Final, &prep, &settings(times)).unwrap();
        let first = &result.points[0];
        for c in &thermal_eigen.channels {
            let got = first.flux_coefficients_s_inv[c.well][c.channel];
            assert!(close(got, c.high_pressure_rate_s_inv, 1e-3, 0.0), "{} at 1e-14 s: {got:e} vs {:e}", c.name, c.high_pressure_rate_s_inv);
        }
        for p in &result.points {
            let mut from_coefficients = 0.0;
            for (w, well) in network.wells.iter().enumerate() {
                let products: f64 = well
                    .channels
                    .iter()
                    .zip(&p.flux_coefficients_s_inv[w])
                    .filter(|(c, _)| matches!(c.destination, ChannelDestination::Products { .. }))
                    .map(|(_, r)| r)
                    .sum();
                from_coefficients += p.well_populations[w] * (products + well.bimolecular_sink_s_inv);
                let total = p.flux_coefficients_s_inv[w].iter().sum::<f64>() + well.bimolecular_sink_s_inv;
                assert!(close(p.loss_coefficients_s_inv[w], total, 1e-14, 0.0), "t = {:e}", p.time_s);
            }
            let fluxes: f64 = p.exit_fluxes.iter().sum();
            assert!(close(from_coefficients, fluxes, 1e-10, 0.0), "t = {:e}: {from_coefficients:e} vs {fluxes:e}", p.time_s);
        }
        let last = result.points.last().unwrap();
        for c in &thermal_eigen.channels {
            let expected = c.thermal_rate_s_inv / thermal_eigen.population_fractions[c.well];
            let got = last.flux_coefficients_s_inv[c.well][c.channel];
            assert!(close(got, expected, 1e-5, 0.0), "{}: {got:e} vs {expected:e}", c.name);
        }
        for (w, d) in thermal_eigen.distributions.iter().enumerate() {
            let mean: f64 = d.iter().enumerate().map(|(i, x)| i as f64 * network.grain_width_cm1 * x).sum();
            assert!(close(last.final_mean_energy_cm1[w], mean, 1e-8, 0.0), "well {w}: {} vs {mean}", last.final_mean_energy_cm1[w]);
        }
    }

    #[test]
    fn the_mean_energy_rate_is_the_time_derivative_of_the_mean_energy_also_with_a_source() {
        // d<E>/dt = (sum E n' - <E> sum n')/C with n' = -J n + sum_a R_a(t) F_a, against a central difference.
        let mut network = two_well_network();
        network.wells.truncate(1);
        network.wells[0].channels.clear();
        let bath = vec![BathSegment { start_s: 0.0, conditions: Conditions { temperature_kelvin: 1000.0, pressure_torr: 760.0 } }];
        let (t, h) = (1e-8, 1e-11);
        let times = vec![t - h, t, t + h];
        let mut prep = pulse(thermal(&network, 0, 300.0).unwrap(), bath);
        prep.channels.push(SourceChannel { name: "feed".into(), distribution: hot_source(&network), profile: TimeProfile::Feed { start_s: 0.0, end_s: None, rate: 1e7 } });
        let r = integrate_preparation(&network, MODEL, SteadyState::Final, &prep, &settings(times)).unwrap();
        let difference = (r.points[2].mean_energy_above_bottom_cm1[0] - r.points[0].mean_energy_above_bottom_cm1[0]) / (2.0 * h);
        let exact = r.points[1].mean_energy_rate_cm1_s[0];
        assert!(close(exact, difference, 1e-5, 0.0), "{exact:e} vs {difference:e}");
        assert!(r.points[1].energy_spread_cm1[0] > 0.0);
    }

    #[test]
    fn in_a_closed_well_e_f_is_the_bath_mean_and_tau_vib_ends_at_the_slowest_relaxation_time() {
        // Barker, King 1995, eqs. 11-12: dE/dt = -(E - E_f)/tau_vib. Without reaction E_f is the Boltzmann mean of the
        // bath; once only the slowest relaxation mode is left, E - E_f decays as exp(-lambda_2 t), so tau_vib = 1/lambda_2.
        let mut network = two_well_network();
        network.wells.truncate(1);
        network.wells[0].channels.clear();
        let hot = Conditions { temperature_kelvin: 1000.0, pressure_torr: 760.0 };
        let eigenvalues = EigenSystem::new(&final_operator(&network, &hot), EigenSolver::FullDecomposition).unwrap().eigenvalues;
        let (lambda_2, lambda_3) = (eigenvalues[1], eigenvalues[2]);
        let bath_mean = thermal(&network, 0, 1000.0).unwrap().mean_energy_above_bottom_cm1(&network)[0];
        let prep = pulse(thermal(&network, 0, 300.0).unwrap(), vec![BathSegment { start_s: 0.0, conditions: hot }]);
        let times: Vec<f64> = [0.1, 4.0, 8.0, 12.0].iter().map(|x| x / lambda_2).collect();
        let precise = TimeIntegrationSettings { relative_tolerance: 1e-11, absolute_tolerance: 1e-20, times_s: times, ..TimeIntegrationSettings::default() };
        let r = integrate_preparation(&network, MODEL, SteadyState::Final, &prep, &precise).unwrap();
        let last = r.points.last().unwrap();
        assert!(close(last.final_mean_energy_cm1[0], bath_mean, 1e-9, 0.0), "{} vs {bath_mean}", last.final_mean_energy_cm1[0]);
        assert!(r.points[0].vibrational_relaxation_time_s[0] > 0.0);
        // The other modes die out as exp(-(lambda_k - lambda_2) t): the error of tau_vib lambda_2 falls with t.
        let errors: Vec<f64> = r.points[1..].iter().map(|p| (p.vibrational_relaxation_time_s[0] * lambda_2 - 1.0).abs()).collect();
        assert!(errors.windows(2).all(|e| e[1] < e[0]) && errors[2] < 1e-3, "{errors:?}; lambda_3/lambda_2 = {}", lambda_3 / lambda_2);
        // No reaction: no incubation time.
        assert!(last.incubation_time_s.is_nan());
    }

    #[test]
    fn the_incubation_time_reaches_the_back_extrapolation_of_the_long_time_decay() {
        // Barker, King 1995, p. 4960 and eq. 9: t_inc is where the extrapolated linear long-time ln(N/N_0) is 0. With
        // N(t) -> a_1 N_0 exp(-lambda_1 t) at late times, t_inc = ln(a_1)/lambda_1; a_1 from a full decomposition.
        let network = one_reactive_well();
        let hot = Conditions { temperature_kelvin: 1000.0, pressure_torr: 760.0 };
        let op = final_operator(&network, &hot);
        let cold = thermal(&network, 0, 300.0).unwrap();
        let symmetrized = crate::masterequation::chemical_activation_steady_state::symmetrize(&op);
        let (values, vectors) = crate::masterequation::chemical_activation_eigen::full_decomposition(&symmetrized.dense(), EigenSolver::FullDecomposition).unwrap();
        let n0 = op.state_populations(&cold.mass);
        let c1: f64 = vectors[0].iter().zip(&n0).zip(&symmetrized.d).map(|((u, n), d)| u * n / d).sum();
        let a1 = c1 * vectors[0].iter().zip(&symmetrized.d).map(|(u, d)| u * d).sum::<f64>();
        let (lambda_1, lambda_2) = (values[0], values[1]);
        let expected = a1.ln() / lambda_1;
        let late = [20.0 / lambda_2, 30.0 / lambda_2];
        assert!(lambda_1 * late[1] < 1.0, "{lambda_1:e} {lambda_2:e}");
        let prep = pulse(cold, vec![BathSegment { start_s: 0.0, conditions: hot }]);
        let r = integrate_preparation(&network, MODEL, SteadyState::Final, &prep, &settings(vec![0.01 / lambda_2, late[0], late[1]])).unwrap();
        for p in &r.points[1..] {
            assert!(close(p.incubation_time_s, expected, 1e-4, 0.0), "t = {:e}: {:e} vs {expected:e}", p.time_s, p.incubation_time_s);
        }
        assert!(expected > 0.0, "a cold start is delayed: {expected:e}");
        // While a continuous source injects, the extrapolation has no reference: NaN.
        let mut fed = pulse(thermal(&network, 0, 300.0).unwrap(), vec![BathSegment { start_s: 0.0, conditions: Conditions { temperature_kelvin: 1000.0, pressure_torr: 760.0 } }]);
        fed.channels.push(SourceChannel { name: "feed".into(), distribution: hot_source(&network), profile: TimeProfile::Feed { start_s: 0.0, end_s: None, rate: 1.0 } });
        let f = integrate_preparation(&network, MODEL, SteadyState::Final, &fed, &settings(vec![late[0]])).unwrap();
        assert!(f.points[0].incubation_time_s.is_nan());
    }

    #[test]
    fn collision_numbers_add_up_over_the_bath_segments() {
        // Z_w(t) = integral of omega_w dt; Eng et al. (2001) give incubation as Z_LJ[M] dt_inc.
        let network = two_well_network();
        let (first, second) = (conditions(), Conditions { temperature_kelvin: 300.0, pressure_torr: 76.0 });
        let omega = |c: &Conditions| -> Vec<f64> { final_operator(&network, c).wells.iter().map(|w| w.collision_frequency_s_inv).collect() };
        let (o1, o2) = (omega(&first), omega(&second));
        let t1 = 1e-7;
        let prep = pulse(hot_source(&network), vec![BathSegment { start_s: 0.0, conditions: first }, BathSegment { start_s: t1, conditions: second }]);
        let r = integrate_preparation(&network, MODEL, SteadyState::Final, &prep, &settings(vec![5e-8, 1e-6])).unwrap();
        for w in 0..2 {
            assert!(close(r.points[0].collision_numbers[w], o1[w] * 5e-8, 1e-12, 0.0));
            assert!(close(r.points[1].collision_numbers[w], o1[w] * t1 + o2[w] * (1e-6 - t1), 1e-12, 0.0));
        }
    }

    #[test]
    fn effective_coefficients_are_exact_for_a_constant_rate_and_converge_with_the_output_density() {
        // BFG 2015 eq. A7: k^e(t_0, t) = (t - t_0)^-1 integral_{t_0}^{t} r dt, t_0 the first output time.
        let times = log_spaced_times(1e-9, 1e-3, 2);
        let constant = effective_coefficient(&times, &vec![7.5; times.len()]);
        assert!(constant[0].is_nan());
        assert!(constant[1..].iter().all(|k| close(*k, 7.5, 1e-14, 0.0)), "{constant:?}");
        // r = a + b exp(-c t): (t - t_0) k^e = a (t - t_0) + (b/c)(exp(-c t_0) - exp(-c t)).
        let (a, b, c) = (2.0, 50.0, 1e6);
        let exact = |t0: f64, t: f64| a + b / c * ((-c * t0).exp() - (-c * t).exp()) / (t - t0);
        let mut errors = Vec::new();
        for per_decade in [2, 4, 8, 16] {
            let times = log_spaced_times(1e-9, 1e-3, per_decade);
            let r: Vec<f64> = times.iter().map(|t| a + b * (-c * t).exp()).collect();
            let k = effective_coefficient(&times, &r);
            errors.push((1..times.len()).map(|i| (k[i] / exact(times[0], times[i]) - 1.0).abs()).fold(0.0, f64::max));
        }
        assert!(errors.windows(2).all(|e| e[1] < 0.4 * e[0]), "second order: {errors:?}");
        assert!(errors[3] < 2e-3, "{errors:?}");
    }

    #[test]
    fn tau_vib_is_not_given_once_the_mean_energy_is_at_e_f_within_the_tolerance() {
        // tau_vib = -(<E> - E_f)/(d<E>/dt) is a ratio of two rounding-level numbers once <E> has reached E_f: NaN when
        // |<E> - E_f| <= 10 rtol E_f.
        let mut network = two_well_network();
        network.wells.truncate(1);
        network.wells[0].channels.clear();
        let hot = Conditions { temperature_kelvin: 1000.0, pressure_torr: 760.0 };
        let lambda_2 = EigenSystem::new(&final_operator(&network, &hot), EigenSolver::FullDecomposition).unwrap().eigenvalues[1];
        let prep = pulse(thermal(&network, 0, 300.0).unwrap(), vec![BathSegment { start_s: 0.0, conditions: hot }]);
        let r = integrate_preparation(&network, MODEL, SteadyState::Final, &prep, &settings(vec![0.1 / lambda_2, 60.0 / lambda_2])).unwrap();
        assert!(r.points[0].vibrational_relaxation_time_s[0] > 0.0);
        assert!(r.points[1].vibrational_relaxation_time_s[0].is_nan(), "{:e}", r.points[1].vibrational_relaxation_time_s[0]);
    }

    #[test]
    fn a_closed_fragment_system_relaxes_to_the_equilibrium_of_its_weights() {
        // C <-> B + A with the partner in excess (Green, Robertson 2014): from a hot start in C, the populations relax to
        // [C]/[B] = K_c [A], the ratio of the summed weights of the operator. [A] = 1e25 cm-3 puts both wells in the same
        // range (at 1e17 cm-3 this shallow C holds 8e-9 of the population at 1000 K).
        use crate::masterequation::chemical_activation_operator::tests::{fragment_network, prior_kernel};
        let mut network = fragment_network(prior_kernel(), 1e25);
        for well in &mut network.wells {
            well.channels.retain(|ch| !matches!(ch.destination, ChannelDestination::Products { .. }));
        }
        let bath = Conditions { temperature_kelvin: 1000.0, pressure_torr: 760.0 };
        let op = final_operator(&network, &bath);
        let max = op.log_boltzmann_weight.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let (mut c, mut b) = (0.0, 0.0);
        for (s, &(w, _)) in op.states.iter().enumerate() {
            let f = (op.log_boltzmann_weight[s] - max).exp();
            if w == 0 { c += f } else { b += f }
        }
        let hot = gaussian(&network, 0, 3600.0, 100.0, Representation::Density, EnergyReference::AboveWellGround).unwrap();
        let prep = pulse(hot, vec![BathSegment { start_s: 0.0, conditions: bath }]);
        let r = integrate_preparation(&network, MODEL, SteadyState::Final, &prep, &settings(vec![1e-9, 1e-2])).unwrap();
        let last = r.points.last().unwrap();
        assert!((last.well_populations.iter().sum::<f64>() - 1.0).abs() < 1e-6, "{:?}", last.well_populations);
        assert!(last.well_populations.iter().all(|c| *c > 0.05), "both populated: {:?}", last.well_populations);
        let ratio = last.well_populations[0] / last.well_populations[1];
        assert!((ratio / (c / b) - 1.0).abs() < 1e-6, "{ratio:e} vs {:e}", c / b);
        // Early on the hot C has formed fragments, which are not yet in equilibrium.
        assert!((r.points[0].well_populations[0] / r.points[0].well_populations[1] / (c / b) - 1.0).abs() > 1e-3);
    }
}
