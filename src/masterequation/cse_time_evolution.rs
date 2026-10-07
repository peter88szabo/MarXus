//! Validity of the CSE description in time for a prepared experiment (reports/nonthermal_sources_design.md, Section 12).
//!
//! Rate coefficients describe the kinetics "for all times significantly beyond that characterizing the slowest decaying
//! energy relaxation process" (Miller et al., J. Phys. Chem. A 120, 306 (2016)); how accurate the phenomenological
//! description is after a disturbance is the question of Barker, Frenklach, Golden, J. Phys. Chem. A 120, 313 (2016).
//! The test here propagates the CSE description of the same preparation and compares it with the direct time
//! integration of the master equation:
//!   dX/dt = K X + sum_a R_a(t) x_a,   dY_nu/dt = sum_i k_(i->nu) X_i + sum_a R_a(t) p_a,nu,
//! with jumps N_a x_a and N_a p_a at the impulses, K from the phenomenological rate coefficients (G13 eqs. 26, 29), x_a
//! the species populations after relaxation and p_a the prompt yields of the source projection (G13 eqs. 13, 14, 24,
//! 37-41; `SourceProjection`). The CSE description assigns the prompt part at the injection; the direct integration
//! forms it during the relaxation, so the two differ at early times and agree once the relaxation is over.
//!
//! G13: Y. Georgievskii, J. A. Miller, M. P. Burke, S. J. Klippenstein, J. Phys. Chem. A 117, 12146 (2013).

use super::chemical_activation_eigen::EigenSolver;
use super::chemical_activation_network::{
    ChannelDestination, ChemicalActivationNetwork, ChemicalActivationOptions, CollisionModel, SteadyState,
};
use super::chemical_activation_operator::assemble_operator;
use super::chemically_significant_eigenvalues::{phenomenological_rate_coefficients_with_sources, CseMerging};
use super::direct_time_integration::TimeIntegrationSettings;
use super::prepared_time_integration::{project, Preparation, TransientResult};
use super::source_profiles::TimeProfile;
use crate::numeric::dense_inverse::invert_dense;
use crate::numeric::integrators::rosenbrock::{integrate, RosenbrockOptions, StiffSystem};

/// Default bound on the relative deviation |CSE - direct| / max(|CSE|, |direct|) for the agreement time.
pub const DEFAULT_CSE_AGREEMENT_TOLERANCE: f64 = 1e-2;

/// The CSE description of a preparation against the direct time integration, at the output times of the latter.
#[derive(Debug, Clone)]
pub struct CseTimeComparison {
    pub times_s: Vec<f64>,
    /// The CSE species compared (wells, or merged wells joined with "+"); wells of a bimolecular group are left out.
    pub species: Vec<String>,
    /// The bimolecular channels and escape sinks of the CSE description.
    pub bimolecular: Vec<String>,
    /// Per output time: the species populations and cumulative bimolecular yields of the CSE description.
    pub cse_species: Vec<Vec<f64>>,
    pub cse_yields: Vec<Vec<f64>>,
    /// The same from the direct integration (wells summed per species, exits per bimolecular channel).
    pub direct_species: Vec<Vec<f64>>,
    pub direct_yields: Vec<Vec<f64>>,
    /// max over species and channels of |CSE - direct| / N_in(t), N_in the amount injected up to t.
    pub deviation: Vec<f64>,
    /// max over species and channels of |CSE - direct| / max(|CSE|, |direct|), over the quantities above 1e-10 N_in: it
    /// shows a delayed onset of products (incubation), whose absolute size is small.
    pub relative_deviation: Vec<f64>,
    /// The species or channel with the largest relative deviation.
    pub largest: Vec<String>,
    pub tolerance: f64,
    /// The first output time from which on every relative deviation is <= `tolerance`; None if there is none.
    pub agreement_time_s: Option<f64>,
    /// sum of the direct yields / N_in at the agreement time (NaN without one).
    pub conversion_at_agreement: f64,
    /// The lowest eigenvalue that is not chemical (s-1): the agreement time in units of 1/lambda is
    /// agreement_time_s * relaxation_eigenvalue_s_inv.
    pub relaxation_eigenvalue_s_inv: f64,
    /// The chemical eigenvalues (s-1).
    pub chemical_eigenvalues_s_inv: Vec<f64>,
    pub notes: Vec<String>,
}

/// dy/dt = A y + sum_a R_a(t) b_a for y = (species populations, bimolecular yields).
struct SpeciesSystem {
    a: Vec<Vec<f64>>,
    sources: Vec<(Vec<f64>, TimeProfile)>,
    inverse: Vec<Vec<f64>>,
}

impl StiffSystem for SpeciesSystem {
    fn dimension(&self) -> usize {
        self.a.len()
    }

    fn rhs(&self, t: f64, y: &[f64], dydt: &mut [f64]) {
        for (d, row) in dydt.iter_mut().zip(&self.a) {
            *d = row.iter().zip(y).map(|(a, y)| a * y).sum();
        }
        for (b, profile) in &self.sources {
            let r = profile.rate(t);
            if r != 0.0 {
                dydt.iter_mut().zip(b).for_each(|(d, b)| *d += r * b);
            }
        }
    }

    fn prepare(&mut self, _t: f64, _y: &[f64], shift: f64) -> Result<(), String> {
        let g: Vec<Vec<f64>> = self
            .a
            .iter()
            .enumerate()
            .map(|(i, row)| row.iter().enumerate().map(|(j, a)| if i == j { shift - a } else { -a }).collect())
            .collect();
        self.inverse = invert_dense(&g)?;
        Ok(())
    }

    fn solve(&self, b: &mut [f64]) -> Result<(), String> {
        let x: Vec<f64> = self.inverse.iter().map(|row| row.iter().zip(b.iter()).map(|(g, b)| g * b).sum()).collect();
        b.copy_from_slice(&x);
        Ok(())
    }
}

/// Propagates the CSE description of `preparation` (one bath segment, final variant) to the output times of `settings`
/// and compares it with `direct`, the result of `integrate_preparation` for the same preparation, settings and
/// collision model with `SteadyState::Final`. `solver`: a full decomposition (all eigenpairs).
#[allow(clippy::too_many_arguments)]
pub fn compare_cse_with_time_integration(
    network: &ChemicalActivationNetwork,
    collision_model: CollisionModel,
    preparation: &Preparation,
    settings: &TimeIntegrationSettings,
    solver: EigenSolver,
    merging: &CseMerging,
    direct: &TransientResult,
    tolerance: f64,
) -> Result<CseTimeComparison, String> {
    if preparation.bath.len() != 1 {
        return Err("CSE comparison in time: only for a constant bath (one bath segment); the rate coefficients change with T and p.".into());
    }
    if direct.exits.iter().any(|e| e.starts_with("stab(")) {
        return Err("CSE comparison in time: the direct integration must use the final variant (no absorbing barrier).".into());
    }
    let times = &settings.times_s;
    if direct.points.len() != times.len() || direct.points.iter().zip(times).any(|(p, t)| p.time_s != *t) {
        return Err("CSE comparison in time: the direct result is not at the output times of the settings.".into());
    }
    let conditions = &preparation.bath[0].conditions;
    let op = assemble_operator(network, conditions, &ChemicalActivationOptions { collision_model, steady_state: SteadyState::Final })?;

    // Sources on the states: the initial population (if any) first, then the channels.
    let mut names = Vec::new();
    let mut on_states = Vec::new();
    if let Some(initial) = &preparation.initial {
        names.push("initial".to_string());
        on_states.push(project(network, &op, &initial.distribution, 1.0, "initial", 0)?.0);
    }
    for c in &preparation.channels {
        names.push(c.name.clone());
        on_states.push(project(network, &op, &c.distribution, 1.0, &c.name, 0)?.0);
    }
    let sources: Vec<(&str, Vec<f64>)> = names.iter().map(|n| n.as_str()).zip(on_states).collect();
    let rates = phenomenological_rate_coefficients_with_sources(network, &op, None, &[], &sources, solver, merging)?;

    // y = (X_species, Y_bimolecular); dX_j/dt = sum_(i != j) k_(i->j) X_i - k_j X_j, dY_nu/dt = sum_i k_(i->nu) X_i.
    let (n, m) = (rates.wells.len(), rates.bimolecular.len());
    let mut a = vec![vec![0.0; n + m]; n + m];
    for i in 0..n {
        for j in 0..n {
            a[j][i] = if i == j { -rates.well_to_well_s_inv[i][i] } else { rates.well_to_well_s_inv[i][j] };
        }
        for nu in 0..m {
            a[n + nu][i] = rates.well_to_bimolecular_s_inv[i][nu];
        }
    }
    let b: Vec<Vec<f64>> = rates
        .source_projections
        .iter()
        .map(|p| p.species_populations.iter().chain(&p.prompt_bimolecular).copied().collect())
        .collect();
    let offset = usize::from(preparation.initial.is_some());
    let channel_b = |k: usize| &b[offset + k];
    let mut y = vec![0.0; n + m];
    if let Some(initial) = &preparation.initial {
        y.iter_mut().zip(&b[0]).for_each(|(y, b)| *y += initial.amount * b);
    }
    let apply_impulses = |t: f64, y: &mut Vec<f64>| {
        for (k, c) in preparation.channels.iter().enumerate() {
            for (tp, amount) in c.profile.impulses() {
                if tp == t {
                    y.iter_mut().zip(channel_b(k)).for_each(|(y, b)| *y += amount * b);
                }
            }
        }
    };
    apply_impulses(0.0, &mut y);
    let mut system = SpeciesSystem {
        a,
        sources: preparation.channels.iter().enumerate().map(|(k, c)| (channel_b(k).clone(), c.profile.clone())).collect(),
        inverse: Vec::new(),
    };
    let t_end = times[times.len() - 1];
    let mut events: Vec<f64> = times.clone();
    for c in &preparation.channels {
        events.extend(c.profile.impulses().iter().map(|&(t, _)| t));
        events.extend(c.profile.breakpoints());
    }
    events.retain(|&t| t > 0.0 && t <= t_end);
    events.sort_by(|x, y| x.total_cmp(y));
    events.dedup();
    let mut cse = Vec::new();
    let mut t = 0.0;
    for &t_event in &events {
        let mut options = RosenbrockOptions::new(settings.method, settings.relative_tolerance, settings.absolute_tolerance);
        options.autonomous = preparation.channels.iter().all(|c| c.profile.is_constant_on(t, t_event));
        integrate(&mut system, &mut y, t, t_event, &options).map_err(|e| format!("CSE comparison in time: {e}"))?;
        t = t_event;
        apply_impulses(t, &mut y);
        if times.contains(&t) {
            cse.push(y.clone());
        }
    }

    // The direct populations per species and yields per bimolecular channel.
    let mut notes: Vec<String> = rates.warnings.iter().map(|w| format!("CSE: {w}")).collect();
    if !rates.bimolecular_group.is_empty() {
        notes.push(format!(
            "wells {} are in equilibrium with bimolecular species (no CSE species of their own) and are not compared",
            rates.bimolecular_group.join(", ")
        ));
    }
    let well_index = |name: &str| network.wells.iter().position(|w| w.name == name);
    let groups: Vec<Vec<usize>> = rates.well_groups.iter().map(|g| g.iter().filter_map(|w| well_index(w)).collect()).collect();
    if groups.iter().any(|g| g.len() > 1) {
        notes.push("merged species: the direct populations are summed over their wells".into());
    }
    let mut exit_to_channel = Vec::new();
    for exit in &direct.exits {
        let target = if exit.starts_with("escape(") {
            exit.clone()
        } else {
            network
                .wells
                .iter()
                .flat_map(|w| w.channels.iter().map(move |c| (w, c)))
                .find_map(|(w, c)| match &c.destination {
                    ChannelDestination::Products { name } if format!("{}->{name}", w.name) == *exit => Some(name.clone()),
                    _ => None,
                })
                .ok_or_else(|| format!("CSE comparison in time: exit '{exit}' has no product channel."))?
        };
        let nu = rates
            .bimolecular
            .iter()
            .position(|b| *b == target)
            .ok_or_else(|| format!("CSE comparison in time: '{target}' is not a bimolecular channel of the CSE description."))?;
        exit_to_channel.push(nu);
    }
    let mut comparison = CseTimeComparison {
        times_s: times.clone(),
        species: rates.wells.clone(),
        bimolecular: rates.bimolecular.clone(),
        cse_species: Vec::new(),
        cse_yields: Vec::new(),
        direct_species: Vec::new(),
        direct_yields: Vec::new(),
        deviation: Vec::new(),
        relative_deviation: Vec::new(),
        largest: Vec::new(),
        tolerance,
        agreement_time_s: None,
        conversion_at_agreement: f64::NAN,
        relaxation_eigenvalue_s_inv: rates.relaxation_eigenvalue_s_inv,
        chemical_eigenvalues_s_inv: rates.chemical_eigenvalues_s_inv.clone(),
        notes,
    };
    let mut conversions = Vec::new();
    for (p, y) in direct.points.iter().zip(&cse) {
        let species: Vec<f64> = groups.iter().map(|g| g.iter().map(|&w| p.well_populations[w]).sum()).collect();
        let mut yields = vec![0.0; m];
        for (x, &nu) in exit_to_channel.iter().enumerate() {
            yields[nu] += p.exit_yields[x];
        }
        let injected = preparation.initial.as_ref().map_or(0.0, |i| i.amount) + p.injected.iter().sum::<f64>();
        let scale = if injected > 0.0 { injected } else { f64::NAN };
        let (mut worst, mut worst_relative, mut name) = (0.0_f64, 0.0, String::new());
        let pairs = y[..n].iter().zip(&species).zip(&comparison.species).chain(y[n..].iter().zip(&yields).zip(&comparison.bimolecular));
        for ((c, d), label) in pairs {
            worst = worst.max((c - d).abs() / scale);
            let size = c.abs().max(d.abs());
            if size > 1e-10 * scale {
                let relative = (c - d).abs() / size;
                if relative > worst_relative || name.is_empty() {
                    (worst_relative, name) = (relative, label.clone());
                }
            }
        }
        conversions.push(yields.iter().sum::<f64>() / scale);
        comparison.cse_species.push(y[..n].to_vec());
        comparison.cse_yields.push(y[n..].to_vec());
        comparison.direct_species.push(species);
        comparison.direct_yields.push(yields);
        comparison.deviation.push(worst);
        comparison.relative_deviation.push(worst_relative);
        comparison.largest.push(name);
    }
    let first_good = (0..times.len()).rev().take_while(|&k| comparison.relative_deviation[k] <= tolerance).last();
    if let Some(k) = first_good {
        comparison.agreement_time_s = Some(times[k]);
        comparison.conversion_at_agreement = conversions[k];
    }
    Ok(comparison)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::masterequation::chemical_activation_network::Conditions;
    use crate::masterequation::chemical_activation_operator::tests::{conditions, two_well_network};
    use crate::masterequation::direct_time_integration::log_spaced_times;
    use crate::masterequation::prepared_distributions::{gaussian, EnergyReference, Representation};
    use crate::masterequation::prepared_time_integration::{integrate_preparation, BathSegment, InitialPopulation, SourceChannel};

    const MODEL: CollisionModel = CollisionModel::ExponentialDown { cutoff_in_mean_down: 10.0 };

    fn settings(times: Vec<f64>) -> TimeIntegrationSettings {
        TimeIntegrationSettings { relative_tolerance: 1e-9, absolute_tolerance: 1e-18, times_s: times, ..TimeIntegrationSettings::default() }
    }

    fn compare(network: &ChemicalActivationNetwork, prep: &Preparation, times: Vec<f64>) -> CseTimeComparison {
        let s = settings(times);
        let direct = integrate_preparation(network, MODEL, SteadyState::Final, prep, &s).unwrap();
        compare_cse_with_time_integration(network, MODEL, prep, &s, EigenSolver::FullDecomposition, &CseMerging::default(), &direct, DEFAULT_CSE_AGREEMENT_TOLERANCE).unwrap()
    }

    #[test]
    fn a_hot_pulse_differs_during_the_relaxation_and_agrees_after_it() {
        // One reactive well, a pulse above its threshold: the CSE description puts the prompt products at t = 0, the
        // master equation forms them during the relaxation; afterwards both follow the slowest eigenvalue.
        let mut network = two_well_network();
        network.wells.truncate(1);
        network.wells[0].channels.retain(|c| matches!(c.destination, ChannelDestination::Products { .. }));
        let hot = gaussian(&network, 0, 3300.0, 150.0, Representation::Density, EnergyReference::AboveWellGround).unwrap();
        let prep = Preparation {
            initial: Some(InitialPopulation { amount: 1.0, distribution: hot }),
            channels: Vec::new(),
            bath: vec![BathSegment { start_s: 0.0, conditions: conditions() }],
        };
        let c = compare(&network, &prep, log_spaced_times(1e-12, 1e-3, 2));
        // At 1e-12 s the master equation has formed almost nothing yet: the absolute deviation is the prompt yield, the
        // relative one about 1.
        let prompt: f64 = c.cse_yields[0].iter().sum();
        assert!(prompt > 1e-2 && (c.deviation[0] / prompt - 1.0).abs() < 1e-2, "{:?}, prompt {prompt:e}", c.deviation);
        assert!(c.relative_deviation[0] > 0.9, "{:?}", c.relative_deviation);
        assert!(*c.deviation.last().unwrap() < 1e-7, "{:?}", c.deviation);
        assert!(*c.relative_deviation.last().unwrap() < 1e-6, "{:?}", c.relative_deviation);
        let t_star = c.agreement_time_s.expect("agreement after the relaxation");
        // The prompt part decays with the relaxation: agreement within 1e-2 comes after a time of the order of
        // 1/lambda_relax (here 0.9/lambda_relax: 0.024 exp(-0.9) = 0.0098).
        let x = t_star * c.relaxation_eigenvalue_s_inv;
        assert!(x > 0.1 && x < 100.0, "t* = {t_star:e}, lambda = {:e}", c.relaxation_eigenvalue_s_inv);
        assert!(c.conversion_at_agreement > 0.0 && c.conversion_at_agreement < 1.0);
        // The CSE description conserves the injected amount.
        for (s, y) in c.cse_species.iter().zip(&c.cse_yields) {
            assert!((s.iter().sum::<f64>() + y.iter().sum::<f64>() - 1.0).abs() < 1e-8);
        }
    }

    #[test]
    fn a_cold_start_has_a_negative_prompt_yield_whose_logarithm_is_the_incubation_time() {
        // Shock heating (Barker, King 1995): a 300 K population in a 1000 K bath. The CSE description decays from t = 0;
        // the master equation first activates the population. For one well the species population after relaxation is
        // x = a_1 = exp(lambda_1 t_inc) (the back-extrapolation of the first-order decay), so ln(x)/lambda_1 is the
        // incubation time of the direct integration, and the prompt yield 1 - x is negative.
        let mut network = two_well_network();
        network.wells.truncate(1);
        network.wells[0].channels.retain(|c| matches!(c.destination, ChannelDestination::Products { .. }));
        let cold = crate::masterequation::prepared_distributions::thermal(&network, 0, 300.0).unwrap();
        let prep = Preparation {
            initial: Some(InitialPopulation { amount: 1.0, distribution: cold }),
            channels: Vec::new(),
            bath: vec![BathSegment { start_s: 0.0, conditions: Conditions { temperature_kelvin: 1000.0, pressure_torr: 760.0 } }],
        };
        let times = log_spaced_times(1e-12, 1e-7, 8);
        let s = settings(times);
        let direct = integrate_preparation(&network, MODEL, SteadyState::Final, &prep, &s).unwrap();
        let c = compare_cse_with_time_integration(&network, MODEL, &prep, &s, EigenSolver::FullDecomposition, &CseMerging::default(), &direct, 1e-2).unwrap();
        let x = c.cse_species[0][0];
        let prompt: f64 = c.cse_yields[0].iter().sum();
        assert!(x > 1.0 && prompt < 0.0, "x = {x}, prompt = {prompt:e}");
        // t_inc is read where the relaxation is over and the population still far above the absolute tolerance:
        // lambda_1 t = 2.
        let lambda_1 = c.chemical_eigenvalues_s_inv[0];
        let at = direct.points.iter().min_by(|a, b| (a.time_s * lambda_1 - 2.0).abs().total_cmp(&(b.time_s * lambda_1 - 2.0).abs())).unwrap();
        assert!(at.time_s * c.relaxation_eigenvalue_s_inv > 30.0, "{:e} {:e}", at.time_s, c.relaxation_eigenvalue_s_inv);
        let t_inc = at.incubation_time_s;
        // x is the CSE population at the first output time t_0; at 0+ it is x exp(lambda_1 t_0).
        let from_projection = x.ln() / lambda_1 + c.times_s[0];
        assert!((from_projection / t_inc - 1.0).abs() < 1e-5, "{from_projection:e} vs {t_inc:e}");
        // The products are delayed: relative deviations of order 1 at first, agreement only after t_inc.
        assert!(c.relative_deviation[0] > 0.5, "{:?}", c.relative_deviation);
        assert!(c.agreement_time_s.expect("agreement") > t_inc, "{:?} {t_inc:e}", c.agreement_time_s);
        assert!(c.notes.iter().all(|n| !n.contains("source projection")), "{:?}", c.notes);
    }

    #[test]
    fn impulses_and_feeds_of_two_wells_agree_after_the_relaxation() {
        let network = two_well_network();
        let hot = gaussian(&network, 0, 3300.0, 150.0, Representation::Density, EnergyReference::AboveWellGround).unwrap();
        let warm = gaussian(&network, 1, 2500.0, 200.0, Representation::Density, EnergyReference::AboveWellGround).unwrap();
        let prep = Preparation {
            initial: None,
            channels: vec![
                SourceChannel { name: "flash".into(), distribution: hot, profile: TimeProfile::Impulse { time_s: 1e-8, amount: 0.5 } },
                SourceChannel { name: "feed".into(), distribution: warm, profile: TimeProfile::Feed { start_s: 0.0, end_s: None, rate: 1e3 } },
            ],
            bath: vec![BathSegment { start_s: 0.0, conditions: Conditions { temperature_kelvin: 300.0, pressure_torr: 760.0 } }],
        };
        let c = compare(&network, &prep, log_spaced_times(1e-10, 1e-3, 2));
        assert!(!c.species.is_empty() && c.bimolecular.contains(&"P".to_string()), "{:?} {:?}", c.species, c.bimolecular);
        // With a continuous feed the relative deviation levels off: the master equation holds molecules in transit in the
        // relaxational modes (about R c/Lambda per mode), which the CSE description books at once as prompt yield or
        // species population. Here 2.5e-6 of the population of A.
        let last = *c.relative_deviation.last().unwrap();
        assert!(last < 1e-4, "{:?} {:?}", c.relative_deviation, c.largest);
        assert!(c.agreement_time_s.is_some());
        // The projections are resolved: no warning of the eigen route.
        assert!(c.notes.iter().all(|n| !n.contains("source projection")), "{:?}", c.notes);
    }

    #[test]
    fn a_bath_history_and_the_intermediate_variant_are_refused() {
        let network = two_well_network();
        let hot = gaussian(&network, 0, 3300.0, 150.0, Representation::Density, EnergyReference::AboveWellGround).unwrap();
        let mut prep = Preparation {
            initial: Some(InitialPopulation { amount: 1.0, distribution: hot }),
            channels: Vec::new(),
            bath: vec![BathSegment { start_s: 0.0, conditions: conditions() }, BathSegment { start_s: 1e-6, conditions: conditions() }],
        };
        let s = settings(vec![1e-7]);
        let direct = integrate_preparation(&network, MODEL, SteadyState::Final, &prep, &s).unwrap();
        let refused = compare_cse_with_time_integration(&network, MODEL, &prep, &s, EigenSolver::FullDecomposition, &CseMerging::default(), &direct, 1e-2);
        assert!(refused.unwrap_err().contains("constant bath"));
        prep.bath.truncate(1);
        let intermediate = integrate_preparation(&network, MODEL, SteadyState::Intermediate { barrier: Default::default() }, &prep, &s).unwrap();
        let refused = compare_cse_with_time_integration(&network, MODEL, &prep, &s, EigenSolver::FullDecomposition, &CseMerging::default(), &intermediate, 1e-2);
        assert!(refused.unwrap_err().contains("final variant"));
    }
}
