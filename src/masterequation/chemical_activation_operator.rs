//! Assembly of the chemical-activation master-equation operator J of a multiwell network.
//!
//! For every well the energy-grained master equation of Olzmann and co-workers reads
//!   dN/dt = R F - J N,   J = omega (I - P) + K + k_c[D] I              (PO14 eq. 2; O02 eqs. 5-6;
//!                                                                      GO10 eq. 7)
//! and wells are coupled by isomerization: a molecule of well w in grain j reacting through a channel
//! to well w' appears in the grain of w' at the same absolute energy. Column j of J therefore holds
//!   J(j|j)  =  omega_w sum_{t != j} P_w(t|j) + sum_r k_r(E_j) + k_c[D]_w       (all losses of grain j)
//!   J(t|j)  = -omega_w P_w(t|j)                                                 (collisions, same well)
//!   J(t'|j) = -k_{w->w'}(E_j)                                                   (isomerization, t' in w')
//! so that the column sum equals the loss out of the network: product formation, the bimolecular
//! sink and, in the intermediate steady state, the flux into the absorbing region.
//!
//! Intermediate steady state: grains below an absorbing barrier are removed from the state space; the
//! collisional (and isomerization) flux into them is the stabilization (GO10 p. 12295; O02 text after
//! eq. 13: "a lower absorbing barrier in the collisional deactivation cascade at energies below the
//! lowest reaction threshold"). By default the barrier lies 10 k_BT below the lowest reaction threshold
//! of the well: "The absorbing boundary is usually placed about 10 k_BT below the reaction threshold"
//! (PR03 Sec. 2.4); "for example, 10 kT below E0" (CD07 p. 125). The collision kernel is always built
//! on the complete grid of the well, so the transition probabilities out of a retained grain remain
//! normalized and the loss into the absorbed grains is counted exactly.
//!
//! Final steady state: no barrier; with completeness of the transition probabilities the final
//! steady-state distribution follows (GO10 p. 12295).
//!
//! Detailed balance: with Boltzmann weights f = rho_w(E) exp(-E/kT) on the ABSOLUTE energy scale, the
//! collision terms obey P_w(t|j) f_j = P_w(j|t) f_t (both kernels, `collision_kernels.rs`) and the
//! isomerization terms obey k_{w->w'}(E) rho_w(E) = k_{w'->w}(E) rho_{w'}(E) = W‡(E - E0)/h (RRKM;
//! PO14 eq. 9), so J(t|j) f_j = J(j|t) f_t and J is similar to a symmetric matrix
//! (`chemical_activation_steady_state.rs`). `isomerization_detailed_balance` reports how well the
//! supplied isomerization rates fulfil this.
//!
//! State ordering: states are sorted by absolute energy and, at equal energy, by well. Collisions
//! couple grains within the kernel band of one well and isomerization couples grains at the same
//! absolute energy, so every non-zero element of J lies within a band of about
//! (number of wells) x (kernel band) around the diagonal.
//!
//! References: see `chemical_activation_network.rs`.

use crate::constants::KB_CM;

use super::chemical_activation_network::{
    AbsorbingBarrier, ChannelDestination, ChemicalActivationNetwork, ChemicalActivationOptions,
    CollisionModel, Conditions, SteadyState, Well,
};
use super::collision_kernels::{exponential_down_kernel, stepladder_kernel, CollisionKernel};
use super::collisional_relaxation::lennard_jones_collision_frequency_s_inv;

/// Collision data of one well at the given conditions.
#[derive(Debug, Clone)]
pub struct WellCollisionData {
    /// Collision frequency omega = Z_LJ [M], s-1 (T77 eqs. 3.1-3.3).
    pub collision_frequency_s_inv: f64,
    /// <dE_down>(T), cm-1.
    pub mean_down_cm1: f64,
    /// Lowest retained grain: the absorbing barrier (intermediate steady state), 0 in the final steady state.
    pub absorbing_barrier_grain: usize,
    /// Exponential down: the grains 0 .. reservoir_grains form the reservoir of the well (the normalization of
    /// Robertson (2019) eq. 4.16 fails at the highest of them; `collision_kernels.rs`), one thermalized
    /// state in J (MESMER manual, Sec. 14.2.1). 0 without a reservoir.
    pub reservoir_grains: usize,
    /// Stepladder: step size in grains.
    pub step_grains: usize,
}

/// The operator J of the steady-state equation J N = F on the retained grains of all wells.
#[derive(Debug, Clone)]
pub struct ChemicalActivationOperator {
    /// (well, grain) of every state, in the order of the rows and columns of J.
    pub states: Vec<(usize, usize)>,
    /// index_of[well][grain] = position of the state, None for absorbed grains.
    pub index_of: Vec<Vec<Option<usize>>>,
    /// J in row form: rows[r] = [(c, J_rc)], columns sorted, (J N)_r = sum_c J_rc N_c. Units s-1.
    pub rows: Vec<Vec<(usize, f64)>>,
    /// ln f of every state, f = rho_w(E) exp(-E/kT) with E on the absolute energy scale.
    pub log_boltzmann_weight: Vec<f64>,
    /// For every state: (well, rate in s-1) of the flux into the absorbed grains of that well.
    pub stabilization: Vec<Vec<(usize, f64)>>,
    pub wells: Vec<WellCollisionData>,
    /// Per well: its reservoir state, if the well has one in J.
    pub reservoirs: Vec<Option<Reservoir>>,
    /// k_B T in cm-1.
    pub kt_cm1: f64,
}

/// The reservoir of a well: its lowest grains, in Boltzmann equilibrium among themselves, represented by
/// one state of J (MESMER manual, Sec. 14.2.1). The state is listed in `states` as (well, 0).
#[derive(Debug, Clone)]
pub struct Reservoir {
    /// Index of the reservoir state.
    pub state: usize,
    /// f_i / Q_res of the grains 0 .. weights.len() (they sum to 1).
    pub weights: Vec<f64>,
    /// ln Q_res = ln sum_i rho_i exp(-E_i/kT), E on the absolute energy scale (the weight of the state).
    pub log_weight: f64,
}

impl ChemicalActivationOperator {
    pub fn dimension(&self) -> usize {
        self.states.len()
    }

    /// The reservoir whose state is `s`, if any.
    fn reservoir_of_state(&self, s: usize) -> Option<&Reservoir> {
        let (w, _) = self.states[s];
        self.reservoirs[w].as_ref().filter(|r| r.state == s)
    }

    /// The rate of a process with grain rates `k` (on the grid of the state's well) out of state `s`: k at
    /// the grain of the state, or for a reservoir state the average sum_i k_i f_i/Q_res over its grains.
    pub fn state_rate(&self, s: usize, k: &[f64]) -> f64 {
        match self.reservoir_of_state(s) {
            Some(r) => r.weights.iter().zip(k).map(|(w, k)| w * k).sum(),
            None => k[self.states[s].1],
        }
    }

    /// Populations on the complete grid of every well (zero in absorbed grains) from a vector over the
    /// states; a reservoir population is spread over its grains with f_i/Q_res.
    pub fn grain_populations(&self, x: &[f64]) -> Vec<Vec<f64>> {
        let mut out: Vec<Vec<f64>> = self.index_of.iter().map(|g| vec![0.0; g.len()]).collect();
        for (s, &(w, i)) in self.states.iter().enumerate() {
            match self.reservoir_of_state(s) {
                Some(r) => r.weights.iter().enumerate().for_each(|(g, f)| out[w][g] = x[s] * f),
                None => out[w][i] = x[s],
            }
        }
        out
    }

    /// The inverse of `grain_populations`: a vector over the states from populations on the grains (the
    /// grains of a reservoir are summed; absorbed grains are left out).
    pub fn state_populations(&self, grains: &[Vec<f64>]) -> Vec<f64> {
        let mut out = vec![0.0; self.dimension()];
        for (w, well) in self.index_of.iter().enumerate() {
            for (i, s) in well.iter().enumerate() {
                if let Some(s) = s {
                    out[*s] += grains[w][i];
                }
            }
        }
        out
    }

    /// (J x) for a vector over the retained states.
    pub fn apply(&self, x: &[f64]) -> Vec<f64> {
        self.rows.iter().map(|row| row.iter().map(|&(c, v)| v * x[c]).sum()).collect()
    }
}

/// Collision kernel of one well at temperature T. It always covers the complete grid of the well, so
/// that the probabilities out of a retained grain stay normalized when the grains below an absorbing
/// barrier are absorbing.
fn well_collision_kernel(
    well: &Well,
    grain_width_cm1: f64,
    temperature_kelvin: f64,
    collision_model: CollisionModel,
) -> Result<CollisionKernel, String> {
    let kt_cm1 = KB_CM * temperature_kelvin;
    let mean_down = well.energy_transfer.mean_down_cm1(temperature_kelvin);
    match collision_model {
        CollisionModel::ExponentialDown { cutoff_in_mean_down } => {
            if !(cutoff_in_mean_down > 0.0) {
                return Err("The exponential-down cutoff must be positive.".into());
            }
            let band = (cutoff_in_mean_down * mean_down / grain_width_cm1).ceil() as usize;
            exponential_down_kernel(&well.density_of_states, grain_width_cm1, mean_down, kt_cm1, band)
        }
        // Step size dE_SL = <dE_down> (GO10, text before eq. 16; PO14: "the step size, dE_SL, which
        // corresponds to the average energy transferred per down collision").
        CollisionModel::Stepladder => stepladder_kernel(&well.density_of_states, grain_width_cm1, mean_down, kt_cm1),
    }
    .map_err(|e| format!("Well '{}': {e}", well.name))
}

/// Low-energy reservoir of one well: its lowest grains, from the grain where the normalization of the
/// exponential-down kernel (Robertson 2019, eq. 4.16) fails down, form one thermalized state (MESMER
/// manual, Sec. 14.2.1; `collision_kernels.rs`).
#[derive(Debug, Clone, PartialEq)]
pub struct LowEnergyReservoir {
    pub well: String,
    /// Number of reservoir grains, counted from the well bottom.
    pub grains: usize,
}

/// The wells with a low-energy reservoir at temperature T (the kernel does not depend on the pressure).
/// Empty when eq. 4.16 holds in every well, and for the stepladder.
pub fn low_energy_reservoirs(
    network: &ChemicalActivationNetwork,
    temperature_kelvin: f64,
    collision_model: CollisionModel,
) -> Result<Vec<LowEnergyReservoir>, String> {
    let mut out = Vec::new();
    for well in &network.wells {
        let kernel = well_collision_kernel(well, network.grain_width_cm1, temperature_kelvin, collision_model)?;
        if kernel.reservoir_grains > 0 {
            out.push(LowEnergyReservoir { well: well.name.clone(), grains: kernel.reservoir_grains });
        }
    }
    Ok(out)
}

/// Assemble J for the network at temperature T and bath-gas pressure p.
pub fn assemble_operator(
    network: &ChemicalActivationNetwork,
    conditions: &Conditions,
    options: &ChemicalActivationOptions,
) -> Result<ChemicalActivationOperator, String> {
    network.validate()?;
    let temperature = conditions.temperature_kelvin;
    if !(temperature > 0.0) || !(conditions.pressure_torr > 0.0) {
        return Err(format!(
            "Temperature ({temperature} K) and pressure ({} Torr) must be positive.",
            conditions.pressure_torr
        ));
    }
    let kt_cm1 = KB_CM * temperature;
    let d_e = network.grain_width_cm1;

    // Per well: collision frequency, <dE_down>(T), absorbing barrier and collision kernel.
    let mut wells = Vec::with_capacity(network.wells.len());
    let mut kernels = Vec::with_capacity(network.wells.len());
    for well in &network.wells {
        let lj = &well.lennard_jones;
        let omega = lennard_jones_collision_frequency_s_inv(
            lj.sigma_angstrom,
            lj.epsilon_kelvin,
            lj.reduced_mass_amu,
            temperature,
            conditions.pressure_torr,
        )
        .map_err(|e| format!("Well '{}': {e}", well.name))?;
        let mean_down = well.energy_transfer.mean_down_cm1(temperature);
        let kernel = well_collision_kernel(well, d_e, temperature, options.collision_model)?;

        let barrier = match &options.steady_state {
            SteadyState::Final => 0,
            SteadyState::Intermediate { barrier } => match barrier {
                AbsorbingBarrier::BelowLowestThreshold { kt_multiple } => {
                    if !(*kt_multiple >= 0.0) {
                        return Err("The absorbing-barrier distance in k_BT must be >= 0.".into());
                    }
                    let threshold = well.lowest_threshold_grain().ok_or_else(|| {
                        format!(
                            "Well '{}' has no open reaction channel, so the absorbing barrier below its \
                             lowest reaction threshold is undefined.",
                            well.name
                        )
                    })?;
                    let distance = (kt_multiple * kt_cm1 / d_e).round() as usize;
                    if distance >= threshold {
                        return Err(format!(
                            "Well '{}': the absorbing barrier {kt_multiple} k_BT below the lowest threshold lies at \
                             or below the bottom of the well (threshold {:.0} cm-1 = {:.1} k_BT above the well \
                             bottom). Nothing can be stabilized and the intermediate steady state is not defined \
                             at {temperature} K; choose a smaller distance (kt_multiple), the final steady state, or the \
                             eigenvalue analysis (k_uni = lambda_1, no absorbing barrier).",
                            well.name,
                            threshold as f64 * d_e,
                            threshold as f64 * d_e / kt_cm1
                        ));
                    }
                    threshold - distance
                }
                AbsorbingBarrier::AtGrains(grains) => {
                    if grains.len() != network.wells.len() {
                        return Err(format!(
                            "{} absorbing-barrier grains given for {} wells.",
                            grains.len(),
                            network.wells.len()
                        ));
                    }
                    let w = wells.len();
                    if grains[w] >= well.grain_count() {
                        return Err(format!(
                            "Well '{}': absorbing barrier at grain {} beyond the grid ({} grains).",
                            well.name,
                            grains[w],
                            well.grain_count()
                        ));
                    }
                    grains[w]
                }
            },
        };
        // An absorbing barrier inside the reservoir would split a thermalized state: refused. A barrier at or
        // above its top absorbs the whole reservoir (no reservoir state).
        if barrier > 0 && barrier < kernel.reservoir_grains {
            return Err(format!(
                "Well '{}': the absorbing barrier at grain {barrier} lies within the {} lowest grains, which form \
                 the low-energy reservoir of the well (the normalization of the collision kernel, Robertson 2019, \
                 eq. 4.16, fails there at {temperature} K; MESMER manual, Sec. 14.2.1). Choose a larger barrier \
                 distance in kT above the bottom (a smaller distance below the threshold), or the final steady state.",
                well.name, kernel.reservoir_grains
            ));
        }

        wells.push(WellCollisionData {
            collision_frequency_s_inv: omega,
            mean_down_cm1: mean_down,
            absorbing_barrier_grain: barrier,
            reservoir_grains: kernel.reservoir_grains,
            step_grains: kernel.step_grains,
        });
        kernels.push(kernel);
    }

    // Retained states, ordered by absolute energy and then by well. A well with a reservoir that is not
    // absorbed has one state (w, 0) for all its reservoir grains.
    let has_reservoir: Vec<bool> =
        wells.iter().map(|d| d.reservoir_grains > 0 && d.absorbing_barrier_grain == 0).collect();
    let mut states: Vec<(usize, usize)> = Vec::new();
    for (w, well) in network.wells.iter().enumerate() {
        let first = if has_reservoir[w] {
            states.push((w, 0));
            wells[w].reservoir_grains
        } else {
            wells[w].absorbing_barrier_grain
        };
        states.extend((first..well.grain_count()).map(|i| (w, i)));
    }
    states.sort_by_key(|&(w, i)| (i as isize + network.wells[w].bottom_offset_grains, w));
    let mut index_of: Vec<Vec<Option<usize>>> =
        network.wells.iter().map(|well| vec![None; well.grain_count()]).collect();
    let mut reservoirs: Vec<Option<Reservoir>> = vec![None; network.wells.len()];
    for (s, &(w, i)) in states.iter().enumerate() {
        if has_reservoir[w] && i == 0 {
            let m = wells[w].reservoir_grains;
            index_of[w][..m].iter_mut().for_each(|x| *x = Some(s));
            // f_i/Q_res, evaluated relative to the largest weight.
            let log_f: Vec<f64> = (0..m)
                .map(|g| network.wells[w].density_of_states[g].ln() - network.absolute_energy_cm1(w, g) / kt_cm1)
                .collect();
            let max = log_f.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let f: Vec<f64> = log_f.iter().map(|l| (l - max).exp()).collect();
            let q: f64 = f.iter().sum();
            reservoirs[w] =
                Some(Reservoir { state: s, weights: f.iter().map(|x| x / q).collect(), log_weight: max + q.ln() });
        } else {
            index_of[w][i] = Some(s);
        }
    }

    // Columns of J (source state c = (w, j)), stored by rows.
    let n = states.len();
    let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    let mut stabilization: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for (c, &(w, j)) in states.iter().enumerate() {
        let well = &network.wells[w];
        let omega = wells[w].collision_frequency_s_inv;
        let reservoir = reservoirs[w].as_ref().filter(|r| r.state == c);
        let mut diagonal = 0.0;

        // Collisions, omega (I - P): the diagonal collects omega sum_{t != j} P(t|j), which equals
        // omega (1 - P(j|j)) of PO14 eq. 2 for a normalized kernel and keeps the column sum exactly
        // equal to the loss out of the network. Transitions into reservoir grains go into the reservoir
        // state (index_of maps all its grains to it).
        match reservoir {
            None => {
                for &(t, p) in &kernels[w].transitions[j] {
                    let rate = omega * p;
                    diagonal += rate;
                    match index_of[w][t] {
                        Some(r) => rows[r].push((c, -rate)),
                        None => stabilization[c].push((w, rate)),
                    }
                }
            }
            // Out of the reservoir: activation into every grain t above it, by detailed balance with the
            // deactivation from t into the reservoir (MESMER manual, Sec. 14.2.1, eqs. 14.15-14.16):
            //   rate(t <- reservoir) = omega sum_{i in reservoir} P(i|t) f_t / Q_res.
            Some(res) => {
                let m = res.weights.len();
                for t in m..well.grain_count() {
                    let down: f64 = kernels[w].transitions[t].iter().filter(|&&(i, _)| i < m).map(|&(_, p)| p).sum();
                    if down == 0.0 {
                        continue;
                    }
                    let log_f_t = well.density_of_states[t].ln() - network.absolute_energy_cm1(w, t) / kt_cm1;
                    let rate = omega * down * (log_f_t - res.log_weight).exp();
                    diagonal += rate;
                    match index_of[w][t] {
                        Some(r) => rows[r].push((c, -rate)),
                        None => stabilization[c].push((w, rate)),
                    }
                }
            }
        }

        // Unimolecular channels, K (PO14 eq. 2): products leave the network, isomerization feeds the
        // grain of the target well at the same absolute energy. Out of a reservoir every grain i reacts
        // with its share f_i/Q_res of the reservoir population (the grains are thermalized).
        let sources: Vec<(usize, f64)> = match reservoir {
            None => vec![(j, 1.0)],
            Some(res) => res.weights.iter().copied().enumerate().collect(),
        };
        for channel in &well.channels {
            for &(g, share) in &sources {
                let k = channel.rate_constant_s_inv[g] * share;
                if k == 0.0 {
                    continue;
                }
                diagonal += k;
                if let ChannelDestination::Well { index: target } = channel.destination {
                    let t = network.aligned_grain(w, g, target);
                    let target_well = &network.wells[target];
                    if t < 0 || t as usize >= target_well.grain_count() {
                        return Err(format!(
                            "Channel '{}' of well '{}' is open at grain {g} ({} cm-1 absolute), outside the grid \
                             of well '{}' (absolute {} to {} cm-1). Extend the energy grids to a common top.",
                            channel.name,
                            well.name,
                            network.absolute_energy_cm1(w, g),
                            target_well.name,
                            target_well.bottom_offset_grains as f64 * d_e,
                            (target_well.bottom_offset_grains + target_well.grain_count() as isize - 1) as f64 * d_e
                        ));
                    }
                    match index_of[target][t as usize] {
                        Some(r) => rows[r].push((c, -k)),
                        None => stabilization[c].push((target, k)),
                    }
                }
            }
        }

        // Bimolecular sink k_c[D] I (PO14 eq. 2).
        diagonal += well.bimolecular_sink_s_inv;
        rows[c].push((c, diagonal));
    }

    // Sort every row by column and merge repeated entries (several channels into the same grain).
    for row in rows.iter_mut() {
        row.sort_by_key(|&(c, _)| c);
        let mut merged: Vec<(usize, f64)> = Vec::with_capacity(row.len());
        for &(c, v) in row.iter() {
            match merged.last_mut() {
                Some(last) if last.0 == c => last.1 += v,
                _ => merged.push((c, v)),
            }
        }
        *row = merged;
    }

    let log_boltzmann_weight = states
        .iter()
        .enumerate()
        .map(|(s, &(w, i))| match reservoirs[w].as_ref().filter(|r| r.state == s) {
            Some(res) => res.log_weight,
            None => network.wells[w].density_of_states[i].ln() - network.absolute_energy_cm1(w, i) / kt_cm1,
        })
        .collect();

    Ok(ChemicalActivationOperator {
        states,
        index_of,
        rows,
        log_boltzmann_weight,
        stabilization,
        wells,
        reservoirs,
        kt_cm1,
    })
}

/// Deviation from microscopic reversibility of the isomerization rates between two wells.
#[derive(Debug, Clone)]
pub struct IsomerizationBalance {
    pub well_a: usize,
    pub well_b: usize,
    /// max over common grains of |rho_a k_ab - rho_b k_ba| / max(rho_a k_ab, rho_b k_ba).
    pub max_relative_deviation: f64,
}

/// Check rho_a(E) sum k_{a->b}(E) = rho_b(E) sum k_{b->a}(E) at every common absolute energy for every
/// pair of wells connected by isomerization (both equal sum W‡/h; PO14 eq. 9).
pub fn isomerization_detailed_balance(network: &ChemicalActivationNetwork) -> Vec<IsomerizationBalance> {
    // Total isomerization rate from well `from` to well `to` in grain i of `from`.
    let total_rate = |from: usize, to: usize, i: usize| -> f64 {
        network.wells[from]
            .channels
            .iter()
            .filter(|ch| ch.destination == ChannelDestination::Well { index: to })
            .map(|ch| ch.rate_constant_s_inv[i])
            .sum()
    };
    let connected = |from: usize, to: usize| {
        network.wells[from]
            .channels
            .iter()
            .any(|ch| ch.destination == ChannelDestination::Well { index: to })
    };

    let mut result = Vec::new();
    for a in 0..network.wells.len() {
        for b in (a + 1)..network.wells.len() {
            if !connected(a, b) && !connected(b, a) {
                continue;
            }
            let mut max_dev: f64 = 0.0;
            for i in 0..network.wells[a].grain_count() {
                let j = network.aligned_grain(a, i, b);
                if j < 0 || j as usize >= network.wells[b].grain_count() {
                    continue;
                }
                let j = j as usize;
                let forward = network.wells[a].density_of_states[i] * total_rate(a, b, i);
                let backward = network.wells[b].density_of_states[j] * total_rate(b, a, j);
                let larger = forward.max(backward);
                if larger > 0.0 {
                    max_dev = max_dev.max((forward - backward).abs() / larger);
                }
            }
            result.push(IsomerizationBalance { well_a: a, well_b: b, max_relative_deviation: max_dev });
        }
    }
    result
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::masterequation::chemical_activation_network::tests::test_well;
    use crate::masterequation::chemical_activation_network::{Channel, Well};

    const D_E: f64 = 10.0;
    const PLANCK_CM1_S: f64 = 3.3356e-11;

    /// Sum of states of an isomerization transition state at absolute grain `ts`.
    fn ts_sum_of_states(absolute_grain: isize, ts: isize) -> f64 {
        if absolute_grain < ts {
            0.0
        } else {
            (1.0 + 0.05 * (absolute_grain - ts) as f64).powi(6)
        }
    }

    /// Two wells: A (bottom at absolute grain 0, 400 grains) and B (bottom 60 grains lower, 460
    /// grains, same top). A <-> B through a transition state at absolute grain 250 with
    /// k = W‡/(h rho) in both directions; A -> products from 300, B -> products from 280.
    pub(crate) fn two_well_network() -> ChemicalActivationNetwork {
        let mut a = test_well("A", 400, 0, 300);
        // At 300 K the lowest grains of A and B form low-energy reservoirs (the exponential-down normalization
        // fails there), below the absorbing barrier 10 k_BT under the lowest threshold (A: grain 41, B: 269).
        a.density_of_states = (0..400).map(|i| (1.0 + 0.02 * i as f64).powi(8)).collect();
        let mut b: Well = test_well("B", 460, -60, 340);
        b.density_of_states = (0..460).map(|i| (1.0 + 0.04 * i as f64).powi(9)).collect();
        b.bimolecular_sink_s_inv = 1.0e5;
        let ts = 250;
        a.channels.push(Channel {
            name: "A->B".into(),
            destination: ChannelDestination::Well { index: 1 },
            threshold_grain: None,
            rate_constant_s_inv: (0..400)
                .map(|i| ts_sum_of_states(i as isize, ts) / (PLANCK_CM1_S * a.density_of_states[i]))
                .collect(),
        });
        b.channels.push(Channel {
            name: "B->A".into(),
            destination: ChannelDestination::Well { index: 0 },
            threshold_grain: None,
            rate_constant_s_inv: (0..460)
                .map(|i| ts_sum_of_states(i as isize - 60, ts) / (PLANCK_CM1_S * b.density_of_states[i]))
                .collect(),
        });
        ChemicalActivationNetwork { grain_width_cm1: D_E, wells: vec![a, b] }
    }

    pub(crate) fn conditions() -> Conditions {
        Conditions { temperature_kelvin: 300.0, pressure_torr: 760.0 }
    }

    pub(crate) fn all_option_combinations() -> Vec<ChemicalActivationOptions> {
        let mut out = Vec::new();
        for collision_model in [CollisionModel::ExponentialDown { cutoff_in_mean_down: 10.0 }, CollisionModel::Stepladder] {
            for steady_state in [SteadyState::Final, SteadyState::Intermediate { barrier: AbsorbingBarrier::default() }] {
                out.push(ChemicalActivationOptions { collision_model, steady_state });
            }
        }
        out
    }

    /// Well "Sparse" (rho = (1 + 0.05 i)^8): at 300 K on 20 cm-1 grains eq. 4.16 fails near its bottom,
    /// so its lowest grains form a low-energy reservoir; well "Smooth" (rho = exp(0.04 i)) has none.
    fn sparse_and_smooth_network() -> ChemicalActivationNetwork {
        let sparse = test_well("Sparse", 200, 0, 150);
        let mut smooth = test_well("Smooth", 200, 0, 150);
        smooth.density_of_states = (0..200).map(|i| (0.04 * i as f64).exp()).collect();
        ChemicalActivationNetwork { grain_width_cm1: 20.0, wells: vec![sparse, smooth] }
    }

    const EXPONENTIAL_DOWN: CollisionModel = CollisionModel::ExponentialDown { cutoff_in_mean_down: 10.0 };

    #[test]
    fn low_energy_reservoirs_list_only_the_wells_where_eq_4_16_fails() {
        let network = sparse_and_smooth_network();
        let reservoirs = low_energy_reservoirs(&network, 300.0, EXPONENTIAL_DOWN).unwrap();
        assert_eq!(reservoirs.len(), 1);
        assert_eq!(reservoirs[0].well, "Sparse");
        assert!(reservoirs[0].grains > 0);
        // The kernel of the operator is the same.
        let options = ChemicalActivationOptions { collision_model: EXPONENTIAL_DOWN, steady_state: SteadyState::Final };
        let op = assemble_operator(&network, &conditions(), &options).unwrap();
        assert_eq!(op.wells[0].reservoir_grains, reservoirs[0].grains);
        assert_eq!(op.wells[1].reservoir_grains, 0);
        assert!(low_energy_reservoirs(&network, 300.0, CollisionModel::Stepladder).unwrap().is_empty());
    }

    #[test]
    fn the_reservoir_is_one_thermalized_state_with_the_boltzmann_weight_of_its_grains() {
        // MESMER manual, Sec. 14.2.1: the reservoir grains are represented by one state in Boltzmann
        // equilibrium; its weight is their partition function, and every process out of reservoir grain i
        // goes with the fraction f_i/Q_res of the reservoir population.
        let network = sparse_and_smooth_network();
        let options = ChemicalActivationOptions { collision_model: EXPONENTIAL_DOWN, steady_state: SteadyState::Final };
        let op = assemble_operator(&network, &conditions(), &options).unwrap();
        let m = op.wells[0].reservoir_grains;
        let r = op.index_of[0][0].unwrap();
        assert!((0..m).all(|i| op.index_of[0][i] == Some(r)));
        assert!((m..200).all(|i| op.index_of[0][i].is_some_and(|s| s != r)));
        assert_eq!(op.dimension(), 400 - m + 1);
        let kt = KB_CM * 300.0;
        let f: Vec<f64> = (0..m)
            .map(|i| network.wells[0].density_of_states[i] * (-network.absolute_energy_cm1(0, i) / kt).exp())
            .collect();
        let q: f64 = f.iter().sum();
        assert!((op.log_boltzmann_weight[r] - q.ln()).abs() < 1e-12);
        // The rate of a channel out of the reservoir is its average over the reservoir grains.
        let k = &network.wells[0].channels[0].rate_constant_s_inv;
        let average: f64 = (0..m).map(|i| k[i] * f[i] / q).sum();
        assert!((op.state_rate(r, k) - average).abs() <= 1e-15 * average.max(1e-300));
        // Populations: the reservoir population is spread over its grains with f_i/Q, and summed back.
        let mut x = vec![0.0; op.dimension()];
        x[r] = 2.0;
        x[op.index_of[0][m].unwrap()] = 1.0;
        let grains = op.grain_populations(&x);
        for i in 0..m {
            assert!((grains[0][i] - 2.0 * f[i] / q).abs() < 1e-15 * 2.0);
        }
        assert_eq!(grains[0][m], 1.0);
        let back = op.state_populations(&grains);
        assert!((back[r] - 2.0).abs() < 1e-14 && back[op.index_of[0][m].unwrap()] == 1.0);
        // The reservoir exchanges population with the grains above it by collisions (activation from
        // detailed balance), and nothing is absorbed in the final steady state.
        assert!(op.rows.iter().any(|row| row.iter().any(|&(c, v)| c == r && v < 0.0)));
        assert!(op.stabilization.iter().all(|s| s.is_empty()));
    }

    #[test]
    fn an_absorbing_barrier_within_the_reservoir_is_an_error() {
        let network = sparse_and_smooth_network();
        let m = low_energy_reservoirs(&network, 300.0, EXPONENTIAL_DOWN).unwrap()[0].grains;
        let options = ChemicalActivationOptions {
            collision_model: EXPONENTIAL_DOWN,
            steady_state: SteadyState::Intermediate { barrier: AbsorbingBarrier::AtGrains(vec![m - 1, 50]) },
        };
        let err = assemble_operator(&network, &conditions(), &options).unwrap_err();
        assert!(err.contains("reservoir"), "{err}");
        // A barrier at or above the top of the reservoir absorbs the whole reservoir: no reservoir state.
        let options = ChemicalActivationOptions {
            collision_model: EXPONENTIAL_DOWN,
            steady_state: SteadyState::Intermediate { barrier: AbsorbingBarrier::AtGrains(vec![m, 50]) },
        };
        let op = assemble_operator(&network, &conditions(), &options).unwrap();
        assert!((0..m).all(|i| op.index_of[0][i].is_none()));
    }

    #[test]
    fn isomerization_into_reservoir_grains_feeds_the_reservoir_with_detailed_balance() {
        // Smooth <-> Sparse at every grain, with microscopic reversibility rho_a k_ab = rho_b k_ba: the flux
        // into reservoir grains of Sparse goes into its reservoir state, the reverse goes out of it with
        // f_i/Q_res, and J stays detailed balanced with the Boltzmann weights of the states.
        let mut network = sparse_and_smooth_network();
        let rho_a = network.wells[1].density_of_states.clone();
        let rho_b = network.wells[0].density_of_states.clone();
        network.wells[1].channels.push(Channel {
            name: "Smooth->Sparse".into(),
            destination: ChannelDestination::Well { index: 0 },
            threshold_grain: None,
            rate_constant_s_inv: (0..200).map(|i| 1.0e3 * rho_b[i] / (rho_a[i] + rho_b[i])).collect(),
        });
        network.wells[0].channels.push(Channel {
            name: "Sparse->Smooth".into(),
            destination: ChannelDestination::Well { index: 1 },
            threshold_grain: None,
            rate_constant_s_inv: (0..200).map(|i| 1.0e3 * rho_a[i] / (rho_a[i] + rho_b[i])).collect(),
        });
        let options = ChemicalActivationOptions { collision_model: EXPONENTIAL_DOWN, steady_state: SteadyState::Final };
        let op = assemble_operator(&network, &conditions(), &options).unwrap();
        let r = op.index_of[0][0].unwrap();
        let from_smooth_0 = op.index_of[1][0].unwrap();
        assert!(element(&op, r, from_smooth_0) < 0.0, "Smooth grain 0 feeds the reservoir of Sparse");
        for (row_index, row) in op.rows.iter().enumerate() {
            for &(c, v) in row.iter().filter(|&&(c, _)| c != row_index) {
                let lhs = v.abs();
                let rhs = element(&op, c, row_index).abs() * (op.log_boltzmann_weight[row_index] - op.log_boltzmann_weight[c]).exp();
                assert!(((lhs - rhs) / lhs.max(rhs)).abs() < 1e-10, "J({row_index}|{c})");
            }
        }
    }

    fn element(op: &ChemicalActivationOperator, r: usize, c: usize) -> f64 {
        op.rows[r].iter().find(|(col, _)| *col == c).map(|(_, v)| *v).unwrap_or(0.0)
    }

    #[test]
    fn column_sums_of_j_equal_the_losses_out_of_the_network() {
        let network = two_well_network();
        for options in all_option_combinations() {
            let op = assemble_operator(&network, &conditions(), &options).unwrap();
            let mut column_sum = vec![0.0; op.dimension()];
            for row in &op.rows {
                for &(c, v) in row {
                    column_sum[c] += v;
                }
            }
            for (c, &(w, i)) in op.states.iter().enumerate() {
                let well = &network.wells[w];
                let products: f64 = well
                    .channels
                    .iter()
                    .filter(|ch| matches!(ch.destination, ChannelDestination::Products { .. }))
                    .map(|ch| op.state_rate(c, &ch.rate_constant_s_inv))
                    .sum();
                let stabilization: f64 = op.stabilization[c].iter().map(|(_, k)| k).sum();
                let loss = products + well.bimolecular_sink_s_inv + stabilization;
                let diagonal = element(&op, c, c);
                assert!(
                    (column_sum[c] - loss).abs() <= 1e-10 * diagonal,
                    "{options:?}: state {c} = ({w},{i}): column sum {:e}, losses {loss:e}",
                    column_sum[c]
                );
            }
        }
    }

    #[test]
    fn j_is_detailed_balanced_with_boltzmann_weights_on_the_absolute_energy_scale() {
        let network = two_well_network();
        for options in all_option_combinations() {
            let op = assemble_operator(&network, &conditions(), &options).unwrap();
            let mut couplings_between_wells = 0;
            for (r, row) in op.rows.iter().enumerate() {
                for &(c, v) in row {
                    if c == r {
                        continue;
                    }
                    if op.states[r].0 != op.states[c].0 {
                        couplings_between_wells += 1;
                    }
                    // J_rc f_c = J_cr f_r
                    let lhs = v.abs();
                    let rhs = element(&op, c, r).abs()
                        * (op.log_boltzmann_weight[r] - op.log_boltzmann_weight[c]).exp();
                    assert!(
                        ((lhs - rhs) / lhs.max(rhs)).abs() < 1e-10,
                        "{options:?}: J({r}|{c}) f_{c} = {lhs:e} vs J({c}|{r}) f_{r} = {rhs:e}"
                    );
                }
            }
            assert!(couplings_between_wells > 0, "the isomerization must couple the wells");
        }
    }

    #[test]
    fn intermediate_steady_state_absorbs_grains_ten_kt_below_the_lowest_threshold() {
        let network = two_well_network();
        let options = ChemicalActivationOptions {
            collision_model: CollisionModel::ExponentialDown { cutoff_in_mean_down: 10.0 },
            steady_state: SteadyState::Intermediate { barrier: AbsorbingBarrier::default() },
        };
        let op = assemble_operator(&network, &conditions(), &options).unwrap();
        let ten_kt_in_grains = (10.0 * KB_CM * 300.0 / D_E).round() as usize;
        // Lowest thresholds: A at grain 250 (isomerization), B at grain 250 + 60 = 310 (isomerization).
        let expected = [250 - ten_kt_in_grains, 310 - ten_kt_in_grains];
        for w in 0..2 {
            assert_eq!(op.wells[w].absorbing_barrier_grain, expected[w]);
            for i in 0..network.wells[w].grain_count() {
                assert_eq!(op.index_of[w][i].is_some(), i >= expected[w], "well {w}, grain {i}");
            }
        }
        assert!(op.stabilization.iter().any(|s| !s.is_empty()));
    }

    #[test]
    fn an_absorbing_barrier_below_the_well_bottom_is_an_error() {
        // At 400 K, 10 kT = 2780 cm-1 exceeds the lowest threshold of well A (2500 cm-1 above its
        // bottom): the barrier would lie below the well, nothing could be stabilized, and the
        // intermediate steady state is not defined.
        let network = two_well_network();
        let options = ChemicalActivationOptions {
            collision_model: CollisionModel::Stepladder,
            steady_state: SteadyState::Intermediate { barrier: AbsorbingBarrier::default() },
        };
        let hot = Conditions { temperature_kelvin: 400.0, pressure_torr: 760.0 };
        let err = assemble_operator(&network, &hot, &options).unwrap_err();
        assert!(err.contains("'A'"), "{err}");
    }

    #[test]
    fn a_smaller_barrier_distance_moves_the_absorbing_barrier_up() {
        // User choice for shallow wells: 5 k_BT instead of 10 k_BT below the lowest threshold.
        let network = two_well_network();
        let options = ChemicalActivationOptions {
            collision_model: CollisionModel::Stepladder,
            steady_state: SteadyState::Intermediate { barrier: AbsorbingBarrier::BelowLowestThreshold { kt_multiple: 5.0 } },
        };
        let op = assemble_operator(&network, &conditions(), &options).unwrap();
        let five_kt = (5.0 * KB_CM * 300.0 / D_E).round() as usize;
        assert_eq!(op.wells[0].absorbing_barrier_grain, 250 - five_kt);
        assert_eq!(op.wells[1].absorbing_barrier_grain, 310 - five_kt);
        // At 400 K the default 10 k_BT lies below the bottom of well A, 5 k_BT does not.
        let hot = Conditions { temperature_kelvin: 400.0, pressure_torr: 760.0 };
        assert!(assemble_operator(&network, &hot, &options).is_ok());
    }

    #[test]
    fn final_steady_state_keeps_every_grain_and_has_no_stabilization() {
        let network = two_well_network();
        let options = ChemicalActivationOptions {
            collision_model: CollisionModel::Stepladder,
            steady_state: SteadyState::Final,
        };
        let op = assemble_operator(&network, &conditions(), &options).unwrap();
        assert_eq!(op.dimension(), 400 + 460);
        assert!(op.stabilization.iter().all(|s| s.is_empty()));
        assert_eq!(op.wells[0].step_grains, 20); // <dE_down>(300 K) = 200 cm-1 = 20 grains
    }

    #[test]
    fn states_are_ordered_by_absolute_energy() {
        let network = two_well_network();
        let options = all_option_combinations().remove(0);
        let op = assemble_operator(&network, &conditions(), &options).unwrap();
        let energies: Vec<f64> = op.states.iter().map(|&(w, i)| network.absolute_energy_cm1(w, i)).collect();
        assert!(energies.windows(2).all(|p| p[0] <= p[1]));
        // Isomerization couples equal energies: neighbours in the ordering.
        for (r, row) in op.rows.iter().enumerate() {
            for &(c, _) in row {
                if op.states[r].0 != op.states[c].0 {
                    assert!(r.abs_diff(c) <= 1, "inter-well coupling {r}-{c} is not at equal energy");
                }
            }
        }
    }

    #[test]
    fn isomerization_detailed_balance_check_flags_inconsistent_reverse_rates() {
        let mut network = two_well_network();
        let balance = isomerization_detailed_balance(&network);
        assert_eq!(balance.len(), 1);
        assert!(balance[0].max_relative_deviation < 1e-12);
        for k in network.wells[1].channels[1].rate_constant_s_inv.iter_mut() {
            *k *= 2.0;
        }
        let balance = isomerization_detailed_balance(&network);
        assert!((balance[0].max_relative_deviation - 0.5).abs() < 1e-12);
    }

    #[test]
    fn isomerization_beyond_the_target_grid_is_an_error() {
        let mut network = two_well_network();
        // Shorten B by 10 grains: A's grains 390..399 now isomerize to energies B does not have.
        let b = &mut network.wells[1];
        b.density_of_states.truncate(450);
        for channel in b.channels.iter_mut() {
            channel.rate_constant_s_inv.truncate(450);
        }
        let options = all_option_combinations().remove(0);
        assert!(assemble_operator(&network, &conditions(), &options).is_err());
    }
}
