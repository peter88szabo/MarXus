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
    CollisionModel, Conditions, SteadyState,
};
use super::collision_kernels::{exponential_down_kernel, stepladder_kernel};
use super::collisional_relaxation::lennard_jones_collision_frequency_s_inv;

/// Collision data of one well at the given conditions.
#[derive(Debug, Clone)]
pub struct WellCollisionData {
    /// Collision frequency omega = Z_LJ [M], s-1 (T77 eqs. 3.1-3.3).
    pub collision_frequency_s_inv: f64,
    /// <dE_down>(T), cm-1.
    pub mean_down_cm1: f64,
    /// Lowest retained grain (0 in the final steady state).
    pub absorbing_barrier_grain: usize,
    /// Exponential down: grains below this carry the low-energy reduction factor.
    pub low_energy_cut_grain: usize,
    /// Exponential down: exponent of the low-energy reduction factor.
    pub reduction_exponent: f64,
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
    /// k_B T in cm-1.
    pub kt_cm1: f64,
}

impl ChemicalActivationOperator {
    pub fn dimension(&self) -> usize {
        self.states.len()
    }

    /// (J x) for a vector over the retained states.
    pub fn apply(&self, x: &[f64]) -> Vec<f64> {
        self.rows.iter().map(|row| row.iter().map(|&(c, v)| v * x[c]).sum()).collect()
    }
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

        // The kernel always covers the complete grid of the well, so that the probabilities out of a
        // retained grain stay normalized when the grains below the barrier are absorbing.
        let kernel = match options.collision_model {
            CollisionModel::ExponentialDown { cutoff_in_mean_down } => {
                if !(cutoff_in_mean_down > 0.0) {
                    return Err("The exponential-down cutoff must be positive.".into());
                }
                let band = (cutoff_in_mean_down * mean_down / d_e).ceil() as usize;
                exponential_down_kernel(&well.density_of_states, d_e, mean_down, kt_cm1, band)
            }
            // Step size dE_SL = <dE_down> (GO10, text before eq. 16; PO14: "the step size, dE_SL, which
            // corresponds to the average energy transferred per down collision").
            CollisionModel::Stepladder => stepladder_kernel(&well.density_of_states, d_e, mean_down, kt_cm1),
        }
        .map_err(|e| format!("Well '{}': {e}", well.name))?;

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
                    threshold.saturating_sub((kt_multiple * kt_cm1 / d_e).round() as usize)
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

        wells.push(WellCollisionData {
            collision_frequency_s_inv: omega,
            mean_down_cm1: mean_down,
            absorbing_barrier_grain: barrier,
            low_energy_cut_grain: kernel.low_energy_cut_grain,
            reduction_exponent: kernel.reduction_exponent,
            step_grains: kernel.step_grains,
        });
        kernels.push(kernel);
    }

    // Retained states, ordered by absolute energy and then by well.
    let mut states: Vec<(usize, usize)> = Vec::new();
    for (w, well) in network.wells.iter().enumerate() {
        states.extend((wells[w].absorbing_barrier_grain..well.grain_count()).map(|i| (w, i)));
    }
    states.sort_by_key(|&(w, i)| (i as isize + network.wells[w].bottom_offset_grains, w));
    let mut index_of: Vec<Vec<Option<usize>>> =
        network.wells.iter().map(|well| vec![None; well.grain_count()]).collect();
    for (s, &(w, i)) in states.iter().enumerate() {
        index_of[w][i] = Some(s);
    }

    // Columns of J (source state c = (w, j)), stored by rows.
    let n = states.len();
    let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    let mut stabilization: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for (c, &(w, j)) in states.iter().enumerate() {
        let well = &network.wells[w];
        let omega = wells[w].collision_frequency_s_inv;

        // Collisions, omega (I - P): the diagonal collects omega sum_{t != j} P(t|j), which equals
        // omega (1 - P(j|j)) of PO14 eq. 2 for a normalized kernel and keeps the column sum exactly
        // equal to the loss out of the network.
        let mut diagonal = 0.0;
        for &(t, p) in &kernels[w].transitions[j] {
            let rate = omega * p;
            diagonal += rate;
            match index_of[w][t] {
                Some(r) => rows[r].push((c, -rate)),
                None => stabilization[c].push((w, rate)),
            }
        }

        // Unimolecular channels, K (PO14 eq. 2): products leave the network, isomerization feeds the
        // grain of the target well at the same absolute energy.
        for channel in &well.channels {
            let k = channel.rate_constant_s_inv[j];
            if k == 0.0 {
                continue;
            }
            diagonal += k;
            if let ChannelDestination::Well { index: target } = channel.destination {
                let t = network.aligned_grain(w, j, target);
                let target_well = &network.wells[target];
                if t < 0 || t as usize >= target_well.grain_count() {
                    return Err(format!(
                        "Channel '{}' of well '{}' is open at grain {j} ({} cm-1 absolute), outside the grid \
                         of well '{}' (absolute {} to {} cm-1). Extend the energy grids to a common top.",
                        channel.name,
                        well.name,
                        network.absolute_energy_cm1(w, j),
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
        .map(|&(w, i)| network.wells[w].density_of_states[i].ln() - network.absolute_energy_cm1(w, i) / kt_cm1)
        .collect();

    Ok(ChemicalActivationOperator {
        states,
        index_of,
        rows,
        log_boltzmann_weight,
        stabilization,
        wells,
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
                    .map(|ch| ch.rate_constant_s_inv[i])
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
