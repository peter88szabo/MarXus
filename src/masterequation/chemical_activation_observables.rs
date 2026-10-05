//! Observables of the steady-state chemical-activation master equation.
//!
//! With a normalized source F and R = 1, N^s = J^-1 F (PO14 eq. 5) and
//!   yield of product channel r          Phi_r = sum_i k_r(E_i) N_i          (O02 eq. 10; O91 eq. 5)
//!   chemical-activation rate coefficient k_r^ca = sum_i k_r(E_i) Ñ_i,  Ñ = N / sum_i N_i over the well
//!                                                                         (GO10 eqs. 8-9; PO14 eq. 6)
//!   bimolecular sink yield              Phi_sink = k_c[D] sum_i N_i        (O02 eq. 5; GO10 eq. 10 gives
//!                                         the relative yield k_c[D]/(sum_r k_r^ca + k_c[D]))
//!   stabilization yield (intermediate steady state)
//!                                       Phi_stab = sum_j N_j (flux into the absorbing region from j)
//!                                                  + source fraction formed below the barrier
//! For isomerization channels sum_i k_r N_i is the gross flux into the other well (internal flux).
//! Every molecule formed leaves the network through a product channel, the sink or (intermediate steady
//! state) the absorbing barrier, so sum Phi_r + sum Phi_stab + sum Phi_sink = 1; this mass balance is
//! returned as a check. In the final steady state without a bimolecular sink there is no net
//! stabilization and the product yields sum to one ("one trivially has R1 = R2 and hence Phi2 = 1",
//! O02, text after eq. 13).
//!
//! References: see `chemical_activation_network.rs`.

use super::chemical_activation_network::{ChannelDestination, ChemicalActivationNetwork};
use super::chemical_activation_operator::ChemicalActivationOperator;
use super::chemical_activation_steady_state::{ProjectedSource, SteadyStateSolution};

/// Results of one unimolecular channel.
#[derive(Debug, Clone)]
pub struct ChannelResult {
    pub well: usize,
    pub channel: usize,
    pub name: String,
    pub destination: ChannelDestination,
    /// sum_i k_r(E_i) N_i: the yield of a product channel, the gross internal flux of an isomerization.
    pub flux: f64,
    /// k_r^ca in s-1 (NaN if the well carries no population).
    pub ca_rate_constant_s_inv: f64,
}

/// Results of one well.
#[derive(Debug, Clone)]
pub struct WellResult {
    pub name: String,
    /// sum_i N_i per unit formation rate (s); the mean residence time in the well.
    pub population: f64,
    /// population / total population of all wells.
    pub population_fraction: f64,
    /// Fraction of the formed molecules stabilized in this well (intermediate steady state).
    pub stabilization_yield: f64,
    /// Fraction of the formed molecules removed by the bimolecular sink k_c[D] of this well.
    pub bimolecular_sink_yield: f64,
    /// Normalized steady-state distribution Ñ on the complete grid of the well (zero below the barrier).
    pub distribution: Vec<f64>,
    /// Mean energy above the well bottom of Ñ, cm-1 (NaN if the well carries no population).
    pub mean_energy_cm1: f64,
}

/// All observables of one steady-state solution.
#[derive(Debug, Clone)]
pub struct ChemicalActivationResult {
    pub channels: Vec<ChannelResult>,
    pub wells: Vec<WellResult>,
    /// sum of all product, stabilization and sink yields (1 up to the solver residual).
    pub mass_balance: f64,
}

impl ChemicalActivationResult {
    /// Sum of the yields of all product channels.
    pub fn total_product_yield(&self) -> f64 {
        self.channels
            .iter()
            .filter(|c| matches!(c.destination, ChannelDestination::Products { .. }))
            .map(|c| c.flux)
            .sum()
    }

    pub fn total_stabilization_yield(&self) -> f64 {
        self.wells.iter().map(|w| w.stabilization_yield).sum()
    }
}

/// Evaluate the observables from the operator, the projected source and the steady-state populations.
pub fn evaluate_observables(
    network: &ChemicalActivationNetwork,
    op: &ChemicalActivationOperator,
    source: &ProjectedSource,
    solution: &SteadyStateSolution,
) -> Result<ChemicalActivationResult, String> {
    let n = op.dimension();
    if solution.population.len() != n || source.on_states.len() != n {
        return Err("Observables: population, source and operator dimensions differ.".into());
    }
    if source.absorbed_per_well.len() != network.wells.len() {
        return Err("Observables: source and network have different numbers of wells.".into());
    }

    // Populations on the complete grid of every well (zero in absorbed grains).
    let mut populations: Vec<Vec<f64>> = network.wells.iter().map(|w| vec![0.0; w.grain_count()]).collect();
    for (s, &(w, i)) in op.states.iter().enumerate() {
        populations[w][i] = solution.population[s];
    }
    let well_population: Vec<f64> = populations.iter().map(|p| p.iter().sum()).collect();
    let total_population: f64 = well_population.iter().sum();

    // Stabilization: flux into the absorbing region of each well, plus direct formation below it.
    let mut stabilization_yield = source.absorbed_per_well.clone();
    for (s, targets) in op.stabilization.iter().enumerate() {
        for &(w, rate) in targets {
            stabilization_yield[w] += rate * solution.population[s];
        }
    }

    let mut channels = Vec::new();
    let mut mass_balance = 0.0;
    for (w, well) in network.wells.iter().enumerate() {
        for (r, channel) in well.channels.iter().enumerate() {
            // sum_i k_r(E_i) N_i (O02 eq. 10) and its average over the normalized distribution (GO10 eq. 9).
            let flux: f64 = channel.rate_constant_s_inv.iter().zip(&populations[w]).map(|(k, p)| k * p).sum();
            if matches!(channel.destination, ChannelDestination::Products { .. }) {
                mass_balance += flux;
            }
            channels.push(ChannelResult {
                well: w,
                channel: r,
                name: channel.name.clone(),
                destination: channel.destination.clone(),
                flux,
                ca_rate_constant_s_inv: if well_population[w] > 0.0 { flux / well_population[w] } else { f64::NAN },
            });
        }
    }

    let mut wells = Vec::with_capacity(network.wells.len());
    for (w, well) in network.wells.iter().enumerate() {
        let population = well_population[w];
        let bimolecular_sink_yield = well.bimolecular_sink_s_inv * population;
        mass_balance += bimolecular_sink_yield + stabilization_yield[w];
        let distribution: Vec<f64> = if population > 0.0 {
            populations[w].iter().map(|p| p / population).collect()
        } else {
            vec![0.0; well.grain_count()]
        };
        let mean_energy_cm1 = if population > 0.0 {
            distribution.iter().enumerate().map(|(i, x)| i as f64 * network.grain_width_cm1 * x).sum()
        } else {
            f64::NAN
        };
        wells.push(WellResult {
            name: well.name.clone(),
            population,
            population_fraction: if total_population > 0.0 { population / total_population } else { 0.0 },
            stabilization_yield: stabilization_yield[w],
            bimolecular_sink_yield,
            distribution,
            mean_energy_cm1,
        });
    }

    Ok(ChemicalActivationResult { channels, wells, mass_balance })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::masterequation::chemical_activation_network::tests::test_well;
    use crate::masterequation::chemical_activation_network::{
        AbsorbingBarrier, Channel, ChemicalActivationOptions, CollisionModel, Conditions, SteadyState,
    };
    use crate::masterequation::chemical_activation_operator::assemble_operator;
    use crate::masterequation::chemical_activation_operator::tests::{all_option_combinations, conditions, two_well_network};
    use crate::masterequation::chemical_activation_steady_state::{project_source, solve_steady_state, LinearSolver};

    fn hot_source_in_first_well(network: &ChemicalActivationNetwork, first: usize) -> Vec<Vec<f64>> {
        network
            .wells
            .iter()
            .enumerate()
            .map(|(w, well)| {
                (0..well.grain_count())
                    .map(|i| if w == 0 && i >= first { (-((i - first) as f64) / 15.0).exp() } else { 0.0 })
                    .collect()
            })
            .collect()
    }

    fn run(
        network: &ChemicalActivationNetwork,
        conditions: &Conditions,
        options: &ChemicalActivationOptions,
        source: &[Vec<f64>],
    ) -> ChemicalActivationResult {
        let op = assemble_operator(network, conditions, options).unwrap();
        let projected = project_source(&op, source).unwrap();
        let solution = solve_steady_state(&op, &projected.on_states, &LinearSolver::BandedCholesky).unwrap();
        evaluate_observables(network, &op, &projected, &solution).unwrap()
    }

    /// Single well with two product channels opening at grains 200 and 260.
    fn single_well_two_channels() -> ChemicalActivationNetwork {
        let mut well = test_well("A", 400, 0, 200);
        well.channels.push(Channel {
            name: "A-second".into(),
            destination: ChannelDestination::Products { name: "Q".into() },
            threshold_grain: None,
            rate_constant_s_inv: (0..400).map(|i| if i >= 260 { 5.0e6 * ((i - 260) as f64 + 1.0).powf(1.5) } else { 0.0 }).collect(),
        });
        ChemicalActivationNetwork { grain_width_cm1: 10.0, wells: vec![well] }
    }

    #[test]
    fn yields_of_products_stabilization_and_sink_add_up_to_one() {
        let network = two_well_network();
        let source = hot_source_in_first_well(&network, 320);
        for options in all_option_combinations() {
            let result = run(&network, &conditions(), &options, &source);
            assert!((result.mass_balance - 1.0).abs() < 1e-9, "{options:?}: mass balance {}", result.mass_balance);
            let sink = &result.wells[1];
            assert!(sink.bimolecular_sink_yield > 0.0);
            assert!((sink.bimolecular_sink_yield - 1.0e5 * sink.population).abs() < 1e-12);
        }
    }

    #[test]
    fn zero_pressure_yields_are_the_rrkm_branching_of_the_nascent_distribution() {
        let network = single_well_two_channels();
        let source = hot_source_in_first_well(&network, 280);
        // 250 K: 10 k_BT (1738 cm-1) lies below the lowest threshold (2000 cm-1).
        let low = Conditions { temperature_kelvin: 250.0, pressure_torr: 1.0e-7 };
        let well = &network.wells[0];
        let f: Vec<f64> = {
            let total: f64 = source[0].iter().sum();
            source[0].iter().map(|x| x / total).collect()
        };
        for options in all_option_combinations() {
            let result = run(&network, &low, &options, &source);
            for r in 0..2 {
                let expected: f64 = (0..well.grain_count())
                    .filter(|&i| f[i] > 0.0)
                    .map(|i| {
                        let total: f64 = well.channels.iter().map(|c| c.rate_constant_s_inv[i]).sum();
                        f[i] * well.channels[r].rate_constant_s_inv[i] / total
                    })
                    .sum();
                let got = result.channels[r].flux;
                assert!(((got - expected) / expected).abs() < 1e-4, "{options:?}: Phi_{r} = {got}, expected {expected}");
            }
        }
    }

    #[test]
    fn high_pressure_intermediate_steady_state_stabilizes_everything() {
        let network = two_well_network();
        let source = hot_source_in_first_well(&network, 320);
        let high = Conditions { temperature_kelvin: 300.0, pressure_torr: 1.0e9 };
        for collision_model in [CollisionModel::ExponentialDown { cutoff_in_mean_down: 10.0 }, CollisionModel::Stepladder] {
            let options = ChemicalActivationOptions {
                collision_model,
                steady_state: SteadyState::Intermediate { barrier: AbsorbingBarrier::default() },
            };
            let result = run(&network, &high, &options, &source);
            let stab = result.total_stabilization_yield();
            assert!(stab > 1.0 - 1e-4, "{collision_model:?}: Phi_stab = {stab}");
        }
    }

    #[test]
    fn final_steady_state_without_a_sink_ends_entirely_in_products() {
        let mut network = two_well_network();
        network.wells[1].bimolecular_sink_s_inv = 0.0;
        let source = hot_source_in_first_well(&network, 320);
        for collision_model in [CollisionModel::ExponentialDown { cutoff_in_mean_down: 10.0 }, CollisionModel::Stepladder] {
            let options = ChemicalActivationOptions { collision_model, steady_state: SteadyState::Final };
            let result = run(&network, &conditions(), &options, &source);
            assert!((result.total_product_yield() - 1.0).abs() < 1e-9, "{collision_model:?}: {}", result.total_product_yield());
            assert_eq!(result.total_stabilization_yield(), 0.0);
        }
    }

    #[test]
    fn final_steady_state_through_the_only_channel_is_the_equilibrium_distribution() {
        // Thermal formation through the reverse of the only channel, no sink: F ∝ K f with
        // f = rho exp(-E/kT), and (I - P) f = 0 for a normalized, detailed-balanced kernel, so J f ∝ F:
        // the final steady state is chemical equilibrium at every pressure, and k^ca is the
        // high-pressure rate coefficient sum k f / sum f (canonical average of k(E)).
        use crate::constants::KB_CM;
        use crate::masterequation::chemical_activation_sources::thermal_source_from_rate;
        let network = ChemicalActivationNetwork { grain_width_cm1: 10.0, wells: vec![test_well("A", 400, 0, 200)] };
        let well = &network.wells[0];
        let t = 600.0;
        let kt = KB_CM * t;
        let source = vec![thermal_source_from_rate(&well.density_of_states, &well.channels[0].rate_constant_s_inv, 10.0, kt).unwrap()];
        let f: Vec<f64> = (0..400).map(|i| well.density_of_states[i] * (-(i as f64) * 10.0 / kt).exp()).collect();
        let f_total: f64 = f.iter().sum();
        let k_inf: f64 = f.iter().zip(&well.channels[0].rate_constant_s_inv).map(|(f, k)| f * k).sum::<f64>() / f_total;
        for collision_model in [CollisionModel::ExponentialDown { cutoff_in_mean_down: 10.0 }, CollisionModel::Stepladder] {
            for pressure_torr in [1.0, 760.0] {
                let options = ChemicalActivationOptions { collision_model, steady_state: SteadyState::Final };
                let result = run(&network, &Conditions { temperature_kelvin: t, pressure_torr }, &options, &source);
                for i in 0..400 {
                    let expected = f[i] / f_total;
                    let got = result.wells[0].distribution[i];
                    assert!((got - expected).abs() <= 1e-8 * expected.max(1e-12), "{collision_model:?}, {pressure_torr} Torr, grain {i}: {got:e} vs {expected:e}");
                }
                assert!((result.channels[0].ca_rate_constant_s_inv / k_inf - 1.0).abs() < 1e-8);
            }
        }
    }

    #[test]
    fn rate_coefficients_average_k_over_the_normalized_well_distribution() {
        let network = two_well_network();
        let source = hot_source_in_first_well(&network, 320);
        let options = all_option_combinations().remove(1);
        let result = run(&network, &conditions(), &options, &source);
        let fractions: f64 = result.wells.iter().map(|w| w.population_fraction).sum();
        assert!((fractions - 1.0).abs() < 1e-14);
        for channel in &result.channels {
            let well = &result.wells[channel.well];
            assert!((well.distribution.iter().sum::<f64>() - 1.0).abs() < 1e-12);
            let k_avg: f64 = network.wells[channel.well].channels[channel.channel]
                .rate_constant_s_inv
                .iter()
                .zip(&well.distribution)
                .map(|(k, n)| k * n)
                .sum();
            assert!((channel.ca_rate_constant_s_inv - k_avg).abs() <= 1e-12 * k_avg);
            assert!((channel.flux - k_avg * well.population).abs() <= 1e-12 * channel.flux);
        }
    }
}
