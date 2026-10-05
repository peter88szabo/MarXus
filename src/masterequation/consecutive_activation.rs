//! Consecutive chemical activation: coupling of master equations (PO14 pp. 236-237).
//!
//! An intermediate C1 formed with a non-thermal distribution reacts in a pseudo-first-order bimolecular
//! step with a thermal partner B (the sink k_c[D] of the first master equation) to a second intermediate
//! C2, which is itself chemically activated. "The output of the first master equation governs the input
//! of a second one. This is what we mean by coupling of master equations" (PO14 p. 237). With the
//! normalized steady-state distribution ñ1^ss(E) of C1 (energy independent k_c, so the molecules that
//! react with B have the distribution ñ1^ss) and the 0 K reaction energy RE of the bimolecular step,
//!   f2(E) = integral_0^{E+RE} ñ_B(e) ñ1^ss(E + RE - e) de                   (PO14 eqs. 10-11)
//! or, neglecting the width of the thermal partner distribution, ñ_B(e) = delta(e - <E_B>),
//!   f2(E) = ñ1^ss(E + RE - <E_B>)                                          (PO14 eqs. 12-13)
//! with <E_B> the average thermal energy of the partner. The energy zero of each master equation is the
//! ground state of its own intermediate (PO14 p. 237).
//!
//! References: see `chemical_activation_network.rs`.

use crate::constants::KB_CM;

use super::chemical_activation_driver::{run_chemical_activation, ChemicalActivationRun, ConditionResult, SourceSpecification};
use super::chemical_activation_network::ChemicalActivationNetwork;
use super::chemical_activation_observables::ChemicalActivationResult;
use super::chemical_activation_sources::{convolution_source, shifted_source, thermal_distribution};

/// Treatment of the internal energy of the thermal partner B.
#[derive(Debug, Clone)]
pub enum PartnerTreatment {
    /// Convolution with the thermal distribution rho_B(E) exp(-E/kT) of the partner (PO14 eqs. 10-11);
    /// rho_B on the grain grid of the networks.
    Convolution { partner_density_of_states: Vec<f64> },
    /// Shift by -RE + <E_B> with the given mean thermal energy of the partner (PO14 eqs. 12-13).
    Shift { partner_mean_energy_cm1: f64 },
}

/// The bimolecular step C1 + B -> C2 linking two master equations.
#[derive(Debug, Clone)]
pub struct ConsecutiveStep {
    /// Well of the first network whose steady-state distribution reacts with B.
    pub from_well: usize,
    /// Well of the second network that is formed.
    pub to_well: usize,
    /// 0 K reaction energy RE of C1 + B -> C2 in cm-1 (negative for an exothermic step).
    pub reaction_energy_cm1: f64,
    pub partner: PartnerTreatment,
}

/// Source of the second master equation (on the grids of all its wells) from the steady state of the first.
pub fn consecutive_source(
    first: &ChemicalActivationResult,
    second: &ChemicalActivationNetwork,
    step: &ConsecutiveStep,
    temperature_kelvin: f64,
) -> Result<Vec<Vec<f64>>, String> {
    let n1 = &first
        .wells
        .get(step.from_well)
        .ok_or_else(|| format!("Consecutive step: well {} not in the first network.", step.from_well))?
        .distribution;
    let receiving = second
        .wells
        .get(step.to_well)
        .ok_or_else(|| format!("Consecutive step: well {} not in the second network.", step.to_well))?;
    if !(step.reaction_energy_cm1 <= 0.0) {
        // PO14 eqs. 10-13 treat an exothermic step with an energy-independent reactive cross section.
        return Err(format!(
            "Consecutive step: reaction energy {} cm-1; the coupling of master equations of PO14 eqs. 10-13 \
             applies to exothermic steps (RE <= 0).",
            step.reaction_energy_cm1
        ));
    }
    let d_e = second.grain_width_cm1;
    let grains = receiving.grain_count();
    let f2 = match &step.partner {
        PartnerTreatment::Convolution { partner_density_of_states } => {
            let n_b = thermal_distribution(partner_density_of_states, d_e, KB_CM * temperature_kelvin)?;
            // E0 = -RE: the energy of C1 + B above the ground state of C2, in whole grains.
            let threshold = (-step.reaction_energy_cm1 / d_e).round() as usize;
            convolution_source(grains, threshold, &n_b, n1)?
        }
        PartnerTreatment::Shift { partner_mean_energy_cm1 } => {
            shifted_source(grains, n1, step.reaction_energy_cm1, *partner_mean_energy_cm1, d_e)?
        }
    };
    Ok(second
        .wells
        .iter()
        .enumerate()
        .map(|(w, well)| if w == step.to_well { f2.clone() } else { vec![0.0; well.grain_count()] })
        .collect())
}

/// Run the first network, then at every (T, p) the second network with the source formed from the
/// first one's steady state. Both runs use the temperatures and pressures of `first_run`; the source
/// specification of `second_run` is replaced.
pub fn run_consecutive_activation(
    first_network: &ChemicalActivationNetwork,
    first_run: &ChemicalActivationRun,
    second_network: &ChemicalActivationNetwork,
    second_run: &ChemicalActivationRun,
    step: &ConsecutiveStep,
) -> Result<Vec<(ConditionResult, ConditionResult)>, String> {
    let first_results = run_chemical_activation(first_network, first_run)?;
    let mut chain = Vec::with_capacity(first_results.len());
    for one in first_results {
        let t = one.conditions.temperature_kelvin;
        let source = consecutive_source(&one.result, second_network, step, t)?;
        let single_condition = ChemicalActivationRun {
            temperatures_kelvin: vec![t],
            pressures_torr: vec![one.conditions.pressure_torr],
            source: SourceSpecification::Fixed(source),
            ..second_run.clone()
        };
        let two = run_chemical_activation(second_network, &single_condition)?.remove(0);
        chain.push((one, two));
    }
    Ok(chain)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::masterequation::chemical_activation_driver::SourceSpecification;
    use crate::masterequation::chemical_activation_network::tests::test_well;
    use crate::masterequation::chemical_activation_network::{
        AbsorbingBarrier, ChemicalActivationOptions, CollisionModel, SteadyState,
    };
    use crate::masterequation::chemical_activation_steady_state::LinearSolver;

    const D_E: f64 = 10.0;

    /// First intermediate: formed above grain 300 by a thermal entrance channel (its product channel
    /// doubles as the back reaction), removed by k_c[D] = 1e6 s-1.
    fn first_network() -> ChemicalActivationNetwork {
        let mut well = test_well("C1", 450, 0, 300);
        well.bimolecular_sink_s_inv = 1.0e6;
        ChemicalActivationNetwork { grain_width_cm1: D_E, wells: vec![well] }
    }

    /// Second intermediate: 1400 grains, products from grain 900.
    fn second_network() -> ChemicalActivationNetwork {
        ChemicalActivationNetwork { grain_width_cm1: D_E, wells: vec![test_well("C2", 1400, 0, 900)] }
    }

    fn run(temperatures: Vec<f64>, pressures: Vec<f64>) -> ChemicalActivationRun {
        ChemicalActivationRun {
            temperatures_kelvin: temperatures,
            pressures_torr: pressures,
            options: ChemicalActivationOptions {
                collision_model: CollisionModel::Stepladder,
                steady_state: SteadyState::Final,
            },
            solver: LinearSolver::BandedCholesky,
            source: SourceSpecification::ThermalEntrance { channels: vec![(0, 0)] },
            tolerance: 1e-8,
        }
    }

    fn first_result() -> ChemicalActivationResult {
        run_chemical_activation(&first_network(), &run(vec![298.0], vec![760.0])).unwrap().remove(0).result
    }

    #[test]
    fn shift_treatment_moves_the_steady_state_distribution_by_minus_re_plus_partner_energy() {
        let first = first_result();
        let step = ConsecutiveStep {
            from_well: 0,
            to_well: 0,
            reaction_energy_cm1: -6630.0,
            partner: PartnerTreatment::Shift { partner_mean_energy_cm1: 210.0 },
        };
        let f2 = consecutive_source(&first, &second_network(), &step, 298.0).unwrap();
        let n1 = &first.wells[0].distribution;
        for (i, &x) in f2[0].iter().enumerate() {
            let expected = if i >= 684 && i - 684 < n1.len() { n1[i - 684] } else { 0.0 };
            assert!((x - expected).abs() <= 1e-14, "grain {i}: {x} vs {expected}");
        }
    }

    #[test]
    fn convolution_treatment_uses_the_thermal_partner_distribution() {
        let first = first_result();
        let rho_b: Vec<f64> = (0..200).map(|i| 1.0 + 0.5 * i as f64).collect();
        let step = ConsecutiveStep {
            from_well: 0,
            to_well: 0,
            reaction_energy_cm1: -6630.0,
            partner: PartnerTreatment::Convolution { partner_density_of_states: rho_b.clone() },
        };
        let f2 = consecutive_source(&first, &second_network(), &step, 298.0).unwrap();
        let n_b = thermal_distribution(&rho_b, D_E, KB_CM * 298.0).unwrap();
        let expected = convolution_source(1400, 663, &n_b, &first.wells[0].distribution).unwrap();
        assert_eq!(f2[0], expected);
    }

    #[test]
    fn the_second_master_equation_is_run_at_every_condition_of_the_first() {
        let step = ConsecutiveStep {
            from_well: 0,
            to_well: 0,
            reaction_energy_cm1: -6630.0,
            partner: PartnerTreatment::Shift { partner_mean_energy_cm1: 210.0 },
        };
        let first_run = run(vec![298.0, 350.0], vec![10.0, 760.0]);
        // The second intermediate has no sink and a deep well: its final steady state is numerically
        // singular at 298 K, so the intermediate steady state is used.
        let second_run = ChemicalActivationRun {
            options: ChemicalActivationOptions {
                collision_model: CollisionModel::Stepladder,
                steady_state: SteadyState::Intermediate { barrier: AbsorbingBarrier::default() },
            },
            ..run(vec![], vec![])
        };
        let chain = run_consecutive_activation(&first_network(), &first_run, &second_network(), &second_run, &step).unwrap();
        assert_eq!(chain.len(), 4);
        for (one, two) in &chain {
            assert_eq!(one.conditions.temperature_kelvin, two.conditions.temperature_kelvin);
            assert_eq!(one.conditions.pressure_torr, two.conditions.pressure_torr);
            // The second run must equal a run with the source built from this condition's first result.
            let source = consecutive_source(&one.result, &second_network(), &step, one.conditions.temperature_kelvin).unwrap();
            let by_hand = run_chemical_activation(
                &second_network(),
                &ChemicalActivationRun {
                    temperatures_kelvin: vec![one.conditions.temperature_kelvin],
                    pressures_torr: vec![one.conditions.pressure_torr],
                    source: SourceSpecification::Fixed(source),
                    ..second_run.clone()
                },
            )
            .unwrap();
            assert_eq!(two.result.channels[0].flux, by_hand[0].result.channels[0].flux);
        }
    }
}
