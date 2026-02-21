//! High-level, struct-based entry points for running master-equation solvers.
//!
//! Design goals:
//! - No explicit lifetimes in the public API.
//! - No heap ownership wrappers (`Box`, `Arc`) forced on the user.
//! - Keep inputs physics-facing (grid + collision model + source + energetics).
//!
//! Microcanonical data providers are passed as function arguments (`&dyn ...`) at call time.
//! This keeps the API simple while avoiding storing borrowed trait objects inside structs.

use crate::masterequation::energy_grained_me::{
    EnergyGrid, MicrocanonicalGridData, OlzmannMasterEquationSettings, SourceConstructionChoice,
    SteadyStateChemicalActivationResult, ThermalDissociationResult,
};
use crate::masterequation::energy_grained_steady_state::OlzmannStepladderMasterEquationSolver;
use crate::masterequation::input_deck::MasterEquationInput;
use crate::masterequation::reaction_network::{
    ChemicalActivationDefinition, MasterEquationSettings, MicrocanonicalProvider, SolutionResults,
    WellDefinition,
};
use crate::masterequation::steady_state_chemical_activation_me::MasterEquationEngine;

// -----------------------------------------------------------------------------
// Multi-well steady-state chemical activation (dense operator)
// -----------------------------------------------------------------------------

pub struct MultiWellChemicalActivationInput {
    pub wells: Vec<WellDefinition>,
    pub settings: MasterEquationSettings,
    pub activation: ChemicalActivationDefinition,
}

pub struct MultiWellFromNetworkInput {
    pub network: MasterEquationInput,
    pub activation: ChemicalActivationDefinition,
}

pub fn run_multiwell_chemical_activation(
    input: MultiWellChemicalActivationInput,
    micro: &dyn MicrocanonicalProvider,
) -> Result<SolutionResults, String> {
    let engine = MasterEquationEngine::new(input.wells, input.settings);
    engine.solve_steady_state_chemical_activation(micro, input.activation)
}

pub fn run_multiwell_from_network(
    input: MultiWellFromNetworkInput,
    micro: &dyn MicrocanonicalProvider,
) -> Result<SolutionResults, String> {
    let (settings, wells) = input.network.build_settings_and_wells()?;
    let engine = MasterEquationEngine::new(wells, settings);
    engine.solve_steady_state_chemical_activation(micro, input.activation)
}

// -----------------------------------------------------------------------------
// Single-well steady-state chemical activation (energy-grained solver)
// -----------------------------------------------------------------------------

pub struct EnergyGrainedSteadyStateInput {
    pub energy_grid: EnergyGrid,
    pub settings: OlzmannMasterEquationSettings,
    pub mode: EnergyGrainedSteadyStateMode,
}

pub enum EnergyGrainedSteadyStateMode {
    ChemicalActivation {
        formation_flux: f64,
        source_choice: SourceConstructionChoice,
    },
    EquilibriumSource {
        formation_flux: f64,
    },
    ThermalDissociation,
}

pub enum EnergyGrainedSteadyStateOutput {
    ChemicalActivation(SteadyStateChemicalActivationResult),
    ThermalDissociation(ThermalDissociationResult),
}

pub fn run_energy_grained_steady_state(
    input: EnergyGrainedSteadyStateInput,
    micro: &dyn MicrocanonicalGridData,
) -> Result<EnergyGrainedSteadyStateOutput, String> {
    let solver = OlzmannStepladderMasterEquationSolver {
        energy_grid: input.energy_grid,
        settings: input.settings,
    };

    match input.mode {
        EnergyGrainedSteadyStateMode::ChemicalActivation {
            formation_flux,
            source_choice,
        } => {
            let result =
                solver.solve_steady_state_with_source_choice(micro, formation_flux, &source_choice)?;
            Ok(EnergyGrainedSteadyStateOutput::ChemicalActivation(result))
        }
        EnergyGrainedSteadyStateMode::EquilibriumSource { formation_flux } => {
            let result =
                solver.solve_steady_state_with_equilibrium_source(micro, formation_flux)?;
            Ok(EnergyGrainedSteadyStateOutput::ChemicalActivation(result))
        }
        EnergyGrainedSteadyStateMode::ThermalDissociation => {
            let result = solver.solve_thermal_dissociation(micro)?;
            Ok(EnergyGrainedSteadyStateOutput::ThermalDissociation(result))
        }
    }
}

/// Physics-facing convenience entry point.
///
/// This is intentionally "boring": you pass one input struct + one microcanonical provider,
/// and you get the requested steady-state (or thermal) result back.
pub fn run_master_equation(
    input: EnergyGrainedSteadyStateInput,
    micro: &dyn MicrocanonicalGridData,
) -> Result<EnergyGrainedSteadyStateOutput, String> {
    run_energy_grained_steady_state(input, micro)
}
