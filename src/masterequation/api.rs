//! High-level, struct-based entry points for running master-equation solvers.
//!
//! Design goals:
//! - No explicit lifetimes in the public API.
//! - No heap ownership wrappers (`Box`, `Arc`) forced on the user.
//! - Keep inputs physics-facing (grid + collision model + source + energetics).
//!
//! Microcanonical data providers are passed as function arguments (`&dyn ...`) at call time.
//! This keeps the API simple while avoiding storing borrowed trait objects inside structs.

use crate::masterequation::input_deck::MasterEquationInput;
use crate::masterequation::reaction_network::{
    ChemicalActivationDefinition, MasterEquationSettings, MicrocanonicalProvider, SolutionResults,
    WellDefinition,
};
pub use crate::masterequation::singlewell_solver::{
    run_energy_grained_steady_state, run_master_equation, EnergyGrainedSteadyStateInput,
    EnergyGrainedSteadyStateMode, EnergyGrainedSteadyStateOutput,
};
use crate::masterequation::steady_state_chemical_activation_me::{
    MasterEquationEngine, MasterEquationSolveDiagnostics,
};

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

pub fn run_multiwell_chemical_activation_with_diagnostics(
    input: MultiWellChemicalActivationInput,
    micro: &dyn MicrocanonicalProvider,
) -> Result<(SolutionResults, MasterEquationSolveDiagnostics), String> {
    let engine = MasterEquationEngine::new(input.wells, input.settings);
    engine.solve_steady_state_chemical_activation_with_diagnostics(micro, input.activation)
}

pub fn run_multiwell_from_network(
    input: MultiWellFromNetworkInput,
    micro: &dyn MicrocanonicalProvider,
) -> Result<SolutionResults, String> {
    let (settings, wells) = input.network.build_settings_and_wells()?;
    let engine = MasterEquationEngine::new(wells, settings);
    engine.solve_steady_state_chemical_activation(micro, input.activation)
}

pub fn run_multiwell_from_network_with_diagnostics(
    input: MultiWellFromNetworkInput,
    micro: &dyn MicrocanonicalProvider,
) -> Result<(SolutionResults, MasterEquationSolveDiagnostics), String> {
    let (settings, wells) = input.network.build_settings_and_wells()?;
    let engine = MasterEquationEngine::new(wells, settings);
    engine.solve_steady_state_chemical_activation_with_diagnostics(micro, input.activation)
}
