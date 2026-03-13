//! Dedicated single-well steady-state master-equation driver.
//!
//! Why this file exists
//! --------------------
//! The low-level numerical/physics pieces for single-well solving already exist in:
//! - `energy_grained_me` (data model and traits),
//! - `chemical_activation_source` (source construction),
//! - `energy_grained_steady_state` (operator construction + linear solves).
//!
//! This module provides one focused, well-commented orchestration layer that decides:
//! - which source mode we are running,
//! - which low-level solver entry point to call,
//! - and how to return a unified output enum.
//!
//! It intentionally does not re-implement physics kernels; it wires existing validated pieces
//! into a single, explicit "single-well solver" API.

use super::energy_grained_me::{
    EnergyGrid, MicrocanonicalGridData, OlzmannMasterEquationSettings, SourceConstructionChoice,
    SteadyStateChemicalActivationResult, ThermalDissociationResult,
};
use super::energy_grained_steady_state::OlzmannStepladderMasterEquationSolver;

/// Input payload for single-well runs.
///
/// This keeps all run configuration together:
/// - the discretized energy grid,
/// - physical/numerical settings,
/// - and the solve mode (CA / equilibrium-source / thermal dissociation).
pub struct EnergyGrainedSteadyStateInput {
    pub energy_grid: EnergyGrid,
    pub settings: OlzmannMasterEquationSettings,
    pub mode: EnergyGrainedSteadyStateMode,
}

/// Supported single-well solve modes.
///
/// - `ChemicalActivation`: driven steady-state with a user-selected source construction.
/// - `EquilibriumSource`: driven steady-state with F(E) ∝ rho(E) exp(-E/kT).
/// - `ThermalDissociation`: no source term; canonical averaging of microcanonical rates.
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

/// Output variants from a single-well run.
pub enum EnergyGrainedSteadyStateOutput {
    ChemicalActivation(SteadyStateChemicalActivationResult),
    ThermalDissociation(ThermalDissociationResult),
}

/// Explicit single-well solver facade.
///
/// This type owns an `OlzmannStepladderMasterEquationSolver`, then routes each mode to the
/// corresponding low-level method. The goal is to keep call sites small and keep orchestration
/// logic in one dedicated module.
pub struct SingleWellSteadyStateSolver {
    inner: OlzmannStepladderMasterEquationSolver,
}

impl SingleWellSteadyStateSolver {
    /// Construct a single-well solver from grid + settings.
    pub fn new(energy_grid: EnergyGrid, settings: OlzmannMasterEquationSettings) -> Self {
        Self {
            inner: OlzmannStepladderMasterEquationSolver {
                energy_grid,
                settings,
            },
        }
    }

    /// Solve the requested mode using the provided microcanonical data.
    ///
    /// Execution flow:
    /// 1. Match the requested mode.
    /// 2. Dispatch to the matching low-level routine in `energy_grained_steady_state`.
    /// 3. Wrap the result in the corresponding output enum variant.
    pub fn solve(
        &self,
        micro: &dyn MicrocanonicalGridData,
        mode: EnergyGrainedSteadyStateMode,
    ) -> Result<EnergyGrainedSteadyStateOutput, String> {
        match mode {
            EnergyGrainedSteadyStateMode::ChemicalActivation {
                formation_flux,
                source_choice,
            } => {
                // Build/normalize source according to `source_choice` and solve J*N = R_form*F.
                let result = self.inner.solve_steady_state_with_source_choice(
                    micro,
                    formation_flux,
                    &source_choice,
                )?;
                Ok(EnergyGrainedSteadyStateOutput::ChemicalActivation(result))
            }
            EnergyGrainedSteadyStateMode::EquilibriumSource { formation_flux } => {
                // Solve the same operator using an equilibrium-like source distribution.
                let result = self
                    .inner
                    .solve_steady_state_with_equilibrium_source(micro, formation_flux)?;
                Ok(EnergyGrainedSteadyStateOutput::ChemicalActivation(result))
            }
            EnergyGrainedSteadyStateMode::ThermalDissociation => {
                // No source term; compute canonical thermal dissociation rates/distribution.
                let result = self.inner.solve_thermal_dissociation(micro)?;
                Ok(EnergyGrainedSteadyStateOutput::ThermalDissociation(result))
            }
        }
    }
}

/// Convenience function: build solver from input and execute one run.
pub fn run_energy_grained_steady_state(
    input: EnergyGrainedSteadyStateInput,
    micro: &dyn MicrocanonicalGridData,
) -> Result<EnergyGrainedSteadyStateOutput, String> {
    let solver = SingleWellSteadyStateSolver::new(input.energy_grid, input.settings);
    solver.solve(micro, input.mode)
}

/// Alias-style convenience entry point kept for compatibility with existing examples/callers.
pub fn run_master_equation(
    input: EnergyGrainedSteadyStateInput,
    micro: &dyn MicrocanonicalGridData,
) -> Result<EnergyGrainedSteadyStateOutput, String> {
    run_energy_grained_steady_state(input, micro)
}
