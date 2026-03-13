//! Unified microcanonical builder for multiwell master-equation networks.
//!
//! Motivation
//! ----------
//! The multiwell master-equation solver (`steady_state_chemical_activation_me`) is written to
//! assemble and solve the linear system once microcanonical inputs are available:
//!   - well densities of states ρ_w(E_i)
//!   - channel microcanonical rates k_{w→ch}(E_i)
//!
//! Historically, those inputs were provided via an external `MicrocanonicalProvider` trait
//! implementation. This module provides a *native* builder that generates those arrays
//! by iterating over all wells/transition states/products in a network and computing:
//!   - ρ(E) for each well (RRHO-style ro-vibrational counting, classical rotor)
//!   - W‡(E) / N‡(E) for each channel transition state (tight RRHO or loose PST core wrapped
//!     in an RRHO-style vibrational/electronic factor)
//!   - k(E) = c * W‡(E - E0) / ρ(E)  (with E in cm^-1, c in cm/s)
//!
//! Notes on units and conventions
//! ------------------------------
//! - Energy grids are in wavenumbers (cm^-1) with uniform bin width ΔE.
//! - ρ(E) is computed as states per cm^-1 (discretized on bins).
//! - W(E) is computed as a cumulative number of states (dimensionless).
//! - With these conventions, RRKM gives:
//!     k(E) [1/s] = c [cm/s] * W‡(E - E0) / ρ(E)
//!   because converting from per-(cm^-1) to per-energy introduces a factor (h c), and
//!   the 1/h in RRKM cancels leaving a factor c.
//! - Symmetry/chirality/electronic degeneracy are handled as multiplicative factors on the
//!   state counts (same placement as MESS RRHO does for weights/states).
//!
//! This module is intentionally straightforward: it does not attempt maximum-entropy
//! smoothing or other statistical post-processing.

use crate::barrierless::phasespace::phase_space_theory::PhaseSpaceTheoryModel;
use crate::masterequation::energy_grained_me::EnergyGrid;
use crate::masterequation::reaction_network::{
    MicrocanonicalProvider, ReactionChannel, WellDefinition,
};
use crate::rrkm::rrkm_rate::rrkm_rate_from_sum_states_and_density;
use crate::rrkm::sum_and_density::get_rovib_WE_or_rhoE;

/// Minimal microcanonical input model for a species (well or tight TS).
///
/// This is meant to be built from your parsed input (or from `MoleculeStruct`),
/// but is independent of any particular parser.
#[derive(Clone, Debug)]
pub struct SpeciesMicroModel {
    pub name: String,
    /// Harmonic vibrational frequencies (cm^-1). Imaginary frequency should be excluded.
    pub vibrational_frequencies_cm1: Vec<f64>,
    /// Rotational constants (cm^-1). Linear rotors may be provided as `[B]` or `[B,B,0]`.
    pub rotational_constants_cm1: Vec<f64>,
    /// Rotational symmetry number σ (dimensionless).
    pub symmetry_number: f64,
    /// Chirality number (enantiomer count), typically 1 or 2.
    pub chirality_number: f64,
    /// Electronic degeneracy factor (dimensionless).
    pub electronic_degeneracy: f64,
}

impl SpeciesMicroModel {
    fn validate(&self) -> Result<(), String> {
        if self
            .vibrational_frequencies_cm1
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(format!(
                "Species '{}' has non-positive or invalid vibrational frequencies.",
                self.name
            ));
        }
        if self.symmetry_number <= 0.0 || !self.symmetry_number.is_finite() {
            return Err(format!(
                "Species '{}' symmetry_number must be positive and finite.",
                self.name
            ));
        }
        if self.chirality_number <= 0.0 || !self.chirality_number.is_finite() {
            return Err(format!(
                "Species '{}' chirality_number must be positive and finite.",
                self.name
            ));
        }
        if self.electronic_degeneracy <= 0.0 || !self.electronic_degeneracy.is_finite() {
            return Err(format!(
                "Species '{}' electronic_degeneracy must be positive and finite.",
                self.name
            ));
        }
        Ok(())
    }

    fn statistical_weight_factor(&self) -> f64 {
        // Same convention used in thermal code: multiply by chirality and divide by σ.
        (self.chirality_number / self.symmetry_number) * self.electronic_degeneracy
    }
}

/// Transition-state model used to build a channel sum-of-states W‡(E).
#[derive(Clone, Debug)]
pub enum TransitionStateModel {
    /// Conventional "tight" transition state: RRHO state counting (rovib + elec + symmetry).
    TightRRHO { species: SpeciesMicroModel },

    /// Loose capture/association transition state: PST core + RRHO wrapper (vib + elec).
    ///
    /// This mirrors the MESS pattern used in your B12 barrier input:
    /// `RRHO { Core PhaseSpaceTheory { ... } Frequencies[...] ElectronicLevels[...] }`
    PhaseSpaceTheoryRRHO {
        pst_core: PhaseSpaceTheoryModel,
        vibrational_frequencies_cm1: Vec<f64>,
        electronic_degeneracy: f64,
    },
}

/// Channel microcanonical definition (per well channel).
#[derive(Clone, Debug)]
pub struct ChannelMicroModel {
    /// Threshold energy E0 for this channel direction, in cm^-1, relative to the well minimum.
    ///
    /// k(E) will be zero for E < E0.
    pub threshold_energy_cm1: f64,

    /// Transition-state sum-of-states model.
    pub transition_state: TransitionStateModel,
}

/// Built microcanonical arrays for a multiwell network.
#[derive(Clone, Debug)]
pub struct MicrocanonicalNetworkData {
    /// Per well: density of states ρ_w(E_i), i = local grain index.
    pub rho_by_well: Vec<Vec<f64>>,

    /// Per well, per channel: microcanonical rate k(E_i) in s^-1.
    ///
    /// This follows the channel ordering in the corresponding `WellDefinition.channels`.
    pub k_by_well_by_channel: Vec<Vec<Vec<f64>>>,

    /// Mirror of the channel bookkeeping (targets) so downstream code can build matrices.
    pub channels_by_well: Vec<Vec<ReactionChannel>>,
}

/// A `MicrocanonicalProvider` backed by the precomputed arrays.
pub struct ArrayMicrocanonicalProvider {
    rho_by_well: Vec<Vec<f64>>,
    k_by_well_by_channel: Vec<Vec<Vec<f64>>>,
}

impl ArrayMicrocanonicalProvider {
    pub fn new(data: &MicrocanonicalNetworkData) -> Self {
        Self {
            rho_by_well: data.rho_by_well.clone(),
            k_by_well_by_channel: data.k_by_well_by_channel.clone(),
        }
    }
}

impl MicrocanonicalProvider for ArrayMicrocanonicalProvider {
    fn density_of_states(&self, well_index: usize, local_grain_index: usize) -> f64 {
        self.rho_by_well
            .get(well_index)
            .and_then(|v| v.get(local_grain_index))
            .copied()
            .unwrap_or(0.0)
    }

    fn microcanonical_rate(
        &self,
        well_index: usize,
        channel_index: usize,
        local_grain_index: usize,
    ) -> f64 {
        self.k_by_well_by_channel
            .get(well_index)
            .and_then(|vv| vv.get(channel_index))
            .and_then(|v| v.get(local_grain_index))
            .copied()
            .unwrap_or(0.0)
    }
}

/// Build ρ(E) and k(E) for a multiwell network.
///
/// Requirements (current implementation):
/// - Each well may have its own ΔE and truncation, but within a well the energy grid is uniform.
/// - `channel_micro_models[well][ch]` must match the well's channel count and ordering.
pub fn build_microcanonical_network_data(
    wells: &[WellDefinition],
    well_models: &[SpeciesMicroModel],
    channel_micro_models: &[Vec<ChannelMicroModel>],
) -> Result<MicrocanonicalNetworkData, String> {
    if wells.len() != well_models.len() || wells.len() != channel_micro_models.len() {
        return Err(
            "wells, well_models, and channel_micro_models must have matching lengths.".into(),
        );
    }

    // Validate wells + models.
    for (idx, (w, m)) in wells.iter().zip(well_models.iter()).enumerate() {
        if w.channels.len() != channel_micro_models[idx].len() {
            return Err(format!(
                "Well '{}' channel count mismatch: network has {}, micro models provide {}.",
                w.well_name,
                w.channels.len(),
                channel_micro_models[idx].len()
            ));
        }
        m.validate()?;
    }

    // 1) Densities of states ρ_w(E)
    let mut rho_by_well: Vec<Vec<f64>> = Vec::with_capacity(wells.len());
    for (well, model) in wells.iter().zip(well_models.iter()) {
        let grid = EnergyGrid {
            number_of_bins: well.one_past_highest_included_grain_index,
            bin_width_wavenumber: well.energy_grain_width_cm1,
            energy_origin_wavenumber: 0.0,
        };
        let rho = compute_rrho_density_of_states(&grid, model)?;
        rho_by_well.push(rho);
    }

    // 2) k(E) per channel, per well
    let mut k_by_well_by_channel: Vec<Vec<Vec<f64>>> = Vec::with_capacity(wells.len());
    for (well_index, well) in wells.iter().enumerate() {
        let model = &well_models[well_index];
        let grid = EnergyGrid {
            number_of_bins: well.one_past_highest_included_grain_index,
            bin_width_wavenumber: well.energy_grain_width_cm1,
            energy_origin_wavenumber: 0.0,
        };

        let rho = &rho_by_well[well_index];
        let mut per_channel: Vec<Vec<f64>> = Vec::with_capacity(well.channels.len());

        for ch in 0..well.channels.len() {
            let ch_model = &channel_micro_models[well_index][ch];
            let w_ts = compute_transition_state_sum_of_states(&grid, &ch_model.transition_state)?;
            let k = rrkm_rate_from_sum_states_and_density(
                grid.bin_width_wavenumber,
                &w_ts,
                rho,
                ch_model.threshold_energy_cm1,
            )?;
            per_channel.push(k);
        }

        k_by_well_by_channel.push(per_channel);
        let _ = model;
    }

    Ok(MicrocanonicalNetworkData {
        rho_by_well,
        k_by_well_by_channel,
        channels_by_well: wells.iter().map(|w| w.channels.clone()).collect(),
    })
}

fn effective_rotational_constants_for_counting(rotational_constants_cm1: &[f64]) -> Vec<f64> {
    // Accept common representations:
    // - atom: [] -> nrot=0
    // - linear: [B] or [B,B,0] -> nrot=2 with [B,B]
    // - nonlinear: [A,B,C] -> nrot=3
    let mut finite: Vec<f64> = rotational_constants_cm1
        .iter()
        .copied()
        .filter(|b| b.is_finite() && *b > 0.0)
        .collect();

    if finite.is_empty() {
        return Vec::new();
    }

    if finite.len() == 1 {
        return vec![finite[0], finite[0]];
    }

    if finite.len() == 2 {
        let b = 0.5 * (finite[0] + finite[1]);
        return vec![b, b];
    }

    finite.truncate(3);
    finite
}

fn compute_rrho_density_of_states(
    grid: &EnergyGrid,
    model: &SpeciesMicroModel,
) -> Result<Vec<f64>, String> {
    let d_e = grid.bin_width_wavenumber;
    let n_ebin = grid.number_of_bins.saturating_sub(1);

    let brot = effective_rotational_constants_for_counting(&model.rotational_constants_cm1);
    let nrot = brot.len();

    let freq_bin: Vec<usize> = model
        .vibrational_frequencies_cm1
        .iter()
        .map(|w| ((*w / d_e) + 0.5).floor().max(1.0) as usize)
        .collect();

    let mut rho = get_rovib_WE_or_rhoE(
        "den".to_string(),
        model.vibrational_frequencies_cm1.len(),
        n_ebin,
        d_e,
        nrot,
        &freq_bin,
        &brot,
    );

    let factor = model.statistical_weight_factor();
    for x in &mut rho {
        *x *= factor;
    }

    // Ensure non-negative and finite.
    for (i, x) in rho.iter().enumerate() {
        if !x.is_finite() || *x < 0.0 {
            return Err(format!(
                "Invalid density of states at bin {} for species '{}'.",
                i, model.name
            ));
        }
    }

    Ok(rho)
}

fn compute_transition_state_sum_of_states(
    grid: &EnergyGrid,
    ts: &TransitionStateModel,
) -> Result<Vec<f64>, String> {
    match ts {
        TransitionStateModel::TightRRHO { species } => compute_rrho_sum_of_states(grid, species),
        TransitionStateModel::PhaseSpaceTheoryRRHO {
            pst_core,
            vibrational_frequencies_cm1,
            electronic_degeneracy,
        } => {
            let d_e = grid.bin_width_wavenumber;
            let n_ebin = grid.number_of_bins.saturating_sub(1);

            // Start from the PST core cumulative states N(E) on the same energy grid.
            let mut w = vec![0.0; n_ebin + 1];
            for i in 1..=n_ebin {
                let e_cm1 = (i as f64) * d_e;
                w[i] = pst_core.cumulative_states_at_energy_cm1(e_cm1)?;
            }

            // RRHO-style electronic degeneracy is a simple multiplication.
            if *electronic_degeneracy <= 0.0 || !electronic_degeneracy.is_finite() {
                return Err("TS electronic_degeneracy must be positive and finite.".into());
            }
            for x in &mut w {
                *x *= *electronic_degeneracy;
            }

            // RRHO-style vibrational convolution: this mirrors the MESS loop:
            //   for each frequency (as integer bins): for e>=bin: W[e] += W[e-bin]
            let freq_bins: Vec<usize> = vibrational_frequencies_cm1
                .iter()
                .filter(|w| w.is_finite() && **w > 0.0)
                .map(|w| ((*w / d_e) + 0.5).floor().max(1.0) as usize)
                .collect();
            convolve_vibrational_sum_states_in_place(&mut w, &freq_bins);

            Ok(w)
        }
    }
}

fn compute_rrho_sum_of_states(
    grid: &EnergyGrid,
    model: &SpeciesMicroModel,
) -> Result<Vec<f64>, String> {
    let d_e = grid.bin_width_wavenumber;
    let n_ebin = grid.number_of_bins.saturating_sub(1);

    let brot = effective_rotational_constants_for_counting(&model.rotational_constants_cm1);
    let nrot = brot.len();

    let freq_bin: Vec<usize> = model
        .vibrational_frequencies_cm1
        .iter()
        .map(|w| ((*w / d_e) + 0.5).floor().max(1.0) as usize)
        .collect();

    let mut w = get_rovib_WE_or_rhoE(
        "sum".to_string(),
        model.vibrational_frequencies_cm1.len(),
        n_ebin,
        d_e,
        nrot,
        &freq_bin,
        &brot,
    );

    let factor = model.statistical_weight_factor();
    for x in &mut w {
        *x *= factor;
    }

    for (i, x) in w.iter().enumerate() {
        if !x.is_finite() || *x < 0.0 {
            return Err(format!(
                "Invalid sum of states at bin {} for species '{}'.",
                i, model.name
            ));
        }
    }

    Ok(w)
}

fn convolve_vibrational_sum_states_in_place(sum_states: &mut [f64], mode_bins: &[usize]) {
    if sum_states.is_empty() {
        return;
    }
    let n = sum_states.len() - 1;
    for &bin in mode_bins {
        if bin == 0 || bin > n {
            continue;
        }
        for e in bin..=n {
            let add = sum_states[e - bin];
            sum_states[e] += add;
        }
    }
}
