//! State counting for RRHO species and transition states on a uniform grain grid.
//!
//! - rho(E) of a species: rovibrational density of states (per cm-1) by direct count of harmonic
//!   vibrations convolved with classical rigid rotors (`rrkm::sum_and_density`).
//! - W‡(E) of a transition state: rovibrational sum of states of a tight transition state, or the
//!   cumulative states of a phase-space-theory core combined with harmonic conserved modes.
//! - Symmetry number, chirality and electronic degeneracy multiply the state counts, so that
//!   k(E) = W‡(E - E0)/(h rho(E)) carries the ratio of these factors (Forst, Theory of Unimolecular
//!   Reactions (1973), Sec. 4.5). With rho per cm-1 and W dimensionless, h is taken in cm-1 s.
//!
//! Grain i lies at E = i dE; rho[i] holds the states in ((i-1) dE, i dE] per dE and rho[0] the ground
//! state (`rrkm::sum_and_density`).

use crate::barrierless::phasespace::phase_space_theory::PhaseSpaceTheoryModel;
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
    /// Input-deck form: `RRHO { Core PhaseSpaceTheory { ... } Frequencies[...] ElectronicLevels[...] }`.
    PhaseSpaceTheoryRRHO {
        pst_core: PhaseSpaceTheoryModel,
        vibrational_frequencies_cm1: Vec<f64>,
        electronic_degeneracy: f64,
    },
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

/// Rovibrational density of states (per cm-1) of an RRHO species on `grains` grains of width
/// `grain_width_cm1`, including the statistical weight chirality g_e/sigma.
pub(crate) fn rrho_density_of_states(
    grains: usize,
    grain_width_cm1: f64,
    model: &SpeciesMicroModel,
) -> Result<Vec<f64>, String> {
    model.validate()?;
    let d_e = grain_width_cm1;
    let n_ebin = grains.saturating_sub(1);

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

/// Sum of states W(E) of a transition state on `grains` grains of width `grain_width_cm1`.
pub(crate) fn transition_state_sum_of_states(
    grains: usize,
    grain_width_cm1: f64,
    ts: &TransitionStateModel,
) -> Result<Vec<f64>, String> {
    match ts {
        TransitionStateModel::TightRRHO { species } => rrho_sum_of_states(grains, grain_width_cm1, species),
        TransitionStateModel::PhaseSpaceTheoryRRHO {
            pst_core,
            vibrational_frequencies_cm1,
            electronic_degeneracy,
        } => {
            let d_e = grain_width_cm1;
            let n_ebin = grains.saturating_sub(1);

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

            // Harmonic conserved modes added to the core sum of states by direct count, for each
            // frequency (in whole grains): W[e] += W[e - bin] for e >= bin (Beyer, Swinehart,
            // Commun. ACM 16, 379 (1973)).
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

/// Rovibrational sum of states of an RRHO species, including chirality g_e/sigma.
pub(crate) fn rrho_sum_of_states(
    grains: usize,
    grain_width_cm1: f64,
    model: &SpeciesMicroModel,
) -> Result<Vec<f64>, String> {
    model.validate()?;
    let d_e = grain_width_cm1;
    let n_ebin = grains.saturating_sub(1);

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
