//! State counting for RRHO species and transition states on a uniform grain grid.
//!
//! - rho(E) of a species: rovibrational density of states (per cm-1) by direct count of harmonic
//!   vibrations convolved with classical rigid rotors (`rrkm::sum_and_density`) and with the quantum
//!   levels of one-dimensional internal rotors (`rrkm::internal_rotor`).
//! - W‡(E) of a transition state: rovibrational sum of states of a tight transition state, or the
//!   cumulative states of a phase-space-theory core combined with harmonic conserved modes.
//! - Excited electronic levels are convolved with the counts: rho(E) = sum_j g_j rho_0(E - eps_j).
//! - Symmetry number, chirality and electronic degeneracy multiply the state counts, so that
//!   k(E) = W‡(E - E0)/(h rho(E)) carries the ratio of these factors (Forst, Theory of Unimolecular
//!   Reactions (1973), Sec. 4.5). With rho per cm-1 and W dimensionless, h is taken in cm-1 s.
//!
//! Grain i lies at E = i dE; rho[i] holds the states in ((i-1) dE, i dE] per dE and rho[0] the ground
//! state (`rrkm::sum_and_density`).

use crate::barrierless::phasespace::phase_space_theory::PhaseSpaceTheoryModel;
use crate::rrkm::internal_rotor::{convolve_rotor_levels, HinderedRotor};
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
    /// Electronic degeneracy factor of the ground level (dimensionless).
    pub electronic_degeneracy: f64,
    /// Excited electronic levels (energy above the ground level in cm-1, degeneracy): the counts become
    /// rho(E) = g_0 rho_0(E) + sum_j g_j rho_0(E - eps_j), the same for W.
    pub excited_electronic_levels: Vec<(f64, f64)>,
    /// One-dimensional hindered or free internal rotors, whose quantum levels are stick-convolved
    /// with the rovibrational counts (`rrkm::internal_rotor`).
    pub internal_rotors: Vec<HinderedRotor>,
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

    /// Rovibrational counts (densities or sums) combined with the internal-rotor levels, each level
    /// counted from the rotor ground level (the species zero energy contains the rotor zero-point
    /// energy). The levels of every rotor must reach the top of the cells.
    fn convolve_internal_rotors(&self, mut counts: Vec<f64>, cell_cm1: f64) -> Result<Vec<f64>, String> {
        let top_cm1 = counts.len().saturating_sub(1) as f64 * cell_cm1;
        for (r, rotor) in self.internal_rotors.iter().enumerate() {
            if rotor.highest_level_above_ground_cm1() < top_cm1 {
                return Err(format!(
                    "Species '{}': the levels of internal rotor {} end at {:.1} cm-1 above its ground, below the top \
                     of the energy grid at {top_cm1:.0} cm-1; a larger Fourier basis is needed.",
                    self.name,
                    r + 1,
                    rotor.highest_level_above_ground_cm1()
                ));
            }
            counts = convolve_rotor_levels(&counts, &rotor.levels_above_ground_cm1, cell_cm1);
        }
        Ok(counts)
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
        /// Excited electronic levels (energy above the ground level in cm-1, degeneracy).
        excited_electronic_levels: Vec<(f64, f64)>,
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
    let rho = model.convolve_internal_rotors(rho, d_e)?;
    let rho = convolve_electronic_levels(rho, model.electronic_degeneracy, &model.excited_electronic_levels, d_e);

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
            excited_electronic_levels,
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

            Ok(convolve_electronic_levels(w, *electronic_degeneracy, excited_electronic_levels, d_e))
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
    let w = model.convolve_internal_rotors(w, d_e)?;
    let w = convolve_electronic_levels(w, model.electronic_degeneracy, &model.excited_electronic_levels, d_e);

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

/// Counts that already carry the ground-level degeneracy g_0, combined with excited electronic levels:
/// c(E) + sum_j (g_j/g_0) c(E - eps_j), each level shifting by ceil(eps_j/cell) cells (as the rotor levels,
/// `rrkm::internal_rotor::convolve_rotor_levels`).
fn convolve_electronic_levels(counts: Vec<f64>, ground_degeneracy: f64, excited: &[(f64, f64)], cell_cm1: f64) -> Vec<f64> {
    if excited.is_empty() {
        return counts;
    }
    let mut out = counts.clone();
    for &(eps, g) in excited {
        let shift = (eps / cell_cm1 - 1e-9).ceil().max(0.0) as usize;
        let weight = g / ground_degeneracy;
        for i in shift..counts.len() {
            out[i] += weight * counts[i - shift];
        }
    }
    out
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rrkm::internal_rotor::{convolve_rotor_levels, TorsionalPotential, DEFAULT_BASIS_SIZE};

    fn species(internal_rotors: Vec<HinderedRotor>) -> SpeciesMicroModel {
        SpeciesMicroModel {
            name: "X".into(),
            vibrational_frequencies_cm1: vec![812.0, 1430.0, 2950.0],
            rotational_constants_cm1: vec![1.2, 0.31, 0.27],
            symmetry_number: 2.0,
            chirality_number: 1.0,
            electronic_degeneracy: 2.0,
            internal_rotors,
            excited_electronic_levels: Vec::new(),
        }
    }

    #[test]
    fn a_hindered_rotor_is_convolved_into_the_density_and_the_sum_of_states() {
        // threefold methyl-like rotor, V = 350 (1 - cos 3 phi) cm-1
        let potential = TorsionalPotential { constant: 350.0, cosine: vec![-350.0], sine: vec![] };
        let rotor = HinderedRotor::new(5.6, 3, potential, DEFAULT_BASIS_SIZE).unwrap();
        let (grains, cell) = (3000, 1.0);

        let rho_without = rrho_density_of_states(grains, cell, &species(vec![])).unwrap();
        let w_without = rrho_sum_of_states(grains, cell, &species(vec![])).unwrap();
        let rho = rrho_density_of_states(grains, cell, &species(vec![rotor.clone()])).unwrap();
        let w = rrho_sum_of_states(grains, cell, &species(vec![rotor.clone()])).unwrap();

        assert_eq!(rho, convolve_rotor_levels(&rho_without, &rotor.levels_above_ground_cm1, cell));
        assert_eq!(w, convolve_rotor_levels(&w_without, &rotor.levels_above_ground_cm1, cell));
    }

    #[test]
    fn a_rotor_whose_levels_end_below_the_top_of_the_grid_is_an_error() {
        // 5 basis functions: levels up to B (2 sigma)^2 = 201.6 cm-1 only
        let potential = TorsionalPotential { constant: 350.0, cosine: vec![-350.0], sine: vec![] };
        let rotor = HinderedRotor::new(5.6, 3, potential, 5).unwrap();
        let err = rrho_density_of_states(3000, 1.0, &species(vec![rotor.clone()])).unwrap_err();
        assert!(err.contains("X") && err.contains("2999"), "{err}");
        assert!(rrho_sum_of_states(3000, 1.0, &species(vec![rotor])).is_err());
    }

    #[test]
    fn excited_electronic_levels_are_convolved_into_the_counts() {
        // rho(E) = sum_j g_j rho_0(E - eps_j) (each level a shift by whole cells), the same for W; then
        // sum_E rho e^(-E/kT) = q_el q_0 with q_el = sum_j g_j e^(-eps_j/kT).
        let ground = species(vec![]);
        let mut excited = species(vec![]);
        excited.excited_electronic_levels = vec![(140.0, 2.0), (1000.0, 4.0)];
        let (grains, cell, g0) = (20_000, 1.0, 2.0);
        for count in [rrho_density_of_states, rrho_sum_of_states] {
            let base = count(grains, cell, &ground).unwrap();
            let with = count(grains, cell, &excited).unwrap();
            for i in [0, 139, 140, 141, 999, 1000, 5000, 19_999] {
                let shifted = |s: usize| if i >= s { base[i - s] } else { 0.0 };
                let expected = base[i] + 2.0 / g0 * shifted(140) + 4.0 / g0 * shifted(1000);
                assert!((with[i] - expected).abs() <= 1e-12 * expected.abs(), "cell {i}: {} vs {expected}", with[i]);
            }
        }
        let kt = crate::constants::KB_CM * 300.0;
        let z = |r: &[f64]| r.iter().enumerate().map(|(i, x)| x * (-(i as f64) / kt).exp()).sum::<f64>();
        let q_el = 2.0 + 2.0 * (-140.0 / kt).exp() + 4.0 * (-1000.0 / kt).exp();
        let (rho0, rho) = (rrho_density_of_states(grains, cell, &ground).unwrap(), rrho_density_of_states(grains, cell, &excited).unwrap());
        assert!((z(&rho) / z(&rho0) - q_el / g0).abs() < 1e-12);
    }

    #[test]
    fn excited_electronic_levels_of_a_phase_space_theory_transition_state_are_convolved() {
        use crate::barrierless::phasespace::phase_space_theory::PhaseSpaceTheoryModel;
        use crate::barrierless::phasespace::types::{CaptureFragment, CaptureFragmentRotorModel, PhaseSpaceTheoryInput, PstTstLevel};
        let pst = || {
            PhaseSpaceTheoryModel::new(PhaseSpaceTheoryInput {
                fragment_a: CaptureFragment { mass_amu: Some(1.0), rotor: CaptureFragmentRotorModel::Atom },
                fragment_b: CaptureFragment { mass_amu: Some(32.0), rotor: CaptureFragmentRotorModel::LinearRigidRotor { rotational_constant_cm1: 1.44 } },
                symmetry_operations: 1.0,
                potential_prefactor_au: 37.0,
                potential_power_exponent: 6.0,
                tst_level: PstTstLevel::E,
            })
            .unwrap()
        };
        let ts = |excited: Vec<(f64, f64)>| TransitionStateModel::PhaseSpaceTheoryRRHO {
            pst_core: pst(),
            vibrational_frequencies_cm1: vec![1580.0],
            electronic_degeneracy: 6.0,
            excited_electronic_levels: excited,
        };
        let base = transition_state_sum_of_states(3000, 1.0, &ts(vec![])).unwrap();
        let with = transition_state_sum_of_states(3000, 1.0, &ts(vec![(250.0, 3.0)])).unwrap();
        for i in [0, 249, 250, 1000, 2999] {
            let expected = base[i] + if i >= 250 { 3.0 / 6.0 * base[i - 250] } else { 0.0 };
            assert!((with[i] - expected).abs() <= 1e-12 * expected.abs(), "cell {i}");
        }
    }
}
