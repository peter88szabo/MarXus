//! Dissociation rate constants k(E) of an energy-selected ion, on the cells of its internal energy.
//!
//! - RRKM (SBB10 eq. 6): k(E) = sigma N(E - E0) / (h rho(E)), with N the sum of states of the transition state and
//!   rho the density of states of the ion, both counted on cells with harmonic vibrations, classical rotors and
//!   internal rotors (`masterequation::microcanonical_builder`). The symmetry numbers, chirality and electronic
//!   degeneracies of the two species enter rho and N; sigma is an additional reaction-path degeneracy (1 when the
//!   symmetry numbers already give it).
//! - Phase space theory and the simplified statistical adiabatic channel model (SBB10 eqs. 7-8) are options of the
//!   input that are refused with an error: they are not implemented yet.
//!
//! Reference: B. Sztáray, A. Bodi, T. Baer, J. Mass Spectrom. 45, 1233 (2010) (SBB10).

use crate::constants::H_PLANCK_CM;
use crate::masterequation::microcanonical_builder::{rrho_density_of_states, rrho_sum_of_states, SpeciesMicroModel};

/// Unimolecular rate theory of an ion dissociation channel (SBB10 pp. 1237-1238).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RateModel {
    /// Rigid activated complex RRKM theory (SBB10 eq. 6).
    Rrkm,
    /// Phase space theory: not available yet.
    PhaseSpaceTheory,
    /// Simplified statistical adiabatic channel model (SBB10 eqs. 7-8): not available yet.
    SimplifiedStatisticalAdiabaticChannel,
}

/// SBB10 eq. 6 on cells: k[i] = sigma N[i - onset] / (h rho[i]) for i >= onset, 0 below the onset and where the ion
/// has no states (s-1; rho per cm-1, N dimensionless). The sum of states must cover the cells above the onset.
pub fn rrkm_rate_constants(rho_ion: &[f64], sum_ts: &[f64], onset_cells: usize, reaction_degeneracy: f64) -> Result<Vec<f64>, String> {
    if sum_ts.len() + onset_cells < rho_ion.len() {
        return Err(format!(
            "RRKM rate constants: the sum of states covers {} cells above the onset at cell {onset_cells}, the density \
             of states {} cells.",
            sum_ts.len(),
            rho_ion.len()
        ));
    }
    Ok(rho_ion
        .iter()
        .enumerate()
        .map(|(i, rho)| {
            if i < onset_cells || !(*rho > 0.0) {
                0.0
            } else {
                reaction_degeneracy * sum_ts[i - onset_cells] / (H_PLANCK_CM * rho)
            }
        })
        .collect())
}

/// k(E) of the ion on `cells` cells of width `cell_cm1` above its ground state, for the dissociation limit (or barrier)
/// at cell `onset_cells`.
pub fn rate_constants(
    model: RateModel,
    ion: &SpeciesMicroModel,
    transition_state: &SpeciesMicroModel,
    onset_cells: usize,
    cells: usize,
    cell_cm1: f64,
    reaction_degeneracy: f64,
) -> Result<Vec<f64>, String> {
    match model {
        RateModel::Rrkm => {
            let rho = rrho_density_of_states(cells, cell_cm1, ion)?;
            let sum = rrho_sum_of_states(cells.saturating_sub(onset_cells).max(1), cell_cm1, transition_state)?;
            rrkm_rate_constants(&rho, &sum, onset_cells, reaction_degeneracy)
        }
        RateModel::PhaseSpaceTheory | RateModel::SimplifiedStatisticalAdiabaticChannel => Err(format!(
            "Rate model {model:?} for the ion '{}' is not available yet; use RRKM.",
            ion.name
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::H_PLANCK_CM;
    use crate::masterequation::microcanonical_builder::{rrho_density_of_states, rrho_sum_of_states, SpeciesMicroModel};

    fn species(name: &str, frequencies: Vec<f64>) -> SpeciesMicroModel {
        SpeciesMicroModel {
            name: name.into(),
            vibrational_frequencies_cm1: frequencies,
            rotational_constants_cm1: vec![0.9, 0.12, 0.11],
            symmetry_number: 1.0,
            chirality_number: 1.0,
            electronic_degeneracy: 2.0,
            internal_rotors: Vec::new(),
        }
    }

    #[test]
    fn rrkm_rate_constants_are_the_transition_state_sum_over_h_times_the_density() {
        // SBB10 eq. 6: k(E) = sigma N(E - E0) / (h rho(E)), zero below the onset
        let rho = vec![0.0, 2.0, 4.0, 5.0, 8.0];
        let sum = vec![1.0, 3.0, 7.0];
        let k = rrkm_rate_constants(&rho, &sum, 2, 2.0).unwrap();
        assert_eq!(&k[..2], &[0.0, 0.0]);
        for (i, w) in [(2, 1.0), (3, 3.0), (4, 7.0)] {
            assert!((k[i] - 2.0 * w / (H_PLANCK_CM * rho[i])).abs() < 1e-6 * k[i]);
        }
        assert!(rrkm_rate_constants(&rho, &sum[..1], 2, 1.0).is_err());
    }

    #[test]
    fn the_rrkm_model_counts_the_states_of_the_ion_and_the_transition_state() {
        let ion = species("ion", vec![300.0, 650.0, 900.0, 1200.0, 1500.0, 3000.0]);
        let ts = species("ts", vec![250.0, 700.0, 1100.0, 1400.0, 2900.0]);
        let (cells, onset) = (12_000, 7000);
        let k = rate_constants(RateModel::Rrkm, &ion, &ts, onset, cells, 1.0, 1.0).unwrap();
        let rho = rrho_density_of_states(cells, 1.0, &ion).unwrap();
        let w = rrho_sum_of_states(cells - onset, 1.0, &ts).unwrap();
        // classical rotors: the transition state has no states at exactly zero energy
        assert_eq!((w[0], k[onset]), (0.0, 0.0));
        for i in [7001, 9000, 11_999] {
            assert!((k[i] - w[i - onset] / (H_PLANCK_CM * rho[i])).abs() < 1e-12 * k[i], "cell {i}: {} vs {}", k[i], w[i - onset] / (H_PLANCK_CM * rho[i]));
        }
        assert!(k[..onset].iter().all(|x| *x == 0.0));
        assert!(k[11_999] > k[9000] && k[9000] > k[7001]);
    }

    #[test]
    fn phase_space_theory_and_the_simplified_adiabatic_channel_model_are_not_available_yet() {
        let ion = species("ion", vec![500.0]);
        for model in [RateModel::PhaseSpaceTheory, RateModel::SimplifiedStatisticalAdiabaticChannel] {
            let err = rate_constants(model, &ion, &ion, 10, 100, 1.0, 1.0).unwrap_err();
            assert!(err.contains("not available yet"), "{err}");
        }
    }
}
