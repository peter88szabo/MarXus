//! Breakdown curves: fractional abundances of the parent and fragment ions against the photon energy.
//!
//! - Fast dissociation (SBB10 eq. 23): every ion above the dissociation limit dissociates,
//!     BD(h nu) = integral_0^{E0 - IE} P_i(E, h nu) dE,
//!   with P_i the normalized ion distribution at h nu (`energy_distributions::ion_energy_distribution`) and E0 the
//!   0 K appearance energy, the dissociation limit of the ion measured from the ground state of the neutral.
//! - Slow dissociation, kinetic shift (SBB10 eq. 24): ions above the limit that do not dissociate within the maximum
//!   flight time tau_max are counted as parent ions,
//!     BD(h nu) = integral_0^{E0 - IE} P_i dE + integral_{E0 - IE}^inf P_i exp(-k(E) tau_max) dE.
//! - Parallel channels within the flight time (SBB10 eq. 21 without isomerization) and sequential fast dissociations
//!   with the statistical product energy distribution of each step (SBB10 eq. 5, p. 1235 steps 4-6).
//! On cells, the ion cells i < onset = round((E0 - IE)/dE) lie below the limit.
//!
//! Reference: B. Sztáray, A. Bodi, T. Baer, J. Mass Spectrom. 45, 1233 (2010) (SBB10).

use super::energy_distributions::{ion_energy_distribution, EV_TO_CM1};
use super::product_energy::daughter_distribution;

/// Rate constants k(E) of the parent ion on its energy cells (s-1) and the maximum flight time (s) within which a
/// dissociation is observed (SBB10 eq. 24).
#[derive(Debug, Clone, Copy)]
pub struct KineticShift<'a> {
    pub rate_constants_s_inv: &'a [f64],
    pub flight_time_s: f64,
}

/// One point of a breakdown curve: fractional abundances of the parent and the fragment ion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BreakdownPoint {
    pub photon_energy_ev: f64,
    pub parent: f64,
    pub fragment: f64,
}

/// SBB10 eqs. 23-24: the fraction of parent ions of the normalized ion distribution `ion` (cells above the ion ground
/// state), with the dissociation limit at cell `onset_cells`.
pub fn parent_fraction(ion: &[f64], onset_cells: usize, kinetic_shift: Option<KineticShift>) -> Result<f64, String> {
    let onset = onset_cells.min(ion.len());
    let below: f64 = ion[..onset].iter().sum();
    let Some(shift) = kinetic_shift else {
        return Ok(below);
    };
    if shift.rate_constants_s_inv.len() < ion.len() {
        return Err(format!(
            "Breakdown curve: rate constants on {} cells, the ion distribution reaches {} cells.",
            shift.rate_constants_s_inv.len(),
            ion.len()
        ));
    }
    if !(shift.flight_time_s >= 0.0) {
        return Err(format!("Breakdown curve: flight time {} s must not be negative.", shift.flight_time_s));
    }
    let surviving: f64 = ion[onset..]
        .iter()
        .zip(&shift.rate_constants_s_inv[onset..])
        .map(|(p, k)| p * (-k * shift.flight_time_s).exp())
        .sum();
    Ok(below + surviving)
}

/// Breakdown curve of a fast dissociation (SBB10 eq. 23) from the normalized thermal distribution of the neutral
/// (`energy_distributions::neutral_thermal_distribution`), the adiabatic ionization energy IE and the 0 K appearance
/// energy E0 (eV), with an optional Gaussian energy resolution (FWHM, cm-1).
pub fn breakdown_curve(
    neutral: &[f64],
    cell_cm1: f64,
    ionization_energy_ev: f64,
    appearance_energy_ev: f64,
    photon_energies_ev: &[f64],
    resolution_fwhm_cm1: Option<f64>,
) -> Result<Vec<BreakdownPoint>, String> {
    if !(appearance_energy_ev >= ionization_energy_ev) {
        return Err(format!(
            "Breakdown curve: the appearance energy {appearance_energy_ev} eV lies below the ionization energy \
             {ionization_energy_ev} eV."
        ));
    }
    let onset = ((appearance_energy_ev - ionization_energy_ev) * EV_TO_CM1 / cell_cm1).round() as usize;
    photon_energies_ev
        .iter()
        .map(|&h_nu| {
            let ion = ion_energy_distribution(neutral, cell_cm1, (h_nu - ionization_energy_ev) * EV_TO_CM1, resolution_fwhm_cm1)?;
            let parent = parent_fraction(&ion, onset, None)?;
            Ok(BreakdownPoint { photon_energy_ev: h_nu, parent, fragment: 1.0 - parent })
        })
        .collect()
}

/// One point of a breakdown curve with parallel dissociation channels of the parent ion.
#[derive(Debug, Clone, PartialEq)]
pub struct ParallelBreakdownPoint {
    pub photon_energy_ev: f64,
    pub parent: f64,
    /// Fractional abundance of the fragment ion of each channel, in the order of the rate constants.
    pub fragments: Vec<f64>,
}

/// Breakdown curve of parallel dissociation channels j with rate constants k_j(E) on the cells of the parent ion
/// (zero below each channel's limit), observed within the flight time tau. Ions at E survive with exp(-k_tot tau) and
/// form fragment j with (k_j / k_tot)(1 - exp(-k_tot tau)), k_tot = sum_j k_j: SBB10 eq. 21 without isomerization
/// (k_1 = k_-1 = 0), integrated over the ion distribution as in eq. 24.
pub fn parallel_breakdown_curve(
    neutral: &[f64],
    cell_cm1: f64,
    ionization_energy_ev: f64,
    rate_constants_s_inv: &[&[f64]],
    flight_time_s: f64,
    photon_energies_ev: &[f64],
    resolution_fwhm_cm1: Option<f64>,
) -> Result<Vec<ParallelBreakdownPoint>, String> {
    if !(flight_time_s >= 0.0) {
        return Err(format!("Breakdown curve: flight time {flight_time_s} s must not be negative."));
    }
    photon_energies_ev
        .iter()
        .map(|&h_nu| {
            let ion = ion_energy_distribution(neutral, cell_cm1, (h_nu - ionization_energy_ev) * EV_TO_CM1, resolution_fwhm_cm1)?;
            if let Some(short) = rate_constants_s_inv.iter().find(|k| k.len() < ion.len()) {
                return Err(format!(
                    "Breakdown curve at {h_nu} eV: rate constants on {} cells, the ion distribution reaches {} cells.",
                    short.len(),
                    ion.len()
                ));
            }
            let mut parent = 0.0;
            let mut fragments = vec![0.0; rate_constants_s_inv.len()];
            for (i, p) in ion.iter().enumerate().filter(|(_, p)| **p != 0.0) {
                let total: f64 = rate_constants_s_inv.iter().map(|k| k[i]).sum();
                let survive = (-total * flight_time_s).exp();
                parent += p * survive;
                if total > 0.0 {
                    for (fragment, k) in fragments.iter_mut().zip(rate_constants_s_inv) {
                        *fragment += p * k[i] / total * (1.0 - survive);
                    }
                }
            }
            Ok(ParallelBreakdownPoint { photon_energy_ev: h_nu, parent, fragments })
        })
        .collect()
}

/// One step of a sequence of fast dissociations: the limit on the energy cells of the dissociating ion, and the
/// densities of states of its fragment ion, the neutral fragment and the relative translation (SBB10 eq. 5).
#[derive(Debug, Clone, Copy)]
pub struct SequentialStep<'a> {
    pub onset_cells: usize,
    pub rho_fragment: &'a [f64],
    pub rho_neutral: &'a [f64],
    pub rho_translation: &'a [f64],
}

/// Fractional abundances [parent, first fragment, ..., last fragment] of a sequence of fast dissociations from the
/// normalized ion distribution `ion` (SBB10 p. 1235, steps 4-6): the ions above a limit dissociate, and the
/// distribution of the fragment ion follows from eq. 5 summed over the dissociating ions
/// (`product_energy::daughter_distribution`).
pub fn sequential_breakdown(ion: &[f64], steps: &[SequentialStep]) -> Result<Vec<f64>, String> {
    let mut abundances = Vec::with_capacity(steps.len() + 1);
    let mut distribution = ion.to_vec();
    for step in steps {
        let onset = step.onset_cells.min(distribution.len());
        abundances.push(distribution[..onset].iter().sum());
        distribution = daughter_distribution(&distribution, step.onset_cells, step.rho_fragment, step.rho_neutral, step.rho_translation)?;
    }
    abundances.push(distribution.iter().sum());
    Ok(abundances)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::photoion::energy_distributions::{ion_energy_distribution, neutral_thermal_distribution, EV_TO_CM1};

    #[test]
    fn fast_dissociation_leaves_the_ions_below_the_dissociation_limit() {
        // SBB10 eq. 23: BD(h nu) = integral_0^{E0 - IE} P_i(E, h nu) dE
        let ion = vec![0.1, 0.2, 0.3, 0.25, 0.15];
        assert!((parent_fraction(&ion, 3, None).unwrap() - 0.6).abs() < 1e-15);
        assert_eq!(parent_fraction(&ion, 0, None).unwrap(), 0.0);
        assert!((parent_fraction(&ion, 10, None).unwrap() - 1.0).abs() < 1e-15);
    }

    #[test]
    fn slow_dissociation_adds_the_ions_that_survive_the_flight_time() {
        // SBB10 eq. 24: BD = integral_0^{E0-IE} P_i dE + integral_{E0-IE}^inf P_i exp(-k(E) tau_max) dE
        let ion = vec![0.1, 0.2, 0.3, 0.25, 0.15];
        let tau: f64 = 2.0e-5;
        let k = vec![0.0, 0.0, 0.0, 1.0e4, 1.0e5];
        let expected = 0.6 + 0.25 * (-1.0e4 * tau).exp() + 0.15 * (-1.0e5 * tau).exp();
        let slow = parent_fraction(&ion, 3, Some(KineticShift { rate_constants_s_inv: &k, flight_time_s: tau })).unwrap();
        assert!((slow - expected).abs() < 1e-15);
        let fast = vec![0.0, 0.0, 0.0, 1.0e12, 1.0e12];
        let limit = parent_fraction(&ion, 3, Some(KineticShift { rate_constants_s_inv: &fast, flight_time_s: tau })).unwrap();
        assert!((limit - 0.6).abs() < 1e-15);
        assert!(parent_fraction(&ion, 3, Some(KineticShift { rate_constants_s_inv: &k[..2], flight_time_s: tau })).is_err());
    }

    #[test]
    fn without_thermal_energy_the_parent_ion_vanishes_in_a_step_at_the_appearance_energy() {
        // Baer, Bodi, Sztáray, PEPICO (Encyclopedia of Spectroscopy and Spectrometry, 2017), p. 641: "in the absence
        // of thermal energy" the onset "would correspond to a sharp step in the ion abundances" at E0.
        let neutral = vec![1.0];
        let (ie, e0) = (10.0, 11.133);
        for (h_nu, parent) in [(11.0, 1.0), (11.132, 1.0), (11.134, 0.0), (11.5, 0.0)] {
            let curve = breakdown_curve(&neutral, 1.0, ie, e0, &[h_nu], None).unwrap();
            assert_eq!(curve[0].parent, parent, "h nu = {h_nu} eV");
            assert!((curve[0].parent + curve[0].fragment - 1.0).abs() < 1e-15);
        }
    }

    #[test]
    fn the_thermal_breakdown_curve_is_the_neutral_distribution_integrated_up_to_the_onset() {
        // Baer, Bodi, Sztáray 2017, eq. 14: parent ion(h nu) = integral_0^{E0 - h nu} P(E, T) dE for h nu < E0, when
        // the whole thermal distribution is ionized (h nu above IE).
        let rho: Vec<f64> = (0..30_000).map(|i| if i % 700 == 0 || i % 1100 == 0 { 1.0 } else { 0.0 }).collect();
        let neutral = neutral_thermal_distribution(&rho, 1.0, 600.0).unwrap();
        let (ie, e0) = (9.0, 10.0);
        let photon_energies = [9.7, 9.85, 9.95, 9.99];
        let curve = breakdown_curve(&neutral, 1.0, ie, e0, &photon_energies, None).unwrap();
        for (point, h_nu) in curve.iter().zip(photon_energies) {
            let below = ((e0 - h_nu) * EV_TO_CM1).round() as usize;
            let expected: f64 = neutral[..below].iter().sum();
            assert!((point.parent - expected).abs() < 1e-12, "h nu = {h_nu}");
            assert_eq!(point.photon_energy_ev, h_nu);
        }
        // the same parent fraction from the ion distribution directly
        let ion = ion_energy_distribution(&neutral, 1.0, (9.95 - ie) * EV_TO_CM1, None).unwrap();
        let onset = ((e0 - ie) * EV_TO_CM1).round() as usize;
        assert!((parent_fraction(&ion, onset, None).unwrap() - curve[2].parent).abs() < 1e-12);
    }

    fn thermal_neutral() -> Vec<f64> {
        let rho: Vec<f64> = (0..30_000).map(|i| if i % 700 == 0 || i % 1100 == 0 { 1.0 } else { 0.0 }).collect();
        neutral_thermal_distribution(&rho, 1.0, 600.0).unwrap()
    }

    /// k(E) = a (E - onset)^2 above the onset, on 40 000 ion cells
    fn rates(onset: usize, a: f64) -> Vec<f64> {
        (0..40_000).map(|i| if i < onset { 0.0 } else { a * ((i - onset) as f64).powi(2) }).collect()
    }

    #[test]
    fn parallel_channels_compete_within_the_flight_time() {
        // Ions at E survive with exp(-k_tot tau) and form fragment j with (k_j / k_tot)(1 - exp(-k_tot tau)); one
        // channel reduces to SBB10 eq. 24.
        let neutral = thermal_neutral();
        let (ie, tau) = (9.0, 2.0e-5);
        let onset = (1.0 * EV_TO_CM1).round() as usize;
        let k1 = rates(onset, 0.5);
        let photon_energies = [9.9, 10.0, 10.2];
        let one = parallel_breakdown_curve(&neutral, 1.0, ie, &[&k1], tau, &photon_energies, None).unwrap();
        for (point, h_nu) in one.iter().zip(photon_energies) {
            let ion = ion_energy_distribution(&neutral, 1.0, (h_nu - ie) * EV_TO_CM1, None).unwrap();
            let expected = parent_fraction(&ion, onset, Some(KineticShift { rate_constants_s_inv: &k1, flight_time_s: tau })).unwrap();
            assert!((point.parent - expected).abs() < 1e-12, "h nu = {h_nu}");
            assert!((point.parent + point.fragments[0] - 1.0).abs() < 1e-12);
        }
        let k2: Vec<f64> = k1.iter().map(|k| 3.0 * k).collect();
        let two = parallel_breakdown_curve(&neutral, 1.0, ie, &[&k1, &k2], tau, &photon_energies, None).unwrap();
        for point in &two {
            assert!((point.fragments[1] - 3.0 * point.fragments[0]).abs() < 1e-12);
            assert!((point.parent + point.fragments.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        }
        let none = vec![0.0; 40_000];
        let still = parallel_breakdown_curve(&neutral, 1.0, ie, &[&none], tau, &[10.2], None).unwrap();
        assert_eq!((still[0].parent, still[0].fragments[0]), (1.0, 0.0));
        assert!(parallel_breakdown_curve(&neutral, 1.0, ie, &[&k1[..100]], tau, &[10.2], None).is_err());
    }

    #[test]
    fn sequential_fast_dissociations_pass_the_daughter_distribution_to_the_next_step() {
        // SBB10 p. 1235, steps 4-6: the internal energy distribution of the daughter ion from eq. 5, summed over the
        // parent distribution above the first limit; the daughter dissociates when above its own limit.
        use crate::photoion::product_energy::{daughter_distribution, translational_density};
        let ion: Vec<f64> = {
            let raw: Vec<f64> = (0..6000).map(|i| ((i as f64 - 3500.0) / 600.0).powi(2)).map(|x| (-x).exp()).collect();
            let total: f64 = raw.iter().sum();
            raw.iter().map(|x| x / total).collect()
        };
        let rho_f: Vec<f64> = (0..6000).map(|i| (i as f64 + 1.0).powf(4.0)).collect();
        let rho_n: Vec<f64> = (0..6000).map(|i| if i % 2100 == 0 { 1.0 } else { 0.0 }).collect();
        let rho_tr = translational_density(6000, 1.0, 2).unwrap();
        let steps = [
            SequentialStep { onset_cells: 2500, rho_fragment: &rho_f, rho_neutral: &rho_n, rho_translation: &rho_tr },
            SequentialStep { onset_cells: 400, rho_fragment: &rho_f, rho_neutral: &rho_n, rho_translation: &rho_tr },
        ];
        let abundances = sequential_breakdown(&ion, &steps).unwrap();
        assert_eq!(abundances.len(), 3);
        assert!((abundances.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        let daughter = daughter_distribution(&ion, 2500, &rho_f, &rho_n, &rho_tr).unwrap();
        let parent: f64 = ion[..2500].iter().sum();
        let daughter_below: f64 = daughter[..400].iter().sum();
        assert!((abundances[0] - parent).abs() < 1e-14);
        assert!((abundances[1] - daughter_below).abs() < 1e-12);
        assert!((abundances[2] - (1.0 - parent - daughter_below)).abs() < 1e-12);
        // a second limit above every daughter energy leaves all daughters intact
        let high = [steps[0], SequentialStep { onset_cells: 10_000, ..steps[1] }];
        let kept = sequential_breakdown(&ion, &high).unwrap();
        assert!((kept[1] - (1.0 - parent)).abs() < 1e-12 && kept[2] == 0.0);
    }
}
