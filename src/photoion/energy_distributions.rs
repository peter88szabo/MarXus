//! Internal energy distributions of a threshold photoionization experiment.
//!
//! - Neutral precursor at the sample temperature (SBB10 eq. 1): P(E) = rho(E) exp(-E/kT) / sum_E rho exp(-E/kT),
//!   with rho the rovibrational density of states (classical rotors, internal rotors and harmonic vibrations
//!   counted on cells, `masterequation::microcanonical_builder`).
//! - Molecular ion at photon energy h nu (SBB10 eq. 2): the neutral distribution is "transposed onto the ion
//!   manifold upon threshold ionization", E_ion = E_neutral + h nu - IE_ad, and convolved with the photon and
//!   electron-analyzer resolution functions (SBB10 p. 1236; approximated here by a normalized Gaussian of given
//!   FWHM). Ions below their ground state are not formed; the ion distribution is normalized.
//!
//! Distributions are probabilities per cell of width dE: cell i holds the energy i dE (rho as in
//! `rrkm::sum_and_density`). Energy shifts are rounded to whole cells.
//!
//! Reference: B. Sztáray, A. Bodi, T. Baer, J. Mass Spectrom. 45, 1233 (2010) (SBB10).

use crate::masterequation::chemical_activation_sources::thermal_distribution;

/// 1 eV / (h c) in cm-1 (CODATA 2018, exact from the defined h, c and e).
pub const EV_TO_CM1: f64 = 8065.543_937_349;

/// SBB10 eq. 1: the normalized thermal distribution of the neutral precursor from its density of states on cells.
pub fn neutral_thermal_distribution(rho_cells: &[f64], cell_cm1: f64, temperature_kelvin: f64) -> Result<Vec<f64>, String> {
    if !(temperature_kelvin > 0.0) || !(cell_cm1 > 0.0) {
        return Err(format!(
            "Neutral energy distribution: temperature ({temperature_kelvin} K) and cell width ({cell_cm1} cm-1) must be positive."
        ));
    }
    thermal_distribution(rho_cells, cell_cm1, crate::constants::KB_CM * temperature_kelvin)
}

/// SBB10 eq. 2: the normalized internal energy distribution of the molecular ion, E_ion = E_neutral + shift with
/// shift = h nu - IE_ad (cm-1), convolved with a normalized Gaussian resolution function of the given FWHM (cm-1).
/// The result runs from the ion ground state to the highest energy reached.
pub fn ion_energy_distribution(
    neutral: &[f64],
    cell_cm1: f64,
    shift_cm1: f64,
    resolution_fwhm_cm1: Option<f64>,
) -> Result<Vec<f64>, String> {
    let kernel = match resolution_fwhm_cm1 {
        None => vec![1.0],
        Some(fwhm) if fwhm > 0.0 => gaussian_kernel(fwhm / cell_cm1),
        Some(fwhm) => return Err(format!("Ion energy distribution: resolution FWHM {fwhm} cm-1 must be positive.")),
    };
    let half = (kernel.len() / 2) as isize;
    let shift = (shift_cm1 / cell_cm1).round() as isize;
    let len = neutral.len() as isize + shift + half;
    if len <= 0 {
        return Err(no_ions(shift_cm1));
    }
    let mut ion = vec![0.0; len as usize];
    for (i, p) in neutral.iter().enumerate().filter(|(_, p)| **p != 0.0) {
        for (k, g) in kernel.iter().enumerate() {
            let target = i as isize + shift + k as isize - half;
            if target >= 0 {
                ion[target as usize] += p * g;
            }
        }
    }
    let total: f64 = ion.iter().sum();
    if !(total > 0.0) {
        return Err(no_ions(shift_cm1));
    }
    ion.iter_mut().for_each(|x| *x /= total);
    Ok(ion)
}

fn no_ions(shift_cm1: f64) -> String {
    format!(
        "Ion energy distribution: no part of the neutral distribution lies above the ionization energy (h nu - IE = \
         {shift_cm1} cm-1)."
    )
}

/// Normalized Gaussian on cells, standard deviation FWHM / (2 sqrt(2 ln 2)), truncated at 6 standard deviations.
fn gaussian_kernel(fwhm_cells: f64) -> Vec<f64> {
    let sigma = fwhm_cells / (2.0 * (2.0 * std::f64::consts::LN_2).sqrt());
    let half = (6.0 * sigma).ceil() as isize;
    let mut kernel: Vec<f64> = (-half..=half).map(|d| (-0.5 * (d as f64 / sigma).powi(2)).exp()).collect();
    let total: f64 = kernel.iter().sum();
    kernel.iter_mut().for_each(|x| *x /= total);
    kernel
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::KB_CM;

    /// One harmonic oscillator of 1000 cm-1 on 1 cm-1 cells: one state at every 1000 cm-1.
    fn oscillator_density(cells: usize) -> Vec<f64> {
        (0..cells).map(|i| if i % 1000 == 0 { 1.0 } else { 0.0 }).collect()
    }

    #[test]
    fn the_neutral_distribution_is_the_boltzmann_weighted_density_of_states() {
        // SBB10 eq. 1: P(E) = rho(E) exp(-E/kT) / sum
        let t = 1000.0;
        let p = neutral_thermal_distribution(&oscillator_density(20_000), 1.0, t).unwrap();
        assert!((p.iter().sum::<f64>() - 1.0).abs() < 1e-14);
        assert!((p[1000] / p[0] - (-1000.0 / (KB_CM * t)).exp()).abs() < 1e-14);
        assert_eq!(p[500], 0.0);
    }

    #[test]
    fn the_ion_distribution_is_the_neutral_one_shifted_by_the_photon_energy_above_the_ionization_energy() {
        // SBB10 eq. 2: E_ion = E_neutral + h nu - IE_ad; ions below their ground state are not formed, and the
        // distribution is renormalized.
        let neutral = neutral_thermal_distribution(&oscillator_density(20_000), 1.0, 2000.0).unwrap();
        let above = ion_energy_distribution(&neutral, 1.0, 500.0, None).unwrap();
        assert_eq!(above.len(), neutral.len() + 500);
        assert!((above[1500] - neutral[1000]).abs() < 1e-15 && above[499] == 0.0);
        let below = ion_energy_distribution(&neutral, 1.0, -1500.0, None).unwrap();
        let kept: f64 = neutral[1500..].iter().sum();
        assert!((below[500] - neutral[2000] / kept).abs() < 1e-15);
        assert!((below.iter().sum::<f64>() - 1.0).abs() < 1e-14);
        assert!(ion_energy_distribution(&neutral, 1.0, -30_000.0, None).is_err());
    }

    #[test]
    fn the_energy_resolution_broadens_the_ion_distribution_by_a_normalized_gaussian() {
        // Convolution with the photon and electron-analyzer resolution (SBB10 p. 1236), approximated by a Gaussian of
        // the given FWHM: the total and the mean are kept, the variance grows by sigma^2 = (FWHM / (2 sqrt(2 ln 2)))^2.
        let mut neutral = vec![0.0; 4000];
        neutral[1000] = 0.5;
        neutral[1300] = 0.5;
        let fwhm = 120.0;
        let sharp = ion_energy_distribution(&neutral, 1.0, 1000.0, None).unwrap();
        let broad = ion_energy_distribution(&neutral, 1.0, 1000.0, Some(fwhm)).unwrap();
        let moment = |p: &[f64], k: i32| p.iter().enumerate().map(|(i, x)| x * (i as f64).powi(k)).sum::<f64>();
        assert!((moment(&broad, 0) - 1.0).abs() < 1e-12);
        assert!((moment(&broad, 1) - moment(&sharp, 1)).abs() < 1e-8);
        let variance = |p: &[f64]| moment(p, 2) - moment(p, 1).powi(2);
        let sigma = fwhm / (2.0 * (2.0 * std::f64::consts::LN_2).sqrt());
        assert!((variance(&broad) - variance(&sharp) - sigma * sigma).abs() < 1e-3 * sigma * sigma);
    }

    #[test]
    fn photon_energies_in_ev_are_converted_with_the_codata_factor() {
        // CODATA 2018: 1 eV / (h c) = 8065.543 937 cm-1
        assert!((EV_TO_CM1 - 8065.543_937).abs() < 1e-6);
    }
}
