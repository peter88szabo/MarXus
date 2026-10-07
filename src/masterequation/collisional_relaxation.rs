use crate::constants::{AMU_TO_KG, ATM_TO_PASCAL, BOLTZMANN_SI, PI};

/// Reduced collision integral Omega(2,2)*(T*) of the Lennard-Jones potential, T* = k_B T/epsilon.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum CollisionIntegral {
    /// Troe, J. Chem. Phys. 66, 4758 (1977), eq. 3.3: [0.636 + 0.567 log10 T*]^-1, accurate to +-7% for
    /// 0.3 <= T* <= 500.
    Troe1977,
    /// Neufeld, Janzen, Aziz, J. Chem. Phys. 57, 1100 (1972): A/T*^B + C exp(-D T*) + E exp(-F T*), the fit to the
    /// tabulated Lennard-Jones integrals for 0.3 <= T* <= 100 (the form used by MESS and MESMER). The default.
    #[default]
    Neufeld1972,
}

/// Omega(2,2)*(T*) of the selected form; an error outside its range of validity.
pub fn reduced_collision_integral(t_star: f64, integral: CollisionIntegral) -> Result<f64, String> {
    match integral {
        CollisionIntegral::Troe1977 => {
            if !(0.3..=500.0).contains(&t_star) {
                return Err(format!(
                    "kT/eps = {t_star:.4} is outside the validity range 0.3..500 of the collision integral approximation \
                     (Troe 1977, eq. 3.3)."
                ));
            }
            Ok((0.636 + 0.567 * t_star.log10()).recip())
        }
        CollisionIntegral::Neufeld1972 => {
            if !(0.3..=100.0).contains(&t_star) {
                return Err(format!(
                    "kT/eps = {t_star:.4} is outside the range 0.3..100 of the collision integral fit (Neufeld, Janzen, \
                     Aziz 1972)."
                ));
            }
            Ok(1.16145 / t_star.powf(0.14874) + 0.52487 * (-0.77320 * t_star).exp() + 2.16178 * (-2.43787 * t_star).exp())
        }
    }
}

/// Lennard-Jones collision frequency Z(T,p) = pi sigma^2 <v_rel> n Omega22*(T*) in s^-1, with the
/// ideal-gas bath-gas number density n = p/(k_B T) and the mean relative speed
/// <v_rel> = sqrt(8 k_B T/(pi mu)), both from physical constants (Troe, J. Chem. Phys. 66, 4758 (1977),
/// eqs. 3.1-3.3), for the combined pair parameters sigma (Å) and epsilon (K), the reduced mass (amu)
/// and the pressure (Torr); Omega22* of the selected form (`CollisionIntegral`).
pub(crate) fn lennard_jones_collision_frequency_s_inv(
    sigma_angstrom: f64,
    epsilon_kelvin: f64,
    reduced_mass_amu: f64,
    temperature_kelvin: f64,
    pressure_torr: f64,
    integral: CollisionIntegral,
) -> Result<f64, String> {
    if sigma_angstrom <= 0.0 || reduced_mass_amu <= 0.0 || epsilon_kelvin <= 0.0 {
        return Err(
            "Invalid collision parameters (sigma, epsilon, reduced mass must be positive)".into(),
        );
    }

    // cross section in cm^2: π σ^2 * 1e-16 (since 1 Å = 1e-8 cm)
    let cross_section_cm2 = PI * sigma_angstrom.powi(2) * 1e-16;

    // mean relative speed in cm/s: sqrt(8 k_B T / (pi mu))
    let reduced_mass_kg = reduced_mass_amu * AMU_TO_KG;
    let mean_speed_cm_s =
        (8.0 * BOLTZMANN_SI * temperature_kelvin / (PI * reduced_mass_kg)).sqrt() * 100.0;

    // bath number density in molecule/cm^3: n = p / (k_B T), 1 Torr = 101325/760 Pa
    let pressure_pa = pressure_torr * ATM_TO_PASCAL / 760.0;
    let bath_number_density = pressure_pa / (BOLTZMANN_SI * temperature_kelvin) * 1.0e-6;

    let omega_22 = reduced_collision_integral(temperature_kelvin / epsilon_kelvin, integral)?;

    Ok(cross_section_cm2 * mean_speed_cm_s * bath_number_density * omega_22)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::BOLTZMANN_SI;

    #[test]
    fn collision_frequency_uses_ideal_gas_density_and_mean_relative_speed() {
        // Z = pi sigma^2 <v_rel> n Omega22*,  n = p/(k_B T),  <v_rel> = sqrt(8 k_B T/(pi mu)),
        // evaluated here independently in SI units (1 Torr = 101325/760 Pa).
        let (t, p_torr) = (300.0, 1.0);
        let z = lennard_jones_collision_frequency_s_inv(3.7, 95.0, 14.0, t, p_torr, CollisionIntegral::Troe1977).unwrap();

        let amu_kg = 1.660_539_066_60e-27;
        let torr_pa = 101_325.0 / 760.0;
        let n_cm3 = p_torr * torr_pa / (BOLTZMANN_SI * t) * 1.0e-6;
        let v_cm_s = (8.0 * BOLTZMANN_SI * t / (PI * 14.0 * amu_kg)).sqrt() * 100.0;
        let omega = (0.636 + 0.567 * (t / 95.0).log10()).recip();
        let z_ref = PI * (3.7e-8_f64).powi(2) * v_cm_s * n_cm3 * omega;
        assert!(((z - z_ref) / z_ref).abs() < 1e-9, "Z = {z:e} s-1, expected {z_ref:e} s-1");
    }

    #[test]
    fn neufeld_collision_integral_is_the_three_term_fit() {
        // Neufeld, Janzen, Aziz, J. Chem. Phys. 57, 1100 (1972): Omega(2,2)* = A/T*^B + C exp(-D T*) + E exp(-F T*),
        // A = 1.16145, B = 0.14874, C = 0.52487, D = 0.77320, E = 2.16178, F = 2.43787, for 0.3 <= T* <= 100.
        for t_star in [0.3_f64, 1.0, 3.0, 10.0, 100.0] {
            let fit = 1.16145 / t_star.powf(0.14874) + 0.52487 * (-0.77320 * t_star).exp() + 2.16178 * (-2.43787 * t_star).exp();
            let omega = reduced_collision_integral(t_star, CollisionIntegral::Neufeld1972).unwrap();
            assert!((omega / fit - 1.0).abs() < 1e-14, "T* = {t_star}: {omega} vs {fit}");
        }
        assert!(reduced_collision_integral(0.2, CollisionIntegral::Neufeld1972).is_err());
        assert!(reduced_collision_integral(150.0, CollisionIntegral::Neufeld1972).is_err());
    }

    #[test]
    fn the_collision_frequency_is_proportional_to_the_selected_collision_integral() {
        let (t, p_torr) = (300.0, 760.0);
        let z_troe = lennard_jones_collision_frequency_s_inv(3.7, 95.0, 14.0, t, p_torr, CollisionIntegral::Troe1977).unwrap();
        let z_neufeld = lennard_jones_collision_frequency_s_inv(3.7, 95.0, 14.0, t, p_torr, CollisionIntegral::Neufeld1972).unwrap();
        let t_star = t / 95.0;
        let ratio = reduced_collision_integral(t_star, CollisionIntegral::Neufeld1972).unwrap()
            / reduced_collision_integral(t_star, CollisionIntegral::Troe1977).unwrap();
        assert!((z_neufeld / z_troe / ratio - 1.0).abs() < 1e-14);
        // The default is the Neufeld fit (the accurate form; Troe 1977 eq. 3.3 is a +-7% approximation).
        assert_eq!(CollisionIntegral::default(), CollisionIntegral::Neufeld1972);
    }

    #[test]
    fn collision_integral_outside_its_validity_range_is_an_error() {
        // Omega22* = [0.636 + 0.567 log10(kT/eps)]^-1 holds for 0.3 <= kT/eps <= 500
        // (Troe, J. Chem. Phys. 66, 4758 (1977), eq. 3.3); kT/eps = 0.1 is outside.
        let z = lennard_jones_collision_frequency_s_inv(3.7, 3000.0, 14.0, 300.0, 1.0, CollisionIntegral::Troe1977);
        assert!(z.is_err(), "kT/eps = 0.1 must be rejected, got {z:?}");
    }
}
