use crate::constants::{ATM_TO_PASCAL, BOLTZMANN_SI, PI};
use super::reaction_network::CollisionModelParams;

/// Computes alpha(T) in cm^-1: alpha(T)=alpha_1000*(T/1000)^exp
pub(crate) fn compute_alpha_cm1(
    alpha_at_1000k: f64,
    alpha_exponent: f64,
    temperature_kelvin: f64,
) -> f64 {
    alpha_at_1000k * (temperature_kelvin / 1000.0).powf(alpha_exponent)
}

/// Lennard-Jones collision frequency Z(T,p) = pi sigma^2 <v_rel> n Omega22*(T*) in s^-1, with the
/// ideal-gas bath-gas number density n = p/(k_B T) and the mean relative speed
/// <v_rel> = sqrt(8 k_B T/(pi mu)), both from physical constants.
///
/// Units expectation:
/// - sigma in Å
/// - reduced mass in amu
/// - epsilon in K
/// - pressure in Torr
pub(crate) fn compute_collision_frequency_s_inv(
    params: &CollisionModelParams,
    temperature_kelvin: f64,
    pressure_torr: f64,
) -> Result<f64, String> {
    if params.lennard_jones_sigma_angstrom <= 0.0
        || params.reduced_mass_amu <= 0.0
        || params.lennard_jones_epsilon_kelvin <= 0.0
    {
        return Err(
            "Invalid collision parameters (sigma, epsilon, reduced mass must be positive)".into(),
        );
    }

    // cross section in cm^2: π σ^2 * 1e-16 (since 1 Å = 1e-8 cm)
    let cross_section_cm2 = PI * params.lennard_jones_sigma_angstrom.powi(2) * 1e-16;

    // mean relative speed in cm/s: sqrt(8 k_B T / (pi mu))
    const AMU_TO_KG: f64 = 1.660_539_066_60e-27;
    let reduced_mass_kg = params.reduced_mass_amu * AMU_TO_KG;
    let mean_speed_cm_s =
        (8.0 * BOLTZMANN_SI * temperature_kelvin / (PI * reduced_mass_kg)).sqrt() * 100.0;

    // bath number density in molecule/cm^3: n = p / (k_B T), 1 Torr = 101325/760 Pa
    let pressure_pa = pressure_torr * ATM_TO_PASCAL / 760.0;
    let bath_number_density = pressure_pa / (BOLTZMANN_SI * temperature_kelvin) * 1.0e-6;

    // collision integral Ω_22(T) (a simple log form)
    // Reduced collision integral, Troe, J. Chem. Phys. 66, 4758 (1977), eq. 3.3:
    //   Omega22* = [0.636 + 0.567 log10(kT/eps)]^-1, accurate to +-7% for 0.3 <= kT/eps <= 500.
    let t_over_eps = temperature_kelvin / params.lennard_jones_epsilon_kelvin;
    if !(0.3..=500.0).contains(&t_over_eps) {
        return Err(format!(
            "kT/eps = {t_over_eps:.4} is outside the validity range 0.3..500 of the collision \
             integral approximation (Troe 1977, eq. 3.3)."
        ));
    }
    let omega_22 = (0.636 + 0.567 * t_over_eps.log10()).recip();

    Ok(cross_section_cm2 * mean_speed_cm_s * bath_number_density * omega_22)
}

/// Returns [min, max_exclusive] indices for banded transitions.
pub(crate) fn band_limits(
    center: usize,
    min_allowed: usize,
    max_exclusive: usize,
    band: usize,
) -> (usize, usize) {
    let min_idx = center.saturating_sub(band).max(min_allowed);
    let max_idx = (center + band + 1).min(max_exclusive);
    (min_idx, max_idx)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::BOLTZMANN_SI;

    #[test]
    fn collision_frequency_uses_ideal_gas_density_and_mean_relative_speed() {
        // Z = pi sigma^2 <v_rel> n Omega22*,  n = p/(k_B T),  <v_rel> = sqrt(8 k_B T/(pi mu)),
        // evaluated here independently in SI units (1 Torr = 101325/760 Pa).
        let params = CollisionModelParams {
            lennard_jones_sigma_angstrom: 3.7,
            lennard_jones_epsilon_kelvin: 95.0,
            reduced_mass_amu: 14.0,
            alpha_at_1000K_cm1: 200.0,
            alpha_temperature_exponent: 0.85,
        };
        let (t, p_torr) = (300.0, 1.0);
        let z = compute_collision_frequency_s_inv(&params, t, p_torr).unwrap();

        let amu_kg = 1.660_539_066_60e-27;
        let torr_pa = 101_325.0 / 760.0;
        let n_cm3 = p_torr * torr_pa / (BOLTZMANN_SI * t) * 1.0e-6;
        let v_cm_s = (8.0 * BOLTZMANN_SI * t / (PI * 14.0 * amu_kg)).sqrt() * 100.0;
        let omega = (0.636 + 0.567 * (t / 95.0).log10()).recip();
        let z_ref = PI * (3.7e-8_f64).powi(2) * v_cm_s * n_cm3 * omega;
        assert!(((z - z_ref) / z_ref).abs() < 1e-9, "Z = {z:e} s-1, expected {z_ref:e} s-1");
    }

    #[test]
    fn collision_integral_outside_its_validity_range_is_an_error() {
        // Omega22* = [0.636 + 0.567 log10(kT/eps)]^-1 holds for 0.3 <= kT/eps <= 500
        // (Troe, J. Chem. Phys. 66, 4758 (1977), eq. 3.3); kT/eps = 0.1 is outside.
        let params = CollisionModelParams {
            lennard_jones_sigma_angstrom: 3.7,
            lennard_jones_epsilon_kelvin: 3000.0,
            reduced_mass_amu: 14.0,
            alpha_at_1000K_cm1: 200.0,
            alpha_temperature_exponent: 0.85,
        };
        let z = compute_collision_frequency_s_inv(&params, 300.0, 1.0);
        assert!(z.is_err(), "kT/eps = 0.1 must be rejected, got {z:?}");
    }
}
