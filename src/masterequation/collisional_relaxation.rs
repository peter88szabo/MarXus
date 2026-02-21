use std::f64::consts::PI;

use super::reaction_network::{
    CollisionModelParams, MasterEquationSettings,
};

/// Computes alpha(T) in cm^-1: alpha(T)=alpha_1000*(T/1000)^exp
pub(crate) fn compute_alpha_cm1(
    alpha_at_1000k: f64,
    alpha_exponent: f64,
    temperature_kelvin: f64,
) -> f64 {
    alpha_at_1000k * (temperature_kelvin / 1000.0).powf(alpha_exponent)
}

/// Collision frequency z(T,p) in s^-1 using a SSUMES-like structure.
///
/// This function is intentionally parameterized to avoid hard-coding your constants.
/// You can wire in your exact prefactors here.
///
/// Units expectation:
/// - sigma in Å
/// - reduced mass in amu
/// - epsilon in K
/// - pressure in Torr
///
/// You must set `bathgas_number_density_prefactor` and `mean_speed_prefactor` in settings
/// to match your unit system.
pub(crate) fn compute_collision_frequency_s_inv(
    params: &CollisionModelParams,
    settings: &MasterEquationSettings,
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

    // mean speed in cm/s: prefactor * sqrt(T / μ)
    let mean_speed_cm_s =
        settings.mean_speed_prefactor * (temperature_kelvin / params.reduced_mass_amu).sqrt();

    // bath number density in molecule/cm^3: p / (prefactor * T)
    let bath_number_density =
        pressure_torr / (settings.bathgas_number_density_prefactor * temperature_kelvin);

    // collision integral Ω_22(T) (a simple log form)
    let t_over_eps = temperature_kelvin / params.lennard_jones_epsilon_kelvin;
    if t_over_eps <= 0.0 {
        return Err("Invalid T/epsilon".into());
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
