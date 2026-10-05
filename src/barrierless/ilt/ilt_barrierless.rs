//! Microcanonical rate coefficients from high-pressure rate coefficients by inverse Laplace
//! transformation (ILT).
//!
//! The high-pressure dissociation rate coefficient is the Laplace transform of k(E) rho(E),
//!   k_inf(T) Q(T) = integral k(E) rho(E) exp(-E/kT) dE,
//! so k(E) rho(E) follows by inversion (Slater; Forst; Davies, Green, Pilling, Chem. Phys. Lett. 126,
//! 373 (1986), eq. 1). For a recombination B + C -> A the reverse rate coefficient is related through
//! the equilibrium constant, k_d = k_r (Q_P/Q_R) exp(-dE0/kT) with Q_P = C' (kT)^3/2 L{N_P}, which gives
//! DGP86 eq. 2 for k_r = A exp(-E/kT). "The Arrhenius form ... need not necessarily be employed - more
//! complex functional forms may be transformed analytically" (DGP86 p. 377); for
//!   k_inf(T) = A (T/T_ref)^n exp(-E_inf/kT) = A beta_ref^n beta^-n exp(-beta E_inf),  beta = 1/kT,
//! the Laplace pair L{x^(nu-1)/Gamma(nu)} = beta^-nu (nu > 0) gives (energies in cm-1):
//!
//!   dissociation, k_inf in s-1:
//!     k(E) rho(E) = A beta_ref^n / Gamma(n) integral_0^{E-E_inf} rho(E - E_inf - x) x^(n-1) dx,  n > 0,
//!     and k(E) rho(E) = A rho(E - E_inf) for n = 0 (Robertson, Comprehensive Chemical Kinetics 43 (2019),
//!     eq. 7.28);
//!   association, k_inf in cm3 s-1, fragments' rovibrational density N_P, reduced mass mu:
//!     k(E) rho(E) = A C'(mu) beta_ref^n / Gamma(n + 3/2)
//!                   integral_0^{E-E_th} N_P(E_P) (E - E_th - E_P)^(n+1/2) dE_P,   n > -3/2,
//!     E_th = dE0 + E_inf above the ground state of A (DGP86 eq. 2 for n = 0),
//!     C'(mu) = (2 pi mu / h^2)^(3/2), so that C'(kT)^(3/2) is the translational partition function per
//!     unit volume of the relative motion.
//! The Gamma-function arguments must be positive: n >= 0 (dissociation) and n > -3/2 (association).
//! E_inf must not be negative: otherwise the inversion "predicts that k(E) is non-zero significantly
//! below the threshold energy" (DGP86 p. 377).
//!
//! Numerics: on the grain grid, rho[m] holds the states in ((m-1) dE, m dE] per dE (rho[0] the ground
//! state, see `rrkm::sum_and_density`). The kernel is integrated analytically over each grain,
//!   w_j = integral_{j dE}^{(j+1) dE} x^(nu-1) dx / Gamma(nu) = [(j+1)^nu - j^nu] dE^nu / Gamma(nu+1),
//! which is finite for every nu > 0 (the integrable singularity of x^(nu-1) at x = 0 for nu < 1 is
//! never evaluated), and the integral becomes the discrete convolution sum_{j=0}^{i} w_j rho[i-j], exact
//! for the piecewise-constant density. For nu -> 0 the weights become w_0 = 1, w_j = 0 (the delta kernel).
//!
//! The result is returned as the RRKM-equivalent sum of states W(e) = h k(E) rho(E) at e = E - E_th =
//! i dE (h in cm-1 s), so that k(E) = W(E - E_th)/(h rho(E)) like for a transition state.

use crate::constants::{AMU_TO_KG, CLIGHT_SI, H_PLANCK_CM, KB_CM, PI, PLANCK_SI};
use crate::numeric::lanczos_gamma::gamma_func;

/// k_inf(T) = A (T/T_ref)^n exp(-E_inf/kT).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModifiedArrhenius {
    /// A: s-1 for a dissociation, cm3 s-1 for an association.
    pub pre_exponential: f64,
    /// n.
    pub temperature_exponent: f64,
    /// T_ref in K.
    pub reference_temperature_kelvin: f64,
    /// E_inf in cm-1.
    pub activation_energy_cm1: f64,
}

impl ModifiedArrhenius {
    pub fn rate(&self, temperature_kelvin: f64) -> f64 {
        self.pre_exponential
            * (temperature_kelvin / self.reference_temperature_kelvin).powf(self.temperature_exponent)
            * (-self.activation_energy_cm1 / (KB_CM * temperature_kelvin)).exp()
    }
}

/// C'(mu) = (2 pi mu / h^2)^(3/2) in cm-3 (cm-1)^(-3/2) for the reduced mass mu in amu:
/// C'(mu) (kT)^(3/2) is the translational partition function per cm3 of the relative motion.
pub fn translational_partition_constant(reduced_mass_amu: f64) -> f64 {
    // Energy unit 1 cm-1 = h c (c in cm/s) joule; (2 pi mu E / h^2)^(3/2) in m-3, times 1e-6 for cm-3.
    let joule_per_cm1 = PLANCK_SI * CLIGHT_SI * 100.0;
    (2.0 * PI * reduced_mass_amu * AMU_TO_KG * joule_per_cm1 / (PLANCK_SI * PLANCK_SI)).powf(1.5) * 1.0e-6
}

fn validate(k_inf: &ModifiedArrhenius, density: &[f64], grain_width_cm1: f64) -> Result<(), String> {
    if !(k_inf.pre_exponential > 0.0) || !k_inf.pre_exponential.is_finite() {
        return Err(format!("ILT: the pre-exponential factor must be positive (got {}).", k_inf.pre_exponential));
    }
    if !(k_inf.reference_temperature_kelvin > 0.0) {
        return Err("ILT: the reference temperature must be positive.".into());
    }
    if !(k_inf.activation_energy_cm1 >= 0.0) {
        return Err(format!(
            "ILT: negative activation energy {} cm-1; the inversion would give k(E) > 0 below the threshold \
             (Davies, Green, Pilling 1986, p. 377).",
            k_inf.activation_energy_cm1
        ));
    }
    if !(grain_width_cm1 > 0.0) {
        return Err("ILT: the grain width must be positive.".into());
    }
    if density.is_empty() || density.iter().any(|r| !(*r >= 0.0) || !r.is_finite()) {
        return Err("ILT: the density of states must be non-empty, finite and >= 0.".into());
    }
    Ok(())
}

/// sum_{j=0}^{i} w_j rho[i-j] with w_j = integral over [j dE, (j+1) dE) of x^(nu-1)/Gamma(nu) dx
/// = [(j+1)^nu - j^nu] dE^nu / Gamma(nu + 1); nu = 0 gives the delta kernel w_0 = 1.
fn grain_integrated_convolution(nu: f64, density: &[f64], grain_width_cm1: f64) -> Vec<f64> {
    let n = density.len();
    let scale = grain_width_cm1.powf(nu) / gamma_func(nu + 1.0);
    // (j+1)^nu - j^nu = j^nu expm1(nu ln(1 + 1/j)) avoids cancellation at large j.
    let w: Vec<f64> = (0..n)
        .map(|j| {
            if j == 0 {
                scale
            } else {
                let jf = j as f64;
                jf.powf(nu) * (nu * (1.0 / jf).ln_1p()).exp_m1() * scale
            }
        })
        .collect();
    (0..n).map(|i| (0..=i).map(|j| w[j] * density[i - j]).sum()).collect()
}

/// W(e) = h k(E) rho(E), e = E - E_inf = i dE, of a dissociation with high-pressure rate coefficient
/// `k_inf` (s-1); `rho_reactant` is the rovibrational density of the dissociating molecule (per cm-1).
pub fn ilt_sum_of_states_dissociation(
    k_inf: &ModifiedArrhenius,
    rho_reactant: &[f64],
    grain_width_cm1: f64,
) -> Result<Vec<f64>, String> {
    validate(k_inf, rho_reactant, grain_width_cm1)?;
    let nu = k_inf.temperature_exponent;
    if !(nu >= 0.0) {
        return Err(format!(
            "ILT of a dissociation: temperature exponent n = {nu} < 0 is outside the domain of the inversion \
             (Gamma(n) requires n >= 0)."
        ));
    }
    let beta_ref = 1.0 / (KB_CM * k_inf.reference_temperature_kelvin);
    let prefactor = H_PLANCK_CM * k_inf.pre_exponential * beta_ref.powf(nu);
    Ok(grain_integrated_convolution(nu, rho_reactant, grain_width_cm1).into_iter().map(|x| prefactor * x).collect())
}

/// W(e) = h k(E) rho(E), e = E - E_th = i dE with E_th = dE0 + E_inf, of the dissociation whose reverse
/// association has the high-pressure rate coefficient `k_inf` (cm3 s-1); `rho_fragments` is the
/// convolved rovibrational density of the two fragments (per cm-1), `reduced_mass_amu` their reduced mass.
pub fn ilt_sum_of_states_association(
    k_inf: &ModifiedArrhenius,
    rho_fragments: &[f64],
    reduced_mass_amu: f64,
    grain_width_cm1: f64,
) -> Result<Vec<f64>, String> {
    validate(k_inf, rho_fragments, grain_width_cm1)?;
    if !(reduced_mass_amu > 0.0) {
        return Err("ILT of an association: the reduced mass must be positive.".into());
    }
    let n = k_inf.temperature_exponent;
    let nu = n + 1.5;
    if !(nu > 0.0) {
        return Err(format!(
            "ILT of an association: temperature exponent n = {n} <= -3/2 is outside the domain of the \
             inversion (Gamma(n + 3/2) requires n > -3/2)."
        ));
    }
    let beta_ref = 1.0 / (KB_CM * k_inf.reference_temperature_kelvin);
    let prefactor = H_PLANCK_CM * k_inf.pre_exponential * translational_partition_constant(reduced_mass_amu) * beta_ref.powf(n);
    Ok(grain_integrated_convolution(nu, rho_fragments, grain_width_cm1).into_iter().map(|x| prefactor * x).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::BOLTZMANN_SI;

    /// Smooth model density (per cm-1) with a ground state, rho[0] = 1/dE.
    fn model_density(n: usize, d_e: f64) -> Vec<f64> {
        (0..n).map(|i| if i == 0 { 1.0 / d_e } else { 0.02 * (i as f64 * d_e).powf(1.5) / 1000.0 + 0.001 }).collect()
    }

    /// k_inf(T) of a dissociation from W: sum_i (W_i/h) exp(-(E_inf + i dE)/kT) dE / Q(T).
    fn forward_dissociation(w: &[f64], rho: &[f64], d_e: f64, e_inf: f64, t: f64) -> f64 {
        let kt = KB_CM * t;
        let flux: f64 = w.iter().enumerate().map(|(i, w)| w / H_PLANCK_CM * (-(e_inf + i as f64 * d_e) / kt).exp() * d_e).sum();
        let q: f64 = rho.iter().enumerate().map(|(i, r)| r * (-(i as f64) * d_e / kt).exp() * d_e).sum();
        flux / q
    }

    #[test]
    fn dissociation_with_n_one_gives_the_sum_of_states_of_the_reactant() {
        // nu = 1: k rho = A beta_ref integral rho = A beta_ref G(E - E_inf), G including the ground state.
        let d_e = 10.0;
        let rho = model_density(300, d_e);
        let k = ModifiedArrhenius { pre_exponential: 1e13, temperature_exponent: 1.0, reference_temperature_kelvin: 300.0, activation_energy_cm1: 2000.0 };
        let w = ilt_sum_of_states_dissociation(&k, &rho, d_e).unwrap();
        let beta_ref = 1.0 / (KB_CM * 300.0);
        let mut g = 0.0;
        for i in 0..rho.len() {
            g += rho[i] * d_e;
            let expected = H_PLANCK_CM * 1e13 * beta_ref * g;
            assert!(((w[i] - expected) / expected).abs() < 1e-12, "W[{i}] = {} vs {expected}", w[i]);
        }
    }

    #[test]
    fn dissociation_with_n_zero_is_the_shifted_density() {
        // Robertson (2019) eq. 7.28: k(E) = A rho(E - E_inf)/rho(E).
        let d_e = 10.0;
        let rho = model_density(300, d_e);
        let k = ModifiedArrhenius { pre_exponential: 3e12, temperature_exponent: 0.0, reference_temperature_kelvin: 298.0, activation_energy_cm1: 1500.0 };
        let w = ilt_sum_of_states_dissociation(&k, &rho, d_e).unwrap();
        for i in 0..rho.len() {
            let expected = H_PLANCK_CM * 3e12 * rho[i];
            assert!(((w[i] - expected) / expected).abs() < 1e-12);
        }
    }

    #[test]
    fn dissociation_inversion_reproduces_k_inf_for_a_fractional_exponent() {
        // n = 0.5: x^(n-1) is singular at x = 0.
        let d_e = 2.0;
        let rho = model_density(20000, d_e);
        let k = ModifiedArrhenius { pre_exponential: 2e13, temperature_exponent: 0.5, reference_temperature_kelvin: 300.0, activation_energy_cm1: 3000.0 };
        let w = ilt_sum_of_states_dissociation(&k, &rho, d_e).unwrap();
        for t in [300.0, 600.0, 1000.0] {
            let got = forward_dissociation(&w, &rho, d_e, 3000.0, t);
            let expected = k.rate(t);
            assert!(((got - expected) / expected).abs() < 1e-2, "T = {t}: {got:e} vs {expected:e}");
        }
    }

    #[test]
    fn association_inversion_reproduces_the_recombination_rate_coefficient() {
        // n = -1: (E - E_th - E_P)^(n+1/2) is singular at the upper limit. Detailed balance:
        //   k_rec(T) C'(mu) (kT)^(3/2) Q_frag(T) = sum_E k(E) rho(E) exp(-E/kT) dE  (E above the asymptote).
        let d_e = 2.0;
        let rho_p = model_density(20000, d_e);
        let mu = 15.0;
        let k = ModifiedArrhenius { pre_exponential: 5e-11, temperature_exponent: -1.0, reference_temperature_kelvin: 298.0, activation_energy_cm1: 100.0 };
        let w = ilt_sum_of_states_association(&k, &rho_p, mu, d_e).unwrap();
        for t in [300.0, 600.0, 1000.0] {
            let kt = KB_CM * t;
            let flux: f64 = w.iter().enumerate().map(|(i, w)| w / H_PLANCK_CM * (-(100.0 + i as f64 * d_e) / kt).exp() * d_e).sum();
            let q_frag: f64 = rho_p.iter().enumerate().map(|(i, r)| r * (-(i as f64) * d_e / kt).exp() * d_e).sum();
            let got = flux / (translational_partition_constant(mu) * kt.powf(1.5) * q_frag);
            let expected = k.rate(t);
            assert!(((got - expected) / expected).abs() < 1e-2, "T = {t}: {got:e} vs {expected:e}");
        }
    }

    #[test]
    fn translational_partition_constant_gives_the_partition_function_per_volume() {
        // (2 pi mu k_B T / h^2)^(3/2) in SI with k_B T in joule, converted to cm-3.
        let (mu, t) = (20.0, 400.0);
        let si = (2.0 * PI * mu * AMU_TO_KG * BOLTZMANN_SI * t / PLANCK_SI.powi(2)).powf(1.5) * 1e-6;
        let ours = translational_partition_constant(mu) * (KB_CM * t).powf(1.5);
        assert!(((ours - si) / si).abs() < 1e-6, "{ours:e} vs {si:e}");
    }

    #[test]
    fn exponents_outside_the_gamma_function_domain_and_negative_activation_energies_are_rejected() {
        let rho = model_density(100, 10.0);
        let base = ModifiedArrhenius { pre_exponential: 1e13, temperature_exponent: 0.0, reference_temperature_kelvin: 300.0, activation_energy_cm1: 0.0 };
        assert!(ilt_sum_of_states_dissociation(&ModifiedArrhenius { temperature_exponent: -0.1, ..base }, &rho, 10.0).is_err());
        assert!(ilt_sum_of_states_association(&ModifiedArrhenius { temperature_exponent: -1.5, ..base }, &rho, 10.0, 10.0).is_err());
        assert!(ilt_sum_of_states_association(&ModifiedArrhenius { temperature_exponent: -1.4, ..base }, &rho, 10.0, 10.0).is_ok());
        assert!(ilt_sum_of_states_dissociation(&ModifiedArrhenius { activation_energy_cm1: -50.0, ..base }, &rho, 10.0).is_err());
        assert!(ilt_sum_of_states_dissociation(&ModifiedArrhenius { pre_exponential: 0.0, ..base }, &rho, 10.0).is_err());
    }
}
