//! MESS-like Phase Space Theory (PST) capture / loose-TS number-of-states model.
//!
//! Purpose:
//! - Provide a compact analytical model for the transition-state number of states N(E)
//!   for barrierless / capture-like association processes.
//! - This is useful when the "TS" is not a tight saddle-point but a loose, long-range
//!   bottleneck dominated by intermolecular motion and long-range attraction.
//!
//! Model summary (as implemented in MESS, mirrored here):
//! - Long-range potential: V(R) = V0 / R^n
//! - Two fragments treated as rigid bodies (linear or nonlinear) + relative motion.
//! - The resulting *cumulative number of states* scales as:
//!     N(E) = states_prefactor * E^power
//!   and the canonical partition-like weight scales as:
//!     Q(T) = weight_prefactor * T^power
//!
//! Notes:
//! - This module works in atomic units internally for consistency with MESS.
//! - Energies passed in cm^-1 are converted to Hartree.

use crate::constants::{AMU_TO_ELECTRON_MASS, CM1_TO_HARTREE, PI};
use crate::numeric::lanczos_gamma::gamma_func;
use crate::utils::atomic_masses::mass_vector_from_symbols_amu;

use super::types::{CaptureFragment, CaptureFragmentRotorModel, PhaseSpaceTheoryInput, PstTstLevel};

/// Resulting analytical PST state-count model:
///   N(E) = states_prefactor * E^power
///   Q(T) = weight_prefactor * T^power
#[derive(Clone, Debug)]
pub struct PhaseSpaceTheoryModel {
    /// Exponent on energy/temperature.
    pub power: f64,

    /// Prefactor for cumulative number of states N(E) in atomic units:
    /// E must be provided in Hartree.
    pub states_prefactor: f64,

    /// Prefactor for canonical weight Q(T) in atomic units:
    /// T must be provided as k_B T in Hartree (i.e. energy units).
    pub weight_prefactor: f64,
}

impl PhaseSpaceTheoryModel {
    /// Build the analytical PST model from fragment + potential parameters.
    ///
    /// This mirrors MESS' `Model::PhaseSpaceTheory` construction logic.
    pub fn new(input: PhaseSpaceTheoryInput) -> Result<Self, String> {
        validate_input(&input)?;

        // Collect fragment masses (atomic units of mass = electron masses).
        let mass_a_amu = fragment_mass_amu(&input.fragment_a)?;
        let mass_b_amu = fragment_mass_amu(&input.fragment_b)?;
        let mass_a_au = mass_a_amu * AMU_TO_ELECTRON_MASS;
        let mass_b_au = mass_b_amu * AMU_TO_ELECTRON_MASS;

        // Collect rotational constants for both fragments in Hartree units.
        // The total list length determines which analytic numerical factor is used,
        // exactly as in MESS.
        let mut rotational_constants_hartree: Vec<f64> = Vec::new();
        push_rotational_constants(&input.fragment_a, &mut rotational_constants_hartree)?;
        push_rotational_constants(&input.fragment_b, &mut rotational_constants_hartree)?;

        // Start with symmetry scaling:
        // MESS uses `_states_factor = 1/symmetry_operations`.
        let mut states_prefactor = 1.0 / input.symmetry_operations;

        let n_rot = rotational_constants_hartree.len();
        let n = input.potential_power_exponent;
        // MESS: power = (number of rotational constants + 2)/2 - 2/n.
        let power = (n_rot as f64 + 2.0) / 2.0 - 2.0 / n;
        if !(2..=6).contains(&n_rot) {
            return Err(format!(
                "Unsupported total number of rotational constants for PST: {n_rot} (expected 2..=6)."
            ));
        }
        // Classical rotor phase-space factor: sqrt(pi) per nonlinear fragment (one nonlinear fragment for
        // 3 or 5 constants, two for 6), as in MESS for the EJ and T levels.
        let nonlinear_rotor_factor = match n_rot {
            3 | 5 => PI.sqrt(),
            6 => PI,
            _ => 1.0,
        };
        match input.tst_level {
            PstTstLevel::E => {
                // MESS E level, microcanonical variational TST (Georgievskii, Klippenstein, J. Chem. Phys.
                // 122, 194103 (2005), eq. 58): c_r d^(2/n) (1 + 1/d)^((r+2)/2), d = (r+2) n/4 - 1.
                states_prefactor *= match n_rot {
                    2 => 1.0,
                    3 => 16.0 / 15.0,
                    4 => 1.0 / 3.0,
                    5 => 32.0 / 105.0,
                    _ => PI / 12.0,
                };
                let internal_dimension_count = n_rot as f64 + 2.0;
                let d_parameter = internal_dimension_count * n / 4.0 - 1.0;
                states_prefactor *= d_parameter.powf(2.0 / n)
                    * (1.0 + 1.0 / d_parameter).powf(internal_dimension_count / 2.0);
            }
            PstTstLevel::EJ => {
                // MESS EJ level (the MESS default), E,J-resolved TST; its canonical capture rate is
                // GK05 eq. 55: 2 ((n-2)/2)^(2/n) Gamma(1 - 2/n) / Gamma(power + 1).
                states_prefactor *= nonlinear_rotor_factor * 2.0 * ((n - 2.0) / 2.0).powf(2.0 / n)
                    * gamma_func(1.0 - 2.0 / n)
                    / gamma_func(power + 1.0);
            }
            PstTstLevel::T => {
                // MESS T level, canonical variational TST: 2 (n/2)^(2/n) exp(2/n) / Gamma(power + 1).
                states_prefactor *= nonlinear_rotor_factor * 2.0 * (n / 2.0).powf(2.0 / n) * (2.0 / n).exp()
                    / gamma_func(power + 1.0);
            }
            PstTstLevel::J0 => {
                return Err("PST TSTLevel J=0 is not supported (as in MESS).".into());
            }
        }

        // Effective mass factor:
        // MESS: dtemp = 1/m1 + 1/m2; states_factor /= dtemp.
        let inverse_reduced_mass_like = (1.0 / mass_a_au) + (1.0 / mass_b_au);
        states_prefactor /= inverse_reduced_mass_like;

        // Rotational constants factor:
        // MESS: states_factor /= sqrt(prod(B_i)).
        let mut rotational_constants_product = 1.0_f64;
        for &b in &rotational_constants_hartree {
            rotational_constants_product *= b;
        }
        states_prefactor /= rotational_constants_product.sqrt();

        // Potential prefactor contribution, MESS: states_factor *= V0^(2/n).
        states_prefactor *= input.potential_prefactor_au.powf(2.0 / n);

        // Canonical weight factor uses Γ(power+1).
        let gamma_factor = gamma_func(power + 1.0);
        let weight_prefactor = states_prefactor * gamma_factor;

        Ok(Self {
            power,
            states_prefactor,
            weight_prefactor,
        })
    }

    /// Cumulative number of states N(E) for energy in cm^-1.
    ///
    /// `energy_wavenumber_cm1` should be the *available energy above threshold*.
    pub fn cumulative_states_at_energy_cm1(
        &self,
        energy_wavenumber_cm1: f64,
    ) -> Result<f64, String> {
        if energy_wavenumber_cm1 <= 0.0 {
            return Ok(0.0);
        }
        if !energy_wavenumber_cm1.is_finite() {
            return Err("Energy must be finite.".into());
        }

        let energy_hartree = energy_wavenumber_cm1 * CM1_TO_HARTREE;
        Ok(self.states_prefactor * energy_hartree.powf(self.power))
    }

    /// Canonical weight Q(T) for temperature in Kelvin.
    ///
    /// This uses k_B*T in energy units. In the MESS formulation, this is consistent because
    /// both N(E) and Q(T) are computed in the same internal energy units.
    ///
    /// Here we accept temperature as an *energy* in Hartree via `k_b_t_hartree`.
    /// If you prefer to pass Kelvin, convert using your chosen k_B convention.
    pub fn canonical_weight_from_kbt_hartree(&self, k_b_t_hartree: f64) -> Result<f64, String> {
        if k_b_t_hartree <= 0.0 {
            return Err("k_B*T must be positive.".into());
        }
        if !k_b_t_hartree.is_finite() {
            return Err("k_B*T must be finite.".into());
        }
        Ok(self.weight_prefactor * k_b_t_hartree.powf(self.power))
    }
}

fn validate_input(input: &PhaseSpaceTheoryInput) -> Result<(), String> {
    if input.symmetry_operations <= 0.0 || !input.symmetry_operations.is_finite() {
        return Err("symmetry_operations must be positive and finite.".into());
    }
    if input.potential_prefactor_au <= 0.0 || !input.potential_prefactor_au.is_finite() {
        return Err("potential_prefactor_au must be positive and finite.".into());
    }
    if input.potential_power_exponent <= 2.0 || !input.potential_power_exponent.is_finite() {
        return Err("potential_power_exponent must be > 2 and finite (as in MESS).".into());
    }
    validate_fragment(&input.fragment_a)?;
    validate_fragment(&input.fragment_b)?;
    Ok(())
}

fn validate_fragment(fragment: &CaptureFragment) -> Result<(), String> {
    match &fragment.rotor {
        CaptureFragmentRotorModel::Atom => Ok(()),
        CaptureFragmentRotorModel::LinearRigidRotor {
            rotational_constant_cm1,
        } => {
            let m = fragment.mass_amu.ok_or_else(|| {
                "Fragment mass_amu is required for non-geometry fragments.".to_string()
            })?;
            if m <= 0.0 || !m.is_finite() {
                return Err("Fragment mass_amu must be positive and finite.".into());
            }
            if *rotational_constant_cm1 <= 0.0 || !rotational_constant_cm1.is_finite() {
                return Err(
                    "Linear rotor rotational_constant_cm1 must be positive and finite.".into(),
                );
            }
            Ok(())
        }
        CaptureFragmentRotorModel::NonlinearRigidRotor {
            rotational_constants_cm1,
        } => {
            let m = fragment.mass_amu.ok_or_else(|| {
                "Fragment mass_amu is required for non-geometry fragments.".to_string()
            })?;
            if m <= 0.0 || !m.is_finite() {
                return Err("Fragment mass_amu must be positive and finite.".into());
            }
            if rotational_constants_cm1
                .iter()
                .any(|x| *x <= 0.0 || !x.is_finite())
            {
                return Err(
                    "Nonlinear rotor rotational_constants_cm1 entries must be positive and finite."
                        .into(),
                );
            }
            Ok(())
        }
        CaptureFragmentRotorModel::GeometryAngstrom {
            symbols,
            coordinates_angstrom,
        } => {
            if symbols.is_empty()
                || coordinates_angstrom.is_empty()
                || symbols.len() != coordinates_angstrom.len()
            {
                return Err(
                    "GeometryAngstrom must have matching non-empty symbols and coordinates.".into(),
                );
            }
            if coordinates_angstrom
                .iter()
                .flat_map(|v| v.iter())
                .any(|x| !x.is_finite())
            {
                return Err("GeometryAngstrom coordinates must be finite.".into());
            }

            // If mass_amu is present, accept it; otherwise we'll infer it from symbols.
            if let Some(m) = fragment.mass_amu {
                if m <= 0.0 || !m.is_finite() {
                    return Err("Fragment mass_amu must be positive and finite.".into());
                }
            }
            // Also validate that we can map all symbols to masses (for inference + inertia).
            let _ = mass_vector_from_symbols_amu(symbols)?;
            Ok(())
        }
    }
}

fn push_rotational_constants(
    fragment: &CaptureFragment,
    out_rotational_constants_hartree: &mut Vec<f64>,
) -> Result<(), String> {
    match &fragment.rotor {
        CaptureFragmentRotorModel::Atom => Ok(()),
        CaptureFragmentRotorModel::LinearRigidRotor {
            rotational_constant_cm1,
        } => {
            out_rotational_constants_hartree.push(rotational_constant_cm1 * CM1_TO_HARTREE);
            out_rotational_constants_hartree.push(rotational_constant_cm1 * CM1_TO_HARTREE);
            Ok(())
        }
        CaptureFragmentRotorModel::NonlinearRigidRotor {
            rotational_constants_cm1,
        } => {
            for b in rotational_constants_cm1 {
                out_rotational_constants_hartree.push(b * CM1_TO_HARTREE);
            }
            Ok(())
        }
        CaptureFragmentRotorModel::GeometryAngstrom {
            symbols,
            coordinates_angstrom,
        } => {
            if symbols.len() == 1 {
                return Ok(());
            }

            // Rotational constants from the principal moments of inertia, as in MESS: B_i = 1/(2 I_i); a
            // fragment with I_min/I_mid < 1e-5 is linear, with B = 1/(2 I_mid) twice.
            let masses_amu = mass_vector_from_symbols_amu(&symbols)?;
            let coords: Vec<[f64; 3]> = coordinates_angstrom.clone();
            let brot = crate::inertia::inertia::get_brot(&coords, &masses_amu);
            // B is proportional to 1/I: sort B descending, i.e. moments ascending (an infinite B is I = 0).
            let mut b_desc: Vec<f64> = brot.into_iter().map(|b| if b.is_finite() { b } else { f64::INFINITY }).collect();
            if b_desc.len() != 3 || b_desc.iter().any(|b| !(*b > 0.0)) {
                return Err("Failed to compute rotational constants from the fragment geometry.".into());
            }
            b_desc.sort_by(|x, y| y.partial_cmp(x).unwrap());
            // I_min/I_mid = B_mid/B_max.
            let linear = !b_desc[0].is_finite() || b_desc[1] / b_desc[0] < 1.0e-5;
            if linear {
                out_rotational_constants_hartree.push(b_desc[1] * CM1_TO_HARTREE);
                out_rotational_constants_hartree.push(b_desc[1] * CM1_TO_HARTREE);
            } else {
                for b in b_desc {
                    out_rotational_constants_hartree.push(b * CM1_TO_HARTREE);
                }
            }
            Ok(())
        }
    }
}

fn fragment_mass_amu(fragment: &CaptureFragment) -> Result<f64, String> {
    if let Some(m) = fragment.mass_amu {
        if m <= 0.0 || !m.is_finite() {
            return Err("Fragment mass_amu must be positive and finite.".into());
        }
        return Ok(m);
    }

    match &fragment.rotor {
        CaptureFragmentRotorModel::GeometryAngstrom { symbols, .. } => {
            let masses = mass_vector_from_symbols_amu(symbols)?;
            Ok(masses.iter().sum::<f64>())
        }
        _ => Err("Fragment mass_amu is required unless GeometryAngstrom is provided.".into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::barrierless::phasespace::types::{CaptureFragment, CaptureFragmentRotorModel, PhaseSpaceTheoryInput, PstTstLevel};

    // Atomic units throughout (hbar = 1, h = 2 pi, k_B = 1): energies and k_B T in Hartree, masses in
    // electron masses, rate coefficients in bohr^3 per atomic time unit.
    const B_LIN_CM1: f64 = 1.44;
    const B_NONLIN_CM1: [f64; 3] = [0.9, 0.31, 0.25];

    fn atom(mass: f64) -> CaptureFragment {
        CaptureFragment { mass_amu: Some(mass), rotor: CaptureFragmentRotorModel::Atom }
    }
    fn linear(mass: f64) -> CaptureFragment {
        CaptureFragment { mass_amu: Some(mass), rotor: CaptureFragmentRotorModel::LinearRigidRotor { rotational_constant_cm1: B_LIN_CM1 } }
    }
    fn nonlinear(mass: f64) -> CaptureFragment {
        CaptureFragment { mass_amu: Some(mass), rotor: CaptureFragmentRotorModel::NonlinearRigidRotor { rotational_constants_cm1: B_NONLIN_CM1 } }
    }

    /// Classical rigid-rotor partition function (symmetry number 1): linear kT/B, nonlinear
    /// pi^(1/2) (kT)^(3/2) / (ABC)^(1/2) (e.g. McQuarrie, Statistical Mechanics (1976), ch. 8).
    fn rotor_partition_function(fragment: &CaptureFragment, kt: f64) -> f64 {
        match &fragment.rotor {
            CaptureFragmentRotorModel::Atom => 1.0,
            CaptureFragmentRotorModel::LinearRigidRotor { rotational_constant_cm1 } => kt / (rotational_constant_cm1 * CM1_TO_HARTREE),
            CaptureFragmentRotorModel::NonlinearRigidRotor { rotational_constants_cm1 } => {
                let abc: f64 = rotational_constants_cm1.iter().map(|b| b * CM1_TO_HARTREE).product();
                PI.sqrt() * kt.powf(1.5) / abc.sqrt()
            }
            CaptureFragmentRotorModel::GeometryAngstrom { .. } => unreachable!(),
        }
    }

    /// Canonical capture rate coefficient from the PST weight Q_PST(T):
    ///   k = (kT/h) Q_PST / ((mu kT / 2 pi)^(3/2) Q_rot,A Q_rot,B)   (atomic units).
    fn capture_rate_from_model(model: &PhaseSpaceTheoryModel, a: &CaptureFragment, b: &CaptureFragment, kt: f64) -> f64 {
        let mu = 1.0 / (1.0 / (a.mass_amu.unwrap() * AMU_TO_ELECTRON_MASS) + 1.0 / (b.mass_amu.unwrap() * AMU_TO_ELECTRON_MASS));
        let q_pst = model.canonical_weight_from_kbt_hartree(kt).unwrap();
        let q_translation = (mu * kt / (2.0 * PI)).powf(1.5);
        kt / (2.0 * PI) * q_pst / (q_translation * rotor_partition_function(a, kt) * rotor_partition_function(b, kt))
    }

    fn input(a: CaptureFragment, b: CaptureFragment, n: f64, level: PstTstLevel) -> PhaseSpaceTheoryInput {
        PhaseSpaceTheoryInput {
            fragment_a: a,
            fragment_b: b,
            symmetry_operations: 1.0,
            potential_prefactor_au: 37.0,
            potential_power_exponent: n,
            tst_level: level,
        }
    }

    fn reduced_mass_au(a: &CaptureFragment, b: &CaptureFragment) -> f64 {
        1.0 / (1.0 / (a.mass_amu.unwrap() * AMU_TO_ELECTRON_MASS) + 1.0 / (b.mass_amu.unwrap() * AMU_TO_ELECTRON_MASS))
    }

    fn fragment_pairs() -> Vec<(CaptureFragment, CaptureFragment)> {
        vec![
            (atom(1.0), linear(32.0)),       // 2 rotational constants
            (atom(1.0), nonlinear(117.0)),   // 3
            (linear(28.0), linear(32.0)),    // 4
            (nonlinear(117.0), linear(32.0)), // 5
            (nonlinear(117.0), nonlinear(33.0)), // 6
        ]
    }

    #[test]
    fn the_default_level_is_ej_as_in_mess() {
        assert_eq!(PstTstLevel::default(), PstTstLevel::EJ);
    }

    #[test]
    fn ej_level_reproduces_the_isotropic_capture_rate_of_georgievskii_klippenstein_eq_55() {
        // k(T) = (8 pi)^(1/2) ((n-2)/2)^(2/n) Gamma(1 - 2/n) mu^(-1/2) V0^(2/n) T^(1/2 - 2/n)
        // (Georgievskii, Klippenstein, J. Chem. Phys. 122, 194103 (2005), eq. 55), independent of the
        // rotational degrees of freedom of the fragments.
        for n in [4.0f64, 6.0] {
            for (a, b) in fragment_pairs() {
                let model = PhaseSpaceTheoryModel::new(input(a.clone(), b.clone(), n, PstTstLevel::EJ)).unwrap();
                for kt in [3.0e-4f64, 1.0e-3, 4.0e-3] {
                    let mu = reduced_mass_au(&a, &b);
                    let expected = (8.0 * PI).sqrt() * ((n - 2.0) / 2.0).powf(2.0 / n) * gamma_func(1.0 - 2.0 / n)
                        * mu.powf(-0.5) * 37.0f64.powf(2.0 / n) * kt.powf(0.5 - 2.0 / n);
                    let got = capture_rate_from_model(&model, &a, &b, kt);
                    assert!((got / expected - 1.0).abs() < 1e-10, "n {n}: {got:e} vs {expected:e}");
                }
            }
        }
    }

    #[test]
    fn ej_level_for_dispersion_gives_the_numerical_coefficient_of_eq_57() {
        // n = 6, V0 = C6: k(T) = 8.55 mu^(-1/2) C6^(1/3) T^(1/6) (GK05 eq. 57).
        let (a, b) = (nonlinear(117.0), linear(32.0));
        let model = PhaseSpaceTheoryModel::new(input(a.clone(), b.clone(), 6.0, PstTstLevel::EJ)).unwrap();
        let kt = 1.0e-3;
        let k = capture_rate_from_model(&model, &a, &b, kt);
        let coefficient = k / (reduced_mass_au(&a, &b).powf(-0.5) * 37.0f64.powf(1.0 / 3.0) * kt.powf(1.0 / 6.0));
        assert!((coefficient - 8.55).abs() < 5e-3, "{coefficient}");
    }

    #[test]
    fn t_level_is_the_canonical_variational_rate_of_an_isotropic_potential() {
        // Canonical variational TST for V = -V0/R^n: the capture flux sqrt(8kT/(pi mu)) pi R^2 exp(V0/(R^n kT))
        // is minimal at R^n = n V0/(2 kT), giving k = sqrt(8kT/(pi mu)) pi (n V0/(2kT))^(2/n) exp(2/n).
        for n in [4.0f64, 6.0] {
            for (a, b) in fragment_pairs() {
                let model = PhaseSpaceTheoryModel::new(input(a.clone(), b.clone(), n, PstTstLevel::T)).unwrap();
                let kt = 1.0e-3;
                let mu = reduced_mass_au(&a, &b);
                let expected = (8.0 * kt / (PI * mu)).sqrt() * PI * (n * 37.0 / (2.0 * kt)).powf(2.0 / n) * (2.0 / n).exp();
                let got = capture_rate_from_model(&model, &a, &b, kt);
                assert!((got / expected - 1.0).abs() < 1e-10, "n {n}: {got:e} vs {expected:e}");
            }
        }
    }

    #[test]
    fn e_level_keeps_the_microcanonical_variational_form_of_eq_58() {
        // E level: prefactor c_r d^(2/n) (1 + 1/d)^((r+2)/2), d = (r+2) n/4 - 1 (GK05 eq. 58), with the
        // rotor factor c_r (1, 16/15, 1/3, 32/105, pi/12 for r = 2..6); relative to the EJ level for r = 5,
        // n = 6 this is 1.125.
        let (a, b) = (nonlinear(117.0), linear(32.0));
        let e = PhaseSpaceTheoryModel::new(input(a.clone(), b.clone(), 6.0, PstTstLevel::E)).unwrap();
        let ej = PhaseSpaceTheoryModel::new(input(a, b, 6.0, PstTstLevel::EJ)).unwrap();
        assert_eq!(e.power, ej.power);
        let d: f64 = 7.0 * 6.0 / 4.0 - 1.0;
        let expected_ratio = (32.0 / 105.0) * d.powf(1.0 / 3.0) * (1.0 + 1.0 / d).powf(3.5)
            / (PI.sqrt() * 2.0 * 2.0f64.powf(1.0 / 3.0) * gamma_func(2.0 / 3.0) / gamma_func(ej.power + 1.0));
        assert!((e.states_prefactor / ej.states_prefactor / expected_ratio - 1.0).abs() < 1e-10);
        assert!((e.states_prefactor / ej.states_prefactor - 1.125).abs() < 1e-3);
    }

    #[test]
    fn a_nearly_linear_geometry_counts_as_linear_with_the_middle_moment() {
        // As in the reference implementation: I_min/I_mid < 1e-5 means linear, B = 1/(2 I_mid).
        let chain = CaptureFragment {
            mass_amu: None,
            rotor: CaptureFragmentRotorModel::GeometryAngstrom {
                symbols: vec!["O".into(), "C".into(), "O".into()],
                coordinates_angstrom: vec![[0.0, 0.0, -1.16], [0.0, 1.0e-7, 0.0], [0.0, 0.0, 1.16]],
            },
        };
        let model = PhaseSpaceTheoryModel::new(input(atom(1.0), chain, 6.0, PstTstLevel::EJ)).unwrap();
        assert!((model.power - (4.0 / 2.0 - 1.0 / 3.0)).abs() < 1e-12, "power {}", model.power);
    }

    #[test]
    fn exponents_not_above_two_and_the_j0_level_are_refused() {
        let (a, b) = (nonlinear(117.0), linear(32.0));
        assert!(PhaseSpaceTheoryModel::new(input(a.clone(), b.clone(), 2.0, PstTstLevel::EJ)).is_err());
        assert!(PhaseSpaceTheoryModel::new(input(a, b, 6.0, PstTstLevel::J0)).is_err());
    }
}
