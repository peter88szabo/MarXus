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

use super::types::{CaptureFragment, CaptureFragmentRotorModel, PhaseSpaceTheoryInput};

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

        // Numerical factor depends on the total number of rotational constants.
        // This encodes different combinations of linear/nonlinear fragments.
        match rotational_constants_hartree.len() {
            2 => {
                // (2 + 0): one linear rotor + one atom
            }
            3 => {
                // (3 + 0): one nonlinear rotor + one atom
                states_prefactor *= 16.0 / 15.0;
            }
            4 => {
                // (2 + 2): two linear rotors
                states_prefactor *= 1.0 / 3.0;
            }
            5 => {
                // (2 + 3): one linear + one nonlinear rotor
                states_prefactor *= 32.0 / 105.0;
            }
            6 => {
                // (3 + 3): two nonlinear rotors
                states_prefactor *= PI / 12.0;
            }
            other => {
                return Err(format!(
                    "Unsupported total number of rotational constants for PST: {} (expected 2..=6).",
                    other
                ));
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

        // Potential / long-range part:
        //
        // MESS uses:
        //   itemp = n_rot_constants + 2
        //   dtemp = itemp * n / 4 - 1
        //   states_factor *= (dtemp*V0)^(2/n) * (1 + 1/dtemp)^(itemp/2)
        //   power = itemp/2 - 2/n
        //
        // where `V0` is in atomic units and n is the potential exponent.
        let internal_dimension_count = (rotational_constants_hartree.len() as f64) + 2.0;
        let n = input.potential_power_exponent;

        let d_parameter = internal_dimension_count * n / 4.0 - 1.0;
        if !(d_parameter > 0.0) || !d_parameter.is_finite() {
            return Err(format!(
                "Invalid PST parameter d = (D*n/4 - 1) = {}. Check potential exponent and fragment model.",
                d_parameter
            ));
        }

        let v0 = input.potential_prefactor_au;
        states_prefactor *= (d_parameter * v0).powf(2.0 / n)
            * (1.0 + 1.0 / d_parameter).powf(internal_dimension_count / 2.0);

        let power = internal_dimension_count / 2.0 - 2.0 / n;

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
    if input.potential_power_exponent <= 1.0 || !input.potential_power_exponent.is_finite() {
        return Err("potential_power_exponent must be > 1 and finite.".into());
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

            // Compute rotational constants from the inertia tensor using the existing inertia module.
            // This provides principal-axis rotational constants in cm^-1.
            let masses_amu = mass_vector_from_symbols_amu(&symbols)?;
            let coords: Vec<[f64; 3]> = coordinates_angstrom.clone();
            let brot = crate::inertia::inertia::get_brot(&coords, &masses_amu);

            // Interpret the inertia-derived constants:
            // - Nonlinear fragments: 3 finite constants.
            // - Linear fragments (incl. diatomics): the inertia tensor has one ~0 principal moment,
            //   which makes one rotational constant blow up. We drop non-finite/huge values and
            //   treat the remaining two as the degenerate linear constant.
            let mut finite: Vec<f64> = brot
                .into_iter()
                .filter(|b| b.is_finite() && *b > 0.0 && *b < 1.0e6)
                .collect();

            if finite.is_empty() {
                return Err("Failed to compute finite rotational constants from geometry.".into());
            }

            if finite.len() == 1 {
                // Should not happen for a real polyatomic; treat as atom-like.
                return Ok(());
            }

            if finite.len() == 2 {
                let b = 0.5 * (finite[0] + finite[1]);
                out_rotational_constants_hartree.push(b * CM1_TO_HARTREE);
                out_rotational_constants_hartree.push(b * CM1_TO_HARTREE);
                return Ok(());
            }

            // 3 or more finite (should be exactly 3). Keep first 3.
            finite.truncate(3);
            for b in finite {
                out_rotational_constants_hartree.push(b * CM1_TO_HARTREE);
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
