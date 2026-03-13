/// Phase-space / capture-style transition-state models.
///
/// This module is intentionally focused on *generic* (model-agnostic) inputs:
/// masses, rotational constants, symmetry operations, and long-range potential parameters.
/// It is meant to be reusable wherever we need a "loose TS" number-of-states model.

/// How a fragment contributes rotational degrees of freedom to the capture complex.
///
/// The PST model needs only:
/// - fragment mass
/// - whether the fragment is linear or nonlinear (to decide if it has 1 or 3 rotational constants)
/// - the rotational constants themselves (in cm^-1)
#[derive(Clone, Debug)]
pub enum CaptureFragmentRotorModel {
    /// Atomic fragment: no rotational constants.
    Atom,

    /// Linear rigid rotor: one independent rotational constant (cm^-1).
    LinearRigidRotor { rotational_constant_cm1: f64 },

    /// Nonlinear rigid rotor: three rotational constants (cm^-1).
    NonlinearRigidRotor { rotational_constants_cm1: [f64; 3] },

    /// Geometry-based fragment definition (Å).
    ///
    /// If you provide a geometry, the PST code can:
    /// - infer the fragment mass from element symbols (if `mass_amu` is omitted), and
    /// - compute effective rotational constants from the inertia tensor.
    ///
    /// This mirrors the convenience behavior of MESS' `PhaseSpaceTheory`, where a fragment
    /// can be specified by `FragmentGeometry` rather than explicitly listing rotational constants.
    GeometryAngstrom {
        symbols: Vec<String>,
        coordinates_angstrom: Vec<[f64; 3]>,
    },
}

/// Fragment description for PST capture.
#[derive(Clone, Debug)]
pub struct CaptureFragment {
    /// Fragment mass in amu.
    ///
    /// If you specify `rotor = GeometryAngstrom { ... }`, this can be omitted and will
    /// be inferred from the element symbols.
    pub mass_amu: Option<f64>,

    /// Rotational model and constants.
    pub rotor: CaptureFragmentRotorModel,
}

/// Inputs for a MESS-like Phase Space Theory (PST) number-of-states model.
///
/// This follows the analytical loose-TS capture model used in MESS (`PhaseSpaceTheory`):
/// a long-range potential V(R) = V0 / R^n with fragment rigid-rotor approximations.
#[derive(Clone, Debug)]
pub struct PhaseSpaceTheoryInput {
    pub fragment_a: CaptureFragment,
    pub fragment_b: CaptureFragment,

    /// Number of symmetry operations (MESS' `SymmetryFactor`).
    /// For example, if there are 2 indistinguishable orientations, set this to 2.
    /// The PST state count is scaled by 1 / symmetry_operations.
    pub symmetry_operations: f64,

    /// Long-range potential prefactor V0 in atomic units.
    ///
    /// In MESS input this is `PotentialPrefactor[au]`.
    pub potential_prefactor_au: f64,

    /// Long-range potential power exponent n in V(R) = V0 / R^n.
    ///
    /// In MESS input this is `PotentialPowerExponent`.
    pub potential_power_exponent: f64,
}
