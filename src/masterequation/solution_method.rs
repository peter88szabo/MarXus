//! Solution method of the master equation and its settings.
//!
//! MarXus has two solution methods:
//! 1. The steady state, J N = F (GO10 eqs. 7 and 8), in two versions (GO10 Sec. 3.2):
//!    - the intermediate steady state, with a lower absorbing barrier below the lowest threshold of each
//!      well ("implemented by introducing a lower absorbing barrier into the master equation");
//!    - the final steady state, without it. It includes the thermal rate coefficients of the same J:
//!      the lowest eigenvalue lambda_1 and the specific rate coefficients averaged over its normalized
//!      eigenvector, the thermal steady-state population (GO10 eq. 12 and the text after it;
//!      `chemical_activation_eigen.rs`).
//! 2. The phenomenological rate coefficients from the chemically significant eigenvalues (CSE; MK06, G13;
//!    `chemically_significant_eigenvalues.rs`).
//!
//! The settings come from the `MarXus ... End` block of the deck header (`mess_input.rs`) and from the
//! command line, which overrides the deck (`SolutionSettings::overridden_by`). `SolutionSettings::resolve`
//! fills in the defaults, refuses impossible combinations and lists the given settings that the selected
//! solution does not use.
//!
//! References:
//! - GO10: G. Gonzalez-Garcia, J. Olzmann, Phys. Chem. Chem. Phys. 12, 12290 (2010).
//! - MK06: J. A. Miller, S. J. Klippenstein, J. Phys. Chem. A 110, 10528 (2006).
//! - G13: Y. Georgievskii, J. A. Miller, M. P. Burke, S. J. Klippenstein, J. Phys. Chem. A 117, 12146 (2013).

use super::chemical_activation_eigen::{EigenSolver, DEFAULT_SUM_RULE_TOLERANCE};

/// Default absorbing barrier of the intermediate steady state: 10 k_BT below the lowest threshold of each well.
pub const DEFAULT_ABSORBING_BARRIER_KT: f64 = 10.0;

/// Solution method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolutionMethod {
    /// Steady state J N = F (GO10 eqs. 7, 8).
    SteadyState,
    /// Phenomenological rate coefficients from the chemically significant eigenvalues (MK06; G13).
    ChemicallySignificantEigenvalues,
}

/// Versions of the steady-state method that are computed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SteadyStateVersions {
    Intermediate,
    Final,
    Both,
}

/// Settings as given in the deck or on the command line; None: not given.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SolutionSettings {
    /// `Method` (`--method`).
    pub method: Option<SolutionMethod>,
    /// `SteadyState` (`--steady-state`).
    pub steady_state: Option<SteadyStateVersions>,
    /// `AbsorbingBarrierBelowThreshold[kT]` (`--barrier-kt`): absorbing barrier of the intermediate steady
    /// state, in k_BT below the lowest threshold of each well.
    pub absorbing_barrier_kt: Option<f64>,
    /// `EigenSolver` (`--eigen-solver`).
    pub eigen_solver: Option<EigenSolver>,
    /// `SumRuleTolerance` (`--sum-rule-tolerance`): relative deviation |lambda_1 - k_uni| / k_uni above which
    /// a warning is given.
    pub sum_rule_tolerance: Option<f64>,
}

/// Eigen-solver and sum-rule tolerance of the thermal rate coefficients of the final steady state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermalEigenSettings {
    pub eigen_solver: EigenSolver,
    pub sum_rule_tolerance: f64,
}

/// The solution to compute, with all its settings.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Solution {
    SteadyState {
        /// Absorbing barrier (k_BT below the lowest threshold) if the intermediate steady state is computed.
        intermediate_absorbing_barrier_kt: Option<f64>,
        /// Settings of the thermal rate coefficients if the final steady state is computed.
        final_steady_state: Option<ThermalEigenSettings>,
    },
    ChemicallySignificantEigenvalues { eigen_solver: EigenSolver },
}

/// The resolved solution and notes on the given settings that it does not use.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedSolution {
    pub solution: Solution,
    pub unused_settings: Vec<String>,
}

/// Keyword value without case, '-' and '_': `SteadyState`, `steady-state` and `steady_state` are equal.
fn normalized(value: &str) -> String {
    value.chars().filter(|c| *c != '-' && *c != '_').collect::<String>().to_lowercase()
}

/// Explanation for a request of the eigenvalue analysis as a method or a version of its own.
const EIGENVALUE_ANALYSIS_IS_PART_OF_THE_FINAL_STEADY_STATE: &str =
    "the thermal eigenvalue analysis is not a method or version of its own: it is part of the final steady \
     state, whose J gives the thermal rate coefficients by its lowest eigenpair (GO10 eq. 12); use Method \
     SteadyState with SteadyState Final or Both (--method steady-state --steady-state final|both)";

impl SolutionMethod {
    /// `SteadyState` (`steady-state`) or `CSE` / `ChemicallySignificantEigenvalues`.
    pub fn from_keyword(value: &str) -> Result<Self, String> {
        match normalized(value).as_str() {
            "steadystate" => Ok(Self::SteadyState),
            "cse" | "chemicallysignificanteigenvalues" => Ok(Self::ChemicallySignificantEigenvalues),
            "eigenvalue" | "eigenvalues" => Err(format!("Method '{value}': {EIGENVALUE_ANALYSIS_IS_PART_OF_THE_FINAL_STEADY_STATE}.")),
            _ => Err(format!("Method '{value}': unknown (SteadyState or CSE).")),
        }
    }
}

impl SteadyStateVersions {
    /// `Intermediate`, `Final` or `Both`.
    pub fn from_keyword(value: &str) -> Result<Self, String> {
        match normalized(value).as_str() {
            "intermediate" => Ok(Self::Intermediate),
            "final" => Ok(Self::Final),
            "both" => Ok(Self::Both),
            "eigenvalue" | "eigenvalues" => {
                Err(format!("SteadyState '{value}': {EIGENVALUE_ANALYSIS_IS_PART_OF_THE_FINAL_STEADY_STATE}."))
            }
            "cse" | "chemicallysignificanteigenvalues" => Err(format!(
                "SteadyState '{value}': the CSE analysis is the other solution method, not a version of the steady \
                 state; use Method CSE (--method cse)."
            )),
            _ => Err(format!("SteadyState '{value}': unknown (Intermediate, Final or Both).")),
        }
    }
}

/// `InverseIteration` (`inverse`), `FullDecomposition` (`full`) or `Lapack`.
pub fn eigen_solver_from_keyword(value: &str) -> Result<EigenSolver, String> {
    match normalized(value).as_str() {
        "inverseiteration" | "inverse" => Ok(EigenSolver::InverseIteration),
        "fulldecomposition" | "full" => Ok(EigenSolver::FullDecomposition),
        "lapack" => Ok(EigenSolver::FullDecompositionLapack),
        _ => Err(format!("EigenSolver '{value}': unknown (InverseIteration, FullDecomposition or Lapack).")),
    }
}

/// Deck keyword and command-line option of each setting, for the messages.
const METHOD: &str = "Method (--method)";
const STEADY_STATE: &str = "SteadyState (--steady-state)";
const ABSORBING_BARRIER: &str = "AbsorbingBarrierBelowThreshold[kT] (--barrier-kt)";
const EIGEN_SOLVER: &str = "EigenSolver (--eigen-solver)";
const SUM_RULE_TOLERANCE: &str = "SumRuleTolerance (--sum-rule-tolerance)";

impl SolutionSettings {
    /// These settings with those that `other` gives replacing them (the command line over the deck).
    pub fn overridden_by(&self, other: &SolutionSettings) -> SolutionSettings {
        SolutionSettings {
            method: other.method.or(self.method),
            steady_state: other.steady_state.or(self.steady_state),
            absorbing_barrier_kt: other.absorbing_barrier_kt.or(self.absorbing_barrier_kt),
            eigen_solver: other.eigen_solver.or(self.eigen_solver),
            sum_rule_tolerance: other.sum_rule_tolerance.or(self.sum_rule_tolerance),
        }
    }

    /// The solution with the defaults filled in. Errors: a non-positive or non-finite number, and the
    /// CSE method with inverse iteration (it needs all eigenpairs). Settings that the selected solution
    /// does not use are listed in `unused_settings`.
    pub fn resolve(&self) -> Result<ResolvedSolution, String> {
        for (name, value) in [(ABSORBING_BARRIER, self.absorbing_barrier_kt), (SUM_RULE_TOLERANCE, self.sum_rule_tolerance)] {
            if let Some(v) = value.filter(|v| !(v.is_finite() && *v > 0.0)) {
                return Err(format!("{name}: {v} is not a positive number."));
            }
        }
        let mut unused_settings = Vec::new();
        let mut note_unused = |name: &str, given: bool, user: &str| {
            if given {
                unused_settings.push(format!("{name} is not used: it applies only to {user}."));
            }
        };
        let method = self.method.unwrap_or(SolutionMethod::SteadyState);
        let solution = match method {
            SolutionMethod::SteadyState => {
                let versions = self.steady_state.unwrap_or(SteadyStateVersions::Both);
                let intermediate = versions != SteadyStateVersions::Final;
                let final_steady_state = versions != SteadyStateVersions::Intermediate;
                let barrier_user = "the intermediate steady state (absorbing barrier)";
                note_unused(ABSORBING_BARRIER, !intermediate && self.absorbing_barrier_kt.is_some(), barrier_user);
                let thermal_user = "the final steady state (thermal eigenpair, GO10 eq. 12) and the CSE method";
                note_unused(EIGEN_SOLVER, !final_steady_state && self.eigen_solver.is_some(), thermal_user);
                let sum_rule_user = "the final steady state (thermal eigenpair, GO10 eq. 12)";
                note_unused(SUM_RULE_TOLERANCE, !final_steady_state && self.sum_rule_tolerance.is_some(), sum_rule_user);
                Solution::SteadyState {
                    intermediate_absorbing_barrier_kt: intermediate
                        .then(|| self.absorbing_barrier_kt.unwrap_or(DEFAULT_ABSORBING_BARRIER_KT)),
                    final_steady_state: final_steady_state.then(|| ThermalEigenSettings {
                        eigen_solver: self.eigen_solver.unwrap_or_default(),
                        sum_rule_tolerance: self.sum_rule_tolerance.unwrap_or(DEFAULT_SUM_RULE_TOLERANCE),
                    }),
                }
            }
            SolutionMethod::ChemicallySignificantEigenvalues => {
                let steady_state_user = "the steady-state method (Method SteadyState)";
                note_unused(STEADY_STATE, self.steady_state.is_some(), steady_state_user);
                note_unused(ABSORBING_BARRIER, self.absorbing_barrier_kt.is_some(), "the intermediate steady state");
                let sum_rule_user = "the final steady state (thermal eigenpair, GO10 eq. 12)";
                note_unused(SUM_RULE_TOLERANCE, self.sum_rule_tolerance.is_some(), sum_rule_user);
                // All eigenpairs are needed (G13 eqs. 21, 25-30): LAPACK unless the in-house decomposition is asked for.
                let eigen_solver = match self.eigen_solver {
                    None => EigenSolver::FullDecompositionLapack,
                    Some(EigenSolver::InverseIteration) => {
                        return Err(format!(
                            "{METHOD} CSE needs a full decomposition (all eigenpairs), which inverse iteration does \
                             not give: set {EIGEN_SOLVER} to Lapack or FullDecomposition."
                        ))
                    }
                    Some(solver) => solver,
                };
                Solution::ChemicallySignificantEigenvalues { eigen_solver }
            }
        };
        Ok(ResolvedSolution { solution, unused_settings })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keywords_of_the_deck_and_of_the_command_line_are_both_accepted() {
        assert_eq!(SolutionMethod::from_keyword("SteadyState").unwrap(), SolutionMethod::SteadyState);
        assert_eq!(SolutionMethod::from_keyword("steady-state").unwrap(), SolutionMethod::SteadyState);
        assert_eq!(SolutionMethod::from_keyword("CSE").unwrap(), SolutionMethod::ChemicallySignificantEigenvalues);
        assert_eq!(
            SolutionMethod::from_keyword("ChemicallySignificantEigenvalues").unwrap(),
            SolutionMethod::ChemicallySignificantEigenvalues
        );
        assert_eq!(SteadyStateVersions::from_keyword("Final").unwrap(), SteadyStateVersions::Final);
        assert_eq!(SteadyStateVersions::from_keyword("intermediate").unwrap(), SteadyStateVersions::Intermediate);
        assert_eq!(SteadyStateVersions::from_keyword("Both").unwrap(), SteadyStateVersions::Both);
        assert_eq!(eigen_solver_from_keyword("InverseIteration").unwrap(), EigenSolver::InverseIteration);
        assert_eq!(eigen_solver_from_keyword("inverse").unwrap(), EigenSolver::InverseIteration);
        assert_eq!(eigen_solver_from_keyword("FullDecomposition").unwrap(), EigenSolver::FullDecomposition);
        assert_eq!(eigen_solver_from_keyword("full").unwrap(), EigenSolver::FullDecomposition);
        assert_eq!(eigen_solver_from_keyword("Lapack").unwrap(), EigenSolver::FullDecompositionLapack);
        assert!(eigen_solver_from_keyword("jacobi").is_err());
    }

    #[test]
    fn the_eigenvalue_analysis_and_cse_are_not_versions_of_the_steady_state() {
        // The thermal eigenpair belongs to the final steady state (GO10 eq. 12); CSE is the other method.
        assert!(SteadyStateVersions::from_keyword("eigenvalue").unwrap_err().contains("final steady state"));
        assert!(SteadyStateVersions::from_keyword("cse").unwrap_err().contains("Method"));
        assert!(SteadyStateVersions::from_keyword("all").is_err());
        assert!(SolutionMethod::from_keyword("eigenvalue").unwrap_err().contains("final steady state"));
    }

    #[test]
    fn the_default_is_the_steady_state_method_in_both_versions() {
        let resolved = SolutionSettings::default().resolve().unwrap();
        assert_eq!(
            resolved.solution,
            Solution::SteadyState {
                intermediate_absorbing_barrier_kt: Some(DEFAULT_ABSORBING_BARRIER_KT),
                final_steady_state: Some(ThermalEigenSettings {
                    eigen_solver: EigenSolver::InverseIteration,
                    sum_rule_tolerance: DEFAULT_SUM_RULE_TOLERANCE,
                }),
            }
        );
        assert!(resolved.unused_settings.is_empty());
    }

    #[test]
    fn the_cse_method_uses_lapack_by_default_and_refuses_inverse_iteration() {
        let cse = SolutionSettings { method: Some(SolutionMethod::ChemicallySignificantEigenvalues), ..Default::default() };
        assert_eq!(
            cse.resolve().unwrap().solution,
            Solution::ChemicallySignificantEigenvalues { eigen_solver: EigenSolver::FullDecompositionLapack }
        );
        let full = SolutionSettings { eigen_solver: Some(EigenSolver::FullDecomposition), ..cse };
        assert_eq!(
            full.resolve().unwrap().solution,
            Solution::ChemicallySignificantEigenvalues { eigen_solver: EigenSolver::FullDecomposition }
        );
        let inverse = SolutionSettings { eigen_solver: Some(EigenSolver::InverseIteration), ..cse };
        assert!(inverse.resolve().unwrap_err().contains("full decomposition"));
    }

    #[test]
    fn settings_the_selected_solution_does_not_use_are_listed_not_refused() {
        // A deck that lists every keyword can be switched to the CSE method.
        let template = SolutionSettings {
            method: Some(SolutionMethod::ChemicallySignificantEigenvalues),
            steady_state: Some(SteadyStateVersions::Both),
            absorbing_barrier_kt: Some(10.0),
            eigen_solver: Some(EigenSolver::FullDecompositionLapack),
            sum_rule_tolerance: Some(1.5e-2),
        };
        let unused = template.resolve().unwrap().unused_settings;
        assert_eq!(unused.len(), 3, "{unused:?}");
        assert!(unused.iter().any(|n| n.contains("SteadyState (--steady-state)")));
        assert!(unused.iter().any(|n| n.contains("AbsorbingBarrierBelowThreshold[kT] (--barrier-kt)")));
        assert!(unused.iter().any(|n| n.contains("SumRuleTolerance (--sum-rule-tolerance)")));
        // The absorbing barrier exists only in the intermediate steady state, the eigenpair only in the final one.
        let final_only = SolutionSettings {
            steady_state: Some(SteadyStateVersions::Final),
            absorbing_barrier_kt: Some(5.0),
            ..Default::default()
        };
        assert_eq!(final_only.resolve().unwrap().unused_settings.len(), 1);
        let intermediate_only = SolutionSettings {
            steady_state: Some(SteadyStateVersions::Intermediate),
            absorbing_barrier_kt: Some(5.0),
            eigen_solver: Some(EigenSolver::FullDecomposition),
            sum_rule_tolerance: Some(0.01),
            ..Default::default()
        };
        let resolved = intermediate_only.resolve().unwrap();
        assert_eq!(
            resolved.solution,
            Solution::SteadyState { intermediate_absorbing_barrier_kt: Some(5.0), final_steady_state: None }
        );
        assert_eq!(resolved.unused_settings.len(), 2);
    }

    #[test]
    fn non_positive_numbers_are_refused() {
        let barrier = SolutionSettings { absorbing_barrier_kt: Some(0.0), ..Default::default() };
        assert!(barrier.resolve().unwrap_err().contains("AbsorbingBarrierBelowThreshold[kT]"));
        let tolerance = SolutionSettings { sum_rule_tolerance: Some(f64::NAN), ..Default::default() };
        assert!(tolerance.resolve().unwrap_err().contains("SumRuleTolerance"));
    }

    #[test]
    fn the_command_line_overrides_the_deck() {
        let deck = SolutionSettings {
            method: Some(SolutionMethod::SteadyState),
            steady_state: Some(SteadyStateVersions::Final),
            sum_rule_tolerance: Some(0.02),
            ..Default::default()
        };
        let command_line = SolutionSettings {
            steady_state: Some(SteadyStateVersions::Both),
            absorbing_barrier_kt: Some(5.0),
            ..Default::default()
        };
        assert_eq!(
            deck.overridden_by(&command_line),
            SolutionSettings {
                method: Some(SolutionMethod::SteadyState),
                steady_state: Some(SteadyStateVersions::Both),
                absorbing_barrier_kt: Some(5.0),
                eigen_solver: None,
                sum_rule_tolerance: Some(0.02),
            }
        );
    }
}
