//! Solution methods of the master equation and their settings.
//!
//! MarXus has four solvers of the same master equation (README, "Master-equation solvers"):
//! 1. `SteadyStateOlzmann`: the final steady state, J N = R F without an absorbing barrier (GO10 eqs. 7, 8,
//!    Sec. 3.2), with the thermal rate coefficients of the same J from its lowest eigenpair (GO10 eq. 12 and
//!    the text after it; `chemical_activation_eigen.rs`) and the thermal fates of the wells.
//! 2. `SteadyStateAbsorbingBarrier`: the intermediate steady state, J N = F with a lower absorbing barrier
//!    below the lowest threshold of each well (GO10 Sec. 3.2: "implemented by introducing a lower absorbing
//!    barrier into the master equation").
//!    The two belong to the steady-state family but are different ways of solving: one run computes one.
//! 3. `CSE`: the phenomenological rate coefficients from the chemically significant eigenvalues (MK06, G13;
//!    `chemically_significant_eigenvalues.rs`).
//! 4. `TimeIntegration`: the direct time integration of the grained populations (`direct_time_integration.rs`).
//!
//! The method must be given (`Method` in the `MarXus ... End` block of the deck header, `mess_input.rs`, or
//! `--method`); there is no default. The command line overrides the deck (`SolutionSettings::overridden_by`).
//! `SolutionSettings::resolve` fills in the defaults of the chosen method, refuses impossible combinations and
//! lists the given settings that the chosen method does not use.
//!
//! References:
//! - GO10: G. Gonzalez-Garcia, J. Olzmann, Phys. Chem. Chem. Phys. 12, 12290 (2010).
//! - MK06: J. A. Miller, S. J. Klippenstein, J. Phys. Chem. A 110, 10528 (2006).
//! - G13: Y. Georgievskii, J. A. Miller, M. P. Burke, S. J. Klippenstein, J. Phys. Chem. A 117, 12146 (2013).

use super::chemical_activation_eigen::{EigenSolver, DEFAULT_SUM_RULE_TOLERANCE};
use super::direct_time_integration::InitialState;
use crate::numeric::integrators::rosenbrock_methods::RosenbrockMethod;

/// Default absorbing barrier of the intermediate steady state: 10 k_BT below the lowest threshold of each well.
pub const DEFAULT_ABSORBING_BARRIER_KT: f64 = 10.0;
/// Default time range of the time integration, s.
pub const DEFAULT_TIME_RANGE_S: (f64, f64) = (1.0e-12, 1.0e2);
/// Default output times per decade.
pub const DEFAULT_TIMES_PER_DECADE: usize = 4;
/// Default relative and absolute tolerance of the time integration.
pub const DEFAULT_INTEGRATION_TOLERANCE: f64 = 1.0e-6;
pub const INTEGRATION_ABSOLUTE_TOLERANCE: f64 = 1.0e-14;

/// The four solvers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolutionMethod {
    /// Final steady state (Olzmann), with its thermal eigenpair and the thermal fates of the wells.
    SteadyStateOlzmann,
    /// Intermediate steady state with an absorbing barrier.
    SteadyStateAbsorbingBarrier,
    /// Phenomenological rate coefficients from the chemically significant eigenvalues (MK06; G13).
    ChemicallySignificantEigenvalues,
    /// Direct time integration of the grained populations.
    TimeIntegration,
}

/// Settings as given in the deck or on the command line; None: not given.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SolutionSettings {
    /// `Method` (`--method`): required.
    pub method: Option<SolutionMethod>,
    /// `AbsorbingBarrierBelowThreshold[kT]` (`--barrier-kt`): absorbing barrier of SteadyStateAbsorbingBarrier,
    /// in k_BT below the lowest threshold of each well.
    pub absorbing_barrier_kt: Option<f64>,
    /// `EigenSolver` (`--eigen-solver`): thermal eigenpair of SteadyStateOlzmann, or all eigenpairs for CSE.
    pub eigen_solver: Option<EigenSolver>,
    /// `SumRuleTolerance` (`--sum-rule-tolerance`): relative deviation |lambda_1 - k_uni| / k_uni above which
    /// a warning is given (SteadyStateOlzmann).
    pub sum_rule_tolerance: Option<f64>,
    /// `Integrator` (`--integrator`): Rosenbrock method of the time integration.
    pub integrator: Option<RosenbrockMethod>,
    /// `InitialState` (`--initial`): a pulse of the source or continuous formation.
    pub initial_state: Option<InitialState>,
    /// `TimeRange[s]` (`--time-range`): first and last output time of the time integration.
    pub time_range_s: Option<(f64, f64)>,
    /// `TimesPerDecade` (`--times-per-decade`): output times per decade.
    pub times_per_decade: Option<usize>,
    /// `IntegrationTolerance` (`--integration-tolerance`): relative tolerance of the time integration.
    pub integration_tolerance: Option<f64>,
    /// `NCores` (`--ncore`): processor cores of the run, for every method. The conditions (T, p) are
    /// computed in batches of up to this many at a time, and each condition's LAPACK calls get the cores
    /// left over (cores / conditions at a time). Default: RAYON_NUM_THREADS, otherwise all logical cores.
    pub cores: Option<usize>,
}

/// Eigen-solver and sum-rule tolerance of the thermal eigenpair of SteadyStateOlzmann.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermalEigenSettings {
    pub eigen_solver: EigenSolver,
    pub sum_rule_tolerance: f64,
}

/// The resolved settings of a time integration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TimeIntegrationPlan {
    pub integrator: RosenbrockMethod,
    pub initial_state: InitialState,
    pub time_range_s: (f64, f64),
    pub times_per_decade: usize,
    pub relative_tolerance: f64,
    pub absolute_tolerance: f64,
}

/// The solver to run, with all its settings.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Solution {
    SteadyStateOlzmann(ThermalEigenSettings),
    SteadyStateAbsorbingBarrier { absorbing_barrier_kt: f64 },
    ChemicallySignificantEigenvalues { eigen_solver: EigenSolver },
    TimeIntegration(TimeIntegrationPlan),
}

/// The resolved solution and notes on the given settings that it does not use.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedSolution {
    pub solution: Solution,
    pub unused_settings: Vec<String>,
}

/// Keyword value without case, '-' and '_': `SteadyStateOlzmann`, `steady-state-olzmann` are equal.
fn normalized(value: &str) -> String {
    value.chars().filter(|c| *c != '-' && *c != '_').collect::<String>().to_lowercase()
}

impl SolutionMethod {
    /// `SteadyStateOlzmann`, `SteadyStateAbsorbingBarrier`, `CSE` (`ChemicallySignificantEigenvalues`) or
    /// `TimeIntegration`; on the command line also `steady-state-olzmann`, ...
    pub fn from_keyword(value: &str) -> Result<Self, String> {
        match normalized(value).as_str() {
            "steadystateolzmann" => Ok(Self::SteadyStateOlzmann),
            "steadystateabsorbingbarrier" => Ok(Self::SteadyStateAbsorbingBarrier),
            "cse" | "chemicallysignificanteigenvalues" => Ok(Self::ChemicallySignificantEigenvalues),
            "timeintegration" => Ok(Self::TimeIntegration),
            "steadystate" => Err(format!("Method '{value}': {STEADY_STATE_KEYWORD_REPLACED}.")),
            "eigenvalue" | "eigenvalues" => Err(format!(
                "Method '{value}': the thermal eigenvalue analysis is not a method of its own: it is part of the final \
                 steady state, whose J gives the thermal rate coefficients by its lowest eigenpair (GO10 eq. 12); use \
                 Method SteadyStateOlzmann (--method steady-state-olzmann)."
            )),
            _ => Err(format!("Method '{value}': unknown; {METHODS}.")),
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

/// `Rodas4` (default), `Rodas3`, `Ros4`, `Ros3` or `Ros2`.
pub fn integrator_from_keyword(value: &str) -> Result<RosenbrockMethod, String> {
    match normalized(value).as_str() {
        "rodas4" => Ok(RosenbrockMethod::Rodas4),
        "rodas3" => Ok(RosenbrockMethod::Rodas3),
        "ros4" => Ok(RosenbrockMethod::Ros4),
        "ros3" => Ok(RosenbrockMethod::Ros3),
        "ros2" => Ok(RosenbrockMethod::Ros2),
        _ => Err(format!(
            "Integrator '{value}': unknown (Rodas4, Rodas3, Ros4, Ros3 or Ros2)."
        )),
    }
}

/// `Pulse` or `Continuous`.
pub fn initial_state_from_keyword(value: &str) -> Result<InitialState, String> {
    match normalized(value).as_str() {
        "pulse" => Ok(InitialState::Pulse),
        "continuous" | "continuousformation" => Ok(InitialState::ContinuousFormation),
        _ => Err(format!(
            "InitialState '{value}': unknown (Pulse or Continuous)."
        )),
    }
}

/// The four methods, for the messages.
const METHODS: &str = "Method (--method) is one of SteadyStateOlzmann (steady-state-olzmann), SteadyStateAbsorbingBarrier \
     (steady-state-absorbing-barrier), CSE (cse) and TimeIntegration (time-integration)";

/// Deck keyword and command-line option of each setting, for the messages.
const ABSORBING_BARRIER: &str = "AbsorbingBarrierBelowThreshold[kT] (--barrier-kt)";
const EIGEN_SOLVER: &str = "EigenSolver (--eigen-solver)";
const SUM_RULE_TOLERANCE: &str = "SumRuleTolerance (--sum-rule-tolerance)";
const INTEGRATOR: &str = "Integrator (--integrator)";
const INITIAL_STATE: &str = "InitialState (--initial)";
const TIME_RANGE: &str = "TimeRange[s] (--time-range)";
const TIMES_PER_DECADE: &str = "TimesPerDecade (--times-per-decade)";
const INTEGRATION_TOLERANCE: &str = "IntegrationTolerance (--integration-tolerance)";

/// Explanation for the removed keyword `SteadyState` (`--steady-state`).
pub const STEADY_STATE_KEYWORD_REPLACED: &str =
    "the steady state is a family of two solvers, chosen by the method: Method SteadyStateOlzmann (final steady \
     state, --method steady-state-olzmann) or Method SteadyStateAbsorbingBarrier (intermediate steady state, \
     --method steady-state-absorbing-barrier); one run computes one of them";

impl SolutionSettings {
    /// These settings with those that `other` gives replacing them (the command line over the deck).
    pub fn overridden_by(&self, other: &SolutionSettings) -> SolutionSettings {
        SolutionSettings {
            method: other.method.or(self.method),
            absorbing_barrier_kt: other.absorbing_barrier_kt.or(self.absorbing_barrier_kt),
            eigen_solver: other.eigen_solver.or(self.eigen_solver),
            sum_rule_tolerance: other.sum_rule_tolerance.or(self.sum_rule_tolerance),
            integrator: other.integrator.or(self.integrator),
            initial_state: other.initial_state.or(self.initial_state),
            time_range_s: other.time_range_s.or(self.time_range_s),
            times_per_decade: other.times_per_decade.or(self.times_per_decade),
            integration_tolerance: other.integration_tolerance.or(self.integration_tolerance),
            cores: other.cores.or(self.cores),
        }
    }

    /// The solver with the defaults of its settings filled in. Errors: no method, a non-positive or
    /// non-finite number, and CSE with inverse iteration (it needs all eigenpairs). Settings that the chosen
    /// method does not use are listed in `unused_settings`.
    pub fn resolve(&self) -> Result<ResolvedSolution, String> {
        let method = self
            .method
            .ok_or_else(|| format!("No solution method given (there is no default): {METHODS}."))?;
        if self.cores == Some(0) {
            return Err("NCores (--ncore): the number of cores must be at least 1.".into());
        }
        for (name, value) in [
            (ABSORBING_BARRIER, self.absorbing_barrier_kt),
            (SUM_RULE_TOLERANCE, self.sum_rule_tolerance),
        ] {
            if let Some(v) = value.filter(|v| !(v.is_finite() && *v > 0.0)) {
                return Err(format!("{name}: {v} is not a positive number."));
            }
        }
        if let Some((t_min, t_max)) = self.time_range_s {
            if !(t_min > 0.0 && t_max > t_min && t_max.is_finite()) {
                return Err(format!(
                    "{TIME_RANGE}: {t_min} .. {t_max} s is not a range 0 < t_min < t_max."
                ));
            }
        }
        if self.times_per_decade == Some(0) {
            return Err(format!(
                "{TIMES_PER_DECADE}: at least one output time per decade."
            ));
        }
        if let Some(r) = self
            .integration_tolerance
            .filter(|r| !(*r > 1e-14 && *r < 1.0))
        {
            return Err(format!(
                "{INTEGRATION_TOLERANCE}: {r} is not a relative tolerance between 1e-14 and 1."
            ));
        }

        // Settings of the other methods are noted, not refused.
        let mut unused_settings = Vec::new();
        let mut note = |name: &str, given: bool, users: &str| {
            if given {
                unused_settings.push(format!("{name} is not used: it applies only to {users}."));
            }
        };
        let olzmann = method == SolutionMethod::SteadyStateOlzmann;
        let barrier = method == SolutionMethod::SteadyStateAbsorbingBarrier;
        let cse = method == SolutionMethod::ChemicallySignificantEigenvalues;
        let time = method == SolutionMethod::TimeIntegration;
        note(
            ABSORBING_BARRIER,
            !barrier && self.absorbing_barrier_kt.is_some(),
            "SteadyStateAbsorbingBarrier",
        );
        note(
            EIGEN_SOLVER,
            !(olzmann || cse) && self.eigen_solver.is_some(),
            "SteadyStateOlzmann and CSE",
        );
        note(
            SUM_RULE_TOLERANCE,
            !olzmann && self.sum_rule_tolerance.is_some(),
            "SteadyStateOlzmann (thermal eigenpair)",
        );
        for (name, given) in [
            (INTEGRATOR, self.integrator.is_some()),
            (INITIAL_STATE, self.initial_state.is_some()),
            (TIME_RANGE, self.time_range_s.is_some()),
            (TIMES_PER_DECADE, self.times_per_decade.is_some()),
            (INTEGRATION_TOLERANCE, self.integration_tolerance.is_some()),
        ] {
            note(name, !time && given, "TimeIntegration");
        }

        let solution = match method {
            SolutionMethod::SteadyStateOlzmann => {
                Solution::SteadyStateOlzmann(ThermalEigenSettings {
                    eigen_solver: self.eigen_solver.unwrap_or_default(),
                    sum_rule_tolerance: self
                        .sum_rule_tolerance
                        .unwrap_or(DEFAULT_SUM_RULE_TOLERANCE),
                })
            }
            SolutionMethod::SteadyStateAbsorbingBarrier => Solution::SteadyStateAbsorbingBarrier {
                absorbing_barrier_kt: self
                    .absorbing_barrier_kt
                    .unwrap_or(DEFAULT_ABSORBING_BARRIER_KT),
            },
            SolutionMethod::ChemicallySignificantEigenvalues => {
                // All eigenpairs are needed (G13 eqs. 21, 25-30): LAPACK unless the in-house decomposition is asked for.
                let eigen_solver = match self.eigen_solver {
                    None => EigenSolver::FullDecompositionLapack,
                    Some(EigenSolver::InverseIteration) => {
                        return Err(format!(
                            "Method CSE needs a full decomposition (all eigenpairs), which inverse iteration does not \
                             give: set {EIGEN_SOLVER} to Lapack or FullDecomposition."
                        ))
                    }
                    Some(solver) => solver,
                };
                Solution::ChemicallySignificantEigenvalues { eigen_solver }
            }
            SolutionMethod::TimeIntegration => Solution::TimeIntegration(TimeIntegrationPlan {
                integrator: self.integrator.unwrap_or(RosenbrockMethod::Rodas4),
                initial_state: self.initial_state.unwrap_or(InitialState::Pulse),
                time_range_s: self.time_range_s.unwrap_or(DEFAULT_TIME_RANGE_S),
                times_per_decade: self.times_per_decade.unwrap_or(DEFAULT_TIMES_PER_DECADE),
                relative_tolerance: self
                    .integration_tolerance
                    .unwrap_or(DEFAULT_INTEGRATION_TOLERANCE),
                absolute_tolerance: INTEGRATION_ABSOLUTE_TOLERANCE,
            }),
        };
        Ok(ResolvedSolution { solution, unused_settings })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn with(method: SolutionMethod) -> SolutionSettings {
        SolutionSettings {
            method: Some(method),
            ..Default::default()
        }
    }

    #[test]
    fn the_four_methods_have_deck_and_command_line_keywords() {
        for (keywords, method) in [
            (
                ["SteadyStateOlzmann", "steady-state-olzmann"],
                SolutionMethod::SteadyStateOlzmann,
            ),
            (
                [
                    "SteadyStateAbsorbingBarrier",
                    "steady-state-absorbing-barrier",
                ],
                SolutionMethod::SteadyStateAbsorbingBarrier,
            ),
            (
                ["CSE", "cse"],
                SolutionMethod::ChemicallySignificantEigenvalues,
            ),
            (
                ["TimeIntegration", "time-integration"],
                SolutionMethod::TimeIntegration,
            ),
        ] {
            for k in keywords {
                assert_eq!(SolutionMethod::from_keyword(k).unwrap(), method, "{k}");
            }
        }
        assert_eq!(
            SolutionMethod::from_keyword("ChemicallySignificantEigenvalues").unwrap(),
            SolutionMethod::ChemicallySignificantEigenvalues
        );
        assert_eq!(
            eigen_solver_from_keyword("inverse").unwrap(),
            EigenSolver::InverseIteration
        );
        assert_eq!(
            eigen_solver_from_keyword("FullDecomposition").unwrap(),
            EigenSolver::FullDecomposition
        );
        assert_eq!(
            eigen_solver_from_keyword("Lapack").unwrap(),
            EigenSolver::FullDecompositionLapack
        );
        assert!(eigen_solver_from_keyword("jacobi").is_err());
    }

    #[test]
    fn the_steady_state_family_needs_its_version_and_the_eigenvalue_analysis_is_no_method() {
        for k in ["SteadyState", "steady-state"] {
            let e = SolutionMethod::from_keyword(k).unwrap_err();
            assert!(
                e.contains("SteadyStateOlzmann") && e.contains("SteadyStateAbsorbingBarrier"),
                "{e}"
            );
        }
        // The thermal eigenpair belongs to the final steady state (GO10 eq. 12).
        assert!(SolutionMethod::from_keyword("eigenvalue")
            .unwrap_err()
            .contains("SteadyStateOlzmann"));
        assert!(SolutionMethod::from_keyword("both").is_err());
    }

    #[test]
    fn a_method_must_be_given() {
        let e = SolutionSettings::default().resolve().unwrap_err();
        for name in [
            "SteadyStateOlzmann",
            "SteadyStateAbsorbingBarrier",
            "CSE",
            "TimeIntegration",
        ] {
            assert!(e.contains(name), "{e}");
        }
    }

    #[test]
    fn each_steady_state_solver_has_its_own_settings() {
        assert_eq!(
            with(SolutionMethod::SteadyStateOlzmann)
                .resolve()
                .unwrap()
                .solution,
            Solution::SteadyStateOlzmann(ThermalEigenSettings {
                eigen_solver: EigenSolver::InverseIteration,
                sum_rule_tolerance: DEFAULT_SUM_RULE_TOLERANCE,
            })
        );
        assert_eq!(
            with(SolutionMethod::SteadyStateAbsorbingBarrier)
                .resolve()
                .unwrap()
                .solution,
            Solution::SteadyStateAbsorbingBarrier {
                absorbing_barrier_kt: DEFAULT_ABSORBING_BARRIER_KT
            }
        );
        let five = SolutionSettings {
            absorbing_barrier_kt: Some(5.0),
            ..with(SolutionMethod::SteadyStateAbsorbingBarrier)
        };
        assert_eq!(
            five.resolve().unwrap().solution,
            Solution::SteadyStateAbsorbingBarrier {
                absorbing_barrier_kt: 5.0
            }
        );
    }

    #[test]
    fn the_cse_method_uses_lapack_by_default_and_refuses_inverse_iteration() {
        let cse = with(SolutionMethod::ChemicallySignificantEigenvalues);
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
    fn settings_the_chosen_method_does_not_use_are_listed_not_refused() {
        // A deck that lists the steady-state keywords can be switched to CSE.
        let template = SolutionSettings {
            absorbing_barrier_kt: Some(10.0),
            eigen_solver: Some(EigenSolver::FullDecompositionLapack),
            sum_rule_tolerance: Some(1.5e-2),
            ..with(SolutionMethod::ChemicallySignificantEigenvalues)
        };
        let unused = template.resolve().unwrap().unused_settings;
        assert_eq!(unused.len(), 2, "{unused:?}");
        assert!(unused
            .iter()
            .any(|n| n.contains("AbsorbingBarrierBelowThreshold[kT] (--barrier-kt)")));
        assert!(unused
            .iter()
            .any(|n| n.contains("SumRuleTolerance (--sum-rule-tolerance)")));
        // The absorbing barrier belongs only to SteadyStateAbsorbingBarrier, the eigenpair only to SteadyStateOlzmann.
        let olzmann = SolutionSettings {
            absorbing_barrier_kt: Some(5.0),
            ..with(SolutionMethod::SteadyStateOlzmann)
        };
        assert_eq!(olzmann.resolve().unwrap().unused_settings.len(), 1);
        let barrier = SolutionSettings {
            eigen_solver: Some(EigenSolver::FullDecomposition),
            sum_rule_tolerance: Some(0.01),
            ..with(SolutionMethod::SteadyStateAbsorbingBarrier)
        };
        assert_eq!(barrier.resolve().unwrap().unused_settings.len(), 2);
    }

    #[test]
    fn the_number_of_cores_is_a_run_setting_of_every_method() {
        // NCores (deck) is overridden by --ncore (command line); it is never an unused setting, because
        // every method runs its conditions on these cores.
        let deck = SolutionSettings { cores: Some(8), ..with(SolutionMethod::ChemicallySignificantEigenvalues) };
        assert_eq!(deck.overridden_by(&SolutionSettings::default()).cores, Some(8));
        let command_line = SolutionSettings { cores: Some(2), ..Default::default() };
        assert_eq!(deck.overridden_by(&command_line).cores, Some(2));
        for method in [
            SolutionMethod::SteadyStateOlzmann,
            SolutionMethod::SteadyStateAbsorbingBarrier,
            SolutionMethod::ChemicallySignificantEigenvalues,
            SolutionMethod::TimeIntegration,
        ] {
            let settings = SolutionSettings { cores: Some(8), ..with(method) };
            assert!(settings.resolve().unwrap().unused_settings.is_empty(), "{method:?}");
        }
        let zero = SolutionSettings { cores: Some(0), ..with(SolutionMethod::SteadyStateOlzmann) };
        assert!(zero.resolve().unwrap_err().contains("NCores"));
    }

    #[test]
    fn non_positive_numbers_are_refused() {
        let barrier = SolutionSettings {
            absorbing_barrier_kt: Some(0.0),
            ..with(SolutionMethod::SteadyStateAbsorbingBarrier)
        };
        assert!(barrier
            .resolve()
            .unwrap_err()
            .contains("AbsorbingBarrierBelowThreshold[kT]"));
        let tolerance = SolutionSettings {
            sum_rule_tolerance: Some(f64::NAN),
            ..with(SolutionMethod::SteadyStateOlzmann)
        };
        assert!(tolerance
            .resolve()
            .unwrap_err()
            .contains("SumRuleTolerance"));
    }

    #[test]
    fn time_integration_is_a_method_with_its_own_settings() {
        assert_eq!(
            integrator_from_keyword("Rodas4").unwrap(),
            RosenbrockMethod::Rodas4
        );
        assert_eq!(
            integrator_from_keyword("ros2").unwrap(),
            RosenbrockMethod::Ros2
        );
        assert!(integrator_from_keyword("Rang3").is_err());
        assert_eq!(
            initial_state_from_keyword("Pulse").unwrap(),
            InitialState::Pulse
        );
        assert_eq!(
            initial_state_from_keyword("continuous").unwrap(),
            InitialState::ContinuousFormation
        );
        assert_eq!(
            with(SolutionMethod::TimeIntegration)
                .resolve()
                .unwrap()
                .solution,
            Solution::TimeIntegration(TimeIntegrationPlan {
                integrator: RosenbrockMethod::Rodas4,
                initial_state: InitialState::Pulse,
                time_range_s: DEFAULT_TIME_RANGE_S,
                times_per_decade: DEFAULT_TIMES_PER_DECADE,
                relative_tolerance: DEFAULT_INTEGRATION_TOLERANCE,
                absolute_tolerance: INTEGRATION_ABSOLUTE_TOLERANCE,
            })
        );
    }

    #[test]
    fn time_integration_settings_are_checked_and_noted_when_unused() {
        let ti = with(SolutionMethod::TimeIntegration);
        assert!(SolutionSettings {
            time_range_s: Some((1.0, 0.5)),
            ..ti
        }
        .resolve()
        .unwrap_err()
        .contains("TimeRange"));
        assert!(SolutionSettings {
            time_range_s: Some((0.0, 1.0)),
            ..ti
        }
        .resolve()
        .unwrap_err()
        .contains("TimeRange"));
        assert!(SolutionSettings {
            times_per_decade: Some(0),
            ..ti
        }
        .resolve()
        .unwrap_err()
        .contains("TimesPerDecade"));
        assert!(SolutionSettings {
            integration_tolerance: Some(2.0),
            ..ti
        }
        .resolve()
        .unwrap_err()
        .contains("IntegrationTolerance"));
        let unused = SolutionSettings {
            absorbing_barrier_kt: Some(5.0),
            ..ti
        }
        .resolve()
        .unwrap()
        .unused_settings;
        assert_eq!(unused.len(), 1, "{unused:?}");
        let olzmann = SolutionSettings {
            integrator: Some(RosenbrockMethod::Ros4),
            times_per_decade: Some(3),
            ..with(SolutionMethod::SteadyStateOlzmann)
        };
        assert_eq!(olzmann.resolve().unwrap().unused_settings.len(), 2);
    }

    #[test]
    fn the_command_line_overrides_the_deck() {
        let deck = SolutionSettings {
            method: Some(SolutionMethod::SteadyStateOlzmann),
            sum_rule_tolerance: Some(0.02),
            ..Default::default()
        };
        let command_line = SolutionSettings {
            method: Some(SolutionMethod::SteadyStateAbsorbingBarrier),
            absorbing_barrier_kt: Some(5.0),
            ..Default::default()
        };
        assert_eq!(
            deck.overridden_by(&command_line),
            SolutionSettings {
                method: Some(SolutionMethod::SteadyStateAbsorbingBarrier),
                absorbing_barrier_kt: Some(5.0),
                sum_rule_tolerance: Some(0.02),
                ..Default::default()
            }
        );
    }
}
