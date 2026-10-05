# The two solution methods of MarXus and their settings in the deck header

**MarXus, 2026-10-05 (evening).** Peter's questions and decisions:
- "why CSE runs with --steady-state? is not it a different method? we either use steady-state method and a linear system, or we use Klippenstein eigenvalue method?"
- "Olzmann eigenvalue analysis is also a step after the steady-state solver not? So we have two methods! eigenvalue à la Klippenstein (cse) and steady-state method with two flavors: Olzmann or the intermediate absorbing-barrier (as in SSUMES)"
- "in the MESS type input file in the header we need also a section where we define these … what type of solver we request"

**My error.** The example program selected the CSE method and the thermal eigenvalue analysis through `--steady-state cse|eigenvalue|all`, which mixed up the method with a version of one method. This should have been recognized when both were built. The thermal eigenpair is computed from the same J as the final steady state, by the same banded Cholesky factor. Peter was right.

## 1. Literature basis

**GO10.** G. Gonzalez-Garcia, J. Olzmann, Phys. Chem. Chem. Phys. 12, 12290 (2010), checked in `/home/peter/Dropbox/Unimolecular_theory/Chemical_Activation/Olzmann_c0cp00284d.pdf`:
- **eq. 7:** the steady-state master equation J·N = F, written for dn(E)/dt = 0. **Eq. 8:** its normalized solution Ñs = J⁻¹F / Σᵢ(J⁻¹F)ᵢ. **Eq. 9:** k^ca = Σ (K_r Ñs).
- **Before eq. 12:** "It is important to note that also the rate coefficient of the thermal (superscript: th) unimolecular decomposition of HSO5 can be obtained from eqn (7) as the lowest eigenvalue λ₁ of the matrix J for [H2O] = 0". **Eq. 12:** k^th = λ₁.
- **After eq. 12:** "Alternatively, the rate coefficient k_i^th (T,P) can also be calculated by an averaging procedure analogous to eqn (9) but with Ñs = Ñs^th being the normalized eigenvector associated with the lowest eigenvalue λ₁." The eigenvector is the **thermal steady-state population** ñs^th(E).
- **Sec. 3.2:** the final ("asymptotic") steady state and the "so called intermediate steady state, which is characterized by the competition between decomposition and stabilization". "In technical terms the intermediate steady state can be implemented by introducing a lower absorbing barrier into the master equation. By omitting this absorbing barrier … the final steady-state solution for ñs(E) is obtained from eqn (8)."

**The CSE method.** MK06: J. A. Miller, S. J. Klippenstein, J. Phys. Chem. A 110, 10528 (2006). G13: Y. Georgievskii, J. A. Miller, M. P. Burke, S. J. Klippenstein, J. Phys. Chem. A 117, 12146 (2013). See `chemically_significant_eigenvalues_method.md`.

## 2. The two methods

| method | equation | versions / parts | code |
|---|---|---|---|
| **1. Steady state** | J·N = F (GO10 eqs. 7, 8) | **intermediate** (absorbing barrier X k_BT below the lowest threshold, as in SSUMES); **final** (no barrier, Olzmann). The final one includes the **thermal rate coefficients** from the lowest eigenpair of the same J (GO10 eq. 12): k_uni = eigenvector average, λ₁ beside it, sum rule | `run_chemical_activation`, `SteadyState::{Intermediate, Final}`; `run_thermal_rate_coefficients`, `chemical_activation_eigen.rs` |
| **2. CSE** | all eigenpairs of J; phenomenological rate coefficients (G13 eqs. 15, 21–30) | none | `run_phenomenological_rates`, `chemically_significant_eigenvalues.rs` |

**Why the thermal eigenpair belongs to the final steady state.**
- **Same operator.** It uses the operator of the final steady state, J without an absorbing barrier (`SteadyState::Final`).
- **Same equation, without the source.** Its eigenvector is the thermal steady-state population, the solution of the same equation with the thermal source replaced by the population itself.
- **The solver shows it.** MarXus's default solver, inverse iteration with the banded Cholesky factor of S + σI, is a sequence of steady-state solves. Each step solves (S + σI)x = u with the previous distribution as the source.
- **Not CSE.** The method is not the CSE method: it needs only the lowest eigenpair and no phenomenological rate-coefficient matrix.

## 3. Settings: the `MarXus` block of the deck header

The new library file is `src/masterequation/solution_method.rs`.

```
MarXus
  Method                              SteadyState        ! SteadyState | CSE
  SteadyState                         Both               ! Intermediate | Final | Both
  AbsorbingBarrierBelowThreshold[kT]  10                 ! intermediate steady state
  EigenSolver                         InverseIteration   ! InverseIteration | FullDecomposition | Lapack
  SumRuleTolerance                    1.5e-2             ! thermal eigenpair of the final steady state
End
```

**Reading the block.** It is a MarXus extension of the MESS format, read by `parse_marxus_header` in `mess_input.rs`, like the `InverseLaplaceTransform` block of a barrier. The deck stores it in `MessGlobal::solution`, a `SolutionSettings` in which every field is optional.

**Errors in the block:**
- an unknown keyword;
- a keyword given twice;
- the block given twice.

**Spelling of values.** Case, '-' and '_' are ignored, so `SteadyState`, `steady-state` and `steady_state` are the same. Short forms are also accepted: `inverse`, `full`, `lapack`, `cse`, and `ChemicallySignificantEigenvalues`.

**Command line.** It overrides the deck, setting by setting (`SolutionSettings::overridden_by`):

| deck keyword | option | default |
|---|---|---|
| `Method` | `--method steady-state\|cse` | SteadyState |
| `SteadyState` | `--steady-state intermediate\|final\|both` | Both |
| `AbsorbingBarrierBelowThreshold[kT]` | `--barrier-kt X` | 10 |
| `EigenSolver` | `--eigen-solver inverse\|full\|lapack` | InverseIteration for the final steady state; Lapack for CSE |
| `SumRuleTolerance` | `--sum-rule-tolerance X` | 1.5e-2 |

**Resolution** (`SolutionSettings::resolve` → `ResolvedSolution { solution: Solution, unused_settings }`):
- **Errors:**
  - `eigenvalue` as a method or as a steady-state version. The message explains that the thermal eigenpair is part of the final steady state (GO10 eq. 12).
  - `cse` as a steady-state version. The message points to `Method CSE` / `--method cse`.
  - CSE with inverse iteration, because CSE needs all eigenpairs.
  - A non-positive or non-finite barrier distance or tolerance.
- **Notes, not errors:** a setting that the selected solution does not use. Examples:
  - the absorbing barrier with `SteadyState Final`;
  - `SumRuleTolerance` with CSE.

  The note is printed as `# note:` in the output and as `note:` on stderr. A template deck that lists every keyword can therefore be switched to the other method from the command line.
- **Why notes and not errors.** This follows the warn-don't-reject line used for the sum rule. A real conflict (CSE with inverse iteration) remains an error.

## 4. Example program `chemical_activation_from_deck`

**Settings.** The command-line settings are combined with the deck's: `deck.global.solution.overridden_by(&command_line).resolve()`.

**Output header.** A new line `# solution method: …` names the resolved method, its versions and solvers. It is followed by any `# note:` lines.

**Order of the output:**
1. intermediate steady state;
2. final steady state;
3. thermal rate coefficients of the final steady state, with the association by detailed balance for one well and one entrance;
4. or, instead of 1–3, the CSE tables.

**Renamed titles** (the values in the tables are unchanged):
- `# eigenvalue analysis (…) of J without absorbing barrier; …` → `# thermal rate coefficients of the final steady state: lowest eigenpair of J (GO10 eq. 12) by …; sum-rule tolerance …`
- `… by detailed balance, eigenvalue analysis` → `… by detailed balance from the thermal k_uni of the final steady state`

**Removed options:** `--steady-state eigenvalue|cse|all`. Each is refused with an explanation.

## 5. Tests (TDD: each written first and seen failing on a stub)

| test | checks |
|---|---|
| `solution_method::keywords_of_the_deck_and_of_the_command_line_are_both_accepted` | deck and command-line spellings of every value |
| `…::the_eigenvalue_analysis_and_cse_are_not_versions_of_the_steady_state` | `eigenvalue`, `cse` and `all` are refused, with the explanation |
| `…::the_default_is_the_steady_state_method_in_both_versions` | defaults: both versions, 10 kT, inverse iteration, 1.5e-2 |
| `…::the_cse_method_uses_lapack_by_default_and_refuses_inverse_iteration` | the CSE solver |
| `…::settings_the_selected_solution_does_not_use_are_listed_not_refused` | a template deck switched to CSE; final only; intermediate only |
| `…::non_positive_numbers_are_refused` | 0 kT, NaN tolerance |
| `…::the_command_line_overrides_the_deck` | merging setting by setting |
| `mess_input::the_marxus_header_block_gives_the_solution_settings` | the block is read; the header lists around it are still read |
| `mess_input::a_deck_without_a_marxus_header_block_leaves_the_solution_settings_unset` | no block: nothing is set (a guard; it passes on the stub by design) |
| `mess_input::the_marxus_header_block_refuses_unknown_keywords_values_and_repetitions` | unknown keyword, the value `Eigenvalue`, repeated keyword, two blocks |

**End-to-end checks of the example**, run on `c2h3_tight_short.inp`:
- The four refused option combinations give the explanatory errors.
- A copy of the deck with a `MarXus` block (`SteadyState Final`, `EigenSolver Lapack`) runs the final steady state with its thermal eigenpair by LAPACK.
- The same deck with `--method cse` runs CSE, with a note that `SteadyState` is not used.

## 6. Validation re-runs (at most 4 cores)

The comparison CSVs were saved before the re-run and compared byte by byte afterwards.

- **`validation/ZZAllyl+O2_Gamma_Case2/`.**
  - `run_marxus.sh` no longer has the two `--steady-state eigenvalue` runs. The thermal block now sits in the `*_steady_states.out` files.
  - `compare_with_mess.py` reads that block.
  - The CSE run uses `--method cse`.
  - **All six comparison CSVs are byte-identical** to the previous ones. The thermal data rows equal those of the former `case2_tstlevel_E_eigenvalue.out`.
  - The obsolete `case2_tstlevel_E_eigenvalue.out` (tracked in git) and `case2_tstlevel_E_mess_eckart_eigenvalue.out` were removed.
- **`validation/c2h3_mess_example_olzmann_eigen/`.** `run_marxus.sh` uses `--method steady-state --steady-state final`, and `plot_comparison.py` reads the renamed thermal title. Result: Re-run after the change: `comparison_table.csv` is byte-identical, and the thermal and association rows of all eight outputs equal the previous ones. The sum-rule warnings are unchanged: 6 for inverse iteration and 9 for LAPACK on the full deck. The final steady-state table now in the outputs is "not available" at 300 and 500 K on the full deck (`c2h3_tight_*`), because J·N = F is numerically singular for the deep well without a sink. The thermal eigenpair is unaffected. `plot_comparison.py` now reads the "not available" and warning lines of the thermal block only, so these lines do not mark the thermal results as missing.
- **`validation/c2h3_mess_example/`** (intermediate steady state, `--steady-state intermediate --barrier-kt X`). The options are unchanged; it was re-run for the new output header. Result: Re-run after the change: `comparison_table.csv` and `barrier_distance_sensitivity.csv` are byte-identical. The default runs (`both`) now also contain the thermal rate coefficients of the final steady state. `plot_comparison.py` recognizes that block: without it, the parser failed with a ValueError on the new table header.

## 7. Files changed

- **New:** `src/masterequation/solution_method.rs`, with a `pub mod` line in `masterequation/mod.rs`.
- **`src/masterequation/mess_input.rs`:** `MessGlobal::solution`, `parse_marxus_header`, `MarXus` as a block starter, 3 tests.
- **`examples/chemical_activation_from_deck.rs`:** options, dispatch on the two methods, output header, titles, module documentation.
- **README.md:** the two methods, the header block and an options table.
- **Validation READMEs and scripts** of the three directories above.
- **Reports:** `olzmann_absorbing_barrier_and_treatment.md`, `chemically_significant_eigenvalues_method.md` (option names).

## 8. Open points

- **Tunneling model in the block.** `--tunneling exact-eckart|mess-eckart` is still a command-line option only. It could become a keyword of the `MarXus` block, so that a deck fully specifies a run. Not done, because it is a model choice rather than a solver choice; this is for Peter to decide.

**Full test suite** (2026-10-05, 4 cores): 176 library tests (166 before, plus the 10 new ones) and 34 binary tests pass.

## 9. README: the solver comparison and the tunneling manual (Peter, 22:45)

Peter's requests:
- "we leave it as optional but not default for users who want more closer comparison with MESS (this should be also in the manual (readme))";
- "the two solvers are not the same and they answer two different questions (and also even the two different steady-state solvers …) so three in total";
- "emphasize that we solve the same master equation, but with different questions to answer";
- "for chemical activation describe what it is and how the different three solvers tell about it and what can provide";
- "give the formulas how the initial flux is prepared for chemical activation from reactants".

**New README section "Master-equation solvers: one master equation, three questions"** (LaTeX math in GitHub's `$$` syntax):
- **The common equation.** PO14 eq. 2; O02 eq. 6; GO10 eq. 7; G13 eq. 1 with Ĝ = J.
- **Chemical activation.** What it is, and the source formulas:
  - the thermal entrance source (PO14 eqs. 7, 9; O02 eq. 11);
  - the multi-entrance weighting on the absolute energy scale;
  - k∞ from the same W‡;
  - the CSE source and its normalization (G13 eqs. 1, 9, 23);
  - non-thermal and shift sources (PO14 eqs. 8, 10–13).
- **One part per solver,** giving its question, its equations, its time window (O02 p. 3618; SN84; GO10 p. 12295; MK06 eq. 19), what it gives and does not give for chemical activation, and its limits.
- **A side-by-side table.**
- **Identities and measured agreement** (C₂H₃ and Case2 validations).

**Checked against the code while writing:**
- The stabilization flux includes isomerization flux that arrives below the barrier of the target well (`chemical_activation_operator.rs`).
- The C₂H₃ barrier-route numbers are at 0.1 atm with 10 kT.

**New README subsection "Tunneling model (`--tunneling`)":**
- the exact Eckart formula (M79 eq. 8) as the default, with its microcanonical and canonical use;
- the MESS model, explicitly optional for comparison with MESS only, with its formula and the Case2 effect.

The Features list now says the same.

**Math rendering fix (Peter, 22:57).** The README formulas rendered with literal commas, e.g. "R , F − J , N".
- **Cause.** Markdown (CommonMark) treats a backslash before ASCII punctuation as an escape. `\,`, `\;` and `\!` therefore became plain `,`, `;` and `!` before the math engine saw them. Commands of backslash plus letters, such as `\frac` and `\mathbf`, are not affected.
- **Fix.** In all README math, `\,` was replaced by `\thinspace`, `\;` by `\thickspace`, and `\;\;` before equation labels by `\quad`; `\!` was removed. No backslash-punctuation is left in any math region (checked by script).

**Math rendering, second fix (Peter, 23:04: "all latex formulas were ugly").**
- **Renderer.** The README Peter saw was the pushed one (commit f91701f) on GitHub.
- **Cause.** GitHub's `$…$` and `$$…$$` math goes through Markdown processing first. This breaks the LaTeX in two ways:
  - backslash escapes turn `\,` into a comma;
  - underscores become emphasis, which mangles the subscripts.
- **Fix.** All 15 display formulas are now in GitHub's ```` ```math ```` fenced blocks, and all 97 inline formulas in GitHub's $\`…\`$ syntax. Both pass the LaTeX unprocessed to the math renderer, so the standard `\,` and `\;` were restored. A script checked that no `$` is left outside math and that no formula inside a table contains `|`.

**README opening paragraph (Peter, 23:05).** It now describes MarXus as a microcanonical rate code and master-equation solver, covering:
- k(E) from direct state counting, RRKM with exact Eckart, PST and ILT;
- the multiwell master equation and its two methods, which answer different questions of the same equation;
- the MESS-format input;
- the thermochemistry.

The status-table row "Olzmann eigenvalue analysis" was renamed to the thermal rate coefficients of the final steady state. The Features line now gives the master equation as dN/dt = R·F − J·N, with J·N = R·F as its steady state.
