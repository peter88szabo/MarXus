# Sum rules and other checks of the validity of a master-equation calculation in MarXus

**MarXus, 2026-10-05.** Peter asked: "write in the report on this sum-rules and other test to check the validity of the master equation".

This report lists every check that tells whether a master-equation (ME) result can be trusted. For each check it gives:
- what is checked and why it must hold, with the literature;
- where MarXus enforces it: a **runtime check**, applied to every calculation, which turns a violation into an error or (sum rule, §4.1) a warning, or a **unit test** (`cargo test`, test name given);
- what the C₂H₃ benchmark shows (`validation/c2h3_mess_example/`, `validation/c2h3_mess_example_olzmann_eigen/` (merged into `validation/c2h3_mess_example/` on 2026-10-06, its Section 4.5)).

## References

| abbreviation | reference |
|---|---|
| O91 | Olzmann, Gebhardt, Scherzer, Int. J. Chem. Kinet. 23, 825 (1991) |
| O02 | Olzmann, Phys. Chem. Chem. Phys. 4, 3614 (2002) |
| GO10 | González-García, Olzmann, Phys. Chem. Chem. Phys. 12, 12290 (2010) |
| PO14 | Pfeifle, Olzmann, Int. J. Chem. Kinet. 46, 231 (2014) |
| R19 | Robertson (ed.), Comprehensive Chemical Kinetics 43 (2019) |
| PR03 | Pilling, Robertson, Annu. Rev. Phys. Chem. 54, 245 (2003) |
| SN84 | Schranz, Nordholm, Chem. Phys. 87, 163 (1984) |
| DGP86 | Davies, Green, Pilling, Chem. Phys. Lett. 126, 373 (1986) |
| M79 | Miller, J. Am. Chem. Soc. 101, 6810 (1979) |
| W12 | Weyl, Math. Ann. 71, 441 (1912) |
| H02 | Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed. (SIAM, 2002), ch. 10 |
| NR92 | Press, Teukolsky, Vetterling, Flannery, *Numerical Recipes in Fortran*, 2nd ed. (1992), §11.2–11.3, §11.7 |

**Notation.** The ME is J·N = R·F with J = ω(I − P) + K + k_c[D]·I (PO14 eq. 2). Here:
- P is the collisional transition-probability matrix;
- K holds the specific rate coefficients (products, isomerization);
- k_c[D] is a pseudo-first-order bimolecular sink;
- f = ρ(E)·e^{−E/kT} is the Boltzmann weight;
- S = D⁻¹JD with D = diag(√f) is the symmetrized operator (R19 eqs. 5.74–5.77);
- λ₁ < λ₂ ≤ … are the eigenvalues of J (equal to those of S).

## 1. The collision operator

**1.1 Completeness (normalization) of P.** Σ_t P(t|j) = 1 for every grain j that is kept in the ME.
- **Why.** Collisions neither create nor destroy molecules. GO10 p. 12296: the final steady state is obtained "by carefully observing the completeness of transition probabilities". A kernel that is not normalized acts as a spurious sink or source. That would change λ₁ and break the sum rule of §3.1.
- **MarXus.** The exponential-down kernel is normalized column by column after its upward part has been fixed by detailed balance. The stepladder is O91 eqs. 16–18.
- **Tests:**
  - `collision_kernels::exponential_down_is_normalized_and_detailed_balanced_for_sparse_low_energy_states`;
  - `collision_kernels::exponential_down_reduces_to_unscaled_back_substitution_for_smooth_densities`;
  - `collision_kernels::stepladder_follows_olzmann_eqs_16_to_18`.

**1.2 Detailed balance of P.** P(t|j)·f_j = P(j|t)·f_t.
- **Why.** At equilibrium, collisions alone must leave the Boltzmann distribution unchanged (stepladder with detailed balance: O91 eqs. 13–18).
- It is also what makes J symmetrizable (§2.2). Symmetrizability is the basis of the Cholesky solver and of the eigenvalue analysis.
- **Tests:** the kernel tests above and `chemical_activation_operator::j_is_detailed_balanced_with_boltzmann_weights_on_the_absolute_energy_scale`.

## 2. The operator J

**2.1 Column sums of J are the loss rates (mass conservation).** Σ_r J_{r,c} = Σ_{product channels} k(E_c) + k_c[D] for every state c.
- **Why.** Collisions (§1.1) and isomerizations move population inside the network. Only product channels and the sink remove it. This identity is the basis of the sum rules of §3.1 and §4.1.
- **Test:** `chemical_activation_operator::column_sums_of_j_equal_the_losses_out_of_the_network`.

**2.2 Symmetrizability.** S = D⁻¹JD must be symmetric. That requires detailed balance of P (§1.2) and of every isomerization pair: k_{w→w'}(E)·ρ_w(E) = k_{w'→w}(E)·ρ_{w'}(E) at the same absolute energy (W‡/h in both directions).
- **Runtime check.** max |S_rc − S_cr| / max(|S_rc|, |S_cr|) is computed for every operator. The banded Cholesky solver and the eigenvalue analysis **refuse** an asymmetry above 10⁻⁸ (`SYMMETRY_TOLERANCE`), with the message that the isomerization rates violate detailed balance.
- **Tests:**
  - `chemical_activation_steady_state::cholesky_refuses_an_operator_without_detailed_balance`;
  - `chemical_activation_operator::isomerization_detailed_balance_check_flags_inconsistent_reverse_rates`;
  - `chemical_activation_from_mess_input::isomerization_rates_obey_detailed_balance_exactly`.

**2.3 Positive definiteness of S.** The eigenvalues of J are positive ("λᵢ > 0", GO10 text before eq. 12) whenever population can leave the network. S is then symmetric positive definite and has a Cholesky factor.
- **Runtime check.** A failed Cholesky factorization ("non-positive pivot") is an error, never a silent result.
- **Eigenvalue analysis.** The inverse iteration factors S + σI with the safety shift σ = n·ε·max Sᵢᵢ (inverse iteration with a shift, NR92 §11.7; Peter's go-ahead 20:18).
  - S + σI has the same eigenvectors as S. Its factor exists whenever the smallest eigenvalue relative to the diagonal exceeds a multiple of n·ε (H02, ch. 10; Demmel, LAPACK Working Note 14 (1989)).
  - With the shift, C₂H₃ at 300 K is solved at all pressures, where the plain factor did not exist at 3 of 5.
  - If the shifted factorization still fails, the message says to use a full decomposition (`--eigen-solver lapack`).
- **λ₁ ≤ 0** is a **warning** (§4.1): rounding noise below the double-precision floor, and the thermal eigenvector, with it k_uni, is still resolved.
- **Tests:**
  - `symmetric_eigen::a_shift_leaves_the_eigenpairs_unchanged`, `a_shift_lets_inverse_iteration_find_the_null_vector_of_a_singular_matrix`, `the_safety_shift_scales_with_dimension_precision_and_the_largest_diagonal_element`;
  - `chemical_activation_eigen::the_shifted_factorization_gives_k_uni_where_the_plain_cholesky_factor_does_not_exist` (125 K deep well, agrees with QL to 10⁻⁵), `the_shift_does_not_change_a_resolved_result`.
- **Why it can fail although the theorem holds.** When λ₁/max(S) is below the double-precision resolution, S is *numerically* indefinite. Example: C₂H₃ at 300 K, where λ₁ ≈ 10⁻¹⁵ s⁻¹, while the diagonal of S contains ω (Lennard-Jones collision frequency, roughly 10⁸–10¹⁰ s⁻¹ over 76–7600 Torr) and k(E) up to the top of the grid. §4.1 explains how this is detected even when the factorization happens to succeed.
- **Tests:** `banded_solvers::banded_cholesky_rejects_an_indefinite_matrix`, `banded_solvers::banded_cholesky_matches_dense_cholesky`.

## 3. Steady-state solutions (absorbing barrier, final steady state with a sink)

**3.1 Mass balance: the yields add up to one.** For a normalized source, Σ Φ_products + Σ Φ_stab + Σ Φ_sink = 1 (O02 eq. 10).
- **Why.** Every molecule formed leaves through exactly one exit.
- **Runtime check.** The driver refuses a condition with |Σ Φ − 1| > tolerance (10⁻⁸ in the example).
- **Test:** `chemical_activation_observables::yields_of_products_stabilization_and_sink_add_up_to_one`.

**3.2 Residual.** ‖F − J·N‖₂/‖F‖₂ is evaluated with the *unsymmetrized* J, independently of the solver.
- **Runtime check.** The driver refuses a residual above the tolerance (10⁻⁸ in the example). In the final steady state without a sink, a large residual means J is numerically singular. The message points to O02 (the final steady state is reached only after about 0.1/k_uni).
- **Test:** `chemical_activation_steady_state::steady_state_satisfies_j_n_equals_f_for_both_solvers_and_all_options`.

**3.3 Two independent solvers agree.** Banded Cholesky on S, and Jacobi-preconditioned BiCGSTAB, which does not use symmetry.
- **Test:** `chemical_activation_steady_state::banded_cholesky_and_bicgstab_give_the_same_populations`.

**3.4 Positivity of the populations.** N ≥ 0 for a non-negative source. J is an M-matrix (positive diagonal, non-positive off-diagonal elements), so J⁻¹ ≥ 0.
- **Test:** the assertion "negative population" inside `steady_state_satisfies_j_n_equals_f_for_both_solvers_and_all_options`.
- **Not yet a runtime check** (open item, §6).

**3.5 Limiting cases** (all unit tests in `chemical_activation_observables`):

| limit | expected (literature) | test |
|---|---|---|
| zero pressure | yields = RRKM branching ratios of the nascent distribution, Φ_r = Σ k_r F/Σ k F | `zero_pressure_yields_are_the_rrkm_branching_of_the_nascent_distribution` |
| infinite pressure, intermediate steady state | everything is stabilized | `high_pressure_intermediate_steady_state_stabilizes_everything` |
| final steady state without a sink | everything ends in products, "trivially R1 = R2 and hence Φ2 = 1" (O02 after eq. 13) | `final_steady_state_without_a_sink_ends_entirely_in_products` |
| equilibrium: final steady state fed thermally through the only channel | N is the Boltzmann distribution, and k^ca is the canonical high-pressure rate coefficient (TST) | `final_steady_state_through_the_only_channel_is_the_equilibrium_distribution` |
| k^ca as an average | k_r^ca = Σ k_r Ñ over the normalized well distribution (GO10 eq. 9) | `rate_coefficients_average_k_over_the_normalized_well_distribution` |

**3.6 Absorbing-barrier plateau.** In the intermediate steady state, the result must not depend on where the absorbing barrier is placed. O02 Fig. 2 shows the plateau of the yield over many decades of the sink rate. The barrier must lie well below the reactive region and well above the thermal distribution (PR03; Carstensen, Dean, Comprehensive Chemical Kinetics 42 (2007)).
- **Runtime check.** A barrier below the well bottom is an error. Test: `chemical_activation_operator::an_absorbing_barrier_below_the_well_bottom_is_an_error`.
- **User control.** The distance is selectable (`--barrier-kt`). Test: `a_smaller_barrier_distance_moves_the_absorbing_barrier_up`.
- **C₂H₃.** `validation/c2h3_mess_example/barrier_distance_sensitivity.csv`. Largest spread of the association between barriers 10, 5 and 3 kT below the threshold:

| T (K) | largest spread |
|---|---|
| 300 | 0.55% |
| 500 | 1.0% |
| 750 | 2.2% |
| 1000 | 4.1% |
| 1250 | 7.4% |
| 1500 | 17% |

  - The plateau is therefore well developed only at ≤ 500 K. It degrades from 750 K and is lost at ≥ 1500 K, where the 10 kT barrier lies inside the thermal distribution of the well: −13% at 1500 K / 0.1 atm, −42% at 1750 K.
  - **A missing plateau means the intermediate-steady-state result is not defined.** The eigenvalue analysis (§4) must then be used. At 1000 K its association agrees with the 10 kT result within 0.02% (§4.6).

**3.7 O02 validity window of a physical sink.** The final steady state with a sink reproduces the absorbing-barrier branching when 0.01·ω > k_c[D] > 10·k_uni (O02, discussion of Fig. 2).
- **MarXus:** `chemical_activation_eigen::steady_state_window`.
- **Test:** `lambda_f_belongs_to_the_spectrum_and_defines_the_window`.
- **Not yet printed** by the driver for each condition (open item, §6).

## 4. Eigenvalue analysis (Olzmann's solution method)

### 4.1 The GO10 eq. 12 sum rule and the reported k_uni (runtime check, warning)

**Statement.** Let Ñ be the thermal eigenvector of λ₁, normalized to unit population. Then, from §2.1,

  λ₁ = 1ᵀJÑ = Σ_j k_j^th + k_c[D] ≡ k_uni,  with k_j^th = Σ_i k_j(E_i)·Ñ_i.

The right-hand side is GO10's "averaging procedure analogous to eqn (9) but with Ñs = Ñs^th being the normalized eigenvector associated with the lowest eigenvalue λ₁" (GO10 after eq. 12). Both forms are GO10's own.

**Peter's decisions (2026-10-05, 20:01):**
- MarXus reports k_uni, the eigenvector average, as the unimolecular rate coefficient, and explains this in the output with the citation.
- The sum rule never rejects a result; it only warns.
- The warning threshold is 1.5%.

**Why the two sides behave differently in double precision.**
- **k_uni is a sum of positive terms**, free of cancellation.
- **λ₁ from any eigen-solver** is an eigenvalue of S + δS, with ‖δS‖ of the order of ε‖S‖:
  - the backward error of the Cholesky factorization (H02, ch. 10);
  - or that of the Householder/QL reduction (Wilkinson, *The Algebraic Eigenvalue Problem*, 1965).

  By Weyl's inequality (W12) it therefore carries an *absolute* error of order ε‖S‖ ≈ 10⁻¹⁶·max(ω + k(E)).
- **The eigenvector is much less sensitive.** Its perturbation is of order ε‖S‖/(λ₂ − λ₁), small because λ₂ ≫ λ₁ in a deep well. This is why k_uni is the reported quantity.
- **So the relative difference measures resolution.** |λ₁ − k_uni|/k_uni grows like ε‖S‖/λ₁ and shows directly whether the eigenpair is resolved in double precision.

**Runtime behaviour.**
- **Fields.** Every result carries `k_uni_s_inv`, `lambda_1_s_inv`, `sum_rule_relative_deviation` and `warning`.
- **Table.** The thermal table has the columns `k_uni[1/s]`, `lambda_1[1/s]`, `sum_rule_deviation`, `lambda_2/k_uni`. Before the header, comment lines explain k_uni and λ₁ with the GO10 citation (`THERMAL_TABLE_EXPLANATION` in `chemical_activation_driver.rs`).
- **Warning.** A deviation above the tolerance (`DEFAULT_SUM_RULE_TOLERANCE` = 1.5·10⁻², option `--sum-rule-tolerance`) produces a `# warning` line before the table, also printed on stderr by the example. It gives:
  - the numbers;
  - the explanation with the citation;
  - **what to do**: use k_uni and confirm it with the inverse iteration, which is the most accurate for the thermal eigenvector; if the solvers give different k_uni, the condition is beyond double precision.
- **λ₁ ≤ 0** (possible with the full decompositions) gets the same warning plus the remark that a non-positive λ₁ is impossible for J ("all eigenvalues positive", GO10 before eq. 12), so it is numerical noise.
- **Precision floor.** ε·max Sᵢᵢ is printed for every condition (column `precision_floor`, field `precision_floor_s_inv`). The warning for λ₁ ≤ 0 (Peter, 20:28) says that λ₁ is **rounding noise below the double-precision floor** (with the floor in s⁻¹), and that this is **not a merging of eigenvalues** (λ₂/k_uni printed). Merging, i.e. thermal decay no longer separated from relaxation at high T, shows as a small λ₂/k_uni, not as a negative λ₁. The table explanation says the same.
- **Errors remain only for:**
  - k_uni ≤ 0 (nothing leaves the network);
  - a Cholesky factorization that fails even with the shift (did not occur in any test or benchmark).

**Tests** (`chemical_activation_eigen`):
- `channel_rates_and_sink_add_up_to_lambda_1`: two wells with isomerization and a sink, deviation < 10⁻⁸.
- `k_uni_is_the_eigenvector_average_and_the_sum_rule_deviation_is_reported`.
- `a_violated_sum_rule_gives_a_warning_and_no_error`: deep well (threshold 6000 cm⁻¹, exponential down), 200 K, 10 Torr. The deviation is about 6·10⁻⁵: no warning at 1.5·10⁻², a warning at 10⁻⁵, with the same k_uni in both cases.
- `a_non_positive_lambda_1_is_a_warning_when_the_eigenvector_average_is_positive`: 150 K, QL decomposition. The warning names the precision floor and "not a merging of eigenvalues"; k_uni agrees with the shifted inverse iteration to 10⁻⁵.
- `a_network_without_losses_is_an_error`.
- Driver tests: `the_thermal_route_gives_falloff_rate_coefficients_at_every_condition` (explanation with citation, columns) and `sum_rule_warnings_are_written_before_the_thermal_table`.

**What the probe and the C₂H₃ benchmark show.** The deep test well has 800 grains and its threshold at grain 600.

**Caveat on this probe:** it used the stepladder model, whose step of 10–14 grains splits J into independent sub-equations (§4.7). The deviations below therefore also contain sub-equation mixing. The unit tests now use exponential down. With exponential down and the shifted inverse iteration, k_uni = 1.61·10⁻¹⁴, 6.13·10⁻¹⁰, 9.47·10⁻⁷ and 2.02·10⁻⁴ s⁻¹ at 125, 150, 175 and 200 K (10 Torr), and the sum-rule deviation is 1.3·10⁶, 66, 1.9·10⁻² and 6.4·10⁻⁵.

In the probe (stepladder, inverse iteration without shift), the deviation grows smoothly as T falls:

| T (K) | sum-rule deviation (1 and 10 Torr) |
|---|---|
| 300 | 5e-11 … 9e-11 |
| 250 | 1e-8 … 3e-8 |
| 225 | 1e-7 … 4e-7 |
| 200 | 1e-5 … 4e-5 |
| 175 | 1.4e-3 … 1.7e-3 |
| 150 | Cholesky factor does not exist |

- **The QL and LAPACK solvers on the same well.** At 150 K both give λ₁ < 0. At 175 K LAPACK gives λ₁ = 8.06e-7 and k_uni = 9.67e-7 s⁻¹, while inverse iteration gives 8.79e-7 and 8.80e-7 s⁻¹. **The dense solvers' eigenvector can therefore also degrade**; hence the advice to confirm with the inverse iteration.
- **Before the check existed,** a run at 400 K and 1 Torr of an even deeper well (threshold 15000 cm⁻¹) returned λ₁ = 7.30e-9 s⁻¹ without comment, against an eigenvector average of 7.02e-9 s⁻¹ (4%). It now reports k_uni = 7.02e-9 s⁻¹ with a warning.

**C₂H₃, full deck** (`validation/c2h3_mess_example_olzmann_eigen/` (merged into `validation/c2h3_mess_example/` on 2026-10-06, its Section 4.5)):
- **≥ 1000 K:** deviation 10⁻¹⁴–10⁻¹¹ (inverse iteration).
- **750 K:** 4·10⁻⁹–7·10⁻⁸.
- **500 K:** 4·10⁻⁴–1.0·10⁻²; no warning.
- **300 K, inverse iteration with the shift.** All five pressures are solved. λ₁ is ±10⁻⁸ … 10⁻⁶ s⁻¹, rounding noise 13 orders below the floor 1.0·10⁻² s⁻¹ (warnings; λ₁ < 0 at three pressures). k_uni is 5.98, 7.49, 8.74, 9.43, 9.81·10⁻¹⁶ s⁻¹.
- **300 K, LAPACK.** k_uni agrees with the inverse iteration to 10⁻⁷ … 10⁻⁴.
- **k_uni is physically right.** Against MESS it is +5.66, +5.35, +5.06, +4.86, +4.73%. The association by detailed balance equals the independent absorbing-barrier route to 0.001 percentage points (+5.828 vs +5.828% at 0.1 atm).

### 4.2 High-pressure limit and fall-off

**Statement.** For ω → ∞ the thermal eigenvector becomes the Boltzmann distribution and λ₁ → k∞ = Σ k f/Σ f. At finite pressure λ₁ < k∞, and λ₁ increases monotonically with pressure.
- **Tests:**
  - `chemical_activation_eigen::at_high_pressure_lambda_1_is_the_boltzmann_average_of_k` (agreement 10⁻³ at 10⁹ Torr; < 0.9 k∞ at 1 Torr);
  - `chemical_activation_driver::the_thermal_route_gives_falloff_rate_coefficients_at_every_condition`.
- **C₂H₃.** k∞ of both directions agrees with MESS within +0.6 … +2.9%, the exact vs semiclassical Eckart difference. Their ratio, the equilibrium constant, agrees within 0.04–0.10%.

### 4.3 Eigenpair residual and solver agreement

**Statement.** J·E₁ = λ₁·E₁ for the returned eigenvector. Inverse iteration, Householder/QL and LAPACK DSYEVD must give the same λ₁ and λ₂.
- **Tests:**
  - `both_solvers_give_the_lowest_eigenpairs_of_j` (λ₁ to 10⁻⁸, λ₂ to 10⁻⁶, residual 10⁻⁸);
  - `the_lapack_decomposition_agrees_with_the_other_solvers`;
  - `symmetric_eigen::*` (5 tests, including `inverse_iteration_resolves_an_eigenvalue_many_orders_below_the_others`);
  - `lapack_interface::divide_and_conquer_agrees_with_the_householder_ql_decomposition`.
- **C₂H₃, short decks (1000 K, 1 atm).** All three solvers give λ₁ = 16004.04 s⁻¹ to 7 digits.

### 4.4 Separation of time scales, λ₂/λ₁

**Statement.** A single phenomenological rate coefficient k_uni (= λ₁) is meaningful only if the thermal decay is much slower than relaxation, λ₂ ≫ λ₁.

The intermediate steady state exists for (0.1·λ_F)⁻¹ < t < (10·k_uni)⁻¹ (SN84; O02, discussion of Fig. 2). λ_F is the eigenvalue "that corresponds to the eigenvector most closely resembling the initial distribution f(E)" (O02).
- **MarXus.** λ₂/λ₁ is printed for every condition (NaN when the deflated inverse iteration cannot resolve λ₂, at separations beyond about 10¹⁶). λ_F and the window come from `EigenSystem::lambda_f` and `steady_state_window`.
- **Test:** `lambda_f_belongs_to_the_spectrum_and_defines_the_window`.
- **C₂H₃ (1 atm):**

| T (K) | λ₂/λ₁ |
|---|---|
| 500 | 3·10¹² |
| 1000 | 1.35·10⁴ |
| 1500 | 72 |
| 2000 | 11.6 |

- **What happens as the separation is lost.** Where λ₂/λ₁ → 10, MESS's pressure-dependent association/dissociation pair departs from MESS's own equilibrium constant: −1.3% at 1500 K, −11% at 2000 K and 1 atm. MarXus imposes k_assoc = λ₁·K. At 2000 K the two codes then differ in opposite directions for the two reactions: dissociation −4.5 … −7.8%, association +3.0 … +7.2%. Small λ₂/λ₁ therefore flags conditions where the phenomenological rate coefficient itself becomes ill defined.

### 4.5 Time-dependent solution, limits

**Statement.** N(t) = R Σᵢ eᵢ (1 − e^{−λᵢt})/λᵢ·Eᵢ (PO14 eqs. 3–4) must satisfy:
- N(t → ∞) = R·J⁻¹F, the steady state;
- N(t → 0) = R·F·t.
- **Test:** `chemical_activation_eigen::time_dependent_population_approaches_the_steady_state`. It requires a full decomposition (`the_time_dependent_solution_needs_a_full_decomposition`).

### 4.6 Agreement of the two independent solution routes

**Statement.** Where both are valid, i.e. where the absorbing-barrier plateau exists (§3.6), they must agree:
- the intermediate steady state with an absorbing barrier, k(R→W) = k∞·Φ_stab (PR03 eq. 44);
- the eigenvalue analysis by detailed balance, k(R→W) = k_uni·k∞,assoc/k∞,diss.

**C₂H₃.** Association deviation from MESS, eigenvalue (k_uni, shifted inverse iteration) vs barrier route (10 kT):

| T (K) | eigenvalue route | barrier route |
|---|---|---|
| 300 | +5.83 / +5.52 / +5.22 / +5.03 / +4.89% | +5.83 / +5.52 / +5.22 / +5.03 / +4.90% |
| 500 | +3.76 … +3.22% | +3.77 … +3.22% |
| 750 | +1.45 … +1.86% | +1.46 … +1.86% |
| 1000 | −0.60 … +0.54% | −0.58 … +0.54% |

- **The routes agree to 0.01% at every (T, p).** With λ₁ instead of k_uni, the difference at 500 K / 1 atm had been 0.34%, exactly the sum-rule deviation of λ₁ there.
- **This confirms that k_uni is the accurate side of the sum rule**, including at 300 K where λ₁ is off by 10⁸.
- **At ≥ 1250 K the routes separate** because the barrier route loses its plateau.

**This is the strongest end-to-end check:** two different formulations of the same ME agree.

### 4.7 Stepladder with a step of several grains: reducible J

**The problem.** With the stepladder model, a collision changes the energy by exactly ±ΔE_SL. If ΔE_SL spans several grains, J splits into ΔE_SL/ΔE independent sub-equations, energetically shifted by one grain. O02 p. 3616 does exactly this on purpose: "19 subequations, which are energetically shifted with respect to each other by 10 cm⁻¹ and which are solved separately … the averaging of the rate coefficients becomes much more precise."

**Consequences for the eigenvalue analysis.**
- Every sub-equation has its own lowest eigenvalue, and these are nearly degenerate.
- "The lowest eigenvalue of J" then belongs to one sub-equation only.
- Inverse iteration separates them only slowly, and not at all once the shift σ exceeds their differences: at ≤ 175 K in the deep test well it did not converge.
- How O02 combines the sub-equations' λ₁ into k₂ᵘⁿⁱ is not stated in the paper.

**Peter's decision (20:25):** exponential down, which keeps J irreducible, is now the default collision model (`CollisionModel::default()`, cutoff 15⟨ΔE_down⟩; test `exponential_down_is_the_default_collision_model`).

**Open:** the treatment of the stepladder in the eigenvalue analysis (per sub-equation, then averaged?) awaits the literature or Peter.

## 5. Input-level checks (the ME is only as good as k(E), ρ(E) and the graining)

| check | literature | MarXus test |
|---|---|---|
| k∞ of the entrance from the ME input reproduces the ILT input | DGP86 | `chemical_activation_from_mess_input::entrance_high_pressure_rate_reproduces_the_ilt_input`; `ilt_barrierless::association_inversion_reproduces_the_recombination_rate_coefficient`, `dissociation_inversion_reproduces_k_inf_for_a_fractional_exponent` |
| k∞ of a tight TS is canonical TST | R19 | `entrance_high_pressure_rate_of_a_tight_transition_state_is_transition_state_theory`; C₂H₃: 1 cm⁻¹ cells → 70 cm⁻¹ grains reproduce TST within 0.08% (`tunneling_ilt_and_energy_graining.md`) |
| graining conserves the cell sums; ρ, W as cell averages | PR03 p. 254 | `energy_graining::grain_averages_conserve_the_cell_sum`, `grain_densities_are_cell_averages` |
| tunneling: the same canonical correction in both directions (detailed balance of k(E) with tunneling) | M79 eqs. 8–9 | `tunneling::canonical_eckart_correction_is_the_same_in_both_directions`, `tunneling_sum_of_states_is_the_stieltjes_convolution_of_miller_eq_9` |
| tunneling below the TS stops at the higher ground state of the two sides | M79 | `tunneling_stops_at_the_highest_ground_state_of_the_two_sides` |

## 6. Open items

1. **Positivity as a runtime check.** N ≥ 0 (§3.4) and a one-signed thermal eigenvector E₁ are guaranteed by the M-matrix structure of J (Perron–Frobenius for the lowest eigenvector of S, whose off-diagonal elements are non-positive). They are tested but not checked at runtime. A sign change in E₁ would be a further, independent signal of lost precision.
2. **O02 window per condition.** The sink window and the intermediate-steady-state time window (§3.7, §4.4) are computed by library functions but not yet printed by the driver.
3. **Double precision for very deep wells.**
   - **Done:** k_uni from the eigenvector average (Peter's decision), and the shifted factorization (Peter's go-ahead).
   - **Remaining:** λ₁ itself (and with it the sum rule as a check) is lost below the floor ε·max Sᵢᵢ.
   - **Planned (Peter, 20:28–20:40):** higher precision.
     - Double-double first (qd + faer), with an MPFR reference path at 192 and 256 bits.
     - The sensitive steps (kernel normalization, detailed balance, diagonal loss sums, symmetrization) must be assembled in that precision; casting a rounded f64 matrix cannot recover what is lost.
     - Symmetric solvers throughout.
     - Convergence of the slow eigenvalues and of the rate coefficients to be checked as the precision increases.
   - **Alternative, not implemented:** a cancellation-free elimination: Grassmann, Taksar, Heyman, Oper. Res. 33, 1107 (1985).
4. **Convergence with respect to the grid** (grain width, top of the grid, `ExcessEnergyOverTemperature`). It was checked for the high-pressure limit (0.08%) but not yet systematically for k_uni and the yields.
5. **Phenomenological rate coefficients at small λ₂/λ₁** (§4.4). To be settled from the literature (Georgievskii et al., J. Phys. Chem. A 117, 12146 (2013), in `papers/ChemAct/`) before any change.
