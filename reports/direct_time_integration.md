# Direct time integration of the master equation: the third solution method

**MarXus, 2026-10-06.** Peter's request: "implement a third method: direct time integration, with an adaptive, L-stable Rosenbrock method as the first choice. You would integrate the actual energy-grained populations, including the early relaxation and subsequent chemistry. This requires neither a steady-state assumption nor a separation between chemical and relaxation modes. Check the KPP software … we gonna adapt their stiff-diff integrators: ROS4 Rosenbrock and others … make a new directory in the numerical modules for the integrators, and a new master equation solver where the other master equation solvers are, and make a new direct_time_integration module."

## 1. Sources

**KPP, the Kinetic PreProcessor** (`/home/peter/Programs/KPP`, int/rosenbrock.f90):
- (C) Adrian Sandu, 2004; revised by P. Miehe and A. Sandu, 2006.
- GNU GPL v3, the license of MarXus, so the adaptation is license-compatible; it is attributed in the source.
- Method references, from the KPP bibliography (docs/source/citations/kpp.bib):

| method | reference |
|---|---|
| Ros2 | Verwer, Spee, Blom, Hundsdorfer, SIAM J. Sci. Comput. 20, 1456 (1999) |
| Ros3, Rodas3 | Sandu, Verwer, Blom, Spee, Carmichael, Potra, Atmos. Environ. 31, 3459 (1997) |
| Ros4, Rodas4 | Hairer, Wanner, Solving ODEs II, Springer (1991; 2nd ed. 1996), Sec. IV.7 |

**Not adopted: KPP's Rang3**, the W-method of Rang and Angermann, BIT 45, 761 (2005).
- **The method itself** is of order 3, as checked: the local error for y' = −y falls with h⁴.
- **Its embedded error estimate vanishes** for a linear problem. For y' = −2y and h = 0.1, 0.5, 2 the estimate is 6·10⁻¹⁶, 1·10⁻¹⁵, 1·10⁻¹⁷, while the true local errors are 3·10⁻⁵, 6·10⁻³, 9·10⁻².
- **Consequence:** the step-size control accepts any step, and the integration of y' = −2y to t = 3 returned −0.033 instead of 0.0025.

## 2. Code

**`src/numeric/integrators/`** (new directory; `mod.rs` holds only `pub mod` lines):

**`rosenbrock_methods.rs`.** `RosenbrockMethod {Ros2, Ros3, Ros4, Rodas3, Rodas4}` with `tableau()`, KPP's coefficients in its row-wise storage of A and C. The method form (KPP, H&W IV.7):
- G = 1/(hγ₁) − ∂f/∂y;
- G K_i = f(t + α_i h, y + Σ A_ij K_j) + Σ (C_ij/h) K_j + hγ_i ∂f/∂t;
- y₁ = y + Σ M_i K_i, with the error estimate Σ E_i K_i.

**`rosenbrock.rs`.** `integrate(system, y, t0, t1, options)`, a port of KPP's `ros_Integrator`:
- f, ∂f/∂t (finite differences unless autonomous) and the preparation of G. A failed factorization halves h, at most 5 times.
- The stages, with stage function values reused where KPP's `ros_NewF` is false.
- The scaled RMS error norm (at least 10⁻¹⁰).
- The step control $h_\mathrm{new} = h\cdot\min(6, \max(0.2,\ 0.9/\mathrm{err}^{1/\mathrm{ELO}}))$; acceptance if err ≤ 1 or h ≤ h_min; no growth after a rejection; ×0.1 after two rejections.
- Tolerance checks as in KPP.
- `StiffSystem` trait: `dimension`, `rhs`, `prepare(t, y, shift)` (factorize shift·I − ∂f/∂y), `solve`.
- `RosenbrockOptions` holds KPP's defaults. One MarXus extension, off by default: `power_of_two_steps` rounds every new step down to a power of two. Steps are then never larger than KPP's, and step sizes repeat, so a constant-Jacobian system can reuse factorizations.

**`src/masterequation/direct_time_integration.rs`:**
- **State:** the populations N of all retained grains and the yield Y_x of every exit: product channels `W->X`, sinks `escape(W)`, and with the absorbing-barrier operator `stab(W)`.
- **Equations:** dN/dt = R F − J N, dY_x/dt = Σ k_x(E) N(E). The source formed directly below an absorbing barrier goes into stab(W).
- **Stage solve.** G is block lower triangular:
  - $(s\,I + J)\,x_N = b_N$ via $x_N = D\,(s\,I + S)^{-1}D^{-1}b_N$, with the banded Cholesky factor of $s\,I + S$ and $D^{-1}b$ formed from logarithms, as in the steady-state solver;
  - $x_Y = (b_Y + K^T x_N)/s$.
- **Cache:** the 4 most recent factors; the number actually computed is reported.
- **`InitialState`:** `Pulse` (N(0) = F, R = 0) or `ContinuousFormation` (N(0) = 0, R = 1 s⁻¹). `log_spaced_times`.
- **`TimeEvolution`:** exits, points (t, well populations, exit yields), final grain populations, statistics, factorizations computed.

**Settings** (`solution_method.rs`):
- `Method TimeIntegration` (`--method time-integration`);
- `Integrator` (Rodas4);
- `InitialState` (Pulse);
- `TimeRange[s] t₁ t₂` (10⁻¹² to 10² s);
- `TimesPerDecade` (4);
- `IntegrationTolerance` (10⁻⁶; absolute 10⁻¹⁴).

The steady-state and CSE settings are noted as unused, and the other way round. The deck block reads all of these keywords (`mess_input.rs`).

**Program** (`chemical_activation_from_deck`):
- Conditions in parallel (`ConditionPool`); the source is the thermal entrance source of each T; the operator is that of the final steady state.
- **Report:** RUN SETTINGS give the method, initial state, time range, integrator and tolerances. Per condition, a table of t against N(W) and the exit yields in %, with the integrator statistics. The groups, in three views: yields at the last time, populations left, yields without the return to R, and the total.
- **Machine-readable:** blocks `# time evolution: T = … K, p = … Torr` in the `--csv` file; the groups in `_tables.csv`.

## 3. Tests (TDD; each seen failing on a stub)

**`numeric::integrators`** (8 tests):

| test | checks |
|---|---|
| `every_tableau_has_consistent_dimensions`, `row_wise_storage_maps_stage_indices` | the transcribed tables |
| `exponential_decay_is_reproduced_by_every_method` | y' = −2y to t = 3, 10⁻⁵ relative |
| `a_stiff_non_autonomous_problem_takes_large_steps` | Prothero–Robinson, λ = 10⁶: y = sin t to 10⁻⁵, < 20 000 steps (Ros2 needs about 6000; an explicit method about 10⁷); tests ∂f/∂t |
| `the_robertson_problem_reaches_its_reference_values` | at t = 40: 0.7158270687, 9.185534764·10⁻⁶, 0.2841637457 to 10⁻⁵; y₁+y₂+y₃ = 1 to 10⁻¹² |
| `every_method_has_its_nominal_order` | y' = −y³ with fixed steps |
| `power_of_two_steps_keep_the_accuracy_and_repeat_the_step_sizes` | Robertson with the extension |
| `unreasonable_tolerances_are_refused` | the tolerance checks |

Observed orders with fixed steps: Ros2 1.88, Ros3 3.00, Ros4 3.86, Rodas3 2.94, Rodas4 4.51. Rodas4's 4.5 comes from higher-order terms at these steps.

**Note on the order test.** It first used y' = −y². Rodas3 integrates that problem exactly, to rounding (a property of this Riccati equation with Rodas3's coefficients), so the ratio of rounding errors suggested "order 1". The test was changed to y' = −y³.

**`masterequation::direct_time_integration`** (6 tests: two wells, entrance, sink; 300 K, 760 Torr):

| test | checks |
|---|---|
| pulse | total N + Y = 1 to 10⁻⁹ at every time; at long times the yields equal the final-steady-state yields to 10⁻⁶; the factorizations are reused (fewer than half the steps) |
| continuous formation | N(t) approaches $N^s = J^{-1}F$ to 10⁻⁶ |
| N(t) | equals the eigenvector expansion (PO14 eqs. 3–4) at 10⁻¹¹ … 10⁻⁵ s to 10⁻⁶ |
| thermalized well | the late decay rate equals λ₁ = k_uni to 10⁻⁴ |
| absorbing barrier | a pulse ends in the stabilization yields of the intermediate steady state |
| `log_spaced_times` | the output times |

**Settings and deck:** 2 tests in `solution_method`, 1 in `mess_input`. **Report:** 1 test in `report_sections`.

**Efficiency.** Before the power-of-two steps and the cache, the six ME tests took 58 s, with about 300 factorizations per integration of the test network. Afterwards they took 14 s; on the four-well network, 120 factorizations for 895 steps.

## 4. Validation: ZZ-allyl + O₂, four wells (`validation/ZZAllyl+O2_Gamma_Case2/`)

**Run.** Pulse, Rodas4, 10⁻¹² … 10² s, 4 per decade, MESS Eckart model, 21 conditions on 4 threads: 102 s.

**Identity.** At t = 100 s the yields of R, P1, P5, P7 and escape equal the final steady state at all 21 conditions, to the 7 printed digits (`time_integration_vs_final_steady_state.csv`). The total stays 100.0%.

**Time scales at 300 K, 760 Torr** (`plots/time_evolution_300K_760torr.png`):
- the nascent G2 redissociates to R within about 10⁻⁹ s (78.4% in the end);
- IEPOX + OH (P5) is complete at about 10⁻⁸ s (0.502%): chemically activated G4;
- escape from G4 builds up between 10⁻¹⁰ and 10⁻⁵ s (21.1%);
- P7 has a prompt part near 10⁻⁹ s and a thermal part near 10⁻⁴ … 10⁻³ s, through stabilized G3 → G6;
- G3 lives longest (until about 10⁻² s).

**Integrator work at 300 K, 760 Torr:** 895 steps, 0 rejected, 5368 function evaluations, 120 factorizations.

## 5. Open points

1. **Cost per factorization.** On multiwell networks the bandwidth of S is set by the isomerization couplings in the well-by-well state order. Sorting the states by energy would reduce it for Case2 from about 840 to about 320 (estimated; for the two-well test network the kernel alone gives 344).
2. **Extensions** for which the time integration is the natural method: time-dependent sources or temperature, and nonlinear (second-order) chemistry. The integrator is general (`StiffSystem`); the ME system is linear today.
