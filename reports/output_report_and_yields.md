# Readable master-equation output, yields in percent, and the CSE / final-steady-state diagnostic

**MarXus, 2026-10-06.** Peter's requests (in order):
- "these agreement between the CSE and Final steady-state calculation should be shown in the report as an important diagnostic"
- "make a proper good looking output for the CSE and also the steady-state variants too … I want a nice looking one similar as in MESS"
- "not the same, but similar output style, temperature vs pressure table, another tables organized by pressure, another by temperature for each reactions"
- "for chemical activation we should report the yields too as T and P functions not just bare rates"
- "we need in the beginning of the core an energetics table in kcal/mol"
- "if there are multiple channels out of a potential well, and/or for reactant --> different bimolecular products, we have to report the yields too in % for different channels (chemical activation + thermal together and also separate)"
- "in the output also must be shown emphasized way what method and what details used for the solvers and run and chemical network"
- "in the output of the validations plot the yields too"
- "instead of the squares use for MESS points some large X"

## 1. The CSE / final-steady-state diagnostic

**Claim checked** (pasted by Peter): from the stored four-well outputs with the MESS tunneling option, the long-time shares reconstructed from MarXus's CSE tables agree with the final steady state.
- IEPOX + OH at 300 K, 760 Torr: 2.32525045% vs 2.32525079%.
- All 21 conditions: within 0.000032% relative.

**Recomputed here** from our own stored outputs (`validation/ZZAllyl+O2_Gamma_Case2/compare_with_mess.py`, section 7):
- The same two numbers.
- Largest relative deviation over the 21 conditions:

| channel | largest relative deviation |
|---|---|
| P5 | 3.14·10⁻⁷ |
| escape | 1.60·10⁻⁸ |
| P1 | 6.13·10⁻⁷ |
| P7 (share about 10⁻⁵) | 1.40·10⁻⁴ |

**Interpretation.** This is an exact identity, not an approximate agreement. With G13 eqs. 21 and 25–30, the absorbing-chain reconstruction gives $\sum_\lambda p^{(x)}_\lambda p^{(R)}_\lambda/(\Lambda_\lambda Q_R)$ over all eigenpairs. That is the spectral form of $k_x^T J^{-1} F$, the final-steady-state yield. The derivation is in `chemically_significant_eigenvalues_method.md`, Section 6.
- It holds independently of the eigenvalue separation.
- It validates the two code paths against each other: the banded Cholesky solve on one side; the eigendecomposition, $M^{-1}$, rate assembly and absorbing chain on the other. It does not validate the CSE approximation.
- The library test `cse_long_time_yields_equal_the_final_steady_state_yields` checks it to 10⁻⁸.

**Pulse identity** (Peter's note, verified). For a normalized pulse $F$, $Y_r(\infty) = \int_0^\infty k_r^T e^{-Jt}F\,dt = k_r^T J^{-1}F$, because the eigenvalues of $J$ are positive. This is the final-steady-state yield per formed adduct. It is in the README ("Where the three must agree").

**Distinction recorded in the README.** The lowest eigenpair belongs to the final-steady-state operator, but its eigenvector is the decay distribution of thermalized molecules without a source. It differs in general from the driven chemical-activation steady-state distribution.

**C₂H₃ correction.** The intermediate steady state (10 kT) and the thermal association of the final steady state agree:

| T | largest deviation |
|---|---|
| 300 K | 0.0007% |
| 500 K | 0.004% |
| 750 K | 0.018% |
| 1000 K | 0.023% |

Above that they separate: 1.2% at 1250 K, 13% at 1500 K, 77% at 1750 K. The previously quoted "0.01%" was corrected in the README, in `validation/c2h3_mess_example_olzmann_eigen/README.md` (two places) and in `olzmann_absorbing_barrier_and_treatment.md`.

## 2. Yields: definitions and code

**Steady state.** All yields are percent of the formed adducts unless stated otherwise.

| quantity | definition | code |
|---|---|---|
| prompt (chemically activated) yields | $\Phi_x$ of the intermediate steady state, with stabilization $\Phi_{\mathrm{stab},w}$ | existing |
| thermal fate of well $w$ | final steady state with the Boltzmann distribution $f^0_w$ of well $w$ as the source: $k_x^T J^{-1} f^0_w$ | `run_thermal_well_fates` (driver), `WellFateConditionResult` |
| all together | final steady state $\Phi_x$ (and without the return to the reactant: % of the net reaction) | existing |
| through the stabilized wells | $\sum_w \Phi_{\mathrm{stab},w}\cdot \mathrm{fate}_w(x)$ | `report_sections::chemical_activation_and_thermal_groups` |
| difference | prompt + through the wells − final steady state (percentage points): how well the separation into a prompt and a thermal stage holds (O02) | same |
| thermal branching | $k^{th}_r/k_{\mathrm{uni}}$ of the thermal eigenpair (products and sinks) | `report_sections::thermal_groups` |

**CSE** (`chemically_significant_eigenvalues::reactant_yields`).

| quantity | definition |
|---|---|
| prompt branching of R | $k_{R\to X}/\sum_{X\ne R} k_{R\to X}$ (wells = stabilization; bimolecular = direct, well-skipping) |
| thermal fate of each well | absorbing chain $B = (I-Q)^{-1}A$ of the CSE well rate coefficients |
| long-time yields | direct $k_{R\to x}$ + through the wells $\sum_i k_{R\to i}B_{ix}$, normalized to the eventual net reaction; equal to the final steady state (Section 1) |

**Case2 example (MESS Eckart model): IEPOX + OH, % of the formed adducts.**

| quantity | value |
|---|---|
| prompt | 0.42–0.70% |
| through the stabilized wells | 0.003–0.036 percentage points |
| final steady state | 0.42–0.72% |
| difference | 6·10⁻⁶ … 8·10⁻⁴ percentage points |

IEPOX + OH comes almost entirely from chemically activated G4, and the prompt + thermal picture holds closely.

## 3. The report (`report_tables.rs`, `report_sections.rs`, example `chemical_activation_from_deck`)

**Header:**
- **RUN SETTINGS:** deck; method and versions; solvers (banded Cholesky, eigen-solver, tolerances); temperatures and pressures; energy grid (grain width, 1 cm⁻¹ cells); collision model with its reference; tunneling model; source (entrance channels, $F \propto \rho k e^{-E/kT}$, PO14 eqs. 7, 9); notes.
- **CHEMICAL NETWORK:**
  - every well with grains, absolute grain range, $\langle\Delta E_{\mathrm{down}}\rangle(T)$, Lennard-Jones σ, ε, reduced mass and escape sink;
  - every channel with from/to and its rate model from the deck (rigid TS, PST and its TST level, ILT; Eckart tunneling), entrances marked.
- **Energetics in kcal/mol relative to the Reactant**, as in the MESS output:
  - wells: G, and D, the lowest barrier out of the well;
  - bimolecular species: G;
  - barriers: H = ZeroEnergy of the barrier in the deck, from/to, model.
  - The first version took the bimolecular asymptote for barrierless barriers. That gave B6P7 at −14.86 instead of MESS's −4.70, and was corrected to the deck ZeroEnergy.

**Views.** Every quantity group is written in three views (`write_groups`):
- tables by temperature (rows: pressure);
- tables by pressure (rows: temperature);
- one temperature–pressure table per quantity (P\T).

Groups wider than 8 quantities are split into several tables. Numbers have six significant digits ("1.04275e+06"); a missing condition is written as `***`.

**Steady-state sections:**
- intermediate and final steady state: yields (%), yields without the return to R (% of the net reaction), k_ca (1/s), and k(R → X) = k∞Φ_X (intermediate);
- thermal rate coefficients of the final steady state: k_uni, λ₁, channel rates; thermal branching (%); k∞ per channel; eigenpair diagnostics; association by detailed balance (single well);
- thermal fates of the wells (%);
- chemical activation and thermal reaction separately and together (Section 2).

**CSE section:**
- species-to-species tables per condition (MESS-like `From\To`), with the chemical eigenvalues, the lowest relaxation eigenvalue, the separation, the capture/return/net of R, and the warnings;
- rate coefficients from every species;
- prompt branching, thermal fates and long-time yields (%);
- diagnostics.

**Machine-readable output.** `--csv FILE` writes the previous machine format (CSV blocks with titles) to FILE. The validation scripts read these files. The new `.csv` files are identical to the previously committed `.out` files except for the solution-method line.

**Removed:** the `println!("Iterative diagonalization is done in N steps.")` in `src/numeric/jacobi_diag.rs`, a debug print of the library that appeared in the report. MarXus only; the sister codes are untouched.

## 4. Tests (TDD, each seen failing on a stub)

| test | checks |
|---|---|
| `report_tables::*` (7) | number format; tables by T and by p; T–P tables; missing values; splitting at 8 columns; long names; labelled tables |
| `report_sections::*` (8) | energetics relative to R; yields sum to 100%; missing conditions; thermal branching and well fates sum to 100%; prompt + via = sum and difference; CSE rate and yield groups and species tables; three views; network summary |
| `chemical_activation_driver::thermal_fates_at_high_pressure_follow_the_high_pressure_branching` | at 10⁶ Torr the thermal fate equals the ratio of the Boltzmann-averaged k(E) to 10⁻³; monotonic approach (1.2·10⁻⁴ at 10⁶, 1.2·10⁻⁵ at 10⁷ Torr) |
| `…::thermal_fates_are_given_for_every_well_and_condition` | order and mass balance |
| `…::cse_long_time_yields_equal_the_final_steady_state_yields` | the identity of Section 1 to 10⁻⁸ |

**Note on the high-pressure test.** It was first written at 10⁹ Torr. There $J\,N = F$ is numerically singular in double precision; the solver reports this correctly. The conditions were moved to 1000 K and 10⁶–10⁷ Torr.

## 5. Validation re-runs (4 cores)

- **ZZ-allyl + O₂ Case 2.**
  - All comparison CSVs are identical to the committed ones.
  - The four machine-readable `.csv` files equal the committed `.out` files except for the method line.
  - `cse_vs_final_steady_state.csv` and `plots/cse_vs_final_steady_state.png` are new.
- **C₂H₃ (both directories).**
  - All comparison CSVs are identical to the committed ones.
  - The 13 machine-readable `.csv` files equal the committed `.out` files except for the method line.
  - **Parser fix:** the default runs of `c2h3_mess_example` (`both`) now also write the thermal block. Its parser recognizes the block; before the fix it failed with a ValueError on the new header.
  - **Parser fix:** the eigen parser reads "not available" and warning lines only inside the thermal block. The final steady state is "not available" at 300 and 500 K on the full deck (J·N = F singular for the deep well without a sink), and those lines must not mark the thermal results as missing.

## 6. Plots

**MESS markers.** In all comparison plots MESS is now drawn as large "x" markers on lines and MarXus as small filled circles, instead of squares over open circles that hid each other.

**Yield plots.**
- **`validation/ZZAllyl+O2_Gamma_Case2/plots/yields.png`** (read from the `*_tables.csv` files):
  - long-time yields of every channel at 760 Torr: MarXus final steady state (exact Eckart, circles; MESS Eckart, triangles) against MESS's long-time fate (x);
  - IEPOX + OH: prompt, through the stabilized wells, and together (final steady state; diamonds);
  - stabilization yields of the wells (intermediate steady state);
  - CSE prompt branching of R at 760 Torr.
- **`validation/c2h3_mess_example/plots/yields.png`:** stabilization and prompt-redissociation yields of the chemically activated C₂H₃ against pressure. MarXus (intermediate steady state) is compared with MESS's k(P1 → W1)/k∞.
- **`validation/ZZAllyl+O2_Gamma_Case2/plots/time_evolution_300K_760torr.png`:** the direct time integration (`direct_time_integration.md`).
- **Tables file.** With `--csv FILE`, the program writes every quantity group of the report to `FILE_tables.csv` (`report_tables::write_groups_csv`, test `groups_are_written_as_titled_csv_blocks`). Each block title carries its section, e.g. "final steady state: Yields (% of the formed adducts)".
