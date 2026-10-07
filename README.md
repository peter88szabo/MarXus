# MarXus

**Author:** Peter Szabo  
**Email:** peter88szabo@gmail.com  

MarXus is a **microcanonical rate code** and **master-equation solver** for gas-phase reaction kinetics, written in Rust. From the molecular data of the reactants, wells and transition states, it computes energy-resolved rate coefficients and solves the energy-grained master equation for pressure- and temperature-dependent kinetics, including chemical activation. The same molecular data also give the thermochemistry.

- **Microcanonical rate coefficients k(E).**
  - Sums and densities of states are counted directly on 1 cm⁻¹ cells (Beyer–Swinehart).
  - Tight transition states are treated by RRKM theory, with exact Eckart tunneling.
  - Barrierless channels are treated by phase space theory (TST levels T, E, EJ) or by the inverse Laplace transform of k∞(T).
  - SACM is in progress.
- **Master equation for multiwell networks.** It includes collisional energy transfer (exponential down or stepladder), isomerization, product channels, bimolecular sinks, and chemically activated formation from bimolecular reactants. It is solved by **four methods in three families**, which answer different questions of the same equation (one method per run; each has its own document):
  - steady state: [SteadyStateOlzmann](docs/methods/steady_state_olzmann.md), the final steady state, with the thermal rate coefficient from the lowest eigenpair of the same matrix;
  - steady state: [SteadyStateAbsorbingBarrier](docs/methods/steady_state_absorbing_barrier.md), the intermediate steady state;
  - eigenvalue: [CSE](docs/methods/chemically_significant_eigenvalues.md), phenomenological rate coefficients for kinetic models (Miller, Klippenstein; Georgievskii et al.);
  - time integration: [TimeIntegration](docs/methods/direct_time_integration.md), the grained populations integrated with adaptive, L-stable Rosenbrock methods (adapted from KPP), without a steady-state assumption or an eigenvalue separation.
- **Input** in the MESS deck format, so existing decks can be used, with MarXus extension blocks for its own settings.
- **Thermochemistry:** partition functions and U, H, F, G, S, Cv, Cp in the RRHO and quasi-RRHO (Grimme) approximations.

The master-equation, tunneling, ILT and numerical code cites the source of each equation (paper and equation number) in its comments. Design notes, derivations and validation results are in `reports/`; validations against reference calculations are in `validation/`.

---

## Current state (2026-10-07)

| area | status |
|---|---|
| Thermochemistry (RRHO, Grimme qRRHO) | implemented |
| Sum and density of states, RRKM k(E), canonical TST | implemented, tested |
| Tunneling (exact Eckart, microcanonical and canonical) | implemented, tested |
| Inverse Laplace transform (ILT) for barrierless channels | implemented, tested |
| Phase space theory (PST) for −C_n/Rⁿ potentials, levels T, E, EJ | implemented, tested, used in the master equation (validated, ZZ-allyl + O₂ Case 2) |
| PST with arbitrary 1D potential, SACM | in progress |
| Multiwell chemical-activation master equation: four methods in three families | |
| – [SteadyStateOlzmann](docs/methods/steady_state_olzmann.md) (final steady state; thermal eigenpair k_uni, λ₁, λ₂; thermal fates of the wells) | implemented, tested, validated (C₂H₃, 300–2000 K; four-well ZZ-allyl + O₂ Case 2) |
| – [SteadyStateAbsorbingBarrier](docs/methods/steady_state_absorbing_barrier.md) (intermediate steady state) | implemented, tested, validated (C₂H₃; Case 2) |
| – [CSE](docs/methods/chemically_significant_eigenvalues.md) (phenomenological rate coefficients, Miller–Klippenstein / Georgievskii et al. 2013) | implemented, tested, validated (four-well MESS run, Case 2; C₂H₃, 300–2000 K) |
| – [TimeIntegration](docs/methods/direct_time_integration.md) (Rosenbrock Ros2–Rodas4, adapted from KPP) | implemented, tested; reproduces the final steady state on the four-well network; late-time decay = k_uni within 4·10⁻⁶ (C₂H₃) |
| Parallel runs over the conditions (T, p) with rayon; cores from `NCores` (deck) or `--ncore` | implemented, tested |
| Exponential-down kernel with the low-energy reservoir state (MESMER) where the normalization fails | implemented, tested, validated on both systems; replaces the former reduction rule, whose integer window caused a step of +3% in R → G4 of Case 2 at 304.7 K (`reports/low_energy_reservoir_state.md`) |
| Lennard-Jones collision integral Ω(2,2)\*: Neufeld, Janzen, Aziz 1972 (default, as MESS and MESMER) or Troe 1977 eq. 3.3 | implemented, tested; with Neufeld, Case 2 agrees with MESS within 1.1% for every reactant rate coefficient (`reports/collision_integral_neufeld.md`) |
| Deck: `RotationalConstants[1/cm]` and `Mass[amu]` in place of a geometry; `Dummy` products (as MESS); ILT barriers without RRHO data | implemented, tested; used to translate MESMER inputs (`reports/deck_rotational_constants_and_dummy_species.md`) |
| Validation against MESMER (acetyl + O₂, two wells, ILT association, Eckart) | done: partition functions within 0.001%, rate coefficients −2.1 … +1.4% of MESMER's double-double run (`validation/acetyl_o2_mesmer/`) |
| Low-energy reservoir validity: warning below 3 k_BT between the reservoir top and the lowest threshold | implemented, tested (`reports/low_energy_reservoir_state.md`) |
| Higher precision (double-double, arbitrary-precision reference) | planned (`reports/higher_precision_decision.md`) |
| Excited electronic levels: partition functions and thermochemistry (`electronic_contributions`), and state counts ρ(E), W(E) of the master equation (all `ElectronicLevels` lines of the deck) | implemented, tested, validated against MESS (`validation/c2h4_ho2_hindered_rotors/`, variant with test levels: the change of every rate coefficient equal in both codes; `reports/electronic_levels_partition_functions.md`) |
| Hindered and free internal rotors (1D), MESS `Rotor Hindered` / `Rotor Free` syntax: quantum levels in a σ-periodic Fourier basis, Kilpatrick–Pitzer reduced moment (default) or bond-axis option (`RotorReducedMoment`, `--rotor-reduced-moment`) | implemented, tested, validated against MESS (`validation/c2h4_ho2_hindered_rotors/`, 7 rotors): rotor levels equal to every printed digit; with the MESS Eckart model, k∞ within 0.07–0.39% and the net reaction within +0.10 … +0.54% at all 30 conditions (`reports/hindered_rotors.md`, `reports/hindered_rotors_validation_mess.md`) |
| Dissociative photoionization: PEPICO breakdown curves by statistical, collision-free modelling (Sztáray, Bodi, Baer, J. Mass Spectrom. 45, 1233 (2010)), `src/photoion` | implemented, tested: thermal and ion energy distributions, fast and RRKM channels with kinetic shift, parallel channels, sequential fast dissociations with statistical product energy partitioning; PST and simplified SACM not available yet (`reports/photoion_statistical_modelling.md`) |
| Equilibrium constants: from the thermochemistry (K_p, K_c from ΔG°) and from the deck species (partition functions on the cells; real equilibrium constants as MESS) | implemented, tested; the two routes agree up to the 1 cm⁻¹ cell counting; against MESS (`validation/c2h4_ho2_hindered_rotors/`) K(W2/P1) is 0.9963–0.9995 of MESS, a difference fully explained by that counting, and the effect of electronic levels is identical (`reports/equilibrium_constants.md`) |
| CSE: rows of the bimolecular products (G13 eqs. 28, 21, 22) and the isomer–bimolecular equilibrium coefficients κ (eq. 34) | implemented, tested, validated against MESS: Case 2 product rows within −0.3 … +0.9% (P7 → G4 −1.4%), κ of the rotor deck within 0.006 (`reports/cse_kappa_and_product_rates.md`) |
| CSE species merging at poor time-scale separation (Georgievskii et al. 2013, Sec. IV; well partition as MESS) | implemented, tested (`reports/cse_species_merging.md`). Chemical eigenvectors by the relaxational projection 1 − F_ne ≤ `ChemicalEigenvalueMax`, as MESS's direct method (default), or by the eigenvalue ratio (`ChemicalSubspaceCriterion EigenvalueRatio`, `--chemical-subspace-criterion`); in `validation/c2h4_ho2_hindered_rotors/` both codes merge at the same 15 of 30 conditions |
| [Prepared experiments](docs/features/prepared_experiments.md) (`Preparation` block of the deck): non-thermal initial populations, time-profiled source channels and piecewise-constant bath histories | implemented, tested. TimeIntegration integrates all of it, with impulses as exact jumps and bath changes as events; the steady-state methods take the shape of the sources; CSE projects every source onto prompt products and post-relaxation species populations (Georgievskii et al. 2013, eqs. 13, 14, 24, 37–41). Diagnostics: flux coefficients r(t) and k^e (Barker, Frenklach, Golden 2015), incubation and vibrational relaxation times (Barker, King 1995), and the CSE description propagated in time against the time integration (`CompareWithCse`). C₂H₃ shock heating: τ_inc = 12.47 ns, converged to 0.09% in the grain width (`validation/c2h3_mess_example/`, Section 4.9) |
| [Fragment wells and lumped reactant states](docs/features/fragment_wells_and_lumped_reactants.md): a dissociation into a fragment well with a partner in excess and an energy-partitioning kernel (prior, modified prior, two-piece Gaussian; Green, Robertson 2014; Shannon et al. 2024); the bimolecular reactant as one thermal state (Miller et al. 2016, SI-VI) | implemented, tested: detailed balance for any kernel, equilibrium [C]/[B] = K_c[A] reached in time, the lumped state's initial loss rate equal to k∞(T)[X] within 7.5·10⁻⁴; OH + glyoxal validation in preparation |

The source contains 390 library unit tests (`cargo test`).

---

## Features

### Thermochemistry
- Molecules built from vibrational frequencies (with scaling), rotational constants, or Cartesian geometries; the moments of inertia are computed from the geometry with isotopic atomic masses (most abundant isotopes, AME2020: Wang et al., Chin. Phys. C 45, 030003 (2021)).
- Partition functions and thermodynamic functions **U, H, F, G, S, Cv, Cp** in the rigid-rotor–harmonic-oscillator approximation.
- **Quasi-RRHO entropy** with free-rotor interpolation for low-frequency modes (Grimme, Chem. Eur. J. 18, 9955 (2012)).
- **Equilibrium constants** $`K_p = e^{-\Delta G^\circ/RT}`$ and $`K_c = K_p (p^\circ/k_BT)^{\Delta n}`$ of any reaction from the standard Gibbs energies (`thermal::equilibrium::equilibrium_constant_from_thermochemistry`).
- **For deck species:** the partition functions from the master-equation cell densities and the real equilibrium constants of the wells and bimolecular species (report section "PARTITION FUNCTIONS AND EQUILIBRIUM CONSTANTS").

### State counting and microcanonical rate theory
- **Sum and density of states** by direct (Beyer–Swinehart) counting, with classical 1D, 2D and 3D rotors (Forst); rovibrational and bimolecular (convolved) states.
- **RRKM / microcanonical TST** specific rate coefficients k(E).
- **Canonical TST** high-pressure rate coefficients, with tunneling corrections.
- **Energy graining** on 1 cm⁻¹ cells averaged to grains (Pilling, Robertson, Annu. Rev. Phys. Chem. 54, 245 (2003)). It reproduces canonical TST within 0.08% for C₂H₃.

### Tunneling
- **Exact Eckart** transmission probability (Miller, J. Am. Chem. Soc. 101, 6810 (1979), eq. 8; Johnston, Heicklen 1962), overflow-free.
- **Microcanonical** tunneling sum of states by convolution with the transition-state sum of states (Miller 1979, eq. 9); the canonical correction κ(T) is the same in both directions.
- Wigner, Bell and Skodje–Truhlar corrections as optional models.
- **MESS Eckart model** (`tunneling/mess_eckart_tunneling.rs`, option `--tunneling mess-eckart`): **optional, not the default**, only for a closer comparison with MESS results (see "Tunneling model" under Usage).
  - It is MESS's semiclassical transmission 1/(1 + e^{−S(E)}) with the Eckart action, cutoff and canonical weight.
  - It reproduces the tunneling factors of a MESS log to 5·10⁻⁵, except below 300 K for κ ≫ 100.
  - MarXus uses the exact Eckart; the MESS model is not the recommended physics.
- `examples/eckart_kappa_from_deck.rs` tabulates κ(T) of every Eckart barrier of a deck, for both models.

### Barrierless reactions
- **Inverse Laplace transform (ILT)** of modified-Arrhenius high-pressure rate coefficients, for association and dissociation, with a grain-integrated kernel (Davies, Green, Pilling, Chem. Phys. Lett. 126, 373 (1986)).
- **Phase space theory for isotropic −C_n/Rⁿ potentials** (`barrierless/phasespace`), from two fragment geometries.
  - TST levels T (canonical variational), E (microcanonical, Georgievskii, Klippenstein, J. Chem. Phys. 122, 194103 (2005), eq. 58) and EJ (default, consistent with eq. 55); input keyword `TSTLevel`, as in MESS.
  - Used for barrierless channels in the master equation (entrance and exit).
- Phase space theory with an arbitrary 1D potential (Troe–Ushakov 2006 form): **in progress**.
- Statistical adiabatic channel model (SACM): **in progress**.

### Master equation: multiwell chemical activation

The energy-grained master equation is $`dN/dt = R\,F - 𝐉 N`$ with $`𝐉 = \omega(𝐈 - 𝐏) + 𝐊 + k_c[\mathrm D]\,𝐈`$ (Pfeifle, Olzmann, Int. J. Chem. Kinet. 46, 231 (2014), eq. 2); its steady state is $`𝐉 N = R\,F`$. See "Master-equation solvers" below for the questions each solver answers.

**Network**
- Any number of wells on a common absolute energy grid, connected by isomerization (exact detailed balance, k = W‡/(hρ) in both directions).
- Product channels, and pseudo-first-order bimolecular sinks k_c[D].

**Collisions**
- **Exponential down** (default), exactly normalized (Robertson (ed.), Comprehensive Chemical Kinetics 43 (2019), eq. 4.16).
  - **Low-energy reservoir.** Where the normalization fails at the sparse bottom of a well (Robertson 2019, p. 294), the grains from the failing one down form one thermalized reservoir state: the reservoir state of MESMER (manual, Sec. 14.2.1).
  - **What it keeps.** It keeps their full Boltzmann weight. Collisions bring molecules into it with the normalized downward probabilities, and activation out of it follows by detailed balance. Partition functions and thermal rate coefficients therefore stay complete.
  - **When it applies.** It forms only where the normalization requires it, and RUN SETTINGS lists it per temperature.
  - **Documentation:** `reports/low_energy_reservoir_state.md` (physical meaning, equations, differences from MESMER's reservoir state, alternatives, validity). How the rule was reached (the temperature step of the former reduction rule, and the population that truncation removes): `reports/low_energy_reduction_temperature_step.md`.
- Olzmann **stepladder** (Olzmann, Gebhardt, Scherzer, Int. J. Chem. Kinet. 23, 825 (1991)).
- Both obey detailed balance exactly; Lennard-Jones collision frequencies (Troe 1977).

**Sources**
- Thermal entrance channels (chemical activation from a bimolecular reactant).
- A given distribution.
- **Consecutive chemical activation:** coupled master equations, the output of one feeding the next (PO14 pp. 236–237).
- **Prepared experiments** (`Preparation` block; [overview](docs/features/prepared_experiments.md)):
  - **Initial population and any number of source channels**, each a normalized distribution over the grains of the wells ([distributions](docs/features/preparation_distributions.md)): thermal at a preparation temperature of its own; Gaussian; single energy; tabulated; the thermal entrance flux; mixtures; shifts.
  - **[Time profiles](docs/features/source_time_profiles.md) of the amplitudes:** impulse, rectangular, Gaussian, constant feed, precursor decay, tabulated, trains.
  - **[Bath history](docs/features/bath_history.md):** a piecewise-constant sequence of T and p, e.g. a shock.
  - **[Transient diagnostics](docs/features/transient_diagnostics.md):** flux coefficients, incubation and vibrational relaxation times.
  - **CSE:** the [source projection](docs/features/cse_source_projection.md), and the [CSE description in time](docs/features/cse_validity_in_time.md) against the master equation.
- **[Fragment wells and lumped reactant states](docs/features/fragment_wells_and_lumped_reactants.md):** a fragment that reacts on is formed with an energy distribution from a kernel (prior, modified prior, two-piece Gaussian); a bimolecular reactant can be one state of the master equation, followed in time.

MarXus has **four solution methods in three families**: [SteadyStateOlzmann](docs/methods/steady_state_olzmann.md) and [SteadyStateAbsorbingBarrier](docs/methods/steady_state_absorbing_barrier.md) (steady state), [CSE](docs/methods/chemically_significant_eigenvalues.md) (eigenvalue) and [TimeIntegration](docs/methods/direct_time_integration.md) (time integration). One is chosen per run, in the `MarXus` block of the deck header or on the command line (see Usage); there is no default.

**Steady-state family, J·N = F (GO10 eqs. 7, 8), two methods (GO10 Sec. 3.2), one per run:**
- [SteadyStateAbsorbingBarrier](docs/methods/steady_state_absorbing_barrier.md): the intermediate steady state, with an absorbing barrier at a user-chosen distance below the threshold (default 10 kT).
- [SteadyStateOlzmann](docs/methods/steady_state_olzmann.md): the final steady state, with physical sinks and no absorbing barrier. Its 𝐉 also gives the thermal rate coefficients (below).
- Solvers: banded Cholesky of the symmetrized operator, or BiCGSTAB.
- Every result is checked for its residual and mass balance.

**Observables**
- Yields of products, stabilization and sink (Olzmann, PCCP 4, 3614 (2002), eq. 10).
- Chemical-activation rate coefficients k^ca (González-García, Olzmann, PCCP 12, 12290 (2010), eq. 9).
- Total loss rate coefficients.
- Bimolecular rate coefficients k(R → X) = k∞Φ_X (PR03 eq. 44).

**Thermal rate coefficients of the final steady state** (part of the final steady state, not a method of its own: the lowest eigenpair of the same J, GO10 eq. 12; its eigenvector is the thermal steady-state population)
- **k_uni**: the average of k(E) over the thermal eigenvector (GO10, after eq. 12). It is reported beside λ₁, the lowest eigenvalue (GO10 eq. 12), with the sum-rule check λ₁ = k_uni. A warning is given above 1.5%, and the double-precision floor is printed.
- Channel rate coefficients, high-pressure limits, λ₂/k_uni (separation of time scales).
- Association rate coefficients by detailed balance.
- Three eigen-solvers:
  - shifted inverse iteration with the banded Cholesky factor (default, most accurate);
  - Householder/QL (EISPACK tred2/tql2, Olzmann's route);
  - LAPACK DSYEVD.
- **Time-dependent populations** N(t) by eigenvalue expansion (PO14 eqs. 3–4), and the validity windows of the steady-state picture (O02).

**Eigenvalue family: [CSE](docs/methods/chemically_significant_eigenvalues.md), phenomenological rate coefficients from the chemically significant eigenvalues** (`Method CSE`, `--method cse`):
- The method is that of Miller, Klippenstein (J. Phys. Chem. A 110, 10528 (2006)) in the formulation of Georgievskii et al. (J. Phys. Chem. A 117, 12146 (2013)).
- It produces species-to-species tables of well ↔ well, well → products, reactant → wells and reactant → products, as MESS does.
- Diagnostics: the eigenvalue separation, the relaxational projections, the loss balance, detailed balance and the precision floor.
- It reproduces a four-well MESS run to a few percent (`validation/ZZAllyl+O2_Gamma_Case2/`).

**Time-integration family: [TimeIntegration](docs/methods/direct_time_integration.md)** (`Method TimeIntegration`, `--method time-integration`):
- The grained populations and the yields of every exit are integrated in time, from a pulse or under continuous formation, with the adaptive, L-stable Rosenbrock methods Ros2–Rodas4 adapted from KPP (`src/numeric/integrators/`).
- At long times the yields equal those of SteadyStateOlzmann exactly ($`k^T 𝐉^{-1} F`$).

**Input**
- Decks in the MESS input format (a subset), so existing decks can be used.
- A MarXus keyword block adds ILT parameters to barriers.
- A `MarXus ... End` block inside a `Well` gives that well its own collision parameters in place of the global model: `Factor[1/cm]`, `Power`, `ReferenceTemperature[K]`, and the pairs `Epsilons[1/cm]`, `Sigmas[angstrom]`, `Masses[amu]` (bath, complex).
- A `MarXus ... End` block in the deck header selects the solution method and its settings (see Usage).
- Exact Eckart tunneling and 1 cm⁻¹ state counting are built from the deck.

**Validation**
- `validation/c2h3_mess_example/` compares H + C₂H₂ ⇌ C₂H₃ at 300–2000 K and 0.1–10 atm with the stored MESS results. It covers all four methods and the three eigen-solvers, in one directory.
- The thermal rate coefficients (k_uni) agree with MESS within −1.4 … +1.8% from 750 to 2000 K. At 300–500 K the offsets of +2.5 … +4.6% are the known exact vs semiclassical Eckart difference.
- At 300–1000 K the associations of all four methods agree (absorbing barrier vs k_uni·K within 0.02%). Above 1000 K the absorbing barrier loses its plateau.
- `validation/acetyl_o2_mesmer/` reproduces the MESMER example AcetylO2 (acetyl + O₂: two wells, ILT association, Eckart tunneling, two sinks). The decks are written from MESMER's XML by `make_deck.py`.
  - MarXus reproduces MESMER's partition functions within 0.001%.
  - The Bartis–Widom rate coefficients agree within −2.1 … +1.4% wherever both codes resolve them.
  - The systematic −1 … −2% is MESMER's 100 cm⁻¹ grains: its canonical rate constants are +1.6 … +2.0% above their exact values, MarXus's +0.01 … +0.24%.
- `validation/ZZAllyl+O2_Gamma_Case2/` reproduces a four-well MESS run: ZZ-allyl + O₂ with two PST channels, Eckart tunneling and an escape sink.
  - The capture rate agrees within 0.5%.
  - Every channel's k∞ agrees within 0.3–0.7% once the tunneling factor is accounted for. MarXus's exact Eckart κ is 2–23% larger than MESS's.
  - The IEPOX + OH share is 1.5–8.0% in MarXus and 1.3–7.2% in MESS.
  - The escape share agrees within 0.93%.
- All validity checks of the master equation (sum rules, detailed balance, limits, solver agreement) are described in `reports/master_equation_validity_checks.md`.

### Dissociative photoionization (`src/photoion/`)
- Breakdown curves of threshold photoelectron photoion coincidence (PEPICO) experiments without collisions, after Sztáray, Bodi, Baer, J. Mass Spectrom. 45, 1233 (2010):
  - the thermal distribution of the neutral (eq. 1), transposed onto the ion at each photon energy (eq. 2) and convolved with a Gaussian energy resolution;
  - fast dissociation (eq. 23), or RRKM rate constants (eq. 6) competing within the flight time (kinetic shift, eqs. 21 and 24);
  - statistical partitioning of the excess energy between the fragment ion, the neutral fragment and their relative translation (1–3 degrees of freedom, eq. 5), passed on to sequential dissociations.
- The deck has MESS-format species blocks (`Species NAME` with RRHO or Atom) and a `Photoionization` block with the channels and their appearance energies.
- Run with `cargo run --release --example photoionization_from_deck -- deck.inp --csv curves.csv`.

### Numerical library (`src/numeric/`)
- Banded Cholesky (factor once, solve many) and LDLᵀ with Bunch–Kaufman pivoting.
- BiCGSTAB and GMRES.
- Symmetric eigensolvers: Householder + implicit QL, Jacobi, inverse iteration with deflation and shift.
- LAPACK interface (DSYEVD, system OpenBLAS); tridiagonal solvers; Lanczos Γ function.

---

## Master-equation solvers: one master equation, four methods, three families

**All methods of MarXus solve the same master equation.** They use the same grains, the same collision operator, the same microcanonical rate coefficients k(E) and the same source of chemically activated adducts. They differ only in the **question** they ask of that equation, that is, which solution of it is computed. Their results are therefore **different quantities, not approximations of one quantity**.

**One method per run.** The two steady-state methods belong to one family, but they are different ways of solving, and each run computes one of them. Each method is described in detail in its own document (follow the links).

| family | method (`Method` / `--method`) | the question it answers |
|---|---|---|
| steady state | [SteadyStateOlzmann](docs/methods/steady_state_olzmann.md) (`steady-state-olzmann`): final steady state, with its thermal eigenpair | Under continuous formation, once the stabilized population has itself reached a steady state (no net stabilization), what are the yields of all products and sinks? And, from the lowest eigenpair of the same 𝐉: how fast does a thermalized adduct react? |
| steady state | [SteadyStateAbsorbingBarrier](docs/methods/steady_state_absorbing_barrier.md) (`steady-state-absorbing-barrier`): intermediate steady state, as in SSUMES | What happens to freshly formed adducts on their first collisional descent? Which fraction redissociates, which decomposes to each product (bimolecular-to-bimolecular), and which is stabilized (bimolecular-to-well)? |
| eigenvalue | [CSE](docs/methods/chemically_significant_eigenvalues.md) (`cse`; Miller, Klippenstein; Georgievskii et al.) | Which species-to-species rate coefficients (reactants, wells, products) define a kinetic model that reproduces the master-equation kinetics after collisional relaxation? |
| time integration | [TimeIntegration](docs/methods/direct_time_integration.md) (`time-integration`; Rosenbrock, adapted from KPP) | How do the populations of all grains and the yields of all exits evolve in time, from a pulse of activated adducts or under continuous formation, through relaxation and chemistry? |

**What each method gives:**
- The two steady-state methods give flux coefficients and yields.
- CSE gives phenomenological rate coefficients.
- TimeIntegration gives the time evolution itself.

**Every method reports the reactant explicitly:**
- **bimolecular-to-bimolecular** rate coefficients and yields, e.g. R → IEPOX + OH (chemical activation);
- where defined, **bimolecular-to-well** (stabilization) rate coefficients and yields, e.g. R → G4. Miller and Klippenstein (MK06, pp. 10529, 10531): "application of the steady-state approximation … is virtually always an attempt to equate a phenomenological rate coefficient to a flux coefficient. Sometimes this is a valid approach, and sometimes it is not."

**References used in this section**

| code | reference |
|---|---|
| O91 | Olzmann, Gebhardt, Scherzer, Int. J. Chem. Kinet. 23, 825 (1991) |
| O02 | Olzmann, Phys. Chem. Chem. Phys. 4, 3614 (2002) |
| GO10 | González-García, Olzmann, Phys. Chem. Chem. Phys. 12, 12290 (2010) |
| PO14 | Pfeifle, Olzmann, Int. J. Chem. Kinet. 46, 231 (2014) |
| SN84 | Schranz, Nordholm, Chem. Phys. 85, 163 (1984) |
| PR03 | Pilling, Robertson, Annu. Rev. Phys. Chem. 54, 245 (2003) |
| MK06 | Miller, Klippenstein, J. Phys. Chem. A 110, 10528 (2006) |
| G13 | Georgievskii, Miller, Burke, Klippenstein, J. Phys. Chem. A 117, 12146 (2013) |
| M79 | Miller, J. Am. Chem. Soc. 101, 6810 (1979) |

### The master equation

The populations $`N_i`$ of the energy grains $`i`$ of all wells obey (PO14 eq. 2; O02 eq. 6; GO10 eq. 7), written per grain:

```math
\frac{dN_i}{dt} = R\,F_i - \sum_j J_{ij}\,N_j .
```

The same equation in matrix notation:

```math
\frac{dN}{dt} = R\,F - 𝐉\,N,
\qquad
𝐉 = \omega\,(𝐈 - 𝐏) + 𝐊 + k_c[\mathrm D]\,𝐈 .
```

**Notation.** **Bold** upper-case letters are matrices; vectors and scalars are in normal type. A product is written by juxtaposition: $`𝐉 N`$ is the matrix $`𝐉`$ times the vector $`N`$, and $`R\,F`$ is the scalar $`R`$ times the vector $`F`$.

| symbol | kind | meaning |
|---|---|---|
| $`N = (N_1, N_2, \dots)`$ | vector | population of every grain, all wells stacked |
| $`R`$ | scalar | total formation rate of the chemically activated adducts, e.g. $`R = k_\infty [\mathrm A][\mathrm B]`$ |
| $`F = (F_1, F_2, \dots)`$ | vector | normalized nascent distribution, $`\sum_i F_i = 1`$ |
| $`R\,F`$ | vector | scalar $`R`$ times vector $`F`$: the formation rate into grain $`i`$ is $`R\,F_i`$ |
| $`𝐉`$ | matrix (s⁻¹) | collisions, reactions and sinks; it has the elements below |
| $`𝐉\,N`$ | vector | matrix–vector product, $`(𝐉\,N)_i = \sum_j J_{ij} N_j`$: the net rate at which population leaves grain $`i`$ |
| $`\omega`$ | scalar (s⁻¹) | collision frequency (Lennard-Jones) |
| $`𝐏`$ | matrix | collisional transition probabilities, $`P_{ij} = P(E_i \leftarrow E_j)`$ (exponential down by default) |
| $`𝐊`$ | matrix (s⁻¹) | microcanonical reactions: products and isomerization |
| $`k_c[\mathrm D]`$ | scalar (s⁻¹) | pseudo-first-order bimolecular sink of a well, e.g. an escape channel |
| $`𝐈`$ | matrix | identity |

**Elements of $`𝐉`$.** Grain $`j`$ belongs to well $`w`$, and grain $`i'`$ is the grain of well $`w'`$ at the same absolute energy:

```math
J_{ij} = \omega\,\big(\delta_{ij} - P_{ij}\big) + \delta_{ij}\,\Big(\sum_r k_r(E_j) + k_c[\mathrm D]\Big) - \delta_{i i'}\,k_{w \to w'}(E_j) .
```

The sum over $`r`$ runs over all channels of well $`w`$: products and isomerizations. The last term puts the isomerization flux into the other well.

**Detailed balance of the collisions.** The probabilities obey $`P(E' \leftarrow E)\,f^0(E) = P(E \leftarrow E')\,f^0(E')`$, with $`f^0(E) = \rho(E)\,e^{-E/k_BT}`$.

**Properties.** $`𝐉`$ has positive eigenvalues when population can leave the network (GO10, text before eq. 12). With $`𝐃 = \mathrm{diag}\big(\sqrt{f^0}\big)`$, the matrix $`𝐒 = 𝐃^{-1}\,𝐉\,𝐃`$ is symmetric, and all solvers work with $`𝐒`$. The CSE literature writes the same equation as $`d|f\rangle/dt = -\hat{𝐆}\,|f\rangle + \sum_\nu s_\nu\,|p^{(\nu)}\rangle`$ (G13 eq. 1), with $`\hat{𝐆} = 𝐉`$.

### Chemical activation

**What it is.** Reactants A + B associate (or react) through an entrance channel to an adduct AB\*. The adduct is born with the energy released by the association plus the thermal energy of the reactants. It therefore lies above the threshold $`E_0`$ of its own redissociation: it is *chemically activated*. Its fates compete:

- redissociation back to A + B;
- decomposition to products, directly or after isomerization to other wells, often before any collision ("well-skipping");
- collisional stabilization into the thermal distribution of a well;
- reaction with a bimolecular partner (the sink $`k_c[\mathrm D]`$).

Pressure (through $`\omega`$) sets the balance between reaction and stabilization. The stabilized adducts react later, thermally and more slowly. **How a solver treats this second stage is the main difference between the three.**

**The initial flux: the source from thermal reactants.** The nascent distribution follows from the reverse (dissociation) channel by detailed balance (PO14 eq. 7; O02 eq. 11):

```math
F(E) = \frac{W^\ddagger(E-E_0)\;e^{-(E-E_0)/k_BT}}{\displaystyle\int_0^\infty W^\ddagger(\varepsilon)\;e^{-\varepsilon/k_BT}\,d\varepsilon},
\qquad E \ge E_0 .
```

Here $`W^\ddagger`$ is the sum of states of the entrance transition state (rigid, phase space theory or inverse Laplace transform). Since $`k(E) = W^\ddagger(E-E_0)/[h\,\rho(E)]`$ (PO14 eq. 9), the source needs only the microcanonical rate coefficient of the reverse reaction:

```math
F(E) \;\propto\; \rho(E)\;k_{\to\mathrm{A+B}}(E)\;e^{-E/k_BT}.
```

**Several entrance channels.** With several entrance channels, into one well or several, the weights are taken on the common absolute energy scale and normalized over all wells (`chemical_activation_sources.rs`). Each channel then contributes in proportion to its thermal flux:

```math
F_w(E) \;\propto \sum_{c\,\in\,\mathrm{entrances\ of\ }w} \rho_w(E)\;k_c(E)\;e^{-E/k_BT},
\qquad \sum_w \sum_E F_w(E) = 1 .
```

**Formation rate.** The formation rate is $`R = k_\infty\,[\mathrm A]\,[\mathrm B]`$. The high-pressure (capture) rate coefficient $`k_\infty`$ comes from the same $`W^\ddagger`$ (`EntranceHighPressureRate`):

```math
k_\infty(T) = \frac{\sum_E W^\ddagger(E)\;e^{-(E-E_{\mathrm{AB}})/k_BT}\,\Delta E}
{h\,\left(2\pi\mu k_BT/h^2\right)^{3/2}\,Q_{\mathrm A}(T)\,Q_{\mathrm B}(T)} ,
```

with $`E_{\mathrm{AB}}`$ the asymptote of A + B, $`\mu`$ the reduced mass, and $`Q_{\mathrm A}`$, $`Q_{\mathrm B}`$ the internal partition functions.

**The source in the CSE method.** The CSE method uses the same source shape, $`s_R\,p^{(R)}(E)`$ with $`p^{(R)}(E) = k_{\to R}(E)\,f^0(E)`$ (G13 eqs. 1, 9). It is normalized by the capture rate coefficient, $`1/Q_R = k_\infty / \sum_E k_{\to R}(E)\,f^0(E)`$ (G13 eq. 23).

**Non-thermal reactants** (consecutive chemical activation; library, steady-state solvers; PO14 eqs. 8, 10–13). With normalized reactant distributions $`n_{\mathrm A}`$, $`n_{\mathrm B}`$:

```math
F(E) = \int_0^{E-E_0} n_{\mathrm A}(\varepsilon)\;n_{\mathrm B}(E-E_0-\varepsilon)\,d\varepsilon .
```

In the shift approximation this becomes $`F(E) = n_{\mathrm A}\big(E + E_R - \langle E_{\mathrm B}\rangle\big)`$, with $`E_R`$ the 0 K reaction energy.

### Side by side

| | [SteadyStateAbsorbingBarrier](docs/methods/steady_state_absorbing_barrier.md) | [SteadyStateOlzmann](docs/methods/steady_state_olzmann.md) | [CSE](docs/methods/chemically_significant_eigenvalues.md) | [TimeIntegration](docs/methods/direct_time_integration.md) |
|---|---|---|---|---|
| master equation | same J, same F | same J, same F | same J, same F | same J, same F |
| mathematical problem | linear system on the grains above the absorbing barriers | linear system on all grains; lowest eigenpair of the same J | all eigenpairs of J | initial-value problem, adaptive Rosenbrock steps |
| stabilized adducts | removed (counted as stabilized) | stay in the well and react thermally | a chemical eigenmode per well | stay in the well; followed in time |
| time scale | $`(0.1\lambda_F)^{-1} < t < (10k_{\mathrm{uni}})^{-1}`$ | $`t \gg 1/k_{\mathrm{uni}}`$ | $`t \gg 1/\Lambda_{n_w+1}`$ | every t |
| output | yields Φ, $`k_\infty\Phi_X`$ | $`k^{ca}`$, yields Φ, sink yields; $`k_{\mathrm{uni}}`$, $`\lambda_1`$, $`k^{th}_r`$ | species-to-species rate coefficients | N(t) of every well, Y(t) of every exit |
| kind of quantity | flux coefficient | flux coefficient; thermal rate coefficient | phenomenological rate coefficient | populations and yields vs time |
| needs | barrier distance (choice) | a sink or exit for a well-conditioned J·N = F | eigenvalue separation | output time range |
| typical use | association falloff, prompt product branching | yields with physical sinks (continuous formation), consecutive activation, thermal $`k_{\mathrm{uni}}`$ | rate coefficients for a kinetic mechanism | time-resolved experiments, checks of the other solvers |

### Where the solvers must agree: identities and measured agreement

**Exact identity 1: pulse and continuous formation have the same long-time yields.** For the linear master equation with a normalized initial distribution $`F`$ (a pulse) and no further source, the population decays as $`N(t) = e^{-𝐉 t} F`$, and the yield of channel $`r`$ accumulated up to $`t \to \infty`$ is

```math
Y_r(\infty) = \int_0^\infty k_r^T\, e^{-𝐉 t} F \, dt = k_r^T\, 𝐉^{-1} F ,
```

since all eigenvalues of $`𝐉`$ are positive. This is exactly the yield of the final steady state per formed adduct with the same $`F`$ (GO10 eqs. 8, 9). A pulsed experiment and a continuously fed experiment therefore share their integrated yields, while their time traces differ. Vereecken et al. found the two "identical for all practical purposes" (J. Chem. Phys. 106, 6564 (1997), Table I).

**Measured for identity 1.** The direct time integration of a pulse (TimeIntegration) on the four-well ZZ-allyl + O₂ network reaches, at t = 100 s, the final-steady-state yields of every exit at all 21 conditions, to the 7 printed digits (`validation/ZZAllyl+O2_Gamma_Case2/time_integration_vs_final_steady_state.csv`, `plots/time_evolution_300K_760torr.png`).

**Exact identity 2: the CSE rate coefficients reproduce the final steady state.**
- **The reconstruction.** Take the long-time yields from the CSE rate coefficients: R forms the wells and the direct products, and each well then ends in a product, the escape or back in R (the absorbing chain of the CSE well rate coefficients, `chemically_significant_eigenvalues::reactant_yields`). These yields equal those of the final steady state, as fractions of the net reaction.
- **Why.** With G13 eqs. 21 and 25–30, both are $`\sum_\lambda p^{(x)}_\lambda\, p^{(R)}_\lambda / (\Lambda_\lambda Q_R)`$ over all eigenpairs, the spectral form of $`k_x^T 𝐉^{-1} F`$. This holds independently of the eigenvalue separation; the individual CSE coefficients lose their meaning without separation, but this sum does not.
- **What it checks.** The two methods reach the same observable through independent code: a banded Cholesky solve on one side; the full eigendecomposition, $`𝐌^{-1}`$, the rate assembly and the absorbing chain on the other.
- **Measured.**
  - The test `cse_long_time_yields_equal_the_final_steady_state_yields` holds to 10⁻⁸.
  - ZZ-allyl + O₂ (`validation/ZZAllyl+O2_Gamma_Case2/cse_vs_final_steady_state.csv` and `plots/cse_vs_final_steady_state.png`): IEPOX + OH at 300 K and 760 Torr is 2.45778% from both the final steady state and the CSE kinetics.
  - Over all 21 conditions the IEPOX + OH share agrees to 6.4·10⁻⁷ relative, the escape to 2.3·10⁻⁸ and P1 to 5.5·10⁻⁷: the precision of the printed digits. P7, at most 0.008% of the reaction and formed through G6, agrees to 1.4·10⁻⁴; the CSE entries on its path (R → G6 ≈ 10⁻²² cm³/s) are at the rounding level.
- **All identities, measured on both systems:** `reports/method_comparison.md` (script `validation/method_comparison.py`; per system `method_comparison.csv` and `plots/method_*.png`).
- **The four methods in figures, against MESS and against each other** (rates, fall-off, yields, chemical activation, thermal activation, total formation yield, time evolution): `reports/four_methods_figures.md` (script `validation/four_methods_figures.py`; per system `four_methods_figures.csv`, `plots/mess_four_methods_*.png` and `plots/internal_four_methods_*.png`).

**Exact identity 3 (one well): the slowest mode of a pulse gives the CSE association.**
- **The statement.** At late times a pulse decays as $`N(t) = A\,e^{-k_{\mathrm{uni}} t}`$. The amplitude extrapolated to t = 0 is $`A = (\sum_E f^{(1)})(\sum_E f^{(1)} k_R)/\sum_E k_R f^0`$, which is $`k(R \to W)/k_\infty`$ of G13 eq. 28.
- **Measured.** TimeIntegration $`k_\infty A`$ equals the CSE association within 1.9·10⁻⁴ at all 40 C₂H₃ conditions (`reports/four_methods_figures.md`, Section 2).

- **One well.** The CSE well → product rate coefficient equals the eigenvector-average $`k_{\mathrm{uni}}`$ of SteadyStateOlzmann (test `a_single_well_gives_the_eigenvector_average_as_its_rate_coefficient`, to 10⁻⁸).
- **Well separation and continuous formation.** Steady-state yields equal time-integrated single-injection yields (Vereecken et al., J. Chem. Phys. 106, 6564 (1997), Table I). The CSE R → product and R → well rate coefficients equal the fractions accumulated on the relaxation time scale (G13 p. 12153).
- **H + C₂H₂ ⇌ C₂H₃, one well** (`validation/c2h3_mess_example/`):
  - **Identities:**
    - CSE's k(W1 → P1) equals SteadyStateOlzmann's k_uni at every condition;
    - the late-time decay of the time-integrated pulse equals k_uni within 4·10⁻⁶;
    - the pulse ends in the final steady state in all printed digits.

    This holds at the 30 conditions where λ₁ is above the double-precision floor.
  - **Up to 1000 K the associations of all methods agree:** CSE (G13 eq. 28), SteadyStateOlzmann's k_uni·K, and SteadyStateAbsorbingBarrier (10 kT, within 0.02%).
  - **Above 1000 K:**
    - the absorbing barrier loses its plateau (0.1 atm: −8% at 1500 K, −29% at 1750 K, no result at 2000 K);
    - CSE's association departs from k_uni·K exactly by its own departure from detailed balance, −14% at 2000 K. MESS's own pair shows the same departure, −13.9%.
  - **Against MESS:**
    - dissociation (CSE = k_uni): +4.4 … +4.6% at 300 K (tunneling model), +0.1 … +1.1% at 1000 K, +0.3 … +0.4% at 2000 K;
    - CSE association: −0.9 … +0.5% at 1250–2000 K.
- **ZZ-allyl + O₂, four wells** (`validation/ZZAllyl+O2_Gamma_Case2/`):
  - CSE reproduces the MESS species tables within 1.3% (entries above the rounding level).
  - The long-time IEPOX + OH share of SteadyStateOlzmann differs from MESS's long-time fate by +12.0 … +17.8% with the exact Eckart tunneling, and by −0.34 … −0.55% with the MESS tunneling model. The escape share agrees within 0.93%.
  - CSE against MESS (MESS tunneling model): R → G4 +0.43 … +0.60%, R → G3 +0.39 … +0.57%, R → G2 +0.90 … +1.13%, R → P5 +0.28 … +0.50%.
  - The prompt $`k(\mathrm R\to \mathrm{IEPOX+OH})`$ of SteadyStateAbsorbingBarrier depends on the tunneling model:
    - exact Eckart tunneling: +19.0 … +20.3% above MESS's R → P5;
    - MESS's Eckart model: +0.44 … +0.74% from MESS, and within 0.44% of CSE's G13 eq. 21.

    The difference from MESS is therefore the tunneling model, not the definition of the quantity (`reports/method_comparison.md`).

---

## Usage example: master equation from a deck

```
cargo run --release --example chemical_activation_from_deck -- deck.inp REACTANT \
    --method steady-state-olzmann --ncore 4 --csv out.csv > out.report
```

The method and its settings can be given in the deck header, in a `MarXus ... End` block (a MarXus extension of the MESS format). `Method` is required (in the deck or on the command line); the other keywords are optional. The command-line option in the second column overrides the deck:

```
MarXus
  Method                              SteadyStateOlzmann ! SteadyStateOlzmann | SteadyStateAbsorbingBarrier | CSE | TimeIntegration
  AbsorbingBarrierBelowThreshold[kT]  10                 ! SteadyStateAbsorbingBarrier
  EigenSolver                         InverseIteration   ! InverseIteration | FullDecomposition | Lapack
  SumRuleTolerance                    1.5e-2             ! thermal eigenpair of the final steady state
  Integrator                          Rodas4             ! time integration: Rodas4 | Rodas3 | Ros4 | Ros3 | Ros2
  InitialState                        Pulse              ! time integration: Pulse | Continuous
  TimeRange[s]                        1e-12  1e2         ! time integration: first and last output time
  TimesPerDecade                      4                  ! time integration
  IntegrationTolerance                1e-6               ! time integration: relative tolerance
  CompareWithCse                      1e-2               ! time integration of a Preparation in a constant bath: CSE in time
  NCores                              8                  ! cores of the run; the (T, p) conditions run in batches of up to 8
  CollisionIntegral                   Neufeld            ! Lennard-Jones Omega(2,2)*: Neufeld (default) | Troe
End
```

| deck keyword | option | meaning (default) |
|---|---|---|
| `Method` | `--method steady-state-olzmann\|steady-state-absorbing-barrier\|cse\|time-integration` | the method of the run (required, no default; see the method documents) |
| `AbsorbingBarrierBelowThreshold[kT]` | `--barrier-kt X` | absorbing barrier of the intermediate steady state, in k_BT below the lowest threshold (10) |
| `EigenSolver` | `--eigen-solver inverse\|full\|lapack` | thermal eigenpair of the final steady state (inverse iteration); CSE needs all eigenpairs (LAPACK; `full` also possible, inverse iteration refused) |
| `SumRuleTolerance` | `--sum-rule-tolerance X` | warning threshold of \|λ₁ − k_uni\|/k_uni (1.5e-2) |
| – | `--csv FILE` | also write the machine-readable tables (CSV with titled blocks) to FILE, and every table of the report to FILE_tables.csv |
| `Integrator`, `InitialState`, `TimeRange[s]`, `TimesPerDecade`, `IntegrationTolerance` | `--integrator`, `--initial`, `--time-range T1 T2`, `--times-per-decade`, `--integration-tolerance` | time integration: Rosenbrock method (Rodas4), pulse or continuous formation (pulse), output times (1e-12 to 1e2 s, 4 per decade), relative tolerance (1e-6) |
| `CompareWithCse` | `--compare-with-cse X` | time integration of a `Preparation` in a constant bath: also propagate the CSE description of the same preparation and report where it agrees with the master equation; X is the tolerance of the relative deviation (no comparison without it) |
| `NCores` | `--ncore N` | number of cores of the run, for every method; `--ncore` overrides `NCores` of the deck, and either may be given alone. The conditions (T, p) are independent and are computed in batches of up to N at a time; LAPACK calls get the cores left over (N / conditions at a time). Default: RAYON_NUM_THREADS, otherwise all logical cores. RUN SETTINGS shows the number and where it came from |
| `ChemicalEigenvalueMax`, `WellProjectionThreshold` (global section of the deck, MESS's keywords) | – | CSE species merging: the chemical eigenvalues are those ≤ ChemicalEigenvalueMax × the lowest relaxation eigenvalue (0 < value < 1; 0.2, where MESS has no default); wells with a projection on the chemical subspace ≥ WellProjectionThreshold are primary wells of the partition (0.2, as MESS) |
| `CollisionIntegral` | `--collision-integral neufeld\|troe` | reduced collision integral Ω(2,2)\* of the Lennard-Jones collision frequency: Neufeld, Janzen, Aziz, J. Chem. Phys. 57, 1100 (1972) (default; also used by MESS and MESMER), or Troe, J. Chem. Phys. 66, 4758 (1977), eq. 3.3 (±7%). With Neufeld the remaining few-percent differences from MESS of the validations disappear (`reports/collision_integral_neufeld.md`) |
| – | `--tunneling exact-eckart\|mess-eckart` | Eckart transmission model: exact Eckart (default); `mess-eckart` only for comparison with MESS (see Tunneling model below) |

A setting that the chosen method does not use is reported as a note in the output, not refused. `SteadyState` / `--steady-state` and `both` no longer exist: the two steady-state methods are chosen by `Method`, one per run. `eigenvalue` is not a method: the thermal eigenpair belongs to SteadyStateOlzmann.

The output tables begin with comment lines that explain the quantities, with their references. The source file `examples/chemical_activation_from_deck.rs` documents all columns.


### Prepared experiments: the `Preparation` block

A top-level `Preparation ... End` block of the deck (a MarXus extension of the MESS format) describes an experiment by three parts: an initial population, source channels and a bath history:

```math
\frac{d n}{dt} = -𝐉[T(t), p(t)]\, n + \sum_a R_a(t)\, F_a, \qquad n(0) = N_0 F_0, \qquad Y_x(t) = \int_0^t k_x^T n\, dt'
```

Every F is normalized over the grains of all wells; the rates R_a(t) carry the amounts.

```
Preparation
  InitialPopulation
    Amount                      1.0
    Distribution Thermal
      Well                      W1
      PreparationTemperature[K] 300
    End
  End
  Source laser
    Profile Impulse
      Time[s]                   1e-9
      Amount                    0.5
    End
    Distribution Gaussian
      Well                      W1
      Centre[kcal/mol]          45.0
      Width[1/cm]               300       ! standard deviation
    End
  End
  Bath
    Segment Start[s] 0      Temperature[K] 300    Pressure[atm] 1
    Segment Start[s] 1e-6   Temperature[K] 1000   Pressure[atm] 1
  End
End
```

| feature | purpose | main controls |
|---|---|---|
| [Prepared experiments](docs/features/prepared_experiments.md) | the block as a whole; what each method does with it | `InitialPopulation`, `Source <name>`, `Bath` |
| [Distributions](docs/features/preparation_distributions.md) | the energy distribution of the prepared molecules | `Thermal`, `Gaussian`, `SingleEnergy`, `Tabulated`, `ThermalEntrance`, `Mixture`; `EnergyReference`, `Shift[unit]` |
| [Time profiles](docs/features/source_time_profiles.md) | when the source molecules are added | `Impulse`, `Rectangular`, `Gaussian`, `Feed`, `PrecursorDecay`, `Tabulated`, `Train` |
| [Bath history](docs/features/bath_history.md) | shocks and other changes of T and p | `Segment Start[s] Temperature[K] Pressure[unit]` |
| [Transient diagnostics](docs/features/transient_diagnostics.md) | what an experiment observes: flux coefficients r, k^e (Barker, Frenklach, Golden 2015), E_f, τ_vib, τ_inc (Barker, King 1995), collision numbers | `TimeRange[s]`, `TimesPerDecade`, `IntegrationTolerance` |
| [CSE source projection](docs/features/cse_source_projection.md) | the species populations after relaxation and the prompt yields of every source (Georgievskii et al. 2013) | `Method CSE` with a `Preparation` |
| [CSE validity in time](docs/features/cse_validity_in_time.md) | when the rate-coefficient description agrees with the master equation for this preparation | `CompareWithCse X`, `--compare-with-cse X` |
| [Fragment wells and lumped reactant states](docs/features/fragment_wells_and_lumped_reactants.md) | a dissociation into a fragment that is a well of the network, with its energy partitioning (Green, Robertson 2014; Shannon et al. 2024); the bimolecular reactant as one state followed in time (Miller et al. 2016, SI-VI) | `FragmentWell`, `FragmentEnergy Prior\|ModifiedPrior\|TwoPieceGaussian`, `LumpedState`, `ExcessFragment`, `PartnerConcentration[molecule/cm^3]` |

Examples: `examples/c2h3_prepared_experiment.inp` (laser pulse, entrance feed and shock in one run); `validation/c2h3_mess_example/input/c2h3_shock_incubation.inp` (shock heating with the CSE comparison).

### Tunneling model (`--tunneling`)

**Default: the exact Eckart transmission probability**, used by MarXus for all `Tunneling Eckart` blocks of a deck (M79 eq. 8; Johnston, Heicklen, J. Phys. Chem. 66, 532 (1962)):

```math
P(E_1) = \frac{\sinh a\,\sinh b}{\sinh^2\big(\tfrac{a+b}{2}\big) + \cosh^2 c},
\qquad
a = \frac{4\pi}{\hbar\omega}\,\frac{\sqrt{E_1 + V_0}}{V_0^{-1/2} + V_1^{-1/2}},
\quad
b = \frac{4\pi}{\hbar\omega}\,\frac{\sqrt{E_1 + V_1}}{V_0^{-1/2} + V_1^{-1/2}},
\quad
c = 2\pi\sqrt{\frac{V_0 V_1}{(\hbar\omega)^2} - \frac{1}{16}} .
```

Here $`E_1`$ is the energy in the reaction coordinate relative to the barrier top, $`V_0`$ and $`V_1`$ are the barrier heights relative to the two sides, and $`\hbar\omega`$ is the magnitude of the imaginary frequency. When $`V_0V_1/(\hbar\omega)^2 < 1/16`$, $`\cosh c`$ becomes $`\cos|c|`$.

**How it enters the rate coefficients.**
- **Microcanonically:** the transition-state sum of states is convolved with $`dP/dE`$ (M79 eq. 9).
- **Canonically:** $`\kappa(T) = \beta\,e^{\beta V_0}\int P(E-V_0)\,e^{-\beta E}\,dE`$, with $`E`$ measured from the forward asymptote and $`V_0`$ the forward barrier height (`examples/eckart_kappa_from_deck.rs`).

**Optional: `--tunneling mess-eckart`**, only for users who want a closer comparison with MESS. It is MESS's semiclassical model (`tunneling/mess_eckart_tunneling.rs`, mirroring the MESS source):

```math
P(E) = \frac{1}{1 + e^{-S(E)}},
\qquad
S(E) = \frac{4\pi}{d_0^{-1/2} + d_1^{-1/2}} \sum_{w=0,1} \Big[\sqrt{\max(E/\omega + d_w,\,0)} - \sqrt{d_w}\Big],
\qquad d_w = V_w/\omega .
```

It also uses MESS's clamps ($`|S| > 100`$), its cutoff at the smaller well depth, and its canonical weight on 0.01 kT steps up to 10 kT.

**What changes with this option.** For deep tunneling it gives smaller factors than the exact Eckart. For the barriers of `validation/ZZAllyl+O2_Gamma_Case2/` at 270–330 K, the exact Eckart κ is 2–23% larger than the MESS model's, most for the H-transfer barriers.
- **Agreement with MESS.** With the option, every high-pressure rate coefficient of that network agrees with MESS within 1.7%.
- **Remaining difference.** The IEPOX + OH share is then 0.34–0.55% below MESS (Neufeld collision integral, the default). With the exact Eckart it is 12–18% above MESS.

**Use the default (exact Eckart) for results.** Use `mess-eckart` only to separate the tunneling difference from other differences when comparing with a MESS run.

---

## Installation and dependencies

MarXus needs **one system library (OpenBLAS)** besides the Rust toolchain. All other dependencies are Rust crates that `cargo` downloads and compiles by itself; nothing else has to be installed or maintained by hand.

### 1. Rust toolchain (required)

Install Rust with rustup (https://rustup.rs):

```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Tested with rustc/cargo 1.98.

### 2. OpenBLAS: BLAS and LAPACK (system library, needed by the default build)

**What uses it.** The default cargo feature `openblas` links the system OpenBLAS. MarXus uses it for the LAPACK symmetric eigensolver DSYEVD, `--eigen-solver lapack` in the master-equation eigenvalue analysis.

**Install it:**

| system | command | status |
|---|---|---|
| Ubuntu / Debian | `sudo apt install libopenblas-dev` | tested (Ubuntu 24.04, libopenblas-dev 0.3.26) |
| Fedora / RHEL | `sudo dnf install openblas-devel` | not tested |
| Arch Linux | `sudo pacman -S openblas` | not tested |
| macOS (Homebrew) | `brew install openblas`, then before building `export RUSTFLAGS="-L $(brew --prefix openblas)/lib"` (Homebrew does not put OpenBLAS on the default library path) | not tested |

**Check that it is found:**

```
ldconfig -p | grep libopenblas        # Linux: should list libopenblas.so
```

**Without OpenBLAS.** If the library is not available, build without it:

```
cargo build --release --no-default-features
```

Everything works except `--eigen-solver lapack`, which then stops with an error naming the missing feature. The default solver (inverse iteration) and the in-house full decomposition (`--eigen-solver full`) need no external library.

**Threads.** OpenBLAS uses all cores by default. To limit it, e.g. to 4 cores:

```
export OPENBLAS_NUM_THREADS=4
```

### 3. Rust crates (automatic)

- `rayon` (pure Rust): the conditions (T, p) of a master-equation run are computed in parallel (`NCores` in the deck, `--ncore N`).

`cargo build` downloads and compiles it by itself; it needs no system library.

The planned higher-precision eigenvalue analysis (`reports/higher_precision_decision.md`) will use pure-Rust crates only:
- `qd`: double-double;
- `faer`: linear algebra;
- `dashu-float`: arbitrary precision.

`cargo build` will download and compile them by itself. They need no system library.

### 4. Python (optional, only for the validation plots)

The plotting scripts in `validation/*/plot_comparison.py` need Python 3 with `numpy` and `matplotlib`:

```
python3 -m venv ~/.venvs/science && source ~/.venvs/science/bin/activate && pip install numpy matplotlib
```

### Build and test

```
cargo build --release -j 4
cargo test -j 4 -- --test-threads=4
```

`-j 4` and `--test-threads=4` limit compilation and tests to 4 cores.

---

---

## To Do (not implemented yet)
- Higher precision for the master equation: double-double assembly and solvers with an arbitrary-precision reference path (planned, `reports/higher_precision_decision.md`).
- Treatment of the stepladder model in the eigenvalue analysis when its step spans several grains (independent sub-equations).
- Photoionization: statistical channels of fragment ions (formation time within the flight time), isomerization (Sztáray et al. 2010, eqs. 12–22), time-of-flight peak shapes, PST and simplified SACM rate models, fitting of appearance energies.
- Microcanonical Variational TST (μVTST).
- Canonical Variational TST (CVTST).
- Submerged barrier with a pre-reaction vdW complex: **μ-canonical, J-resolved 2-TST treatment**.

---

## License
This project is licensed under the **GNU General Public License v3.0 (GPL-3.0)**.
