# MarXus

**Molecular Statistical Physics for Kinetics and Thermochemistry**

**Author:** Peter Szabo  
**Email:** peter88szabo@gmail.com  

MarXus is a **microcanonical rate code** and **master-equation solver** for gas-phase reaction kinetics, written in Rust. From the molecular data of the reactants, wells and transition states, it computes energy-resolved rate coefficients and solves the energy-grained master equation for pressure- and temperature-dependent kinetics, including chemical activation. The same molecular data also give the thermochemistry.

- **Microcanonical rate coefficients k(E).**
  - Sums and densities of states are counted directly on 1 cm⁻¹ cells (Beyer–Swinehart).
  - Tight transition states are treated by RRKM theory, with exact Eckart tunneling.
  - Barrierless channels are treated by phase space theory (TST levels T, E, EJ) or by the inverse Laplace transform of k∞(T).
  - SACM is in progress.
- **Master equation for multiwell networks.** It includes collisional energy transfer (exponential down or stepladder), isomerization, product channels, bimolecular sinks, and chemically activated formation from bimolecular reactants. It is solved by **two methods that answer different questions of the same equation**:
  - the **steady state**, either intermediate (absorbing barrier) or final (Olzmann). The final one includes the thermal rate coefficient from the lowest eigenpair of the same matrix.
  - the **chemically significant eigenvalues**, which give phenomenological rate coefficients for kinetic models (Miller, Klippenstein; Georgievskii et al.).
- **Input** in the MESS deck format, so existing decks can be used, with MarXus extension blocks for its own settings.
- **Thermochemistry:** partition functions and U, H, F, G, S, Cv, Cp in the RRHO and quasi-RRHO (Grimme) approximations.

The master-equation, tunneling, ILT and numerical code cites the source of each equation (paper and equation number) in its comments. Design notes, derivations and validation results are in `reports/`; validations against reference calculations are in `validation/`.

---

## Current state (2026-10-05)

| area | status |
|---|---|
| Thermochemistry (RRHO, Grimme qRRHO) | implemented |
| Sum and density of states, RRKM k(E), canonical TST | implemented, tested |
| Tunneling (exact Eckart, microcanonical and canonical) | implemented, tested |
| Inverse Laplace transform (ILT) for barrierless channels | implemented, tested |
| Phase space theory (PST) for −C_n/Rⁿ potentials, levels T, E, EJ | implemented, tested, used in the master equation (validated, ZZ-allyl + O₂ Case 2) |
| PST with arbitrary 1D potential, SACM | in progress |
| Multiwell chemical-activation master equation (steady states) | implemented, tested, validated (C₂H₃; four-well ZZ-allyl + O₂ Case 2) |
| Thermal rate coefficients of the final steady state (lowest eigenpair of J: k_uni, λ₁, λ₂; N(t)) | implemented, tested, validated (C₂H₃, 300–2000 K) |
| CSE method (phenomenological rate coefficients, Miller–Klippenstein / Georgievskii et al. 2013) | implemented, tested, validated against a four-well MESS run |
| Higher precision (double-double, arbitrary-precision reference) | planned (`reports/higher_precision_decision.md`) |

The source contains 176 library unit tests (`cargo test`).

---

## Features

### Thermochemistry
- Molecules built from vibrational frequencies (with scaling), rotational constants, or Cartesian geometries; the moments of inertia are computed from the geometry with isotopic atomic masses (most abundant isotopes, AME2020: Wang et al., Chin. Phys. C 45, 030003 (2021)).
- Partition functions and thermodynamic functions **U, H, F, G, S, Cv, Cp** in the rigid-rotor–harmonic-oscillator approximation.
- **Quasi-RRHO entropy** with free-rotor interpolation for low-frequency modes (Grimme, Chem. Eur. J. 18, 9955 (2012)).
- Equilibrium constants: used inside the master equation (k∞,assoc/k∞,diss for detailed balance); a general thermochemistry routine is not available yet (see To Do).

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

The energy-grained master equation is dN/dt = R·F − J·N with J = ω(I − P) + K + k_c[D]·I (Pfeifle, Olzmann, Int. J. Chem. Kinet. 46, 231 (2014), eq. 2); its steady state is J·N = R·F. See "Master-equation solvers" below for the questions each solver answers.

**Network**
- Any number of wells on a common absolute energy grid, connected by isomerization (exact detailed balance, k = W‡/(hρ) in both directions).
- Product channels, and pseudo-first-order bimolecular sinks k_c[D].

**Collisions**
- **Exponential down** (default), exactly normalized (Robertson (ed.), Comprehensive Chemical Kinetics 43 (2019), eq. 4.16).
- Olzmann **stepladder** (Olzmann, Gebhardt, Scherzer, Int. J. Chem. Kinet. 23, 825 (1991)).
- Both obey detailed balance exactly; Lennard-Jones collision frequencies (Troe 1977).

**Sources**
- Thermal entrance channels (chemical activation from a bimolecular reactant).
- A given distribution.
- **Consecutive chemical activation:** coupled master equations, the output of one feeding the next (PO14 pp. 236–237).

MarXus has **two solution methods**: the steady state and the chemically significant eigenvalues (CSE). They are chosen in the `MarXus` block of the deck header or on the command line (see Usage).

**Solution method 1: steady state, J·N = F (GO10 eqs. 7, 8), in two versions (GO10 Sec. 3.2)**
- **Intermediate steady state** with an absorbing barrier, at a user-chosen distance below the threshold (default 10 kT).
- **Final steady state** (Olzmann) with physical sinks and no absorbing barrier. Its J also gives the thermal rate coefficients (below).
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

**Solution method 2: phenomenological rate coefficients from the chemically significant eigenvalues** (`Method CSE`, `--method cse`):
- The method is that of Miller, Klippenstein (J. Phys. Chem. A 110, 10528 (2006)) in the formulation of Georgievskii et al. (J. Phys. Chem. A 117, 12146 (2013)).
- It produces species-to-species tables of well ↔ well, well → products, reactant → wells and reactant → products, as MESS does.
- Diagnostics: the eigenvalue separation, the relaxational projections, the loss balance, detailed balance and the precision floor.
- It reproduces a four-well MESS run to a few percent (`validation/ZZAllyl+O2_Gamma_Case2/`).

**Input**
- Decks in the MESS input format (a subset), so existing decks can be used.
- A MarXus keyword block adds ILT parameters to barriers.
- A `MarXus ... End` block in the deck header selects the solution method and its settings (see Usage).
- Exact Eckart tunneling and 1 cm⁻¹ state counting are built from the deck.

**Validation**
- `validation/c2h3_mess_example/` (intermediate steady state) and `validation/c2h3_mess_example_olzmann_eigen/` (thermal rate coefficients of the final steady state) compare H + C₂H₂ ⇌ C₂H₃ at 300–2000 K and 0.1–10 atm with the stored MESS results.
- The thermal rate coefficients agree within ±2.5% from 750 to 1750 K. At 300–1000 K the residual offsets are the known exact vs semiclassical Eckart difference.
- Where both are valid, the intermediate steady state and the thermal rate coefficients of the final steady state agree with each other to 0.01%.
- `validation/ZZAllyl+O2_Gamma_Case2/` reproduces a four-well MESS run: ZZ-allyl + O₂ with two PST channels, Eckart tunneling and an escape sink.
  - The capture rate agrees within 0.5%.
  - Every channel's k∞ agrees within 0.3–0.7% once the tunneling factor is accounted for. MarXus's exact Eckart κ is 2–23% larger than MESS's.
  - The IEPOX + OH share is 1.4–7.7% in MarXus and 1.3–7.2% in MESS.
  - The escape share agrees within 0.5%.
- All validity checks of the master equation (sum rules, detailed balance, limits, solver agreement) are described in `reports/master_equation_validity_checks.md`.

### Numerical library (`src/numeric/`)
- Banded Cholesky (factor once, solve many) and LDLᵀ with Bunch–Kaufman pivoting.
- BiCGSTAB and GMRES.
- Symmetric eigensolvers: Householder + implicit QL, Jacobi, inverse iteration with deflation and shift.
- LAPACK interface (DSYEVD, system OpenBLAS); tridiagonal solvers; Lanczos Γ function.

---

## Master-equation solvers: one master equation, three questions

**All solvers of MarXus solve the same master equation.** They use the same grains, the same collision operator, the same microcanonical rate coefficients k(E) and the same source of chemically activated adducts. They differ only in the **question** they ask of that equation, that is, which solution of it is computed.

There are two solution methods and three solvers in total: the steady-state method with two solvers (intermediate and final steady state) and the CSE method. Their results are therefore **different quantities, not three approximations of one quantity**:

| solver | the question it answers |
|---|---|
| **1. Intermediate steady state** (absorbing barrier, as in SSUMES) | What happens to freshly formed adducts on their first collisional descent? Which fraction redissociates, which decomposes to each product, and which is stabilized? |
| **2. Final steady state** (Olzmann), with its thermal eigenpair | Under continuous formation, once the stabilized population has itself reached a steady state (no net stabilization), what are the yields of all products and sinks? And, from the lowest eigenpair of the same J: how fast does a thermalized adduct react? |
| **3. Chemically significant eigenvalues** (CSE; Miller, Klippenstein; Georgievskii et al.) | Which species-to-species rate coefficients (reactants, wells, products) define a kinetic model that reproduces the master-equation kinetics after collisional relaxation? |

Solvers 1 and 2 give flux coefficients and yields; solver 3 gives phenomenological rate coefficients. Miller and Klippenstein (MK06, pp. 10529, 10531): "application of the steady-state approximation … is virtually always an attempt to equate a phenomenological rate coefficient to a flux coefficient. Sometimes this is a valid approach, and sometimes it is not."

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

The master equation for the grained populations $`\mathbf N`$ of all wells (PO14 eq. 2; O02 eq. 6; GO10 eq. 7) is

```math
\frac{d\mathbf N}{dt} = R\,\mathbf F - \mathbf J\,\mathbf N,
\qquad
\mathbf J = \omega\,(\mathbf I - \mathbf P) + \mathbf K + k_c[\mathrm D]\,\mathbf I .
```

- $`\omega`$ is the collision frequency (Lennard-Jones) and $`\mathbf P`$ the matrix of collisional transition probabilities (exponential down by default). They obey detailed balance, $`P(E' \leftarrow E)\,f^0(E) = P(E \leftarrow E')\,f^0(E')`$, with $`f^0(E) = \rho(E)\,e^{-E/k_BT}`$.
- $`\mathbf K`$ is the microcanonical reaction: $`\sum_r k_r(E)`$ on the diagonal of every well, and the isomerization couplings between wells.
- $`k_c[\mathrm D]`$ is the pseudo-first-order bimolecular sink of a well, for example an escape channel.
- $`R\,\mathbf F`$ is the formation of chemically activated adducts: the rate $`R`$ times the normalized nascent distribution $`\mathbf F`$.

**Properties.** $`\mathbf J`$ has positive eigenvalues when population can leave the network (GO10, text before eq. 12). With $`\mathbf D = \mathrm{diag}\big(\sqrt{f^0}\big)`$, the matrix $`\mathbf S = \mathbf D^{-1}\mathbf J\,\mathbf D`$ is symmetric, and all solvers work with $`\mathbf S`$. The CSE literature writes the same equation as $`d|f\rangle/dt = -\hat G\,|f\rangle + \sum_\nu s_\nu\,|p^{(\nu)}\rangle`$ (G13 eq. 1), with $`\hat G = \mathbf J`$.

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

### Solver 1: intermediate steady state (absorbing barrier)

**Question.** What is the fate of the nascent adducts during their first collisional descent?

**Equation.** Each well gets an absorbing barrier at $`E_{\mathrm{abs}} = E_{\mathrm{thr,min}} - X\,k_BT`$ (default $`X = 10`$, `--barrier-kt`). A molecule transferred below it counts as stabilized and is removed. With $`\mathbf J_{\mathrm{abs}}`$, the operator on the grains above the barriers, and $`R = 1`$:

```math
\mathbf J_{\mathrm{abs}}\,\mathbf N^{s} = \mathbf F,
\qquad
\Phi_r = \sum_E k_r(E)\,N^s(E),
\qquad
\Phi_{\mathrm{stab}} = \sum_{E \ge E_{\mathrm{abs}}} \omega \sum_{E' < E_{\mathrm{abs}}} P(E' \leftarrow E)\,N^s(E)
```

($`\Phi_{\mathrm{stab}}`$ also includes the isomerization flux that arrives below the barrier of the target well, and the part of the source formed below a barrier; `chemical_activation_operator.rs`). The yields obey $`\sum_r \Phi_r + \Phi_{\mathrm{stab}} + \Phi_{\mathrm{sink}} = 1`$, and the apparent bimolecular rate coefficients are $`k(\mathrm{A+B} \to X) = k_\infty\,\Phi_X`$ (PR03 eq. 44).

**Time window.** The solution holds for $`(0.1\,\lambda_F)^{-1} < t < (10\,k_{\mathrm{uni}})^{-1}`$ (O02 p. 3618; SN84), with $`\lambda_F`$ the eigenvalue whose eigenvector has the largest weight in $`\mathbf F`$. In this window the activated population has relaxed, and the stabilized adducts have not yet reacted thermally.

**What it gives for chemical activation:**

- the prompt branching of the nascent adducts: redissociation, chemically activated products, stabilization;
- the falloff of the association, $`k(\mathrm{A+B}\to W) = k_\infty\,\Phi_{\mathrm{stab},W}`$, which is "the rate into the absorbing barrier" (PR03);
- the chemically activated product channels, $`k(\mathrm{A+B}\to P) = k_\infty\,\Phi_P`$.

**What it does not give:**

- the later thermal reaction of the stabilized adducts;
- the thermal rate coefficient.

**Limits:**

- The barrier position is a choice, and for shallow wells the results depend on it. A barrier below the well bottom is refused.
- With a physical bimolecular sink the barrier is "artificial", and "a too low product yield would be predicted" (O02 p. 3617).

### Solver 2: final steady state (Olzmann) and its thermal eigenpair

**Question.** Under continuous formation, after the stabilized population has itself reached a steady state, where does every formed molecule end? And how fast does a thermalized adduct react?

**Equation.** The same $`\mathbf J`$ is used, without a barrier (GO10 eq. 8; PO14 eq. 5):

```math
\mathbf J\,\mathbf N^{s} = R\,\mathbf F,
\qquad
\tilde{\mathbf N}^{s} = \frac{\mathbf J^{-1}\mathbf F}{\sum_i \big(\mathbf J^{-1}\mathbf F\big)_i},
\qquad
k^{ca}_r = \sum_E k_r(E)\,\tilde N^s(E) \quad \text{(GO10 eq. 9)},
```

```math
\Phi_r = \sum_E k_r(E)\,N^s(E) \quad \text{(O02 eq. 10)},
\qquad
\Phi_{\mathrm{sink}} = k_c[\mathrm D] \sum_E N^s(E) .
```

**The thermal eigenpair of the same $`\mathbf J`$** (GO10 eq. 12 and the text after it) is part of this solver, not a method of its own. Its eigenvector is the thermal steady-state population:

```math
\mathbf J\,\tilde{\mathbf n}^{th} = \lambda_1\,\tilde{\mathbf n}^{th},
\qquad
k^{th} = \lambda_1,
\qquad
k^{th}_r = \sum_E k_r(E)\,\tilde n^{th}(E),
\qquad
k_{\mathrm{uni}} = \sum_r k^{th}_r + k_c[\mathrm D] .
```

**How it is computed and reported:**

- MarXus reports the eigenvector average $`k_{\mathrm{uni}}`$, with $`\lambda_1`$ beside it as the sum-rule check.
- The default solver, inverse iteration, factors $`\mathbf S + \sigma\mathbf I`$ and solves $`(\mathbf S+\sigma\mathbf I)\,\mathbf x = \mathbf u`$ repeatedly. Each step is a steady-state solve with the previous distribution as the source.
- For one well with one entrance channel, the association follows by detailed balance: $`k(\mathrm{A+B}\to W) = k_{\mathrm{uni}}\;k_{\infty,\mathrm{assoc}}/k_{\infty,\mathrm{diss}}`$.

**Time scale.** The solution is reached when the experimental time is distinctly longer than $`1/k_{\mathrm{uni}}`$ (O02 p. 3618). At that point "there is no more net stabilization; the stabilization reservoir is filled up, and time-independent energy distributions have been established" (GO10 p. 12295).

**What it gives for chemical activation:**

- the complete fate of all formed molecules, including those that react after stabilization;
- with physical sinks (for example O₂ addition or an escape channel), the yields that a continuously fed system shows;
- non-thermal sources and consecutive activation (PO14);
- $`k^{ca}`$ (activated adduct) and $`k^{th}`$ (thermalized adduct) from one and the same $`\mathbf J`$, directly comparable (GO10).

**Link to solver 1.** A physical sink in the window $`0.01\,\omega > k_c[\mathrm D] > 10\,k_{\mathrm{uni}}`$ gives the absorbing-barrier yield within 10% (O02 p. 3618). The sink acts as a physically defined absorbing barrier.

**Limits:**

- Without a sink, and with a single exit such as back to A + B, every molecule eventually leaves through it. In O02's words (text after eq. 13), "one trivially has … $`\Phi_2 = 1`$".
- For deep wells at low T, $`\mathbf J\,\mathbf N = \mathbf F`$ then becomes numerically singular, for example C₂H₃ at 300–500 K. The thermal eigenpair is still obtained (shifted inverse iteration). For the stabilization, use solver 1 or 3.

### Solver 3: chemically significant eigenvalues (CSE)

**Question.** Which phenomenological rate coefficients between species reproduce the kinetics after relaxation?

**Equations** (G13; MK06). All eigenpairs of $`\hat G = \mathbf J`$ are computed, $`\hat G\,f^{(\lambda)} = \Lambda_\lambda\,f^{(\lambda)}`$. For $`N`$ wells, the $`N`$ lowest eigenvalues are chemically significant and must be well separated from the relaxation eigenvalues, $`\Lambda_N \ll \Lambda_{N+1}`$ (MK06 eq. 19). With $`Q_i = \sum_{E \in i} f^0(E)`$ and $`p^{(\nu)}_\lambda = \sum_E f^{(\lambda)}(E)\,k_{\to\nu}(E)`$ (eq. 15):

```math
M_{i\lambda} = Q_i^{-1/2} \sum_{E\in i} f^{(\lambda)}(E) \quad \text{(eq. 25)},
\qquad
k_{j\to i} = -\sqrt{Q_i/Q_j}\;\big(M\Lambda M^{-1}\big)_{ij} \quad \text{(eq. 27)},
```

```math
k_{i\to\nu} = Q_i^{-1/2} \sum_{\lambda \le N} \big(M^{-1}\big)_{\lambda i}\,p^{(\nu)}_\lambda \quad \text{(eq. 30)},
\qquad
k_{R\to i} = \frac{\sqrt{Q_i}}{Q_R} \sum_{\lambda \le N} M_{i\lambda}\,p^{(R)}_\lambda \quad \text{(eq. 28)},
```

```math
k_{R\to\mu} = \frac{1}{Q_R} \sum_{\lambda > N} \frac{p^{(\mu)}_\lambda\,p^{(R)}_\lambda}{\Lambda_\lambda} \quad \text{(eq. 21)}.
```

**Validity.** The rate coefficients describe the kinetics for $`t \gg 1/\Lambda_{N+1}`$, after relaxation. They exist only while the chemically significant eigenvalues stay separated from the relaxation eigenvalues. Otherwise G13 (Sec. IV) merges species. MarXus warns when $`\Lambda_N/\Lambda_{N+1} > 0.1`$; merging is not implemented.

**What it gives for chemical activation:**

- $`k_{R\to i}`$ (eq. 28): stabilization into well $`i`$;
- $`k_{R\to\mu}`$ (eq. 21): the chemically activated, well-skipping products, carried by the relaxation modes;
- the thermal well → well and well → product rate coefficients;
- together, a complete mechanism for kinetic models.

**What it does not give:**

- yields of a continuously fed system;
- non-thermal sources;
- results when the separation fails.

### Side by side

| | 1. intermediate steady state | 2. final steady state (+ thermal eigenpair) | 3. CSE |
|---|---|---|---|
| master equation | same J, same F | same J, same F | same J, same F |
| mathematical problem | linear system on the grains above the absorbing barriers | linear system on all grains; lowest eigenpair of the same J | all eigenpairs of J |
| stabilized adducts | removed (counted as stabilized) | stay in the well and react thermally | a chemical eigenmode per well |
| time scale | $`(0.1\lambda_F)^{-1} < t < (10k_{\mathrm{uni}})^{-1}`$ | $`t \gg 1/k_{\mathrm{uni}}`$ | $`t \gg 1/\Lambda_{N+1}`$ |
| output | yields Φ, $`k_\infty\Phi_X`$ | $`k^{ca}`$, yields Φ, sink yields; $`k_{\mathrm{uni}}`$, $`\lambda_1`$, $`k^{th}_r`$ | species-to-species rate coefficients |
| kind of quantity | flux coefficient | flux coefficient; thermal rate coefficient | phenomenological rate coefficient |
| needs | barrier distance (choice) | a sink or exit for a well-conditioned J·N = F | eigenvalue separation |
| typical use | association falloff, prompt product branching | yields with physical sinks (continuous formation), consecutive activation, thermal $`k_{\mathrm{uni}}`$ | rate coefficients for a kinetic mechanism |

### Where the three must agree: identities and measured agreement

- **One well.** The CSE well → product rate coefficient equals the eigenvector-average $`k_{\mathrm{uni}}`$ of solver 2 (test `a_single_well_gives_the_eigenvector_average_as_its_rate_coefficient`, to 10⁻⁸).
- **Well separation and continuous formation.** Steady-state yields equal time-integrated single-injection yields (Vereecken et al., J. Chem. Phys. 106, 6564 (1997), Table I). The CSE R → product and R → well rate coefficients equal the fractions accumulated on the relaxation time scale (G13 p. 12153).
- **H + C₂H₂ ⇌ C₂H₃, one well** (`validation/c2h3_mess_example*`):
  - At 300–1000 K, solver 1 (10 kT) and the thermal association of solver 2 agree to 0.01%.
  - Above about 1500 K, solver 1 fails because of the barrier distance (0.1 atm, 10 kT): −13% at 1500 K, −42% at 1750 K, no result at 2000 K.
  - Solver 2 stays within ±2.5% of the MESS (CSE) result from 750 to 1750 K.
- **ZZ-allyl + O₂, four wells** (`validation/ZZAllyl+O2_Gamma_Case2/`):
  - Solver 3 reproduces the MESS species tables to a few percent.
  - The long-time IEPOX + OH share of solver 2 differs from MESS's long-time fate by +7 … +12% with the exact Eckart tunneling, and by −5.1 … −5.9% with the MESS tunneling model. The escape share agrees within 0.5%.
  - The apparent $`k(\mathrm R\to \mathrm{IEPOX+OH})`$ of solver 1 is 15% above MESS's R → P5. These are different quantities: the prompt formation is assigned by the barrier in one and by the eigenvalue splitting in the other.

---

## Usage example: master equation from a deck

```
cargo run --release --example chemical_activation_from_deck -- deck.inp REACTANT \
    --method steady-state --steady-state both --eigen-solver inverse
```

The solution method and its settings can be given in the deck header, in a `MarXus ... End` block (a MarXus extension of the MESS format). Every keyword is optional; the command-line option in the last column overrides it:

```
MarXus
  Method                              SteadyState        ! SteadyState | CSE
  SteadyState                         Both               ! Intermediate | Final | Both
  AbsorbingBarrierBelowThreshold[kT]  10                 ! intermediate steady state
  EigenSolver                         InverseIteration   ! InverseIteration | FullDecomposition | Lapack
  SumRuleTolerance                    1.5e-2             ! thermal eigenpair of the final steady state
End
```

| deck keyword | option | meaning (default) |
|---|---|---|
| `Method` | `--method steady-state\|cse` | solution method (steady state) |
| `SteadyState` | `--steady-state intermediate\|final\|both` | versions of the steady-state method (both); the final one includes the thermal rate coefficients |
| `AbsorbingBarrierBelowThreshold[kT]` | `--barrier-kt X` | absorbing barrier of the intermediate steady state, in k_BT below the lowest threshold (10) |
| `EigenSolver` | `--eigen-solver inverse\|full\|lapack` | thermal eigenpair of the final steady state (inverse iteration); CSE needs all eigenpairs (LAPACK; `full` also possible, inverse iteration refused) |
| `SumRuleTolerance` | `--sum-rule-tolerance X` | warning threshold of \|λ₁ − k_uni\|/k_uni (1.5e-2) |
| – | `--tunneling exact-eckart\|mess-eckart` | Eckart transmission model: exact Eckart (default); `mess-eckart` only for comparison with MESS (see Tunneling model below) |

A setting that the selected solution does not use is reported as a note in the output, not refused. `eigenvalue` is neither a method nor a steady-state version: the thermal eigenpair belongs to the final steady state.

The output tables begin with comment lines that explain the quantities, with their references. The source file `examples/chemical_activation_from_deck.rs` documents all columns.


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
- **Remaining difference.** The IEPOX + OH share is then 5.1–5.9% below MESS. With the exact Eckart it is 7–12% above MESS. The remainder is collisional, not tunneling.

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

The current version has no crate dependencies.

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
- Well merging in the CSE method when chemical eigenvalues approach the relaxation ones (Georgievskii et al. 2013, Sec. IV).
- Treatment of the stepladder model in the eigenvalue analysis when its step spans several grains (independent sub-equations).
- Excited electronic states in the partition functions and state counts.
- A general equilibrium-constant routine (thermochemistry).
- Microcanonical Variational TST (μVTST).
- Canonical Variational TST (CVTST).
- Submerged barrier with a pre-reaction vdW complex: **μ-canonical, J-resolved 2-TST treatment**.

---

## License
This project is licensed under the **GNU General Public License v3.0 (GPL-3.0)**.
