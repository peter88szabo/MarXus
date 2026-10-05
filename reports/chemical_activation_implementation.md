# Multiwell chemical-activation solver: implementation report

**MarXus**, 2026-10-05. Everything described here is **uncommitted**.  
This report implements `chemical_activation_code_design.md`. Background:
`chemical_activation_three_approaches.md` and `collision_kernel_detailed_balance_and_normalization.md`.

**Test status at the time of writing:** `cargo test` passes 89 library tests and 18 binary tests,
and all targets build. The production code of the new files was written test-first: each test was
seen to fail before its implementation. The one exception is noted under §2.1.

---

## 0. Abbreviations of the cited literature (also used in the code comments)

| Abbreviation | Reference |
|---|---|
| O91 | Olzmann, Gebhardt, Scherzer, Int. J. Chem. Kinet. 23, 825 (1991) |
| O02 | Olzmann, Phys. Chem. Chem. Phys. 4, 3614 (2002) |
| GO10 | González-García, Olzmann, Phys. Chem. Chem. Phys. 12, 12290 (2010) |
| PO14 | Pfeifle, Olzmann, Int. J. Chem. Kinet. 46, 231 (2014) |
| T77 | Troe, J. Chem. Phys. 66, 4758 (1977) |
| R19 | Robertson (ed.), Comprehensive Chemical Kinetics 43 (2019) |
| PR03 | Pilling, Robertson, Annu. Rev. Phys. Chem. 54, 245 (2003) |
| CD07 | Carstensen, Dean, Comprehensive Chemical Kinetics 42 (2007) |
| DGP86 | Davies, Green, Pilling, Chem. Phys. Lett. 126, 373 (1986) |
| F73 | Forst, Theory of Unimolecular Reactions (1973) |

### Citations checked against the papers in this session

The places were checked in the text extracts of the PDFs.

| Statement | Place |
|---|---|
| ME **J = ω(I − P) + K + k_c[D]·I** | PO14 eq. 2; O02 eq. 6; GO10 eq. 7 |
| Steady state **N = R·J⁻¹F** | PO14 eq. 5; O02 eq. 7 |
| Normalized distribution | GO10 eq. 8; O02 eq. 8 |
| k^ca = Σ(K_r Ñ) | GO10 eq. 9; PO14 eq. 6 |
| Yield Φ = Σ(K J⁻¹F) | O02 eq. 10 |
| Relative sink yield | GO10 eq. 10 |
| k^th = λ₁ | GO10 eq. 12 |
| Thermal source W(E−E₀)·e^(−(E−E₀)/kT) | PO14 eq. 7; O02 eq. 11 |
| k = W/(hρ) | PO14 eq. 9 |
| Convolution source | PO14 eq. 8 |
| Consecutive steps | PO14 eqs. 10–11 (convolution), 12–13 (shift) |
| "energy-independent reactive cross section ⇔ vibrational PST without angular-momentum restrictions" | PO14 p. 235 |
| "the output of the first master equation governs the input of a second one" | PO14 p. 237 |
| Intermediate steady state "implemented by introducing a lower absorbing barrier"; final steady state = "no more net stabilization" | GO10 p. 12295 |
| "lower absorbing barrier in the collisional deactivation cascade at energies below the lowest reaction threshold" | O02, text after eq. 13 |
| Final steady state after ≈ 0.1/k_uni; "R1 = R2 and hence Φ2 = 1" | O02, text after eq. 13 |
| "The absorbing boundary is usually placed about 10 k_BT below the reaction threshold" | PR03, Sec. 2.4 |
| "for example, 10 kT below E0" | CD07 p. 125 |
| Stepladder step ΔE_SL "represents the average amount of energy transferred in down collisions"; eq. 16 relates it to ⟨ΔE⟩ = ΔE_SL·tanh(ΔE_SL/2F_E·kT) | GO10, text **before** eq. 16 (the code comments were corrected from "eq. 16 and text") |
| Symmetrization S_ij = M_ij (b_j/b_i)^½, S = F⁻¹MF, similarity transform | R19 eqs. 5.74–5.77, p. 319 |
| ILT, plain-Arrhenius form: k(E) = A∞·ρ(E−E∞)/ρ(E) | R19 eq. 7.28, p. 428 |
| ILT, generalized Arrhenius form credited to DGP86 | R19 p. 428 (refs [24] and [69] = DGP86) |
| Recombination ILT; "more complex functional forms may be transformed analytically"; with E∞ < 0, k(E) is non-zero below the threshold | DGP86 eq. 2 and p. 377 |

---

## 1. New and changed files

`mod.rs` contains only `pub mod` lines.

### 1.1 New production files (`src/masterequation/`)

| File | Content |
|---|---|
| `collision_kernels.rs` | **Exponential down**: Robertson eqs. 4.4, 4.7, 4.11, 4.16, back substitution from the top grain, with low-energy reduction factors on the lower grain of each pair. **Stepladder**: O91 eqs. 13–18. Both exactly detailed-balanced. (Earlier phase.) |
| `chemical_activation_network.rs` | Data model: `Well`, `Channel`, `ChannelDestination::{Products, Well}`, `LennardJonesPair`, `EnergyTransferParameters`, `CollisionModel`, `SteadyState`, `AbsorbingBarrier`, `Conditions`, `ChemicalActivationOptions`, `validate()`. ⟨ΔE_down⟩(T) = ⟨ΔE_down⟩(T_ref)·(T/T_ref)ⁿ with an explicit T_ref. |
| `chemical_activation_operator.rs` | Assembly of J. Ordering, Boltzmann weights, stabilization bookkeeping and the isomerization check: §3. |
| `chemical_activation_steady_state.rs` | `project_source`: normalizes F; source weight below the barrier counts as directly stabilized. `solve_steady_state`: §3.3. |
| `chemical_activation_sources.rs` | Thermal source from k_reverse (PO14 eqs. 7 and 9); multi-channel thermal entrance on the absolute energy scale; thermal source from W‡ (PO14 eq. 7, O02 eq. 11); single grain; convolution (PO14 eq. 8); shift (PO14 eqs. 12–13); thermal distribution; mean energy. Weight falling outside the receiving grid is an error above 10⁻¹⁰. |
| `chemical_activation_observables.rs` | Φ_r, k_r^ca, Φ_stab, Φ_sink, population fractions, normalized per-well distributions, ⟨E⟩, mass balance. |
| `chemical_activation_driver.rs` | Loop over T (outer) and p (inner); source per T; residual and mass-balance checks; results table. |
| `consecutive_activation.rs` | ME₁ → ñ₁ˢˢ → f₂ by convolution (PO14 eqs. 10–11) or shift (PO14 eqs. 12–13) → ME₂ at every (T, p). |
| `chemical_activation_from_mess_input.rs` | Input deck (MESS format) → network: §5. |

### 1.2 Changed files

- **`numeric/banded_solvers.rs`**: new `SymmetricBandMatrix` and `BandedCholeskyFactor`. These are a row-oriented Cholesky–Crout with factor-once / solve-many, plus 2 tests. The old function there had no tests.
- **`barrierless/ilt/ilt_barrierless.rs`**: rewritten; see §4.
- **`mess_input.rs`**:
  - lists for T and p;
  - `ExponentCutoff`, `Escape`, well order, `InverseLaplaceTransform`, tunneling flag, kJ/mol;
  - legacy builder removed;
  - 6 new tests.
- **`microcanonical_builder.rs`**: state counting on a plain (grains, ΔE) grid (`rrho_density_of_states`, `rrho_sum_of_states`, `transition_state_sum_of_states`). The legacy network builder and provider were removed. Comments naming another program were replaced by literature (Beyer–Swinehart; F73 §4.5).
- **`collisional_relaxation.rs`**: `lennard_jones_collision_frequency_s_inv(σ, ε, μ, T, p)` (T77 eqs. 3.1–3.3). The legacy wrapper, α(1000 K) helper and band helper were removed; the tests were ported.
- **`constants.rs`**: added `AMU_TO_KG` (CODATA 2018), now shared with `collisional_relaxation.rs`.

### 1.3 Examples

- New `examples/chemical_activation_from_deck.rs`: runs any deck through the new code, intermediate and final steady state, and prints the table.
- New `examples/c2h3_chemical_activation.inp`: the data of the former C₂H₃ single-well example, as an input deck. See §6.
- `examples/mess_parse_zzallyl_o2.rs`: now only lists what the parser read. It shows the barrier flags core / ILT / tunneling.

---

## 2. Data model and options (`chemical_activation_network.rs`)

- **Grid.** All wells share one grain width ΔE. Well w has an integer `bottom_offset_grains`: grain i lies at the absolute energy (i + offset)·ΔE. `aligned_grain(from, i, to)` maps equal absolute energies exactly, with no rounding.
- **Collision model.**
  - `ExponentialDown{cutoff_in_mean_down}`: band = ⌈cutoff·⟨ΔE_down⟩/ΔE⌉.
  - `Stepladder`: ΔE_SL = ⟨ΔE_down⟩(T), rounded to whole grains.
- **Steady state.**
  - `Final`: no barrier.
  - `Intermediate{barrier}`: `BelowLowestThreshold{kt_multiple}` (default 10, from PR03 and CD07) or `AtGrains(Vec)`.
- **Sink.** `bimolecular_sink_s_inv` is k_c[D] per well. It is energy independent, as in PO14 p. 234 ("because little is known about this energy dependence … the thermal rate coefficient is used").

### 2.1 TDD note

This file's data model was written before its 4 tests. The tests were then checked by mutation: the alignment sign, the self-isomerization check, the threshold detection and the sign of the temperature exponent were each broken, and all 4 tests failed. All other files were developed RED → GREEN.

---

## 3. Operator and solution

### 3.1 Assembly of J (`chemical_activation_operator.rs`)

For a retained state c = (w, j), the column of J is:

- **diagonal**: ω_w·Σ_{t≠j} P_w(t|j) + Σ_r k_r(E_j) + k_c[D]_w;
- **collisions**: −ω_w·P_w(t|j) into retained grains t of the same well;
- **isomerization**: −k_{w→w'}(E_j) into the grain of w' at the same absolute energy;
- **absorbed targets**: flux into an absorbed grain, by collision or by isomerization, is recorded as stabilization into that well.

The diagonal is the sum of the stored transitions, which equals ω(1 − P(j|j)). It therefore makes the **column sum equal the loss out of the network exactly**, which gives an exact mass balance.

**Kernel on the full grid.** The kernel is always computed on the complete grid of the well. Probabilities out of a retained grain thus stay normalized when grains below the barrier are removed.

**Barrier.** The default is grain(lowest threshold) − round(10·kT/ΔE). The lowest threshold is the first grain with any open channel, isomerization included.

**Errors** (never silent):
- isomerization that targets energies outside the target grid;
- a well without open channels when the barrier is defined relative to its threshold;
- an explicit barrier beyond the grid.

**Ordering.** States are ordered by absolute energy, then by well. All non-zero elements then lie within about n_wells × (kernel band) of the diagonal. Isomerization couples neighbours at equal energy (tested).

**Weights.** The log Boltzmann weights ln f = ln ρ_w − E_abs/kT live on the **absolute** energy scale, which is required for detailed balance between wells.

### 3.2 Detailed balance

- Collisions: P(t|j)·f_j = P(j|t)·f_t for both kernels.
- Isomerization: ρ_a·k_ab = ρ_b·k_ba = W‡/h.

Together these give J_rc·f_c = J_cr·f_r, which is tested for all four combinations of options. `isomerization_detailed_balance(network)` reports, for each pair of connected wells, the largest relative deviation of ρ_a·Σk_ab from ρ_b·Σk_ba.

### 3.3 Steady state (`chemical_activation_steady_state.rs`)

- **Symmetrization.** S = D⁻¹JD with D = diag(√f), i.e. S_rc = J_rc·(f_c/f_r)^½ (R19 eq. 5.75). The weights are taken relative to their maximum and handled in log form, so they neither underflow nor overflow. The right-hand side F/D is also formed from logarithms.
- **Asymmetry check.** max |S_rc − S_cr| / max(|S_rc|, |S_cr|) is reported. The Cholesky path **refuses** asymmetry above 10⁻⁸ and names the cause: isomerization rates that violate detailed balance.
- **`BandedCholesky`.** Bandwidth = max |r − c| over the non-zeros. Lower band = ½(S_rc + S_cr). Factor once, solve, then up to 3 steps of iterative refinement on the **unsymmetrized** J·N = F, kept only while the residual decreases.
- **`BiCgStab{tol, max_iter}`.** Jacobi-preconditioned BiCGSTAB on D⁻¹JD (`numeric/krylov.rs`). It works without symmetry.
- **Residual.** ‖F − J·N‖₂/‖F‖₂ is evaluated independently with the unsymmetrized J.
- **Driver checks.** The driver refuses a row when the residual or |mass balance − 1| exceeds the tolerance.

### 3.4 Observables (`chemical_activation_observables.rs`)

F is normalized and R = 1.

- **Yield** of a product channel: Φ_r = Σ_i k_r·N_i (O02 eq. 10).
- **Rate coefficient**: k_r^ca = Σ_i k_r·Ñ_i, with Ñ normalized over the well (GO10 eqs. 8–9). It is NaN if the well holds no population.
- **Internal flux**: for an isomerization channel, Σ_i k·N_i is the gross internal flux.
- **Sink**: Φ_sink = k_c[D]·Σ_i N_i.
- **Stabilization**: Φ_stab = Σ_s N_s·(rate into the absorbing region) plus the source fraction formed below the barrier.
- **Mass balance**: Σ Φ_products + Σ Φ_stab + Σ Φ_sink. It is returned and checked.
- **Total loss rate**: k_tot = 1/Σ_w Σ_i N_i is the total loss rate coefficient of the intermediates (column `k_tot`).

### 3.5 Results table (`write_results_table`)

The columns are:

`T[K], P[Torr], Phi(well:channel)…, Phi_stab(well), Phi_sink(well)…, k_ca(well:channel)[1/s]…, k_tot[1/s], fpop(well)…, <E>(well)[cm-1]…, mass_balance, residual`.

Commas in names are replaced by `;`. The table carries the same information as SSUMES carate's CSV (fractions, stabilization, rates, k_tot, internal rates, populations) with Olzmann's definitions.

### 3.6 Consecutive activation (`consecutive_activation.rs`)

- **Source.** `consecutive_source(first_result, second_network, step, T)` takes ñ₁ˢˢ of `from_well` and produces f₂ in `to_well`, by:
  - `Convolution{partner_density_of_states}`: thermal partner distribution ⊗ ñ₁ˢˢ, placed at E₀ = −RE (PO14 eqs. 10–11);
  - `Shift{partner_mean_energy_cm1}`: ñ₁ˢˢ(E + RE − ⟨E_B⟩) (PO14 eqs. 12–13).
- **Restriction.** RE ≤ 0 is required: PO14 treats exothermic steps with an energy-independent cross section.
- **Chaining.** `run_consecutive_activation` runs ME₁ for all (T, p), then ME₂ at each condition with that condition's f₂.

---

## 4. ILT check and rewrite (`barrierless/ilt/ilt_barrierless.rs`)

### 4.1 Theory

**Dissociation**, k∞ = A·(T/T_ref)ⁿ·e^(−E∞/kT) in s⁻¹:

k(E)ρ(E) = A·β_refⁿ/Γ(n) ∫₀^{E−E∞} ρ(E−E∞−x)·x^(n−1) dx, valid for n > 0.

The case n = 0 is the δ kernel, k = A·ρ(E−E∞)/ρ(E) (R19 eq. 7.28).

**Association**, k∞ in cm³ s⁻¹, with N_P the convolved rovibrational density of the fragments:

k(E)ρ(E) = A·C′(μ)·β_refⁿ/Γ(n+3/2) ∫₀^{E−E_th} N_P(E_P)·(E−E_th−E_P)^(n+1/2) dE_P, valid for n > −3/2,

with E_th = ΔE₀ + E∞. This is DGP86 eq. 2 for n = 0, generalized through the Laplace pair L{x^(ν−1)/Γ(ν)} = β^(−ν), as DGP86 p. 377 allows.

C′(μ) = (2πμ/h²)^{3/2}, so that C′·(kT)^{3/2} is the translational partition function per cm³. With energies in cm⁻¹, C′ = 3.2433×10²⁰·μ^{3/2} cm⁻³ (μ in amu). It is computed from `PLANCK_SI`, `CLIGHT_SI` and `AMU_TO_KG`.

**Validity.** n ≥ 0 for dissociation, n > −3/2 for association, and **E∞ ≥ 0** (DGP86 p. 377). Violations are errors.

### 4.2 Problems found in the former code

1. **The kernel was evaluated at x = 0 by trapezoid quadrature.** x^(n−1) (dissociation, 0 < n < 1) and x^(n+1/2) (association, n < −1/2, i.e. typical recombination exponents such as −0.5 or −1) are infinite there, so **W(E) = ∞**. For dissociation with n ≤ 0, Γ(n) gave NaN.
2. **ρ sampling.** ρ was sampled at the nearest grain under a 4000-point trapezoid for every energy. This was slow and inconsistent with the grain grid.
3. **No checks** on n, E∞ or A.
4. **C′ was not computed**: it was a caller-supplied `bimol_prefactor_cm`. No code called the functions, so nothing produced wrong numbers in practice.

### 4.3 New numerics

The method follows MESMER's ILT, read on your request. No code was copied and the program is not named in the code.

- **Grain-integrated kernel.** w_j = ∫ over [jΔE, (j+1)ΔE) of x^(ν−1)/Γ(ν) dx = [(j+1)^ν − j^ν]·ΔE^ν/Γ(ν+1). This is finite for every ν > 0, and ν → 0 gives the δ kernel automatically (w₀ = 1). For j ≥ 1 it is evaluated as j^ν·expm1(ν·ln1p(1/j)) to avoid cancellation.
- **Discrete convolution.** Σ_{j=0}^{i} w_j·ρ[i−j]. This is **exact for MarXus's piecewise-constant density** convention: ρ[m] holds the states in ((m−1)ΔE, mΔE], and the term j = i captures the ground state ρ[0].
- **Output.** W(ε) = h·k·ρ at ε = E − E_th. A channel then gets k(E) = W(E−E_th)/(hρ(E)), like a transition state.

### 4.4 Tests

| Test | Result |
|---|---|
| Dissociation, n = 1 | exact: W = h·A·β_ref·G(E−E∞) |
| Dissociation, n = 0 | exact: R19 eq. 7.28 |
| Dissociation, n = 0.5 (singular kernel); forward Laplace transform | within 1% of A(T/T_ref)ⁿe^(−E∞/kT) at 300, 600, 1000 K |
| Association, n = −1 (singular kernel); forward transform via detailed balance with C′(kT)^{3/2}Q_frag | within 1% at 300, 600, 1000 K |
| C′(μ)·(kT)^{3/2} vs (2πμk_BT/h²)^{3/2} in SI | agree to 10⁻⁶ |
| n < 0 (dissociation), n ≤ −3/2 (association), E∞ < 0, A ≤ 0 | rejected |

The residual 1% of the forward checks is the test's own sampling of e^(−E/kT) at grid points. It is of order βΔE with ΔE = 2 cm⁻¹; the inversion itself is exact for the piecewise density.

---

## 5. Input decks (MESS format) → network

### 5.1 Parser extensions (`mess_input.rs`)

- **T and p lists.** `TemperatureList[K]` is now a list. `PressureList[torr|atm|bar]` is a list converted to Torr; 1 bar = 10⁵·760/101325 Torr. Both were single values before.
- **Exponential down.** `ExponentCutoff` is read.
- **Sink.** The well block `Escape … PseudoFirstOrderRateConstant[1/sec]` is read and becomes k_c[D] (`well_escape_rate_s_inv`).
- **Well order.** The order of the deck is kept (`well_order`); it was a `HashMap` before, which made the order nondeterministic. A duplicate well is an error.
- **Tunneling.** `has_tunneling` is set when a barrier has a `Tunneling` block.
- **Units.** Energies in kJ/mol are accepted (1 kcal = 4.184 kJ).
- **ILT keyword.** New MarXus keyword block inside the barrier's **RRHO** block (your decision: a keyword in the deck, both directions):

```
    InverseLaplaceTransform
      Direction                    Association      ! or Dissociation
      PreExponential[cm^3/s]       6.0e-12          ! [1/s] for Dissociation (unit must match)
      TemperatureExponent          -0.5
      ReferenceTemperature[K]      298.0
      ActivationEnergy[kcal/mol]   0.0              ! [1/cm], [kJ/mol] also accepted
    End
```

The block must sit inside RRHO, because a barrier ends with the `End` of its RRHO block. `InverseLaplaceTransform` is registered as a block starter. A unit that does not match the direction, an unknown keyword or a missing value is an error.

### 5.2 Adapter rules (`chemical_activation_from_mess_input.rs`)

**Grid**
- **Grain.** ΔE comes from the settings, otherwise EnergyStepOverTemperature·k_B·(lowest T of the deck), the finest grain.
- **Top.** The common top is from the settings, otherwise ModelEnergyLimit.
- **Rounding.** **Every deck energy is rounded once to absolute grains**, and thresholds are integer differences. The two directions of an isomerization therefore use the same TS grain, and ρ_a·k_ab = ρ_b·k_ba holds **exactly** (tested < 10⁻¹²).
- **Why the old builder broke this.** It rounded E_TS − E_well separately for each well, so the two sides could be one grain apart and detailed balance failed.
- **Grid start.** Each well's grid starts at the **first grain that contains states**. With classical rotors G(0) = Q′·0^{r/2} = 0, so the ground-state grain is empty. An empty grain would break the ρ_t/ρ_j ratios of the kernels.

**Rate coefficients**
- **Tight TS.** k = W‡(E−E₀)/(hρ). The symmetry factor and ground electronic degeneracy enter W‡ and ρ as g_e/σ (F73 §4.5).
- **ILT association.** ρ_AB(E) = Σ ρ_A(E′)ρ_B(E−E′)·ΔE from the two fragments of the Bimolecular species. μ comes from the fragment geometries' atomic masses. E_th = E(asymptote) + E∞.
- **ILT dissociation.** The well's own density is used, with E_th = E(well) + E∞. A threshold below the asymptote is an error.

**Refusals** (errors, never silent)
- **PST core without an ILT block**: the barrierless module is not yet connected (your instruction).
- **Barrier with `Tunneling`**: refused unless `ignore_tunneling` is set (§7).
- ILT between two wells.
- A barrier between two non-wells.
- Missing ExponentCutoff, Factor, Power or LJ data.

**Collisions**
- Lennard-Jones: σ = (σ₁+σ₂)/2, ε = √(ε₁ε₂)·1.438776877 K, μ = m₁m₂/(m₁+m₂) (T77 Sec. III).
- Energy transfer: ⟨ΔE_down⟩(T) = Factor·(T/300 K)^Power.
- Collision model: exponential down with the deck's ExponentCutoff.

**Entrance.** Every channel to the deck's `Reactant` becomes an entrance channel. Several of them are weighted by W‡·e^(−E_abs/kT) on the absolute scale.

---

## 6. First run: C₂H₃ (`examples/c2h3_chemical_activation.inp`)

**System.** H + C₂H₂ → C₂H₃* through the tight TS B1 (+4.42 kcal/mol); the well lies at −34.44 kcal/mol.

**Conditions.** He-like bath: ε = 6.95 / 292 cm⁻¹, σ = 2.55 / 4.36 Å, m = 4 / 75 amu, Factor 200 cm⁻¹, Power 0.85, cutoff 15. T = 1000 K; ΔE = 0.1·kT = 69.5 cm⁻¹; top at 400 kcal/mol, giving 2185 grains. Runtime ≈ 3 s in release mode.

**Provenance of the data.** The data are those of the former example: geometries, frequencies, energies and the collision model, the latter recovered from `HEAD:examples/steadystate_me_c2h3_input.txt`. The **C₂H₂ geometry (r_CC = 1.203 Å, r_CH = 1.062 Å) is my input**, because the deck format derives B from a geometry. It gives B ≈ 1.18 cm⁻¹; the old code had hard-coded 1.176 cm⁻¹.

Results, from `cargo run --release --example chemical_activation_from_deck`:

| Steady state | P (Torr) | Φ(back, B1) | Φ_stab | k^ca (s⁻¹) | k_tot (s⁻¹) | ⟨E⟩ (cm⁻¹) | residual |
|---|---|---|---|---|---|---|---|
| intermediate | 76 | 0.98834 | 0.01166 | 9.312e8 | 9.422e8 | 11343 | 7e-17 |
| intermediate | 760 | 0.93911 | 0.06089 | 1.600e9 | 1.704e9 | 11604 | 3e-16 |
| intermediate | 7600 | 0.77269 | 0.22731 | 3.284e9 | 4.250e9 | 11959 | 1e-15 |
| final | 76 / 760 / 7600 | 1 | 0 | **2.567245e5 (all)** | 2.567245e5 | **2905.4 (all)** | ≤ 2e-11 |

**Interpretation of the final steady state.** When the entrance channel is the only exit and the reactants are thermal, F ∝ K·f with f = ρ·e^(−E/kT). For a normalized, detailed-balanced kernel (I − P)·f = 0, so J·f = K·f ∝ F. **The final steady state is chemical equilibrium at every pressure**, and k^ca equals the canonical average of k(E), i.e. the high-pressure rate coefficient.

This is O02's "R1 = R2, Φ2 = 1" limit. The pressure independence to 7 digits confirms that both kernels are exactly normalized and detailed-balanced. It is now a test for both kernels at 1 and 760 Torr (`final_steady_state_through_the_only_channel_is_the_equilibrium_distribution`).

---

## 7. Open items (decisions or papers needed)

> **Update (same day):** items 1–3 below are resolved: Eckart tunneling, the B12 ILT block, and graining on 1 cm⁻¹ cells. Details are in `tunneling_ilt_and_energy_graining.md`. The C₂H₃ numbers of §6 were also updated there (cell-based: Φ_stab −2%; k∞ now within 0.08% of canonical TST instead of +6.3%).

1. **Tunneling in k(E).** 8 of the 9 barriers of the Case1 deck have `Tunneling Eckart`; these are H-shifts with imaginary frequencies up to 2659 cm⁻¹. At 300 K tunneling changes k(E) by orders of magnitude.
   - MarXus has Eckart transmission probabilities (`tunneling.rs`: `tunprop1`, `tunprop2`), but only inside a canonical κ(T) integral.
   - The microcanonical form, W‡_tun(E) = ∫ ρ‡(ε)·P_tun(E−ε) dε, is Miller, J. Am. Chem. Soc. 101, 6810 (1979) (in `papers/`). The Eckart P(E) is Johnston & Heicklen (1962) (in `papers/Tunneling/`).
   - Until it is implemented, such decks run only with `ignore_tunneling = true`.
   - **Decision needed:** implement it next?
2. **The Case1 deck needs an `InverseLaplaceTransform` block** in barrier B12 (R → G2, PST core) **with your k_rec∞(T) parameters.** I did not invent any.
3. **Grain size.** MarXus's ρ[i] counts states in ((i−1)ΔE, iΔE]. Discrete convolutions such as ρ_AB = ρ_A ⊗ ρ_B therefore represent energy ≈ (i−1)ΔE at index i, half a grain from the centre. Thresholds and Boltzmann factors carry errors of order βΔE.
   - Example: in the adapter test, the ILT association channel opens 2 grains above the asymptote, because both fragments' classical-rotor densities vanish in their first grain.
   - With the deck default ΔE = 0.2·kT, this can be ~10% in thermal averages.
   - The usual remedy is to count states on fine cells (~1 cm⁻¹) and average into ME grains.
   - **Literature needed** (grain averaging of ρ and k(E); e.g. Gilbert & Smith 1990, or Holbrook–Pilling–Robertson 1996). The `books/gilbert_smith.txt` extract is empty: the PDF text layer may be missing.
4. **Final steady state for deep wells without a sink.** With k_th/ω ≲ 10⁻¹⁴, J is numerically singular in double precision. The driver refuses such rows (residual check) with an explanation, and the intermediate steady state should be used there. Example: the consecutive-activation test, second well 9000 cm⁻¹ deep, 298 K. This is physical: the final steady state is reached only after ≈ 0.1/k_uni (O02).
5. **Excited electronic levels** in `ElectronicLevels` are ignored; only the ground-level degeneracy is used.
6. **Barrierless (PST/SACM) module for chemical activation.** Its connection is deferred (your instruction). `TransitionStateModel::PhaseSpaceTheoryRRHO` and its W(E) are kept in `microcanonical_builder.rs` for that.
7. **Output noise.** `inertia::get_brot` prints "Iterative diagonalization is done in 1 steps." to stdout for every geometry. It is noise in the example output and was not changed.
8. **Pre-existing warnings** (not from the new code): an unused `ReactantType` import in `pst_troe_ushakov2006.rs`; unused `AU_TO_KCAL` and `CM1_TO_KCAL` imports in `thermal/thermofuncs.rs`.

---

## 8. Removed files (legacy master-equation engine; all recoverable from git HEAD)

**`src/masterequation/`**
- `api.rs`
- `chemical_activation_source.rs` (old)
- `chemical_network_builder.rs`
- `collisional_energy_transfer.rs`
- `energy_grained_me.rs`, `energy_grained_steady_state.rs`
- `graph_utils.rs`
- `input_deck.rs`
- `matrix_physics_assembly.rs`
- `network_builder.rs`
- `reaction_network.rs`
- `report.rs`
- `singlewell_solver.rs`
- `state_index.rs`
- `steady_state_chemical_activation_me.rs`
- `text_input.rs`

The earlier phase removed `current_session.txt`, `mastereq_ssumes.txt`, `example.rs` and `placeholder_microcanonical.rs`.

Also removed: inside `mess_input.rs`, `build_multiwell_network`, `MessBuildOptions` and `MessBuiltNetwork`, together with their 2 tests (now covered by `collision_parameters_and_sink_come_from_the_deck`). Inside `microcanonical_builder.rs`, `build_microcanonical_network_data`, `MicrocanonicalNetworkData`, `ArrayMicrocanonicalProvider` and `ChannelMicroModel`.

**`examples/`**
- `mess_build_network_and_rates.rs`
- `mess_case1_zzallyl_o2_run.rs`
- `steadystate_me_c2h3.rs` and `steadystate_me_c2h3_input.txt`
- `zzallyloh_o2_mess_mirror.rs` and `zzallyloh_o2_multiwell_input.txt`
- `cases/` (`c2h3_molecules.rs`, `c2h3_rrkm.rs`, `mod.rs`): the C₂H₃ data moved to `examples/c2h3_chemical_activation.inp`

The earlier phase removed `masterequation_input.rs` and `.txt`, `run_masterequation.rs`, and `steadystate_me.rs` with its input.

**Kept:** `high_pressure_limit.rs` (independent TST/Eyring k∞), `mess_input.rs`, `microcanonical_builder.rs`, `collisional_relaxation.rs`, the barrierless examples (`co_oh*.rs`, `phasespace_co_oh_capture.rs`, `zzallyloh_o2.rs`), the Case1 deck and the excerpt deck.

---

## 9. Test inventory of this phase (all passing)

**`chemical_activation_network`**
- aligned_grain_maps_equal_absolute_energies
- lowest_threshold_is_the_first_grain_with_an_open_channel
- mean_down_follows_the_power_law_from_its_reference_temperature
- validation_rejects_self_isomerization_and_rate_arrays_of_wrong_length

**`chemical_activation_operator`**
- column_sums_of_j_equal_the_losses_out_of_the_network (4 option combinations)
- j_is_detailed_balanced_with_boltzmann_weights_on_the_absolute_energy_scale
- intermediate_steady_state_absorbs_grains_ten_kt_below_the_lowest_threshold
- final_steady_state_keeps_every_grain_and_has_no_stabilization
- states_are_ordered_by_absolute_energy
- isomerization_detailed_balance_check_flags_inconsistent_reverse_rates
- isomerization_beyond_the_target_grid_is_an_error

**`chemical_activation_steady_state`**
- steady_state_satisfies_j_n_equals_f_for_both_solvers_and_all_options
- banded_cholesky_and_bicgstab_give_the_same_populations
- cholesky_refuses_an_operator_without_detailed_balance
- source_in_absorbed_grains_counts_as_directly_stabilized
- project_source_rejects_negative_or_empty_distributions

**`chemical_activation_sources`**
- source_from_the_reverse_rate_equals_the_transition_state_sum_of_states_form
- entrance_channels_into_several_wells_are_weighted_by_their_thermal_fluxes
- convolution_of_two_single_grains_is_a_single_grain_above_the_threshold
- convolution_with_a_sharp_partner_reduces_to_the_shift_approximation
- thermal_distribution_is_normalized_boltzmann
- sources_that_do_not_fit_on_the_receiving_grid_are_rejected

**`chemical_activation_observables`**
- yields_of_products_stabilization_and_sink_add_up_to_one
- zero_pressure_yields_are_the_rrkm_branching_of_the_nascent_distribution (10⁻⁷ Torr, within 10⁻⁴)
- high_pressure_intermediate_steady_state_stabilizes_everything (10⁹ Torr, Φ_stab > 0.9999)
- final_steady_state_without_a_sink_ends_entirely_in_products
- final_steady_state_through_the_only_channel_is_the_equilibrium_distribution
- rate_coefficients_average_k_over_the_normalized_well_distribution

**`chemical_activation_driver`**
- every_temperature_and_pressure_is_solved_in_order
- stabilization_grows_with_pressure
- thermal_entrance_source_is_rebuilt_at_every_temperature
- results_table_has_a_header_and_one_row_per_condition
- a_source_with_the_wrong_number_of_wells_is_an_error

**`consecutive_activation`**
- shift_treatment_moves_the_steady_state_distribution_by_minus_re_plus_partner_energy
- convolution_treatment_uses_the_thermal_partner_distribution
- the_second_master_equation_is_run_at_every_condition_of_the_first

**`chemical_activation_from_mess_input`**
- wells_share_one_grid_up_to_the_model_energy_limit
- isomerization_rates_obey_detailed_balance_exactly
- channels_follow_the_barriers_of_the_deck
- association_ilt_forms_the_entrance_channel_at_the_asymptote
- collision_parameters_and_sink_come_from_the_deck
- a_phase_space_barrier_without_an_ilt_block_is_refused
- tunneling_must_be_ignored_explicitly
- the_deck_runs_through_the_chemical_activation_driver

**`ilt_barrierless`**: 6 tests (§4.4).

**`banded_solvers`**
- banded_cholesky_matches_dense_cholesky
- banded_cholesky_rejects_an_indefinite_matrix

**`mess_input`**, new tests
- temperature_and_pressure_lists_are_read_in_torr
- exponent_cutoff_and_well_escape_rate_are_read
- wells_keep_the_order_of_the_input_deck
- inverse_laplace_transform_block_is_read_inside_a_barrier
- tunneling_blocks_are_detected
- inverse_laplace_transform_units_must_match_the_direction

**`collisional_relaxation`**: 2 tests, ported to the new function.

**`collision_kernels`**: 4 tests (earlier phase).

**Test-expectation corrections made during this phase**, with reasons; none was a code change:
- The consecutive-activation fixture grid was too short. The code correctly refused to truncate a 5·10⁻⁵ tail.
- The second consecutive ME was run in the final steady state for a deep well without a sink, which is numerically singular (item 7.4). The test now uses the intermediate steady state.
- The adapter's channel-opening grains: with classical rotors W‡(0) = 0, so the isomerization opens one grain above the TS grain. The convolved fragment density vanishes in two grains (item 7.3).
