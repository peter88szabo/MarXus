# How Olzmann treats stabilization, the absorbing barrier and the whole chemical-activation problem, and what MarXus implements

**MarXus, 2026-10-05.** Answer to Peter's question: "what does Olzmann do with this absorption barrier and treatment of the whole problem? did we implement that?"

The quotations below come from the text extracts of the papers in `papers/ChemAct/`. Abbreviations are those of `chemical_activation_implementation.md`.

## 1. What Olzmann does

### 1.1 The physical sink replaces the absorbing barrier

**O91** (Olzmann, Gebhardt, Scherzer 1991, 2-butyl + H, eqs. 1–6).
- **The sink is in the master equation.** The master equation of the chemically activated radical contains the physical sink, i.e. reactive collisions with the excess H atoms: "Every collision, reactive or nonreactive, removes a radical from level i (ω, second term) but only nonreactive collisions contribute to transitions j → i (ω_nr, third term)". This gives J = ωI − ω_nr P + K (eq. 2).
- **Radicals that survive react, they are not stabilized:** "all radicals that do not decompose, react with hydrogen atoms and hence do not contribute to stabilization S", so R = D₁ + R′ (eq. 6).
- **The product of that reaction is a second chemically activated species.** Its source is the convolution of the radical's steady-state distribution with the H-atom distribution ("the positive part of a one-dimensional Maxwellian for hydrogen atoms"), shifted by D₀ (eq. 8).
- **Second master equation:** J′ = ωI − ωP′ + K′ (eq. 9) with R′ = D₂ + S (eq. 11).
- **Stepladder with upward transitions**, which keeps detailed balance (eqs. 13–18). O91 discusses why ⟨ΔE⟩ obtained with and without upward collisions differ.

**O02** (Olzmann 2002, s-C₄H₉ + H, the central paper on this question; text after eq. 13).

*Two time regimes (Schranz & Nordholm 1984, O02 ref. 12):*

- **Intermediate steady state.** Decomposition competes with stabilization, D/S is time independent, and "the traditional chemical activation formalism applies". "In the respective steady-state master equation the quasi-irreversible stabilization has to be accounted for by either neglecting upward collisions [ref. 7] or, if detailed balancing is included, by introducing a lower absorbing barrier in the collisional deactivation cascade at energies below the lowest reaction threshold [ref. 9 = Holbrook, Pilling, Robertson, *Unimolecular Reactions* (1996)]."
- **Final steady state**, reached after ≈ 0.1/k_uni: "the stabilization reservoir is filled up … the ratio of decomposition over net stabilization … loses its physical sense, and one trivially has R1 = R2 and hence Φ2 = 1. It should be emphasized that a steady-state master equation with upward collisions included corresponds to this case if no absorbing barrier exists."

*The key statement.* With the consecutive bimolecular reactions, "these reactions form a physically reasonable sink for the intermediates, and the artificial condition of a lower absorbing barrier in the master equation must be dropped". The master equation then becomes J = ω(I − P) + K₂ + k₄[H]·I (eq. 6) and is solved in the **final steady state, without a barrier**.

*Validity window, shown numerically:*
- For 10⁻⁷ < k₄[H] < 10⁵ s⁻¹ the yield is Φ₂ ≈ 0.43, which "agrees with the result from a master-equation analysis neglecting reactions (4) and (4′) and using an absorbing barrier instead". That is the plateau in Fig. 2.
- The sink and the barrier pictures agree within 10% for **0.01 ω > k₄[H] > 10 k₂uni**.
- Outside that window, either excited radicals react (k₄[H] ≈ ω), or the thermal decomposition competes and Φ₂ → 1.
- Compared with the time window of the intermediate steady state, (0.1 λ_F)⁻¹ < t < (10 k_uni)⁻¹ (Schranz & Nordholm): "(k₄[H])⁻¹ can be loosely interpreted as an average reaction time for the chemically activated population".

*Thermal rate coefficient:* "k₂uni = 4.1 × 10⁻¹⁰ s⁻¹ as the **lowest eigenvalue of the matrix J** (routine tql1 … applied after symmetrization of J)".

**GO10** (González-García & Olzmann 2010, HSO₅, p. 12295).
- "In technical terms the intermediate steady state can be implemented by introducing a lower absorbing barrier into the master equation. **By omitting this absorbing barrier, that is by carefully observing the completeness of transition probabilities, the final steady-state solution** for ñ_s(E) is obtained."
- The calculation is done in the final steady state with the physical sink k₄[H₂O]:
  - k^ca = Σ (K_r Ñ_s)_i (eq. 9);
  - yield f₄ = k₄[H₂O]/(k^ca_2a + k^ca_2b + k₄[H₂O]) (eq. 10);
  - **k^th = λ₁, the lowest eigenvalue of J for [H₂O] = 0** (eq. 12), or the average of k(E) over the corresponding eigenvector.

**PO14** (Pfeifle & Olzmann 2014, isoprene + OH + O₂).
- **Every well is solved in the final steady state with its physical sink**, the bimolecular capture by O₂: "It is these bimolecular capture reactions that prevent the molecular population from becoming thermalized even in the long-time limit".
- **Coupled master equations** (consecutive activation; eqs. 8, 10–13).
- **Time-dependent solutions** N(t) = R Σ e_i (1 − e^{−λᵢt})/λᵢ · E_i (eqs. 3–4) show how the steady state is approached. The eigenvalues and eigenvectors come from "EISPACK routine tql2 after symmetrization". The steady state itself is "most conveniently obtained by directly solving R_a F = J N_ss" with a band solver.
- **No absorbing barrier is used.** For a well *without* a capture reaction (Fig. 3), the time-dependent population is shown to approach the thermal distribution.

### 1.2 Olzmann's approach in summary

1. **Physical sinks.** Write the master equation with the physical sinks (pseudo-first-order bimolecular reactions, reactive collisions), and use the **final steady state without an absorbing barrier**. The absorbing barrier is regarded as an *artificial* device of the classical D/S treatment, needed only when no physical sink exists.
2. **Observables.** Rate coefficients k^ca as averages of k(E) over the normalized steady-state distribution; yields; competition with the sink.
3. **Thermal rate coefficient** k^th = λ₁, the lowest eigenvalue of the symmetrized J.
4. **Time dependence** by eigenvalue expansion, when needed (PO14 eq. 3).
5. **Consecutive steps** by coupling master equations through convolution or shift (O91 eq. 8; PO14 eqs. 8, 10–13).
6. **Validity of the absorbing-barrier (D/S) picture.** It holds only in the window 0.01 ω > k_c[D] > 10 k_uni (O02).

## 2. What MarXus implements

| Olzmann element | MarXus | Where |
|---|---|---|
| J = ω(I − P) + K + k_c[D]·I with a physical sink, final steady state, no barrier | **yes** (`SteadyState::Final`; sink `Well::bimolecular_sink_s_inv`, from the deck `Escape`) | `chemical_activation_operator.rs` |
| Reactive collisions as part of the collisions (O91: ω vs ω_nr) | equivalent: our sink k_c[D] is added to ω(I − P), which corresponds to O91's ω = ω_nr + ω_r | — |
| Intermediate steady state with a lower absorbing barrier (O02, ref. HPR 1996) | **yes** (`SteadyState::Intermediate`); default 10 kT below the threshold (PR03, CD07), distance user-selectable | `chemical_activation_operator.rs` |
| k^ca = Σ k Ñ (GO10 eq. 9), yields Φ (O02 eq. 10), sink yield | **yes** | `chemical_activation_observables.rs` |
| Stepladder with detailed balance (O91 eqs. 13–18) | **yes** (with exponential down as a second model) | `collision_kernels.rs` |
| Coupled master equations: convolution and shift (PO14 eqs. 8, 10–13) | **yes** | `consecutive_activation.rs` |
| Partner H atom as the "positive part of a 1-D Maxwellian" (O91 eq. 8) | **no** special form; any partner distribution can be given to the convolution | — |
| **k^th = λ₁, lowest eigenvalue of J** (GO10 eq. 12; O02), and k^th by the eigenvector average (GO10 after eq. 12) | **yes** (update below): `thermal_rate_coefficients` reports k_uni = eigenvector average, with λ₁ beside it; solvers inverse iteration (default), Householder/QL, LAPACK DSYEVD; sum rule λ₁ = k_uni checked at runtime (warning above 1.5%) | `chemical_activation_eigen.rs`, `numeric/symmetric_eigen.rs`, `numeric/lapack_interface.rs` |
| **Time-dependent solution** by eigenvalue expansion (PO14 eqs. 3–4) | **yes** (library: `EigenSystem::time_dependent_population`; not yet in the example output) | `chemical_activation_eigen.rs` |
| **Validity window 0.01 ω > k_c[D] > 10 k_uni** and (0.1 λ_F)⁻¹ < t < (10 k_uni)⁻¹ (O02; SN84) | **yes** (library: `steady_state_window`, `EigenSystem::lambda_f`; not yet printed per condition) | `chemical_activation_eigen.rs` |

## 3. Consequence for the shallow-well / C₂H₃ problem

**The case.** The C₂H₃ benchmark (`validation/c2h3_mess_example/`) is a pure fall-off reaction without a physical sink. In Olzmann's framework there is no absorbing barrier to place for it. Its pressure-dependent thermal rate coefficient is

**k_uni(T, p) = λ₁(J)**, the lowest eigenvalue of J of the final steady state without a sink (GO10 eq. 12; O02),

and the association follows by detailed balance: k_assoc(T, p) = λ₁·K_eq, with K_eq = k∞_assoc/k∞_diss.

**Why this fits the high-temperature points.** It needs no barrier distance at all. It would therefore treat the 1500–2000 K conditions, where the absorbing barrier lies inside the thermal distribution, on the same footing as low T, provided λ₁ is still separated from the relaxation eigenvalues.

**A proposal that follows Olzmann's own numerics** (tql1/tql2 after symmetrization; implemented since, see §4):
- compute λ₁ of the symmetrized, banded J by inverse iteration with the banded Cholesky factor MarXus already has;
- inverse iteration converges fastest exactly where J is nearly singular, which is the case that makes the final steady state fail in double precision (deep well, low T);
- this would be the first part of the "eigenvalue route", and also the basis for O02's validity check of the absorbing-barrier picture.

Not implemented: Peter decided that the eigenvalue route comes later.

## 4. Update (2026-10-05, later the same day): the eigenvalue machinery is implemented

Additional abbreviations:
- SN84 = Schranz, Nordholm, Chem. Phys. 87, 163 (1984);
- NR92 = Press, Teukolsky, Vetterling, Flannery, *Numerical Recipes in Fortran*, 2nd ed. (1992).

**Why.** Peter: "we want Olzmann formulation and solution method … keep the current model, but we need to build Olzmann machinery … so we can avoid this stupid absorbing barrier"; "implement both approaches and the user can chose, and make the better one and more stable one as default".

**Implemented** (all TDD, all tests green):

- **`numeric/symmetric_eigen.rs`.**
  - Householder tridiagonalization (tred2) + QL with implicit shifts (tql2), Olzmann's EISPACK route (PO14 ref. 39; NR92 §11.2–11.3).
  - Inverse iteration with the banded Cholesky factor for the lowest eigenpairs, with deflation for λ₂ (NR92 §11.7).
- **`numeric/lapack_interface.rs`.** LAPACK DSYEVD for the full spectrum: Behemoth's binding, system OpenBLAS, feature `openblas` (default). See `eigen_solvers_lapack_and_davidson_assessment.md`.
- **`chemical_activation_eigen.rs`.**
  - `thermal_rate_coefficients`: k_uni = Σ_j k_j^th + k_c[D] (GO10's "averaging procedure", the reported rate coefficient), λ₁, λ₂, channel rates k_j^th, k∞ per channel, sink rates, populations and distributions.
  - It includes the **sum-rule check** λ₁ = k_uni (GO10 eq. 12). A deviation above 1.5% gives a **warning** with the explanation, the citation and what to do. Nothing is rejected (Peter's decision, 20:01).
  - `EigenSystem`: all eigenpairs, N(t) (PO14 eqs. 3–4), λ_F.
  - `steady_state_window` (O02).
- **Driver and example.**
  - `run_thermal_rate_coefficients`, `write_thermal_table`.
  - Example options `--steady-state eigenvalue|all`, `--eigen-solver inverse|full|lapack`, `--sum-rule-tolerance`.
  - For one well with one entrance, the association by detailed balance, k = λ₁·k∞,assoc/k∞,diss.

**Result for C₂H₃** (`validation/c2h3_mess_example_olzmann_eigen/README.md`).
- **No absorbing barrier and no barrier distance are needed.** The association agrees with MESS within ±2.5% from 750 to 1750 K. The barrier route gave −13% at 1500 K and −42% at 1750 K with 10 kT.
- **At 500–1000 K both routes agree** within 0.00–0.34%.
- **At 300 K λ₁ is below the double-precision resolution** (off by 10⁸; warnings). k_uni from the thermal eigenvector is still correct: +4.7 … +5.7% vs MESS. The association equals the absorbing-barrier route to 0.01% at 300–1000 K.
- **With the shifted factorization S + σI** (Peter's go-ahead, 20:18), the inverse iteration solves all 300 K conditions. The association equals the barrier route to 0.001 percentage points.
- **Exponential down is the default collision model** (Peter, 20:25). The stepladder on grains finer than its step splits J into O02's sub-equations (O02 p. 3616), whose eigenvalue treatment is open.

**Checks.** All validity checks, sum rules included, are collected in `master_equation_validity_checks.md`.
