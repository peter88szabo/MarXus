# Multiwell chemical-activation master-equation code: design

**MarXus**: a steady-state chemical-activation solver for multiwell networks. It follows Olzmann's
formulation and uses SSUMES' `carate` program as its template for workflow and outputs.  
**Date:** 2026-10-05  
**Status:** implemented (test-driven, uncommitted); details, results and open items in `chemical_activation_implementation.md`.

Background and literature: `chemical_activation_three_approaches.md` and
`collision_kernel_detailed_balance_and_normalization.md`.

## Decisions (Peter, 2026-10-05)

1. **Collision models.** Both are available:
   - Olzmann's **stepladder** (Olzmann 1991, eqs. 13–18);
   - **exponential down**, as used in SSUMES and MESS.
2. **Low-energy treatment for exponential down.** SSUMES-style reduction factors on the lower grain of
   each pair. These keep normalization and detailed balance exact. The UNIMOL E₀/2 rule is replaced.
3. **Steady state.** The user chooses either:
   - **final** (no absorbing barrier; GO10 p. 12295, PO14 p. 238), or
   - **intermediate** (absorbing barrier, whose flux is stabilization; O02, SSUMES `truncate`).
4. **Barrier position.** By default **10 k_BT below the lowest reaction threshold** of each well
   (Pilling & Robertson 2003; Carstensen & Dean 2007, p. 125). An explicit grain can be given instead.
5. **Entrance channel.** Barrierless entrance channels use the inverse Laplace transform (ILT, Green & Pilling
   1986) until the barrierless module is connected. The ILT is checked before use.

## Workflow (template: SSUMES `carate`)

For each (T, p):

1. compute the collision frequency and ⟨ΔE_down⟩(T) of each well;
2. build the collision kernel per well (`collision_kernels.rs`);
3. assemble J;
4. build the source F;
5. solve J·N = F;
6. compute the observables;
7. write a table row.

## Equations

| Quantity | Equation | Source |
|---|---|---|
| Master equation | dN/dt = R·F − J·N, J = ω(I − P) + K + k_c[D]·I | PO14 eq. 2; O02 eqs. 5–6 |
| Steady state | J·N^s = R·F, N^s = R·J⁻¹F | PO14 eq. 5; O91 eqs. 2–3 |
| Collision frequency | ω = Z_LJ[M] = πσ²⟨v⟩n·Ω^(2,2)*, σ_AM = (σ_A+σ_M)/2, ε_AM = √(ε_Aε_M) | Troe 1977 eqs. 3.1–3.3 |
| Exponential down | Robertson CCK 43 eqs. 4.7, 4.11, 4.16 with low-energy reduction factors | `collision_kernels.rs` |
| Stepladder | P(i+n\|i) = A/(1+A), P(i\|i+n) = 1 − P(i+n\|i), A = (ρ_{i+n}/ρ_i)·e^(−ΔE_SL/kT) | O91 eqs. 15–18 |
| Thermal source | F(E) ∝ W‡(E−E₀)·e^(−(E−E₀)/kT), equivalently ρ(E)·k_entrance(E)·e^(−E/kT) | GO10 eq. 15; PO14 eq. 7 |
| Non-thermal source | F(E) = ∫ ñ_A(ε) ñ_B(E−E₀−ε) dε | PO14 eq. 8 |
| Consecutive activation | partner convolution, or the shift f₂(E) = ñ₁^ss(E + RE − ⟨E⟩_partner) | PO14 eqs. 10–13 |
| Yield of channel r | Φ_r = Σ_i k_r(E_i) N_i (normalized F) | O91 eq. 5; O02 eq. 10 |
| CA rate coefficient | k_r^ca = Σ_i k_r(E_i) Ñ_i, Ñ = N/Σ N | GO10 eqs. 8–9; PO14 eq. 6 |
| Stabilization (intermediate) | Φ_stab = Σ_j N_j · ω Σ_{t < barrier} P(t\|j), plus isomerization flux landing below the barrier | O02 p. 3616; SSUMES k_stab |
| Bimolecular sink yield | Φ_sink = k_c[D]·Σ_i N_i | O02 eq. 5; GO10 eq. 10 |
| Thermal rate coefficient | k^th = lowest eigenvalue of J (final steady state) | GO10 eq. 12 |
| Mass balance (check) | Σ_r Φ_r + Σ_w Φ_stab,w + Σ_w Φ_sink,w = 1 | — |

## Files

All in `src/masterequation/`. `mod.rs` receives only the module declarations.

| File | Content |
|---|---|
| `collision_kernels.rs` | exponential-down and stepladder kernels (pure functions, tested) — **done** |
| `chemical_activation_from_mess_input.rs` | input deck (MESS format) → network; ILT keyword for barrierless channels — **done** |
| `chemical_activation_network.rs` | wells, channels (to products or to another well), collision parameters, bimolecular sink, conditions, collision-model and steady-state options — **done** |
| `chemical_activation_operator.rs` | assembly of J on the common grain grid, stabilization (leakage) rates, isomerization coupling with a detailed-balance check — **done** |
| `chemical_activation_sources.rs` | Olzmann sources: thermal (entrance channel or W‡ array, e.g. from ILT), arbitrary distribution, single grain, convolution, shift approximation — **done** |
| `chemical_activation_steady_state.rs` | symmetrization W⁻¹JW, Cholesky or BiCGSTAB, residual check — **done** |
| `chemical_activation_observables.rs` | yields, k^ca, stabilization, sink, population fractions, normalized distributions per well — **done** |
| `chemical_activation_driver.rs` | carate-like loop over T and p, results table (CSV) — **done** |
| `consecutive_activation.rs` | chaining ME₁ → ñ₁^ss → f₂ → ME₂ — **done** |

**Grid rules**

- All wells share one grain width ΔE.
- Each well has an integer bottom offset on the common absolute grid, so inter-well alignment is
  exact and needs no rounding.
- Grain 0 of each well is its well bottom.

## Validation (tests)

- **Mass balance:** the yields sum to 1 in both steady-state modes.
- **Zero pressure:** Φ_r → Σ_i F_i·k_r(E_i)/k_tot(E_i), RRKM branching averaged over the nascent
  distribution.
- **High pressure, intermediate mode:** Φ_stab → 1.
- **Residual:** ‖J·N − F‖/‖F‖ small for every solver path.
- **Detailed balance of the collision operator:** exact for both kernels.
- **Consecutive activation:** shift and convolution against analytic cases.

## Code removed after the new code replaced it (done; list in `chemical_activation_implementation.md` §8)

- **Legacy multiwell engine:** `steady_state_chemical_activation_me.rs`, `matrix_physics_assembly.rs`,
  `api.rs`, `network_builder.rs`, `chemical_network_builder.rs`, `input_deck.rs`, `report.rs`.
- **Single-well solver:** `singlewell_solver.rs`, `energy_grained_me.rs`, `energy_grained_steady_state.rs`,
  `collisional_energy_transfer.rs`, `text_input.rs`, and the old `chemical_activation_source.rs` (its
  sources move to `chemical_activation_sources.rs`).
- **Kept as input front-ends:** `mess_input.rs` and `microcanonical_builder.rs`, adapted to produce the new
  network. Also kept: `collisional_relaxation.rs` (collision frequency) and `high_pressure_limit.rs`.
- **Already removed:** the AI session transcript, the AI-written SSUMES notes, the fake-data provider and
  the examples that used it.
