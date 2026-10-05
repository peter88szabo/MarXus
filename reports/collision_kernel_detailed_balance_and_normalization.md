# Collision kernel: detailed balance and normalization

**MarXus** — energy-grained master equation, multiwell collision operator  
**Date:** 2026-10-05  
**Status:** fixed in the working tree (not committed); all 50 unit tests pass

---

## 1. Summary

The collisional energy-transfer kernel of a master equation must satisfy two independent physical
constraints (Robertson 2019, eqs. 4.4 and 4.6; Pilling & Robertson 2003, eqs. 15 and 19):

1. **Normalization.** For every source energy, the transition probabilities over all final energies sum to 1.
2. **Detailed balance.** At thermal equilibrium, every pair of energies exchanges equal fluxes in both directions.

Before the fix, neither MarXus multiwell kernel satisfied both:

| Kernel (before) | Normalization | Detailed balance |
|---|---|---|
| `CollisionKernelImplementation::Mess` | exact | **violated** (up to 85 % in the test well) |
| `CollisionKernelImplementation::Spd` | **not imposed** (one global scale factor) | exact |

The `Mess` kernel now uses the exponential-down model with normalization coefficients obtained by back
substitution from the top grain (Robertson 2019, eq. 4.16). This satisfies **both** constraints exactly.
Below half the lowest reaction threshold, the low-energy rule of Section 5.2 applies. The `Spd` kernel is
unchanged; its status is discussed in Section 7.

> **Note on MESS.** The former MarXus `Mess` kernel was documented as a "MESS-style per-source
> normalization", and it violated detailed balance. MESS itself, however, does **not** violate detailed
> balance. Its source code scales both elements of each grain pair by one common factor, which keeps
> detailed balance exact. What MESS does not satisfy exactly is the **normalization**: in the test well
> its total transition probability per collision reaches 4.7 (Section 4). The former MarXus kernel was
> therefore not a faithful reproduction of MESS, and the detailed-balance error was MarXus's own.
> MarXus does not copy MESS; it only reads MESS-format input decks.

---

## 2. Physical requirements

Notation:

| Symbol | Meaning |
|---|---|
| s, t | source and target grains |
| P(t\|s) | probability per collision of the transition s → t |
| ω | collision frequency |
| R_ts = ω P(t\|s) | rate of s → t; the master-equation matrix is stored as R[target, source] |
| ρ_i | density of states of grain i |
| E_i | energy of grain i |
| f_i = ρ_i exp(−E_i/kT) | Boltzmann distribution (unnormalized) |

| Requirement | Equation | Reference |
|---|---|---|
| Normalization | Σ_t P(t\|s) = 1 for every s | Robertson 2019 eq. 4.4; Pilling & Robertson 2003 eq. 15 |
| Detailed balance | P(t\|s) f_s = P(s\|t) f_t | Robertson 2019 eq. 4.6; Pilling & Robertson 2003 eq. 19 |
| Exponential down (deactivating, t < s) | P(t\|s) = A_s exp(−(E_s − E_t)/α) | Robertson 2019 eq. 4.7 |
| Activating (t > s), from detailed balance | P(t\|s) = A_t (ρ_t/ρ_s) exp(−(E_t − E_s)(1/α + 1/kT)) | Robertson 2019 eq. 4.11 |

The normalization coefficient A belongs to the **source** of the deactivating step. An activating step s → t
is the detailed-balance partner of the deactivating step t → s, so it carries **A_t, the coefficient of
the target**.

Inserting both forms into the normalization condition gives, for every grain i (Robertson 2019, eq. 4.16):

```
A_i · Σ_{j ≤ i} exp(−(E_i − E_j)/α)  +  Σ_{j > i} A_j (ρ_j/ρ_i) exp(−(E_j − E_i)(1/α + 1/kT))  =  1
```

This system is upper triangular. It is solved by **back substitution from the highest grain**, where
only deactivating collisions exist (a reflecting upper boundary; Robertson 2019, p. 278).

Detailed balance is **not** a numerical convenience. It is a symmetry of the underlying dynamics
(Robertson 2019, p. 275). If the kernel violates it:

- the Boltzmann distribution is no longer a stationary solution;
- the similarity-transformed operator W⁻¹RW, with W = √f, is no longer symmetric, so symmetric
  solvers (Cholesky, LDLᵀ) are not applicable.

---

## 3. What was wrong in MarXus

### 3.1 The former `Mess` kernel: detailed balance violated

In `src/masterequation/matrix_physics_assembly.rs` the kernel was built as follows:

1. Deactivating weights: exp(−ΔE/α).
2. Activating weights: exp(−ΔE/α)·(ρ_t/ρ_s)·exp(−ΔE/kT), the detailed-balance ratio applied to
   **unnormalized** weights.
3. Each **source column** divided by its own sum N_s.

For a pair of grains s < t:

```
R_ts f_s = ω (g ρ_t/ρ_s e^{−ΔE/kT} / N_s) ρ_s e^{−E_s/kT} = ω g ρ_t e^{−E_t/kT} / N_s
R_st f_t = ω (g / N_t)                    ρ_t e^{−E_t/kT} = ω g ρ_t e^{−E_t/kT} / N_t
```

Detailed balance therefore held only where N_s = N_t. That condition fails near the boundaries of the
grid and wherever ρ varies rapidly.

| Measurement (test well, Section 6) | Violation |
|---|---|
| Lowest grain pair | 17 % (7.10·10⁵ vs 8.59·10⁵) |
| Maximum over all pairs | 85 % |

### 3.2 Consequences that hid the error

Three other defects in the same code path made the detailed-balance error hard to see:

- **Inverted similarity transform.** The code built W·L·W⁻¹ instead of W⁻¹·L·W. With W = √f,
  only W⁻¹LW is symmetric for a matrix stored as R[target, source].
- **No symmetry check before Cholesky.** The "Direct" solver ran Cholesky, which reads only the
  lower triangle, on a non-symmetric matrix. It returned a wrong solution without any error.
- **Collision frequency effectively zero.** The bath-gas number density was computed as
  p/(3.262·10¹⁶·T) instead of p/(k_B·T), about 10³⁶ times too small. The multiwell master equation
  therefore ran in the collisionless limit, and the kernel barely mattered.

All three are fixed (Section 5.4).

### 3.3 The `Spd` kernel: normalization not imposed

The `Spd` kernel keeps the unnormalized weights (A ≡ 1), which preserves detailed balance exactly.
It then scales **all** columns by one global factor, so that no source exceeds a total probability
of 1. Its normalization is therefore not the exponential-down model's. In the test well the
total transition probability per collision ranges from 0.03 to 1.00, and it deviates from the exact
model by up to 97 %. This kernel has not been changed (Section 7).

---

## 4. Comparison with MESS (for reference only)

MESS builds the kernel in `MasterEquation::Well::_set_kernel`, `src/libmess/mess.cc`, lines 6566–6614
(2026 source; this is the file compiled into the `mess` executable, and `new_mess.cc`, lines 392–478, contains an
identical but uncompiled copy). The grid index increases toward **lower** energy, and the
corresponding Boltzmann factor is `exp(+e·ΔE/kT)`.

1. For every pair (higher-energy grain h, lower grain l), MESS sets the deactivating element to
   exp(−ΔE/α) and the activating element by detailed balance.
2. **Both elements of the pair are divided by the same factor** `nfac(h)`. This factor is computed for
   the higher grain h and sums the elastic term, the deactivating steps from h, and the activating
   steps into h.
3. The diagonal is fixed from conservation.

A common factor per pair preserves the detailed-balance ratio, so **MESS satisfies detailed balance
exactly**. A source's deactivating probabilities, however, carry the factor of the source, while its
activating probabilities carry the factors of the target grains. The normalization of eq. 4.16 is
therefore not satisfied.

The comparison in Section 6 is a simplified single-exponential re-evaluation of these lines on the MarXus
test well. It was done for analysis only and is not part of MarXus.

---

## 5. The fix: what MarXus uses now

### 5.1 Exponential down with exact normalization (Robertson 2019, eq. 4.16)

`CollisionKernelImplementation::Mess` in `src/masterequation/matrix_physics_assembly.rs`:

1. ⟨ΔE_down⟩(T) = α(T) from the energy-transfer parameters of the well.
2. Normalization coefficients A_i by back substitution of eq. 4.16 from the highest grain of the well,
   restricted to the collision band.
3. Off-diagonal rates:
   - deactivating: R_ts = ω A_s exp(−(E_s − E_t)/α)
   - activating: R_ts = ω A_t (ρ_t/ρ_s) exp(−(E_t − E_s)(1/α + 1/kT))
4. Diagonal: R_ss = −Σ_{t≠s} R_ts. Population is conserved exactly, and the elastic probability is A_s.

Detailed balance holds exactly by construction: the activating element is the detailed-balance partner
of the deactivating one, with the same coefficient. Normalization holds exactly wherever eq. 4.16 is
used.

### 5.2 Low-energy rule (E < E₀/2)

For sparse low-energy states, the back substitution can return non-positive coefficients. Robertson (2019,
p. 294) shows that this is caused by the sparsity of the states and not by numerical error. MarXus applies
the rule of Gilbert's UNIMOL master-equation code (`mas55c3.f`, NREACT/NCUT with NCUT = 2; the file is shipped in
the `unimol/` directory of the SSUMES distribution):

- E₀ is the lowest reaction threshold of the well: the first grain where any of its channels,
  dissociation or isomerization, has k(E) > 0.
- Grains below **E₀/2** all share the normalization coefficient of the grain above them. Their
  populations are always at equilibrium, so the precise kernel there has no physical effect on the
  rate coefficients.
- Detailed balance remains exact below E₀/2, because the activating probabilities use the same
  coefficients.
- Normalization below E₀/2 is approximate. In the test well the total transition probability per
  collision there is 1.03–2.20.

**Errors instead of silent repairs.** UNIMOL silently copies a neighbouring coefficient if the back
substitution still fails above the cutoff. MarXus does not. A breakdown **above E₀/2**, or in a well
**without any reactive channel** (where E₀ is undefined), stops with an error that names the well and
the grain. In such a case a different transition model is needed (Robertson 2019, p. 294, eqs. 4.60–4.61).

### 5.3 Energy-transfer parameters from MESS-format decks

In the input format, `Factor[1/cm]` is ⟨ΔE_down⟩ at T₀ = 300 K and
⟨ΔE_down⟩(T) = Factor·(T/T₀)^Power (MESS manual). The reader previously stored `Factor` as MarXus's
α(1000 K), so α(300 K) was 72 cm⁻¹ instead of 200 cm⁻¹ for Factor = 200 and Power = 0.85.

- It is now converted exactly: α(1000 K) = Factor·(1000/300)^Power.
- A missing `Factor` or `Power` is now an error. It used to default silently to 200 cm⁻¹ and 0.85.

### 5.4 Related fixes in the same code path

| Defect | Fix | Reference |
|---|---|---|
| Similarity transform W L W⁻¹ (multiwell and single-well stepladder) | W⁻¹ L W | W = √f, R[target, source] |
| Cholesky on a non-symmetric matrix | Cholesky/LDLᵀ only if the relative asymmetry ≤ 10⁻¹⁰, otherwise BiCGSTAB | — |
| LDLᵀ: wrong 2×2-pivot solves, missing Bunch–Kaufman test | Fixed | Bunch & Kaufman 1977 |
| Collision frequency ~10⁻³⁶ too small | n = p/(k_B T), ⟨v⟩ = √(8k_B T/πμ) from physical constants | Troe 1977 eq. 3.1 |
| Lennard-Jones parameters: species only | σ_AM = (σ_A+σ_M)/2, ε_AM = √(ε_A ε_M) | Troe 1977, Sec. III |
| Ω(2,2)* used outside its validity range | Error outside 0.3 ≤ kT/ε ≤ 500 | Troe 1977 eq. 3.3 |

---

## 6. Verification

### 6.1 Test well

All numbers in this report refer to the following synthetic well, which has a steeply rising density of
states:

| Parameter | Value |
|---|---|
| Grains | 60 |
| ΔE | 20 cm⁻¹ |
| T | 500 K |
| ⟨ΔE_down⟩ | 166.4 cm⁻¹ |
| Collision band | ±20 grains |
| ρ_i | (1 + 0.05 i)¹⁰ |
| Reaction threshold | grain 40 (so E₀/2 is grain 20) |

### 6.2 Comparison of kernel constructions

| Construction | Max. detailed-balance violation | Total transition probability per collision Q(s) = Σ_{t≠s} P(t\|s) | Max. deviation of Q from eq. 4.16 (grains ≥ E₀/2) |
|---|---|---|---|
| Former MarXus `Mess` (per-column normalization) | **8.5·10⁻¹** | 0.88 – 1.00 | 2.0·10⁻² |
| MESS (pairwise normalization) | 7·10⁻¹⁶ | **0.33 – 4.71** | 6.3·10⁻¹ |
| MarXus `Spd` (no normalization, global cap) | 7·10⁻¹⁶ | 0.03 – 1.00 | 9.7·10⁻¹ |
| **New MarXus `Mess` (eq. 4.16 + E₀/2 rule)** | **1·10⁻¹⁵** | 0.88 – 0.99 (≥ E₀/2); 1.03 – 2.20 (< E₀/2) | **0** |

Q(s) > 1 means more than one transition per collision, which is unphysical. These numbers characterize
this synthetic well only. Real molecules have different densities of states, but the qualitative
behaviour of each construction is the same.

### 6.3 Unit tests

| Test | What it checks |
|---|---|
| `matrix_physics_assembly::tests::exponential_down_kernel_is_normalized_and_obeys_detailed_balance` | Every column conserves population, and R_ij f_j = R_ji f_i to 10⁻¹⁰ |
| `matrix_physics_assembly::tests::normalization_breakdown_without_reaction_threshold_is_an_error` | A breakdown without E₀ is reported, not repaired |
| `steady_state_chemical_activation_me::tests::steady_state_populations_solve_the_master_equation` | Returned populations satisfy L p + s = 0 (‖·‖/‖s‖ < 10⁻⁶) for both kernels and all three solvers |
| `energy_grained_steady_state::tests::stepladder_operator_is_the_symmetrized_master_equation` | The single-well stepladder operator equals W⁻¹JW |
| `mess_input::tests::energy_transfer_factor_is_the_value_at_300_kelvin` | α(300 K) = Factor |
| `mess_input::tests::lennard_jones_parameters_combine_species_and_bath_gas` | Combining rules |
| `collisional_relaxation::tests::collision_frequency_uses_ideal_gas_density_and_mean_relative_speed` | Z = πσ²⟨v⟩nΩ, checked against an independent SI evaluation |
| `collisional_relaxation::tests::collision_integral_outside_its_validity_range_is_an_error` | Range check of Ω(2,2)* |

Each test was first observed to **fail** on the old code and then to **pass** after the fix.

---

## 7. Open points (decisions needed)

1. **`Spd` kernel.** It keeps detailed balance but does not normalize, and it deviates by up to 97 % from
   the exponential-down model. Options: keep it as a documented alternative, or remove it.
2. **Name of the variant `Mess`.** It now implements Robertson's eq. 4.16 rather than any MESS
   procedure. Renaming it (for example to `ExponentialDown`) would also change the input-deck keyword
   `mess`.
3. **Collision band.** The multiwell band is fixed at ±20 grains, and the deck keyword
   `ExponentCutoff` is ignored. A band of 20·ΔE ≈ 2.4 ⟨ΔE_down⟩ in the test well truncates the
   exponential.
4. **Wells without reactive channels.** These use strict eq. 4.16 everywhere. A breakdown is an error.
5. **Inter-well detailed balance.** Well offsets and transition-state thresholds are rounded to the grain
   separately, which breaks microscopic reversibility between wells. The code audit estimated up to
   13 % in the Case1 deck; this figure has not yet been re-verified. This is not addressed yet.
6. **Multiple exponentials** (`Fraction` in the input format) are not supported.

---

## 8. References

- S. H. Robertson, *Unimolecular Kinetics, Parts 2 and 3: Collisional Energy Transfer and the Master
  Equation*, Comprehensive Chemical Kinetics **43**, Elsevier (2019). Ch. 4: eqs. 4.4, 4.6, 4.7, 4.11,
  4.16 (pp. 275–278); breakdown of the normalization for sparse states, p. 294.
- M. J. Pilling, S. H. Robertson, *Master equation models for chemical reactions of importance in
  combustion*, Annu. Rev. Phys. Chem. **54**, 245 (2003). Eqs. 14–21.
- J. Troe, *Theory of thermal unimolecular reactions at low pressures. II*, J. Chem. Phys. **66**, 4758
  (1977). Sec. III, eqs. 3.1–3.3.
- J. R. Bunch, L. Kaufman, *Some stable methods for calculating inertia and solving symmetric linear
  systems*, Math. Comp. **31**, 163 (1977).
- UNIMOL master-equation code (R. G. Gilbert), `mas55c3.f` as shipped in the SSUMES distribution
  (normalization vector ANORM; low-energy rule NREACT/NCUT).
- SSUMES (A. Miyoshi), `source/ssulibc.cc`, `umolProb::setETProb` (lines 1015–1064): its own
  chemical-activation solver uses a *different* low-energy remedy, reduction factors
  redfac_i = (ρ_{i+upref}/ρ_i)^m on the lower grain of each pair (detailed balance also exact). Which rule
  MarXus should finally use is an open decision (see `chemical_activation_three_approaches.md`).
- MESS source (2026), `src/libmess/mess.cc` lines 6566–6614 (compiled; identical copy in the uncompiled
  `new_mess.cc` lines 392–478), and the MESS manual (Exponential model:
  `Factor`, `Power`, `ExponentCutoff`). Used only to interpret the input format and for the comparison in
  Section 4.
