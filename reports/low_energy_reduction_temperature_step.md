# Temperature step from the low-energy reduction of the exponential-down kernel

**Date:** 2026-10-06. **Status:** cause confirmed. Nothing is changed in the code; the choice of the low-energy rule is open and needs Peter's decision (Section 6).

## 1. Finding

In the ZZ-allyl + O₂ Case 2 validation, several MarXus results jump between 300 and 310 K. MESS is smooth there.

The jump appears in all four methods, because they share the operator 𝐉. It also appears with both tunneling models, and the high-pressure rate coefficients k∞ are smooth. It is therefore in the collisional part of 𝐉.

**Step in MarXus between 300 and 310 K, 760 Torr** (CSE, MESS Eckart model; `validation/ZZAllyl+O2_Gamma_Case2/cse_comparison.csv`):

| quantity | 290 → 300 K | 300 → 310 K | 310 → 320 K | MESS 300 → 310 K |
|---|---|---|---|---|
| k(R → G4) | ×0.891 | **×0.910** | ×0.874 | ×0.881 |
| k(R → G2) | ×0.844 | **×0.830** | ×0.832 | ×0.836 |

As a result, CSE's k(R → G4) deviates from MESS by +1.4 … +2.2% at 270–300 K, but by +4.9 … +6.1% at 310–330 K. The absorbing barrier's k(R → G4) shows the same ratio (×0.908).

## 2. Cause

`src/masterequation/collision_kernels.rs`, `exponential_down_kernel` (lines 74–120): the low-energy reduction.

**The rule.** Below the grain `low_cut`, every transition probability involving a grain i < low_cut is divided by

```
redfac_i = max(1, (ρ_{i+n_ref} / ρ_i)^m),   n_ref = floor(1.5 <ΔE_down> / ΔE) + 1   (line 84)
```

`low_cut` is the highest grain where ρ_{i+n_ref}/ρ_i > 3.

The factor is attached to the lower grain of each pair, so detailed balance stays exact. The exponent m starts at 1.0 and grows in steps of 0.1, up to 3.05, until the back substitution of Robertson (2019) eq. 4.16 succeeds.

**Origin.** The rule is variant V2 of `reports/chemical_activation_three_approaches.md`, after SSUMES `umolProb::setETProb`. In `reports/collision_kernel_detailed_balance_and_normalization.md` (Section 8) its choice is recorded as an open decision.

**Why it gives a step.** n_ref is an integer, but ⟨ΔE_down⟩(T) = 200 cm⁻¹ (T/300 K)^0.85 is continuous. In Case 2 the grains are 38 cm⁻¹, and 1.5⟨ΔE_down⟩/ΔE crosses 8 at ⟨ΔE_down⟩ = 202.67 cm⁻¹, that is at **T = 304.7 K**:

| T (K) | ⟨ΔE_down⟩ (cm⁻¹) | 1.5⟨ΔE_down⟩/ΔE | n_ref |
|---|---|---|---|
| 270 | 182.6 | 7.21 | 8 |
| 300 | 200.0 | 7.89 | 8 |
| 304 | 202.27 | 7.984 | 8 |
| 305 | 202.83 | 8.006 | **9** |
| 310 | 205.7 | 8.12 | 9 |
| 330 | 216.8 | 8.56 | 9 |

When n_ref changes by one grain, both `low_cut` and every factor redfac_i change. So the whole low-energy kernel of every well changes at once.

The exponent m can also jump by 0.1 at some temperature. This was not observed here.

**The rule is always active, not only as a repair.** The reduction applies whenever the density of states rises by more than a factor of 3 over n_ref grains, which is the case at the bottom of every well. It is not limited to conditions where the plain normalization fails. Its parameters (`low_cut`, n_ref, m) appear in no output.

## 3. Confirmation: 1 K scan

**Prediction:** if n_ref is the cause, the step must lie between 304 and 305 K and nowhere else.

**Setup.** The deck `marxus_input/case2_tstlevel_E.inp` with `TemperatureList 270 300 301 … 310 330` and `PressureList 760`, everything else unchanged. Keeping 270 K and 330 K in the list keeps the grid identical: grain 38 cm⁻¹, same top. Run with `--method cse --tunneling mess-eckart --threads 4`; the deck and output are in the session scratchpad.

| T (K) | lowest relaxation eigenvalue (1/s) | k(R → G2) | k(R → G3) | k(R → G4) (cm³/s) | ratio of k(R → G4) to the previous T |
|---|---|---|---|---|---|
| 300 | 3.2577e8 | 6.4000e-12 | 5.2614e-14 | 1.9133e-12 | |
| 301 | 3.2467e8 | 6.2901e-12 | 5.2074e-14 | 1.8903e-12 | 0.9879 |
| 302 | 3.2359e8 | 6.1817e-12 | 5.1536e-14 | 1.8673e-12 | 0.9878 |
| 303 | 3.2251e8 | 6.0747e-12 | 5.1000e-14 | 1.8444e-12 | 0.9878 |
| 304 | 3.2145e8 | 5.9692e-12 | 5.0467e-14 | 1.8217e-12 | 0.9877 |
| **305** | **3.0388e8** | **5.8101e-12** | **5.1182e-14** | **1.8560e-12** | **1.0188** |
| 306 | 3.0270e8 | 5.7079e-12 | 5.0638e-14 | 1.8328e-12 | 0.9875 |
| 307 | 3.0153e8 | 5.6071e-12 | 5.0096e-14 | 1.8097e-12 | 0.9874 |
| 308 | 3.0037e8 | 5.5077e-12 | 4.9556e-14 | 1.7868e-12 | 0.9873 |
| 309 | 2.9922e8 | 5.4097e-12 | 4.9020e-14 | 1.7639e-12 | 0.9872 |
| 310 | 2.9808e8 | 5.3131e-12 | 4.8485e-14 | 1.7412e-12 | 0.9871 |

The step lies exactly between 304 and 305 K; every other 1 K interval is smooth. Size of the step, after removing the smooth trend per kelvin:

| quantity | step at 304 → 305 K |
|---|---|
| lowest relaxation eigenvalue | −5.1% |
| k(R → G4) | +3.2% |
| k(R → G3) | +2.5% |
| k(R → G2) | −0.9% |
| k(R → escape) (a small negative entry, rounding level) | +15% |
| thermal loss of G2, G2 → R | −0.25%, −1% |
| thermal loss of G6 | −1% |
| k(R → IEPOX + OH) | about +0.2% |

## 4. Where it matters

- **Every method and every network** with the exponential-down kernel. Steps occur wherever the T grid crosses an integer value of 1.5⟨ΔE_down⟩/ΔE.
  - In C₂H₃ (21 cm⁻¹ grains, ⟨ΔE_down⟩ from 200 to 1003 cm⁻¹ between 300 and 2000 K), n_ref runs from 15 to 72. On a 250 K grid the steps are not visible as such.
- **The size of the low-energy effect.** It is not small in Case 2: 5% in the lowest relaxation eigenvalue, 3% in R → G4. So the low-energy kernel is not physically irrelevant for the bimolecular-to-well rate coefficients of these wells. This also bears on the few-% residual between MarXus and MESS that was attributed to "the collisional part" (Case 2 README, Section 4.4).
- **Thermal rate coefficients.** They change by up to 1% at the step (G2, G6).

## 5. What is not yet known

- Whether the plain normalization of Robertson eq. 4.16, without any reduction, succeeds for the Case 2 and C₂H₃ wells. If it does, the reduction is not needed there at all.
- How large the effect of the reduction itself is: with it, against without it, where both exist.
- Which of the two sides of the step (n_ref = 8 or 9) is closer to MESS cannot decide the question. R → G4 is closer with n_ref = 8, R → G2 slightly closer with n_ref = 9.

## 6. Decision needed (Peter)

The low-energy rule is a physics choice. The options found in the literature and code already reviewed:

1. **Keep the current rule** (SSUMES V2). Document the step, and choose temperature grids accordingly.
2. **Apply a reduction only where the plain back substitution fails.** Start without any reduction, i.e. plain Robertson eq. 4.16. Where no failure occurs there is then no step. Where the reduction does switch on, the switch is still discontinuous.
3. **The UNIMOL rule used in the older multiwell kernel path** (`matrix_physics_assembly.rs`; `reports/collision_kernel_detailed_balance_and_normalization.md`, Section 5.2). Grains below E₀/2 share the normalization coefficient of the grain above them. Its switch point is fixed by E₀, not by T, so it gives no T step.
4. **A different transition model at low energies** (Robertson 2019, p. 294, eqs. 4.60–4.61).

A continuous variant of rule 1, e.g. an interpolated n_ref, would be a MarXus invention without a paper. It is therefore not proposed.

**Diagnostic that can be added independently of the decision.** Report `low_cut`, n_ref and m per well and condition in the output (RUN SETTINGS or the network block), so that a switch is visible.

## 7. References

- S. H. Robertson, Comprehensive Chemical Kinetics 43 (2019): eq. 4.16 (normalization), p. 294 (breakdown for sparse states, eqs. 4.60–4.61).
- `reports/chemical_activation_three_approaches.md`: variant V2 (low-energy reduction factors) and the row "Low-energy treatment" of the comparison table.
- `reports/collision_kernel_detailed_balance_and_normalization.md`: Section 5.2 (UNIMOL E₀/2 rule) and Section 8 (the open decision).
- `validation/ZZAllyl+O2_Gamma_Case2/cse_comparison.csv`, `four_methods_comparison.csv`.
