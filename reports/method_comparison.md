# Comparison of the four methods: identities, agreements, diagnostics

**Date:** 2026-10-06. **Request (Peter):** "find and make analysis for the comparison of the methods (as you wrote in the main README that they have to fulfill some identities and measured agreement (or not, just deviate), also try to find other information to extract and if possible plot".

**Script:** `validation/method_comparison.py`. One script for both validation systems. It reads the outputs of their run scripts, runs no master equation, and writes into each system's directory:

| output | content |
|---|---|
| `method_comparison.csv` | every compared quantity: system, check, kind (identity / approximate / MESS), T, p, quantity, the two values and their names, relative deviation, extra (e.g. the CSE separation), and whether λ₁ is below the double-precision floor at that condition |
| `plots/method_identities.png` | the exact identities, \|a/b − 1\| per condition and quantity, log scale |
| `plots/method_agreement.png` | left: the approximate relations against the CSE separation Λ_N/Λ_{N+1}; right: each method against MESS versus T (band: range over the pressures) |
| `plots/method_diagnostics.png` | CSE separation, sum rule of the thermal eigenpair, reservoir grains per well, run time per method |

**Run:** `source ~/.venvs/science/bin/activate && python3 validation/method_comparison.py`. It is also the last step of the validation re-runs.

**Data.** All numbers below are from the re-run of 2026-10-06 with the low-energy reservoir state (`reports/low_energy_reservoir_state.md`).

## 1. What is checked

**Exact identities** (must hold to the printed precision):

| | identity | source |
|---|---|---|
| I1 | pulse (TimeIntegration at the last time) = final steady state: Y_r(∞) = ∫k_rᵀe^{−𝐉t}F dt = k_rᵀ𝐉⁻¹F | README, "Exact identity 1" |
| I2 | long-time yields from the CSE rate coefficients = final steady state (% of the net reaction) | README, "Exact identity 2"; G13 eqs. 21, 25–30 |
| I3 | one well: CSE k(W → P) = k_uni of SteadyStateOlzmann (the eigenvector average) | GO10 after eq. 12 |
| I4 | late-time decay of the pulse, −d ln N/dt from the last two output times with 10⁻⁸ < N < 10⁻³, = k_uni | the thermal eigenpair; the time integration uses no eigenvector |
| I5 | conservation: populations + yields of the pulse = 100% at the last time | column sums of 𝐉 |
| I6 | CSE loss balance (G13 eq. 29); capture balance: R → wells + R → products = capture − return (G13 eq. 22) | G13 |
| I7 | sum rule λ₁ = k_uni of SteadyStateOlzmann (MarXus's own full-precision deviation) | GO10 eq. 12 |
| E | eigen-solvers: k_uni from inverse iteration = from LAPACK (C₂H₃) | — |

**Approximate relations** (hold where the chemical and relaxation time scales separate):

| | relation | why approximate |
|---|---|---|
| A1 | prompt R → P: SteadyStateAbsorbingBarrier k∞Φ_P against CSE (G13 eq. 21) | flux definition (barrier) vs eigenmode definition (CSE) |
| A2 | R → W: barrier k∞Φ_stab,W against CSE (G13 eq. 28), each well and the total; for one well also against SteadyStateOlzmann's k_uni·K | where the barrier sits in the well's thermal distribution; the separation |
| A3 | prompt yield + Σ_W stabilization × thermal fate of W = final steady state | assumes the stabilized molecules are thermalized before they react |
| A4 | detailed balance of the association/dissociation pair, [k(P → W)/k(W → P)]/K − 1, for CSE and for MESS (one well) | CSE's rate coefficients satisfy it only with separated time scales |

**Against MESS:** each method's quantity that MESS also gives (R → P5, R → G2, R → G3, R → G4 for ZZ-allyl + O₂; association and dissociation for C₂H₃).

**Filters, so that nothing meaningless is compared:**
- **Rounding-level CSE entries.** CSE entries below 10⁻⁶ of the largest of their row are rounding noise (and may be negative, as in MESS), e.g. R → G6 ≈ 10⁻²² cm³/s. They are not compared.
- **Precision floor.** Conditions where the precision floor ε·max(S_ii) exceeds 10⁻² k_uni, i.e. λ₁ is not resolved in double precision (C₂H₃ at 300 and 500 K), are flagged in the CSV, plotted as grey crosses, and left out of the maxima of the identities.
- **MESS duplicate column name.** MESS's species table names the escape column of G4 "G4" as well. The reader renames the second occurrence `escape(G4)`.

## 2. Results: ZZ-allyl + O₂, Case 2 (four wells; MESS Eckart model; 21 conditions)

**Exact identities** (`plots/method_identities.png`):

| identity | n | max \|a/b − 1\| |
|---|---|---|
| I1 pulse = final steady state (every exit) | 81 | equal in all printed digits |
| I5 pulse conservation | 21 | equal in all printed digits |
| I2 CSE long-time = final steady state | 80 | 5.3·10⁻⁴, only P7; P5, P1 and the escape equal in all printed digits |
| I4 pulse decay = k_uni | 21 | 2.1·10⁻⁶ |
| I6 CSE loss balance | 21 | 1.1·10⁻⁷ |
| I6 CSE capture balance | 21 | 5.7·10⁻⁷ |
| I7 sum rule | 21 | 5.3·10⁻⁹ |

**P7** is at most 0.008% of the reaction and is formed through G6. The CSE entries on that path (R → G6 ≈ 10⁻²² cm³/s) are at the rounding level, which limits its reconstruction.

**Approximate relations** (`plots/method_agreement.png`, left):

| relation | range |
|---|---|
| A1 prompt R → P, barrier vs CSE | P1 0.0 … 0.1%, P5 0.0 … 0.44%, P7 0.0 … 0.2% |
| A2 R → W, barrier vs CSE | G2 −5.7 … −0.6%, G3 +3.2 … +10.7%, G4 −18.3 … −9.0% |
| A2 total stabilization, barrier vs CSE | −9.2 … −2.8% |
| A3 prompt + stabilization × fate = final | max 1.8·10⁻³ |

**Reading.**
- **Prompt chemical activation is the same quantity in both methods,** to within 0.44%.
- **Stabilization is not.** The absorbing barrier counts the flux 10 k_BT below the lowest threshold of each well. Part of a well's thermal distribution can lie above that barrier; this is not quantified here, because the thermal distributions above the barrier were not evaluated. So the barrier counts less stabilization (−3 … −9% in total) and splits it differently among the wells. CSE counts the population of the chemical eigenmodes.

**Against MESS** (right panel, all 21 conditions):

| quantity | CSE | absorbing barrier | Olzmann / TimeIntegration (overall) |
|---|---|---|---|
| R → P5 | −3.6 … −3.0% | −3.5 … −2.6% | −2.6 … +3.8% (includes the thermal formation through the wells; MESS's R → P5 is prompt only) |
| R → G2 | +3.0 … +5.2% | −0.7 … +2.4% | — |
| R → G3 | −1.2 … −0.7% | | — |
| R → G4 | −1.0 … −0.3% | −19.0 … −9.5% | — |

**Diagnostics** (`plots/method_diagnostics.png`):
- **CSE separation** Λ_N/Λ_{N+1}: 0.049 … 0.10 at 760 Torr, up to 0.137 at 330 K and 500 Torr.
- **Sum rule:** 10⁻⁹ … 5·10⁻⁹.
- **Reservoir grains:** 4 … 10 per well.
- **Run time on 4 cores** (21 conditions): absorbing barrier 1–2 s, SteadyStateOlzmann 6 s, CSE 8 s, TimeIntegration 111 s.

## 3. Results: H + C₂H₂ ⇌ C₂H₃ (one well; exact Eckart; 40 conditions)

**Exact identities**, at the 30 conditions where λ₁ is resolved (750–2000 K):

| identity | n | max \|a/b − 1\| |
|---|---|---|
| I1 pulse = final steady state | 28 | equal in all printed digits |
| I5 pulse conservation | 30 | equal in all printed digits |
| I3 one well: CSE k(W1 → P1) = k_uni | 30 | equal in all printed digits |
| E inverse iteration = LAPACK (k_uni) | 30 | equal in all printed digits |
| I6 CSE capture balance | 30 | equal in all printed digits |
| I4 pulse decay = k_uni | 30 | 3.1·10⁻⁶ |
| I6 CSE loss balance | 30 | 3.0·10⁻⁵ |
| I7 sum rule | 30 | 3.2·10⁻⁸ |

**At 300 and 500 K, below the double-precision floor** (floor/k_uni = 40 … 10¹³; grey crosses):
- pulse conservation 2·10⁻⁴ … 6·10⁻⁴;
- CSE loss balance up to 4;
- sum rule up to 10⁹;
- CSE k(W → P) vs k_uni and inverse iteration vs LAPACK 5·10⁻⁵.

This is the known double-precision limit for the deep well without a sink (`reports/eigen_solvers_lapack_and_davidson_assessment.md`), not a method error.

**Approximate relations, against the separation** (`plots/method_agreement.png`, left). The deviations grow with Λ₁/Λ₂ and reach about 10⁻¹ at Λ₁/Λ₂ ≈ 0.1:

| T (K) | Λ₁/Λ₂ | CSE association vs Olzmann k_uni·K | CSE pair vs K (A4) | MESS pair vs K (A4) | barrier association vs Olzmann |
|---|---|---|---|---|---|
| 300–750 | ≤ 10⁻⁶ | 0.00% | 0.00% | 0.00 … −0.04% | 0.00 … +0.02% |
| 1000 | 3·10⁻⁵ … 1.4·10⁻⁴ | −0.01% | −0.01% | −0.02 … −0.15% | 0.00 … +0.02% |
| 1250 | 10⁻³ … 3.4·10⁻³ | −0.05 … −0.21% | −0.05 … −0.21% | −0.13 … −0.58% | −0.4 … −1.2% |
| 1500 | 8·10⁻³ … 0.019 | −0.5 … −1.7% | −0.5 … −1.7% | −0.7 … −2.1% | −4.2 … −9.1% |
| 1750 | 0.029 … 0.054 | −2.6 … −6.1% | −2.6 … −6.1% | −2.6 … −6.1% | −20 … −34% |
| 2000 | 0.065 … 0.10 | −7.2 … −14.0% | −7.2 … −14.0% | −7.2 … −13.9% | (not defined) |

**Reading:**
- **For one well, the CSE association and Olzmann's k_uni·K differ exactly by CSE's departure from detailed balance:** the A2 and A4 columns are identical.
- **MESS's own pair departs from K by the same amount as MarXus's CSE pair at 1750–2000 K.** The departure is a property of the CSE rate coefficients when the time scales do not separate, not an error of either code. This settles the question left open in the earlier eigen-solver validation, whether MESS's high-T pair was in error.
- **The absorbing barrier fails first.** Its distance from the well bottom shrinks with T: 9.6 k_BT at 1000 K, 3 k_BT at 1500 K, 1.2 k_BT at 1750 K, below the bottom at 2000 K.

**Against MESS** (right panel; range over the pressures):

| T (K) | dissociation (CSE = Olzmann) | association CSE | association Olzmann (k_uni·K) | association barrier |
|---|---|---|---|---|
| 300 | +4.7 … +5.7% | +4.9 … +5.8% | +4.9 … +5.8% | +4.9 … +5.8% |
| 750 | +1.3 … +1.8% | +1.4 … +1.9% | +1.4 … +1.9% | +1.5 … +1.9% |
| 1000 | −0.8 … +0.5% | −0.6 … +0.5% | −0.6 … +0.5% | −0.6 … +0.5% |
| 1500 | −3.5 … −1.9% | −3.1 … −1.7% | −1.4 … −1.2% | −10.4 … −5.3% |
| 1750 | −2.5 … −1.8% | −2.4 … −1.8% | +0.8 … +3.9% | −31 … −20% |
| 2000 | −2.6 … −2.1% | −2.6 … −2.1% | +5.5 … +13.2% | — |

**Reading:**
- **300–750 K:** the offsets are the known exact- vs semiclassical-Eckart difference of the high-pressure limits (`validation/c2h3_mess_example/README.md`, §5.2).
- **High T:** MarXus's CSE agrees with MESS (which is a CSE code) within −3.5 … −0.8% at 1250–2000 K (association −3.1 … −0.8%, dissociation −3.5 … −1.0%).

**Diagnostics:**
- **Reservoir grains:** 2 (300 K) … 44 (2000 K).
- **Run time on 4 cores** (40 conditions): absorbing barrier 7 s, CSE 50 s, SteadyStateOlzmann with LAPACK 68 s, TimeIntegration 658 s. The collision band reaches 716 grains at 2000 K, which makes the time integration the costliest.

## 4. What was found along the way

- **The temperature step at 304.7 K** in all methods of ZZ-allyl + O₂: the integer window of the former low-energy reduction rule. Solved by the reservoir state (`reports/low_energy_reduction_temperature_step.md`, `reports/low_energy_reservoir_state.md`).
- **A misleading README explanation.** The README said the absorbing barrier's R → IEPOX + OH is 15% above MESS's R → P5, because "these are different quantities".
  - The 15% is correct for the exact-Eckart runs (+14.8 … +15.7%).
  - With MESS's Eckart model the same quantity is −2.6 … −3.5% from MESS and agrees with CSE (G13 eq. 21) within 0.44%.
  - So the difference is the tunneling model, not the definition. Corrected in the README.
- **Wrong numbers in the method documents,** e.g. CSE R → G4 "+1.0 … +1.7%". Corrected and now re-derived from the CSV files.

## 5. References

- G13: Y. Georgievskii, J. A. Miller, M. P. Burke, S. J. Klippenstein, J. Phys. Chem. A 117, 12146 (2013).
- GO10: G. González-García, J. Olzmann, Phys. Chem. Chem. Phys. 12, 12290 (2010).
- PO14: M. Pfeifle, J. Olzmann, Int. J. Chem. Kinet. 46, 231 (2014).
- MK06: J. A. Miller, S. J. Klippenstein, J. Phys. Chem. A 110, 10528 (2006).
