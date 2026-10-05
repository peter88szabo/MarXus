# Validation: H + C₂H₂ ⇌ C₂H₃, MarXus versus stored MESS results

**Date:** 2026-10-05. The MarXus results come from the current, uncommitted working tree.

**Purpose.** Reproduce a published master-equation example with MarXus's steady-state chemical-activation solver:
- Olzmann-type steady state;
- 1 cm⁻¹ cells averaged into grains;
- exact Eckart tunneling (Miller 1979);
- exponential-down kernel.

**Reference.** The MESS example set ("Examples_From_Argon"), case `c2h3`: one well and one tight transition state, with and without Eckart tunneling. It is the model of Miller & Klippenstein (2004), H + C₂H₂ (+M) ⇌ C₂H₃ (+M); the PDF `miller2004.pdf` is in the original example directory.

## 1. Contents of this directory

| Path | Content |
|---|---|
| `input/c2h3_tight.inp` | deck run by both codes: 8 T (300–2000 K) × 5 p (0.1–10 atm), Eckart tunneling |
| `input/c2h3_tight_short.inp` | 1000 K, 1 atm, Eckart tunneling |
| `input/c2h3_tight_short_notunneling.inp` | 1000 K, 1 atm, no tunneling |
| `reference_mess_output/*.out` | stored MESS results for the same decks (unchanged copies) |
| `marxus_output/*.out` | MarXus results, from `run_marxus.sh`; `c2h3_tight_barrier_{5,3}kT.out` are the barrier-distance runs |
| `run_marxus.sh` | builds the MarXus example program and runs the three decks with reactant P1 |
| `plot_comparison.py` | reads all outputs, writes the figures and `comparison_table.csv` |
| `comparison_table.csv` | every compared number of the full deck |
| `barrier_distance_sensitivity.csv` | association for barrier distances 10, 5 and 3 kT, all conditions |
| `plots/*.png` | the figures below |

The decks are byte-identical copies of `MESS_kinetics/Examples_From_Argon/examples/c2h3/*.inp`, and MarXus reads them unchanged. The decks contain no `Reactant` line, so the reactant P1 (H + C₂H₂) is given on the command line.

To reproduce:

    ./run_marxus.sh                                                       # about 1.5 min
    source ~/.venvs/science/bin/activate && python3 plot_comparison.py

## 2. System and settings

**System.** Zero-point-corrected energies in kcal/mol:
- C₂H₃ (W1): −34.44;
- TS B1: +4.42, Eckart imaginary frequency 872.23 cm⁻¹, well depths 4.42 and 38.86;
- C₂H₂ + H (P1): 0.

**Collisions.**
- Bath gas, Lennard-Jones: ε = 6.95 / 292 cm⁻¹, σ = 2.55 / 4.36 Å, m = 4 / 75 amu.
- Energy transfer: ⟨ΔE_down⟩ = 200 (T/300 K)^0.85 cm⁻¹, ExponentCutoff 15.

| | MESS (stored results) | MarXus |
|---|---|---|
| Energy grid | 0.1 kT per temperature (EnergyStepOverTemperature 0.1). States are counted on 1 cm⁻¹, splined, and evaluated at the grid nodes. | States and all convolutions are on 1 cm⁻¹ cells. Grains are 0.1 kT(T_min) wide: 70 cm⁻¹ for the 1000 K decks, 21 cm⁻¹ for the full deck at all T. Grain values are cell averages. |
| Grid top | highest barrier + 30 kT (ExcessEnergyOverTemperature) | the same, taken at the highest temperature of the deck |
| Tunneling | semiclassical P = 1/(1+e^(−S)), S = WKB action of the Eckart potential | **exact Eckart** transmission (Miller 1979 eq. 8), convolved with the TS states (Miller eq. 9) |
| Rate coefficients | eigenvalue / chemically-significant-eigenvalue analysis | steady state (two rows below) |
| k∞ of the association | from partition functions | canonical sum of the cell numbers of states W(E) e^(−E/kT) over the reactant partition functions (tested against TST to within 0.5%) |

How MarXus obtains the rate coefficients:
- **Association.** k(P1→W1) = k∞ Φ_stab. Φ_stab comes from the **intermediate steady state** with an absorbing barrier 10 kT below the classical threshold (Pilling & Robertson 2003, eq. 44).
- **Dissociation.** k(W1→P1) = k∞,d Φ_stab, by detailed balance. k∞,d comes from the **final steady state**, which is equilibrium for thermal formation through the only channel.

## 3. Potential-energy surface

![PES](plots/pes.png)

The figure shows the stationary points of the deck, with the TS and its Eckart parameters. The dashed lines are the absorbing barriers of the intermediate steady state, 10 kT below the classical threshold:

| T | absorbing barrier |
|---|---|
| 1000 K | 9.6 kT above the well bottom |
| 1500 K | only 3 kT above the bottom |
| 2000 K | below the bottom |

## 4. Results

### 4.1 1000 K, 1 atm, with and without tunneling

![1000 K](plots/short_decks_1000K.png)

| | MESS | MarXus | Δ |
|---|---|---|---|
| **no tunneling:** k∞(P1→W1), cm³ s⁻¹ | 3.88e-11 | 3.876e-11 | −0.1% |
| k(P1→W1, 1 atm) | 2.34e-12 | 2.310e-12 | −1.3% |
| k∞(W1→P1), s⁻¹ | 2.41e5 | 2.413e5 | +0.1% |
| k(W1→P1, 1 atm) | 1.46e4 | 1.438e4 | −1.5% |
| **Eckart:** k∞(P1→W1) | 4.0860e-11 | 4.145e-11 | +1.4% |
| k(P1→W1, 1 atm) | 2.5705e-12 | 2.571e-12 | 0.0% |
| k∞(W1→P1) | 2.5451e5 | 2.580e5 | +1.4% |
| k(W1→P1, 1 atm) | 1.6022e4 | 1.600e4 | −0.1% |

The MESS values without tunneling are printed with three digits only.

### 4.2 High-pressure limits (full deck, Eckart)

![high-pressure limits](plots/high_pressure_limits.png)

| T (K) | MESS k∞(P1→W1) | MarXus | Δ | MESS k∞(W1→P1) | MarXus | Δ |
|---|---|---|---|---|---|---|
| 300 | 2.4433e-13 | 2.5602e-13 | +4.8% | 9.5862e-16 | — | — |
| 500 | 3.5805e-12 | 3.6855e-12 | +2.9% | 2.9835e-04 | — | — |
| 750 | 1.6985e-11 | 1.7316e-11 | +1.9% | 2.4393e+02 | 2.4851e+02 | +1.9% |
| 1000 | 4.0861e-11 | 4.1451e-11 | +1.4% | 2.5451e+05 | 2.5804e+05 | +1.4% |
| 1250 | 7.3151e-11 | 7.3985e-11 | +1.1% | 1.7464e+07 | 1.7655e+07 | +1.1% |
| 1500 | 1.1180e-10 | 1.1285e-10 | +0.9% | 3.0096e+08 | 3.0364e+08 | +0.9% |
| 1750 | 1.5523e-10 | 1.5646e-10 | +0.8% | 2.3332e+09 | 2.3507e+09 | +0.7% |
| 2000 | 2.0230e-10 | — | — | 1.0932e+10 | 1.1002e+10 | +0.6% |

**Why some MarXus values are missing.**
- **Dissociation at 300 and 500 K.** The final steady state of the 34 kcal/mol well without a sink is numerically singular in double precision (the lines marked "not available" in `marxus_output/c2h3_tight.out`). This is the postponed double-precision topic.
- **Association at 2000 K.** The example program prints k∞ only next to an intermediate steady state, and that steady state is not defined at 2000 K (§5.3).

**Both directions deviate by the same amount**, as detailed balance requires: both codes use the same equilibrium constant. The deviation grows from +0.6% at 2000 K to +4.8% at 300 K. This is the exact Eckart transmission (MarXus) against the semiclassical one (MESS); see §5.2.

### 4.3 Fall-off (full deck, Eckart)

![fall-off](plots/falloff_P1_W1.png)

![deviation](plots/deviation.png)

Association: k∞(T) Φ_stab, from the intermediate steady state. Dissociation: k∞,d(T) Φ_stab. k∞,d is pressure independent, so it is taken from any temperature at which the MarXus final steady state exists (750–2000 K).

| T (K) | p (atm) | MESS k(P1→W1) | MarXus | Δ | MESS k(W1→P1) | MarXus | Δ |
|---|---|---|---|---|---|---|---|
| 300 | 0.1 | 1.4421e-13 | 1.5261e-13 | +5.8% | 5.6580e-16 | — | — |
| 300 | 0.3 | 1.8112e-13 | 1.9113e-13 | +5.5% | 7.1063e-16 | — | — |
| 300 | 1 | 2.1202e-13 | 2.2310e-13 | +5.2% | 8.3187e-16 | — | — |
| 300 | 3 | 2.2930e-13 | 2.4083e-13 | +5.0% | 8.9965e-16 | — | — |
| 300 | 10 | 2.3870e-13 | 2.5040e-13 | +4.9% | 9.3652e-16 | — | — |
| 500 | 0.1 | 6.9740e-13 | 7.2366e-13 | +3.8% | 5.8116e-05 | — | — |
| 500 | 0.3 | 1.1408e-12 | 1.1826e-12 | +3.7% | 9.5065e-05 | — | — |
| 500 | 1 | 1.7523e-12 | 1.8139e-12 | +3.5% | 1.4601e-04 | — | — |
| 500 | 3 | 2.3395e-12 | 2.4184e-12 | +3.4% | 1.9494e-04 | — | — |
| 500 | 10 | 2.8901e-12 | 2.9832e-12 | +3.2% | 2.4082e-04 | — | — |
| 750 | 0.1 | 7.9550e-13 | 8.0714e-13 | +1.5% | 1.1429e+01 | 1.1584e+01 | +1.4% |
| 750 | 0.3 | 1.5805e-12 | 1.6055e-12 | +1.6% | 2.2704e+01 | 2.3042e+01 | +1.5% |
| 750 | 1 | 3.0644e-12 | 3.1164e-12 | +1.7% | 4.4016e+01 | 4.4726e+01 | +1.6% |
| 750 | 3 | 5.1099e-12 | 5.2010e-12 | +1.8% | 7.3391e+01 | 7.4643e+01 | +1.7% |
| 750 | 10 | 8.0055e-12 | 8.1542e-12 | +1.9% | 1.1497e+02 | 1.1703e+02 | +1.8% |
| 1000 | 0.1 | 5.1285e-13 | 5.0987e-13 | −0.6% | 3.1993e+03 | 3.1741e+03 | −0.8% |
| 1000 | 0.3 | 1.1451e-12 | 1.1415e-12 | −0.3% | 7.1401e+03 | 7.1064e+03 | −0.5% |
| 1000 | 1 | 2.5737e-12 | 2.5735e-12 | −0.0% | 1.6042e+04 | 1.6020e+04 | −0.1% |
| 1000 | 3 | 4.9976e-12 | 5.0107e-12 | +0.3% | 3.1142e+04 | 3.1193e+04 | +0.2% |
| 1000 | 10 | 9.3954e-12 | 9.4462e-12 | +0.5% | 5.8533e+04 | 5.8805e+04 | +0.5% |
| 1250 | 0.1 | 2.7659e-13 | 2.6677e-13 | −3.6% | 6.6421e+04 | 6.3658e+04 | −4.2% |
| 1250 | 0.3 | 6.6713e-13 | 6.4696e-13 | −3.0% | 1.6000e+05 | 1.5438e+05 | −3.5% |
| 1250 | 1 | 1.6573e-12 | 1.6173e-12 | −2.4% | 3.9695e+05 | 3.8593e+05 | −2.8% |
| 1250 | 3 | 3.5798e-12 | 3.5140e-12 | −1.8% | 8.5656e+05 | 8.3854e+05 | −2.1% |
| 1250 | 10 | 7.6924e-12 | 7.5985e-12 | −1.2% | 1.8389e+06 | 1.8132e+06 | −1.4% |
| 1500 | 0.1 | 1.4204e-13 | 1.2299e-13 | −13.4% | 3.9042e+05 | 3.3094e+05 | −15.2% |
| 1500 | 0.3 | 3.6165e-13 | 3.1837e-13 | −12.0% | 9.9060e+05 | 8.5664e+05 | −13.5% |
| 1500 | 1 | 9.6459e-13 | 8.6600e-13 | −10.2% | 2.6318e+06 | 2.3301e+06 | −11.5% |
| 1500 | 3 | 2.2504e-12 | 2.0588e-12 | −8.5% | 6.1189e+06 | 5.5396e+06 | −9.5% |
| 1500 | 10 | 5.3391e-12 | 4.9869e-12 | −6.6% | 1.4469e+07 | 1.3418e+07 | −7.3% |
| 1750 | 0.1 | 7.4261e-14 | 4.2763e-14 | −42.4% | 1.1890e+06 | 6.4249e+05 | −46.0% |
| 1750 | 0.3 | 1.9658e-13 | 1.1847e-13 | −39.7% | 3.1213e+06 | 1.7800e+06 | −43.0% |
| 1750 | 1 | 5.5199e-13 | 3.5240e-13 | −36.2% | 8.6792e+06 | 5.2945e+06 | −39.0% |
| 1750 | 3 | 1.3623e-12 | 9.2275e-13 | −32.3% | 2.1226e+07 | 1.3864e+07 | −34.7% |
| 1750 | 10 | 3.4783e-12 | 2.5279e-12 | −27.3% | 5.3670e+07 | 3.7979e+07 | −29.2% |
| 2000 | 0.1 | 4.1158e-14 | — | — | 2.5837e+06 | — | — |
| 2000 | 0.3 | 1.1207e-13 | — | — | 6.9238e+06 | — | — |
| 2000 | 1 | 3.2672e-13 | — | — | 1.9798e+07 | — | — |
| 2000 | 3 | 8.4038e-13 | — | — | 4.9970e+07 | — | — |
| 2000 | 10 | 2.2655e-12 | — | — | 1.3188e+08 | — | — |

## 5. Interpretation

### 5.1 750–1250 K: agreement

| | deviation from MESS |
|---|---|
| association | −3.6% … +1.9% |
| dissociation | −4.2% … +1.8% |
| both, at 1000 K | ±0.8% |

The two codes differ in graining (cell-averaged grains versus nodes), in method (steady state versus eigenvalue analysis) and in tunneling (exact versus semiclassical), yet still agree this closely.

### 5.2 300–500 K: +3% to +6% in the association, from the tunneling model

The deviation is already present in k∞ (+4.8% at 300 K, +2.9% at 500 K). The fall-off ratio k/k∞ agrees to within 1%: at 300 K and 0.1 atm it is 0.590 in MESS and 0.596 in MarXus.

MarXus uses the exact Eckart transmission, as decided. MESS's semiclassical form gives a lower transmission at low T. For a Case1 H-shift we found 13% at 300 K; here the barrier is low and broad (4.42 kcal/mol, 872 cm⁻¹), so the difference is smaller.

### 5.3 1500 K and above: the absorbing-barrier picture breaks down

| T | well depth below the TS (13 591 cm⁻¹) | absorbing barrier above the well bottom |
|---|---|---|
| 1500 K | 13.0 kT | 3.0 kT |
| 1750 K | 11.2 kT | 1.2 kT |
| 2000 K | 9.8 kT | below the bottom |

The thermal distribution of C₂H₃ has ⟨E⟩ ≈ 4 kT (2935 cm⁻¹ at 1000 K, MarXus final steady state). At these temperatures the barrier therefore cuts through the thermal distribution. Thermalized molecules above the barrier are not counted as stabilized, so k(P1→W1) and k(W1→P1) come out too low:
- 1500 K: −7% to −15%;
- 1750 K: −27% to −46%.

MESS obtains these rate coefficients from the eigenvalue analysis, which needs no barrier.

At 2000 K MarXus now refuses the intermediate steady state with an explanatory error. It previously clamped the barrier to the well bottom without warning and returned zero stabilization.

### 5.4 Decision for shallow wells (Peter, 2026-10-05)

- **(b) The user chooses the absorbing-barrier distance.** The default stays 10 kT below the lowest threshold (Pilling & Robertson 2003; Carstensen & Dean 2007). The library accepts any distance (`AbsorbingBarrier::BelowLowestThreshold { kt_multiple }`) or explicit grains (`AtGrains`). The example program takes `--barrier-kt X`. If the barrier would lie below the well bottom, the error message suggests a smaller distance.
- **Eigenvalue route.** It is an optional choice, `SteadyState::EigenvalueAnalysis` (`--steady-state eigenvalue`), but **it is not available yet**. Selecting it is reported as "not available yet" and is never replaced by another method. MarXus currently has no eigenvalue analysis of J; there is only a general Jacobi diagonalizer in `numeric/jacobi_diag.rs`.

### 5.5 Sensitivity to the absorbing-barrier distance

![barrier distance](plots/barrier_distance_sensitivity.png)

`run_marxus.sh` also runs the full deck with the barrier 5 and 3 kT below the threshold (`marxus_output/c2h3_tight_barrier_5kT.out`, `..._3kT.out`). All 40 conditions are in `barrier_distance_sensitivity.csv`.

The 1 atm rows (MarXus k(P1→W1) relative to MESS):

| T (K) | MESS k(P1→W1) | 10 kT (default) | 5 kT | 3 kT |
|---|---|---|---|---|
| 300 | 2.1203e-13 | +5.2% | +5.2% | +5.3% |
| 500 | 1.7523e-12 | +3.5% | +3.6% | +3.8% |
| 750 | 3.0644e-12 | +1.7% | +1.9% | +2.5% |
| 1000 | 2.5737e-12 | −0.0% | +0.6% | +2.0% |
| 1250 | 1.6573e-12 | −2.4% | −0.8% | +2.0% |
| 1500 | 9.6459e-13 | −10.2% | −3.4% | +1.5% |
| 1750 | 5.5199e-13 | −36.2% | −9.9% | −1.8% |
| 2000 | 3.2672e-13 | not defined | −26.8% | −10.0% |

**At low T (≤ 750 K) the result does not depend on the barrier distance.** The steady state lies on a plateau: the change is less than 1% between 10 and 3 kT, so the remaining +2…+5% is the tunneling-model difference (§5.2).

**At high T a smaller distance recovers MESS.** The barrier then lies above the thermal distribution of the shallow well. With 3 kT the agreement is within 2% up to 1750 K.

At 1000–1250 K, 3 kT overshoots slightly (+2%). The distance should therefore be chosen per condition: as large as the well allows while staying above its thermal distribution, and checked for a plateau by varying it.

At 2000 K (well depth 9.8 kT) even 3 kT remains 10% low. That is the regime of the eigenvalue route, which is not available yet.

### 5.6 Still open

- **Final steady state of deep wells without a sink at low T** (double precision): postponed by Peter.
- **Eigenvalue route:** to be implemented later.

## 6. Code changes made for this validation (MarXus, uncommitted)

- **Input reader.**
  - Reads `Atom` fragments.
  - Accepts `GroundEnergy` in any unit.
  - Bug fix: block keywords are now matched exactly. `WellCutoff 10` had been read as a Well block named "10" that swallowed the whole model.
- **Entrance high-pressure rate coefficient** (`EntranceHighPressureRate`). Tested against the ILT input and against canonical TST.
- **Example program** `chemical_activation_from_deck`:
  - optional reactant argument;
  - bimolecular rate coefficients k(R→X) = k∞ Φ_X;
  - every (T, p) is solved separately; a condition without a valid steady state is reported as "not available";
  - `--barrier-kt X` sets the absorbing-barrier distance;
  - `--steady-state intermediate|final|eigenvalue|both` selects the solution.
- **`SteadyState::EigenvalueAnalysis`**: a selectable option that is reported as not available yet.
- **Absorbing barrier below the well bottom.** Now an error instead of a silent clamp.

Full test suite: 117 library and 18 binary tests pass.
