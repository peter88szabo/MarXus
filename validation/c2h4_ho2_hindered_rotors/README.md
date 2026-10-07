# Hindered rotors against MESS: C₂H₄ + HO₂ → CH₂CH₂OOH → OH + oxirane

**Date:** 2026-10-07. MarXus results come from the current working tree, run on at most 4 cores.

## 1. Purpose and system

**Purpose.** This validates the one-dimensional hindered rotors of MarXus (`reports/hindered_rotors.md`) against a MESS run of the same deck.

**System.** HO₂ + C₂H₄ (P1) → CH₂CH₂OOH (W2) → OH + oxirane (P2), the HO₂ + C₂H₄ part of the C₂H₅ + O₂ example of the MESS distribution (`Examples_From_Argon/examples/test_c2h5o2/test_c2h5o2.inp`).
- The well, both barriers and the fragments are copied unchanged.
- W2 has three hindered rotors (CH₂ with σ = 2, OOH, OH), and each barrier has two (OOH, OH). All potentials are given as 12 equidistant points.
- B4's OH rotor lists 14 values for `Potential[kcal/mol] 12`. Both codes read the first 12 and ignore the rest of the line.
- Both barriers have Eckart tunneling; P2 is a `Dummy` product.

**Deck changes.**
- Temperatures reduced to 300, 500, 700, 1000, 1500, 2000 K.
- Pressures reduced to 0.01, 0.1, 1, 10, 100 bar.
- `Reactant P1`.
- `OutPrecision 6` and `LogPrecision 6`, for six printed digits.

**Energetics** (kcal/mol, relative to P1):
- W2 is at −2.36.
- The entrance barrier B3 is at 13.25 and the exit barrier B4 at 11.32.
- The well is shallow: its lowest barrier is 13.68 kcal/mol above its bottom.

**Runs.** MESS is the static binary of the 2026 source (`run_mess.sh`, 2 s). MarXus runs all four methods with two tunneling models (`run_marxus.sh`):
- `default`: the exact Eckart transmission, the MarXus default;
- `mess_eckart`: the MESS Eckart model (`--tunneling mess-eckart`). This variant isolates the rotor and state-counting treatment from the tunneling model.

## 2. Rotor levels

The MESS log prints, for every rotor, the effective rotational constant, the ground level, the number of levels and the nine lowest levels above the ground. MarXus prints the same in the `Rotors:` lines of its report (`rotor_comparison.csv`):

| species | rotor | B MESS (cm⁻¹) | B MarXus (cm⁻¹) | E₀ MESS = MarXus (kcal/mol) | levels | lowest nine levels, max \|Δ\| (kcal/mol) |
|---|---|---|---|---|---|---|
| W2 | CH₂, σ = 2 | 9.93852 | 9.93840 | 0.155905 | 997 = 997 | 4·10⁻⁶ |
| W2 | OOH | 1.94297 | 1.94294 | 0.236392 | 997 = 997 | 4·10⁻⁶ |
| W2 | OH | 19.5404 | 19.5402 | 0.300510 | 997 = 997 | 4·10⁻⁶ |
| B3 | OOH | 1.96105 | 1.96102 | 0.128274 | 997 = 997 | 4·10⁻⁶ |
| B3 | OH | 19.6833 | 19.6830 | 0.601209 | 997 = 997 | 5·10⁻⁶ |
| B4 | OOH | 1.77550 | 1.77547 | 0.138660 | 991 = 991 | 5·10⁻⁶ |
| B4 | OH | 19.9134 | 19.9132 | 0.116406 | 997 = 997 | 5·10⁻⁶ |

- **Ground levels** are identical to all printed digits.
- **The nine lowest levels** agree within half a unit of the last printed digit.
- **The rotational constants** differ by 0.9–1.5·10⁻⁵ (relative). MESS converts B to cm⁻¹ for printing with rounded atomic-unit constants (B·I = 16.8578262 cm⁻¹ amu Å², against 16.8576304 in MarXus); the identical levels show that its internal B is the same.
- **Number of levels** below the top of the basis: 997 of 999 functions, and 991 for the OOH rotor of B4 with its 37 kcal/mol barrier.

## 3. High-pressure rate coefficients

Deviation of MarXus from MESS (`high_pressure_comparison.csv`, `plots/high_pressure_deviation.png`). The first number is with the exact Eckart transmission; the second, in brackets, is with the MESS Eckart model.

| T (K) | W2 → P1 (s⁻¹) | W2 → P2 (s⁻¹) | P1 → W2 (cm³ s⁻¹) |
|---|---|---|---|
| 300 | +2.02% (+0.19%) | +1.57% (+0.13%) | +1.43% (−0.39%) |
| 500 | +1.24% (+0.14%) | +0.99% (+0.12%) | +0.82% (−0.27%) |
| 700 | +0.90% (+0.11%) | +0.73% (+0.11%) | +0.59% (−0.20%) |
| 1000 | +0.64% (+0.09%) | +0.54% (+0.11%) | +0.42% (−0.12%) |
| 1500 | +0.44% (+0.07%) | +0.40% (+0.11%) | +0.31% (−0.06%) |
| 2000 | +0.34% (+0.07%) | +0.33% (+0.12%) | +0.26% (−0.02%) |

- **The same tunneling model.** These are canonical averages of the microcanonical k(E), so they test the rotor-convolved densities of the well, the sums of states of both barriers and the fragment states. With the MESS Eckart model, MarXus agrees with MESS within 0.07–0.39% from 300 to 2000 K.
- **Exact Eckart.** Its larger transmission adds 0.3–2.0%, most at 300 K, as in the other MESS validations (`validation/ZZAllyl+O2_Gamma_Case2/`).

## 4. Phenomenological rate coefficients

**What is compared.** MarXus CSE (chemically significant eigenvalues) against the MESS species-species tables at every (T, p) (`rate_comparison.csv`, `plots/rate_deviation.png`). In the MESS tables, W2 is a species of its own only at:
- 300 and 500 K;
- 700 K above 0.01 bar;
- 1000 K at 100 bar.

Elsewhere MESS merges W2 with the reactant P1 (`ChemicalEigenvalueMax 0.2`). Since 2026-10-07 MarXus applies the same merging criterion as MESS's direct method (Section 5), and merges W2 at exactly the same 15 conditions.

**The net rate coefficient.** The net P1 → P2 (capture minus return; the long-time formation of P2) is compared at every condition. It is the P1 → P2 entry of MESS where W2 is merged, and the P1 → P1 diagonal of MESS otherwise.

Deviation of MarXus from MESS with the MESS Eckart model:

| T (K) | p (bar) | W2 | net P1 → P2 | W2 → P1 | W2 → P2 | P1 → W2 | P1 → P2 |
|---|---|---|---|---|---|---|---|
| 300 | 0.01 | both keep | +0.14% | +0.48% | +0.35% | −0.06% | +0.16% |
| 300 | 0.1–100 | both keep | +0.14% | +0.71 … +0.73% | +0.64 … +0.70% | +0.13 … +0.15% | +0.07 … +0.15% |
| 500 | 0.01 | both keep | +0.10% | −0.23% | −0.67% | −0.27% | +0.11% |
| 500 | 0.1–100 | both keep | +0.11 … +0.14% | +0.36 … +0.65% | +0.06 … +0.55% | +0.09 … +0.26% | −0.33 … +0.11% |
| 700 | 0.01 | both merge | +0.10% | | | | |
| 700 | 0.1 | both keep | +0.11% | +3.13% | +2.61% | +3.20% | +0.07% |
| 700 | 1 | both keep | +0.15% | +1.94% | +1.60% | +1.75% | −0.04% |
| 700 | 10 | both keep | +0.22% | +1.15% | +0.96% | +0.86% | −0.28% |
| 700 | 100 | both keep | +0.20% | +0.67% | +0.59% | +0.36% | −0.60% |
| 1000 | 0.01 | both merge | +0.11% | | | | |
| 1000 | 0.1 / 1 / 10 | both merge | +0.12% / +0.12% / +0.12% | | | | |
| 1000 | 100 | both keep | +0.54% | +1.67% | +1.49% | +1.69% | −0.86% |
| 1500, 2000 | all | both merge | +0.14 … +0.16% | | | | |

With the exact Eckart transmission all values are higher by the tunneling difference (net P1 → P2 +0.4 … +2.0%; `comparison_summary.txt`).

**Reading.**
1. **The net reaction** agrees within +0.10 … +0.22% at all conditions except 1000 K and 100 bar (+0.54%). This includes the 15 conditions where both codes merge W2 (+0.10 … +0.16%).
2. **The individual rate coefficients** of the well agree within 0.7% where the eigenvalue separation is good (300 K; 500 K).
   - They deviate by up to +3.2% at 700 K and 0.1 bar, where the separation λ₁/λ₂ is 0.075 in MarXus.
   - The deviation falls with pressure as the separation improves (+0.7% at 100 bar, separation 0.006).
   - In the CSE these coefficients depend on how the chemical and relaxational eigenvectors are separated, and they are sensitive to the relaxation eigenvalue (next point). The net reaction does not depend on this split.
3. **Merged conditions.** With the eigenvalue-ratio criterion (the MarXus rule until 2026-10-07, now the option `EigenvalueRatio`), MarXus kept W2 at 700 K and 0.01 bar and at 1000 K and 0.1–10 bar.
   - There the net P1 → P2 of MESS refers to the merged pool P1 + W2 and differed by up to +6.25% (1000 K, 10 bar); MESS's own value goes from 7.25·10⁻¹⁶ (10 bar, merged) to 9.77·10⁻¹⁶ cm³ s⁻¹ (100 bar, kept).
   - With the same criterion as MESS the difference is +0.12%.

## 5. Why MESS merges W2 under more conditions

### 5.1 The criterion of MESS's direct method, adopted by MarXus

**MESS** (2026 source, `MasterEquation::direct_diagonalization_method`, `mess.cc` lines 12105–12123; `CalculationMethod direct`, as in this deck) interprets `ChemicalEigenvalueMax` (`chemical_threshold`) by its range:
- **0 < value < 1:** "relaxation projection threshold". An eigenvector is chemical while its projection on the relaxational subspace, 1 − F_ne (column *P of the log), is at most the value.
- **value > 1:** "absolute eigenvalue threshold", Λ_relax/Λ ≥ value.
- **value < −1:** "relative eigenvalue threshold".

**MarXus.** Until 2026-10-07 MarXus counted the chemical eigenvalues as Λ ≤ value × Λ_{N+1}. That is the rule of MESS's reaction-complex code (`ReactiveComplex::there_are_bound_groups`, `mess.cc` line 662; `ReactiveComplex::well_reduction_method`, around line 1148), not of the direct method.

**The projections agree, and the projection rule reproduces the MESS decisions.** 1 − F_ne of the lowest eigenvector:

| T (K) | p (bar) | MESS *P | MarXus | MESS decision |
|---|---|---|---|---|
| 700 | 0.01 | 0.236 | 0.230 | merged (> 0.2) |
| 700 | 0.1 | 0.142 | 0.139 | kept |
| 700 | 1 | 0.060 | 0.058 | kept |
| 700 | 10 | 0.0137 | 0.0134 | kept |
| 700 | 100 | 0.00117 | 0.00115 | kept |
| 1000 | 0.01 | 0.772 | (merged) | merged |
| 1000 | 0.1 | 0.645 | 0.580 | merged |
| 1000 | 1 | 0.442 | 0.404 | merged |
| 1000 | 10 | 0.215 | 0.203 | merged |
| 1000 | 100 | 0.056 | 0.054 | kept |

With the MESS criterion, the MarXus projections give every merging decision of MESS in this deck, including the borderline 1000 K, 10 bar case (0.203 > 0.2).

**Adopted (2026-10-07).** MarXus now uses the relaxational projection by default (`ChemicalSubspaceCriterion RelaxationProjection`, `chemical_projection_count`); the eigenvalue ratio is the option `ChemicalSubspaceCriterion EigenvalueRatio` / `--chemical-subspace-criterion eigenvalue-ratio`. In this deck, both codes now merge W2 at the same 15 conditions (`comparison_summary.txt`).

### 5.2 Eigenvalues

Both codes divide by the same relaxation limit, the (N+1)-th eigenvalue (MESS: `min_relax_eval = eigenval[well_size()]`). λ₁/λ₂ of the lowest eigenvalue:

| p (bar) | 0.01 | 0.1 | 1 | 10 | 100 |
|---|---|---|---|---|---|
| 700 K, MESS | 0.121 | 0.086 | 0.052 | 0.023 | 0.0059 |
| 700 K, MarXus | 0.098 | 0.075 | 0.047 | 0.022 | 0.0057 |
| 1000 K, MESS | 0.465 | 0.393 | 0.285 | 0.164 | 0.074 |
| 1000 K, MarXus | (merged) | 0.197 | 0.162 | 0.119 | 0.066 |

At 1000 K and 1 bar:
- λ₁ is 1.355·10⁸ s⁻¹ in MESS and 1.505·10⁸ s⁻¹ in MarXus (+11%);
- λ₂ is 4.76·10⁸ s⁻¹ in MESS and 9.27·10⁸ s⁻¹ in MarXus.

**Grain width is excluded.** A MarXus run at 1000 K alone, with MESS's grain of 139 cm⁻¹ (0.2 kT at 1000 K instead of 0.2 kT at 300 K), gives the same ratios (0.197, 0.163, 0.119, 0.066).

**A candidate, not tested here: the bottom of the well.** There the exponential-down kernel cannot be normalized.
- MarXus lumps these grains into one thermalized reservoir state (MESMER manual, Sec. 14.2.1; `reports/low_energy_reservoir_state.md`): 15 grains (630 cm⁻¹) at 700 K and 28 grains (1176 cm⁻¹) at 1000 K, in a well 4790 cm⁻¹ deep.
- A thermalized reservoir has no slow relaxation inside it, which raises λ₂.
- MESS treats the well bottom differently.

**Effect on the results.** The difference grows with T and falls with p. It is the likely source of the deviations of the individual CSE rate coefficients at poor separation (Section 4), whereas the net reaction agrees within +0.10 … +0.22%. It concerns the collision treatment, not the rotors.

## 6. Conclusion

- **Rotor levels.** The hindered-rotor treatment of MarXus (Fourier basis of period 2π/σ, potential interpolated through the deck points, Kilpatrick–Pitzer reduced moment, levels convolved with the other degrees of freedom) reproduces the MESS levels of all seven rotors to every printed digit.
- **Rate coefficients.** With the same tunneling model, the high-pressure rate coefficients agree within 0.07–0.39%, and the net P1 → P2 rate coefficient within +0.10 … +0.22% (+0.54% at 1000 K and 100 bar).
- **Individual CSE coefficients** of the shallow well deviate by up to 3% at poor eigenvalue separation. λ₂ is about twice MESS's at 1000 K; the likely origin is the treatment of the well bottom (Section 5.2), not the rotors.
- **Merging.** MESS's direct method merges by the relaxational projection 1 − F_ne. MarXus does the same by default since 2026-10-07 (Section 5.1), and both codes merge the well at the same 15 conditions. There the net P1 → P2 agrees within +0.10 … +0.16%.

## 7. Files and how to run

| file | content |
|---|---|
| `input/c2h4_ho2.inp` | the deck (MESS format; read unchanged by both codes) |
| `run_mess.sh` | MESS run; output in `reference_mess/` (`.out` rate tables, `.log` rotor data and eigenvalues) |
| `run_marxus.sh` | MarXus, four methods × two tunneling models; output in `marxus_output/<variant>_<method>.{out,csv,_tables.csv,err}` |
| `compare_with_mess.py` | comparison; writes `rotor_comparison.csv`, `high_pressure_comparison.csv`, `rate_comparison.csv`, `plots/*.png`, and prints the summary in `comparison_summary.txt` |

```
bash run_mess.sh
bash run_marxus.sh
source ~/.venvs/science/bin/activate && python compare_with_mess.py > comparison_summary.txt
```
