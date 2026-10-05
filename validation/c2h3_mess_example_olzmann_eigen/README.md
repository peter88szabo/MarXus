# H + C₂H₂ ⇌ C₂H₃: Olzmann's eigenvalue analysis in MarXus compared with MESS

**MarXus, 2026-10-05**. Updated after Peter's decisions of 20:01 (k_uni = eigenvector average; the sum rule warns and never rejects; threshold 1.5%) and of 20:18 (shifted factorization in the inverse iteration). All runs on at most 4 cores.

This directory repeats the comparison of `../c2h3_mess_example/` with the **eigenvalue analysis** (Olzmann's solution method). That earlier comparison used the absorbing-barrier route (intermediate steady state) and is kept unchanged.

The pressure-dependent unimolecular rate coefficient comes from the thermal eigenvector of the relaxation matrix J of the final steady state, **without an absorbing barrier**. The association follows by detailed balance. Same decks, same MESS reference output, same MarXus adapter as before; only the solution method differs.

## 1. Files

| file | content |
|---|---|
| `input/*.inp` | the three decks (copies of `../c2h3_mess_example/input`) |
| `reference_mess_output/*.out` | the stored MESS results (copies) |
| `run_marxus.sh` | runs `chemical_activation_from_deck --steady-state eigenvalue` with the eigen-solvers |
| `marxus_output/<deck>_<solver>.out` | MarXus results, solver = `inverse` (default), `lapack`, `full` |
| `plot_comparison.py` | reads all of the above (and, read-only, `../c2h3_mess_example/comparison_table.csv`, `../c2h3_mess_example/barrier_distance_sensitivity.csv`), writes the plots and `comparison_table.csv` |
| `comparison_table.csv` | every compared number (Section 4) |
| `plots/pes.png` | stationary points of the deck; no absorbing barrier in this route |
| `plots/falloff_W1_P1.png` | k_uni(T, p) against MESS k(W1→P1), with the high-pressure limits |
| `plots/deviation.png` | MarXus/MESS − 1: dissociation (k_uni) and association (k_uni·K) |
| `plots/eigen_vs_absorbing_barrier.png` | association deviation: eigenvalue route versus absorbing barrier 10, 5, 3 kT |
| `plots/sum_rule.png` | sum-rule deviation of inverse iteration and LAPACK, and λ₂/k_uni |
| `plots/short_decks_1000K.png` | 1000 K, 1 atm, with and without tunneling, all three solvers |

Reproduce:

```
bash run_marxus.sh                     # ~2.5 min (LAPACK on the full deck dominates; 2634 grains)
source ~/.venvs/science/bin/activate && python3 plot_comparison.py
```

## 2. Method

**References.**
- GO10: González-García, Olzmann, Phys. Chem. Chem. Phys. 12, 12290 (2010).
- PO14: Pfeifle, Olzmann, Int. J. Chem. Kinet. 46, 231 (2014).
- O02: Olzmann, Phys. Chem. Chem. Phys. 4, 3614 (2002).

**Operator.** J = ω(I − P) + K is the final-steady-state operator of PO14 eq. 2, with exact detailed balance of P and no sink, since this system has none. It is symmetrized: S = D⁻¹JD, D = diag(√f), f = ρ e^{−E/kT}.

**Reported rate coefficient.** k_uni = Σᵢ k(Eᵢ) Ñᵢ: the specific rate coefficients averaged over Ñ, the normalized eigenvector of the lowest eigenvalue λ₁. GO10 calls this the averaging procedure "analogous to eqn (9) but with Ñs = Ñs^th being the normalized eigenvector associated with the lowest eigenvalue λ₁" (text after eq. 12).

**λ₁ and the sum rule.** λ₁, the lowest eigenvalue of J (GO10 eq. 12), is printed beside k_uni. The column sums of J are the loss rates, so λ₁ = k_uni exactly.
- **λ₁ degrades first.** Numerically λ₁ has an absolute error of order ε‖S‖ and is lost first when it lies many orders of magnitude below the collision frequency.
- **The eigenvector is much less sensitive**, because its error scales with ε‖S‖/λ₂.
- **The output explains this** with the reference, in comment lines before the table.
- **Precision floor:** ε·max Sᵢᵢ (column `precision_floor`, 1.0·10⁻² s⁻¹ for this deck) is the order of the absolute rounding error of λ₁ (Weyl 1912). A λ₁ near or below it is rounding noise, even a negative one. This is not a merging of eigenvalues, which would show as a small λ₂/k_uni.
- **Sum-rule deviation:** |λ₁ − k_uni|/k_uni is printed for every condition. Above 1.5% (option `--sum-rule-tolerance`) a `# warning` line, also on stderr, gives the numbers and what to do. Nothing is rejected.

**Association by detailed balance.** k(P1→W1, T, p) = k_uni · k∞,assoc/k∞,diss:
- k∞,assoc is the canonical capture rate of the entrance channel;
- k∞,diss is the Boltzmann average of k(E) over the well.

**Solvers** (option `--eigen-solver`):

| solver | method | cost |
|---|---|---|
| `inverse` (default) | inverse iteration with the banded Cholesky factor of S + σI, σ = n·ε·max Sᵢᵢ (inverse iteration with a shift, Numerical Recipes §11.7; same eigenvectors, eigenvalues θ − σ); λ₂ by deflation | O(n·bw²) |
| `lapack` | LAPACK DSYEVD: Householder reduction + divide and conquer, all eigenpairs (`src/numeric/lapack_interface.rs`, system OpenBLAS) | O(n³), blocked |
| `full` | in-house Householder (tred2) + implicit QL (tql2), all eigenpairs, Olzmann's route | O(n³), unblocked |

**What MarXus tells the user to do** (in the messages themselves):

| situation | message | what to do |
|---|---|---|
| sum-rule deviation > 1.5% | `# warning`; result kept | use k_uni; confirm with the inverse iteration (most accurate for the thermal eigenvector); if the solvers give different k_uni, the condition is beyond double precision |
| λ₁ ≤ 0 | `# warning`, saying a non-positive λ₁ is impossible for J (GO10 before eq. 12): rounding noise below the double-precision floor (printed); not a merging of eigenvalues (λ₂/k_uni printed); result kept | as above |
| inverse iteration: Cholesky factor of S + σI does not exist (did not occur here) | `# not available` | use a full decomposition, `--eigen-solver lapack`, which does not need the factor |
| k_uni ≤ 0 | error | nothing leaves the network (no product channel or sink), or the eigenvector is not resolved |

## 3. Results in brief

1. **1000 K is reproduced.** k_uni agrees with MESS within −0.8 … +0.5% and the association within −0.6 … +0.5%. The short deck at 1 atm gives:
   - k_uni = 16004.0 s⁻¹ against MESS 16022.4 (−0.11%);
   - k(P1→W1) = 2.57086e-12 against 2.57053e-12 cm³ s⁻¹ (+0.01%).

   All three solvers give the same 7 digits.
2. **The absorbing-barrier problem at high T is gone.** For the association (`plots/eigen_vs_absorbing_barrier.png`):
   - The barrier route needed a user-chosen barrier distance and still failed above ~1500 K: −13% at 1500 K / 0.1 atm with 10 kT, −42% at 1750 K, no result at 2000 K.
   - The eigenvalue route needs no choice. It is within ±2.5% from 750 to 1750 K at all pressures, and +3 … +7% at 2000 K.
3. **At 300–1000 K both solution methods give the same association to 0.01%.** Eigenvalue route (k_uni·K) vs absorbing barrier at 10 kT, at every pressure:

| T (K) | eigenvalue route | barrier route |
|---|---|---|
| 300 | +5.83 / +5.52 / +5.22 / +5.03 / +4.89% | +5.83 / +5.52 / +5.22 / +5.03 / +4.90% |
| 500 | +3.76 … +3.22% | +3.77 … +3.22% |
| 750 | +1.45 … +1.86% | +1.46 … +1.86% |
| 1000 | −0.60 … +0.54% | −0.58 … +0.54% |

   The offset from MESS at low T is the known exact- vs semiclassical-Eckart difference of the high-pressure limits (`../c2h3_mess_example/README.md`).
4. **300 K: complete results from the thermal eigenvector, while λ₁ is noise.**
   - **λ₁ is far below the resolution.** It is about 10⁻¹⁵ s⁻¹, 13 orders below the double-precision floor ε·max Sᵢᵢ = 1.0·10⁻² s⁻¹.
   - **λ₁ is noise.** Its computed values are ±10⁻⁸ … 10⁻⁶ s⁻¹ and their sign changes between solvers. The warnings say so.
   - **The plain Cholesky factor of S does not exist** at 0.1, 3 and 10 atm. With the shift σ = n·ε·max Sᵢᵢ = 26 s⁻¹, the inverse iteration gives all five pressures.
   - **k_uni agrees with LAPACK** to 1·10⁻⁷ … 1·10⁻⁴.
   - **Against MESS** it is +5.66, +5.35, +5.06, +4.86, +4.73%. The association equals the barrier route to 0.001 percentage points (point 3).

5. **High temperature: λ₂/k_uni approaches 10.** Above ~1500 K the dissociation (−2 … −8%) and association (−2 … +7%) deviate from MESS with opposite trends. MarXus's equilibrium constants agree with MESS to 0.1% at every T (k∞,assoc/k∞,diss), but MESS's own pressure-dependent pair departs from its own high-pressure ratio (Section 6).

## 4. Full table

All rows: inverse iteration with the shifted factorization. "dev." is MarXus/MESS − 1. Column "absorbing barrier 10 kT" is the association deviation of the earlier route (`../c2h3_mess_example/comparison_table.csv`).

| T (K) | p (atm) | MESS k(W1→P1) (s⁻¹) | MarXus k_uni | dev. | λ₁ | sum-rule dev. | warning | λ₂/k_uni | MESS k(P1→W1) (cm³ s⁻¹) | MarXus k_uni·K | dev. | absorbing barrier 10 kT dev. |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 300 | 0.1 | 5.658e-16 | 5.978e-16 | +5.66% | -1.69e-07 | 2.8e+08 | yes | 1.45e+23 | 1.442e-13 | 1.526e-13 | +5.83% | +5.83% |
| 300 | 0.3 | 7.106e-16 | 7.487e-16 | +5.35% | 2.35e-07 | 3.1e+08 | yes | 3.46e+23 | 1.811e-13 | 1.911e-13 | +5.52% | +5.52% |
| 300 | 1 | 8.319e-16 | 8.739e-16 | +5.06% | 3.16e-07 | 3.6e+08 | yes | 9.89e+23 | 2.12e-13 | 2.231e-13 | +5.22% | +5.22% |
| 300 | 3 | 8.996e-16 | 9.434e-16 | +4.86% | -3.51e-08 | 3.7e+07 | yes | 2.75e+24 | 2.293e-13 | 2.408e-13 | +5.03% | +5.03% |
| 300 | 10 | 9.365e-16 | 9.809e-16 | +4.73% | -2.84e-06 | 2.9e+09 | yes | 8.81e+24 | 2.387e-13 | 2.504e-13 | +4.90% | +4.90% |
| 500 | 0.1 | 5.812e-05 | 6.024e-05 | +3.65% | 6.03e-05 | 1.1e-03 | no | 7.46e+11 | 6.974e-13 | 7.236e-13 | +3.76% | +3.77% |
| 500 | 0.3 | 9.506e-05 | 9.844e-05 | +3.56% | 9.86e-05 | 1.7e-03 | no | 1.37e+12 | 1.141e-12 | 1.183e-12 | +3.66% | +3.66% |
| 500 | 1 | 0.000146 | 0.000151 | +3.41% | 0.000152 | 6.0e-03 | no | 2.98e+12 | 1.752e-12 | 1.814e-12 | +3.52% | +3.52% |
| 500 | 3 | 0.0001949 | 0.0002013 | +3.27% | 0.000202 | 3.8e-03 | no | 6.7e+12 | 2.339e-12 | 2.418e-12 | +3.37% | +3.37% |
| 500 | 10 | 0.0002408 | 0.0002483 | +3.12% | 0.000246 | 1.1e-02 | no | 1.81e+13 | 2.89e-12 | 2.983e-12 | +3.22% | +3.22% |
| 750 | 0.1 | 11.43 | 11.58 | +1.33% | 11.6 | 1.0e-09 | no | 2.44e+06 | 7.955e-13 | 8.07e-13 | +1.45% | +1.46% |
| 750 | 0.3 | 22.7 | 23.04 | +1.48% | 23 | 1.1e-09 | no | 3.67e+06 | 1.58e-12 | 1.605e-12 | +1.57% | +1.58% |
| 750 | 1 | 44.02 | 44.72 | +1.61% | 44.7 | 9.2e-09 | no | 6.3e+06 | 3.064e-12 | 3.116e-12 | +1.69% | +1.70% |
| 750 | 3 | 73.39 | 74.64 | +1.70% | 74.6 | 2.6e-08 | no | 1.13e+07 | 5.11e-12 | 5.201e-12 | +1.78% | +1.78% |
| 750 | 10 | 115 | 117 | +1.78% | 117 | 4.8e-08 | no | 2.41e+07 | 8.006e-12 | 8.154e-12 | +1.86% | +1.86% |
| 1000 | 0.1 | 3199 | 3173 | -0.81% | 3.17e+03 | 5.4e-12 | no | 6.93e+03 | 5.128e-13 | 5.098e-13 | -0.60% | -0.58% |
| 1000 | 0.3 | 7140 | 7105 | -0.49% | 7.11e+03 | 3.5e-12 | no | 9.22e+03 | 1.145e-12 | 1.141e-12 | -0.32% | -0.31% |
| 1000 | 1 | 1.604e+04 | 1.602e+04 | -0.14% | 1.6e+04 | 1.6e-11 | no | 1.35e+04 | 2.574e-12 | 2.573e-12 | -0.02% | -0.01% |
| 1000 | 3 | 3.114e+04 | 3.119e+04 | +0.16% | 3.12e+04 | 1.6e-11 | no | 2.07e+04 | 4.998e-12 | 5.01e-12 | +0.26% | +0.26% |
| 1000 | 10 | 5.853e+04 | 5.88e+04 | +0.46% | 5.88e+04 | 6.3e-11 | no | 3.63e+04 | 9.395e-12 | 9.446e-12 | +0.54% | +0.54% |
| 1250 | 0.1 | 6.642e+04 | 6.441e+04 | -3.03% | 6.44e+04 | 2.4e-13 | no | 298 | 2.766e-13 | 2.699e-13 | -2.41% | -3.55% |
| 1250 | 0.3 | 1.6e+05 | 1.559e+05 | -2.56% | 1.56e+05 | 1.0e-13 | no | 364 | 6.671e-13 | 6.534e-13 | -2.06% | -3.02% |
| 1250 | 1 | 3.969e+05 | 3.889e+05 | -2.02% | 3.89e+05 | 6.0e-13 | no | 477 | 1.657e-12 | 1.63e-12 | -1.65% | -2.41% |
| 1250 | 3 | 8.566e+05 | 8.435e+05 | -1.53% | 8.43e+05 | 1.3e-12 | no | 648 | 3.58e-12 | 3.535e-12 | -1.26% | -1.84% |
| 1250 | 10 | 1.839e+06 | 1.821e+06 | -1.00% | 1.82e+06 | 2.9e-13 | no | 978 | 7.692e-12 | 7.629e-12 | -0.82% | -1.22% |
| 1500 | 0.1 | 3.904e+05 | 3.745e+05 | -4.09% | 3.74e+05 | 4.6e-14 | no | 51.9 | 1.42e-13 | 1.392e-13 | -2.01% | -13.41% |
| 1500 | 0.3 | 9.906e+05 | 9.544e+05 | -3.66% | 9.54e+05 | 1.8e-13 | no | 59.7 | 3.617e-13 | 3.547e-13 | -1.92% | -11.97% |
| 1500 | 1 | 2.632e+06 | 2.549e+06 | -3.15% | 2.55e+06 | 1.0e-14 | no | 72.4 | 9.646e-13 | 9.473e-13 | -1.79% | -10.22% |
| 1500 | 3 | 6.119e+06 | 5.957e+06 | -2.65% | 5.96e+06 | 9.7e-14 | no | 90.1 | 2.25e-12 | 2.214e-12 | -1.63% | -8.51% |
| 1500 | 10 | 1.447e+07 | 1.417e+07 | -2.09% | 1.42e+07 | 3.5e-14 | no | 122 | 5.339e-12 | 5.265e-12 | -1.39% | -6.60% |
| 1750 | 0.1 | 1.189e+06 | 1.136e+06 | -4.42% | 1.14e+06 | 7.8e-13 | no | 18.5 | 7.426e-14 | 7.564e-14 | +1.86% | -42.41% |
| 1750 | 0.3 | 3.121e+06 | 2.996e+06 | -4.02% | 3e+06 | 3.0e-13 | no | 20.4 | 1.966e-13 | 1.994e-13 | +1.43% | -39.73% |
| 1750 | 1 | 8.679e+06 | 8.37e+06 | -3.56% | 8.37e+06 | 2.0e-14 | no | 23.5 | 5.52e-13 | 5.571e-13 | +0.93% | -36.16% |
| 1750 | 3 | 2.123e+07 | 2.057e+07 | -3.11% | 2.06e+07 | 5.0e-13 | no | 27.6 | 1.362e-12 | 1.369e-12 | +0.48% | -32.27% |
| 1750 | 10 | 5.367e+07 | 5.227e+07 | -2.61% | 5.23e+07 | 1.5e-14 | no | 34.5 | 3.478e-12 | 3.479e-12 | +0.02% | -27.32% |
| 2000 | 0.1 | 2.584e+06 | 2.382e+06 | -7.81% | 2.38e+06 | 3.6e-13 | no | 9.81 | 4.116e-14 | 4.41e-14 | +7.15% | – |
| 2000 | 0.3 | 6.924e+06 | 6.434e+06 | -7.07% | 6.43e+06 | 1.7e-12 | no | 10.5 | 1.121e-13 | 1.191e-13 | +6.29% | – |
| 2000 | 1 | 1.98e+07 | 1.857e+07 | -6.20% | 1.86e+07 | 5.4e-13 | no | 11.6 | 3.267e-13 | 3.438e-13 | +5.23% | – |
| 2000 | 3 | 4.997e+07 | 4.728e+07 | -5.38% | 4.73e+07 | 9.6e-14 | no | 13 | 8.404e-13 | 8.754e-13 | +4.16% | – |
| 2000 | 10 | 1.319e+08 | 1.26e+08 | -4.47% | 1.26e+08 | 3.0e-13 | no | 15.3 | 2.265e-12 | 2.332e-12 | +2.95% | – |

**Short decks, 1000 K, 1 atm (all three solvers give these digits):**

| deck | quantity | MarXus | MESS | dev. | sum-rule dev. (inverse / LAPACK / QL) |
|---|---|---|---|---|---|
| Eckart tunneling | k_uni (s⁻¹) | 16004.04 | 16022.4 | −0.11% | 1.0e-11 / 1.4e-10 / 1.7e-10 |
| Eckart tunneling | k(P1→W1) (cm³ s⁻¹) | 2.570858e-12 | 2.57053e-12 | +0.01% | |
| no tunneling | k_uni (s⁻¹) | 14379.14 | 1.46e4 (3 digits printed) | −1.5% | 1.1e-11 / 2.2e-10 / 2.6e-10 |
| no tunneling | k(P1→W1) (cm³ s⁻¹) | 2.309837e-12 | 2.34e-12 (3 digits printed) | −1.3% | |

## 5. Numerical behaviour of the solvers

**Sum-rule deviation |λ₁ − k_uni|/k_uni** (`plots/sum_rule.png`):

| T (K) | inverse iteration | LAPACK DSYEVD |
|---|---|---|
| 750 | 4e-9 … 7e-8 | 5e-8 … 1e-5 |
| 1000 | 3e-12 … 4e-11 | 2e-10 … 7e-8 |
| ≥ 1250 | 2e-14 … 2e-12 | 2e-11 … 9e-10 |
| 500 | 1.1e-3 … 1.1e-2 (no warning) | 1.3e-3 … 11 (warnings at 4 pressures; λ₁ < 0 at 2) |
| 300 | 4e7 … 3e9, all five pressures (warnings; λ₁ < 0 at 3) | 8e7 … 1e11 (warnings; λ₁ < 0 at 2) |

**Why λ₁ is most accurate from inverse iteration.**
- A dense decomposition has a *normwise* backward error ε‖S‖, so the absolute error of λ₁ is of order ε·max(ω + k(E)) (Weyl's inequality).
- The banded Cholesky factorization has a *componentwise* one (Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., ch. 10).

**Why k_uni is robust.** The eigenvector error is of order ε‖S‖/(λ₂ − λ₁), and λ₂ ≫ λ₁. For C₂H₃, k_uni agrees between inverse iteration and LAPACK to 7 digits at 300 and 500 K, even where LAPACK's λ₁ is negative.

**This robustness is not guaranteed in general.** On an extreme synthetic test well (threshold 6000 cm⁻¹, 175 K, unit test `a_violated_sum_rule_gives_a_warning_and_no_error`), LAPACK's k_uni differed from the inverse iteration's by 10%. That is why the warning asks for a cross-check with inverse iteration.

**Timing.** Full deck (2634 grains), two conditions:

| solver | time | memory |
|---|---|---|
| inverse | 6.3 s | 137 MB |
| LAPACK | 5.0 s | 391 MB |
| QL | 238 s | 269 MB |

**Shifted factorization (implemented after Peter's go-ahead, 20:18).**
- **What it does.** The inverse iteration factors S + σI with σ = n·ε·max Sᵢᵢ (26 s⁻¹ here). This is inverse iteration with a shift, Numerical Recipes §11.7. The matrix has the same eigenvectors, and its Cholesky factor exists whenever the smallest eigenvalue relative to the diagonal exceeds a multiple of n·ε (Higham 2002, ch. 10; Demmel 1989).
- **Convergence is unaffected.** σ ≪ λ₂ (λ₂ ≈ 10⁸ s⁻¹ at 300 K), so the ratio (λ₁ + σ)/(λ₂ + σ) stays tiny.
- **k_uni is unchanged** to 7 digits wherever it existed before. Only λ₁, which is noise there, changes at 300 and 500 K.
- **Limitation.** If several lowest eigenvalues lie within σ of each other, the iteration cannot separate them. This happens with the stepladder model when the step is several grains wide: J then splits into ΔE_SL/ΔE independent sub-equations (O02 p. 3616: "solved separately"). Exponential down, the default collision model since 20:25, is not affected.

## 6. High temperature: MESS's pressure-dependent pair and detailed balance

**MarXus's equilibrium constants agree with MESS.** MarXus's high-pressure limits are +0.6 … +2.9% from MESS, larger at low T because of exact vs semiclassical Eckart tunneling. Their ratio k∞,assoc/k∞,diss agrees with MESS's to 0.04–0.10% at every T. MarXus imposes k(P1→W1) = k_uni·K.

**MESS's own ratio drifts from its K.** In the MESS output, [k(P1→W1)/k(W1→P1)] / [k∞(P1→W1)/k∞(W1→P1)] − 1 is:

| T (K) | 0.1 atm | 1 atm | 10 atm | MarXus λ₂/k_uni (1 atm) |
|---|---|---|---|---|
| 1000 | −0.15% | −0.07% | −0.02% | 1.35e4 |
| 1250 | −0.58% | −0.32% | −0.13% | 477 |
| 1500 | −2.07% | −1.34% | −0.67% | 72 |
| 1750 | −6.12% | −4.41% | −2.59% | 24 |
| 2000 | −13.9% | −10.8% | −7.2% | 11.6 |

**The two deviations bracket MESS.** Where the separation λ₂/k_uni shrinks towards 10, the two codes partly differ in how a single phenomenological rate coefficient is extracted. MarXus's k_uni lies between MESS's two directions:
- at 2000 K and 0.1 atm: dissociation −7.8%, association +7.2%;
- at 1750 K: −4.4% and +1.9%.

Which extraction is appropriate at λ₂/λ₁ ≈ 10 has to be settled from the literature, starting with Georgievskii et al., J. Phys. Chem. A 117, 12146 (2013), in `papers/ChemAct/`. That has not been done yet. No claim is made about the reason for MESS's departure.

## 7. Conclusions

1. **Olzmann's eigenvalue route reproduces MESS for C₂H₃ from 300 to 1750 K** within the model differences already known from the barrier route (exact vs semiclassical Eckart), **without any absorbing barrier**. Where the barrier route is valid (300–1000 K), both routes agree to 0.01%, at 300 K to 0.001 percentage points.
2. **k_uni (GO10's eigenvector average) is robust where λ₁ is not.** The sum-rule deviation shows where λ₁ is lost, and the messages tell the user what to do.
3. **Inverse iteration with the shifted factorization is the most accurate solver, stays the default, and covers all 40 conditions.** LAPACK DSYEVD gives the full spectrum, about 50 times faster than the in-house QL, and serves as the independent check.
4. **λ₁ itself is resolved only down to the double-precision floor** (≈ 10⁻² s⁻¹ here). Resolving it at 300 K needs higher precision in the assembly and the solver; this is planned (double-double with an MPFR reference path).
