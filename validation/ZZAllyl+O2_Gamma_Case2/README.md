# ZZ-allyl + O₂, Gamma Case 2: MarXus reproduction of the reference MESS run

**MarXus, 2026-10-05.** Peter's request: reproduce his MESS run of this multiwell network. The key quantity is the yield of P5 (IEPOX + OH), a few percent of the total. Most of the reaction ends in the escape channel after stabilization of Gamma-4.

## 1. Files

| file | content |
|---|---|
| `Gamma-Case2_..._12.7kcal.inp/.log/.out` | Peter's MESS input, log and output (2025-09-29), **untouched** |
| `marxus_input/case2_tstlevel_E.inp` | the same deck with `TSTLevel E` added to the two phase-space-theory cores (Section 3.1); otherwise identical (checked with `diff`) |
| `run_marxus.sh` | the MarXus runs (at most 4 cores) |
| `marxus_output/case2_tstlevel_E_steady_states.out` | intermediate (absorbing barrier) and final steady state |
| `marxus_output/case2_tstlevel_E_eigenvalue.out` | eigenvalue analysis: k_uni and k∞ of every channel |
| `marxus_output/case2_default_EJ_steady_states.out` | the original deck as is (PST cores at the EJ level), for comparison |
| `compare_with_mess.py` | all comparisons; writes the CSV tables and `plots/` |
| `capture_comparison.csv`, `high_pressure_comparison.csv`, `net_yields_comparison.csv`, `apparent_rates_comparison.csv` | the compared numbers |
| `plots/pes.png` | the network |
| `plots/p5_share.png` | P5 share of the net reaction, MESS vs MarXus, and the deviations |
| `plots/high_pressure_deviation.png` | k∞ of every channel, MarXus/MESS − 1 |
| `plots/apparent_rates_760torr.png` | k(R → X) at 760 Torr |

Reproduce:

```
bash run_marxus.sh
source ~/.venvs/science/bin/activate && python3 compare_with_mess.py
```

## 2. The network (`plots/pes.png`)

Energies in kcal/mol relative to R = ZZ-allyl + O₂, from the MESS output header.

**Entrance and wells.**
- R → G2 (−16.3) through B12 (0.0), barrierless, with a phase-space-theory core (−C₆/R⁶, C₆ = 2.4 au).
- Wells G3 (−18.7), G4 (−19.4) and G6 (−19.7) are connected by B23 (0.0), B24 (−4.6), B34 (+0.4) and B36 (−1.4). All isomerizations carry Eckart tunneling; B23, B24 and B34 are H-transfers with imaginary frequencies of 2200–2700 cm⁻¹.

**Products and sink.**
- G4 → P1 (HPALD + HO₂, −5.4) through B4P1 (+1.6).
- G4 → P5 (IEPOX + OH, −27.4) through B4P5 (−3.6).
- G6 → P7 (+ HO₂, −14.9) through B6P7 (−4.7), phase-space theory, C₆ = 8.5 au.
- G4 escape: pseudo-first-order sink 2.5·10⁷ s⁻¹.

**Conditions.** 270–330 K, 500 / 600 / 760 Torr, exponential down (200 cm⁻¹ (T/300)^0.85), Lennard-Jones collisions.

## 3. What was needed in MarXus

### 3.1 Phase-space-theory barriers in the chemical-activation adapter

The PST core already existed (`src/barrierless/phasespace/`, used by `microcanonical_builder::PhaseSpaceTheoryRRHO`). What was missing:
- **The adapter wiring.** `chemical_activation_from_mess_input.rs` refused PST barriers without an ILT block.
- **The `TSTLevel` keyword.**

**Comparison with the MESS source** (`model.cc`, `Model::PhaseSpaceTheory`). The MarXus module was MESS's **E level**, term by term:
- the exponent (r+2)/2 − 2/n;
- 1/σ;
- the rotor factors 1, 16/15, 1/3, 32/105, π/12;
- d^{2/n}(1+1/d)^{(r+2)/2} (Georgievskii, Klippenstein, J. Chem. Phys. 122, 194103 (2005), eq. 58);
- V₀^{2/n}, μ and 1/√(ΠBᵢ).

The current MESS source defaults to the **EJ** level, which the module did not have.

**Now implemented, as in MESS:**
- `TSTLevel` T, E and EJ, default EJ; J=0 refused.
- The linear-fragment criterion I_min/I_mid < 10⁻⁵, with B from the middle moment.
- n > 2.

**Tests** (`phase_space_theory.rs`):
- the EJ level reproduces GK05 **eq. 55** to 10⁻¹⁰ for all fragment combinations at n = 4 and 6, and eq. 57's coefficient 8.55;
- the T level is the canonical variational capture rate;
- the E level is unchanged;
- the parser reads `TSTLevel`;
- adapter test `a_phase_space_barrier_without_an_ilt_block_uses_the_phase_space_core`: the entrance k∞ equals eq. 55 within 0.5%.

**Which level the reference run used.** MESS's log says **"TST level: E"** for B12 and B6P7: the 2025 MESS version defaulted to E. With the deck as is, MarXus's EJ capture rate is 10.6% below MESS (`capture_comparison.csv`). With `TSTLevel E` it is **+0.4 … +0.5%**. All results below use `marxus_input/case2_tstlevel_E.inp`.

### 3.2 Nothing else

Everything else in the deck was already supported:
- wells and barriers with Eckart tunneling (WellDepth);
- the escape sink;
- exponential down, Lennard-Jones;
- the 1 cm⁻¹ cells.

The MESS-only output options (HotEnergies, TimeEvolution, PED, eigenvector output) are ignored by the reader.

## 4. Results

### 4.1 High-pressure rate coefficients of every channel (`high_pressure_comparison.csv`, `plots/high_pressure_deviation.png`)

| channel | MarXus/MESS − 1, 270–330 K |
|---|---|
| B12 (G2 → R), PST | +0.3% |
| B6P7 (G6 → P7), PST | +0.7% |
| B36 | +1.5 … +2.3% |
| B4P1 | +3.0 … +3.5% |
| B4P5 | +5.3 … +6.1% |
| B34 | +16.8 … +17.9% |
| B23 | +20.3 … +21.5% |
| B24 | +21.4 … +22.7% |

**The differences are entirely the Eckart tunneling factor.** The canonical κ(300 K) of MarXus's exact Eckart transmission (`tunneling::eckart`; Miller, J. Am. Chem. Soc. 101, 6810 (1979), eq. 8) compares with the factors in the MESS log ("tunneling partition function correction factors") as follows:

| barrier | imaginary frequency (cm⁻¹) | κ MESS | κ MarXus | ratio | k∞ deviation at 300 K |
|---|---|---|---|---|---|
| B23 | 2742.1 | 2.357·10⁴ | 2.820·10⁴ | 1.197 | +21.0% |
| B24 | 2191.7 | 304.9 | 367.9 | 1.207 | +21.8% |
| B34 | 2666.8 | 4.567·10⁴ | 5.331·10⁴ | 1.167 | +16.9% |
| B36 | 602.6 | 1.436 | 1.460 | 1.017 | +1.7% |
| B4P1 | 586.3 | 1.390 | 1.423 | 1.024 | +3.2% |
| B4P5 | 988.9 | 2.962 | 3.087 | 1.042 | +5.7% |

After κ, every channel agrees within 0.3–0.7%. This includes the PST channels, which have no tunneling. Partition functions, symmetry factors and energies therefore agree, and the remaining difference is how the Eckart tunneling factor is evaluated.

**Interpretation.** MarXus uses the exact Eckart transmission probability, as decided (`reports/tunneling_ilt_and_energy_graining.md`). MESS's κ is lower: by 2–4% for the moderate barriers, and by 17–21% in the deep-tunneling H-transfers, where κ ≈ 10²–10⁵. This is the same direction as in the C₂H₃ benchmark (`validation/c2h3_mess_example/`). How MESS evaluates its Eckart factor was not examined here.

### 4.2 The P5 yield (`net_yields_comparison.csv`, `plots/p5_share.png`)

**How the shares are computed.**
- *Normalization.* All shares are fractions of the net reaction, P1 + P5 + P7 + escape = 1. MarXus normalizes to the capture rate k∞, which includes the prompt redissociation of hot G2 to R (65–89% of k∞); MESS's phenomenological k(R → X) do not.
- *MarXus.* Final steady state (Olzmann, PCCP 4, 3614 (2002)) with the G4 escape as its physical sink. This gives the long-time yields directly.
- *MESS.* The long-time fate computed from its (T, p) rate tables: R forms each well, and each well then ends in a product or the escape (absorption probabilities of the well network).

| T (K) | P5 share MESS, 500 / 600 / 760 Torr | P5 share MarXus | deviation |
|---|---|---|---|
| 270 | 2.22 / 1.76 / 1.29% | 2.48 / 1.97 / 1.45% | +11.5 … +12.0% |
| 300 | 4.03 / 3.27 / 2.47% | 4.40 / 3.57 / 2.70% | +9.1 … +9.5% |
| 330 | 7.17 / 5.99 / 4.72% | 7.66 / 6.41 / 5.04% | +6.9% |

**Other channels.**
- **Escape** agrees within 0.2–0.5%: 92–99% of the net reaction.
- **P1** is +7 … +11%, and **P7** −9 … −14%. Both are tiny, 10⁻⁶–10⁻⁵ of the net reaction.

**The P5 offset follows the tunneling difference.** It falls with T exactly as the κ ratios do: B4P5 itself is +4% at 300 K, and the isomerization network that feeds G4 is +17–22%.

### 4.3 Apparent rate coefficients (`apparent_rates_comparison.csv`, `plots/apparent_rates_760torr.png`)

MarXus intermediate steady state (absorbing barrier 10 kT below the lowest threshold of each well), k(R → X) = k∞Φ_X, against the MESS R row. At 270 K and 500 Torr:

| channel | MESS (cm³/s) | MarXus (cm³/s) |
|---|---|---|
| R → G2 | 8.01e-12 | 7.74e-12 |
| R → G4 | 3.05e-12 | 2.93e-12 |
| R → P5 | 2.48e-13 | 2.87e-13 |

The MESS escape entry R → G4-escape is **negative**, −1.3·10⁻¹³ cm³/s, and several other MESS entries are negative or tiny (e.g. G2 → G4-escape −2.9·10⁻¹⁰ s⁻¹). The phenomenological splitting of MESS's chemically significant eigenvalues is therefore not resolved for the escape channel at these conditions. The long-time shares of Section 4.2 do not depend on this splitting.

## 5. Conclusions

1. **MarXus reproduces the reference MESS run of this four-well network with two phase-space-theory channels and a physical sink.**
   - It needs no change to the deck except stating the TST level that the reference run actually used (E).
   - The capture rate agrees within 0.5%, and every channel's k∞ within 0.3–0.7% once the tunneling factor is accounted for.
2. **The P5 (IEPOX + OH) share** is 1.3–7.2% (MESS) and 1.4–7.7% (MarXus), rising with T and falling with p. MarXus is 7–12% higher, entirely through its exact Eckart tunneling factors, which are larger by 2–21% than MESS's.
3. **The escape after stabilization of G4** is 92–99% in both codes, agreeing within 0.5%.

## 6. Open points

1. **MESS's Eckart factor.** Whether to examine how MESS evaluates its Eckart factor. MarXus keeps the exact Eckart transmission by decision.
2. **Default PST level.** Whether MarXus's default PST level should be E, as the 2025 MESS version that produced most existing decks, or EJ, as the current MESS source. It is EJ now; the deck copy states E explicitly.
3. **Atomic masses.** MarXus's atomic-mass table mixes isotopic (C 12.0) and standard (H 1.00784, O 15.999) masses; MESS uses isotopic masses. The effect here is about 0.03%.
