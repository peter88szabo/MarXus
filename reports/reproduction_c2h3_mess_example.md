# Reproduction of a reference example: H + C₂H₂ ⇌ C₂H₃ (MESS example set)

**MarXus, 2026-10-05.** Everything here is uncommitted.

> **Validation folder:** `validation/c2h3_mess_example/` holds the input decks, the stored MESS results, the MarXus results, `run_marxus.sh`, `plot_comparison.py`, the plots (PES, fall-off, deviations, high-pressure limits, 1000 K summary) and its own README with all tables. The reverse direction W1→P1 is also compared there.

## 1. Choice of the example

Peter asked to "find a good example to reproduce from Mesmer or Mess example files". Candidates:

| Example | Location | Content | Reference output | Suitability |
|---|---|---|---|---|
| **c2h3** (chosen) | `MESS_kinetics/Examples_From_Argon/examples/c2h3/` | C₂H₃ ⇌ C₂H₂ + H, one well, tight TS; three variants (`c2h3_tight_short_notunneling`, `c2h3_tight_short` with Eckart, `c2h3_tight` with Eckart at 8 T × 5 p); paper `miller2004.pdf` | `.out` files | **ideal first target**: input format MarXus reads, physics MarXus has (RRHO, tight TS, Eckart, exponential down, LJ) |
| AcetylO2 | `Mesmer7.1-source/examples/AcetylO2`, `examples/Tunnelling`, `MesmerQA/Acetyl O2 association` | CH₃CO + O₂ → CH₃C(O)OO* → wells/products; ILT association, tight TSs, Eckart variants | `baselines/Linux64/*.test` | **next target** (multiwell, ILT, Eckart); the XML input must be translated into a deck |
| c2h5o2, c2h3o2, hco, ho2, c2h6 | MESS examples from Argon | multiwell/barrierless (VRC-TST "rotd" files) | `.out` files | need variational/ROTD barrierless models that MarXus lacks |

## 2. What was added to run the example (test-driven)

- **Parser**:
  - `Atom` fragments (`Mass[amu]`, `ElectronicLevels`), with `Atom` registered as a block opener;
  - Bimolecular `GroundEnergy` in any energy unit (only `[kcal/mol]` was accepted before);
  - **bug fix:** block keywords are now matched exactly. `WellCutoff 10` used to be read as a Well block named "10" that swallowed the whole model (`line.starts_with("Well")`). The same matching was applied to Barrier, Bimolecular and Fragment.
- **High-pressure rate coefficient of the entrance** (`EntranceHighPressureRate` in `chemical_activation_from_mess_input.rs`):

  k∞(T) = Σ_E W(E)·e^{−(E−E_AB)/kT}·δ / (h·C′(μ)·(kT)^{3/2}·Q_A·Q_B)

  It uses the same 1 cm⁻¹ cell numbers of states as k(E), tunneling included. Tests:
  - the ILT entrance gives back the ILT input A within 1%;
  - the C₂H₃ tight TS gives canonical TST from partition functions within 0.5%.
- **Bimolecular rate coefficients** k(R → X) = k∞·Φ_X in the example output. For stabilization this is "the rate into the absorbing barrier" (Pilling & Robertson 2003, eq. 44).
- **Example program:** an optional second argument names the reactant (`… <deck> <reactant>`). The reference decks have no `Reactant` line.

## 3. Results

**Settings.** MarXus used its default graining: 1 cm⁻¹ cells. The grain is EnergyStepOverTemperature (0.1)·k_B·T_min, i.e. 70 cm⁻¹ for the 1000 K decks and 21 cm⁻¹ for the 300–2000 K deck. MESS used 0.1·kT per temperature.

The comparison was made with the intermediate steady state, with the absorbing barrier 10 kT below the classical threshold.

### 3.1 1000 K, 1 atm (`c2h3_tight_short_notunneling`, `c2h3_tight_short`)

| | MESS | MarXus | Δ |
|---|---|---|---|
| **no tunneling:** k∞(P1→W1), cm³ s⁻¹ | 3.88e-11 | 3.876e-11 | −0.1% |
| k(P1→W1, 1 atm) | 2.34e-12 | 2.310e-12 | −1.3% |
| k∞(W1→P1), s⁻¹ | 2.41e5 | 2.413e5 (final steady state) | +0.1% |
| **Eckart:** k∞(P1→W1) | 4.0860e-11 | 4.145e-11 | +1.4% |
| k(P1→W1, 1 atm) | 2.5705e-12 | 2.571e-12 | 0.0% |
| k∞(W1→P1) | 2.54513e5 | 2.580e5 | +1.4% |

### 3.2 Full falloff, `c2h3_tight` (Eckart, k(P1 → W1) in cm³ s⁻¹; MarXus run ≈ 19 s)

| T (K) | p (atm) | MESS | MarXus | Δ |
|---|---|---|---|---|
| 300 | 0.1 | 1.4421e-13 | 1.5261e-13 | +5.8% |
| 300 | 0.3 | 1.8112e-13 | 1.9113e-13 | +5.5% |
| 300 | 1 | 2.1202e-13 | 2.2310e-13 | +5.2% |
| 300 | 3 | 2.2930e-13 | 2.4083e-13 | +5.0% |
| 300 | 10 | 2.3870e-13 | 2.5040e-13 | +4.9% |
| 500 | 0.1 | 6.9740e-13 | 7.2366e-13 | +3.8% |
| 500 | 0.3 | 1.1408e-12 | 1.1826e-12 | +3.7% |
| 500 | 1 | 1.7523e-12 | 1.8140e-12 | +3.5% |
| 500 | 3 | 2.3395e-12 | 2.4184e-12 | +3.4% |
| 500 | 10 | 2.8901e-12 | 2.9832e-12 | +3.2% |
| 750 | 0.1 | 7.9550e-13 | 8.0714e-13 | +1.5% |
| 750 | 0.3 | 1.5805e-12 | 1.6055e-12 | +1.6% |
| 750 | 1 | 3.0644e-12 | 3.1164e-12 | +1.7% |
| 750 | 3 | 5.1099e-12 | 5.2010e-12 | +1.8% |
| 750 | 10 | 8.0055e-12 | 8.1542e-12 | +1.9% |
| 1000 | 0.1 | 5.1285e-13 | 5.0987e-13 | −0.6% |
| 1000 | 0.3 | 1.1451e-12 | 1.1415e-12 | −0.3% |
| 1000 | 1 | 2.5737e-12 | 2.5735e-12 | −0.0% |
| 1000 | 3 | 4.9976e-12 | 5.0107e-12 | +0.3% |
| 1000 | 10 | 9.3954e-12 | 9.4462e-12 | +0.5% |
| 1250 | 0.1 | 2.7659e-13 | 2.6676e-13 | −3.6% |
| 1250 | 0.3 | 6.6713e-13 | 6.4696e-13 | −3.0% |
| 1250 | 1 | 1.6573e-12 | 1.6173e-12 | −2.4% |
| 1250 | 3 | 3.5798e-12 | 3.5140e-12 | −1.8% |
| 1250 | 10 | 7.6924e-12 | 7.5985e-12 | −1.2% |
| 1500 | 0.1 | 1.4204e-13 | 1.2299e-13 | −13.4% |
| 1500 | 0.3 | 3.6165e-13 | 3.1837e-13 | −12.0% |
| 1500 | 1 | 9.6459e-13 | 8.6600e-13 | −10.2% |
| 1500 | 3 | 2.2504e-12 | 2.0588e-12 | −8.5% |
| 1500 | 10 | 5.3391e-12 | 4.9869e-12 | −6.6% |
| 1750 | 0.1 | 7.4261e-14 | 4.2763e-14 | −42.4% |
| 1750 | 0.3 | 1.9658e-13 | 1.1847e-13 | −39.7% |
| 1750 | 1 | 5.5199e-13 | 3.5240e-13 | −36.2% |
| 1750 | 3 | 1.3623e-12 | 9.2275e-13 | −32.3% |
| 1750 | 10 | 3.4783e-12 | 2.5279e-12 | −27.3% |
| 2000 | all | 4.1e-14 … 2.3e-12 | — (now an error, §4.2) | |

## 4. Interpretation

### 4.1 Agreement range

At **750–1250 K** MarXus and MESS agree within **−3.6% … +1.9%**, and within ±0.6% at 1000 K. These are two independent codes: different graining, different solvers (steady state versus eigenvalue/CSE), and different tunneling forms.

### 4.2 Low temperature (300–500 K): +3% … +6% from the tunneling model

The deviation is already present in k∞ (300 K: MESS 2.443e-13 versus MarXus 2.560e-13, +4.8%). The falloff ratio agrees within 1%:
- MESS: 1.442e-13 / 2.443e-13 = 0.590;
- MarXus: 1.526e-13 / 2.560e-13 = 0.596.

MarXus uses the **exact Eckart** transmission, as Peter decided. MESS uses the semiclassical 1/(1+e^{−S}), which is known to give a lower κ at low T (see `tunneling_ilt_and_energy_graining.md` §5.1). The difference is therefore the tunneling model, as intended.

### 4.3 High temperature (≥ 1500 K): the absorbing-barrier picture breaks down

The C₂H₃ well is 13 591 cm⁻¹ deep below the TS. With the barrier 10 kT below the threshold, the barrier lies above the well bottom by:

| T | well depth | barrier above the well bottom |
|---|---|---|
| 1000 K | 19.6 kT | 9.6 kT |
| 1250 K | 15.6 kT | 5.6 kT |
| 1500 K | 13.0 kT | 3.0 kT |
| 1750 K | 11.2 kT | 1.2 kT |
| 2000 K | 9.8 kT | below the bottom |

The thermal distribution of C₂H₃ has ⟨E⟩ ≈ 4 kT (2920 cm⁻¹ at 1000 K, final steady state). From 1500 K on, the barrier therefore cuts through the thermal distribution: thermalized molecules above the barrier are counted as not stabilized, and k(R → W) is underestimated.

MESS obtains these rate coefficients from the eigenvalue (CSE) analysis, which needs no barrier. The steady-state absorbing-barrier picture assumes that the well is deep compared with its thermal distribution (PR03: "usually placed about 10 k_BT below the reaction threshold"; CD07 p. 125).

**Fix made.** At 2000 K the barrier would lie below the well bottom. The old code silently clamped it to the bottom and returned zero stabilization. It is now an **error** that explains this (test `an_absorbing_barrier_below_the_well_bottom_is_an_error`).

The same silent clamp had made 5 existing tests pass for the wrong reason:
- driver tests run the two-well test network at 400 K;
- the zero-pressure test ran at 300 K with a 2000 cm⁻¹ threshold.

They now run at 250/300 K and 250 K, where the barrier lies inside the well.

**Open (decision for Peter).** How should MarXus report rate coefficients where the well is not deep enough compared with 10 kT plus its thermal width?
- (a) Report the thermal fraction of each well below the absorbing barrier as a validity diagnostic, and refuse below a threshold.
- (b) Let the user lower the barrier distance; the result then depends on it.
- (c) Add an eigenvalue/CSE route for such conditions.

### 4.4 Final steady state at 300 K

Final-steady-state rows fail at 300 K. The 34 kcal/mol well has no sink, so J is singular in double precision (Cholesky pivot −2.5e9). The driver refuses this with an explanation. This is the double-precision topic Peter postponed.

## 5. How to rerun

```
cargo run --release --example chemical_activation_from_deck -- \
  /home/peter/Dropbox/Research_Leuven/MESS_kinetics/Examples_From_Argon/examples/c2h3/c2h3_tight.inp P1
```

The bimolecular rate-coefficient table follows the intermediate-steady-state results.

## 6. Test status

115 library and 18 binary tests pass. New in this step:
- `atom_fragments_and_ground_energy_in_wavenumbers_are_read`
- `keywords_that_begin_like_block_names_are_not_blocks`
- `entrance_high_pressure_rate_reproduces_the_ilt_input`
- `entrance_high_pressure_rate_of_a_tight_transition_state_is_transition_state_theory`
- `an_absorbing_barrier_below_the_well_bottom_is_an_error`
- `deep_tunneling_keeps_the_eckart_correction` (Wigner fallback removed)
