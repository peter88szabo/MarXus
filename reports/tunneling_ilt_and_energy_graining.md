# Tunneling in k(E), ILT estimate for the Case1 entrance, energy graining

**MarXus, 2026-10-05.** Research results and decisions. The implementation status is in §6.

Sources read in this session:
- **Papers:** Miller 1979 (PDF page rendered and read); Johnston & Heicklen 1962; Georgievskii & Klippenstein 2005; Pilling & Robertson 2003.
- **Code sources:** MESS, MESMER 7.1 and TUMME 2023, read by four read-only research agents. Their file:line references were kept (§5).

---

## 1. Decisions (Peter, 2026-10-05)

1. **Tunneling**
   - Use the **exact Eckart** transmission probability only: Miller 1979 eq. 8, which is identical to Johnston & Heicklen 1962 eq. 13.
   - Do **not** use the semiclassical P = 1/(1+e^{−S}) form that MESS uses. Peter: "here we want just Eckart not general semiclassical".
   - Keep **one** formula. Peter: "make clear which tunprop formula to use, and remove the other!". Done in §6.
   - Convolve P with the TS states as in Miller 1979.
2. **B12 ILT block.** Make it for Case1 B12 with estimated recombination parameters (§3).
3. **Energy graining.** Chosen convention: "Cells + common centres" (§4).
   - 1 cm⁻¹ cells for all counting and convolutions.
   - Grains are cell sums with flux-weighted k (PR03/MESMER).
   - Boltzmann factors and the collision kernel use the grain centre on the common absolute grid, so detailed balance between wells stays exact.
4. **Double-precision limits** of the final steady state: postponed ("later").

---

## 2. Tunneling: theory used (Miller, J. Am. Chem. Soc. 101, 6810 (1979))

The equations below were checked against the rendered PDF page 6811.

**Rate coefficient and tunneling states.**
- eq. 1: k(E) = N(E)/(2πħ ρ(E)).
- eq. 6: N_QM(E) = Σ_n P(E − ε_n‡), with ε_n‡ = V₀ + Σᵢ ħωᵢ‡(nᵢ + ½). The sum runs over the TS levels, without the reaction coordinate.
- eq. 7: k_QM(E) = Σ_n P(E − ε_n‡) / (2πħ N₀′(E)).
- eq. 9: N_QM(E) = ∫_{−V₀}^{E−V₀} dE₁ P(E₁) N′(E−E₁) = ∫ dE₁ P′(E₁) N(E−E₁). This is "a convolution of the classical approximation … and the tunneling probability".
- Miller's warning: in the threshold region, do not replace the TS level sum by a classical closed form (eq. 10). Keep the discrete count (eq. 11).

**Eckart transmission, eq. 8.**

P(E₁) = sinh a · sinh b / [sinh²((a+b)/2) + cosh² c], with

- a = (4π/ħω_b)·√(E₁+V₀)·(V₀^{−½} + V₁^{−½})⁻¹
- b = (4π/ħω_b)·√(E₁+V₁)·(V₀^{−½} + V₁^{−½})⁻¹
- c = 2π·√(V₀V₁/(ħω_b)² − 1/16)

Here E₁ is the reaction-coordinate energy relative to the barrier top, and V₀ and V₁ are the barrier heights from the two sides.

- **Symmetry.** P is symmetric in V₀ ↔ V₁, so one N_QM serves both directions and detailed balance holds.
- **Lower limit.** P = 0 for E₁ ≤ −min(V₀, V₁).
- **Imaginary c.** When V₀V₁/(ħω)² < 1/16, c is imaginary and cosh² c becomes cos² |c| (Johnston & Heicklen, J. Phys. Chem. 66, 532 (1962), text after eq. 13).

**Overflow-free evaluation used in MarXus.** The formula is rewritten as

P = [cosh(a+b) − cosh(a−b)] / [cosh(a+b) + cosh 2c]

(the Johnston–Heicklen form). All terms are scaled by e^{−max(a+b, 2c)}, and expm1 is used for small a and b.

**Thermal factor.** Γ(T) = ∫_{−V₀}^{∞} dE₁ P(E₁) e^{−E₁/kT}/kT (Miller p. 6811). It is the same in both directions.

---

## 3. ILT estimate for Case1 barrier B12 (R = C₅H₉O₃ + O₂ → G2)

**The deck's B12 model** is a PST core with:
- isotropic potential V = −C₆/R⁶, C₆ = `PotentialPrefactor` = 2.4 au, n = 6;
- core SymmetryFactor 2;
- fragments C₅H₉O₃ (doublet, σ = 1) and O₂ (triplet, σ = 2);
- TS ElectronicLevels: doublet;
- asymptote GroundEnergy = 16.3 kcal/mol, which is the deck's energy reference ("the energy is measured from here").

**Capture rate for an isotropic potential.** Georgievskii & Klippenstein, J. Chem. Phys. 122, 194103 (2005), **eq. 55**:

k(T) = (8π)^{½}·((n−2)/2)^{2/n}·Γ(1−2/n)·μ^{−½}·V₀^{2/n}·T^{½−2/n}   (atomic units).

For n = 6 this is 8.55·μ^{−½}·C₆^{1/3}·T^{1/6} (their eq. 57). It agrees with my independent derivation from σ(E) = (3π/2)(2C/E)^{1/3}: k = πΓ(2/3)·√(8kT/πμ)·(2C/kT)^{1/3}.

GK2005 (text after eq. 55): for isotropic potentials, additional uncoupled modes give "a simple product … which cancels with the corresponding partition function for the reactants".

**Statistical factors.**
- Electronic: g‡/(g_A·g_B) = 2/(2·3) = 1/3.
- Symmetry: the core SymmetryFactor 2 equals σ_A·σ_B = 1·2, so the symmetry factors cancel.
- Hence k_PST = k_cap/3.

**Numbers.** μ = 117.06·31.998/149.06 = 25.13 amu (MarXus atomic-mass table).

| T (K) | k_cap (cm³ s⁻¹) | k_PST = k_cap/3 (cm³ s⁻¹) |
|---|---|---|
| 200 | 9.61e-11 | 3.20e-11 |
| 298 | 1.027e-10 | **3.42e-11** |
| 500 | 1.119e-10 | 3.73e-11 |
| 1000 | 1.256e-10 | 4.19e-11 |
| 1500 | 1.344e-10 | 4.48e-11 |

**ILT parameters** (association): **A = 3.42×10⁻¹¹ cm³ s⁻¹, n = 1/6, T_ref = 298 K, E∞ = 0.**
- The inversion is valid because ν = n + 3/2 = 5/3 > 0.
- The parameters reproduce the deck's own barrierless model at the canonical level.
- **Caveat:** isotropic PST capture is an upper bound for radical + O₂ association. Real anisotropy reduces it, and measured RO₂-forming rate coefficients are often ~10⁻¹²–10⁻¹¹ cm³ s⁻¹. **Replace with literature values if available.**

---

## 4. Energy graining: what the codes do; decision

| | State counting | Mapping to master-equation grains | Grain energy | Grain size |
|---|---|---|---|---|
| **MESS** | per species on `InterpolationEnergyStep` = 1 cm⁻¹ (model.cc:25587); classical rotor sampled, oscillators by direct count; natural cubic spline of ρ or N (not log) | point evaluation at the grid nodes, no averaging (mess.cc:6334-6350, 7001-7020) | node E_i = E_ref − i·dE | dE = EnergyStepOverTemperature·kT, recomputed per T; E_ref = highest barrier ground + ExcessEnergyOverTemperature·kT (mess_driver.cc:1851-1855); manual gives no convergence advice |
| **TUMME** | DE = 1 cm⁻¹ (default), Beyer–Swinehart, up to EMAX = 400 kcal/mol | nearest fine-grid value at the grain energy, no averaging (species_lib.py:282-287) | point E_i = E_top − i·ΔE | ΔE = 0.1·kT (default); E_top = max E₀ + 30·kT; manual: decrease ESOT at low T |
| **MESMER** | **cells** of 1 cm⁻¹ (System.cpp:194-196): classical rotors exact per cell, Beyer–Swinehart for vibrations | **grain DOS = Σ cell DOS** (MesmerTools.cpp:54-100); **grain flux = Σ cell flux**, so k_grain = Σ F/Σ ρ (Reaction.cpp:251-265) | DOS-weighted mean of the cell centres (MesmerTools.cpp:79-87) | GrainSize 100 cm⁻¹ (default); EMax = highest ZPE + 25·kT |

**Pilling & Robertson, Annu. Rev. Phys. Chem. 54, 245 (2003), p. 254:** grains are "contiguous intervals … for which mean values of energy, microcanonical rate coefficient and density of states are assigned. The densities of states are usually calculated by a direct-count algorithm, such as that developed by Beyer & Swinehart. A δ-function is then used to represent each grain and is located at the mean energy."

**Decision ("Cells + common centres").**
- All counting and convolutions on 1 cm⁻¹ cells: ρ, W‡, Eckart, ILT, fragment pairs.
- Grain ρ = average of the cell ρ over the grain (states per cm⁻¹).
- Grain W = average of the cell N over the grain, so k_g = W_g/(hρ_g) = Σ cell flux / Σ cell states (PR03 and MESMER).
- The grain sits at its centre on the common absolute grid, so all wells share the same Boltzmann factor per grain. This keeps inter-well detailed balance exact.
- **Deviation from PR03:** the "δ at the mean energy" differs only in the lowest, partially filled grain of each well. A per-well mean energy would make the inter-well couplings asymmetric by about (ΔE/kT)·ΔE·ρ′/ρ ≈ 10⁻³–10⁻⁴ and break the Cholesky path.
- **Grain width** is rounded to whole cells.

**What this removes.** MarXus's current counting on master-equation grains has the following O(ΔE) biases:
- A discrete convolution of two grain densities represents energy (i−1)·ΔE at index i: half a grain off the grain centre.
- Thresholds are rounded to grains.
- Classical-rotor grains are empty at the ground state.

With 1 cm⁻¹ cells these become 1 cm⁻¹ effects.

**Top energy.** The adapter used ModelEnergyLimit as the master-equation top (300–400 kcal/mol in the decks). MESS uses ModelEnergyLimit only as the top of state counting. Its master-equation top is E_ref = highest barrier ground + ExcessEnergyOverTemperature·kT. MarXus will follow that, with T = the highest temperature of the deck, because its network is built once.

---

## 5. Research details per code (with file:line)

### 5.1 MESS: tunneling

Paths are relative to `MESS_kinetics/Source_from_2026/MESS/src/`. The compiled files are libmess/model.cc, libmess/mess.cc and mess_driver.cc.

**Transmission** (model.cc:5665-5684 and 6430-6476): P(ε) = 1/(1+e^{−S}) with the **WKB action of the Eckart potential**:
- x = ε/ħω, d_{a,b} = V_{a,b}/ħω, F = 4π/(d_a^{−½} + d_b^{−½});
- S = F·[√(x+d_a) − √d_a + √(x+d_b) − √d_b], with each √ argument clamped at ≥ 0;
- S → 2πε/ħω at the top, the parabolic limit.

This is **not** the exact Eckart formula (it drops cosh 2π(α−β) and approximates cosh). The agent compared both for 2658.84 cm⁻¹, 19.5 / 21.3 kcal/mol:
- at the top, P_MESS = 0.500 versus exact 0.537;
- thermal κ_MESS = 0.865 × exact at 200–300 K, 0.92 at 1000 K, 0.96 at 2000 K.

MESS cites no literature for tunneling.

**WellDepth.**
- Exactly two values are required (model.cc:6369-6373) and they are sorted (6376-6378), because the formula is symmetric.
- Units: kcal/mol, 1/cm or kJ/mol. ImaginaryFrequency is compulsory.

**Convolution** (model.cc:5460-5512, `Tunnel::convolute`). On the RRHO 1 cm⁻¹ grid:

N_tun[j] = N_j·P_0 + Σ_{i=1}^{min(j,n−1)} N_{j−i}·(P_i − P_{i−1}) + [j≥n]·N_{j−n}·(1 − P_{n−1}),

a Riemann–Stieltjes sum of ∫N(E−ε)dP(ε) with ε_i = −c + iΔ and n = ceil((c + 2ħω)/Δ).
- **Cutoff c:** CutoffEnergy if given, otherwise the smaller WellDepth (6382-6384). It is reduced if −S(−c) > TunnelActionMax (default 100).
- **Range:** P is taken as 1 above +2ħω.

**Placement and k(E).**
- ε = 0 is the TS ground state (ZeroEnergy).
- The tunneling ground is real_ground − c, clamped at the highest connected-well ground (26434-26468).
- k(E) = N_tun/(2πħρ).
- One N_tun is used for both directions (mess.cc:11030-11040), so detailed balance holds.

**Other models:** Harmonic, Quartic, Read and Expansion. All use the same P = 1/(1+e^{−S}).

**Manual versus code.** The manual's "TunnelingActionCutoff, default 40" is `TunnelActionMax`, default 100, in the code.

### 5.2 MESS: energy grid

- dE = EnergyStepOverTemperature·T; the value must be in (0, 1) (mess_driver.cc:724-731).
- Grid nodes E_i = E_ref − i·dE (mess.cc:219).
- **Wells:** ρ(E_i) is point-evaluated from the spline (6334-6350). The grid ends above E_base, which is the well bottom, raised by WellCutoff or GlobalCutoff if given.
- **Barriers:** N(E_i) point-evaluated (7001-7020).
- **Matrix elements:**
  - off-diagonal −N/(2π√(ρ₁ρ₂)) (11030-11038);
  - bimolecular source N·e^{(E_ref−E)/T}/(2π) (11120-11129).
- **ModelEnergyLimit** is the top of state counting (model.cc:26482-26502); "Typically … 400 kcal/mol" (manual p. 9).
- **States above the spline range** are extrapolated with a power law (26717-26721).

### 5.3 MESMER 7.1: tunneling

Paths are relative to `Mesmer7.1-source/src`.

**Eckart** (plugins/EckartCoefficients.cpp:122-133):
- Formula: **exactly Miller's eq. 8**, with E measured from the TS ZPE. The code comment cites Miller 1979.
- Default V₀ and V₁ are *classical* barriers (ZPE removed); `me:useZPE` gives ZPE-based barriers.
- P is evaluated at the lower edge of each 1 cm⁻¹ cell, with no cell integration.
- P = 0 below max(ZPE_rct, ZPE_pdt).
- NaN is set to 0.
- **The double-precision formula overflows.** P becomes 0 for a+b ≳ 711, or at every energy for c ≳ 355 (√(V₀V₁)/ν̃ ≳ 56). MarXus's scaled form avoids this.

**Convolution** (plugins/RRKM.cpp:69-103):
- N[i] = Σ_{j≤i} ρ_TS[j]·P[i−j], computed by FFT.
- The flux bottom is the ZPE of the higher well.
- One flux is used for both directions (IsomerizationReaction.cpp:276-281), so detailed balance holds.

**Other options:** WKB (Garrett–Truhlar), user-defined coefficients, spin crossing.

### 5.4 MESMER 7.1: cells and grains

- **Sizes.** CellSize 1 cm⁻¹; GrainSize default 100 cm⁻¹ (defaults.xml:55).
- **Energy range.**
  - EMin = lowest ZPE.
  - EMax = highest ZPE + energyAboveTheTopHill·kT, default 25.
- **State counting.** Classical rotors give exact states per cell (ClassicalRotor.cpp:91-146); vibrations are added by Beyer–Swinehart.
- **Grain DOS.** The sum of the cell DOS. The first, partial grain holds cellPerGrain − offset cells (MesmerTools.cpp:54-100).
- **Grain flux.** The sum of the cell fluxes, so k_g = Σ F/Σ ρ.
- **Grain energy.** The DOS-weighted mean of the cell centres. It is used in the exponential-down kernel and in ρ_g·e^{−βE_g}.
- **Literature.** None is cited in the code. The manual §14.1 refers to Robertson 2007, Pilling & Robertson 2003, Holbrook–Pilling–Robertson 1996 and Robertson 2019.

### 5.5 TUMME 2023

**Grid.**
- ME grain ΔE = 0.1·kT (default); E_top = max E₀ + 30·kT.
- Fine grid DE = 1 cm⁻¹ for ρ (Beyer–Swinehart), up to EMAX = 400 kcal/mol.
- Grain values are the nearest fine-grid point; there is no averaging (species_lib.py:282-287).

**Tunneling.**
- **Eckart (Johnston–Heicklen):**
  - ε is measured from the higher asymptote, with the two momenta √ε and √(ε + |V₁ − V₂|).
  - This uses √|α₁α₂ − π²/4| in cosh. That is wrong when the argument is negative, where cos should be used.
- **SCT/ZCT/LCT/μOMT:** P(E) tables are read from Polyrate.
- **Convolution:** N(E_η) = DE·Σ P[k]·ρ‡[j−k] from max(E₀,R, E₀,P) (reaction_lib.py:392-478). It has the form of Miller's but Miller is not cited.

**Chemical activation.**
- No separate steady-state solver. The association source is B_η ∝ N·e^{−βE}·ΔE.
- With all CSE modes merged, k_{R→P} = Σ k̂·W⁻¹·Δk̂, which is the steady-state CA expression (FD2022 eq. 25).

**Bugs noted (TUMME, not ours).**
- Electronic levels in the DOS are overwritten instead of added (molecule_lib.py:453).
- The Eckart d term above.

---

## 6. Implementation status (MarXus)

### 6.1 Done

**Tunneling functions** (`src/tunneling/tunneling.rs`), developed test-first.

- **`eckart_transmission_probability(E₁, V₀, V₁, ħω)`.** Miller eq. 8 in the overflow-free form, with the cos branch for imaginary c. Tests:
  - agrees with the literal sinh/cosh form to 10⁻¹²;
  - symmetric in V₀ ↔ V₁;
  - P = 0 below −min(V₀, V₁) and P → 1 far above the barrier;
  - parabolic limit 1/(1+e^{−2πE₁/ħω}) near the top of a wide barrier;
  - finite for narrow, high barriers (a, b, c of several hundred);
  - imaginary-c branch.
- **`eckart_tunneling_sum_of_states(N‡, δ, V₀, V₁, ħω)`.** Miller eq. 9 as a Stieltjes sum.
  - Cells of width δ are centred at E₁ = jδ, each carrying ΔP_j = P((j+½)δ) − P((j−½)δ).
  - N_QM(kδ) = Σ_j ΔP_j·N‡((k−j)δ).
  - For a step P this gives N_QM = N‡ exactly (tested).
  - The states extend m = floor(min(V₀,V₁)/δ + ½) cells below the top.
  - N‡ must be supplied m cells beyond the highest energy wanted. That is the documented contract, so nothing is truncated silently.
- **Canonical `eckart(β, ω, V_f, V_b, dE, E_max)`.** Now returns one κ and uses only the exact transmission. Tested: κ is the same in both directions.

**Removed:** `tunprop1` and `tunprop2`, as Peter asked.

**Bug fixed in the old code.** The canonical κ used `tunprop1`, which had `|α₁ − α₂|` while E was measured from the *forward* asymptote.
- **Effect:** for the endothermic direction (V_f > V_b) the transmission was wrong. It was non-zero below the product asymptote: 1.26e-6 where it must be 0. At E = 4000 cm⁻¹ it was 0.98 instead of 0.012.
- **Who was affected:** `high_pressure_limit.rs:242` (canonical TST with Eckart).
- **Where else this formula appears:** TUMME uses the same |α₁ − α₂| form *correctly*, because it measures ε from the higher asymptote.

`tunprop2` had a related defect: sqrt of a negative number below the product asymptote gave NaN. Its κ₂ was never used.

**Parser** (`mess_input.rs`):
- `TunnelingSpecification::Eckart{imaginary_frequency_cm1, well_depths_cm1: [_;2]}`.
- Other models are recorded as `Unsupported{model}`.
- Two positive WellDepth values are required.
- Tests: values read in cm⁻¹; other models recorded; a missing WellDepth is an error.

### 6.2 Wigner fallback removed (Peter: "we do not need fallback, Eckart works in general")

`high_pressure_limit::eckart_tunneling_kappa` used to replace the Eckart κ silently by **Wigner** whenever κ ≥ 10⁴ or κ was non-finite. That guard dated from the overflow-prone old routine. It is removed. The function now returns the exact-Eckart κ, or an error if κ is not a positive finite number.

Size of the effect, for Case1 B23 (ħω = 2658.84 cm⁻¹, 19.5 / 21.3 kcal/mol). κ_Eckart was evaluated independently by a trapezoid on 1 cm⁻¹:

| T | κ_Eckart | κ_Wigner |
|---|---|---|
| 200 K | 8.30e10 | 16.2 |
| 300 K | 7.39e4 | 7.8 |
| 500 K | 23.8 | 3.4 |

At 300 K the old guard would have made the TST rate about 10⁴ times too small.

The new test `deep_tunneling_keeps_the_eckart_correction` (300 K, κ = 7.39e4 within 1%) failed with the guard and passes without it.

Wigner, Bell and Skodje–Truhlar (`wigner`, `bell`, `skodje_truhlar`, `skodje_truhlar_exact` in `tunneling.rs`) remain **optional canonical κ(T) models** that the user can choose; the TST routines take `tunneling_kappa` as an input. None of them is used as a substitute for another (Peter: "Wigner is just optional model or Skodje formulas or others as well").

### 6.3 Energy graining: implemented

**Module `src/masterequation/energy_graining.rs`** (new, 5 tests):
- `GrainGrid{cell_width_cm1, cells_per_grain}`. The grain width is rounded to whole cells.
- Grain g collects the absolute cells g·n − ⌊n/2⌋ … g·n − ⌊n/2⌋ + n − 1. It is centred at g·ΔE, exactly for odd n and half a cell higher for even n.
- `average_over_grains`: grain values are cell averages. Cells below a species' first cell count as 0, which gives the partial lowest grain. It is an error if the cell values do not reach the top grain.
- Tests:
  - grain centres and the cell ↔ grain mapping;
  - rounding of the grain width;
  - the sum over the grains equals the sum over the cells;
  - the partially filled lowest grain;
  - the coverage error.

**`Channel.threshold_grain: Option<usize>`** (classical threshold).
- `Well::lowest_threshold_grain` uses it when given, otherwise the first grain with k > 0.
- The absorbing barrier is therefore placed 10 kT below the *classical* threshold, also when tunneling makes k(E) > 0 below it.
- Validation: the threshold must lie on the grid. Tested.

**Adapter (`chemical_activation_from_mess_input.rs`) rebuilt on cells.**

*Grid*
- **Cells:** `MessNetworkSettings.cell_width_cm1`, default 1 cm⁻¹. The grain is EnergyStepOverTemperature·k_B·T_min, rounded to whole cells; 41.7 → 42 cm⁻¹ for Case1.
- **Top:** the highest barrier ZeroEnergy or bimolecular GroundEnergy + ExcessEnergyOverTemperature·k_B·T_max (the input-format convention), otherwise ModelEnergyLimit.
  - `ExcessEnergyOverTemperature` is now parsed (test).
  - For Case1 the top is 19.5 kcal/mol + 40·kT(300 K) ≈ 43 kcal/mol instead of 300 kcal/mol.

*Well densities*
- ρ is counted on the cells from the ground state to the top and averaged over the grains.
- Leading empty grains are dropped. With classical rotors the ground-state cell is empty.

*Channels*
- **Common rule:** k_g = ⟨W⟩_g/(h·⟨ρ⟩_g). The isomerization uses the same ⟨W⟩_g in both directions, so detailed balance is exact.
- **Tight TS:** N‡ on the cells.
- **With `Tunneling Eckart`:**
  - N‡ is counted m = ⌊min(WellDepth)/δ + ½⌋ cells beyond the top;
  - then `eckart_tunneling_sum_of_states` (Miller eqs. 8–9) is applied;
  - the cells below the **higher ground state of the two sides** are dropped, so tunneling only reaches energies where both sides have states. MESS clamps its cutoff the same way (model.cc:26442-26456).
- **ILT association:** fragment densities on the cells, convolved as ρ_AB = ρ_A ⊗ ρ_B·δ, then the ILT on the cells.
- **ILT dissociation:** on the well's cell density.
- **Thresholds:** `threshold_grain` is the grain of the TS (or ILT threshold) cell.
- **Unsupported tunneling models** are an error unless `ignore_tunneling` is set.

*Adapter tests (13)*
1. Shared grid and top.
2. Grain densities are cell averages (number of states conserved to 10⁻¹²).
3. Exact isomerization detailed balance, with and without tunneling.
4. Channels and classical thresholds.
5. ILT entrance at the asymptote.
6. Collision parameters and sink.
7. A PST core without ILT is refused.
8. Eckart tunneling opens the isomerization more than 10 grains below the TS, increases k above it, and keeps the threshold.
9. Tunneling stops at the higher ground state when WellDepth exceeds the deck's energy differences.
10. `ignore_tunneling` reproduces the classical rates exactly.
11. Unsupported models are refused.
12. The deck runs through the driver, with and without tunneling.
13. **C₂H₃: the final steady state reproduces canonical TST within 0.5%** (next paragraph).

### 6.4 Validation: master equation versus canonical TST (C₂H₃, 1000 K)

With thermal formation through the only channel, the final steady state is equilibrium, and k^ca must be the canonical TST rate coefficient.

Independent evaluation from partition functions (quantum harmonic vibrations from the zero-point level, classical rigid rotors, the deck's geometries, E₀ = 38.86 kcal/mol): **k∞(1000 K) = 2.4150×10⁵ s⁻¹**.

| Master equation | k^ca (final steady state) | Deviation |
|---|---|---|
| old counting on 70 cm⁻¹ grains | 2.5672×10⁵ s⁻¹ | +6.3% |
| **1 cm⁻¹ cells → 70 cm⁻¹ grains** | **2.4130×10⁵ s⁻¹** | **−0.08%** |

The intermediate steady state at 1000 K changed by about −2% in Φ_stab:

| P (Torr) | Φ_stab, old | Φ_stab, cells |
|---|---|---|
| 76 | 0.01166 | 0.01140 |
| 760 | 0.0609 | 0.0596 |
| 7600 | 0.2273 | 0.2233 |

### 6.5 Case1 with ILT entrance and Eckart tunneling

**Deck.** `examples/Case1_ZZ_allyl+O2_2025_12_02_ilt.inp`. It is a copy of Peter's deck with the B12 `InverseLaplaceTransform` block of §3; **the original deck is unchanged**.

**Run** (`cargo run --release --example chemical_activation_from_deck -- <deck>`, ≈ 4 s):
- 5 wells. Grain 42 cm⁻¹ on 1 cm⁻¹ cells. Grains per well: G2 362, G3 347, G4 416, G6 404, G13 380.
- 8 Eckart barriers, the ILT entrance, and the G4 escape (2.5×10⁷ s⁻¹).

Yields at 300 K, 760 Torr:

| | Back to R (B12) | P1 (B4P1) | P7 (B6P7) | P5 (B4P5) | Stabilized | Sink (G4) | Mass balance; residual |
|---|---|---|---|---|---|---|---|
| intermediate | 0.43734 | 3.53e-5 | 1.33e-3 | 1.57e-3 | G2 0.0667, G3 2.4e-4, G4 0.4051, G6 9.0e-6, G13 4.5e-5 | 0.0876 | 1; 3e-15 |
| final | 0.43736 | 3.53e-5 | 1.33e-3 | 1.57e-3 | 0 | 0.5597 | 1 (2e-11); 1.5e-10 |

Population fractions in the final steady state: **G6 0.998**. G6 is a dead end; its only exit, B36, lies 17.9 kcal/mol above it. Its population is large and drains slowly through G3 → G4 → escape.

This is the near-singular situation of the postponed double-precision topic. The residual (1.5e-10) is still acceptable here.

### 6.6 Data questions on the Case1 deck (for Peter)

1. **WellDepth versus ZeroEnergy.** Some WellDepth values do not equal the ZeroEnergy differences:
   - B23: 21.3 kcal/mol from G3, versus 19.5 − 1.8 = 17.7;
   - B34: 20.8 from G3, versus 19.0 − 1.8 = 17.2.
   
   The Eckart shape uses the WellDepths as given (as MESS does). The tunneling range is clamped to the higher ground state.
2. **B6P7 connects G4 → P7**, while its comment says "TS-6-P7". It may be a typo for G6.
3. **The B12 ILT parameters** are the isotropic-PST estimate (§3). Please replace them with literature values when available.

### 6.7 Still open

- **Final steady state, double precision:** postponed by Peter.
- **ILT and fragment convolutions are O(n²)** on the cells: about 10⁸ operations for Case1, < 1 s in release. An FFT would be needed for much larger grids.
- **Excited electronic levels** are still ignored (ground-level degeneracy only).
- **`inertia::get_brot` prints** "Iterative diagonalization is done …" for every geometry.

## Update (2026-10-05): the MESS Eckart model as an option

**Peter's request:** "build a new tunneling model call it mess_eckart_tunneling, mimic it" (MESS is under the Apache License 2.0).

**What MESS calls Eckart tunneling.** It is a semiclassical model, not the exact Eckart transmission of Miller (1979), eq. 8. Source: MESS `src/libmess/model.cc`, `Model::Tunnel` and `Model::EckartTunnel`.
- **Transmission:** P(E) = 1/(1 + e^{−S(E)}), with S(E) = 4π/(d₀^{−½} + d₁^{−½})·Σ_w[√(max(E/ω + d_w, 0)) − √d_w], d_w = V_w/ω. Near the top this is the parabolic barrier.
- **Clamps:** P = 1 for S > 100, P = 0 for S < −100.
- **Cutoff:** E_c is the smaller well depth, lowered by bisection if −S(−E_c) > 100.
- **Convolution:** of the number of states up to 2ω above the top.
- **Canonical factor:** a rectangle sum on 0.01 kT up to 10 kT.

**In MarXus:**
- `src/tunneling/mess_eckart_tunneling.rs`, with the same interface as `eckart_tunneling_sum_of_states`.
- Selection through `MessNetworkSettings::eckart_tunneling` (`EckartTunnelingModel::{Exact, Mess}`, default Exact) or `--tunneling mess-eckart`.
- `examples/eckart_kappa_from_deck.rs` prints both factors.

**Tests:**
- the parabolic limit, P(0) = ½, the clamps;
- the cutoff rule;
- the Stieltjes property of the convolution;
- the reproduction of the factors printed in the MESS log of `validation/ZZAllyl+O2_Gamma_Case2`;
- the adapter scaling of k∞ by κ_MESS/κ_exact.

**Agreement with the MESS log.** ≤ 5·10⁻⁵ for 15 of 18 values. For deep tunneling below 300 K (B23, B24, B34), +0.06 … +0.4% at 300 K and up to +2.5% at 200 K: MESS's ground-state bookkeeping lowers the cutoff by about 60 cm⁻¹ there.

**The exact model is larger.** The exact Eckart factor exceeds the MESS model by 17–23% for deep H-transfer tunneling (κ ≈ 10²–10⁵) and by 2–6% otherwise. This explains the systematic MarXus/MESS offsets of the validations (C₂H₃: +3–6% at 300–500 K; ZZ-allyl + O₂ Case 2: +17–22% for B23, B24, B34).

**The exact Eckart remains the MarXus default** (Peter's decision on exact Eckart tunneling).

