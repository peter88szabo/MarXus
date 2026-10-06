# Low-energy treatment of the exponential-down kernel: from a temperature step to the MESMER reservoir state

**Date:** 2026-10-06. **Status:**
- The temperature step is found and confirmed (Sections 1–5).
- The SSUMES-type reduction is replaced, by Peter's decisions (Sections 7, 11): first by MESS's truncation (Sections 8–10), then by **MESMER's reservoir state**, which is the current rule (Sections 11–13).

**Summary of why** (each step is documented below):
1. **SSUMES-type reduction factors.** The integer window n_ref = ⌊1.5⟨ΔE_down⟩/ΔE⌋ + 1 jumps with T. That gave a step of +3% in R → G4 of ZZ-allyl + O₂ at 304.7 K, in every method.
2. **"As MESS or MESMER".** Both use the plain normalization of eq. 4.16 from the top.
   - It fails in the lowest grains of every well of both validation systems, so option (b), reduction only where it fails, still gave the step.
   - MESS's answer to that failure, truncation of the well (kernel mode `down`), removed the step in Case 2.
3. **Truncation removes states.** In a small molecule those states are thermally populated: C₂H₃ loses 2 grains at 300 K, which raised k_uni by 2.5% and broke the detailed balance of CSE's own pair by 2.4%. MESS's own reference runs never truncate (default kernel mode).
4. **MESMER's reservoir state** (manual, Sec. 14.2.1). The failing grains are kept as one thermalized state, with their full Boltzmann weight: no step, and no lost population.

**The reference document for the current rule** (physics, equations, differences from MESMER, validity, tests, results) is `reports/low_energy_reservoir_state.md`. This report records how the rule was reached.

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

## 6. Options considered

The low-energy rule is a physics choice. The options found in the literature and code already reviewed:

1. **Keep the current rule** (SSUMES V2). Document the step, and choose temperature grids accordingly.
2. **Apply a reduction only where the plain back substitution fails.** Start without any reduction, i.e. plain Robertson eq. 4.16. Where no failure occurs there is then no step. Where the reduction does switch on, the switch is still discontinuous.
3. **The UNIMOL rule used in the older multiwell kernel path** (`matrix_physics_assembly.rs`; `reports/collision_kernel_detailed_balance_and_normalization.md`, Section 5.2). Grains below E₀/2 share the normalization coefficient of the grain above them. Its switch point is fixed by E₀, not by T, so it gives no T step.
4. **A different transition model at low energies** (Robertson 2019, p. 294, eqs. 4.60–4.61).

A continuous variant of rule 1, e.g. an interpolated n_ref, would be a MarXus invention without a paper. It is therefore not proposed.

**Diagnostic that can be added independently of the decision.** Report `low_cut`, n_ref and m per well and condition in the output (RUN SETTINGS or the network block), so that a switch is visible.

## 7. Decision (Peter, 2026-10-06)

"Use the best rule as MESS or MESMER, or if there is nothing then do (b)."

**What the two codes do** (sources read):
- **MESMER 7.1** (`src/TMatrix.h`, `normalizeProbabilityMatrix`; `src/gWellProperties.h`, `collisionOperator`). Plain back substitution from the top grain (column sums 1), i.e. Robertson eq. 4.16, with no low-energy reduction and no check: where it fails, the coefficients come out negative. The Gaussian kernel plug-in mentions these "negative normalization coefficients".
- **MESS 2026** (`src/libmess/mess.cc`, `MasterEquation::Well::_set_kernel`):
  - **Kernel mode `down`.** The same plain normalization from the top. Where it fails ("cannot satisfy the constant collision frequency"), MESS truncates the well: the failing grain and all grains below are removed, and the kernel is built again (do-while loop). With `notruncation`, the failing grain's deactivating transitions are dropped instead.
  - **Default mode.** A per-pair common factor, not normalized; rejected earlier for MarXus (`collision_kernel_detailed_balance_and_normalization.md`, Section 4).
- **Neither code has a low-energy reduction like SSUMES.**

**Option (b) alone is not enough.**
- **Implemented first:** plain eq. 4.16 wherever it holds, the reduction only as a fallback.
- **Result:** in Case 2 and C₂H₃ the plain normalization fails at every temperature, in the lowest 3–9 grains of every well (Case 2: 114–342 cm⁻¹ above the bottom of about 430 grains). The fallback was therefore always active, and it reduced the transitions of all grains below 82–104 (3100–4000 cm⁻¹), with its window still jumping at 304.7 K.

**Peter's choice:** MESS's truncation.

## 8. Implementation: MESS truncation (superseded by Section 12)

- **`collision_kernels.rs`, `exponential_down_kernel`.**
  - Eq. 4.16 is solved from the top. At the first grain j where the activating probabilities alone reach 1, grains 0..=j are removed, and eq. 4.16 is solved again from the top on the remaining grid. The grains just above the new bottom lose deactivating targets, so their normalization changes. This repeats until it holds everywhere (MESS's do-while loop).
  - `CollisionKernel::truncated_grains` reports the number removed. The SSUMES-type reduction (n_ref, redfac, m) is removed.
- **`chemical_activation_operator.rs`.**
  - **States.** Truncated grains are not states (`WellCollisionData::truncated_grains`).
  - **Absorbing barrier.** A barrier at or below the truncated grains is an error with an explanation: nothing could be stabilized.
  - **Isomerization.** Isomerization at energies whose target grain is truncated is left out in both directions. That is what MESS does: an inner barrier spans `min(well(i1).size(), well(i2).size())` grains (`mess.cc`, around line 2737). In Case 2 this concerns the deep Eckart tunneling of B34/B36 into the lowest grains of G3.
  - **Reporting.** `low_energy_truncations(network, T, model)` lists the truncated wells per temperature.
- **`chemical_activation_steady_state.rs`, `project_source`.** Source weight on truncated grains is left out, and F is normalized over the existing grains. Before this change, the weight was counted as stabilized even in the final steady state, which a test exposed.
- **Output.** RUN SETTINGS, "Collisions:", lists per temperature the wells and the number of truncated grains (and their energy above the bottom).
- **Case 2, 1 K scan.** Truncated: 4–12 grains (152–456 cm⁻¹), for example G2 9 grains at 300 K and 10 at 308 K.
- **Tests:**
  - `exponential_down_truncates_the_well_below_the_grain_where_eq_4_16_fails`, against an independent implementation;
  - `exponential_down_uses_plain_back_substitution_wherever_it_succeeds`;
  - `exponential_down_kernel_is_continuous_in_the_mean_energy_transfer`, red with the old rule: P(0|0) jumped from 0.460 to 0.493;
  - `low_energy_truncations_list_only_the_wells_where_eq_4_16_fails`;
  - `truncated_grains_are_not_states_and_no_flux_reaches_them`;
  - `an_absorbing_barrier_within_the_truncated_grains_is_an_error`;
  - `isomerization_into_truncated_grains_is_omitted_in_both_directions`;
  - `source_on_truncated_grains_is_left_out_of_the_normalization`.
- **Fixtures and references adapted, with the reasons in comments:**
  - the fixtures `two_well_network` (well A: ρ = (1 + 0.02i)⁸) and `single_well_two_channels` (same ρ), so that their absorbing barriers lie above the truncated grains;
  - the high-pressure λ₁ reference: the Boltzmann average over the existing grains;
  - the equilibrium-distribution reference;
  - the temperatures of the two double-precision tests (λ₁ ≤ 0 at 120 K; no plain Cholesky factor at 130 K).

## 9. Result: ZZ-allyl + O₂ Case 2 re-run with the truncation (MESS Eckart model, CSE; superseded by Section 13)

**The step is gone.** The second differences of ln k(T) at 760 Torr follow MESS's smooth trend:
- R → G4: −0.0105, −0.0099, −0.0097, −0.0093, −0.0089 at 280–320 K (MESS: −0.0107 … −0.0090). Before, the values were +0.0216 and −0.0404 at 300 and 310 K.
- R → G2, R → G3 and R → P5 behave the same way.

**CSE against MESS after the change**, all 21 conditions (`cse_comparison.csv`; before → after):

| entry | before (SSUMES-type reduction) | after (MESS truncation) |
|---|---|---|
| R → P5 | −2.7 … −3.5% | −3.0 … −3.6% |
| R → G2 | +2.1 … +3.3% | +3.0 … +5.2% |
| R → G3 | +0.9 … +3.6% | −1.2 … −0.8% |
| R → G4 | +1.4 … +6.1% (step at 304.7 K) | −1.0 … −0.3% |
| R → escape | (rounding level) | −4.8 … −5.6% at 760 Torr |

**Truncated grains** (RUN SETTINGS, 38 cm⁻¹ grains):

| T (K) | G2 | G3 | G4 | G6 |
|---|---|---|---|---|
| 270 | 6 | 7 | 6 | 4 |
| 300 | 9 | 9 | 8 | 8 |
| 330 | 12 | 11 | 11 | 8 |

That is 152–456 cm⁻¹ above the well bottoms.

**Residual.** The thermal losses of the wells still wobble slightly, because the truncation changes by one grain at some temperatures: about ±0.4% for G6 (4, 5, 5, 8 grains at 270–300 K) and ±0.15% for G2. MESS's reference run used its default kernel mode, which does not truncate, and is smooth.

## 10. Second finding: truncation removes thermally populated grains (C₂H₃)

**C₂H₃ re-run with the MESS truncation** (`validation/c2h3_mess_example/`, `method_comparison.csv`, 2026-10-06):
- **Truncated grains** (21 cm⁻¹ grains): 2 at 300 K (42 cm⁻¹), 3 at 500 K, 5 at 1000 K, 45–53 at 1500–2000 K (945–1113 cm⁻¹).
- **Dissociation k_uni** (equal in SteadyStateOlzmann and CSE, the one-well identity): +7.3 … +8.3% from MESS at 300 K. Before it was +4.7 … +5.7%.
- **Association:** CSE (G13 eq. 28) and the absorbing barrier stay at +4.9 … +5.8%. SteadyStateOlzmann's k_uni·K moves to +7.5 … +8.5%.
- **CSE's own pair** k(P1 → W1)/k(W1 → P1) departs from K = k∞,a/k∞,d by −2.4% at 300 K, −2.1% at 500 K, −1.2% at 750 K and −0.6% at 1000 K. Before the truncation, CSE and Olzmann agreed to 7·10⁻⁵ there.

**Reading.** The truncated grains carry part of the thermal population of a small molecule. With them removed, the thermal distribution is renormalized over fewer states, and thermal rate coefficients rise by Q_all/Q_kept.
- In Case 2 (four large wells), the 4–12 truncated grains carry almost no population, so the effect is negligible.
- MESS's reference runs use the default kernel mode, which never truncates. The only "truncating" lines in the Case 2 MESS log concern the states of barrier B36, not a well.

**Decision needed again.** Options given to Peter: keep the truncation; MESMER's reservoir state; or back to (b).

## 11. Decision (Peter, 2026-10-06): MESMER's reservoir state

**Source read:** MESMER 7.1 manual, Sec. 14.2.1 "The Reservoir State Approximation", and `src/gWellProperties.h`, `constructReservoir`.
- **What it assumes:** "This method assumes that significant portions of low energy molecular phase space are in a Boltzmann distribution throughout the course of the reaction. It is usually appropriate for grains which are more than a few kT below the lowest reaction threshold".
- **How it works:** "the bimolecular source term represents a collection of grains that are represented with one grain because we assume that these grains are always thermalized … The reservoir state can be formulated by analogy".
- **Equations:**
  - Deactivation into the reservoir from grain E: k_d(E) = Σ over the reservoir grains of the normalized downward probabilities P(i|E).
  - Activation out of it follows from detailed balance: k_a x_B = k_d x_C (eq. 14.15), with k_d = Σ_E k_d(E) f(E) (eq. 14.16).
  - In the source: "k_a = k_d(E) * f(E) / x_r"; "upward transitions are determined as part of symmetrization".

**Why it removes both problems:**
- **No temperature step.** There is no reduction factor and no integer window: the normalization above the reservoir is the plain eq. 4.16.
- **No lost population.** The reservoir keeps its Boltzmann weight Q_res = Σ f_i, so partition functions, thermal distributions and k_uni stay complete.

### 11.1 What the reservoir state means physically, and why MarXus uses it

**Why there is a problem at the bottom of a well.** The exponential-down model fixes how much energy a deactivating collision removes, P(E' ← E) ∝ exp(−(E − E')/⟨ΔE_down⟩). The activating probabilities then follow from detailed balance, P(E ← E')/P(E' ← E) = ρ(E)/ρ(E') e^{−(E − E')/kT}.
- Near the bottom of a well the density of states is sparse and rises steeply.
- From such a low grain, the activating probabilities into the many states above can add up to more than 1.
- Normalization (Robertson 2019, eq. 4.16) then has no positive solution.
- Robertson (2019, p. 294) traces this to the sparsity of the states, not to numerical error: there the exponential-down model and detailed balance cannot both hold grain by grain.

**What the reservoir assumes.** Physically, the lowest part of a well is the region of thermalized molecules.
- There collisions exchange energy many times before any molecule can react, because the reaction thresholds lie many k_BT higher.
- So the populations of these grains stay in Boltzmann equilibrium with each other at all times. MESMER's manual: "significant portions of low energy molecular phase space are in a Boltzmann distribution throughout the course of the reaction".
- In that region the master equation does not need the grain-to-grain detail that the exponential-down model fails to give. It only needs the total population and its exchange with the grains above.

**What the reservoir state is.** These grains are represented by one state, a population N_res distributed as f_i/Q_res, with Q_res = Σ_i ρ_i e^{−E_i/kT}. It is a "stabilized, thermalized" pool at the bottom of the well:
- **In:** molecules arrive by collisional deactivation from the grains above, with the normalized downward probabilities, which do exist there.
- **Out:** they are activated back by detailed balance, so the equilibrium is exact.
- **Reactions:** any process from its grains (in MarXus also the deep tunneling reactions) goes with the thermal share f_i/Q_res.

**What it changes, and what not.**
- **Above the reservoir:** every grain keeps the plain exponential-down model, exactly normalized. Chemical activation, falloff and the competition between reaction and stabilization are computed in full detail.
- **The partition function of the well stays complete,** because the reservoir carries the full Boltzmann weight of its grains. Equilibrium constants, thermal distributions, thermal rate coefficients (k_uni) and CSE's detailed balance are therefore unaffected.
- **Truncation was different.** It deleted these grains, which in a small molecule such as C₂H₃ carry thermal population: k_uni rose by 2.5% at 300 K.
- **Only the internal energy-transfer detail is lost,** among the lowest grains themselves. That does not matter where the grains are thermalized anyway.

**When it is valid.** MESMER's manual: "usually appropriate for grains which are more than a few kT below the lowest reaction threshold, and when the rate of collisional deactivation is faster than the rate of reaction (which is usually the case at moderate pressures)".
- In MarXus the reservoir is only as large as the normalization requires. Its size is set by where eq. 4.16 fails, not chosen by the user as in MESMER.
- In the validations it is far below every threshold:
  - ZZ-allyl + O₂: 4–12 grains, 152–456 cm⁻¹ above the bottoms of wells whose thresholds lie thousands of cm⁻¹ higher;
  - C₂H₃: at 2000 K the reservoir top is 0.8 k_BT above the bottom, against the threshold at 9.8 k_BT.
- **Not implemented, only proposed to Peter:** a warning if a reservoir ever reaches within a few k_BT of a well's lowest threshold.

**Why it is the rule in MarXus.** Where eq. 4.16 holds, MarXus uses it plainly, with no reservoir. Where it fails, the reservoir is the only treatment considered here that keeps all three:
1. the exponential-down model exactly normalized above it;
2. detailed balance exactly;
3. the full thermal population.

**The alternatives:**
- MESMER without a reservoir gives negative probabilities;
- MESS's truncation deletes thermally populated states;
- the SSUMES-type reduction changes the kernel with a temperature-dependent integer window, which caused the step.

## 12. Implementation: the reservoir state

**Kernel** (`collision_kernels.rs`, `exponential_down_kernel`).
- **Boundary.** Eq. 4.16 is solved from the top grain down. At the first grain g where the activating probabilities alone reach 1 the back substitution stops. Grains 0..=g form the reservoir (`CollisionKernel::reservoir_grains = g + 1`).
- **Above the boundary.** Every grain keeps its normalization from the top, together with its downward transitions into reservoir grains. Its normalization needs only the higher grains, so it is exact (Σ_t P(t|j) = 1).
- **Unlike MESS's truncation** there is no renormalization loop. The reservoir boundary is MESS's failure criterion; the treatment below it is MESMER's.
- **The reservoir size is automatic.** In MESMER it is a user input (`me:reservoirSize`).

**Operator** (`chemical_activation_operator.rs`).
- **One state per reservoir.** A well with a reservoir has one state (w, 0) for all its reservoir grains (`index_of` maps them all to it), with the weight ln Q_res, where Q_res = Σ_{i≤g} ρ_i e^{−E_i/kT} on the absolute scale (`Reservoir { state, weights = f_i/Q_res, log_weight }`).
- **Into the reservoir.** Collisions from a grain t above it are the summed downward probabilities ω Σ_{i≤g} P(i|t); `index_of` sends them to the reservoir state. Isomerization into a reservoir grain of another well also goes into the reservoir state.
- **Out of the reservoir:**
  - activation into grain t: ω Σ_{i≤g} P(i|t) f_t/Q_res (detailed balance, MESMER eqs. 14.15–14.16);
  - reactions and isomerization: every reservoir grain i reacts with its share f_i/Q_res of the reservoir population. *This is a MarXus extension of the thermalized-reservoir assumption.* MESMER requires the reservoir to lie below the reaction thresholds. Here the deep Eckart tunneling of B34/B36 in ZZ-allyl + O₂ reaches the lowest grains, so these rates are not exactly zero.
  - the bimolecular sink k_c[D].
- **Detailed balance.** J stays exactly detailed balanced with the weights f_t and Q_res. For isomerization this needs ρ_a k_ab = ρ_b k_ba at every energy.
- **Absorbing barrier** (SteadyStateAbsorbingBarrier):
  - a barrier inside the reservoir would split a thermalized state, so it is an error, with an explanation;
  - a barrier at or above the top of the reservoir absorbs all of it, so there is no reservoir state.
- **Helpers:**
  - `state_rate(s, k)`: k at the grain, or Σ k_i f_i/Q_res for a reservoir;
  - `grain_populations(x)`: the reservoir population spread over its grains with f_i/Q_res;
  - `state_populations(grains)`: the inverse.

**Consumers** use the helpers, so every grain-resolved quantity (Σ_E k(E) N(E), distributions, ⟨E⟩) is exact:
- the observables of the steady states and the thermal eigenpair (`grain_populations`);
- the CSE product rates p^(ν) and the time-integration exits (`state_rate`);
- the source projection (weights on reservoir grains add up in the reservoir state; nothing is lost).

**Output.** RUN SETTINGS, "Collisions:", lists per temperature the wells with a reservoir and the number of its grains (with their energy above the bottom).

**Tests.** All were written before the code. Only the kernel test was run red first. The operator and source tests compiled for the first time together with the implementation, because the field and the functions did not exist before.
- `exponential_down_lumps_the_grains_from_the_first_failure_down_into_a_reservoir`: against an independent eq. 4.16 stopped at the first failure. It was red on the truncation code: 22 truncated vs 19 reservoir grains.
- `low_energy_reservoirs_list_only_the_wells_where_eq_4_16_fails`.
- `the_reservoir_is_one_thermalized_state_with_the_boltzmann_weight_of_its_grains`: the weight ln Q_res, `state_rate` as the Boltzmann average, `grain_populations`/`state_populations`, activation out of the reservoir, no absorption.
- `an_absorbing_barrier_within_the_reservoir_is_an_error`, and a barrier at its top absorbs it.
- `isomerization_into_reservoir_grains_feeds_the_reservoir_with_detailed_balance`.
- `source_on_reservoir_grains_goes_into_the_reservoir_state`.
- **Restored to the complete well** (the truncation had made them partial):
  - the equilibrium-distribution test of the final steady state;
  - λ₁ → canonical k∞ at high pressure;
  - the CSE high-pressure and detailed-balance tests.
- **Double-precision tests moved to temperatures where their regime holds with the reservoir:**
  - sum rule between 10⁻⁵ and 1.5·10⁻²: 190 K (1.7·10⁻⁴);
  - λ₁ < 0 from the full decomposition: 160 K;
  - no plain Cholesky factor: still 130 K.

Full suite: 229 library and 43 binary tests pass.

## 13. Results with the reservoir state

In `reports/low_energy_reservoir_state.md`, Section 9:
- the 1 K scan is smooth, also where the reservoirs change by whole grains;
- the Case 2 MESS deviations are as with the truncation, and the thermal losses are now smooth as well;
- C₂H₃: the thermal population is complete again (CSE detailed balance 0.00% at 300–1000 K, dissociation +4.7 … +5.7% at 300 K), and the high-T results improved compared with the former reduction rule (2000 K: −2.6 … −2.1% instead of −7.8 … −4.5%).

## 14. References

- S. H. Robertson, Comprehensive Chemical Kinetics 43 (2019): eq. 4.16 (normalization), p. 294 (breakdown for sparse states, eqs. 4.60–4.61).
- `reports/chemical_activation_three_approaches.md`: variant V2 (low-energy reduction factors) and the row "Low-energy treatment" of the comparison table.
- `reports/collision_kernel_detailed_balance_and_normalization.md`: Section 5.2 (UNIMOL E₀/2 rule) and Section 8 (the open decision).
- `validation/ZZAllyl+O2_Gamma_Case2/cse_comparison.csv`, `four_methods_comparison.csv`.
- MESMER 7.1: manual, Sec. 14.2.1 "The Reservoir State Approximation" (eqs. 14.15–14.16); source `src/gWellProperties.h` (`constructReservoir`), `src/TMatrix.h` (`normalizeProbabilityMatrix`).
- MESS 2026: `src/libmess/mess.cc`, `MasterEquation::Well::_set_kernel` (kernel modes; truncation in mode `down`); inner-barrier sizes around line 2737.
