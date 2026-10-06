# The low-energy reservoir state of the exponential-down kernel

**Date:** 2026-10-06. **Status:** implemented, tested (229 library and 43 binary tests pass) and validated on both systems (Section 9). Decided by Peter on 2026-10-06. How this rule was reached, step by step, is in `reports/low_energy_reduction_temperature_step.md`.

## Contents

1. Summary
2. The problem: the normalization fails at the bottom of a well
3. Physical meaning of the reservoir state
4. Formulation
5. Differences from MESMER's reservoir state
6. Why this rule: the alternatives
7. Validity and limits
8. Implementation and tests
9. Validation results
10. References

## 1. Summary

**What happens.**
- **Where the normalization works.** MarXus's exponential-down kernel is normalized exactly by back substitution from the top grain (Robertson 2019, eq. 4.16). Where that works, nothing else is done.
- **Where it fails.** In the sparse lowest grains of a well it fails. There the grain at which it fails and all grains below it form **one thermalized state, the reservoir state**, as in MESMER (manual, Sec. 14.2.1).
- **What the reservoir keeps.** It keeps the full Boltzmann weight of its grains. Collisions bring molecules into it with the exactly normalized downward probabilities of the grains above, and take them out by detailed balance.

**Why this rule.** Of the treatments considered (Section 6), it is the only one that keeps all three:
1. the exponential-down model, exactly normalized, at every grain above the reservoir;
2. detailed balance exactly;
3. the complete thermal population of the well.

**What it fixed:**
- **Case 2.** A step of +3% in R → G4 of ZZ-allyl + O₂ at 304.7 K, in every method. It came from the former reduction rule.
- **C₂H₃.** A rise of k_uni by 2.5% at 300 K from MESS-type truncation, which removed thermally populated grains.

## 2. The problem: the normalization fails at the bottom of a well

**The model** (Robertson 2019, eqs. 4.7, 4.11):
- Deactivating collisions: P(t ← j) = A_j exp(−(E_j − E_t)/⟨ΔE_down⟩), E_t ≤ E_j.
- Activating collisions, by detailed balance: P(t ← j) = A_t (ρ_t/ρ_j) exp(−(E_t − E_j)(1/⟨ΔE_down⟩ + 1/kT)), E_t > E_j.

**The normalization** Σ_t P(t ← j) = 1 (eq. 4.16) is solved for A_j from the top grain down. It needs only the A_t of higher grains.

**Why it fails.** For a grain j at the sparse bottom of a well, ρ_t/ρ_j is large for the grains above: the density of states rises steeply from the zero-point level. The activating probabilities alone can then add up to 1 or more, and no positive A_j exists. Robertson (2019, p. 294) attributes this to the sparsity of the states, not to numerical error: the exponential-down form and detailed balance cannot both hold grain by grain there.

**How often it happens.** In both validation systems it fails at every temperature. At the first failure, from the top:
- **ZZ-allyl + O₂:** in the lowest 3–9 grains of 38 cm⁻¹, of about 430 per well.
- **C₂H₃:** in the lowest grains of 21 cm⁻¹, more of them at high T.

## 3. Physical meaning of the reservoir state

**The bottom of a well is the region of thermalized molecules.**
- Its energies lie many k_BT below the reaction thresholds.
- A molecule there is "stabilized": before it can react it undergoes many collisions, which exchange energy among the low grains far faster than any reaction removes molecules.
- The populations of these grains therefore stay in Boltzmann equilibrium with each other at all times. MESMER's manual: "significant portions of low energy molecular phase space are in a Boltzmann distribution throughout the course of the reaction".

**What the master equation needs from this region.** It needs the total population there and its exchange with the grains above, not the energy transfer between the low grains themselves. That internal detail is exactly what the exponential-down model fails to provide there (Section 2).

**The reservoir state.** One population N_res, distributed over its grains as f_i/Q_res, with f_i = ρ_i e^{−E_i/kT} and Q_res = Σ_i f_i. It acts as a pool of thermalized molecules at the bottom of the well:
- **Into the pool:** collisional deactivation from the grains above, with their normalized downward probabilities (which exist).
- **Out of the pool:** activation back up, by detailed balance, so that the equilibrium between the pool and the grains above is exact.
- **Reactions:** a process from a reservoir grain i goes with the share f_i/Q_res of the pool (in MarXus; MESMER has no reactions out of the reservoir, Section 5).

**What it changes, and what it does not:**
- **The grains above** keep the full, exactly normalized exponential-down model. Chemical activation, falloff, and the competition of reaction and stabilization are unchanged.
- **The partition function of the well is complete.** The pool carries the Boltzmann weight of all its grains. Equilibrium constants, thermal distributions, thermal rate coefficients (k_uni, λ₁) and the detailed balance of CSE's rate coefficients are not affected.
- **Lost:** only the energy-transfer detail among the reservoir grains themselves. That does not matter where these grains are thermalized.

## 4. Formulation

**Notation.** Grains 0 … g form the reservoir, where g is the first grain from the top at which eq. 4.16 fails. Grains t > g are "active", ω is the collision frequency, and f_i = ρ_i e^{−E_i/kT} on the absolute energy scale.

| term of J | value |
|---|---|
| weight of the reservoir state (symmetrization) | Q_res = Σ_{i≤g} f_i |
| collision, active t → reservoir | ω Σ_{i≤g} P(i ← t) (the normalized downward probabilities of t) |
| collision, reservoir → active t | ω Σ_{i≤g} P(i ← t) · f_t/Q_res (detailed balance) |
| reaction or isomerization out of the reservoir | Σ_{i≤g} k(E_i) f_i/Q_res, isomerization going to the grain of the target well at E_i |
| isomerization into a reservoir grain of another well | into the reservoir state |
| bimolecular sink | k_c[D], as for every grain |
| collisions within the reservoir | none (internal equilibrium) |

**Detailed balance.** J_{t,res} Q_res = J_{res,t} f_t exactly. For isomerization it holds whenever ρ_a k_ab = ρ_b k_ba at every energy (as MarXus builds the isomerization rates).

**Column sums.** The column sums of J are the losses out of the network, so the yields still add up to 1.

**Observables.**
- **Grain-resolved quantities** (Σ_E k(E) N(E), distributions, ⟨E⟩) use the reservoir population spread over its grains with f_i/Q_res. All methods therefore see the same well.
- **Sources** on reservoir grains (e.g. the Boltzmann source of the thermal fates) go into the reservoir state, so no weight is lost.

**SteadyStateAbsorbingBarrier.**
- An absorbing barrier inside the reservoir would split a thermalized state, so it is an error, with an explanation.
- A barrier at or above the top of the reservoir absorbs all of it, and the well then has no reservoir state.

## 5. Differences from MESMER's reservoir state

**Source read:** MESMER 7.1 manual, Sec. 14.2.1, and the manual's `me:reservoirSize` entry; source files `src/gWellProperties.cpp` (reservoir size), `src/gWellProperties.h` (`constructReservoir`, symmetrization), `src/TMatrix.h` (normalization), `src/IsomerizationReaction.cpp` and `src/IrreversibleUnimolecularReaction.cpp` (reaction terms).

| aspect | MESMER | MarXus |
|---|---|---|
| **Purpose** | efficiency: "it significantly truncates the size of the matrix that must be diagonalized … up to a factor of 30 faster" | the normalization of eq. 4.16 fails in the sparse lowest grains; a valid kernel is needed there |
| **When used** | optional, only if the user gives `me:reservoirSize` | automatically, only where eq. 4.16 fails; never otherwise |
| **Size** | set by the user as an energy, from the well bottom up, or a given distance below the lowest threshold; capped at the lowest barrier ("corrected according to the lowest barrier height") | the first failing grain from the top and all grains below it: as small as the normalization requires |
| **Temperature dependence of the size** | fixed energy, the same grains at every T | follows the failure point: e.g. C₂H₃ 2 grains (42 cm⁻¹) at 300 K, more at higher T |
| **Wells** | "applies only to isomer wells" | every well |
| **Normalization of the grains above** | back substitution from the top over the complete grid (`normalizeProbabilityMatrix`), without a failure check | the same back substitution from the top, stopped at the first failure; for the grains above the reservoir the coefficients are identical, because they need only higher grains |
| **Deactivation into the reservoir** | sum of the normalized downward probabilities into the reservoir grains (`constructReservoir`) | the same |
| **Activation out of the reservoir** | "determined as part of symmetrization", with the reservoir weight = Σ f of its grains | explicitly, ω Σ P(i ← t) f_t/Q_res: the same detailed-balance value |
| **Reactions out of reservoir grains** | none: the reaction terms start at the reaction threshold, and the reservoir is capped below the lowest barrier | each reservoir grain reacts with its share f_i/Q_res (MarXus extension, below) |
| **Isomerization into reservoir grains** | not reached (reservoir below the thresholds) | into the reservoir state, with detailed balance |
| **Absorbing barrier** | no such method in MESMER | a barrier inside the reservoir is refused; above it, it absorbs the reservoir |
| **Reporting** | the reservoir size in the input | RUN SETTINGS lists per temperature the wells with a reservoir and their grains (number and cm⁻¹) |

**The MarXus extension: reactions out of the reservoir.**
- **Why it is needed.** MESMER needs no rule for reactions inside the reservoir because it keeps the reservoir below the thresholds. In MarXus the reservoir size is not chosen; it follows from the normalization. Exact Eckart tunneling also gives small but nonzero k(E) far below the classical threshold, down to the higher of the two well bottoms (ZZ-allyl + O₂: B34 and B36 into the lowest grains of G3). So a rule is needed.
- **The rule.** Each reservoir grain reacts with its thermal share f_i/Q_res. That is the consequence of the reservoir's own assumption (Boltzmann equilibrium within it).
- **Status.** This rule is not in MESMER; it is documented here as MarXus's.
- **Size of the effect.** In the validations these rates are many orders of magnitude below the collision frequency.

## 6. Why this rule: the alternatives

**The question.** What to do where eq. 4.16 fails. The alternatives considered and their measured consequences (details in `low_energy_reduction_temperature_step.md`):

| treatment | where | consequence |
|---|---|---|
| reduction factors (ρ_{i+n}/ρ_i)^m below a cut, n = ⌊1.5⟨ΔE_down⟩/ΔE⌋ + 1 | SSUMES (`umolProb::setETProb`); MarXus until 2026-10-06 | the integer window n jumps with T: a step of +3% in R → G4 (ZZ-allyl + O₂, 304.7 K), −5% in the lowest relaxation eigenvalue |
| the same reduction, only where eq. 4.16 fails ("option b") | — | eq. 4.16 fails in every well of both systems, so the same step remains |
| no treatment; coefficients become negative | MESMER default | negative transition probabilities |
| per-pair common factor, not normalized | MESS default kernel mode | eq. 4.16 not satisfied (rejected earlier, `collision_kernel_detailed_balance_and_normalization.md`) |
| truncation: delete the failing grain and all below, renormalize | MESS kernel mode `down` | removes thermally populated grains in small molecules: C₂H₃ k_uni +2.5% at 300 K, CSE detailed balance −2.4% |
| grains below E₀/2 share the coefficient of the grain above | UNIMOL (older MarXus multiwell path) | normalization only approximate below E₀/2 |
| **reservoir state** | **MESMER (optional); MarXus (automatic)** | **no step; complete partition function; exact normalization above, exact detailed balance** |

## 7. Validity and limits

**The approximation assumes** that the reservoir grains are in Boltzmann equilibrium among themselves. MESMER's manual: "usually appropriate for grains which are more than a few kT below the lowest reaction threshold, and when the rate of collisional deactivation is faster than the rate of reaction (which is usually the case at moderate pressures)".

**In the validations** the reservoirs lie far below every threshold:
- ZZ-allyl + O₂: a few to a dozen grains of 38 cm⁻¹, against thresholds thousands of cm⁻¹ higher;
- C₂H₃: the reservoir top at 2000 K is below 1 k_BT above the bottom, against the threshold at 9.8 k_BT.

**Possible check, not implemented** (proposed to Peter): a warning when a reservoir reaches within a few k_BT of a well's lowest threshold.

**The size moves with T.** It changes by whole grains where the failure point moves. Inside a thermalized reservoir this should not matter. It is checked with the 1 K temperature scan of ZZ-allyl + O₂ (Section 9).

## 8. Implementation and tests

**Code:**
- `src/masterequation/collision_kernels.rs`, `exponential_down_kernel`: back substitution from the top, stopped at the first failure; `CollisionKernel::reservoir_grains`.
- `src/masterequation/chemical_activation_operator.rs`:
  - `Reservoir { state, weights = f_i/Q_res, log_weight = ln Q_res }`;
  - the reservoir column of J and `index_of` mapping all reservoir grains to the reservoir state;
  - the helpers `state_rate`, `grain_populations` and `state_populations`;
  - `low_energy_reservoirs(network, T, model)`;
  - the absorbing-barrier rule.
- **Consumers:**
  - `chemical_activation_steady_state.rs` (`project_source`: reservoir weights add up);
  - `chemical_activation_observables.rs` and `chemical_activation_eigen.rs` (`grain_populations`);
  - `chemically_significant_eigenvalues.rs` and `direct_time_integration.rs` (`state_rate`).
- `examples/chemical_activation_from_deck.rs`: RUN SETTINGS, "Collisions:" lines.

**Tests:**
- **Kernel:** `exponential_down_lumps_the_grains_from_the_first_failure_down_into_a_reservoir`, against an independent eq. 4.16 stopped at the first failure (first run red on the earlier truncation code: 22 vs 19 grains).
- **Operator:**
  - `low_energy_reservoirs_list_only_the_wells_where_eq_4_16_fails`;
  - `the_reservoir_is_one_thermalized_state_with_the_boltzmann_weight_of_its_grains`;
  - `an_absorbing_barrier_within_the_reservoir_is_an_error`;
  - `isomerization_into_reservoir_grains_feeds_the_reservoir_with_detailed_balance`;
  - the existing `j_is_detailed_balanced_with_boltzmann_weights_on_the_absolute_energy_scale` and `column_sums_of_j_equal_the_losses_out_of_the_network`, now with reservoirs in both test wells.
- **Sources:** `source_on_reservoir_grains_goes_into_the_reservoir_state`.
- **Physics restored to the complete well:**
  - the final steady state through the only channel is the equilibrium distribution on all grains;
  - λ₁ → canonical k∞ over all grains at high pressure;
  - the CSE high-pressure limits and detailed balance with the full partition functions.
- **Test-first record.** All tests were written before the code. The kernel test was run red first; the operator and source tests compiled for the first time together with the implementation.

## 9. Validation results

**Runs.** Both validation directories were re-run with the reservoir state (2026-10-06; 4 cores):
- ZZ-allyl + O₂ Case 2: 4 min 18 s, all four methods, both tunneling models;
- C₂H₃: 13 min 46 s, all four methods and the three eigen-solvers.

### 9.1 Reservoir sizes (RUN SETTINGS, "Collisions:")

**ZZ-allyl + O₂ Case 2** (38 cm⁻¹ grains):

| T (K) | G2 | G3 | G4 | G6 |
|---|---|---|---|---|
| 270 | 6 | 6 | 5 | 4 |
| 300 | 8 | 8 | 6 | 5 |
| 330 | 10 | 10 | 9 | 8 |

That is 152–380 cm⁻¹ above the well bottoms.

**C₂H₃** (21 cm⁻¹ grains):

| T (K) | 300 | 500 | 750 | 1000 | 1250 | 1500 | 1750 | 2000 |
|---|---|---|---|---|---|---|---|---|
| grains | 2 | 2 | 3 | 4 | 5 | 36 | 40 | 44 |
| top (cm⁻¹ above the bottom) | 42 | 42 | 63 | 84 | 105 | 756 | 840 | 924 |

The 2000 K top is 0.67 k_BT; the threshold lies at 9.8 k_BT.

### 9.2 No temperature step: 1 K scan of ZZ-allyl + O₂ (CSE, MESS Eckart model, 760 Torr, 300–310 K)

**Setup.** The same deck and grid as in `low_energy_reduction_temperature_step.md`, Section 3. In this interval the reservoirs change by whole grains: G4 from 6 to 8 grains at 301–302 K, G2 from 8 to 9 at 302 K, G6 from 5 to 7 at 306 K.

**Result.** The ratio of each quantity to its value 1 K lower is smooth everywhere:

| quantity | ratios 301 … 310 K | spread |
|---|---|---|
| lowest relaxation eigenvalue | 0.9947 … 0.9948 | 0.0002 |
| k(R → G2) | 0.9829 → 0.9823 (linear) | 0.0006 |
| k(R → G3) | 0.9898 → 0.9892 | 0.0006 |
| k(R → G4) | 0.9879 → 0.9871 | 0.0008 |
| k(R → P5) | 1.0038 → 1.0028 | 0.0010 |
| thermal loss of G2 | 1.0216 → 1.0218 | 0.0002 |
| thermal loss of G6 | 1.0645 → 1.0588 (linear) | 0.0057 |

**With the former reduction rule** the same scan had steps at 304 → 305 K: k(R → G4) ratio 1.0188 against a trend of 0.9877, and the relaxation eigenvalue 0.9454 against 0.9967.

### 9.3 ZZ-allyl + O₂ Case 2: CSE against MESS (all 21 conditions; `cse_comparison.csv`)

| entry | former reduction rule | MESS truncation | **reservoir state** |
|---|---|---|---|
| R → P5 (IEPOX + OH) | −2.7 … −3.5% | −3.0 … −3.6% | **−3.0 … −3.6%** |
| R → G2 | +2.1 … +3.3% | +3.0 … +5.2% | **+3.0 … +5.2%** |
| R → G3 | +0.9 … +3.6% | −1.2 … −0.8% | **−1.2 … −0.7%** |
| R → G4 | +1.4 … +6.1% (step) | −1.0 … −0.3% | **−1.0 … −0.3%** |
| R → P1 / R → P7 | | | **−2.3 … −1.5% / −5.1 … −3.9%** |
| thermal losses G2 / G3 / G4 / G6 | | | **+0.8 … +1.1% / −0.04 … +0.4% / +0.6% / +2.9 … +4.7%** |

**Smoothness.** The second differences of ln k(T) at 760 Torr follow MESS, now also for the thermal losses of the wells. G6: −0.0775, −0.0698, −0.0629 against MESS's −0.0774, −0.0699, −0.0627. With the truncation they still wobbled by ±0.4%.

### 9.4 C₂H₃: the thermal population is complete again (`method_comparison.csv`)

| quantity | MESS truncation | **reservoir state** |
|---|---|---|
| CSE's own pair against K: [k(P1 → W1)/k(W1 → P1)]/K − 1, 300–1000 K | −2.4 … −0.6% | **0.00 … −0.01%** |
| dissociation k_uni (= CSE k(W1 → P1)) vs MESS, 300 K | +7.3 … +8.3% | **+4.7 … +5.7%** (the known exact vs semiclassical Eckart offset) |
| association, SteadyStateOlzmann k_uni·K vs MESS, 300 K | +7.5 … +8.5% | **+4.9 … +5.8%**, equal to CSE and to the absorbing barrier |

**C₂H₃ by temperature**, against MESS. Each cell gives the range over the five pressures:

| T (K) | dissociation (k_uni = CSE) | association CSE | association Olzmann (k_uni·K) | association absorbing barrier |
|---|---|---|---|---|
| 300 | +4.7 … +5.7% | +4.9 … +5.8% | +4.9 … +5.8% | +4.9 … +5.8% |
| 1000 | −0.8 … +0.5% | −0.6 … +0.5% | −0.6 … +0.5% | −0.6 … +0.5% |
| 1500 | −3.5 … −1.9% | −3.1 … −1.7% | −1.4 … −1.2% | −10.4 … −5.3% |
| 2000 | −2.6 … −2.1% | −2.6 … −2.1% | +5.5 … +13.2% | not defined |

**At 2000 K with the former reduction rule:** dissociation −7.8 … −4.5%, CSE association −9.1 … −4.9%. The reduction factors had acted on 300 grains (6300 cm⁻¹) of C₂H₃ at 2000 K. The reservoir now holds 44 grains.

**Detailed balance of the CSE pair at high temperature.** CSE's own pair departs from K as MESS's own pair does:

| T (K) | MarXus CSE | MESS |
|---|---|---|
| 1250 | −0.05 … −0.21% | −0.13 … −0.58% |
| 1500 | −0.5 … −1.7% | −0.7 … −2.1% |
| 1750 | −2.6 … −6.1% | −2.6 … −6.1% |
| 2000 | −7.2 … −14.0% | −7.2 … −13.9% |

This is a property of the CSE rate coefficients when the separation of the time scales is poor (Λ₁/Λ₂ = 0.065–0.10 at 2000 K), not an error of either code. SteadyStateOlzmann imposes K exactly, k(P1 → W1) = k_uni·K. Its association therefore lies above MESS's at 2000 K, while the dissociation agrees within −2.6 … −2.1%.

### 9.5 Identities between the methods (both systems; `method_comparison.csv`, `reports/method_comparison.md`)

**All exact identities hold** where λ₁ is above the double-precision floor:
- the time integration reaches the final steady state in all printed digits (81 exits in Case 2, 28 conditions in C₂H₃);
- the late-time decay of the pulse equals k_uni within 3·10⁻⁶;
- the CSE capture and loss balances hold within 6·10⁻⁷;
- the one-well CSE k(W → P) equals k_uni within 5·10⁻⁵.

**The exceptions** are the C₂H₃ conditions at 300–500 K, where λ₁ lies 40 to 10¹³ times below the floor (known since the eigen-solver study). There the conservation of the pulse is 2·10⁻⁴ … 6·10⁻⁴, and CSE's loss balance and λ₁ are noise.

## 10. References

- S. H. Robertson, Comprehensive Chemical Kinetics 43 (2019): eqs. 4.7, 4.11, 4.16 (exponential down and its normalization), p. 294 (breakdown for sparse states).
- MESMER 7.1 manual: Sec. 14.2.1 "The Reservoir State Approximation" (eqs. 14.15, 14.16) and the `me:reservoirSize` entry. MESMER 7.1 source: `src/gWellProperties.cpp`, `src/gWellProperties.h` (`constructReservoir`, `collisionOperator`), `src/TMatrix.h` (`normalizeProbabilityMatrix`), `src/IsomerizationReaction.cpp`, `src/IrreversibleUnimolecularReaction.cpp`.
- MESS 2026 source: `src/libmess/mess.cc`, `MasterEquation::Well::_set_kernel` (kernel modes: default per-pair factor; `down`; truncation).
- `reports/low_energy_reduction_temperature_step.md`: the finding of the temperature step and the decisions.
