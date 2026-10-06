# CSE: merging of species that are not kinetically distinct (plan)

**Date:** 2026-10-06. **Status:** sources read; plan, not yet implemented.

**Request (Peter, 2026-10-06):**
- "when CSE is used and we have similar time-scales and chemical activation rates similar to relaxation rates, then some merging might happen and also when some species are not distinguishable … we also should use the merging process described in the 2013 master equation paper of Georgievskii and Klippenstein and also how they use it in MESS".
- "if such case happens when the users want to run the CSE module, then the output must warn about this, also it does not make sense to distinguish some species that are involved in this".

**Current MarXus behaviour.**
- `chemically_significant_eigenvalues.rs` always takes N = number of wells chemical eigenvalues.
- When Λ_N/Λ_{N+1} > 0.1 it only warns: "wells in fast equilibrium should be merged (Georgievskii et al. 2013, Sec. IV)". Merging is not implemented.
- In the validations the warning appears for ZZ-allyl + O₂ at 310–330 K (Λ_N/Λ_{N+1} up to 0.137) and for C₂H₃ at 2000 K and 76 Torr (0.10).

## 1. The paper: Georgievskii, Miller, Burke, Klippenstein, J. Phys. Chem. A 117, 12146 (2013), Sec. IV

- **When.** "Under the conditions of high temperatures and/or low pressures, some of the chemical eigenvalues may approach the energy relaxation limit. When this happens, the projection of the corresponding eigenvectors onto the relaxational subspace increases and eventually this projection becomes the dominant one."
- **What it means.** "The situation where some of the chemical eigenvalues approach the energy relaxation limit implies that the equilibration between the corresponding groups of species occurs so rapidly that it cannot be separated from the energy relaxation processes. The corresponding eigenstates should then be considered as effectively relaxational ones and the effective number of species should be reduced. Those species that are in equilibrium with each other should be united and treated as one and the dimensionality of the effective chemical subspace should be reduced as well."
- **Merged species.** "Then, the expressions for the phenomenological rate coefficients given by eqs 27, 28, and 30 remain valid if one interprets the ith isomer to mean the united species g", with
  - |g⟩ = f_j⁽⁰⁾(E)/√Q_g for j ∈ g, 0 otherwise (eq. 31);
  - Q_g = Σ_{j∈g} Q_j (eq. 32);
  - "The matrix M, eq 25, should be modified accordingly."
- **Eq. 21** (bimolecular-to-bimolecular) "remains formally unchanged". Its sum over the relaxational eigenstates now includes the eigenstates that left the chemical subspace.
- **Wells in equilibrium with bimolecular species** (eqs. 33–34): n_i/(n_A n_B) = κ_iν Q_i/Q_ν, with κ_iν = (1/Q_i) Σ_λ ⟨i|λ⟩ p_λ^(ν)/Λ_λ over the relaxational eigenstates. κ_iν is "close to unity if the ith isomer is in equilibrium with the νth bimolecular species and to zero otherwise".

## 2. MESS (2026 source, `src/libmess/mess.cc`; keywords in `src/mess_driver.cc`)

**Number of chemical eigenvalues** (around line 1138):
- **Keyword.** `ChemicalEigenvalueMax` (and `ChemicalThreshold`) set `MasterEquation::chemical_threshold` (mess_driver.cc, around line 1113). Both validation decks use `ChemicalEigenvalueMax 0.2`.
- **0 < threshold < 1 (absolute):** the chemical eigenvalues are those with Λ ≤ threshold × `relax_eval_min`, the minimal relaxation eigenvalue. MESS computes that one per well from the relaxation problem alone; its log prints "minimal relaxation eigenvalue" per well.
- **threshold > 1:** by the relaxational projection of the eigenvectors.
- **−1 < threshold < 0:** by the ratio of successive eigenvalues.

**Partition of the wells** (`threshold_well_partition`, around line 4400), used when chem_size < number of wells:
1. **pop_chem.** pop_chem(w, λ) is the population of well w in chemical eigenvector λ, G13's M_wλ (eq. 25).
2. **Projections.** The projection of a group g onto the chemical subspace is Σ_λ (Σ_{w∈g} √Q_w pop_chem(w, λ))²/Q_g = |P_chem |g⟩|², with |g⟩ of eq. 31 (`projection`, around line 4653).
3. **Primary wells.** Each single well's projection is computed. Taken in decreasing order, these are the primary wells: at least chem_size of them, plus every further well with projection ≥ `WellProjectionThreshold` (default 0.2).
4. **Best partition.** Every partition of the primary wells into chem_size groups is tried (`PartitionGenerator`); the one with the largest total projection is kept.
5. **Remaining wells.** Each is added greedily, the best (well, group) pair first, to the group whose projection it increases most, while the increase is positive.
6. **Bimolecular group.** Wells that are left over form the "bimolecular group": in equilibrium with bimolecular species, with no rate coefficients of their own.
7. **Log.** MESS prints "well partition", "bimolecular group" and "partition projection error" = chem_size − total projection.

**Rate coefficients.** The chemical eigenvectors are converted to the group basis (`basis(well_partition)`), and G13 eqs. 27, 28 and 30 are evaluated for the groups. The output names a merged species by its wells joined with "+".

## 3. Plan for MarXus

1. **Number of chemical eigenvalues, as MESS:** Λ_λ ≤ `ChemicalEigenvalueMax` × the minimal relaxation eigenvalue.
   - The deck keyword `ChemicalEigenvalueMax` is read; its default must be checked in mess_driver.cc.
   - The minimal relaxation eigenvalue must be computed as MESS does. *To check before implementing:* how `relax_eval_min` is obtained (per-well relaxation problem, around line 941).
2. **Partition, as MESS:** `threshold_well_partition`, with `WellProjectionThreshold` (default 0.2) as a deck keyword.
3. **Rate coefficients of the groups:** G13 eqs. 27, 28 and 30 in the group basis (eqs. 31–32), and eq. 21 over all non-chemical eigenstates.
4. **Wells in the bimolecular group:** no rate coefficients. Their κ_iν (eq. 34) are reported as the equilibrium with the bimolecular species.
5. **Output:**
   - a **warning** per condition where species are merged, with the eigenvalues, the threshold, the partition and the projection error;
   - the tables use the **merged names** (e.g. `G2+G4`), with no separate rate coefficients for merged wells, because "it does not make sense to distinguish" them (Peter);
   - the yields of merged species are those of the group.
6. **Tests (TDD):**
   - a two-well network in fast equilibrium (low isomerization barrier) merges into one group, and its rate coefficients equal those of the one-well network with the combined density of states;
   - separated wells are not merged;
   - a well in fast equilibrium with the products goes into the bimolecular group;
   - the identities of `method_comparison.py` (capture balance, long-time yields = final steady state) still hold after merging.

**Found after the first reading:**
- **relax_eval_min.** In `direct` it is `eigenval[well_size()]` (mess.cc, lines 662 and 941): the (N+1)-th eigenvalue of the full kinetic matrix, Λ_{N+1}, the "lowest relaxation eigenvalue" that MarXus already prints. MESS's rule is therefore Λ_λ ≤ `ChemicalEigenvalueMax` × Λ_{N+1}.
- **Consequence for the validations.** With `ChemicalEigenvalueMax 0.2`, as in both validation decks, neither system merges: the largest Λ_N/Λ_{N+1} is 0.137 (ZZ-allyl + O₂, 330 K, 500 Torr) and 0.10 (C₂H₃, 2000 K, 76 Torr). MarXus's present warning threshold, 0.1, is stricter than this criterion.
- **Default.** The initial value in mess.cc is −2, which leads to the error "chemical threshold has not been initialized properly" unless a keyword is given (the deck must set it).

**Open before implementing (to read):**
- MESS's default `ChemicalEigenvalueMax` and its `relax_eval_min`;
- how MESS's output tables name merged species;
- whether MESS also merges a well into the reactant (bimolecular group) in the bimolecular-to-well rates.
