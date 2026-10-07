# CSE: phenomenological rate coefficients from the chemically significant eigenvalues

[← README](../../README.md) · the four methods: [SteadyStateOlzmann](steady_state_olzmann.md) · [SteadyStateAbsorbingBarrier](steady_state_absorbing_barrier.md) · **CSE** · [TimeIntegration](direct_time_integration.md)

**Family:** eigenvalue methods. **Run with:** `Method CSE`, or `--method cse`.

## 1. The question it answers

Which species-to-species rate coefficients define a kinetic model that reproduces the master-equation kinetics after collisional relaxation? The species are reactants, wells, products and sinks.

The result is a set of phenomenological rate coefficients for a kinetic mechanism, not flux coefficients. MK06 (pp. 10529, 10531): "application of the steady-state approximation … is virtually always an attempt to equate a phenomenological rate coefficient to a flux coefficient. Sometimes this is a valid approach, and sometimes it is not."

## 2. Equations

The formulation of G13 (bimolecular reactant as a thermal source, eq. 1; MK06 for the concept):

```math
\frac{d|f\rangle}{dt} = -\hat{𝐆}\,|f\rangle + \sum_\nu s_\nu\,|p^{(\nu)}\rangle, \qquad \hat{𝐆} = 𝐉, \qquad \hat{𝐆}\,f^{(\lambda)} = \Lambda_\lambda\,f^{(\lambda)} .
```

**Chemically significant eigenvalues.** For $`n_w`$ wells the $`n_w`$ lowest eigenvalues are chemically significant when they are separated from the relaxation eigenvalues, $`\Lambda_{n_w} \ll \Lambda_{n_w+1}`$ (MK06 eq. 19). As in MESS, the chemical eigenvalues are those, counted from the lowest, with

```math
\Lambda_\lambda \le c_{\max}\,\Lambda_{n_w+1}, \qquad c_{\max} = \texttt{ChemicalEigenvalueMax}\ (0 < c_{\max} < 1) .
```

Their number $`n \le n_w`$ is the number of kinetic species.

**Definitions:**
- $`Q_i = \sum_{E\in i} f^0(E)`$;
- $`p^{(\nu)}_\lambda = \sum_E f^{(\lambda)}(E)\,k_{\to\nu}(E)`$ (eq. 15);
- $`M_{i\lambda} = Q_i^{-1/2}\sum_{E\in i} f^{(\lambda)}(E)`$ (eq. 25), the elements of $`𝐌`$;
- $`𝚲 = \mathrm{diag}(\Lambda_1,\dots,\Lambda_{n_w})`$.

```math
k_{j\to i} = -\sqrt{Q_i/Q_j}\,\big(𝐌𝚲𝐌^{-1}\big)_{ij} \;\; \text{(eq. 27)}, \qquad k_{i\to\nu} = Q_i^{-1/2}\sum_{\lambda\le n_w}\big(𝐌^{-1}\big)_{\lambda i}\,p^{(\nu)}_\lambda \;\; \text{(eq. 30)},
```

```math
k_{R\to i} = \frac{\sqrt{Q_i}}{Q_R}\sum_{\lambda\le n_w} M_{i\lambda}\,p^{(R)}_\lambda \;\; \text{(eq. 28, bimolecular-to-well)}, \qquad k_{R\to\mu} = \frac{1}{Q_R}\sum_{\lambda>n_w}\frac{p^{(\mu)}_\lambda\,p^{(R)}_\lambda}{\Lambda_\lambda} \;\; \text{(eq. 21, bimolecular-to-bimolecular)} .
```

**Reactant normalization** (eq. 23): $`1/Q_R = k_\infty / \sum_E k_{\to R}(E)\,f^0(E)`$, with $`k_\infty`$ the capture rate coefficient. No partition function of the reactants is needed.

**Product rows.** Eqs. 28, 21 and 22 hold for every bimolecular species, not only the reactant.
- For a product P, $`1/Q_P = k^{(c)}_P / \sum_E k_{\to P}(E)\,f^0(E)`$ (eq. 23), with $`k^{(c)}_P`$ the high-pressure rate coefficient of the reverse association P → wells.
- This gives k(P → W), k(P → P′) and the return k(P → P).

**Isomer–bimolecular equilibrium coefficients** (eqs. 33–34), per well $`w`$ and bimolecular channel $`\nu`$:

```math
\frac{n_w^{(\nu)}}{n_A n_B} = \kappa_{w\nu}\,\frac{Q_w}{Q_\nu}, \qquad \kappa_{w\nu} = \frac{1}{\sqrt{Q_w}}\sum_{\lambda > n}\frac{M_{w\lambda}\,p^{(\nu)}_\lambda}{\Lambda_\lambda} .
```

- κ is close to 1 when the well is in equilibrium with ν, and close to 0 otherwise.
- Because the collisions conserve $`f^0`$, $`(1/\sqrt{Q_w})\sum_{\text{all }\lambda}M_{w\lambda}\sum_\nu p^{(\nu)}_\lambda/\Lambda_\lambda = 1`$ over all loss channels. Its deviation is reported as a check.

**Species merging** (G13 Sec. IV), when $`n < n_w`$. Wells that equilibrate on the time scale of energy relaxation "should be united and treated as one". The wells are partitioned into $`n`$ groups $`g`$, the merged species, with

```math
|g\rangle = \sum_{w\in g}\sqrt{Q_w/Q_g}\;|w\rangle \;\; \text{(eq. 31)}, \qquad Q_g = \sum_{w\in g} Q_w \;\; \text{(eq. 32)}, \qquad M_{g\lambda} = \sum_{w\in g}\sqrt{Q_w/Q_g}\;M_{w\lambda} .
```

Eqs. 27, 28 and 30 above hold for the groups with $`\lambda \le n`$; eq. 21 runs over all other eigenstates, $`\lambda > n`$. Wells in no group are in equilibrium with the bimolecular species (MESS's "bimolecular group") and have no rate coefficients of their own.

**Partition**, as MESS (`threshold_well_partition`). The projection of a group on the chemical subspace is $`P_g = \sum_{\lambda\le n}\big(\sum_{w\in g}\sqrt{Q_w}\,M_{w\lambda}\big)^2/Q_g`$.
1. **Primary wells:** in decreasing order of their own projection, at least $`n`$ of them, plus every well with projection ≥ `WellProjectionThreshold`.
2. **Best partition:** of all partitions of the primary wells into $`n`$ groups, the one with the largest $`\sum_g P_g`$.
3. **Other wells:** one at a time, the best (well, group) pair first, into the group whose projection they raise most, while the raise is not negative. The rest form the bimolecular group.
4. **Partition projection error:** $`n - \sum_g P_g`$.

**Yields from the rate coefficients** (`reactant_yields`):
- **Prompt branching of R:** bimolecular-to-bimolecular + bimolecular-to-well, % of the net reaction.
- **Thermal fate of each well:** the absorbing chain $`B = (𝐈 - 𝐐)^{-1} 𝐀`$ of the well rate coefficients.
- **Long-time yields:** direct + through the wells. They equal the final steady state exactly (Section 7).

## 3. Algorithm

**Eigenpairs.** All eigenpairs of the symmetrized $`𝐒`$ come from LAPACK DSYEVD (default) or the in-house Householder/QL. Inverse iteration is refused, because it gives only the lowest eigenpair.

**Inverse of 𝐌.** $`𝐌^{-1}`$ is computed by Gauss–Jordan elimination with partial pivoting.

**Diagnostics:**
- the chemical eigenvalues;
- the number of species and, where wells are merged, a warning with the eigenvalues, the threshold, the species, the bimolecular group and the partition projection error;
- the lowest non-chemical eigenvalue $`\Lambda_{n+1}`$ and the separation $`\Lambda_n/\Lambda_{n+1}`$, with a warning above 0.1. This warning is stricter than the merging criterion and says so;
- the relaxational projections;
- the loss balance and detailed balance;
- the double-precision floor.

## 4. Settings

| deck keyword | option | default |
|---|---|---|
| `Method CSE` | `--method cse` | required (no default method) |
| `EigenSolver` | `--eigen-solver lapack\|full` | Lapack (inverse iteration refused) |
| `ChemicalEigenvalueMax` (global section, MESS's keyword) | – | 0.2 (0 < value < 1; MESS has no default and requires the keyword) |
| `WellProjectionThreshold` (global section, MESS's keyword) | – | 0.2 (as MESS) |

## 5. Output

**Species-to-species tables,** one per (T, p): `From\To`.
- **Rows:** the species (1/s), the reactant and every bimolecular product that is not Dummy (cm³/s).
- **Diagonal:** the total loss of a species; for a bimolecular species, capture − return (as MESS).
- **With each table:** its eigenvalues, separation, capture/return/net and warnings, then the κ matrix (wells × bimolecular channels) with the sum-rule deviation.

**Merged species** are named by their wells joined with "+" (e.g. `G2+G4`). Their wells get no separate rate coefficients, because they cannot be distinguished at that condition. MESS writes the rates of a group under the name of its deepest well (`Group::group_index`) and `***` for the others; the "+" names appear only in its log.

**Rate coefficients from every well** (1/s). Then, for the reactant:
- **Bimolecular-to-bimolecular rate coefficients** $`k(R\to P)`$ (cm³/s, eq. 21), with **yields** (% of the net reaction of R);
- **Bimolecular-to-well rate coefficients** $`k(R\to W)`$ (cm³/s, eq. 28), with **yields**;
- **Capture, return and net reaction of R**;
- **Thermal fate of each well** (%);
- **Long-time yields:** direct, through the wells, total (% of the eventual net reaction).
- **Rate coefficients from every product** P: P → species, P → other channels, capture, return, net (cm³/s).
- **κ** of every well and bimolecular channel.

**Prepared distributions** (`Preparation` block of the deck; [CSE source projection](../features/cse_source_projection.md)): every source F is projected onto the eigenvectors (G13 eqs. 13, 14, 24, 37–41).
- **Projection:** c_λ = ⟨f^λ|F⟩.
- **Species populations after the relaxation:** n_g = √Q_g Σ_chem M_gλ c_λ.
- **Prompt yields:** Σ_relax p_λ^ν c_λ/Λ_λ, formed during the relaxation.
- **Identity:** prompt + Σ_g n_g × (fate of g) = k_x^T 𝐉⁻¹F.
- **Precision check:** Σ n + Σ prompt = ΣF holds exactly with all eigenpairs; a violation is a warning (Frankcombe, Smith, J. Theor. Comput. Chem. 2, 179 (2003)).
- **Negative prompt yields** are physical: a preparation colder than the steady state reacts later than the species decaying from t = 0. For one well, ln(n)/λ₁ is the incubation time.
- **CSE validity in time** (`CompareWithCse`, with TimeIntegration): the projected description is propagated in time and compared with the master equation ([CSE validity in time](../features/cse_validity_in_time.md)).

All of these come in three views: by temperature, by pressure, and temperature–pressure. When the species differ between conditions, every species of every condition is a row, empty where it does not exist. The diagnostics table gives the number of species per condition.

## 6. Validity and limits

**Time scale.** The rate coefficients describe the kinetics for $`t \gg 1/\Lambda_{n_w+1}`$. They exist only while the eigenvalues are separated. Otherwise the wells are merged into species (G13 Sec. IV; Section 2), and the rate coefficients are those of the merged species.

**Rounding noise.** Entries many orders of magnitude below the largest of their row are rounding noise, and may be negative. MESS shows the same.

**What it does not give:** the time dependence before the relaxation is over, time profiles of sources and bath histories (TimeIntegration gives them). Non-thermal sources are projected onto prompt yields and post-relaxation species populations (Section 5). The long-time yields of a thermal source come out exactly; see Section 7.

## 7. Validation and identities

**Identity with the final steady state.** With all eigenpairs, the long-time yields reconstructed from the CSE rate coefficients equal $`k_x^T\,𝐉^{-1} F`$ of [SteadyStateOlzmann](steady_state_olzmann.md).
- Both are $`\sum_\lambda p^{(x)}_\lambda p^{(R)}_\lambda/(\Lambda_\lambda Q_R)`$ over all eigenpairs, independently of the separation.
- Test `cse_long_time_yields_equal_the_final_steady_state_yields`: 10⁻⁸.
- Case 2: IEPOX + OH at 300 K, 760 Torr is 2.45778% from both. Over all 21 conditions P5, the escape and P1 agree to within 6·10⁻⁷, the printed precision. P7, at most 0.008% of the reaction, agrees to 1·10⁻⁴: it is formed through G6, whose CSE entries are at the rounding level.

**One well.** The CSE well → product rate coefficient equals the eigenvector-average $`k_{\mathrm{uni}}`$ (test, 10⁻⁸).

**ZZ-allyl + O₂, four wells** (MESS Eckart model; `validation/ZZAllyl+O2_Gamma_Case2/`):
- The species tables of MESS are reproduced within 1.3% (entries above the rounding level; Neufeld collision integral, `reports/collision_integral_neufeld.md`).
- **Against MESS** (all 21 conditions): R → IEPOX + OH +0.28 … +0.50%; R → G2 +0.90 … +1.13%, R → G3 +0.39 … +0.57%, R → G4 +0.43 … +0.60%; R → P1 and P7 +0.17 … +0.46%; well → well −0.3 … +1.3%.
- **Low-energy reservoir state.** These numbers are with it (`reports/low_energy_reservoir_state.md`). The former reduction rule of the collision kernel had a step at 304.7 K, R → G4 +1.4 … +6.1%.
- **MESS's negative R → escape entry is reproduced** (−0.67 … +0.12%).
- **Identities:** capture balance within 4·10⁻⁷, loss balance within 1.3·10⁻⁷ (`reports/method_comparison.md`).

**Product rows, captures and κ** (`reports/cse_kappa_and_product_rates.md`):
- **Case 2,** MESS Eckart model, 21 conditions:
  - the rows of P1, P5 and P7 agree with MESS within −0.3 … +0.4% (P1, P5) and ±0.9% (P7; P7 → G4 −1.4%) for the entries above the rounding level;
  - the capture rate coefficients agree with MESS's high-pressure rows within +0.13 … +0.84%.
- **C₂H₄ + HO₂ hindered-rotor deck,** at the 16 conditions where MESS prints κ ≥ 0.05: MarXus agrees within 0.006, and κ(W2,P1) + κ(W2,P2) = 1 in both codes.

**H + C₂H₂ ⇌ C₂H₃, one well** (`validation/c2h3_mess_example/`, §4.4; 40 conditions, 300–2000 K):
- **CSE's k(W1 → P1) equals SteadyStateOlzmann's k_uni** in all printed digits wherever λ₁ is resolved (750–2000 K).
- **Up to 1000 K** the association k(P1 → W1) (eq. 28) equals k_uni·K within 0.02%.
- **Above that, CSE's own pair departs from detailed balance** as the separation Λ₁/Λ₂ grows: −0.5 … −1.7% at 1500 K (Λ₁/Λ₂ ≤ 0.015), −7.2 … −13.9% at 2000 K (0.045–0.067). **MESS's own pair departs by the same amount** (−7.2 … −13.9% at 2000 K). This is a property of the CSE rate coefficients at poor separation, not an error.
- **TimeIntegration gives the same association:** k_∞ times the slowest-mode amplitude of a pulse equals eq. 28 within 1.9·10⁻⁴ (`reports/four_methods_figures.md`).
- **Deviation from MESS:**
  - dissociation +4.4 … +4.6% at 300 K (tunneling model), +0.1 … +1.1% at 1000 K, +0.3 … +0.4% at 2000 K;
  - association −0.9 … +0.5% at 1250–2000 K.
- **No species are merged** with `ChemicalEigenvalueMax 0.2` (the value of both validation decks). The largest $`\Lambda_N/\Lambda_{N+1}`$ is 0.067 (C₂H₃, 2000 K) and 0.107 (ZZ-allyl + O₂, 330 K, 500 Torr).

**Species merging** (`reports/cse_species_merging.md`):
- tests on a two-well network with a low isomerization barrier: A and B merge into `A+B` when the threshold lies between $`\Lambda_1`$ and $`\Lambda_2`$; $`Q_{A+B} = Q_A + Q_B`$; the merged species decays with $`k_{\mathrm{uni}}`$ of the thermal eigenvector (10⁻⁸);
- with one species and with none (both wells in the bimolecular group), the long-time yields still equal the final steady state (10⁻⁸): the identity holds for any invertible group matrix;
- separated wells are not merged with the default threshold;
- the partition (three tests): wells in equilibrium form one group and a decoupled well joins the bimolecular group; the partition with the largest projection is chosen; a well with a small projection joins the group whose projection it raises.

## 8. Code

- `src/masterequation/chemically_significant_eigenvalues.rs` (`phenomenological_rate_coefficients`, `reactant_yields`, `CseMerging`, `chemical_eigenvalue_count`, `partition_wells`), `chemical_activation_driver.rs` (`run_phenomenological_rates`), `numeric/lapack_interface.rs`, `numeric/dense_inverse.rs`, `report_sections.rs` (`cse_groups`, `write_cse_species_tables`).
- `chemical_activation_from_mess_input.rs`: `bimolecular_high_pressure_rates`, the capture rate coefficients of all bimolecular species.
- `phenomenological_rate_coefficients_with_sources` and `projection_warnings` (source projections), `cse_time_evolution.rs` (`compare_cse_with_time_integration`), `report_sections.rs` (`write_cse_source_projections`, `write_cse_time_comparison`).
- Reports: `reports/chemically_significant_eigenvalues_method.md`, `reports/cse_species_merging.md`, `reports/cse_kappa_and_product_rates.md`, `reports/nonthermal_sources_design.md` (Sections 6, 10.3, 12).

## 9. References

- MK06: J. A. Miller, S. J. Klippenstein, J. Phys. Chem. A 110, 10528 (2006).
- G13: Y. Georgievskii, J. A. Miller, M. P. Burke, S. J. Klippenstein, J. Phys. Chem. A 117, 12146 (2013).
- BW74: J. T. Bartis, B. Widom, J. Chem. Phys. 60, 3474 (1974).
