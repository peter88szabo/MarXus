# Phenomenological rate coefficients from the chemically significant eigenvalues (CSE method)

**MarXus, 2026-10-05.** Peter's request: "make the eigenvalue method of Klippenstein as a new method" (also: "check /home/peter/Programs/pTDME … this is also a CSE code").

## 1. Literature

| abbreviation | reference | used for |
|---|---|---|
| MK06 | J. A. Miller, S. J. Klippenstein, J. Phys. Chem. A 110, 10528 (2006) | the CSE concept: S − 1 chemically significant eigenmodes (eq. 19), separation \|λ_N\| ≪ \|λ_N+1\|, initial-rate and long-time methods (eqs. 24, 25), infinite sinks (eqs. 26–29) |
| G13 | Y. Georgievskii, J. A. Miller, M. P. Burke, S. J. Klippenstein, J. Phys. Chem. A 117, 12146 (2013) | the formulation implemented: bimolecular reactants as decoupled thermal sources (eq. 1), all rate-coefficient expressions (eqs. 15, 21–30) |
| BW74 | J. T. Bartis, B. Widom, J. Chem. Phys. 60, 3474 (1974) | origin of the long-time/eigenvector method (cited by MK06 and G13) |
| NS20 | T. L. Nguyen, J. F. Stanton, J. Phys. Chem. A 124, 2907 (2020) | pTDME (`/home/peter/Programs/pTDME`) |

**How pTDME compares.** pTDME is a fixed-J, time-dependent (E,J) master-equation code:
- it uses ARPACK for the CSEs, or LAPACK for all eigenpairs;
- its outputs are time-dependent species populations ("population.i") and time-dependent flux rate coefficients ("rate.i");
- stabilized intermediates are collected in sinks.

It does not extract the phenomenological rate-coefficient matrix of G13. That makes it a different eigenvalue approach, close in spirit to Olzmann's N(t) expansion. It was not implemented here.

## 2. The method as implemented (`src/masterequation/chemically_significant_eigenvalues.rs`)

**Operator.** J is the relaxation operator of the final steady state (`assemble_operator`, `SteadyState::Final`): all wells, no absorbing barrier, the product channels and escape sinks as losses. Its symmetrized form is S = D⁻¹JD with D = diag(√f⁰), where f⁰ = ρ e^{−E/kT} on the common absolute energy scale (`symmetrize`).

**Eigenpairs.** All eigenpairs (Λ_λ, u_λ) come from a full decomposition: LAPACK DSYEVD, or the in-house Householder/QL. The eigenvectors of G are f^(λ) = D·u_λ, normalized under the scalar product of G13 eq. 12.

**Chemical eigenstates.** For N wells, the N lowest eigenpairs are the chemical eigenstates (G13 Sec. III; MK06 eq. 19).

**Expressions** (G13):

| quantity | expression | G13 eq. |
|---|---|---|
| Q_i | Σ_(grains of i) f⁰ | 17 |
| M_(i,λ) | (1/√Q_i) Σ_(grains of i) f^(λ) | 25 |
| p_λ^(ν) | Σ_(all grains) f^(λ)(E)·k_(→ν)(E) | 15 |
| k_(j→i) | −√(Q_i/Q_j)·(MΛM⁻¹)_(i,j) | 27 |
| k_i (total loss) | (MΛM⁻¹)_(i,i) | 29 |
| k_(i→ν) | (1/√Q_i)·Σ_(λ chem) (M⁻¹)_(λ,i)·p_λ^(ν) | 30 |
| k_(R→i) | (√Q_i/Q_R)·Σ_(λ chem) M_(i,λ)·p_λ^(R) | 28 |
| k_(R→μ) | (1/Q_R)·Σ_(λ relax) p_λ^(μ)·p_λ^(R)/Λ_λ | 21 |
| k_(R→R) | k^(c) − Σ_μ k_(R→μ) − Σ_i k_(R→i) | 22 |

**Bimolecular channels ν.** They are the product names of the channels, the reactant included, plus one "escape(W)" for every well with a pseudo-first-order sink.

**Normalization of the reactant rates.** The reactant enters only through 1/Q_R = k^(c) / Σ_(all grains) k_(→R)·f⁰ (G13 eq. 23). k^(c) is the capture (high-pressure association) rate coefficient of the entrance channels, already computed by the deck adapter (`EntranceHighPressureRate`), so no partition function of the reactants is needed.

**Dense inverse.** M⁻¹ is computed by Gauss–Jordan elimination with partial pivoting (`numeric/dense_inverse.rs`; Numerical Recipes in Fortran, 2nd ed. (1992), §2.1).

**Diagnostics, printed with every table:**
- the chemical eigenvalues;
- the lowest relaxation eigenvalue and the separation Λ_N/Λ_(N+1). A warning is given above 0.1, because MK06 requires \|λ_N\| ≪ \|λ_N+1\| and G13 Sec. IV says to merge wells otherwise. Merging is not implemented.
- the relaxational projection 1 − Σ_i M_(i,λ)² of every chemical eigenvector (G13 Fig. 2);
- the loss balance of eq. 29;
- the largest detailed-balance deviation k_(i→j)Q_i vs k_(j→i)Q_j;
- the double-precision floor ε·max S_ii, with a warning when a chemical eigenvalue is within a factor 100 of it.

**Use.** `cargo run --release --example chemical_activation_from_deck -- deck.inp R --method cse`, or `Method CSE` in the `MarXus` block of the deck header (`solution_methods_and_deck_settings.md`). CSE is the second solution method of MarXus, beside the steady state:
- **Eigen-solver:** LAPACK unless `--eigen-solver full`; `--eigen-solver inverse` is refused.
- **Output:** MESS-style species tables per (T, p): rows from, columns to, wells in s⁻¹, the reactant row in cm³ s⁻¹. The diagonal holds the total loss of each well, and for the reactant the net reaction (capture minus return), as in MESS.
- **Explanation:** comment lines with the references precede the tables.

## 3. Tests

| test | checks |
|---|---|
| `chemically_significant_eigenvalues::a_single_well_gives_the_eigenvector_average_as_its_rate_coefficient` | one chemical eigenstate: k_(W→P) = Olzmann's eigenvector average k_uni (GO10 after eq. 12) to 10⁻⁸; k_W = Λ₁ |
| `…::well_rates_satisfy_the_loss_balance_and_detailed_balance` | two wells, products, escape sink: eq. 29 to 10⁻⁸, which is an identity because the column sums of J are the losses; detailed balance of the isomerization to 10⁻⁴; separation |
| `…::at_high_pressure_the_rates_approach_the_high_pressure_limits` | 10⁹ Torr: k_(A→B) and k_(A→P) equal the Boltzmann-averaged k(E) to 10⁻³ |
| `…::reactant_rates_satisfy_detailed_balance_with_the_reverse_dissociation` | k_(R→i)/k_(i→R) = Q_i/Q_R: 10⁻³ for the directly connected well (3·10⁻⁵ found); 10⁻² for the indirect well (0.2% found, see below); capture balance; non-negative return |
| `numeric::dense_inverse::*` (3) | inverse × matrix = I; pivoting; singular and non-square refused |
| `chemical_activation_driver::the_cse_route_writes_a_species_table_per_condition` | the table writer |

**Detailed balance of indirect pairs.** G13 Sec. IV states that their expressions satisfy detailed balance to the extent that M is orthogonal. For a well reached from R only through another well, k_(R→B) is a small difference of the two chemical modes, which amplifies the non-orthogonality (about 3·10⁻⁵ here) to 0.2%.

## 4. Validation: ZZ-allyl + O₂, Gamma Case 2 (`validation/ZZAllyl+O2_Gamma_Case2/`)

**Setup.** MarXus CSE (`--method cse --tunneling mess-eckart`, deck with `TSTLevel E`) against Peter's MESS run: 21 conditions, 270–330 K, 500–760 Torr. Results are in `cse_comparison.csv` and `plots/cse_vs_mess.png`.

**Example: 270 K, 500 Torr.**

| entry | MESS | MarXus CSE |
|---|---|---|
| G2 → G3 / G4 / R / P5 | 3.612e4 / 1.004e6 / 2201 / 166.3 s⁻¹ | 3.645e4 / 1.010e6 / 2251 / 162.0 s⁻¹ |
| G4 → G2 / P5 / escape | 3182 / 71.9 / 2.486e7 s⁻¹ | 3213 / 72.5 / 2.500e7 s⁻¹ |
| R → G2 / G3 / G4 | 8.005e-12 / 8.14e-14 / 3.052e-12 cm³/s | 8.207e-12 / 8.24e-14 / 3.110e-12 cm³/s |
| R → P5 (IEPOX + OH) | 2.483e-13 cm³/s | 2.406e-13 cm³/s |
| R → escape(G4) | −1.32e-13 cm³/s | −1.47e-13 cm³/s |

**All 21 conditions, significant entries:**
- R → G2: +2.1 … +3.3%; R → G3: +0.9 … +3.6%; R → G4: +1.4 … +6.1%.
- **R → P5: −2.8 … −3.5%**; R → P1: −1.4 … −2.3%; R → P7: −3.9 … −5.0%.
- Well → well: within ±1%, except G4 → G3, +1.0 … +2.5%.
- Well → R: +0.2 … +6%.
- Well → products: within −3.4 … +3.8%.
- G4 → escape: +0.6%.
- Entries of 10⁻⁹ to 10⁻²³ (from G6, to G6, and to escape from G2, G3, G6) are rounding noise in both codes, and are partly negative in MESS.
- **The negative R → escape(G4) entry of MESS is reproduced** (MarXus +10 … +37% on this small mixing term). It is therefore a property of the G13 bimolecular-to-bimolecular expression (eq. 21) under these conditions, not a MESS artefact.

**Separation.** The escape mode of G4 (2.5·10⁷ s⁻¹) has Λ_N/Λ_(N+1) ≈ 0.10–0.11. A warning is printed; MESS's deck allows chemical eigenvalues up to 0.2 (`ChemicalEigenvalueMax`).

**Remaining differences of a few percent.** These are consistent with the collisional part. Every k∞ agrees within 1.7% with the MESS Eckart model, and the residual of the long-time IEPOX + OH share was attributed to the collision kernel normalization, the collision frequency or the graining (README Section 4.4).

## 5. Open points

1. **Well merging** (G13 Sec. IV, eq. 31) when the chemical eigenvalues approach the relaxation ones: not implemented; only a warning.
2. **Product → well rate coefficients** (bimolecular-to-isomer for the products). They need the products' partition functions (Q_ν), which the adapter does not yet provide for non-reactant bimolecular species.
3. **Precision.** Chemical eigenvalues near the double-precision floor are only warned about. This is the planned higher-precision work (`higher_precision_decision.md`).
