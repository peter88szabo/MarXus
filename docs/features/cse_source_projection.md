# CSE source projection: a prepared distribution on the phenomenological description

[← README](../../README.md) · prepared experiments: [Prepared experiments](prepared_experiments.md) · [Distributions](preparation_distributions.md) · [Time profiles](source_time_profiles.md) · [Bath history](bath_history.md) · [Transient diagnostics](transient_diagnostics.md) · **CSE source projection** · [CSE validity in time](cse_validity_in_time.md) · [Fragment wells](fragment_wells_and_lumped_reactants.md)

**Method:** [CSE](../methods/chemically_significant_eigenvalues.md) with a [`Preparation` block](prepared_experiments.md). **Output:** the section "source projections" of the CSE report.

## 1. The question it answers

A kinetic model with phenomenological rate coefficients needs initial concentrations of its species. For a prepared, non-thermal distribution two questions follow:
- **Species:** how much of each species is present once the internal-energy relaxation is over?
- **Prompt products:** how much has already reacted during that relaxation?

A thermal entrance rate coefficient cannot simply be reused for a laser, photolysis or beam source.

## 2. Equations

**Decomposition.** The pulse F is decomposed on the eigenvectors of the symmetrized operator (Georgievskii et al. 2013, eqs. 12–14):

```math
c_\lambda = \sum_k u_{\lambda k}\,\frac{F_k}{d_k},
```

evaluated in logarithms, since d spans many orders of magnitude.

**The chemical modes** give the species populations after the relaxation (eqs. 24, 37):

```math
n_g = \sqrt{Q_g}\sum_{\lambda\ \mathrm{chem}} M_{g\lambda}\,c_\lambda .
```

**The relaxational modes** give the prompt yields formed during the relaxation (eqs. 38, 41):

```math
Y^\mathrm{prompt}_\nu = \sum_{\lambda\ \mathrm{relax}} \frac{p^{(\nu)}_\lambda\,c_\lambda}{\Lambda_\lambda} .
```

**Identities:**
- **Long-time yields:** $`Y^\mathrm{prompt}_x + \sum_g n_g\,(\text{fate of } g \text{ in } x) = k_x^T 𝐉^{-1} F`$, the long-time yield of the pulse.
- **Mass:** $`\sum_g n_g + \sum_\nu Y^\mathrm{prompt}_\nu = \sum F`$, exact with all eigenpairs (without a bimolecular group).

## 3. How to read it

**Hot preparation** (above the thresholds): positive prompt yields, which are products formed before the molecules relax.

**Cold preparation** (colder than the steady state, e.g. before a shock): a negative prompt yield and species populations above the injected amount.
- **Why it is negative:** the species of the CSE description decay from t = 0, while the master equation first has to activate the population.
- **Meaning:** this is physical, the incubation. For one well, $`n = e^{\lambda_1 \tau_\mathrm{inc}}`$, so $`\ln(n)/\lambda_1`$ is the incubation time of the time integration.

**Warning.** If the mass identity is violated by more than 10⁻⁶ (relative), the eigenpairs are not resolved in double precision. Frankcombe and Smith (2003) describe this failure of eigenvector routes. The direct time integration should then be used.

## 4. Settings and output

- **Settings:** every source of the `Preparation` block is projected at every condition (T, p) of the deck. The CSE settings apply: `EigenSolver` (Lapack or FullDecomposition), `ChemicalEigenvalueMax`, `WellProjectionThreshold`, `ChemicalSubspaceCriterion`.
- **Output:** one table per condition, with one row per source. The columns are the species populations n(species), the prompt yields of every bimolecular channel, and their sum.

## 5. Validation

**Unit tests** (`chemically_significant_eigenvalues.rs`):
- the projection of the reactant's thermal entrance shape gives the reactant rows per capture, to 10⁻⁹;
- prompt + species × thermal fates = $`k_x^T 𝐉^{-1} F`$ of the final steady state for a hot source, to 10⁻⁷;
- the mass check.

**C₂H₃ shock heating** (`validation/c2h3_mess_example/`, Section 4.9): n = 1.000181 and prompt −1.81·10⁻⁴. ln(n)/k_uni = 12.47 ns, the incubation time of the time integration to five digits.

## 6. Code and references

**Code:**
- `src/masterequation/chemically_significant_eigenvalues.rs`: `phenomenological_rate_coefficients_with_sources`, `SourceProjection`, `projection_warnings`;
- `chemical_activation_driver.rs`: `run_phenomenological_rates_with_sources`;
- `report_sections.rs`: `write_cse_source_projections`.

**References:**
- Y. Georgievskii, J. A. Miller, M. P. Burke, S. J. Klippenstein, J. Phys. Chem. A 117, 12146 (2013).
- T. J. Frankcombe, S. C. Smith, J. Theor. Comput. Chem. 2, 179 (2003).
