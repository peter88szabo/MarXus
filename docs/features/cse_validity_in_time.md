# CSE validity in time: when do rate coefficients describe the experiment?

[← README](../../README.md) · prepared experiments: [Prepared experiments](prepared_experiments.md) · [Distributions](preparation_distributions.md) · [Time profiles](source_time_profiles.md) · [Bath history](bath_history.md) · [Transient diagnostics](transient_diagnostics.md) · [CSE source projection](cse_source_projection.md) · **CSE validity in time** · [Fragment wells](fragment_wells_and_lumped_reactants.md)

**Run with:** `Method TimeIntegration`, a [`Preparation` block](prepared_experiments.md) in a constant bath, and `CompareWithCse X` (deck) or `--compare-with-cse X`.

## 1. The question it answers

Rate coefficients describe the kinetics "for all times significantly beyond that characterizing the slowest decaying energy relaxation process" (Miller et al. 2016). How long that takes after a disturbance, how much reacts before, and how wrong a rate-coefficient model is in the meantime is the question of Barker, Frenklach and Golden (2016).

The comparison answers it for the actual preparation. The CSE description of the same experiment is propagated in time and compared with the master equation.

## 2. Equations

**Species and yields:**

```math
\frac{dX}{dt} = 𝐊 X + \sum_a R_a(t)\,x_a, \qquad \frac{dY_\nu}{dt} = \sum_i k_{i\to\nu} X_i + \sum_a R_a(t)\,p_{a,\nu},
```

**Ingredients:**
- **Jumps:** $`N_a x_a`$ and $`N_a p_a`$ at the impulses, and $`N_0 x_0`$, $`N_0 p_0`$ at t = 0.
- **𝐊:** from the phenomenological rate coefficients (G13 eqs. 26, 29).
- **$`x_a`$, $`p_a`$:** the [source projection](cse_source_projection.md) of every source.

**Integrator.** The same Rosenbrock integrator as the master equation, with the same events.

**Comparison per output time:**
- **absolute deviation:** $`\max |X^\mathrm{CSE} - X^\mathrm{ME}| / N_\mathrm{in}(t)`$ over the species (wells summed per species) and the bimolecular channels (exits summed per product);
- **relative deviation:** $`\max |X^\mathrm{CSE} - X^\mathrm{ME}| / \max(|X^\mathrm{CSE}|, |X^\mathrm{ME}|)`$ over the quantities above 10⁻¹⁰ $`N_\mathrm{in}`$.

The relative deviation shows a delayed onset of products, whose absolute size is small.

**Agreement time t\*:** the first output time after which the relative deviation stays ≤ X.

## 3. Settings

| deck keyword | option | meaning |
|---|---|---|
| `CompareWithCse X` | `--compare-with-cse X` | tolerance of the relative deviation; without it, no comparison |
| `EigenSolver`, `ChemicalEigenvalueMax`, `WellProjectionThreshold`, `ChemicalSubspaceCriterion` | as for CSE | the CSE description that is propagated |
| `TimeRange[s]`, `TimesPerDecade`, `IntegrationTolerance` | as for TimeIntegration | the output times, which are also the resolution of t* |

**Requirements:**
- a constant bath (no `Bath` block, or one segment), since the rate coefficients change with T and p;
- the final variant, without an absorbing barrier;
- wells in equilibrium with a bimolecular species have no CSE species and are not compared (noted).

## 4. Output

**Report** (section "CSE description in time"):
- the chemical eigenvalues;
- the lowest relaxation eigenvalue $`\lambda_\mathrm{relax}`$;
- t\*, $`t^*\lambda_\mathrm{relax}`$ (t\* in relaxation times) and the conversion at t\*;
- notes and warnings;
- per output time: the absolute and relative deviation, the species or channel where it is largest, and $`X`$, $`Y`$ of both descriptions.

**CSV:** the block `# CSE description in time …`.

## 5. How to read it

- **Hot pulse:** the CSE description has the prompt products at t = 0. The master equation forms them during the relaxation, so the relative deviation starts near 1 and falls with the relaxation.
- **Cold start (shock):** the CSE yield is negative until about τ_inc; it decays too early. Agreement comes after a few incubation times.
- **Continuous feed:** the relative deviation levels off at the share of molecules in transit through the relaxational modes, about $`R\,c_\lambda/\Lambda_\lambda`$ per mode. That is a physical difference, not an error.

## 6. Validation

**Unit tests** (`cse_time_evolution.rs`):

| case | result |
|---|---|
| hot pulse | the early absolute deviation equals the prompt yield; late agreement below 10⁻⁷ |
| cold start | ln(n)/λ₁ = the incubation time of the time integration (10⁻⁵); agreement only after τ_inc |
| two wells with an impulse and a feed | late relative deviation below 10⁻⁴ |
| refusals | a bath history and the intermediate variant |

**C₂H₃ shock heating** (`validation/c2h3_mess_example/`, Section 4.9): agreement within 1% from t\* = 31.6 ns, which is 2.5 τ_inc or 7.4 relaxation times, at a conversion of 2.8·10⁻⁴.

## 7. Code and references

**Code:**
- `src/masterequation/cse_time_evolution.rs`: `compare_cse_with_time_integration`, `CseTimeComparison`;
- `report_sections.rs`: `write_cse_time_comparison`;
- `solution_method.rs` and `mess_input.rs`: `CompareWithCse`.

**References:**
- J. A. Miller, S. J. Klippenstein, S. H. Robertson, M. J. Pilling, R. Shannon, J. Zádor, A. W. Jasper, C. F. Goldsmith, M. P. Burke, J. Phys. Chem. A 120, 306 (2016).
- J. R. Barker, M. Frenklach, D. M. Golden, J. Phys. Chem. A 120, 313 (2016).
- Y. Georgievskii, J. A. Miller, M. P. Burke, S. J. Klippenstein, J. Phys. Chem. A 117, 12146 (2013) (G13).
