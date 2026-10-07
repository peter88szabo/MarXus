# Transient diagnostics: flux coefficients, relaxation and incubation

[← README](../../README.md) · prepared experiments: [Prepared experiments](prepared_experiments.md) · [Distributions](preparation_distributions.md) · [Time profiles](source_time_profiles.md) · [Bath history](bath_history.md) · **Transient diagnostics** · [CSE source projection](cse_source_projection.md) · [CSE validity in time](cse_validity_in_time.md) · [Fragment wells](fragment_wells_and_lumped_reactants.md)

**Method:** TimeIntegration with a [`Preparation` block](prepared_experiments.md). **Output:** the section "PREPARED EXPERIMENT: DIRECT TIME INTEGRATION" and the CSV block `# prepared time integration: …` (`--csv`).

## 1. The question it answers

What does an experiment observe, and when is the system in a steady energy distribution? The diagnostics are those used in the literature on non-steady-state energy distributions and on shock-tube incubation and relaxation:
- Barker, Frenklach, Golden 2015;
- Barker, King 1995;
- Eng et al. 2001.

## 2. Quantities

One table per group, with the output times as rows.

| group | quantity | definition | purpose |
|---|---|---|---|
| populations and yields | $`C_w`$, $`Y_x`$ | sum over the grains of well w; cumulative yield of exit x | the measured concentrations and product yields |
| fluxes | $`q_x`$ | $`k_x^T n`$ (1/s) | product formation rates, e.g. a detector signal |
| injected amounts | $`N_{a,\mathrm{in}}`$ | per source | the denominator of yields |
| energies | ⟨E⟩, tail | mean energy above the well bottom; fraction above the lowest threshold | distance from the bath distribution; the reactive tail |
| hazard and balance | $`k_\mathrm{inst}`$, balance | $`\sum q/\sum C`$; $`(\sum n + \sum Y - N_\mathrm{in})/N_\mathrm{in}`$ | overall decay rate; numerical check (round-off level) |
| flux coefficients | $`r_{wc}`$, $`R_w`$ | $`r_{wc} = \sum_i k_{wc}(E_i)\,n_{wi}/C_w`$ for every channel, isomerization included; $`R_w = \sum_c r_{wc} + k_c[\mathrm D]_w`$ (BFG eqs. A5, 1) | time-dependent rate coefficients: they start at the high-pressure (Boltzmann) value for a thermal start and end at the thermal rate coefficient of the final steady state |
| effective coefficients | $`k^e_{wc}(t_1, t)`$ | $`(t - t_1)^{-1}\int_{t_1}^{t} r_{wc}\,dt'`$, $`t_1`$ the first output time (BFG eq. A7) | the rate an experiment would deduce over a time window |
| energy spread and relaxation | σ_E, d⟨E⟩/dt, E_f, τ_vib | standard deviation of E; d⟨E⟩/dt exact from $`\dot n = -𝐉n + \sum_a R_a F_a`$; E_f the mean energy of the lowest eigenvector of 𝐉 (the final steady state, depleted at high E by reaction); $`\tau_\mathrm{vib} = -(\langle E\rangle - E_f)/(d\langle E\rangle/dt)`$ (Barker, King eq. 11) | vibrational relaxation (laser schlieren measures d⟨E⟩/dt) |
| incubation and collisions | τ_inc, Z_w | $`\tau_\mathrm{inc}(t) = (t - s) + \ln(N/N_\mathrm{ref})/k_\mathrm{inst}`$ from the start s of the bath segment, with $`N_\mathrm{ref}`$ the population after s plus the impulses since; $`Z_w = \int_0^t \omega_w\,dt'`$ | the delay of the reaction after a shock (Barker, King eq. 9), in s and in collisions (Eng et al.: $`Z_\mathrm{LJ}[M]\,\Delta t_\mathrm{inc}`$) |

## 3. How to read them

**τ_inc:**
- **Plateau.** τ_inc(t) reaches a plateau once the decay is first order: $`N = N_\mathrm{ref}\,e^{-k(t - s - \tau_\mathrm{inc})}`$, the back-extrapolation of Barker and King. Read it where $`k_\mathrm{uni}\,t \approx 1`$.
- **Late times.** At $`t \gg \tau_\mathrm{inc}`$ it is a small difference of large numbers, with a relative error of about $`t/\tau_\mathrm{inc}`$ times that of $`k_\mathrm{inst}`$.
- **Not defined** (`***`) while a continuous source injects in the segment.

**τ_vib:**
- **Sign.** It is negative while ⟨E⟩ moves away from E_f, e.g. during an injection.
- **Not given** (`***`) once |⟨E⟩ − E_f| ≤ 10 × the relative integration tolerance × E_f. There both are rounding and integration error.
- **E_f** is not available in the intermediate (absorbing-barrier) variant.

**k^e** is integrated over the output times, with r linear in ln t between them. It is exact for a constant r and second order otherwise, so its accuracy is set by `TimesPerDecade`.

## 4. Settings

| deck keyword | option | meaning (default) |
|---|---|---|
| `TimeRange[s]` | `--time-range T1 T2` | first and last output time (1e-12 … 1e2 s) |
| `TimesPerDecade` | `--times-per-decade N` | output density (4); also the accuracy of k^e |
| `IntegrationTolerance` | `--integration-tolerance X` | relative tolerance (1e-6); also the resolution of τ_vib |

## 5. Validation

**Unit tests** (`prepared_time_integration.rs`):
- $`\sum_w C_w(\sum_\mathrm{products} r + k_c[\mathrm D]) = \sum_x q_x`$ at every time;
- r from a Boltzmann start equals the high-pressure value at 10⁻¹⁴ s (to 10⁻³) and the eigenvector average at late times (to 10⁻⁵);
- E_f equals the thermal-eigenvector mean, and the Boltzmann mean for a closed well;
- d⟨E⟩/dt equals a central difference (10⁻⁵), also with a source;
- in a closed well τ_vib → 1/λ₂, the slowest relaxation time;
- the τ_inc plateau equals $`\ln(a_1)/\lambda_1`$ from a full decomposition (10⁻⁴);
- collision numbers add up over segments;
- k^e is exact for a constant rate and converges at second order.

**C₂H₃ shock heating** (300 K population into 1000 K, 1 atm; `validation/c2h3_mess_example/`, Section 4.9):
- τ_inc = 12.47 ns = 62.6 collisions;
- τ_vib = 4.08 ns at τ_inc;
- the late flux coefficient equals k_uni of SteadyStateOlzmann to six digits.

**Grain convergence** of τ_inc: 12.78, 12.51, 12.47, 12.46 ns for grains of 278, 139, 70, 35 cm⁻¹.
- The states are counted exactly on 1 cm⁻¹ cells and summed into grains, so coarse grains do not smooth the low-energy density of states.
- Eng et al. found that such smoothing shortens incubation times drastically.

## 6. Code

- `src/masterequation/prepared_time_integration.rs`: `TransientPoint`, `effective_coefficient`.
- `chemical_activation_eigen.rs`: `lowest_eigenvector_grain_populations`.
- `report_sections.rs`: `write_transient_tables`, `channel_labels`.

## 7. References

- J. R. Barker, M. Frenklach, D. M. Golden, J. Phys. Chem. A 119, 7451 (2015), eqs. 1, A5, A7.
- J. R. Barker, K. D. King, J. Chem. Phys. 103, 4953 (1995), eqs. 9, 11–12.
- C. Eng, A. Gebert, E. Goos, H. Hippler, C. Kachiani, Phys. Chem. Chem. Phys. 3, 2258 (2001).
- J. H. Kiefer, G. C. Buzyna, A. Dib, K. P. Sundaram, J. Chem. Phys. 113, 48 (2000).
