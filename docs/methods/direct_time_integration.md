# TimeIntegration: direct time integration of the master equation

[← README](../../README.md) · the four methods: [SteadyStateOlzmann](steady_state_olzmann.md) · [SteadyStateAbsorbingBarrier](steady_state_absorbing_barrier.md) · [CSE](chemically_significant_eigenvalues.md) · **TimeIntegration**

**Family:** time integration. **Run with:** `Method TimeIntegration`, or `--method time-integration`.

## 1. The question it answers

How do the populations of all grains and the yields of all exits evolve in time, through the early relaxation and the later chemistry? The start is either a pulse of chemically activated adducts or continuous formation.

No steady state is assumed, and no separation of chemical and relaxation modes is needed.

## 2. Equations

The populations $`N`$ of all grains (operator of the final steady state, no absorbing barrier) and the yield $`Y_x`$ accumulated in every exit $`x`$ (product channels and escape sinks):

```math
\frac{dN}{dt} = R\,F - 𝐉\,N, \qquad \frac{dY_x}{dt} = \sum_E k_x(E)\,N(E) .
```

**Initial states:**
- **Pulse:** $`N(0) = F`$, $`R = 0`$.
- **Continuous formation:** $`N(0) = 0`$, $`R = 1\ \mathrm{s^{-1}}`$.

The column sums of $`𝐉`$ are the losses out of the network. For a pulse, therefore, $`\sum N + \sum Y = 1`$ at all times.

**Long-time limit of a pulse:**

```math
Y_x(\infty) = \int_0^\infty k_x^T\,e^{-𝐉 t} F\,dt = k_x^T\,𝐉^{-1} F ,
```

the yield of the final steady state ([SteadyStateOlzmann](steady_state_olzmann.md)) with the same $`F`$.

## 3. Algorithm

**Rosenbrock methods** (`src/numeric/integrators/`). Adapted from KPP, the Kinetic PreProcessor (int/rosenbrock.f90; (C) A. Sandu 2004, revised by P. Miehe and A. Sandu 2006; GPL-3.0, as MarXus). One step from $`t`$ to $`t + h`$ (Hairer, Wanner, Sec. IV.7):

```math
𝐆 = \frac{1}{h\gamma_1}\,𝐈 - \frac{\partial f}{\partial y}, \qquad 𝐆\,K_i = f\Big(t + \alpha_i h,\ y + \sum_{j<i} A_{ij} K_j\Big) + \sum_{j<i} \frac{C_{ij}}{h}\,K_j, \qquad y_1 = y + \sum_i M_i K_i ,
```

- **Error estimate:** $`\sum_i E_i K_i`$ in the scaled RMS norm.
- **Step control:** $`h_{\mathrm{new}} = h \min\big(6, \max(0.2,\ 0.9/\mathrm{err}^{1/\mathrm{ELO}})\big)`$.

| method | stages | order | stability | reference |
|---|---|---|---|---|
| Rodas4 (default) | 6 | 4 | stiffly accurate | Hairer, Wanner (1991, 1996) |
| Rodas3 | 4 | 3 | stiffly accurate | Sandu et al. (1997) |
| Ros4 | 4 | 4 | L-stable | Hairer, Wanner (1991, 1996) |
| Ros3 | 3 | 3 | L-stable | Sandu et al. (1997) |
| Ros2 | 2 | 2 | L-stable | Verwer et al. (1999) |

**Not adopted: KPP's Rang3.** Its embedded error estimate vanishes for linear problems, so the step control accepts any step (`reports/direct_time_integration.md`).

**Stage solve for the master equation.** The Jacobian is constant: $`-𝐉`$ for $`N`$, $`k_x^T`$ for $`Y`$. $`𝐆`$ is block lower triangular:

```math
(s\,𝐈 + 𝐉)\,x_N = b_N \;\Rightarrow\; x_N = 𝐃\,(s\,𝐈 + 𝐒)^{-1}\,𝐃^{-1} b_N, \qquad x_Y = \frac{b_Y + K^T x_N}{s}, \qquad s = \frac{1}{h\gamma_1} .
```

**Factorization reuse.** The banded Cholesky factors of $`s\,𝐈 + 𝐒`$ are cached. Steps are rounded down to powers of two (a MarXus option of the integrator, off in plain KPP), so step sizes repeat and factorizations are reused. Case 2: 120 factorizations for 895 steps.

## 4. Settings

| deck keyword | option | default |
|---|---|---|
| `Method TimeIntegration` | `--method time-integration` | required (no default method) |
| `Integrator` | `--integrator rodas4\|rodas3\|ros4\|ros3\|ros2` | Rodas4 |
| `InitialState` | `--initial pulse\|continuous` | Pulse |
| `TimeRange[s]` | `--time-range T1 T2` | 1e-12 1e2 |
| `TimesPerDecade` | `--times-per-decade N` | 4 |
| `IntegrationTolerance` | `--integration-tolerance X` | 1e-6 (relative; absolute 1e-14) |

## 5. Output

**Per condition:** a table of time against the population of every well and the yield of every exit (% of the formed adducts), with the integrator work (steps, rejected steps, function evaluations, factorizations).

**At the last output time,** in three views (by temperature, by pressure, temperature–pressure):
- the yields;
- the populations left;
- the yields without the return to the reactant;
- the **bimolecular-to-bimolecular rate coefficients, overall**, $`k(R \to P) = k_\infty Y_P`$, with their **yields**;
- the conserved total.

**Machine-readable:** blocks `# time evolution: T = … K, p = … Torr` with t, N(W), exit yields; `FILE_tables.csv`.

## 6. Validity and limits

The method has no physical restriction beyond the master equation itself. The output time range must cover the slowest process of interest. Thermal decay of deep wells at low T (k_uni ≪ 1 s⁻¹) is slower than any practical range; at the last time those populations are reported as left in the wells.

## 7. Validation

**Unit tests (8, integrator):**
- exponential decay;
- Prothero–Robinson (λ = 10⁶, non-autonomous);
- Robertson (reference values at t = 40 to 10⁻⁵, linear invariant to 10⁻¹²);
- observed orders (Ros2 1.9, Ros3 3.0, Ros4 3.9, Rodas3 2.9, Rodas4 4.5);
- power-of-two steps;
- tolerance checks.

**Unit tests (6, master equation):**
- **pulse:** total = 1 to 10⁻⁹, long-time yields = final steady state to 10⁻⁶;
- **continuous formation:** approaches $`N^s`$ to 10⁻⁶;
- **eigenvector expansion:** $`N(t)`$ equals PO14 eqs. 3–4 to 10⁻⁶;
- **thermalized well:** late decay rate = $`\lambda_1`$ to 10⁻⁴;
- **absorbing-barrier operator:** a pulse ends in the stabilization yields of the intermediate steady state.

**ZZ-allyl + O₂, four wells** (pulse, Rodas4, 10⁻¹² … 10² s, MESS Eckart model, 21 conditions; `validation/ZZAllyl+O2_Gamma_Case2/`):
- At t = 100 s every yield equals SteadyStateOlzmann to all 7 printed digits; the total stays 100.0%.
- Time scales at 300 K, 760 Torr (`plots/time_evolution_300K_760torr.png`):
  - redissociation to R within 10⁻⁹ s (78.4%);
  - IEPOX + OH from chemically activated G4 by 10⁻⁸ s (0.502%);
  - escape from G4 between 10⁻¹⁰ and 10⁻⁵ s (21.1%);
  - a thermal P7 stage near 10⁻⁴ … 10⁻³ s.
- The integrator needed 851–938 steps per condition, none rejected, and 117–132 factorizations (111 s for 21 conditions on 4 cores).

**H + C₂H₂ ⇌ C₂H₃, one well** (`validation/c2h3_mess_example/`, §4.4; 40 conditions):
- **Late-time decay.** The decay rate of the pulse, −d ln N/dt from the last two output times with 10⁻⁸ < N < 10⁻³, equals SteadyStateOlzmann's thermal k_uni within 3.1·10⁻⁶. That holds at all 30 conditions that decay inside the window, 750–2000 K (`time_integration_decay_vs_k_uni.csv`). The time integration computes no eigenvector.
- **Association from the slowest mode.** The amplitude A of the late single-exponential decay, extrapolated to t = 0, times k_∞ equals the CSE association k(R → W) (G13 eq. 28) within 4.7·10⁻⁵ at all 40 conditions. For one well both are (Σ f⁽¹⁾)(Σ f⁽¹⁾ k_R)/Σ k_R f⁰ (`reports/four_methods_figures.md`, Section 2).
- **Conservation.** Populations + yields stay 100% in all printed digits at 750–2000 K. At 300–500 K, where λ₁ lies below the double-precision floor, the total drifts by up to 6·10⁻⁴ over 100 s (`reports/method_comparison.md`).
- **At 300 K and 1 atm,** 12.9% redissociates within about 10⁻⁹ s. The remaining 87.1% stays as C₂H₃ until 100 s (k_uni = 8.7·10⁻¹⁶ s⁻¹).
- **Integrator work:** 620–700 steps per condition and about 115–122 factorizations. The full deck took 11 min on 4 cores: the collision band reaches 716 grains at 2000 K.

## 8. Code

- `src/numeric/integrators/rosenbrock_methods.rs` (coefficients), `src/numeric/integrators/rosenbrock.rs` (`integrate`, `StiffSystem`), `src/masterequation/direct_time_integration.rs` (`integrate_master_equation`), `report_sections.rs` (`write_time_evolution_tables`, `time_integration_groups`).
- Report: `reports/direct_time_integration.md`.

## 9. References

- E. Hairer, G. Wanner, Solving Ordinary Differential Equations II, Springer (1991; 2nd ed. 1996), Sec. IV.7.
- A. Sandu, J. G. Verwer, J. G. Blom, E. J. Spee, G. R. Carmichael, F. A. Potra, Atmos. Environ. 31, 3459 (1997).
- J. G. Verwer, E. J. Spee, J. G. Blom, W. Hundsdorfer, SIAM J. Sci. Comput. 20, 1456 (1999).
- KPP, the Kinetic PreProcessor (GPL-3.0).
- PO14: M. Pfeifle, J. Olzmann, Int. J. Chem. Kinet. 46, 231 (2014), eqs. 3–4.
