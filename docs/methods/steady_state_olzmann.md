# SteadyStateOlzmann: the final steady state

[← README](../../README.md) · the four methods: **SteadyStateOlzmann** · [SteadyStateAbsorbingBarrier](steady_state_absorbing_barrier.md) · [CSE](chemically_significant_eigenvalues.md) · [TimeIntegration](direct_time_integration.md)

**Family:** steady state. **Run with:** `Method SteadyStateOlzmann` in the `MarXus` block of the deck, or `--method steady-state-olzmann`.

## 1. The question it answers

Reactants form chemically activated adducts continuously. The adducts decompose, isomerize or are stabilized, and the stabilized ones react thermally in turn. Once the stabilized population has itself reached a steady state, there is no net stabilization any more.

**Questions:**
- Where does every formed molecule end in this state?
- How fast does a thermalized adduct react?

GO10 (p. 12295): the final steady state is reached when "there is no more net stabilization; the stabilization reservoir is filled up, and time-independent energy distributions have been established".

## 2. Equations

**The master equation** (PO14 eq. 2; GO10 eq. 7) for the populations $`N_i`$ of the grains of all wells:

```math
\frac{dN}{dt} = R\,F - 𝐉\,N, \qquad 𝐉 = \omega\,(𝐈 - 𝐏) + 𝐊 + k_c[\mathrm D]\,𝐈 .
```

- **Collisions:** $`\omega`$ is the collision frequency and $`𝐏`$ the collisional transition probabilities.
- **Reactions:** $`𝐊`$ holds the microcanonical reactions (product channels and isomerization).
- **Sink:** $`k_c[\mathrm D]`$ is the pseudo-first-order bimolecular sink of a well.
- **Formation:** $`R\,F`$ is the formation of chemically activated adducts. Its source is the normalized nascent distribution $`F(E) \propto \rho(E)\,k_{\to\mathrm{A+B}}(E)\,e^{-E/k_BT}`$ of the entrance channels (PO14 eqs. 7, 9).

**Final steady state:** no absorbing barrier, $`dN/dt = 0`$ (GO10 eq. 8; PO14 eq. 5):

```math
𝐉\,N^{s} = R\,F, \qquad \tilde N^{s} = \frac{𝐉^{-1} F}{\sum_i (𝐉^{-1} F)_i}, \qquad k^{ca}_r = \sum_E k_r(E)\,\tilde N^s(E) \quad \text{(GO10 eq. 9)},
```

```math
\Phi_r = \sum_E k_r(E)\,N^s(E) \quad \text{(O02 eq. 10)}, \qquad \Phi_{\mathrm{sink}} = k_c[\mathrm D] \sum_E N^s(E), \qquad \sum_r \Phi_r + \Phi_{\mathrm{sink}} = 1 .
```

**Thermal rate coefficients of the same $`𝐉`$** (GO10 eq. 12 and the text after it):

```math
𝐉\,\tilde n^{th} = \lambda_1\,\tilde n^{th}, \qquad k^{th}_r = \sum_E k_r(E)\,\tilde n^{th}(E), \qquad k_{\mathrm{uni}} = \sum_r k^{th}_r + k_c[\mathrm D] = \lambda_1 .
```

The eigenvector $`\tilde n^{th}`$ is the decay distribution of thermalized adducts without a source, what GO10 calls the thermal steady-state population. In general it differs from the driven distribution $`\tilde N^s`$.

**Thermal fate of each well:** the final steady state with the Boltzmann distribution $`f^0_w`$ of one well as the source, $`Y_x(w) = k_x^T\,𝐉^{-1} f^0_w`$. It is the probability that a molecule thermalized in well $`w`$ ends in exit $`x`$.

## 3. Algorithm

**Linear systems.** $`𝐒 = 𝐃^{-1}\,𝐉\,𝐃`$ with $`𝐃 = \mathrm{diag}(\sqrt{f^0})`$ is symmetric, by detailed balance. $`𝐒\,y = 𝐃^{-1} F`$ is solved by banded Cholesky factorization, followed by iterative refinement on $`𝐉\,N = F`$. The relative residual must be at most 10⁻⁸.

**Thermal eigenpair:**
- **Default:** inverse iteration with the banded Cholesky factor of $`𝐒 + \sigma𝐈`$, with $`\sigma = n\,\varepsilon\,\max S_{ii}`$. It also works where $`\lambda_1`$ lies below the double-precision floor $`\varepsilon\,\max S_{ii}`$.
- **Alternatives:** the full Householder/QL decomposition (`FullDecomposition`) or LAPACK DSYEVD (`Lapack`).
- **Reported value:** $`k_{\mathrm{uni}}`$ is the eigenvector average. $`\lambda_1`$ is printed beside it, and the sum rule $`|\lambda_1 - k_{\mathrm{uni}}|/k_{\mathrm{uni}}`$ warns above 1.5%, without rejecting.

## 4. Settings

| deck keyword | option | default |
|---|---|---|
| `Method SteadyStateOlzmann` | `--method steady-state-olzmann` | required (no default method) |
| `EigenSolver` | `--eigen-solver inverse\|full\|lapack` | InverseIteration |
| `SumRuleTolerance` | `--sum-rule-tolerance X` | 1.5e-2 |

## 5. Output

**Report:** RUN SETTINGS, CHEMICAL NETWORK, energetics. The tables below come in three views: by temperature, by pressure, and temperature–pressure.
- **Yields:** % of the formed adducts, and without the return to the reactant (% of the net reaction).
- **Chemical-activation rate coefficients** $`k^{ca}`$ (1/s).
- **Bimolecular-to-bimolecular rate coefficients, overall** (chemical activation + thermal reaction of the stabilized adducts): $`k(R \to P) = k_\infty \Phi_P`$ (cm³/s), summed over the channels that form P. The sinks are `R->escape(W)`.
- **Bimolecular-to-bimolecular yields, overall** (% of the net reaction).
- **Capture, return and net reaction of R** (cm³/s).
- **Thermal rate coefficients:** $`k_{\mathrm{uni}}`$, $`\lambda_1`$ and the channels (1/s); thermal branching (%); high-pressure rate coefficients; eigenpair diagnostics.
- **Thermal fate of every well** (%).

There is no bimolecular-to-well (stabilization) rate here: in the final steady state there is no net stabilization. Those rates come from [SteadyStateAbsorbingBarrier](steady_state_absorbing_barrier.md) and [CSE](chemically_significant_eigenvalues.md).

**Machine-readable** (`--csv FILE`): the blocks `# final steady state` and `# thermal rate coefficients of the final steady state`; `FILE_tables.csv` holds every table of the report.

## 6. Validity and limits

**Time scale.** The solution holds when the experimental time is distinctly longer than $`1/k_{\mathrm{uni}}`$ (O02 p. 3618).

**Physical sinks** (e.g. O₂ addition, an escape channel) give the yields of a continuously fed system. A sink in the window $`0.01\,\omega > k_c[\mathrm D] > 10\,k_{\mathrm{uni}}`$ reproduces the absorbing-barrier yield within 10% (O02 p. 3618).

**Without a sink and with a single exit,** every molecule leaves through that exit ($`\Phi = 1`$, O02 after eq. 13). For deep wells at low T, $`𝐉\,N = F`$ is then numerically singular, e.g. C₂H₃ at 300–500 K; the condition is reported as "not available". The thermal eigenpair is still obtained.

## 7. Validation

**ZZ-allyl + O₂, four wells** (`validation/ZZAllyl+O2_Gamma_Case2/`):
- The long-time yields equal CSE (P5, escape, P1 to 6·10⁻⁷) and TimeIntegration (in all printed digits) at all 21 conditions.
- The overall k(R → IEPOX + OH) = 1.943·10⁻¹³ cm³/s at 300 K and 760 Torr, MESS Eckart model.
- With the exact Eckart tunneling, the IEPOX + OH share is 6.6–12.0% above MESS's long-time fate; with the MESS Eckart model it is 5.1–5.9% below.
- The escape share agrees within 0.5%.

**H + C₂H₂ ⇌ C₂H₃** (`validation/c2h3_mess_example/`, Section 4.5):
- **Thermal k_uni against MESS:** +4.7 … +5.7% at 300 K (tunneling model), −2.9 … +1.8% from 750 to 2000 K. It equals CSE's k(W1 → P1) in all printed digits.
- **The three eigen-solvers** agree to 7 digits at 1000 K. Inverse iteration and LAPACK agree in all printed digits at 750–2000 K, and to 5·10⁻⁵ at 300 K, where λ₁ is 13 orders below the double-precision floor.
- **Association by detailed balance, k_uni·K, all 40 conditions:**
  - 300–1000 K: equal to the absorbing-barrier association within 0.02% and to CSE within 0.01%;
  - 300–500 K: +3.2 … +5.8% from MESS (tunneling model). The final steady state itself is singular there, but the thermal eigenpair is not;
  - 750–1750 K: −2.3 … +3.9%;
  - 2000 K: +5.5 … +13.2%. MESS's pair departs from K there by −7 … −14% (CSE's by the same amount), while k_uni·K imposes K.
- **Low-energy reservoir state** at the bottom of the well (`reports/low_energy_reservoir_state.md`): 2 grains at 300 K, 44 at 2000 K. It keeps the thermal population complete, so k_uni is not biased.
- **Thermal k_uni** = the late-time decay rate of [TimeIntegration](direct_time_integration.md) within 3·10⁻⁶ (30 conditions, 750–2000 K). The time integration uses no eigenvector.

## 8. Code

- `src/masterequation/chemical_activation_operator.rs` (𝐉), `chemical_activation_steady_state.rs` (solve), `chemical_activation_observables.rs` (yields, $`k^{ca}`$), `chemical_activation_eigen.rs` (thermal eigenpair), `chemical_activation_driver.rs` (`run_chemical_activation`, `run_thermal_rate_coefficients`, `run_thermal_well_fates`), `report_sections.rs` (report groups).
- Reports: `reports/olzmann_absorbing_barrier_and_treatment.md`, `reports/master_equation_validity_checks.md`, `reports/output_report_and_yields.md`.

## 9. References

- GO10: G. González-García, J. Olzmann, Phys. Chem. Chem. Phys. 12, 12290 (2010).
- O02: J. Olzmann, Phys. Chem. Chem. Phys. 4, 3614 (2002).
- PO14: M. Pfeifle, J. Olzmann, Int. J. Chem. Kinet. 46, 231 (2014).
- PR03: M. J. Pilling, S. H. Robertson, Annu. Rev. Phys. Chem. 54, 245 (2003).
