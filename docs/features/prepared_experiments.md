# Prepared experiments: initial populations, sources and bath histories

[← README](../../README.md) · prepared experiments: **Prepared experiments** · [Distributions](preparation_distributions.md) · [Time profiles](source_time_profiles.md) · [Bath history](bath_history.md) · [Transient diagnostics](transient_diagnostics.md) · [CSE source projection](cse_source_projection.md) · [CSE validity in time](cse_validity_in_time.md) · [Fragment wells](fragment_wells_and_lumped_reactants.md)

**Deck:** a top-level `Preparation ... End` block (a MarXus extension of the MESS format). **Methods:** all four; the time integration uses all of it.

## 1. The question it answers

What happens to molecules that are **not** formed by a thermal entrance flux in a constant bath? Experiments prepare them in other ways:
- **Laser excitation:** a narrow energy band, at a given time.
- **Photolysis of a precursor:** a pulse of finite length, or formation with the precursor's decay.
- **Shock heating:** a cold population put into a hot bath.
- **Beams and flows:** continuous formation.

Such molecules start far from the bath distribution and react before or while they relax. The `Preparation` block describes the experiment in four independent parts:

| part | question | deck | page |
|---|---|---|---|
| initial population | which molecules are present at t = 0, with which energy distribution? | `InitialPopulation` | [Distributions](preparation_distributions.md) |
| source channels | which molecules are added later, with which distribution and time profile? | `Source <name>`, any number | [Distributions](preparation_distributions.md), [Time profiles](source_time_profiles.md) |
| bath history | what are T and p, and when do they change? | `Bath` | [Bath history](bath_history.md) |
| observation | what is measured, and when does a rate-coefficient description apply? | output tables; `CompareWithCse` | [Transient diagnostics](transient_diagnostics.md), [CSE validity in time](cse_validity_in_time.md) |

## 2. Equations

```math
\frac{dn}{dt} = -𝐉[T(t), p(t)]\,n + \sum_a R_a(t)\,F_a, \qquad n(0) = N_0 F_0, \qquad Y_x(t) = \int_0^t k_x^T n\,dt' .
```

**Notation:**
- $`F_0`$ and $`F_a`$ are energy distributions over the grains of all wells, each normalized to 1;
- $`N_0`$ is the initial amount and $`R_a(t)`$ the rate of channel $`a`$;
- $`Y_x`$ is the yield of exit $`x`$ (product channel or escape sink).

**Balance,** checked at every output time: $`\sum n + \sum Y = N_0 + \sum_a N_{a,\mathrm{in}}(t)`$.

## 3. The deck

```
Preparation
  InitialPopulation
    Amount                          1.0
    Distribution Thermal
      Well                          W1
      PreparationTemperature[K]     300
    End
  End
  Source laser
    Profile Impulse
      Time[s]                       1e-9
      Amount                        0.5
    End
    Distribution Gaussian
      Well                          W1
      Centre[kcal/mol]              45.0
      Width[1/cm]                   300
    End
  End
  Source entrance
    Profile Feed
      Start[s]                      0
      Rate[1/s]                     1e3
    End
    Distribution ThermalEntrance
    End
  End
  Bath
    Segment Start[s] 0      Temperature[K] 300    Pressure[atm] 1
    Segment Start[s] 1e-6   Temperature[K] 1000   Pressure[atm] 1
  End
End
```

**Rules:**
- **Optional parts:** all three are optional, but something must be injected.
- **Amounts:** `Amount` of the initial population defaults to 1. Amounts are relative numbers, and yields are in the same units.
- **Reactant:** with a `Preparation` block the deck needs no `Reactant` with an entrance barrier, unless `ThermalEntrance` is used.

## 4. What each method does with it

| method | uses | gives |
|---|---|---|
| [TimeIntegration](../methods/direct_time_integration.md) | everything: impulses as exact jumps, profile edges and bath changes as events of the integrator | the experiment in time ([Transient diagnostics](transient_diagnostics.md)); optionally the CSE description in time ([CSE validity in time](cse_validity_in_time.md)) |
| [SteadyStateOlzmann](../methods/steady_state_olzmann.md), [SteadyStateAbsorbingBarrier](../methods/steady_state_absorbing_barrier.md) | one source shape: the rate-weighted shape of the open-ended feeds, otherwise the amount-weighted shape of everything injected; the deck's conditions | steady-state yields; equal to the infinite-time yields of the whole preparation in a constant bath |
| [CSE](../methods/chemically_significant_eigenvalues.md) | every source distribution; the deck's conditions | the species populations after relaxation and the prompt yields of every source ([CSE source projection](cse_source_projection.md)) |

## 5. Validation

- **`examples/c2h3_prepared_experiment.inp`** (C₂H₃: 300 K population, laser pulse, entrance feed, shock 300 → 1000 K at 1 µs):
  - population balance ≤ 5·10⁻¹² at every output time;
  - 86% of the laser-excited molecules dissociate promptly;
  - after the shock, the mean energy and the decay rate are those of the 1000 K bath;
  - at late times the feed equals the loss.
- **`validation/c2h3_mess_example/`, Section 4.9:** shock heating; incubation, relaxation, CSE validity and grain convergence.
- **Unit tests** (`src/masterequation/prepared_time_integration.rs`):
  - an initial population equals an impulse at t = 0 and the pulse of the original driver (10⁻¹²);
  - a feed from t = 0 equals continuous formation;
  - superposition of channels;
  - a narrowing rectangular pulse converges to the impulse;
  - balance;
  - absorbing barrier;
  - closed-well relaxation to the bath distribution;
  - bath segments;
  - the two-grain analytical solution;
  - cold against hot preparations.

## 6. Code

- **Distributions:** `src/masterequation/prepared_distributions.rs`.
- **Time profiles:** `source_profiles.rs`.
- **Integration and diagnostics:** `prepared_time_integration.rs`.
- **Deck:** `preparation_input.rs`.
- **CSE in time:** `cse_time_evolution.rs`.
- **Projections:** `chemically_significant_eigenvalues.rs`.
- **Report sections:** `report_sections.rs` (`write_preparation_summary`, `write_transient_tables`, `write_cse_source_projections`, `write_cse_time_comparison`).

## 7. References

- J. R. Barker, M. Frenklach, D. M. Golden, J. Phys. Chem. A 119, 7451 (2015); reply in J. Phys. Chem. A 120, 313 (2016).
- J. A. Miller, S. J. Klippenstein, S. H. Robertson, M. J. Pilling, R. Shannon, J. Zádor, A. W. Jasper, C. F. Goldsmith, M. P. Burke, J. Phys. Chem. A 120, 306 (2016).
- S. H. Robertson, "Foundations of the master equation", in Comprehensive Chemical Kinetics 43 (2019), Ch. 5.
- Y. Georgievskii, J. A. Miller, M. P. Burke, S. J. Klippenstein, J. Phys. Chem. A 117, 12146 (2013).
