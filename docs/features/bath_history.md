# Bath history: shocks and other changes of T and p

[← README](../../README.md) · prepared experiments: [Prepared experiments](prepared_experiments.md) · [Distributions](preparation_distributions.md) · [Time profiles](source_time_profiles.md) · **Bath history** · [Transient diagnostics](transient_diagnostics.md) · [CSE source projection](cse_source_projection.md) · [CSE validity in time](cse_validity_in_time.md) · [Fragment wells](fragment_wells_and_lumped_reactants.md)

**Deck:** a `Bath ... End` block in the [`Preparation` block](prepared_experiments.md). **Method:** TimeIntegration.

## 1. The question it answers

What happens when the temperature or pressure of the bath changes during the experiment? Examples:
- the passage of an incident and a reflected shock;
- a temperature jump;
- a pressure change.

## 2. The deck

```
Bath
  Segment Start[s] 0      Temperature[K] 300    Pressure[atm] 1
  Segment Start[s] 1e-6   Temperature[K] 1000   Pressure[atm] 1
End
```

- **Segments:** piecewise-constant conditions. The first segment starts at 0 s; the starts increase.
- **Pressure units:** `[torr]`, `[atm]`, `[bar]`.
- **Without a `Bath` block:** the time integration runs once per (T, p) of the deck's `TemperatureList` and `PressureList`.

## 3. What happens at a segment start

1. **Operator.** The operator 𝐉 is assembled for the new conditions on the same grain grid: collision frequency, ⟨ΔE_down⟩(T), the kernel and its low-energy reservoir.
2. **Populations.** They are carried over through the grains, so the vector is continuous across the jump and the total is conserved. A reservoir state of the old bath is first spread with its Boltzmann weights; the grains are then summed into the partition of the new bath.
3. **Sources.** Every source distribution is projected again onto the states of the new segment, and the report lists that projection.
4. **Integrator.** It restarts with its own factorization cache, so no factor of the earlier operator is reused.

## 4. Output

**Report:**
- the segments (start, T, p);
- the projection notes per segment;
- one set of time tables for the whole history.

**Diagnostics per segment:**
- the incubation time is measured from the start of the current segment ([Transient diagnostics](transient_diagnostics.md));
- E_f refers to the operator of the segment;
- the collision numbers add up over the segments.

## 5. Validity and limits

- **What changes.** Only T and p change; the energy release of the reaction does not change T (an isothermal bath).
- **Step changes only.** A continuously varying history (tabulated T(t), p(t)) is not available yet; it can be approximated by short segments.
- **Other methods.** The steady-state and CSE methods use the deck's conditions, not the history.

## 6. Validation and code

**Tests** (`prepared_time_integration.rs`):
- a constant history split into segments equals one segment;
- a jump 300 → 1500 K conserves the population;
- the collision numbers add up over two segments.

**Example:** `examples/c2h3_prepared_experiment.inp`.
