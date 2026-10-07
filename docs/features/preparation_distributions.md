# Distributions: which molecules are prepared, and at which energies

[← README](../../README.md) · prepared experiments: [Prepared experiments](prepared_experiments.md) · **Distributions** · [Time profiles](source_time_profiles.md) · [Bath history](bath_history.md) · [Transient diagnostics](transient_diagnostics.md) · [CSE source projection](cse_source_projection.md) · [CSE validity in time](cse_validity_in_time.md) · [Fragment wells](fragment_wells_and_lumped_reactants.md)

**Deck:** `Distribution <kind> ... End` inside `InitialPopulation` or `Source <name>` of the [`Preparation` block](prepared_experiments.md).

## 1. The question it answers

How is the internal energy of the prepared molecules distributed, and over which wells? The distribution sets everything that follows:
- whether molecules react before they relax (hot);
- whether they are first activated (cold);
- how the products branch.

## 2. Kinds

| kind | purpose | controls |
|---|---|---|
| `Thermal` | Boltzmann equilibrium at a temperature of its own: the gas before a shock, a cold or hot preparation | `Well`, `PreparationTemperature[K]` |
| `Gaussian` | a band of energies: laser or photolysis preparation with a finite energy spread | `Well`, `Centre[unit]`, `Width[unit]` (standard deviation σ), `Representation` |
| `SingleEnergy` | all molecules at one energy (δ preparation) | `Well`, `Energy[unit]` |
| `Tabulated` | any distribution from a file: trajectories, another master equation, an experiment | `Well`, `File`, `Representation`, `EnergyUnit`, `Support` |
| `ThermalEntrance` | the thermal chemical-activation flux of the deck's `Reactant` at the temperature of the first bath segment | – |
| `Mixture` | several populations at once (two wells, two energy bands) | sub-blocks `Component <weight>`, one distribution each; the weights sum to 1 and set the well fractions |

**Options of every distribution:**
- `EnergyReference AboveWellGround` (default): energies above the ZeroEnergy of the well.
- `EnergyReference Absolute`: the energy scale of the deck.
- `Shift[unit]`: a shift of the whole distribution, e.g. by a photon energy.

Energy units: `[1/cm]`, `[kcal/mol]`, `[kJ/mol]`.

## 3. Equations

**Gaussian:**

```math
f(E) \propto \exp\!\left[-\frac{(E - E_c)^2}{2\sigma^2}\right] .
```

**Representations** (how a value is turned into a grain mass):
- **`Density`** (default for Gaussian): the density is integrated over each grain [E_i − ΔE/2, E_i + ΔE/2], with erf/erfc accurate in both tails, then normalized on the grid. The mass outside the grid is reported.
- **`PerStateWeight`:** the value is a weight per state, so the grain mass is $`w(E_i)\,\rho(E_i)\,\Delta E`$, as for a Boltzmann factor.
- **`BinMass`** (default for `Tabulated`): the value is the mass of the bin.

**Tabulated** files have lines "lower upper value" (whitespace or commas, `#` comments).
- **Rebinning:** the bins are rebinned onto the grains in proportion to their overlap, which conserves the mass; a zero-width bin goes to the grain that contains it.
- **`Support Complete`** (default): mass outside the grid of the well is an error.
- **`Support Truncate`:** that mass is removed and reported.

**Shifts** are rebinned the same way.

**Normalization:** every distribution is normalized over all wells together. The well fractions of a `Mixture` are its component weights.

## 4. Output

For the initial population and every source, the report gives:
- its description;
- its fraction in every well;
- the mean energy above the well bottom;
- the probability that lay outside the grid.

**Low-energy reservoirs.** Where a well has a reservoir state (exponential down, where the normalization fails at the sparse well bottom), the part of a distribution inside the reservoir grains takes the reservoir's Boltzmann shape at the bath temperature. A cold preparation in a hot bath is therefore partly heated at once.

The report lists, per bath segment and distribution:
- the fraction in a reservoir;
- the fraction below an absorbing barrier (intermediate steady state; counted as stabilized at once);
- the mean energy before and after the projection.

## 5. Validity and limits

- **Energy only.** The distributions depend on energy only, as the master equation does (no angular-momentum resolution).
- **Grid.** Grains refer to the well bottom, which is placed on the common grain grid. Mean energies at coarse grains can therefore differ by up to half a grain.

## 6. Code

`src/masterequation/prepared_distributions.rs`:
- `thermal`, `gaussian`, `tabulated`, `mixture`, `shifted`, `rebin`;
- `GrainDistribution`.

`preparation_input.rs`: `build_distribution`, `read_distribution_file`.

Tests:
- conservative rebinning;
- thermal at its own temperature;
- the Gaussian integral and its tail accuracy;
- the per-state weight;
- tabulated reference scales;
- support errors;
- mixtures;
- shifts;
- the text of the descriptions.
