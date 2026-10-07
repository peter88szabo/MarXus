# Fragment wells and lumped reactant states

[← README](../../README.md) · related: [Prepared experiments](prepared_experiments.md) · [Distributions](preparation_distributions.md) · [Transient diagnostics](transient_diagnostics.md)

**Deck:** MarXus keywords in a `Bimolecular` block (`FragmentWell`, `LumpedState`, `ExcessFragment`, `PartnerConcentration[molecule/cm^3]`) and a `FragmentEnergy` block in the barrier to it. **Methods:** all; time integration is the main use.

## 1. The question it answers

In a standard master equation a bimolecular product is an exit: once formed, its molecules leave the network. Two kinds of experiments need more.

1. **A fragment that reacts on.** A complex C dissociates into a fragment B and a co-fragment A, and B reacts further: it dissociates again, or adds another species. B is then formed with a non-thermal energy distribution, and how much reacts before it is collisionally relaxed depends on how the energy was shared at the dissociation. Example: OH + glyoxal → HC(O)CO + H₂O, after which hot HC(O)CO decomposes promptly or adds O₂ (Shannon, Blitz, Seakins, J. Phys. Chem. A 128, 1501 (2024)).
2. **A reactant that is consumed in time.** In chemical activation the bimolecular reactant A + X is used up while it forms the well, and molecules that redissociate return to it. Following the reactant as a dependent variable replaces the instantaneous "nascent" start (Miller et al., J. Phys. Chem. A 120, 306 (2016), SI-VI).

## 2. Equations

**Fragment wells** (Green, Robertson, Chem. Phys. Lett. 605–606, 44 (2014), eq. 15). C at energy $`E_i`$ has $`X_i = E_i - (E_{B0} + E_{A0})`$ above the pair asymptote.

Forward rate and its reverse:

```math
G_{j \leftarrow i} = k_C(E_i)\,P(\varepsilon_j \mid X_i), \qquad \sum_j P(\varepsilon_j \mid X_i) = 1, \qquad
G_{i \leftarrow j} = G_{j \leftarrow i}\,\frac{f_C(i)}{f_B(j)} ,
```

with the partner A, in excess, in the weight of the fragment states:

```math
f_B(j) = \rho_B(\varepsilon_j)\,e^{-(E_{B0}+\varepsilon_j)/kT}\,\phi_A, \qquad
\phi_A = \frac{Q_{A,\mathrm{int}}\,e^{-E_{A0}/kT}\,C'(\mu)(kT)^{3/2}}{[A]} ,
```

**Why this is consistent:**
- In equilibrium [C]/[B] = K_c[A].
- Detailed balance holds for any kernel P, so all solvers apply.
- Every well gets a weight offset along its couplings, and a cycle that gives a well two different offsets is refused.

**Kernels** (`FragmentEnergy <kind>`):

| kind | P(ε \| X) | source |
|---|---|---|
| `Prior` | $`\propto \rho_B(\varepsilon)\,[\rho_A \otimes \rho_t](X - \varepsilon)`$, ρ_t ∝ E^{1/2} | Green, Robertson 2014, eq. 18 |
| `ModifiedPrior` | $`\propto \rho_B(\varepsilon)^n\,[\rho_A \otimes \rho_t](X - \varepsilon)`$, $`n = n_0 (T/T_\mathrm{ref})^m`$ | Shannon et al. 2024 and their ref. 5 |
| `TwoPieceGaussian` | $`\exp[-(\varepsilon-\mu)^2/2\sigma_L^2]`$ for ε ≤ μ, $`\exp[-(\varepsilon-\mu)^2/2\sigma_R^2]`$ for ε > μ; μ, σ_L, σ_R linear in X | Shannon et al. 2024, eq. 4 |

**Grid handling of the kernels:**
- They are evaluated on the 1 cm⁻¹ cells, or integrated over the grain edges for the Gaussian.
- They are truncated to [0, X] and normalized over the fragment grid.

**Lumped reactant state** (Robertson, Comprehensive Chemical Kinetics 43 (2019), eq. 5.187; Green, Plane, Robertson, Chem. Phys. Lett. 865, 141943 (2025), eq. 2).
- **The state:** the pair A + X is one thermal state of weight $`Q_{AX}(T)\,e^{-E_0/kT}/[X]`$, implemented as a one-grain well at the asymptote.
- **The coupling:** its channels into the wells are the entrance channels of the deck. Their reverse, the association, follows from detailed balance.
- **Consequence:** the initial loss rate of the state is the high-pressure capture rate k∞(T)[X].

## 3. Deck

**Fragment well:** the fragment is both a `Well` and a `Fragment` of the pair, with the same name and the same molecular data.

```
Bimolecular  P_HCOCO_H2O
  Fragment HCOCO ... End
  Fragment H2O ... End
  GroundEnergy[kcal/mol]               0.0
  FragmentWell                         HCOCO
  PartnerConcentration[molecule/cm^3]  1e14
End
Barrier R3 C2 P_HCOCO_H2O
  RRHO
    ...
    FragmentEnergy TwoPieceGaussian
      Mu[1/cm]          -2685.46  0.50158      ! intercept [unit], gradient in X
      SigmaLeft[1/cm]   -132.008  0.053485
      SigmaRight[1/cm]  -1122.11  0.13401
    End
  End
```

`FragmentEnergy Prior` (no keywords) and `FragmentEnergy ModifiedPrior` (`Order`, `TemperatureExponent`, `ReferenceTemperature[K]`) are the other kinds.

**Lumped reactant state:**

```
Bimolecular  R
  Fragment OH ... End
  Fragment GLYOXAL ... End
  GroundEnergy[kcal/mol]               0.0
  LumpedState
  ExcessFragment                       GLYOXAL
  PartnerConcentration[molecule/cm^3]  1e15
End
```

The state appears as a well named after the pair (here `R`). An initial population is put into it with the `Preparation` block, e.g. `Distribution Thermal` with `Well R`; any distribution in a one-grain well is its single state.

| keyword | where | meaning |
|---|---|---|
| `FragmentWell <name>` | Bimolecular | this fragment is a well of the network; the other fragment is the partner in excess |
| `LumpedState` | Bimolecular | the pair is one thermal state of the master equation |
| `ExcessFragment <name>` | Bimolecular, with `LumpedState` | the fragment in excess |
| `PartnerConcentration[molecule/cm^3]` | Bimolecular | concentration of the partner, or of the excess fragment |
| `FragmentEnergy <kind> ... End` | barrier to a pair with a fragment well | the energy partitioning (required) |

**Checks:**
- **Same fragment data:** the fragment's density of states as a `Fragment` must equal the one as a `Well`, otherwise the deck is refused.
- **Required blocks:** a missing `FragmentEnergy` block, or a missing `ExcessFragment` or concentration, is refused.

## 4. Output

- **Fragment wells** are wells like any other: populations, energies, flux coefficients and diagnostics in time, and species in CSE.
- **Fragment channels** appear in the network summary as `C -> B+A` and in the flux-coefficient tables as `C->B`.
- **A lumped state** is reported as a well. Its population is the reactant concentration in units of the initial amount.

## 5. Validity and limits

- **Pseudo-first-order:** the partner, or the excess fragment, is in excess at a fixed concentration.
- **Tabulated kernels** P(ε | X_k), for example from trajectories, are not available yet.
- **No absorbing barrier for a lumped state:** it has no reaction threshold, so the intermediate steady state is not defined for it. Use the final variant or TimeIntegration.

## 6. Validation

**Unit tests:**
- **Kernels:**
  - the two-piece Gaussian with the coefficients of the OH + glyoxal deck gives the prompt fractions 0.332, 0.705, 0.900, 0.996 above TS2 at 5, 15, 25, 50 kJ/mol above TS1;
  - the prior equals the statistical partitioning of the photoion module grain by grain;
  - the modified prior of order 1 is the prior.
- **Operator:**
  - detailed balance and conservation for both kernels, all collision models and both steady-state variants;
  - a closed C ⇌ B system is stationary, and its ratio doubles with [A];
  - from a hot start, the time integration relaxes to [C]/[B] = K_c[A] (10⁻⁶);
  - an inconsistent cycle is refused.
- **Lumped state:** the initial loss rate of the lumped HCO + O₂ state equals the high-pressure capture rate k∞(T)[O₂] within 7.5·10⁻⁴ at 300 and 500 K.

**In preparation:** OH + glyoxal (Shannon et al. 2024), with 91 measured OH yields against [O₂] at 212–295 K and 5–80 Torr.

## 7. Code and references

**Code:**
- `src/masterequation/fragment_partition.rs`;
- `chemical_activation_network.rs` (`ChannelDestination::Fragment`);
- `chemical_activation_operator.rs` (`well_log_offsets`, fragment couplings);
- `mess_input.rs` and `chemical_activation_from_mess_input.rs` (deck).

**References:**
- N. J. B. Green, S. H. Robertson, Chem. Phys. Lett. 605–606, 44 (2014).
- R. J. Shannon, M. A. Blitz, P. W. Seakins, J. Phys. Chem. A 128, 1501 (2024).
- S. H. Robertson, "Foundations of the master equation", Comprehensive Chemical Kinetics 43 (2019), Ch. 5.
- J. A. Miller, S. J. Klippenstein, S. H. Robertson, M. J. Pilling, R. Shannon, J. Zádor, A. W. Jasper, C. F. Goldsmith, M. P. Burke, J. Phys. Chem. A 120, 306 (2016).
- N. J. B. Green, J. M. C. Plane, S. H. Robertson, Chem. Phys. Lett. 865, 141943 (2025).
- B. Sztáray, A. Bodi, T. Baer, J. Mass Spectrom. 45, 1233 (2010) (the statistical partitioning reused for the prior).
