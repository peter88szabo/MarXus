# SteadyStateAbsorbingBarrier: the intermediate steady state

[← README](../../README.md) · the four methods: [SteadyStateOlzmann](steady_state_olzmann.md) · **SteadyStateAbsorbingBarrier** · [CSE](chemically_significant_eigenvalues.md) · [TimeIntegration](direct_time_integration.md)

**Family:** steady state. **Run with:** `Method SteadyStateAbsorbingBarrier`, or `--method steady-state-absorbing-barrier`.

## 1. The question it answers

What happens to the freshly formed, chemically activated adducts on their first collisional descent? Which fraction redissociates, which decomposes to each product (chemical activation), and which is stabilized in each well?

**The intermediate steady state** (SN84; GO10 Sec. 3.2) is "characterized by the competition between decomposition and stabilization of the chemically activated adduct". GO10 states that it "can be implemented by introducing a lower absorbing barrier into the master equation". This is also the approach of SSUMES.

## 2. Equations

**Operator.** The same master equation as [SteadyStateOlzmann](steady_state_olzmann.md). Each well gets an absorbing barrier at

```math
E_{\mathrm{abs}} = E_{\mathrm{thr,min}} - X\,k_BT \qquad (X = 10 \text{ by default}).
```

The grains below the barrier are removed. A molecule transferred below it by a collision, or by isomerization into the grains below the barrier of the target well, counts as stabilized. $`𝐉_{\mathrm{abs}}`$ is the operator on the grains above the barriers. With $`R = 1`$:

```math
𝐉_{\mathrm{abs}}\,N^{s} = F, \qquad \Phi_r = \sum_E k_r(E)\,N^s(E), \qquad \Phi_{\mathrm{stab},w} = \sum_{E \ge E_{\mathrm{abs}}} \omega \sum_{E' < E_{\mathrm{abs}}} P(E' \leftarrow E)\,N^s(E) + \ldots
```

**Stabilization yield.** It also includes the isomerization flux into the absorbed grains, and the part of the source formed below a barrier.

**Balance.** $`\sum_r \Phi_r + \sum_w \Phi_{\mathrm{stab},w} + \Phi_{\mathrm{sink}} = 1`$.

**Bimolecular rate coefficients of the reactant** (PR03 eq. 44):

```math
k(R \to P) = k_\infty\,\Phi_P \quad \text{(bimolecular-to-bimolecular, chemical activation)}, \qquad k(R \to W) = k_\infty\,\Phi_{\mathrm{stab},W} \quad \text{(bimolecular-to-well, stabilization)},
```

with $`k_\infty`$ the capture (high-pressure association) rate coefficient of the entrance channels, from the same $`W^\ddagger`$ as the source.

## 3. Algorithm

The same banded Cholesky solve with iterative refinement as [SteadyStateOlzmann](steady_state_olzmann.md), on the grains above the absorbing barriers.

## 4. Settings

| deck keyword | option | default |
|---|---|---|
| `Method SteadyStateAbsorbingBarrier` | `--method steady-state-absorbing-barrier` | required (no default method) |
| `AbsorbingBarrierBelowThreshold[kT]` | `--barrier-kt X` | 10 |

## 5. Output

The tables come in three views: by temperature, by pressure, and temperature–pressure.
- **Yields:** % of the formed adducts (products, `stab(W)`, `escape(W)`), and without the return to the reactant (% of the net reaction).
- **Chemical-activation rate coefficients** $`k^{ca}`$ (1/s).
- **Bimolecular-to-bimolecular rate coefficients (chemical activation):** $`k(R \to P) = k_\infty \Phi_P`$ (cm³/s), summed over the channels that form P, and `R->escape(W)`. With the matching **yields** (% of the net reaction).
- **Bimolecular-to-well rate coefficients (stabilization):** $`k(R \to W)`$ (cm³/s). With the matching **yields**.
- **Capture, return and net reaction of R** (cm³/s).

**Machine-readable:** the blocks `# intermediate steady state` and `# bimolecular rate coefficients of R [cm3/s], intermediate steady state`; `FILE_tables.csv`.

## 6. Validity and limits

**Time window.** $`(0.1\,\lambda_F)^{-1} < t < (10\,k_{\mathrm{uni}})^{-1}`$ (O02 p. 3618; SN84), with $`\lambda_F`$ the eigenvalue whose eigenvector has the largest weight in $`F`$. In this window the activated population has relaxed and the stabilized adducts have not yet reacted.

**The barrier position is a choice.** For shallow wells, or at high T where the well depth approaches 10 k_BT, the results depend on it. A barrier below the well bottom is refused with an explanation.

**With a physical bimolecular sink** the barrier is "artificial", and "a too low product yield would be predicted" (O02 p. 3617). [SteadyStateOlzmann](steady_state_olzmann.md) is then the right solver.

**What it does not give:** the later thermal reaction of the stabilized adducts, and thermal rate coefficients.

## 7. Validation

**ZZ-allyl + O₂, four wells** (MESS Eckart model, 760 Torr, 270–330 K):
- **Prompt products.** k(R → IEPOX + OH) = k∞Φ_P5 agrees with CSE's G13 eq. 21 within 0.44% (P1 0.11%, P7 0.18%), and is 0.44–0.74% above MESS.
  - With the exact Eckart tunneling it is +19.0 … +20.3% above MESS: the difference is the tunneling model.
- **Stabilization, against MESS:** k(R → G4) −18.8 … −9.5%, k(R → G2) −4.8 … +0.2%.
  - Total stabilization is 3.0–9.9% below CSE's.
  - The barrier counts the flux 10 kT below the lowest threshold, a different definition from CSE's chemical eigenmode.
- **Decomposition.** The IEPOX + OH prompt yield, plus the stabilization yields × the thermal fates of the wells (from SteadyStateOlzmann), equals the SteadyStateOlzmann total within 1.8·10⁻³ relative (`reports/method_comparison.md`).

**H + C₂H₂ ⇌ C₂H₃** (`validation/c2h3_mess_example/`):
- **300–1000 K.** The association falloff is within a few % of MESS, with the residual from the tunneling model. It equals SteadyStateOlzmann's k_uni·K within 0.02%.
- **Above about 1500 K** the result depends on the barrier distance: −8.3% at 1500 K and −29.3% at 1750 K (0.1 atm, 10 kT), and the barrier falls below the well bottom at 2000 K. With 3 kT the agreement at 1 atm is within 3.6% at 1000–1750 K.

## 8. Code

- `chemical_activation_network.rs` (`AbsorbingBarrier`), `chemical_activation_operator.rs` (absorbed grains, stabilization flux), `chemical_activation_steady_state.rs`, `chemical_activation_observables.rs`, `report_sections.rs`.
- Report: `reports/olzmann_absorbing_barrier_and_treatment.md`; validation: `validation/c2h3_mess_example/README.md`, §5 (barrier-distance sensitivity).

## 9. References

- SN84: K. Schranz, S. Nordholm, Chem. Phys. 85, 163 (1984).
- O02: J. Olzmann, Phys. Chem. Chem. Phys. 4, 3614 (2002).
- GO10: G. González-García, J. Olzmann, Phys. Chem. Chem. Phys. 12, 12290 (2010).
- PR03: M. J. Pilling, S. H. Robertson, Annu. Rev. Phys. Chem. 54, 245 (2003).
- CD07: J. Carstensen, A. M. Dean, Comprehensive Chemical Kinetics 42 (2007), the 10 k_BT convention.
