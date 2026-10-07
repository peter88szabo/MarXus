# Time profiles: when the source molecules are added

[← README](../../README.md) · prepared experiments: [Prepared experiments](prepared_experiments.md) · [Distributions](preparation_distributions.md) · **Time profiles** · [Bath history](bath_history.md) · [Transient diagnostics](transient_diagnostics.md) · [CSE source projection](cse_source_projection.md) · [CSE validity in time](cse_validity_in_time.md) · [Fragment wells](fragment_wells_and_lumped_reactants.md)

**Deck:** `Profile <kind> ... End` inside `Source <name>` of the [`Preparation` block](prepared_experiments.md).

## 1. The question it answers

At which rate R(t) does a source add molecules? This covers short or long laser pulses, photolysis with the decay of a precursor, continuous flows, and measured rate histories.

## 2. Kinds

| kind | R(t) | controls |
|---|---|---|
| `Impulse` | $`N_p\,\delta(t - t_p)`$: an exact population jump | `Time[s]`, `Amount` |
| `Rectangular` | $`N_p/(t_2 - t_1)`$ on $`[t_1, t_2)`$ | `Start[s]`, `Stop[s]`, `Amount` |
| `Gaussian` | $`N_p\,e^{-(t - t_c)^2/2\sigma^2}/(\sqrt{2\pi}\,\sigma)`$ for t ≥ 0, renormalized if it starts before 0 | `Centre[s]`, `Width[s]` (σ), `Amount` |
| `Feed` | $`R_0`$ on $`[t_1, t_2)`$, open-ended without `Stop` | `Start[s]`, optional `Stop[s]`, `Rate[1/s]` |
| `PrecursorDecay` | $`k_f N_\mathrm{prec}\,e^{-k_\mathrm{tot}(t - t_1)}`$ for t ≥ t₁: formation by a first-order precursor with competing losses ($`k_\mathrm{tot} \ge k_f`$) | `Start[s]`, `FormationRate[1/s]`, `TotalLossRate[1/s]`, `PrecursorAmount` |
| `Tabulated` | measured (t, R), linear between the points, zero outside | `File` (lines "time rate", s and 1/s) |
| `Train` | the sum of several profiles: repetitive pulses, a pulse with a background | sub-blocks `Profile <kind>` |

## 3. How the integrator treats them

**Impulses** are not rates. They are applied as exact jumps $`n(t_p^+) = n(t_p^-) + N_p F`$.

**Events.** The integrator stops at every impulse, at every rate edge (rectangular and feed limits), at every bath change and at every output time, and restarts its step control there.

**Between events:**
- it runs autonomously where all rates are constant;
- otherwise it runs non-autonomously, with $`\partial f/\partial t = \sum_a \dot R_a F_a`$.

**Injected amounts.** $`N_{a,\mathrm{in}}(t) = \int_0^t R_a\,dt' + \sum_{t_p \le t} N_p`$ is given in closed form for every profile and reported at every output time. It is the denominator of the yields and of the balance.

## 4. Output

- **Run settings:** each source with its profile in words, e.g. "impulse at 1e-9 s, amount 0.5".
- **Report:** the total amount, or "open-ended" for a feed without `Stop`.
- **Time tables:** the injected amount of every source.

## 5. Validity and limits

- **Output times.** An output time that coincides with an impulse is reported after the jump.
- **Narrow pulses.** A rectangular or Gaussian pulse much shorter than the relaxation behaves as an impulse with the same amount (tested).

## 6. Code

- `src/masterequation/source_profiles.rs`: `TimeProfile` with `rate`, `impulses`, `breakpoints`, `injected_until`, `total_amount`; `Display`.
- `prepared_time_integration.rs`: events and jumps.

Tests:
- impulses are jumps, not rates;
- rectangular and Gaussian integrals;
- the injected amount of a feed;
- precursor decay;
- tabulated profiles;
- trains;
- invalid profiles;
- the text of the descriptions;
- a narrowing rectangular pulse approaches the impulse.
