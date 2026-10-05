# Chemical activation in the energy-grained master equation: Olzmann, SSUMES and MESS compared

**MarXus** — background for implementing chemical activation (CA)  
**Date:** 2026-10-05  
**Status:** literature and source-code review. No MarXus code was changed for this report.

---

## 0. How this report was made

- **Papers.** Every paper in `MarXus/papers/ChemAct/` was read, plus Miller & Klippenstein 2006,
  Pilling & Robertson 2003, Caralp et al. 2008 and Johnson & Green 2022 from `papers/Master_Equiton/`.
- **Source code.** SSUMES (`/home/peter/Programs/SSUMES/ssumes/source/`) and MESS 2026
  (`/home/peter/Dropbox/Research_Leuven/MESS_kinetics/Source_from_2026/MESS/src/`).
- **Citations.** Every statement below cites a paper (equation or page) or a source file and line.
  Quotations used for the main conclusions were checked against the PDF text or the source code.
- **Own derivations** are marked **[derived]** and are not taken from a publication.
- **Earlier internal notes** (LaTeX notes on SSUMES, Olzmann and MESS) were treated as hypotheses only.
  Section 6 lists which of their claims the literature supports and which it contradicts.

**Abbreviations**

| Code | Reference |
|---|---|
| O91 | Olzmann, Gebhardt, Scherzer, *Int. J. Chem. Kinet.* **23**, 825 (1991) |
| O02 | Olzmann, *Phys. Chem. Chem. Phys.* **4**, 3614 (2002) |
| GO10 | González-García, Olzmann, *Phys. Chem. Chem. Phys.* **12**, 12290 (2010) |
| PO14 | Pfeifle, Olzmann, *Int. J. Chem. Kinet.* **46**, 231 (2014) |
| MK06 | Miller, Klippenstein, *J. Phys. Chem. A* **110**, 10528 (2006) |
| G13 | Georgievskii, Miller, Burke, Klippenstein, *J. Phys. Chem. A* **117**, 12146 (2013) |
| SN84 | Schranz, Nordholm, *Chem. Phys.* **85**, 163 (1984) |
| TR77 | Tardy, Rabinovitch, *Chem. Rev.* **77**, 369 (1977) |
| S89 | Smith, McEwan, Gilbert, *J. Chem. Phys.* **90**, 4265 (1989) |
| V97 | Vereecken, Huyberechts, Peeters, *J. Chem. Phys.* **106**, 6564 (1997) |
| PR03 | Pilling, Robertson, *Annu. Rev. Phys. Chem.* **54**, 245 (2003) |
| CD07 | Carstensen, Dean, *The Kinetics of Pressure-Dependent Reactions*, Comprehensive Chemical Kinetics **42** (2007) |
| C08 | Caralp, Forst, Bergeat, *Phys. Chem. Chem. Phys.* **10**, 5746 (2008) |
| JG22 | Johnson, Green, *Faraday Discuss.* **238**, 380 (2022) |
| Z22 | Zhang, Chen, Truhlar, Xu, *Faraday Discuss.* **238**, 431 (2022) |
| B01 | Barker, Ortiz, *Int. J. Chem. Kinet.* **33**, 246 (2001) |

---

## 1. Summary

1. **One equation, three ways of reading it.** All three approaches assemble the same kind of linear
   energy-grained master equation: collisional energy transfer plus microcanonical reaction sinks plus a
   source of chemically activated adducts. They differ in five places:
   - which **observable** is extracted (flux coefficients or yields vs. phenomenological rate coefficients);
   - which **steady state** is meant (intermediate or final);
   - how **stabilization** is defined;
   - how **multi-step (consecutive) activation** is coupled;
   - the **collision model** and its low-energy treatment.

2. **The chemical-activation source is the same for thermal reactants.**
   - Olzmann: F(E) ∝ W‡(E−E₀)·e^(−(E−E₀)/kT) (GO10 eq 15, PO14 eq 7).
   - SSUMES: r(E) ∝ ρ(E)·e^(−E/kT)·k_back(E) (`ssulibc.cc:1174`).
   - MESS: N‡(E)·e^(−E/kT)/2π (`mess.cc:11127`; G13 eq 9).
   - Because k(E) = W‡(E−E₀)/(hρ(E)) (GO10 eq 13), these are identical after normalization.
   - Olzmann additionally provides **non-thermal sources**: convolution and the shift approximation
     (PO14 eqs 8, 10–13). SSUMES and MESS do not.

3. **The decisive physical difference is stabilization: intermediate vs. final steady state.**
   - SN84 (abstract, p. 174) describes three time scales after the source is switched on: a transient,
     an *intermediate steady state* (activated adducts decompose or are stabilized), and an *asymptotic
     (final) steady state*. In the final state there is no more net stabilization, because the stabilized
     adducts react thermally.
   - **SSUMES** computes the intermediate steady state by an absorbing barrier: grains below a
     user-chosen `truncate` index leak out as k_stab (`ssulibc.cc:1069–1073`).
   - **Olzmann** used this only in O02 and called it "artificial" when a bimolecular sink exists:
     "If the absorbing barrier would be retained, an irreversible loss … would result and as a
     consequence a too low product yield would be predicted" (O02 p. 3617). GO10 (p. 12295) and PO14
     (p. 238) use the **final** steady state without a barrier.
   - **MESS** represents stabilization by the slow chemical eigenvector, i.e. the R → well rate
     coefficient (G13 eq 28, p. 12153).

4. **The eigenvalue (MESS) formulation does not lose the formation physics.**
   - The nascent distribution enters explicitly (G13 eqs 9, 21). Prompt, well-skipping products come from
     the relaxational eigenmodes (G13 eq 21; MK06 worked well-skipping examples).
   - Its condition of validity is **separation of the chemically significant eigenvalues from the
     energy-relaxation eigenvalues, plus a source that varies slowly on the relaxation time scale**
     (G13 p. 12147, eq 19; JG22 eqs 29, 40). It is *not* simply "thermalization faster than reaction".
   - When the separation fails, phenomenological rate coefficients do not exist. MESS then merges species
     or drops wells into a bimolecular group (`mess.cc:12840–12846`). The steady-state and time-dependent
     master-equation solutions remain well defined in that situation.

5. **Rate coefficients and flux coefficients are different quantities.** MK06 (pp. 10529, 10531):
   "application of the steady-state approximation … is virtually always an attempt to equate a
   phenomenological rate coefficient to a flux coefficient. Sometimes this is a valid approach, and
   sometimes it is not." The steady-state CA quantities of Olzmann and SSUMES (k^ca, yields) are
   flux-type quantities. They are the natural observables of a continuously fed experiment.
   MESS's quantities are phenomenological rate coefficients for kinetic models.

6. **Is the steady-state (SSUMES/Olzmann) picture "more real" than the eigenvalue picture?** My
   assessment, based on the sources above:
   - **No, not in general.** The two answer different questions.
   - **Yes, for specific cases:** continuous formation, removal of the energized or stabilized adduct by
     a bimolecular reaction (O2 capture, PO14), consecutive activation with a non-thermal parent
     distribution, and conditions where eigenvalues do not separate. There the Olzmann **final**
     steady-state formulation is the appropriate one, and it is more complete than SSUMES's absorbing
     barrier.
   - Where the separation holds, the three approaches agree:
     - steady-state yields equal time-integrated single-injection yields (V97, Table I: "identical for
       all practical purposes");
     - MESS's R → P and R → well rate coefficients equal the fractional populations accumulated on the
       relaxation time scale (G13 p. 12153).

---

## 2. The common equation and its sign conventions

| | Olzmann | SSUMES | MESS (G13) |
|---|---|---|---|
| Master equation | dN/dt = R·F − J·N, with J = ω(I − P) + K + k_c[D]·I (PO14 eq 2; O02 eq 6; GO10 eq 7) | dg/dt = J·g + k_in·r (manual eq 8) | d\|f⟩/dt = −Ĝ\|f⟩ + Σ_ν s_ν\|p^(ν)⟩ (G13 eq 1) |
| Sign of the operator | J positive definite | J = −(Olzmann J) | Ĝ positive definite |
| Steady state | N^s = R·J⁻¹·F (PO14 eq 5) | g = −k_in·J⁻¹·r (manual eq 9) | f_λ = Σ_ν s_ν p_λ/Λ_λ for relaxational modes (G13 eq 19) |
| Matrix storage | — | `b[source][target]` (`ssulibc.cc:1051, 1332`) | rows/columns symmetrized |
| Symmetrization | after detailed balance (O02 p. 3616; tql2 after symmetrization) | b[source][target]·s[source]/s[target], i.e. **W⁻¹MW** in target-row form (`ssulibc.cc:1285`) | kernel·√f_i/√f_j (`mess.cc:11085–11087`) |

**Note on the former MarXus bug.** SSUMES's documented transform "S·B·S⁻¹" is correct only for its
`b[source][target]` storage. MarXus stores R[target, source], for which the same formula is inverted.
That is the likely origin of the inverted transform fixed today (see
`collision_kernel_detailed_balance_and_normalization.md`).

---

## 3. Comparison table

| Aspect | **Olzmann** (O91, O02, GO10, PO14) | **SSUMES** (source code + manual) | **MESS** (G13, MK06, source code) |
|---|---|---|---|
| **Purpose / observable** | CA rate coefficients k_r^ca = Σ_i (K_r Ñ^s)_i (GO10 eq 9) and yields Φ_r = Σ_i (K_r J⁻¹F)_i (O91 eq 5; O02 eq 10) | population-averaged k_ch = Σ g_i k_ch(E_i)/G and fractions k/k_out,tot (`ssulibc.cc:1412–1443`); yields per molecule formed | phenomenological rate coefficients: well↔well, well↔bimolecular, bimolecular→bimolecular (G13 eqs 21, 27, 28, 30; `mess.cc:12262–12965`) |
| **State space** | one energized intermediate per ME; several MEs chained (PO14) | several wells in one matrix; bimolecular species are sinks | wells only; bimolecular species are sources and sinks, not states (G13 abstract; `mess.cc:10968`) |
| **Source (thermal reactants)** | F ∝ W‡(E−E₀)·e^(−(E−E₀)/kT) (GO10 eq 15, PO14 eq 7; O91 eq 12) | r_i ∝ ρ_i γ^i k_back(E_i), normalized on [lowest, upbound) (`ssulibc.cc:1174–1178`) | N‡(E)·e^(−E/T)/2π (G13 eq 9; `mess.cc:11127`) |
| **Source (non-thermal)** | convolution of reactant distributions (PO14 eq 8; O91 eq 8); shift approximation (PO14 eqs 12–13; O91 eq 19) | only an externally computed vector (`excitArbit`) or a single grain (`excitGrain`) | not provided (`TimeEvolution` allows a thermal start at another temperature, `mess.cc:11726`) |
| **Collision model** | stepladder with detailed balance; tridiagonal J (O91 eqs 13–18; O02 eqs 13–15) | exponential down, top-down normalization, banded (`ssulibc.cc:1015–1064`) | exponential down (multi-exponential possible), pairwise normalization (`mess.cc:6566–6614`) |
| **Normalization and detailed balance** | completeness and detailed balance imposed (O91 eqs 13–14; GO10: "carefully observing the completeness of transition probabilities") | detailed balance exact; normalization exact above `ncutlow` (back substitution, `norm[jc] = dnorm/(1−usum)`) | detailed balance exact (pairwise common factor); per-source normalization approximate |
| **Low-energy treatment** | — (stepladder) | reduction factors redfac_i = (ρ_{i+upref}/ρ_i)^m ≥ 1 below `ncutlow`, m increased automatically from 1.0 to 3.05 (`ssulibc.cc:1025–1032, 991–1005`) | none |
| **Stabilization** | final steady state, no absorbing barrier (GO10 p. 12295, PO14 p. 238). O02 used a barrier only as the "intermediate steady state", which it calls artificial with a bimolecular sink. O91: S = R′ − D₂ (eq 11) | absorbing barrier: leakage into grains below `lowest` gives k_stab (`ssulibc.cc:1069–1073, 1414–1422`). Default `truncate 0` means no stabilization | chemical eigenvector, i.e. the R → well rate coefficient (G13 eq 28, p. 12153) |
| **Bimolecular sink of the adduct** | k_c[D]·I on the diagonal (O02 eq 5; GO10 eq 6; PO14 eqs 1–2) | not available | `Escape` pseudo-first-order rate added to the diagonal (`model.cc:28679`; `mess.cc:11060`) |
| **Consecutive activation** | yes: solve ME₁, normalize ñ₁^ss, build f₂ by convolution or shift, solve ME₂ (PO14 eqs 10–13) | no (exactly one source, `ssulibc.cc:309–312`) | no |
| **Solution method** | banded LU / tridiagonal (O91 p. 831; GO10 bandec/banbks); eigenvalues by tql2 after symmetrization; 128-bit word length in O02 | Cholesky DPOSV on the symmetrized −J with symmetry check (`ssulibc.cc:740, 1287–1302`); DLSODES for time dependence; DSYEVR for the eigenvalue mode | full diagonalization (default `direct`), or `low-eigenvalue` / `well-reduction` (`mess_driver.cc:931–936`) |
| **Precision** | 128-bit word length (O02) | double | double by default; extended precision only when compiled with WITH_MPACK, which CMake does not set (`mess_driver.cc:1046–1050`) |
| **Thermal rate coefficient** | k^th = lowest eigenvalue of J (GO10 eq 12) | `diseig` (largest negative eigenvalue) or iterative `dislit` | chemically significant eigenvalues |
| **Time dependence** | N(t) = R Σ e_i (1 − e^(−λ_i t))/λ_i E_i (PO14 eq 3) | `catime`, DLSODES (`ssulibc.cc:803`) | `TimeEvolution`: analytic from the eigen-decomposition, constant source (`mess.cc:11597–11808`) |
| **Validity statements** | intermediate steady state exists for (0.1·λ_F)⁻¹ < t < (10·k_uni)⁻¹; with a bimolecular reaction of the adduct it "remains an adequate assumption as long as 10·k_uni < k[B] < 0.1·γ_c·ω" (O02 p. 3618); final steady state is reached when the experimental time scale is distinctly larger than 1/k_uni (O02 p. 3618; GO10 p. 12295) | "least-negative eigenvalue … is not always the solution needed" (manual) | eigenvalue separation; merging otherwise (G13 eqs 31–34; `mess.cc:12092–12141`) |

---

## 4. Discussion of the differences

### 4.1 Source term

For thermal reactants the three sources are the same distribution **[derived]**:

```
ρ(E)·k_back(E)·e^(−E/kT) = (W‡(E−E₀)/h)·e^(−E₀/kT)·e^(−(E−E₀)/kT)  ∝  W‡(E−E₀)·e^(−(E−E₀)/kT)
```

This follows from k(E) = W‡(E−E₀)/(hρ(E)) (GO10 eq 13). The constant cancels on normalization. B01
eq 9 uses the same ρ·k form; V97 (p. 6567) uses the G‡ form.

Minor caveats:

- O02 approximates W‡ by the sum of states of the adduct.
- SSUMES places no source below its truncation grain or above `upbound`.
- MESS uses one vector for both entrance and exit (`mess.cc:11127`).

The **non-thermal** sources are what distinguish Olzmann:

- **Convolution** of two non-thermal reactant distributions (PO14 eq 8). This assumes a reactive cross
  section independent of the internal energies, equivalent to vibrational PST (PO14).
- **Partner convolution** for consecutive activation (PO14 eqs 10–11), with the thermal partner (O2)
  distribution.
- **Shift approximation** (PO14 eqs 12–13): f₂(E) = ñ₁^ss(E + RE − ⟨E⟩_O2), with RE < 0. The
  distribution is "shifted to higher energies by the sum of the reaction energy at 0 K and the average
  thermal energy of O2". O91 eq 19 is the variant without the partner's energy (a δ-function at zero).
- The energy zero is the rovibrational ground state of the intermediate (PO14 p. 234; O02 p. 3615).

MarXus already implements these three Olzmann sources in `chemical_activation_source.rs`, single well
only. The shift-approximation sign follows PO14 eqs 12–13.

### 4.2 Collision model, normalization and the low-energy region

| | Model | Detailed balance | Normalization | Low energies |
|---|---|---|---|---|
| Olzmann | stepladder, P_{i+1,i} = A/(1+A), A = (ρ_{i+1}/ρ_i)e^(−ΔE_SL/kT) (O91 eqs 15–18); ⟨ΔE⟩ = ΔE_SL tanh(ΔE_SL/2F_E kT) (O91 eq 22) | exact | exact | not needed (tridiagonal) |
| SSUMES | exponential down | exact | exact above `ncutlow` | reduction factors on the lower grain of each pair |
| MESS | exponential down | exact (pairwise factor) | approximate | none |
| MarXus now | exponential down, Robertson CCK 43 eq 4.16 | exact | exact above E₀/2 | shared normalization below E₀/2 (UNIMOL rule, Gilbert's `mas55c3.f`) |

The low-energy problem is physical. For sparse states the exponential-down normalization cannot be
satisfied (Robertson CCK 43, p. 294). SSUMES and UNIMOL repair it in different ways. Both keep detailed
balance, both affect only grains far below the reaction threshold, and both are justified by the levels
there staying at equilibrium. Neither repair comes from a model in the literature. MarXus must choose one
(Section 7).

### 4.3 Stabilization: intermediate vs. final steady state

SN84 identifies three time scales after a constant source is switched on: an initial transient, an
**intermediate** steady state, and an **asymptotic** steady state. "The intermediate steady state is well
defined only if the eigenvalues of the corresponding thermal rate matrix can be separated into two groups
one of which contains eigenvalues of substantially smaller magnitude" (SN84 abstract).

- **Intermediate steady state.** Activated adducts either decompose or are stabilized, and the stabilized
  ones are lost from the problem. Technically this means an absorbing barrier (GO10 p. 12295).
  - SSUMES `truncate` does exactly this. The barrier position is user-chosen.
  - The literature places it about 10 kT below the threshold (PR03; CD07 p. 125), or low enough that
    results are insensitive (V97 p. 6568).
  - With weak collisions and a barrier at E₀, errors of 10–30 % occur (SN84 Tables 9–10). Neglecting
    activating collisions from below the threshold "can lead to significant error" (S89 p. 4272).
- **Final steady state.** "There is no more net stabilization; the stabilization reservoir is filled
  up, and time-independent energy distributions have been established. This allows the use of rate
  coefficients instead of eigenvalue-eigenvector expansions" (GO10 p. 12295).
  - Without a bimolecular sink the stabilization flux is zero and the distribution is non-thermal
    (C08 p. 5749). C08 notes that earlier steady-state treatments "implicitly assumed … a reactant sink
    slightly below threshold" (C08 p. 5752).
  - With a sufficiently fast bimolecular removal of the stabilized adducts, "the final steady state
    becomes identical to the intermediate steady state" (PO14 p. 238).
- **Eigenvalue picture.** The R → well rate coefficient is "related to the part of the activated
  reactive complex population that is associated with the chemical eigenstate" (G13 p. 12153). This is
  the eigenvalue analogue of the intermediate-steady-state stabilization, valid under separation.

**Consequence for MarXus.** The choice between intermediate and final steady state is not numerical. It
depends on the experiment, i.e. on the time scale of observation compared with the thermal lifetime of
the adduct (O02 abstract: "true only for reaction times shorter than the thermal lifetime").

### 4.4 Observables: flux coefficients, yields and rate coefficients

- **Olzmann.** k_r^ca = Σ_i (K_r Ñ^s)_i with Ñ^s = J⁻¹F/Σ(J⁻¹F)_i (GO10 eqs 8–9), a population
  average. Yields: Φ_r = Σ_i (K_r J⁻¹F)_i for normalized F (O91 eq 5; O02 eq 10). Branching including a
  bimolecular sink: φ₄ = k₄[H₂O]/(Σk^ca + k₄[H₂O]) (GO10 eq 10).
- **SSUMES.** The same population average k_ch = Σ g_i k_ch(E_i)/G over all wells (`ssulibc.cc:1419`),
  with fractions f = k/k_out,tot. Since the source sums to 1, the fractions are yields per molecule
  formed. The SSUMES manual: these k "are NOT the bimolecular rate coefficients" (Quick Start 1).
  Bimolecular values need k_add,∞ × f, computed outside SSUMES.
- **MESS.** Phenomenological rate coefficients: bimolecular→bimolecular
  k_{ν→μ} = (1/Q_ν) Σ_relax p_λ^(μ) p_λ^(ν)/Λ_λ (G13 eq 21); isomerization, entrance and exit coefficients
  (G13 eqs 27, 28, 30).
- **Relation.** "The rate coefficients for the R → well and R → P reactions are then equal to fractional
  macroscopic populations accumulated in the well and in the products" when integrated over the
  relaxation time, "as long as the chemical eigenvalue is much smaller in magnitude than the ones
  describing energy relaxation" (G13 p. 12153). The flux vs. rate-coefficient distinction (MK06
  pp. 10529, 10531) applies whenever this separation is incomplete.

**Algebraic link [derived by one of the readers from G13 eqs 21, 27, 28, 30; not stated by the authors
and not yet verified numerically]:** the energy-resolved steady-state solve with reversible wells equals
the eigenvalue rate coefficients combined with a macroscopic steady-state approximation on the wells.
MESS itself once computed the bimolecular→bimolecular coefficients by a steady-state Cholesky solve; that
line is commented out and replaced by the eigen-sum (`mess.cc:12193, 12257`).

### 4.5 When do steady-state and eigenvalue results agree?

| Condition | Steady state (Olzmann, SSUMES) | Eigenvalue (MESS) |
|---|---|---|
| Eigenvalues separated, slowly varying source | well defined; equals time-integrated single-injection yields (V97 Table I) | well defined; agrees with time-dependent solution (G13 p. 12153; Z22 Fig. S3; JG22 at long times) |
| Eigenvalues not separated (high T) | final steady state still well defined; intermediate steady state ill defined (SN84 abstract, Table 6) | phenomenological rate coefficients do not exist; species merged (G13 eqs 31–34; MESMER manual p. 64) |
| Source or sink varying on the relaxation time scale | time-dependent ME needed (O02 p. 3615) | residual error (JG22 eq 40, Figs 11–14) |
| Weak collisions, absorbing barrier at E₀ | 10–30 % error (SN84 Tables 9–10) | — |

CD07 (p. 125) states the common requirement: "Both methods essentially rely on reaching steady-state
concentrations … the smallest negative eigenvalues must be clearly separated from those that describe
energy relaxation."

### 4.6 Consecutive (multi-step) activation

Only Olzmann treats it explicitly: "the output of the first master equation governs the input of a
second one" (PO14 p. 237). The coupling is one-way and uses the normalized ñ₁^ss.

In MESS, bimolecular reactants enter as thermal species. G13 removes them from the state space, so a
non-thermal reactant distribution has no place in the formulation. MESMER has a related device: a
non-thermal deficient reactant handled by "pseudo-isomerization" with a Prior fragmentation distribution
(MESMER manual §11.2.6, eqs 11.17–11.19).

### 4.7 Numerics

- **Olzmann:** tridiagonal or banded LU; 128-bit arithmetic in O02.
- **SSUMES:** Cholesky on the symmetrized matrix. It refuses to proceed for asymmetry above 5 % and warns
  above 1 % (`ssulibc.cc:1287–1302`).
- **MESS:** double precision by default. Extended precision requires a compile flag that is not set in
  the CMake build (`mess_driver.cc:1046–1050`).
- **Eigen-decomposition at low T:** MK06 (p. 10535) reports that "it can be difficult numerically to
  obtain accurate eigenvalues and eigenvectors" and recommends quadruple precision or direct time
  integration. JG22 reports matrix condition numbers of 10¹³–10²¹ (Fig. 10). Note that this is a property
  of the matrix itself, so it also affects a linear solve. No paper reviewed here shows that the
  steady-state linear solve is better conditioned than the eigenvalue problem.

---

## 5. Additional findings from the source code

**MESS 2026**

- `new_mess.cc` is not compiled into `mess`. The compiled master equation is in `mess.cc`
  (`CMakeLists.txt:314, 351`).
- The global keyword `Reactant` is only an energy reference: "bimolecular species to use as an energy
  reference" (`model.hh:2425`).
- `ExcessReactantConcentration` exists only inside the `TimeEvolution` block (`model.cc:1227`).
- `TimeEvolution` evaluates populations analytically from the eigen-decomposition with a constant,
  non-depleting source (`mess.cc:11657–11685`).
- The default number of chemical eigenvalues is decided by a gap ratio (`chemical_threshold = −2`,
  `mess.cc:114, 12132`).
- The low-eigenvalue fallback described in the manual for `ChemicalEigenvalueMin` is commented out in
  the direct method (`mess.cc:11321–11565`).

**SSUMES**

- **Untracked loss at the top of the grid.** Upward probabilities into grains at or above `upbound`
  count in the normalization but are not stored (`ssulibc.cc:1047–1051`). That loss appears neither in
  k_stab nor in k_out.
- **Possible out-of-bounds write** in `checkEqkE` (`ssulibc.cc:450–454`).
- **Memory leak.** `delete[] norm, redfac;` frees only `norm` (`ssulibc.cc:1058, 1079`).
- **Manual out of date.** The intro manual names DSYSV; the code uses DPOSV (`ssulibc.cc:740`).
- The SSUMES Fortran directory `unimol/` (UNIMOL's `mas55c3.f`) is not used by the CA solver.

---

## 6. Assessment of the earlier internal notes

| Claim in the notes | Verdict | Evidence |
|---|---|---|
| MESS row-normalizes its kernel, P_ij = w_ij/Σ_m w_mj | **Contradicted** | pairwise common factor, `mess.cc:6566–6614` |
| SSUMES source differs physically from Olzmann's F(E) | **Contradicted** (thermal case) | identical after normalization, Section 4.1 |
| The steady state is a residence-time/flux distribution, not a population | **Contradicted as worded** | all papers call N^s a population (O91; PO14 eqs 3→5; S89 eq 19). It has units of time per formation flux (PO14 Fig. 2), so the residence-time reading is an interpretation, not a different object |
| Consecutive activation by chained MEs with convolution or shift | **Supported** | PO14 eqs 8, 10–13; O91 eqs 8, 19 |
| SSUMES stabilization by leakage below `lowest`; population-averaged k | **Supported** | `ssulibc.cc:1069–1073, 1412–1422` |
| Eigenvalue method contains no formation physics and cannot describe prompt reaction | **Contradicted** | G13 eqs 9, 21, 38, 41; MK06 well-skipping examples |
| Steady state and eigenvalue agree only for rapid thermalization | **Contradicted** (the condition is wrong) | the condition is eigenvalue separation plus a slowly varying source (G13; SN84; CD07 p. 125) |
| The Klippenstein reformulation enlarges the state space and the steady state is the null space | **Contradicted for G13** | G13 removes the reactants from the state space; the enlarged pseudo-first-order form is MK06 eqs 12–13 |
| k_{ν→μ} = Σ_relax p_λ p_λ/Λ_λ | **Supported, with prefactor 1/Q_ν** | G13 eq 21 |
| Eigenvalue extraction is ill-conditioned; steady-state solve is better conditioned | **First half partly supported; second half not supported** | MK06 p. 10535; JG22 Fig. 10 (condition number of the matrix itself) |
| Steady-state yields equal time-integrated single-injection yields | **Supported, conditional** | V97 Table I (requires absorbing stabilized adducts) |
| Stabilization defined by an absorbing barrier; results depend on its position | **Partly supported** | used by TR77, V97, PR03, CD07, Z22; insensitive if low enough (V97 p. 6568); rejected as arbitrary by SN84 (eq 5) |

---

## 7. Proposal: Olzmann's method in MarXus, in combination with SSUMES

MarXus already has the pieces: the corrected master-equation operator (exact detailed balance and
normalization, correct similarity transform, collision frequency from physical constants), the Olzmann
sources (single well), and an SSUMES-style multiwell steady-state solver. The proposal follows Olzmann
where the three differ and keeps the SSUMES features Olzmann does not have.

**Implementation steps**

1. **Sources (Olzmann).** Thermal TS sum of states (PO14 eq 7), non-thermal convolution (PO14 eq 8),
   partner convolution and shift approximation for consecutive activation (PO14 eqs 10–13), extended
   from single-well to multiwell.
2. **Bimolecular sink of the adduct (Olzmann).** Add k_c[D] on the diagonal (O02 eq 5; PO14 eq 2) as an
   explicit channel that is reported in the outputs.
3. **Stabilization: both steady states, explicitly selectable.**
   - **Final steady state** (no absorbing barrier; GO10, PO14). Proposed default.
   - **Intermediate steady state** with an absorbing barrier (SSUMES `truncate`; O02). The barrier
     position is a user input, and MarXus reports the sensitivity of the yields to it.
4. **Observables.**
   - k_r^ca (GO10 eq 9), yields Φ_r (O91 eq 5; O02 eq 10), branching including the sink (GO10 eq 10).
   - SSUMES-style k_stab and fractions in the intermediate-steady-state mode.
   - The thermal k^th from the lowest eigenvalue (GO10 eq 12).
   - Diagnostics: the ratio of the smallest relaxation eigenvalue to the chemical eigenvalue.
5. **Consecutive activation.** Chain MEs: ME₁ → normalized ñ₁^ss → f₂ → ME₂ (PO14).
6. **Time dependence (verification tool).** N(t) from PO14 eq 3, to check the approach to the
   steady state.

**Validation**

1. **Single-injection identity.** Steady-state yields = time-integrated single-injection yields (V97).
2. **Pressure limits.** As P → 0, ñ^s → f(E); as P → ∞, both distributions become thermal and give the
   same k (GO10 p. 12296).
3. **Reproduce published numbers.** The plateau Φ₂ ≈ 0.43 (O02); the PO14 results for the isoprene/OH/O2
   system.
4. **Cross-check with MESS** only where eigenvalues are separated.

**Decisions needed from Peter**

1. **Collision model.** Olzmann uses the stepladder; SSUMES and MESS use exponential down. Offer both,
   and choose which is the default.
2. **Low-energy normalization for exponential down.**
   - (a) UNIMOL rule: shared normalization below E₀/2. This is what MarXus does now.
   - (b) SSUMES reduction factors.

   Both keep detailed balance exactly; neither comes from a published model.
3. **Default steady state.** Final or intermediate.
4. **Absorbing barrier.** Its default position in intermediate mode, e.g. about 10 kT below the lowest
   threshold (PR03; CD07 p. 125).

---

# Appendices: detailed findings

The quotations used in Sections 1–7 were re-checked against the PDF text or the source code. The appendix
entries below come from the reading notes. Each keeps its page, equation or file:line citation for checking;
not every entry was individually re-checked.

## Appendix A — Per-paper extractions (papers/ChemAct and papers/Master_Equiton)

### A.1 Olzmann, Gebhardt, Scherzer, Int. J. Chem. Kinet. 23, 825 (1991) — O91

Two coupled steady-state MEs: 2-butyl, then butane.

**Master equation**
- Eq 1 (p. 827): dn_i/dt = R f_i − ω n_i + ω_nr Σ_j P_ij n_j − k_i n_i.
- ω counts all collisions and ω_nr only the non-reactive ones; ω_nr = ω − ¼ ω_H/C4H9 (eq 21).
- Eq 2: R F = (ωI − ω_nr P + K) N^s = J N^s; eq 3: N^s = R J⁻¹F.

**Collisions**
- Stepladder with upward steps: completeness (eq 13) and detailed balance (eq 14), approximated by eq 15.
- P_{i+1,i} = A/(1+A) with A = (ρ_{i+1}/ρ_i) e^(−ΔE_SL/kT) (eqs 16–18).
- ⟨ΔE⟩ = −ΔE_SL tanh(ΔE_SL/2F_E kT) (eq 22); F_E defined in eq 23.

**Sources**
- Thermal (eq 12, p. 830): f(E) = G‡(E−E₀) e^(−E/kT)/∫G‡ e^(−E/kT) dE. G‡ is the TS sum of states of the
  reverse reaction 2-C4H9 → C4H8 + H.
- Convolution (eq 8): f′(E) = ∫₀^{E−D₀} g_H^eq(ε) n̂^s(E−D₀−ε) dε, with g_H^eq "the positive part of a
  one-dimensional Maxwellian for hydrogen atoms".
- Shift (eq 19): f′(E) = n̂^s_C4H9(E−D₀), D₀ = 34625 cm⁻¹. A δ-function at zero H-atom energy, so no
  mean-energy term.

**Observables**
- Decomposition D1/R = Σ_i (K J⁻¹F)_i for normalized F (eqs 4–5); R = D1 + R′ (eq 6).
- Normalized distribution N̂^s = J⁻¹F/Σ_i (J⁻¹F)_i (eq 7).
- Second step: D2/R′ = Σ(K′J′⁻¹F′)_i (eq 10); R′ = D2 + S (eq 11).
- Stabilization S is therefore formation minus decomposition, a flux balance, not computed as a
  collisional flux. Where the non-decomposing butane flux leaves the grid is not stated; the Fig. 1 caption
  mentions only "lower break-off of the n̂′ curve is due to limited computer memory".

**Other**
- Stationary solution only.
- Validity: "under steady-state conditions, it seems not possible to suppress the second way completely
  by more frequent and/or more efficient collisions" (pp. 833–834).
- Numerics: "Gaussian elimination algorithm for tridiagonal matrices" (p. 831).

### A.2 Olzmann, Phys. Chem. Chem. Phys. 4, 3614 (2002) — O02

`olzmann2002.pdf` and `Olzmann_b203244a.pdf` are the same paper; the checksums differ only through the
download watermark.

**Master equation**
- Eq 5: dn_i/dt = R1 f_i − ω n_i + ω Σ P_ij n_j − k2i n_i − k4[H] n_i.
- Eq 6: R1 F = {ω(I−P) + K2 + k4[H] I} N^s ≡ J N^s.
- Energy "counted from the vibrational ground state of s-C4H9" (p. 3615).

**Collisions**
- Stepladder (eqs 13–14); ⟨ΔE⟩ ≈ ΔE_SL tanh(ΔE_SL/2F_E kT) (eq 15).
- ΔE_SL = 190 cm⁻¹ (⟨ΔE⟩ ≈ 60 cm⁻¹), F_E ≈ 1.33 at 300 K, on a 10 cm⁻¹ grain.
- The ME is split into 19 energetically shifted sub-equations (p. 3616).

**Source and observables**
- Source (eq 11): f(E) = N W(E−E₀(−1)) exp[−(E−E₀(−1))/kT]. W is approximated by the sum of states of
  s-C4H9.
- Normalization (eq 8); R2 = Σ(K2 N^s)_i (eq 9); Φ2 = R2/R1 = Σ(K2 J⁻¹F)_i (eq 10).

**Stabilization**
- Treated "either neglecting upward collisions or … introducing a lower absorbing barrier … below the
  lowest reaction threshold" (p. 3616).
- With a bimolecular sink the "artificial condition of a lower absorbing barrier in the master equation
  must be dropped. If the absorbing barrier would be retained, … a too low product yield would be
  predicted" (p. 3617).

**Steady state vs. eigenvalues**
- "a steady-state master equation with upward collisions included corresponds to this case [final steady
  state, Φ2 = 1] if no absorbing barrier exists" (p. 3616).
- k2^uni is taken as the lowest eigenvalue of J. The plateau Φ2 ≈ 0.43 agrees with an absorbing-barrier
  calculation.

**Validity**
- "true only for reaction times shorter than the thermal lifetime" (abstract).
- Intermediate steady state for (0.1λ_F)⁻¹ < t < (10k_uni)⁻¹. With a bimolecular reaction it remains
  adequate "as long as 10k_uni < k[B] < 0.1γ_cω" (p. 3618).
- The time-dependent ME is needed "at elevated temperatures where the rates of chemically and thermally
  activated reactions … become comparable" (p. 3615).

**Numerics:** bandec/banbks; tql1 after symmetrization; "128 bit word length".

### A.3 González-García, Olzmann, Phys. Chem. Chem. Phys. 12, 12290 (2010) — GO10 (HSO5)

- **Master equation:** sinks k_{−2a}, k_{2b} and k4[H2O] (eq 6); R2a F = [ω(I−P) + K_{−2a} + K_{2b} +
  k4[H2O] I] N^s (eq 7).
- **Collisions:** stepladder "obeying detailed balancing"; ΔE_SL "represents the average amount of energy
  transferred in down collisions" (eq 16).
- **Source:** f(E) = W_{−2a}(E−E₀(−2a)) exp[−(E−E₀(−2a))/kT]/∫₀^∞ W_{−2a}(ε) e^{−ε/kT} dε (eq 15).
- **Microcanonical rates:** k_r(E) = W_r(E−E₀(r))/hρ(E) (eq 13).
- **Observables:**
  - Ñ^s = J⁻¹F/Σ(J⁻¹F)_i (eq 8); **k_r^ca = Σ_i (K_r Ñ^s)_i** (eq 9).
  - φ4 = k4[H2O]/(k^ca_{−2a} + k^ca_{2b} + k4[H2O]) (eq 10); k^ca = k^ca_{−2a} + k^ca_{2b} (eq 11).
- **Thermal rate:** k^th = λ1, the lowest eigenvalue (eq 12). "Alternatively … averaging procedure
  analogous to eqn (9) but with Ñ^s = Ñ^s_th [lowest eigenvector]".
- **Steady state:** the intermediate steady state uses a lower absorbing barrier; "by omitting this
  absorbing barrier … the final steady-state solution … is obtained" (p. 12295). The paper uses the final
  steady state.
- **Pressure limits:** as P → 0, ñ^s → f(E); as P → ∞ both distributions become thermal and give the same
  high-pressure k (p. 12296).
- **Validity:** the final steady state is justified because reactant concentrations do not change on the
  ~2 ms thermal lifetime. The binary-collision model fails above 100 bar.
- **Numerics:** banbks/bandec; tql2 after symmetrization; 10 cm⁻¹ bins.

### A.4 Pfeifle, Olzmann, Int. J. Chem. Kinet. 46, 231 (2014) — PO14 (coupled MEs)

- **Master equation:**
  - Continuous form (eq 1, p. 234): ω∫P(E,ε)n(ε,t)dε − {k_{−a} + k_b + k_c[D]}n.
  - Discretized form (eq 2): dN/dt = R_a F − JN, J = ω(I−P) + K + k_c[D] I.
  - Energy zero: "rovibrational ground state of the intermediate C".
  - Stepladder with detailed balance; 10 cm⁻¹ grains.
- **Sources:**
  - Thermal (eq 7), same form as GO10 eq 15, for E ≥ E₀(−a).
  - Non-thermal A and B (eq 8): f(E) = ∫₀^{E−E₀(−a)} ñ_A(ε) ñ_B(E−E₀(−a)−ε) dε, assuming an
    internal-energy-independent cross section (vibrational PST).
  - Consecutive, partner convolution (eqs 10–11): f2A(E) = ∫₀^{E+RE2a} ñ_O2(ε) ñ1^ss(E+RE2a−ε) dε, with
    RE2a = −79.3 kJ/mol.
  - Shift (eqs 12–13): f2A(E) = ñ1^ss(E + RE2a − ⟨E⟩_O2), "shifted to higher energies by the sum of the
    reaction energy at 0 K and the average thermal energy of O2". Which degrees of freedom ⟨E⟩_O2 covers is
    not specified.
- **Steady state and observables:** N^ss = R_a J⁻¹F, "which can be normalized if required" (eq 5);
  k_j^(ca) = ∫ k_j(E) ñ^ss dE (eq 6). No explicit stabilization observable.
- **Time dependence:**
  - N(t) = R_a Σ e_i [(1−e^{−λ_i t})/λ_i] E_i (eq 3), with F = Σ e_i E_i (eq 4).
  - "The sum of these rate coefficients k_j^(th) … is identical to the lowest eigenvalue of the matrix J".
  - The intermediate→final transition occurs "at times t ∼ 0.1 × τ_therm".
- **Statements:**
  - "the final steady state becomes identical to the intermediate steady state, if a sufficiently fast
    bimolecular reaction occurs" (p. 238).
  - "a thermal master equation approach does not yield the same results … neglects the fact that the
    forming reactions … constantly populate the high-energy levels" (p. 240).
- **Numerics:** tql2 after symmetrization; banbks/bandec.

### A.5 Rabinovitch & Diesen, J. Chem. Phys. 30, 735 (1959)

- **Nascent distribution:** f(E)dE = k_E′K(E)dE/∫_{Emin}^∞ k_E′K(E)dE (eq d). k_E′ is the reverse C–H
  rupture rate and K(E) a Maxwell–Boltzmann factor. Energy from the butyl ground state, Emin = ΔE₀′ + E_H
  (p. 735).
- **Strong collisions:** "each collision leads to deactivation; i.e., k3 = ω" (fn. 4).
- **Rate coefficient:** k_a ≡ ωD/S = ω∫[k_E/(k_E+ω)] f dE/∫[ω/(k_E+ω)] f dE (eq a). Limits ⟨k_E⟩ (eq b) and
  1/⟨1/k_E⟩ (eq c).
- Analytic steady-state competition. The radicals are "closely monoenergetic".

### A.6 Schranz & Nordholm, Chem. Phys. 85, 163 (1984) — SN84

- **Formulation:** ME with a source term (eq 1). The flux enters at E ≥ E_F > E₀ (eq 2). D(t) is eq 3.
- **Stabilization:** S(t) = ∫₀^∞ ∂p/∂t dE (eq 5), i.e. "all reactant molecules" count as stabilized. The
  traditional sink method (eqs 7–9, 17) makes the levels below E₀ absorbing.
- **Collision models:** strong collision (SCA), exponential (EXP), stepladder (SL).
- **Solution:** exact eigen-solution with a constant source switched on at t = 0:
  q_l(t) = (f_l/ωλ_l)[1 − e^{−ωλ_l t}] (eq 30). D(t) and S(t) are eqs 34–35.
- **Intermediate steady state** (eqs 44–50): S ≈ r1 f1 (eq 50); |q_SS⟩ ≈ ω⁻¹S⁻¹(|f̂⟩ − |q_S⟩) (eq 52).
- **Three time scales.** In the final steady state S → 0 (eq 43). The intermediate steady state is "well
  defined only if the eigenvalues … can be separated" (abstract), needing about 3–4 orders of magnitude
  of separation (p. 174). It disappears at high T (Table 6).
- **Sink method errors:** 10–30 % for weak colliders at low P (Tables 9–10, p. 177).
- A δ-pulse of flux is mentioned only as a possible extension (p. 176).

### A.7 Tardy & Rabinovitch, Chem. Rev. 77, 369 (1977) — TR77

- ME with an external input f_i (eqs 8, 11); steady state N^ss = (ωI + k − ωP)⁻¹ f (eq 14).
- Chemical activation: "Collisional quenching to levels below Eo results in stabilization"; k_a = β_c ωD/S
  (p. 396).
- The techniques reviewed are "all of which are steady state" (p. 377). Quasi-steady state is assumed above
  "Eo − Δ" (p. 377).
- The ⟨ΔE⟩ obtained depends on the assumed form of P and on the cross section (p. 396).

### A.8 Barker & Ortiz, Int. J. Chem. Kinet. 33, 246 (2001) — B01 (MultiWell)

- **Nascent distribution:** y0^(ca,i)(E) ∝ k_i(E)ρ(E) exp(−E/k_B T_vib), normalized over E ≥ E₀ (eq 9).
  This is the same ρ·k form as SSUMES.
- **Stochastic solution:** "an exact stochastic algorithm, which Gillespie showed to be mathematically
  equivalent to … the master equation", on a hybrid grained/continuum grid (p. 252).
- Each CA trial runs 0.001 s, "more than enough for collisional stabilization" (p. 253). The reported
  yields are therefore time-integrated after a single injection.
- Converged for ΔE_grain ≤ 50 cm⁻¹.
- `Barker_Multiwell_testcase_ftp.pdf` has identical extracted text.

### A.9 Vereecken, Huyberechts, Peeters, J. Chem. Phys. 106, 6564 (1997) — V97

- **Formation distribution:** P_form(E) ∝ G‡(E−E‡) exp(−(E−E‡)/kT) for E ≥ E‡ (p. 6567). Troe
  biexponential collision model.
- **Sink:** placed 40 kJ/mol below TS2. Lower-energy molecules "can then be collected in a 'sink' without
  significantly altering the calculated product distribution" (p. 6568).
- **Three methods:**
  - ESM: Gillespie stochastic simulation.
  - DCPD: absorbing Markov chain, Fraction_Y = Σ P_form(X(E)) · P_{X(E)→Y}.
  - CSSPI: steady state from the λ = 0 eigenvector, with products re-injected.
- **Result:** "The results of all three methods are identical for all practical purposes" (p. 6573,
  Table I). DCPD is "only applicable in steady-state conditions or for the cumulative product distribution
  of a completed reaction" (p. 6569).
- At 1500 K / 100 atm, thermal reactions of the sink "should be taken into account" (p. 6571).

### A.10 Smith, McEwan, Gilbert, J. Chem. Phys. 90, 4265 (1989) — S89

- **Source term:** k⁻¹(E−ΔH₀¹) f_r(E−ΔH₀¹) A(t)B(t) (eqs 2–4), with microscopic reversibility (eq 3).
  Weak collisions, all states included.
- **Stabilization:** defined spectrally through the lowest eigenvector x1 (g^s, eq 18), not by an energy
  cutoff.
- **Rate coefficients:**
  - k_s = K_eq k¹_uni f_ne (eq 23);
  - k_d^i = k_ss^i − k_s k^i_uni/k_uni (eq 25);
  - Σ k_d^i + k_s = k_cap = k_rec^∞ (eq 26).
- **Eigen-expansion:** ∫ds e^{λ_i(t−s)} AB ≈ AB/|λ_i| (eq 12); Σ_i (q_i/|λ_i|) ψ_i = −B⁻¹u (eq 14).
- g*(t) is "the steady-state population distribution of the collision complex" (eq 19, App. A3).
- Assumes the higher eigenvalues relax very quickly (p. 4267).
- "neglect of activating collisions can lead to significant errors, particularly at lower pressures"
  (p. 4265; also p. 4272).

### A.11 Carstensen & Dean, Comprehensive Chemical Kinetics 42, 105–187 (2007) — CD07

- **CA master equation** (eqs 69/73/76): h(E) = k2(E)f(E)/Σ_{E>E₀} k2(E)f(E) (eqs 74, 106). Energy is
  relative to non-activated AB (p. 115).
- **Modified-strong-collision steady-state apparent rate constants** k_stab,A, k_prodA, k_stab,B, k_prodB
  (eqs 105–109); limits in eqs 110–117.
- **Absorbing boundary:** "sufficiently below the barrier (for example, 10 kT below E₀)" (p. 125).
- **Methods covered:** numerical integration, absorbing-boundary steady state, eigenvalue methods
  (eqs 122–130), and MultiWell.
- **Statements:**
  - "Both methods essentially rely on reaching steady-state concentrations … the smallest negative
    eigenvalues must be clearly separated from those that describe energy relaxation" (p. 125).
  - "Waiting for the energy relaxation to be completed means to assume steady-state conditions" (p. 139).
  - Ethoxy example: "20 % of the decay is not captured" (p. 168).
  - Prompt vs. delayed products are justified by "two completely separated time domains" (p. 176).

### A.12 Pilling & Robertson, Annu. Rev. Phys. Chem. 54, 245 (2003) — PR03

- **Source:** g = k_a,∞[B][C]φ(E), with φ = k f/∫k f dE (eqs 39–42).
- **Absorbing boundary:** "usually placed about 10 k_BT below the reaction threshold".
- **Steady-state density:** ρ_s = −R M̂⁻¹φ (eq 43), "steady-state … in the sense that the relative
  population of each energy state is constant while the absolute value decays". Association coefficient
  in eq 44.
- **Collision kernel:** normalization and detailed balance as in eqs 14–21 (see the detailed-balance
  report).
- **Mixed eigenvalues:** eigenvalue mixing (pp. 262–263). For C2H5 + O2 at 575–750 K, "prompt" vs. thermal
  HO2 with mixed eigenvalues (p. 265).

### A.13 Caralp, Forst, Bergeat, Phys. Chem. Chem. Phys. 10, 5746 (2008) — C08

- Incoming flux F(E,t) = ℜ f(E,t), constant in time, with N(0) = 0 (eqs 12–14). Eigen-expansion (eq 15),
  D(t) (eq 16), k(t) (eq 17).
- Three regimes. The long-time steady state has "a stationary energy distribution … and a zero
  stabilisation flux" (p. 5749), and that distribution is non-thermal (Figs 9–10).
- Earlier steady-state versions "implicitly assumed … a reactant sink slightly below threshold" (p. 5752).

### A.14 Miller & Klippenstein, J. Phys. Chem. A 110, 10528 (2006) — MK06

- **Master equation:** eqs 3a–3c plus a rate equation for n_R (eq 5), under the pseudo-first-order
  condition n_B ≫ n_m ≫ n_R (eq 2). Products are infinite sinks.
- **Enlarged symmetric matrix:** G (eq 12) acting on a vector that includes (n_m/Q_Rm δE)^{1/2} X_R
  (eq 13).
- **Kernel:** P(E,E′) = exp(−ΔE/α)/C_N(E′) for E ≤ E′ (eq 6); "activating wing … from detailed balance".
- **Rate coefficients:**
  - Number of chemical eigenvalues N_chem = S − 1 (eq 19).
  - Initial-rate method (eq 24) needs |λ_Nchem| ≪ |λ_Nchem+1|; long-time method (eq 25) needs only "<"
    (p. 10534).
  - Worked well-skipping examples: C2H5 + O2 → C2H4 + HO2 and allene → propyne.
- **Flux vs. rate coefficient:** "Considerable confusion exists … failure to make a distinction between a
  rate coefficient and what might best be called a 'flux coefficient'" (p. 10529). "Application of the
  steady-state approximation … is virtually always an attempt to equate a phenomenological rate coefficient
  to a flux coefficient. Sometimes this is a valid approach, and sometimes it is not" (p. 10531).
- **Numerics:**
  - At high T: "combine the two (or more) species being equilibrated into one compound species".
  - At low T: "it can be difficult numerically to obtain accurate eigenvalues and eigenvectors …
    (1) quadruple-precision … (2) integrating the ME directly in time" (pp. 10535–10536).

### A.15 Georgievskii, Miller, Burke, Klippenstein, J. Phys. Chem. A 117, 12146 (2013) — G13

- **Master equation:** d|f⟩/dt = −Ĝ|f⟩ + Σ_ν s_ν(t)|p^(ν)⟩ (eq 1), with s_ν = n_A n_B/Q_ν (eq 2) and
  p_i^(ν)(E) = N#_{i,ν}(E) e^{−E/T}/2π (eq 9).
- **State space:** "the dynamical phase space consists of only the microscopic populations of the various
  isomers … bimolecular reactants and products are treated equally as sources and sinks" (abstract).
- **Operator:** self-adjoint under the weight 1/f^(0) (eqs 11–12); detailed balance (eqs 6, 10). The total
  collision rate ω_i(E) is assumed constant (eq 5).
- **Rate coefficients:**
  - Relaxational modes are at quasi-steady state: f_λ = Σ_ν s_ν p_λ^(ν)/Λ_λ (eq 19).
  - Bimolecular→bimolecular: k_{ν→μ} = (1/Q_ν) Σ_relax p_λ^(μ) p_λ^(ν)/Λ_λ (eq 21).
  - Isomerization: k_{j→i} = −√(Q_i/Q_j)(MΛM⁻¹)_ij (eq 27).
  - Bimolecular→isomer: k_{ν→i} = (√Q_i/Q_ν) Σ_CSE M_iλ p_λ^(ν) (eq 28).
  - Isomer→bimolecular: k_{i→ν} = (1/√Q_i) Σ_CSE M⁻¹_λi p_λ^(ν) (eq 30).
  - Capture balance: eqs 22–23.
- **Change from 2006:** does "not require the pseudo-first order assumption … especially important for
  self-reactions" (p. 12147). The current approach "formally corresponds to the limit of zero-concentration
  for the excess bimolecular reactant" (p. 12153).
- **Loss of separation:** merge into "united species" (eqs 31–34); detailed balance holds approximately via
  M⁻¹ ≈ Mᵀ (eq 36).
- **Chemical activation:** branching for arbitrary non-thermal initial distributions (eqs 37, 38, 41;
  Fig. 4). R → well and R → P rate coefficients equal the fractional populations accumulated over the
  relaxation time (p. 12153).
- "results from this code are fully consistent with those derived from our earlier implementation"
  (p. 12153).

### A.16 Johnson & Green, Faraday Discuss. 238, 380 (2022) — JG22

- **Forms:** source form (eq 1) vs. pseudo-first-order augmented form (eq 2). Eq 10 is G13 eq 21 with
  negative eigenvalues, "obtained by explicitly assuming that the relaxational eigenstates are at steady
  state" (p. 383).
- **Accuracy:**
  - "At low temperatures, the accuracy of the eigendecomposition becomes problematic" (p. 401).
  - Matrix condition number about 10¹³–10²¹ (Fig. 10). Neglecting high-energy collisions reduces it, but
    gave no consistent accuracy gain (§3.4).
  - The CSE methods "do not necessarily satisfy equilibrium" (p. 392).
  - No robust method exists for separating eigenvalues that do not separate (p. 402).
  - The cse_g (G13) formulation is more robust at low T than Allen's cse (Fig. 2).
- **Comparison with the time-dependent solution:**
  - With a constant source, the CSE flux error decays on the energy-transfer time scale (eq 29).
  - With a decaying source a residual error remains (eq 40), "can be significant for strongly
    chemically-activated cases".
  - Errors exceed a factor of 10 at 800–1500 K and low P (Fig. 3). With a lowered barrier, CSE fails at most
    conditions (Fig. 4).

### A.17 Zhang, Chen, Truhlar, Xu, Faraday Discuss. 238, 431 (2022) — Z22 (TUMME) and its SI

- **Formulation:**
  - dy/dt = −Wy + B n n (eq 2), symmetrized G = F⁻¹WF (eq 3).
  - Exponential-down kernel normalized by A(E_η′) (eq 13). Quadruple precision.
- **Stabilized fraction:** f(t) = Σ_{E ≤ E₀‡} y_RC/n_OH(0) (eq 27), a time-dependent energy cutoff.
- **Rate constants:**
  - N_CSE = S − δ_MS (eq 22).
  - k_{R→P} = steady-state W⁻¹ term minus the chemically significant mode term (eq 24). When merged, the
    pure W⁻¹ expression applies (eq 25).
  - Merging test P = 1 − EPCS > 0.2 (eq 26), with an "ambiguous zone" for 0.2 < P < 0.6.
  - Table 5: k_exp^SSA/k_exp = 1.00–1.15.
- **SI:**
  - TUMME 3.0 input: N2 bath; E_down 200 cm⁻¹ at 300 K, exponent 0.85; T = 20–1800 K;
    P = 10⁻²–10⁷ Torr; grain 0.2 cm⁻¹.
  - Fig. S3: time profiles from rate constants "match perfectly well" with the full solution.
  - Tables S1/S2: the equilibrium-constant check drifts when the separation fails (e.g. 300 K, 10⁻² Torr:
    1.40·10⁻²² vs 2.84·10⁻²²).
- **Pressure behaviour:** a single-exponential k_exp fails at high p, where the decay is biexponential
  (p. 451).

### A.18 Alexander, Hall, Dagdigian, J. Chem. Educ. 88, 1538 (2011)

Non-reactive ME only: symmetrized (eq 13), with one zero eigenvalue whose eigenvector is the Boltzmann
distribution, and all others negative. Chemical activation is not addressed.

---

## Appendix B — MESMER (manual, version 5.0, 2017) as a fourth reference point

- **Master equation:** one-dimensional in E only, no J resolution (§13.1, p. 112); eq 13.1, p. 113. Products
  are infinite sinks.
- **Kernel:** "only the exponential down model is implemented": P(E|E′) = A(E′) exp(−(E′−E)/⟨ΔE⟩_d)
  (§11.2.2, pp. 93–94). Activating probabilities come from detailed balance. How A(E′) is computed, and how
  normalization is handled at low energy, is **not described**.
- **Temperature dependence:** ⟨ΔE⟩_d(T) = ⟨ΔE⟩_d,ref (T/T_ref)^n, with T_ref = 298 K (eq 11.5).
  ⟨ΔE⟩_d defaults to 130 cm⁻¹ "with warning".
- **Bimolecular source:** a deficient reactant represented as a single thermal grain that **depletes** inside
  the eigenproblem. Block form [[ω(P−I) − K, k_f,∞[B]φ], [k, −k_f,∞[B]]], where "φ is the chemical
  activation distribution" (§13.1.1, pp. 118–120). This is not a constant source flux.
- **Solution methods:**
  - Eigen-decomposition time evolution, n(t) = U e^{Λt} U⁻¹ n(0) (eq 9.4).
  - Yields P(t) (eqs 9.5–9.6).
  - Bartis–Widom rate coefficients.
  - **No steady-state chemical-activation solver; no stabilization fraction; no "prompt" or well-skipping
    treatment.**
- **Validity:** warns when chemically significant and relaxation eigenvalues are "not well separated by more
  than an order of magnitude" (p. 64).
- **Precision and size reduction:** double-double / quad-double arithmetic (p. 41, §8.4). Reservoir states
  lump the low-energy grains (§13.2.1).
- **Non-thermal starts:** Boltzmann at another temperature, or Prior (§11.2.6). Consecutive activation via
  pseudo-isomerization with a Prior fragmentation distribution (eqs 11.17–11.19).

---

## Appendix C — SSUMES source code, claim by claim

`source/` = `/home/peter/Programs/SSUMES/ssumes/source/`. The CA path is `carate.cc` → `umolProb` in
`ssulibc.cc`. The Fortran in `unimol/` (UNIMOL) is not used by the CA solver.

| # | Item | Code |
|---|---|---|
| V1 | α(T) = alpha1000·(T/1000)^alpTex; Z = πσ²·10⁻¹⁶ · CVMEAN√(T/μ) · p/(RTCMMLC·T) · Ω22, Ω22 = (0.636 + 0.567 log10(T/ε))⁻¹. **CVMEAN = 14550.8069** cm s⁻¹ K⁻½ amu½ (= 100√(8k_B/π·amu)); **RTCMMLC = 1.03557154·10⁻¹⁹** Torr cm³ K⁻¹ (= k_B) | `ssulibc.cc:172–177`; `ssulibb.h:65, 69` |
| V2 | Exponential down, normalized top-down: norm[jc] = dnorm/(1 − usum); up probabilities use the target's norm; reduction factors redfac below ncutlow (ρ-gradient threshold thRhoGrad = 3); redExpon automatically 1.0 → 3.05 in steps of 0.1; band noff = ⌊ln(err3)/ln β⌋ + 1, elements outside the band zeroed. **Detailed balance exact** (redfac always on the lower grain) | `ssulibc.cc:1015–1064, 991–1005, 1351–1358` |
| V3 | Stabilization leakage d_j = Z Σ_{i<lowest} P(i←j); k_stab = Σ g_i d_i/G. Default `truncate 0` = no stabilization | `ssulibc.cc:1069–1073, 1414–1422` |
| V4 | Diagonal loss of all channels; internal channels as off-diagonal terms; kEthOut zeroes small outgoing k(E); kEthInt zeroes both directions; detailed-balance correction corf = √((ρ_u/ρ_w)/(k_wu/k_uw)) if the deviation exceeds 1 % (warning above 10 %, stop above 60 %) | `ssulibc.cc:1093–1095, 1150, 504–517, 450–478` |
| V5 | Storage b[source][target]; symmetMat b[jc][ir] *= s[jc]/s[ir], i.e. W⁻¹MW in target-row form; symmetry check: error above 5 %, warning above 1 % | `ssulibc.cc:1285–1302` |
| V6 | CA source r_i = ρ_i γ^i k_back(E_i), normalized on [lowest, upbound); `recombChan` is the adduct's back-dissociation channel; alternatives `excitGrain` and `excitArbit` | `ssulibc.cc:1174–1190`; `carate.cc:87` |
| V7 | carate: DPOSV (Cholesky); dislit: DPOTRS iterations; catime: DLSODES (mf = 121, rtol = 1e−6, atol = 1e−10); diseig: DSYEVR (top 20 % of eigenvalues) | `ssulibc.cc:740, 758, 772–803, 697` |
| V8 | k_ch = Σ g_i k_ch(E_i)/G with G summed over all wells; k_out,tot = outgoing channels + Σ k_stab; fractions k/k_out,tot are printed. The per-well fractions are computed but never printed | `ssulibc.cc:1387–1443`; `carate.cc:269–271` |

- **Modes:**
  - `dislit` iterates a source built from the normalized steady population until k_out,tot changes by
    less than 1e−5.
  - `diseig` uses the eigenvector of the largest negative eigenvalue, requiring non-negative populations
    and a reactant-well population ≥ 0.5.
- **Not present:** re-formation (re-association), consecutive activation, more than one source
  (`ssulibc.cc:309–312, 335–336`).
- **Problems found:**
  - untracked top loss (`:1047–1051`);
  - possible out-of-bounds write in `checkEqkE` (`:450–454`);
  - memory leak `delete[] norm, redfac;` (`:1058, 1079`);
  - the manual's sizeLSDSRW default (0.3) differs from the code's (0.1) (`:533`).

---

## Appendix D — MESS 2026 source code: modes and rate-coefficient extraction

Paths relative to `/home/peter/Dropbox/Research_Leuven/MESS_kinetics/Source_from_2026/MESS/`.

- **Build:** `mess` is built from `src/mess_driver.cc` + `src/libmess/mess.cc` (`CMakeLists.txt:314,
  351`). `new_mess.cc` is not compiled.
- **Methods** (`mess_driver.cc:931–936`, default direct at `:140`):
  - **direct:** full diagonalization of the symmetrized well-only matrix (`mess.cc:10867, 11224–11234`).
  - **low-eigenvalue:** Schur complement onto the thermal well vectors (`mess.cc:7384–7390, 7399`).
  - **well-reduction:** fast "horizontal relaxation" modes are removed per energy bin and used only for
    bimolecular→bimolecular contributions (`mess.cc:8304–8306, 8323, 8476–8571`; threshold 10 × collision
    frequency, `mess.cc:122`).
- **HotEnergies:** branching fractions of hot isomers at given energies (`mess_driver.cc:1245–1360`;
  `mess.cc:6227, 13369–13429`).
- **TimeEvolution:**
  - Analytic from the eigen-decomposition; no numerical integration.
  - The bimolecular reactant is a constant, non-depleting source scaled by `ExcessReactantConcentration`
    (`mess.cc:11597–11808, 11657–11685`).
  - Alternatively, a bound reactant starts thermal at `EffectiveTemperature` (`mess.cc:11726`).
- **Bimolecular entrance and exit:** one vector per well–bimolecular barrier,
  N‡(E) exp(−(E−E_ref)/kT)/2π/√(ρe^{−E/kT}) (`mess.cc:11127`), used for both entrance and exit.
- **Number of chemical eigenvalues:** `chemical_threshold = −2` (gap ratio) by default
  (`mess.cc:114, 12092–12141`). If no gap is found, only bimolecular→bimolecular rates are produced.
- **Rate formulas:**
  - bimolecular→bimolecular: Σ_relax ⟨b_i|l⟩⟨l|b_j⟩/λ_l (`mess.cc:12262`);
  - well→well: from M Λ M⁻¹ (`mess.cc:12934`);
  - well→bimolecular and bimolecular→well (`mess.cc:12955, 12965`).
- **Merging:** wells merge via `threshold_well_partition` (`mess.cc:146, 12842`). Wells left in the
  "bimolecular group" are treated as equilibrated with the products and printed as `***`.
- **Steady-state history:** the steady-state Cholesky version of bb_rate is commented out (`mess.cc:12193,
  12257`). The fallback described in the manual for `ChemicalEigenvalueMin` is commented out
  (`mess.cc:11321–11565`).
- **Kernel:** pairwise normalization, detailed balance exact (`mess.cc:6566–6614`). Optional
  `EnergyRelaxationFlags` switch to a row-wise constant-collision-rate normalization
  (`mess.cc:6444–6549`).
- **Precision:** double by default (`mpack.cc:18`). Multiple precision needs WITH_MPACK, which is not set in
  CMake (`mess_driver.cc:1046–1050`), plus `UseMultiPrecision`. Results are cast back to double
  (`mpack.hh:182, 190`).
- **Escape:** pseudo-first-order loss of hot isomers on the diagonal (`model.cc:28679`; `mess.cc:11060`).
  The global `Reactant` keyword is only an energy reference (`model.hh:2425`).
- **Outputs:** "prompt isomerization/dissociation" table (`mess.cc:13436–13556`); PEDs of bimolecular
  products are called "Product steady-state energy distributions … in a chemically activated process"
  (manual §5.5).

---

## Appendix E — Chemical-activation observables defined in the literature

What a code should be able to output (source: paper and equation):

1. **Nascent distribution** f(E) / h(E) / φ(E): Rabinovitch eq d; O91 eq 12; GO10 eq 15; PO14 eqs 7–8;
   B01 eq 9; CD07 eqs 74, 106; PR03 eq 42; V97 P_form; S89 eq 4.
2. **D/S and k_a = ωD/S** (strong collisions, with β_c for weak collisions): Rabinovitch eqs a–c; TR77
   p. 396; SN84 eq 40.
3. **Time-resolved fluxes** D(t), S(t), with total flux = D + S: SN84 eqs 3–5, 34–35. Intermediate steady
   state S ≈ r1 f1: SN84 eq 50. C08 D(t): eq 16.
4. **Yields per molecule formed** Φ_r = Σ_i (K_r J⁻¹F)_i: O91 eq 5; O02 eq 10. Fractions k/k_out,tot:
   SSUMES.
5. **Chemically activated rate coefficients** k_r^ca = Σ_i (K_r Ñ^s)_i: GO10 eq 9; PO14 eq 6. Branching
   including a bimolecular sink: GO10 eq 10.
6. **Apparent CA rate constants** k_stab,i, k_prod,i: CD07 eqs 105–109. k_s, k_d^i, k_cap: S89 eqs 23,
   25, 26.
7. **Absorbing-boundary steady-state density** and association coefficient: PR03 eqs 43–44.
8. **Energy-resolved branching** P_{X(E)→Y} and Fraction_Y: V97 Sec. III.
9. **Thermal rate coefficient** k^th = lowest eigenvalue of J: GO10 eq 12.
10. **Phenomenological rate coefficients including well-skipping:** G13 eqs 21, 27, 28, 30; Z22 eqs
    23–25; CD07 eq 129; PR03 eq 47.
11. **Time-dependent stabilized fraction** f(t): Z22 eq 27.
12. **Separation diagnostics:** eigenvalue ratios, flux overlaps r_l f_l/F (SN84 Tables 2–5), merge
    criterion P = 1 − EPCS (Z22 eq 26).
13. **Prompt fraction:** C08 p. 5752; CD07 p. 176.
