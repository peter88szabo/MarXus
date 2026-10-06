# pTDME: what MarXus could use as alternative options

**Date:** 2026-10-06. **Request (Peter):** "check if we can use anything from this: /home/peter/Programs/pTDME/pTDME-To-Multiwell-Aug-30-2022 (other collision models, other anything that can be useful as alternative option)".

**How this was done.** A read-only survey of the source, its manual (`/home/peter/Programs/pTDME/Manual-pTDME-Feb-1-2023-V2.pdf`) and its presentation (`pTDME-presentation-Feb-2-2023.pptx`). Nothing was modified. File paths below are relative to the code root.

## 1. What pTDME is

**Code.** Fortran 90 with MPI, by T. L. Nguyen and J. F. Stanton (main header dated 2021-02-04).
- **Method.** An E,J-resolved master equation in the fixed-J approximation: for every total angular momentum J an independent 1D energy-grained ME. J is conserved in collisions.
- **Results.** Combined as k = Σ_J w_J λ(J), with w_J from the formation distribution F(E,J).

**References in the package.** The manual asks users to cite J. R. Barker et al., MultiWell (2023), and T. L. Nguyen, J. F. Stanton, "Pragmatic Solution for a Fully E,J-Master Equation", J. Phys. Chem. A 124, 2907 (2020).
- In the code itself there are almost no references: "Reid–Prausnitz–Sherwood, 1977" for the collision integral (`CalCF.f90:53`), and ARPACK's technical report.
- **No equation numbers are cited anywhere.**

**License: none found.** There is no LICENSE/COPYING file and no license text in the headers; only "Purpose: Education and Research should be encouraged" (`pTDME.f90:18`). The bundled ARPACK has no license text either.
- **Consequence for MarXus:** no pTDME code may be copied. Anything adopted must be implemented from the papers. The authors could be asked (addresses in the manual).

## 2. Collision model

**Only one kernel: single exponential down.** It acts on grain-index differences, with one constant ⟨ΔE_down⟩ for all wells, independent of T and E (`Pdij.f90:43`: `exp(-(J-I)*dE/AL)/CN(J)`).
- Up transitions come by detailed balance with the J-resolved densities (`Puij.f90:49`).
- There is no biexponential, Gaussian, stepladder, Troe-type, ⟨ΔE⟩_all or α(E) form.
- **Nothing new for MarXus,** which has exponential down with ⟨ΔE_down⟩(T) and the stepladder.

**Normalization: the same top-down back substitution as MarXus** (Robertson 2019, eq. 4.16): `CalCN.f90:46-66`, CN(J) = su1/(1 − su2).

**Where it fails** (su2 > 1), pTDME has a fourth low-energy rule, besides the reduction, truncation and reservoir rules considered for MarXus:
- `IF(CN(J).LT.0.0d0) CN(J) = 0.0d0`: that grain loses all its downward and elastic transitions, and keeps its upward ones.
- The diagonal is then built from the actual column sums (`CalA.f90:41-55`; residuals subtracted in `CALRATE_ASYM.f90:440-451`), so population is conserved.
- The clamped grain gets an effective collision rate ω·su2 > ω.
- su2 = 1 exactly is not guarded (division by zero).

Assessment for MarXus: it keeps all states, but it changes the collision rate of the affected grains, and nothing documents it physically. It is not better than the reservoir state (`reports/low_energy_reservoir_state.md`). **Not recommended.**

**A J-dependent ⟨ΔE_down⟩ exists in the code but is disabled** (`pTDME.f90:340-352`): ⟨ΔE⟩_down = ⟨ΔE⟩_vib + ⟨B⟩·ΔJ·(2J + 1 − ΔJ) ("spherical top"). Its coefficient is never read (NJJ = 0).
- **Possible later option** for a J-resolved MarXus master equation, but only from a paper with the derivation. The code cites none.

## 3. Collision frequency: a possible alternative option

`CalFC` in `CalCF.f90` offers two collision integrals Ω(2,2)*:
- **TROE:** 1/(0.636 + 0.246 ln T*). This is numerically MarXus's current formula, 1/(0.636 + 0.567 log₁₀ T*).
- **RPS** (the default): 1.16145/T*^0.14874 + 0.52487 e^{−0.7732 T*} + 2.16178 e^{−2.43787 T*}. The code cites only "Reid–Prausnitz–Sherwood, 1977".
  - The survey recognized these coefficients as the Neufeld–Janzen–Aziz fit. That attribution is ours, not the code's, and must be checked in the paper before use.

**Candidate option for MarXus: the RPS/Neufeld form as an alternative collision integral.** It is more accurate than Troe's two-parameter fit over a wide T* range; Troe's fit is stated to be good to ±7% (MarXus `collisional_relaxation.rs`). **Needed before implementing:** the paper (Neufeld, Janzen, Aziz, J. Chem. Phys. 57, 1100 (1972), to be confirmed) or the Reid–Prausnitz–Sherwood book page.

## 4. Solution methods

**The variants:**
- **Full eigen-expansion in time** (ASYM + LAPACK, `CALRATE_ASYM.f90:522-698`): DGEEV of the non-symmetric matrix, then C(t) = Σ_j v_j e^{λ_j t}(V⁻¹C₀)_j. It has eigenvalue clean-ups by thresholds (|λ| below thresholds set to 0, positive λ set to 0).
- **Shift-invert ARPACK** for a few eigenvalues near a guess (`Call_ANEV.f`).
- **A symmetric variant** with DSYEV.

**None of these is new for MarXus**, which already has the full symmetric decomposition (LAPACK DSYEVD and Householder/QL), inverse iteration for the thermal eigenpair, and direct time integration with Rosenbrock methods. The eigenvector expansion of N(t) is used in a MarXus test (PO14 eqs. 3–4).

**Not present in pTDME:** steady-state solvers, the absorbing-barrier steady state, ODE integrators, extended precision (its `qp` kind is ordinary double precision).

## 5. Densities and sums of states, rotations, tunneling

**Vibrations.** Harmonic direct count with every oscillator level rounded to a grain (`HO.f90`). There is no anharmonicity and no hindered rotors. **Nothing new:** MarXus counts states on 1 cm⁻¹ cells.

**J-resolved densities** with a symmetric-top rotational energy E_rot = √(BC)·J(J+1) + (A − √(BC))K² (`Er.f90`, `cdsj*.f90`).
- The mode is hard-wired to "INACTIVE" (`pTDME.f90:175`).
- **This is relevant only if MarXus gets an E,J-resolved master equation** (Section 7).

**Variational RRKM** along a splined bond-length path (`driver_VRRKM.f90`): minimum of N(E,J) over the path points, with an optional free rotor. The linear and nonlinear versions differ by a factor 2 in the free-rotor term, an inconsistency in the code. MarXus treats barrierless channels with ILT and PST (and SACM in progress).

**Tunneling:**
- **Eckart** with the Johnston–Heicklen-type transmission (no citation in the code). MarXus already has the exact Eckart transmission (Miller 1979, eq. 8).
- **A Kemble/WKB transmission** on a splined 1D path. Its input reading looks inconsistent and untested (`pTDME.f90:510-523`).
- **Nothing to adopt.** MarXus decided on exact Eckart only (no semiclassical P).

## 6. Chemical activation

**CAT1** ("Barker"): start in the well with F(E,J) = (2J+1) N‡(E,J) e^{−βE}; k = k∞ × yield. This is the MarXus chemical-activation route.

**CAT2** ("Pilling"): a pseudo-first-order reactant state feeds the well grains (`CALRATE_ASYM.f90:423-436`), and k = −λ₁/[C]. In pTDME:
- [C] is hard-coded to 10¹⁵ cm⁻³;
- detailed balance between the formation and the redissociation is not enforced;
- the entrance must be "product 1 of well 1".

**Assessment:** the idea (the reactant as a state of the master equation, the bimolecular eigenvalue route) is G13's formulation. MarXus's CSE uses it with exact detailed balance. **Nothing to adopt from this implementation.**

## 7. What is useful for MarXus

| item | use | status |
|---|---|---|
| **Benchmark examples with stored outputs** (`EXAMPLES/THERMAL_DECOMPOSITION/{C2H5, C2H5_CSE, NH3, NH3_CSE, CH3O, CH3O_2wells}`, `EXAMPLES/CHEMICAL_ACTIVATION/{H+C2H4 CAT1/CAT2/CAT2_CSE, OH+CO CAT1/CAT2, OH+HNO3 CAT1/CAT2}`, with `pTDME.out`, populations and rates per pressure) | further validation systems for MarXus: thermal decomposition (one and two wells) and chemical activation; e.g. C₂H₅ at 1000 K, 760 Torr: λ = −1.319·10⁴ s⁻¹ | **recommended**; their decks must be translated to the MESS format MarXus reads, and the E,J (fixed-J) vs E-only difference must be kept in mind |
| Neufeld-type collision integral Ω(2,2)* | alternative collision-frequency option | candidate; needs the paper |
| E,J-resolved master equation (fixed-J approximation, Nguyen & Stanton 2020) | a possible later MarXus extension | needs the paper (JPCA 124, 2907) and a decision |
| J-dependent ⟨ΔE_down⟩ | only with a J-resolved ME | disabled in pTDME; needs a paper |
| CN clamp at low energy | — | not recommended (Section 2) |
| everything else (single exponential down, eigen-expansion, ARPACK, harmonic counting, Eckart, WKB, CAT1/CAT2) | — | already in MarXus in an equivalent or better form, or decided against |

**Quirks of pTDME not to copy:**
- Time-array sizing can write one point past the end.
- `RKinf` is used before being set (harmless).
- F(E,J) is uninitialized in one LOOSE branch.
- The external F(E,J) path reads an uninitialized index.

## 8. Recommendation

1. **Use the pTDME examples as additional validation cases** (thermal decomposition C₂H₅, CH₃O with two wells; chemical activation H + C₂H₄), after translating the decks.
2. **Consider the Neufeld-type Ω(2,2)\* as an option of the collision frequency**, from the paper.
3. **E,J-resolved master equation:** a later decision, from Nguyen & Stanton (2020).

No code is copied; there is no license.
