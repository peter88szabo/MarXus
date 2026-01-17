# SACM Port Audit

## Findings vs `barrierless/sacm/SACM_mod.f`

### 1) Missing Faminf baseline for angular momentum corrections (major)
- **Fortran:** `WPST(I)` is initialized by `WCOUNT/SCOUNT` (Faminf), then blended with `FAME`:
  `WPST(I)=WPST(I)+(FAMES-WPST(I))*EXP(-C3*ALBEEF)`.
- **Rust:** `apply_corrections` starts from `1.0` and blends `FAME` against that, so Faminf is skipped entirely.
- **Impact:** SACM correction factors are too simplified; results differ when `alpha_over_beta > 0`.

### 2) Reactant density corrections missing (major)
- **Fortran:** Reactant density `RHO(E)` is multiplied by `FANHRO` and `RRHIND` before the J‑resolved densities are built.
- **Rust:** No use of `anharmonic_reactant_rhoE` or `hindered_rotor_density_factor` on the denominator.
- **Impact:** `k(E,J)` denominator is incorrect compared to the Fortran logic.

### 3) Thermal rate expression differs (major)
- **Fortran:** High‑pressure rates use QCENT, QSTARP, QVIB/QIR/QEL partition functions, and `QSTAR` products.
- **Rust:** Uses a simple Boltzmann average of `k(E)` with `rho(E)` weights.
- **Impact:** Thermal rates are not equivalent to the Fortran outputs.

### 4) Fam energy shift inputs are manual (medium)
- **Fortran:** `EIF = E + AWR(2E/SEPSE, BETWR) * SEPSE/2` with `SEPSE`/`BETWR` derived from reactant vibrational data.
- **Rust:** Requires `reactant_freq_sum` and `wr_beta` as manual inputs.
- **Impact:** Extra manual input and possible mismatch if not derived consistently.

### 5) PST construction simplified (medium)
- **Fortran:** `SCOUNT/WCOUNT` handle disappearing/conserved degrees of freedom and K‑rotor treatment.
- **Rust:** Uses convolution of RRKM states for two independent fragments.
- **Impact:** PST baseline differs in detail from Fortran unless the simplified model is intentionally accepted.

## Notes
- These findings reflect the current Rust implementation as of the most recent changes in `src/barrierless/sacm`.
- The differences may be acceptable if the simplified PST and correction models are desired; otherwise the missing pieces should be implemented.
