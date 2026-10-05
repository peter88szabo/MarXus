# Eigen-solvers for the master equation: LAPACK import from Behemoth, and assessment of Behemoth's Davidson solvers

**MarXus, 2026-10-05.**

## The request

Peter asked:
- "I do not know that we need all the eigenvalues with eigensolver, if not then a Lanczos Krylov type solver maybe better, for those check Behemoth quantum chemical code source code, there is a regular for BSE and TDDFT and also a Pulay variant with smaller memory footprint), import those if you can use them";
- "if you need the Lapack/Blas interface then use that one from Behemoth too ... and put it into our numerical module".

Behemoth was only read; nothing in Behemoth was changed.

## 1. What the eigenvalue analysis needs

| quantity | eigenpairs needed | used for |
|---|---|---|
| k_uni = λ₁ (GO10 eq. 12), k_j^th, sum rule | lowest one, with high *relative* accuracy although λ₁ can lie many orders of magnitude below the largest elements of S (ω + k(E) at the top of the grid) | every condition |
| λ₂/λ₁ (separation of time scales) | second lowest, moderate accuracy | diagnostic |
| N(t) = R Σ eᵢ(1 − e^{−λᵢt})/λᵢ Eᵢ (PO14 eqs. 3–4), λ_F (O02) | **all** eigenpairs that overlap the source F; a chemically activated F at high energy overlaps essentially the whole spectrum | time dependence, O02 window |

The matrices are symmetric (after S = D⁻¹JD) and banded, with n = 500–3000 grains per well in the C₂H₃ decks and a bandwidth of a few hundred grains.

## 2. Measured cost and accuracy (C₂H₃, full deck, n = 2634)

Two conditions (300 K and 2000 K at 1 atm), release build:

| solver | time | memory | sum-rule deviation at 2000 K |
|---|---|---|---|
| inverse iteration, banded Cholesky (default) | 6.3 s | 137 MB | 5.4·10⁻¹³ |
| LAPACK DSYEVD (new) | 5.0 s | 391 MB | 2.4·10⁻¹¹ |
| in-house Householder/QL (tred2/tql2) | 238 s | 269 MB | 2.1·10⁻⁹ |

Over all 40 conditions of the deck (`validation/c2h3_mess_example_olzmann_eigen/`):

| T (K) | sum-rule deviation, inverse iteration | sum-rule deviation, LAPACK |
|---|---|---|
| 750 | 4·10⁻⁹–7·10⁻⁸ | 5·10⁻⁸–10⁻⁵ |
| 1000 | 3·10⁻¹²–4·10⁻¹¹ | 2·10⁻¹⁰–7·10⁻⁸ |
| 500 | sum rule within 1.5% at 5 of 5 pressures | within 1.5% at 1 of 5 pressures (λ₁ < 0 at 10 atm) |

**Why.**
- The dense decomposition has a normwise backward error of order ε‖S‖.
- The banded Cholesky factor has a componentwise one (Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., ch. 10), which respects the grading of S.

**Conclusions.**
- **Inverse iteration remains the default for λ₁.**
- **LAPACK DSYEVD is the route to the full spectrum.** It is about 50 times faster than the in-house QL at n = 2634.

**Update (Peter's decision, 20:01): k_uni is now the eigenvector average, and the sum rule only warns (threshold 1.5%).**
- For C₂H₃, k_uni agrees between inverse iteration and LAPACK to 7 digits at 300 and 500 K, even where LAPACK's λ₁ is wrong or negative.
- LAPACK therefore also covers the conditions where the inverse iteration has no Cholesky factor (300 K at 0.1, 3, 10 atm). The error message of the inverse iteration says so.
- On an extreme synthetic well (175 K, threshold 6000 cm⁻¹, stepladder), however, LAPACK's k_uni differed by 10% from the inverse iteration. Part of this is the reducible stepladder J (validity report §4.7). Inverse iteration stays the reference.

**Update (Peter's go-ahead, 20:18): shifted factorization.**
- **What changed.** The inverse iteration now factors S + σI, σ = n·ε·max Sᵢᵢ (`cholesky_safety_shift`; NR92 §11.7, Higham ch. 10, Demmel LAWN 14).
- **Result.** It solves all 40 C₂H₃ conditions, including 300 K, where the plain factor did not exist at 3 pressures. k_uni agrees with LAPACK to 10⁻⁷ … 10⁻⁴ there.
- **Limitation.** It cannot separate eigenvalues closer than σ; this happens with the stepladder on fine grains, where J is reducible.
- **Resources.** All runs now use at most 4 cores (`OPENBLAS_NUM_THREADS=4`, `cargo -j 4`; Peter's rule). The validation run takes 178 s wall and 394 s CPU.

## 3. What was imported: the LAPACK interface

**New file `src/numeric/lapack_interface.rs`** (declared in `numeric/mod.rs`, declaration line only).
- **Function.** `symmetric_eigen_lapack(matrix: &[Vec<f64>]) -> Result<(Vec<f64>, Vec<Vec<f64>>), String>`. It has the same contract as `symmetric_eigen::symmetric_eigen`: ascending eigenvalues, eigenvectors as rows, only the lower triangle read.
- **The binding is Behemoth's** (`Behemoth/src/numeric/linalg.rs`, `symmetric_eigen_blas`). It is a direct call of the Fortran symbol `dsyevd_` through `extern "C"`:
  - with a workspace query (LWORK = LIWORK = −1);
  - with the row-major/column-major identity of a symmetric matrix (UPLO = 'U' on the row-major data reads the lower triangle);
  - with explicit `info < 0` / `info > 0` errors.

  It is adapted from ndarray to `Vec<Vec<f64>>` and `String` errors, the conventions of MarXus.
- **Algorithm reference.** DSYEVD = Householder reduction + divide and conquer: Cuppen, Numer. Math. 36, 177 (1981); Gu, Eisenstat, SIAM J. Matrix Anal. Appl. 16, 172 (1995); LAPACK Users' Guide, 3rd ed. (1999), §2.4.4.

**Linking.**
- **The library.** The system OpenBLAS is used: `/usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so`, which exports `dsyevd_`, `dsbevx_`, `dpbtrf_` and `dstebz_` (checked with `nm -D`).
- **How it is linked.** `#[link(name = "openblas")]` on the `extern` block, behind a cargo feature `openblas`, which is on by default. `Cargo.toml` now has a `[features]` section (`default = ["openblas"]`, `openblas = []`) and **still no dependencies**.
- **Why not the openblas-src crate.**
  - Behemoth links through ndarray-linalg → openblas-src 0.10.13 (system).
  - The current openblas-src (0.10.16) requires a TLS backend and pulls about 48 crates for *downloading and building* OpenBLAS from source. They are not in the offline cargo cache, and are not needed for a system library.
  - Pinning 0.10.13 still pulled newer transitive crates.
  - With the `system` feature, openblas-src only emits `rustc-link-lib=openblas` on Linux, which the attribute does directly.
- **Why an attribute and not a build script.** A build-script link line reaches only the library target. `src/main.rs` compiles its own copy of the `numeric` tree, so its tests failed to link (`undefined symbol: dsyevd_`). The attribute applies wherever the module is compiled.
- **Without the feature** (`cargo build --no-default-features`), `symmetric_eigen_lapack` returns an error naming the feature. There is **no silent fallback** to another algorithm.

**Wiring.**
- `EigenSolver::FullDecompositionLapack` in `chemical_activation_eigen.rs`.
- `EigenSystem::new(op, solver)` now takes the solver (`FullDecomposition` or `FullDecompositionLapack`). `InverseIteration` is an error there, because N(t) needs all eigenpairs.
- Example option `--eigen-solver lapack`.

**Tests** (TDD; all green in both build modes):
- `lapack_interface::divide_and_conquer_agrees_with_the_householder_ql_decomposition` (n = 1, 2, 7, 60, 150; eigenvalues 10⁻¹², A v = λ v, overlap with the QL vectors);
- `lapack_interface::only_the_lower_triangle_is_read`;
- `lapack_interface::a_non_square_matrix_is_an_error`;
- `lapack_interface::without_the_openblas_feature_the_lapack_route_is_an_error` (compiled only without the feature);
- `chemical_activation_eigen::the_lapack_decomposition_agrees_with_the_other_solvers`;
- `chemical_activation_eigen::the_time_dependent_solution_needs_a_full_decomposition`.

**Full suite:**

| build | lib tests | bin tests |
|---|---|---|
| default (OpenBLAS) | 134 | 26 |
| `--no-default-features` | 132 | 25 |

**Not imported from Behemoth's `linalg.rs` (1779 lines).**
- The ndarray-based `LinAlg` backend (Backend::{Auto, Blas, Pure}, matmul, Cholesky, triangular solves).
- `FactorizedGeneralizedEigh` (DSYGST + DSYEVD + DTRSM, for SCF).
- The SCF commutator helpers.

MarXus works on `Vec` data and has no dense SCF-type problems. These would bring ndarray, ndarray-linalg, cblas-sys and anyhow as dependencies for no present use. MarXus's `numeric/krylov.rs` already contains the BiCGSTAB/GMRES solvers.

## 4. Behemoth's Davidson solvers: assessed, not imported

**What they are.**
- **`Behemoth/src/numeric/davidson.rs`** (664 lines). Davidson's method (E. R. Davidson, J. Comput. Phys. 17, 87 (1975)) for the lowest eigenpairs of a large symmetric operator:
  - the operator enters through the `DavidsonOperator` trait (dimension, diagonal, apply, precondition, apply_block);
  - corrections are the diagonally preconditioned residual or Olsen's correction;
  - the subspace eigenproblem is solved with Behemoth's `LinAlg`.
- **`Behemoth/src/numeric/davidson_lowmem.rs`** (335 lines). The space-saving variant of van Lenthe and Pulay (J. Comput. Chem. 11, 1164 (1990)):
  - three vectors per root;
  - roots one at a time with deflation;
  - overlap floor 10⁻¹².
- **Both depend on** ndarray, anyhow and `LinAlg`.

**Why they are not imported now.**

1. **λ₁ needs relative accuracy that Davidson cannot give.**
   - Davidson stops on the residual ‖Sx − θx‖. The Ritz value θ then has an absolute error of at least order ε‖S‖, the same floor as any backward-stable method.
   - For C₂H₃ at 500 K, λ₁ ≈ 10⁻⁴ s⁻¹ is already at the resolution limit: the sum-rule deviation is 10⁻³ even with the most accurate solver. ‖S‖ itself was not evaluated.
   - Inverse iteration with the exact banded Cholesky factor converges per step with the ratio λ₁/λ₂: ≤ 1.5·10⁻⁴ at ≤ 1000 K, ≤ 0.1 at 2000 K. It has the more favourable componentwise error (§2).
   - Davidson's convergence for the lowest ME mode with a diagonal preconditioner was **not tested**. No claim is made about its rate.
2. **The full spectrum is cheap at ME sizes.** Davidson and Lanczos pay off when the dimension forbids a dense decomposition, as for TDDFT/BSE response matrices (10⁵–10⁶). For the ME (n ≤ a few thousand per well, banded), LAPACK DSYEVD gives all eigenpairs in about 2.5 s per condition at n = 2634. N(t) and λ_F need essentially all of them.
3. **λ₂ is only a diagnostic.** The deflated inverse iteration returns NaN when λ₂/λ₁ exceeds about 10¹⁶ (300 K). That is reported, and the full decompositions give λ₂ whenever needed.

**When they would become useful.** In a large multiwell analysis (many wells, n ~ 10⁴–10⁵), only the lowest N_wells + 1 eigenpairs and the gap are needed (chemically significant eigenvalues). An iterative solver for a few eigenpairs would then be appropriate.

Because of point 1 it should work on S⁻¹, applied through the banded Cholesky factor:
- shift-and-invert Lanczos;
- block inverse iteration;
- or Davidson with S⁻¹ as preconditioner.

Davidson on S itself would not do. At that point the Davidson code (DavidsonOperator trait, van Lenthe–Pulay storage scheme) can be ported to `Vec` data in `src/numeric/`, with the operator applying S⁻¹.

**A LAPACK alternative for selected eigenpairs.** The banded routines DSBEVX (selected eigenpairs of a band matrix: DSBTRD reduction + DSTEBZ bisection + DSTEIN) and DPBTRF (banded Cholesky) are present in the system OpenBLAS. They would give λ₁…λ_k of the band matrix without dense storage. They are not bound yet.

## 5. Files changed

- **New:** `src/numeric/lapack_interface.rs`, this report.
- **Modified:**
  - `src/numeric/mod.rs` (one `pub mod` line);
  - `Cargo.toml` (`[features]`);
  - `src/masterequation/chemical_activation_eigen.rs` (`FullDecompositionLapack`, `EigenSystem::new(op, solver)`, sum-rule check);
  - `src/masterequation/chemical_activation_driver.rs` (tolerance argument, table columns);
  - `examples/chemical_activation_from_deck.rs` (`--eigen-solver lapack`, `--sum-rule-tolerance`).

Nothing is committed.
