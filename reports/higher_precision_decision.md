# Higher precision for the master-equation eigenvalue analysis: decision and plan

**MarXus, 2026-10-05.**

Peter's requirements (20:28–20:45):
- higher precision than double, starting with double-double, with a higher-precision reference;
- assemble the sensitive matrix operations in the selected precision;
- symmetric solvers where detailed balance permits;
- convergence checks as the precision increases;
- no built-in Rust f128 in production;
- **easy installation**: only crates that cargo installs by itself, or third-party sources vendored inside MarXus and compiled with it. The user never downloads, builds or maintains a dependency by hand. The system OpenBLAS is an accepted exception.

## 1. Why double precision is not enough

In the C₂H₃ benchmark the double-precision floor ε·max Sᵢᵢ is 1.0·10⁻² s⁻¹. At 300 K the thermal eigenvalue is λ₁ ≈ 10⁻¹⁵ s⁻¹, so the computed λ₁ is rounding noise, sometimes negative.

The reported rate coefficient k_uni, the eigenvector average of GO10, is still correct to 10⁻⁷ … 10⁻⁴ between solvers. It also agrees with the absorbing-barrier route to 0.001 percentage points. But:
- the sum rule λ₁ = k_uni can no longer serve as a check at low T;
- λ₂/λ₁ cannot be resolved there;
- the final steady state without a sink is numerically singular at low T.

A backward-stable symmetric eigensolver has an absolute eigenvalue error of order u‖S‖, so a small chemical eigenvalue can carry a large relative error even when the residual looks excellent (Peter).

**Precision needed.** Double-double (u ≈ 10⁻³²) lowers the floor to about 10⁻¹⁸ s⁻¹ for this deck, which resolves λ₁ at 300 K. Deeper wells or lower T need more, hence the arbitrary-precision reference.

## 2. Candidates checked (crates fetched and inspected, not built into MarXus)

| candidate | precision | installation | verdict |
|---|---|---|---|
| **qd 0.8** (`qd::Quad`) | double-double, 31 digits (`DIGITS = 31`), f64 exponent range | pure Rust (deps bytemuck, libm); cargo only | **chosen**: production precision |
| **faer 0.24** | f64 and `fx128 = qd::Quad` (faer-traits 0.24 implements `RealField` for it) | pure Rust; cargo only | **chosen**: dense symmetric eigendecomposition (`self_adjoint_eigen`) in f64 and double-double |
| **dashu-float 0.4** (`FBig`) | arbitrary (e.g. 192, 256 bits), arbitrary exponent range | pure Rust (dashu-int, dashu-base); no build script, no native library | **chosen**: reference path |
| astro-float 0.9 | arbitrary | pure Rust | alternative to dashu-float; every operation needs explicit precision, rounding mode and constants cache, which makes generic code clumsier |
| rug 1.30 (GMP/MPFR) | arbitrary, correctly rounded | without system headers (absent here: no gmp.h, mpfr.h) it compiles bundled GMP/MPFR C sources and needs a C compiler, make and m4 on the user's machine | **rejected** by the installation rule |
| MPLAPACK (MPFR) via C/C++ wrapper | arbitrary, full LAPACK | large C++ build | **rejected** by the installation rule |
| Rust built-in f128 | binary128 | experimental, incomplete platform support | rejected (Peter) |
| f128 crate (libquadmath) | binary128 | needs GCC's libquadmath; maintenance mode | rejected |

**Naming trap.** The Rust `qd::Quad` is *double-double* (about 31 digits). The C++ QD library also offers quad-double (about 62 digits). The Rust crate keeps the f64 exponent range, so it improves precision but not underflow (Peter). MarXus forms Boltzmann factors in log form, and √f stays ≥ 10⁻⁵⁷ on the C₂H₃ grids, so the f64 range suffices there. dashu-float has no exponent limit.

## 3. Decision

1. **Production precision: double-double (`qd::Quad`)** for the eigenvalue analysis, with f64 still selectable for speed.
2. **Symmetric solvers throughout** (detailed balance makes S = D⁻¹JD symmetric):
   - **λ₁, k_uni, λ₂:** MarXus's own banded Cholesky and shifted inverse iteration, made generic over the scalar type. They cost O(n·bw²), a few seconds in double-double for 2634 grains.
   - **Full spectrum** (N(t), λ_F): faer's dense `self_adjoint_eigen` in f64 or double-double.
3. **Reference path: dashu-float at 192 and 256 bits**, through the same generic banded inverse iteration. It is used on the difficult cases to check convergence with precision:
   - C₂H₃ at 300 and 500 K;
   - the deep test well at 125–200 K.
4. **Assembly in the selected precision.** Casting a rounded f64 matrix to higher precision cannot recover what was lost (Peter). The following are therefore computed in the selected type:
   - the collision-kernel normalization;
   - the upward transitions from detailed balance;
   - the diagonal of J as the sum of the off-diagonal loss rates plus k(E), so the column sums are exactly the losses in that precision;
   - the Boltzmann weights from logarithms;
   - the symmetrization factors.

   ρ(E) and k(E) remain f64 inputs: their relative rounding perturbs the model by ε, not the conservation.
5. **OpenBLAS stays a default feature.** Peter (20:47): "openblas is fine to have as a system library". The f64 dense symmetric eigensolver remains LAPACK DSYEVD (`--eigen-solver lapack`). faer provides the double-double dense eigensolver, which LAPACK does not have.
6. **Dependencies:**
   - qd, faer and dashu-float, all fetched and built by cargo, with versions pinned in `Cargo.lock`;
   - optionally `cargo vendor` into `MarXus/third_party/` with a `.cargo/config.toml` source replacement, so that MarXus builds fully offline from sources inside its own tree.
7. **Builds and runs on at most 4 cores** (Peter's rule).
8. **System dependencies are documented** in the root `README.md`, section "Installation and dependencies" (Peter, 20:47). It covers what each library is used for, the install commands per OS (tested ones marked), how to check it, and how to build without it. Currently that is OpenBLAS only; the planned crates are pure Rust and need no entry beyond the list of crates.

## 4. Work plan (test-driven, in this order)

1. `src/numeric/real_scalar.rs`: a `RealScalar` trait (field operations, sqrt, exp, ln, abs, ε, conversion from and to f64), implemented for f64, `qd::Quad` and a fixed-precision dashu `FBig`. Tests: the identities, and that ε and the precision of each type are what they claim.
2. Generic banded Cholesky and shifted inverse iteration (`banded_solvers.rs`, `symmetric_eigen.rs`). Tests: the f64 results unchanged; double-double resolves a 10⁻²⁰ eigenvalue that f64 cannot.
3. Generic assembly of the symmetrized operator for the eigenvalue analysis (kernel, diagonal loss sums, symmetrization in the selected type). Tests: detailed balance and column sums exact to the precision of the type.
4. `Precision::{Double, DoubleDouble, Arbitrary { bits }}` in `thermal_rate_coefficients` and the example (`--precision`).
5. faer for the dense symmetric eigendecomposition in double-double (f64 stays LAPACK/OpenBLAS, default).
6. Validation directory `validation/precision_convergence/`:
   - λ₁, λ₂, k_uni and the sum-rule deviation for f64 → double-double → 192 → 256 bits on the difficult cases;
   - plots and README.
7. Later, possibly: the steady-state solvers (final steady state at low T) in the selected precision.
