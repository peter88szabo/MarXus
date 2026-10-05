# MarXus

**Molecular Statistical Physics for Kinetics and Thermochemistry**

**Author:** Peter Szabo  
**Email:** peter88szabo@gmail.com  

MarXus is a molecular statistical-physics toolkit for **thermochemistry** and **chemical kinetics**, written in Rust. It covers partition functions, state counting, statistical rate theories (TST, RRKM, ILT, PST, SACM) and an energy-grained **master equation** for gas-phase reactions. The master equation follows Olzmann's formulation and solution method: the steady states with physical sinks, the eigenvalue analysis, and consecutive chemical activation.

The master-equation, tunneling, ILT and numerical code cites the source of each equation (paper and equation number) in its comments. Design notes, derivations and validation results are in `reports/`; validations against reference calculations are in `validation/`.

---

## Current state (2026-10-05)

| area | status |
|---|---|
| Thermochemistry (RRHO, Grimme qRRHO) | implemented |
| Sum and density of states, RRKM k(E), canonical TST | implemented, tested |
| Tunneling (exact Eckart, microcanonical and canonical) | implemented, tested |
| Inverse Laplace transform (ILT) for barrierless channels | implemented, tested |
| Phase space theory (PST) with arbitrary 1D potential, SACM | in progress |
| Multiwell chemical-activation master equation (steady states) | implemented, tested, validated (C₂H₃) |
| Olzmann eigenvalue analysis (k_uni, λ₁, λ₂, N(t)) | implemented, tested, validated (C₂H₃, 300–2000 K) |
| Higher precision (double-double, arbitrary-precision reference) | planned (`reports/higher_precision_decision.md`) |

The source contains 144 unit tests (`cargo test`).

---

## Features

### Thermochemistry
- Molecules built from vibrational frequencies (with scaling), rotational constants, or Cartesian geometries; the moments of inertia are computed from the geometry.
- Partition functions and thermodynamic functions **U, H, F, G, S, Cv, Cp** in the rigid-rotor–harmonic-oscillator approximation.
- **Quasi-RRHO entropy** with free-rotor interpolation for low-frequency modes (Grimme, Chem. Eur. J. 18, 9955 (2012)).
- Equilibrium constants: used inside the master equation (k∞,assoc/k∞,diss for detailed balance); a general thermochemistry routine is not available yet (see To Do).

### State counting and microcanonical rate theory
- **Sum and density of states** by direct (Beyer–Swinehart) counting, with classical 1D, 2D and 3D rotors (Forst); rovibrational and bimolecular (convolved) states.
- **RRKM / microcanonical TST** specific rate coefficients k(E).
- **Canonical TST** high-pressure rate coefficients, with tunneling corrections.
- **Energy graining** on 1 cm⁻¹ cells averaged to grains (Pilling, Robertson, Annu. Rev. Phys. Chem. 54, 245 (2003)). It reproduces canonical TST within 0.08% for C₂H₃.

### Tunneling
- **Exact Eckart** transmission probability (Miller, J. Am. Chem. Soc. 101, 6810 (1979), eq. 8; Johnston, Heicklen 1962), overflow-free.
- **Microcanonical** tunneling sum of states by convolution with the transition-state sum of states (Miller 1979, eq. 9); the canonical correction κ(T) is the same in both directions.
- Wigner, Bell and Skodje–Truhlar corrections as optional models.

### Barrierless reactions
- **Inverse Laplace transform (ILT)** of modified-Arrhenius high-pressure rate coefficients, for association and dissociation, with a grain-integrated kernel (Davies, Green, Pilling, Chem. Phys. Lett. 126, 373 (1986)).
- Phase space theory (Troe–Ushakov 2006 form; capture models): **in progress**.
- Statistical adiabatic channel model (SACM): **in progress**.

### Master equation: multiwell chemical activation

The energy-grained master equation is J·N = R·F with J = ω(I − P) + K + k_c[D]·I (Pfeifle, Olzmann, Int. J. Chem. Kinet. 46, 231 (2014), eq. 2).

**Network**
- Any number of wells on a common absolute energy grid, connected by isomerization (exact detailed balance, k = W‡/(hρ) in both directions).
- Product channels, and pseudo-first-order bimolecular sinks k_c[D].

**Collisions**
- **Exponential down** (default), exactly normalized (Robertson (ed.), Comprehensive Chemical Kinetics 43 (2019), eq. 4.16).
- Olzmann **stepladder** (Olzmann, Gebhardt, Scherzer, Int. J. Chem. Kinet. 23, 825 (1991)).
- Both obey detailed balance exactly; Lennard-Jones collision frequencies (Troe 1977).

**Sources**
- Thermal entrance channels (chemical activation from a bimolecular reactant).
- A given distribution.
- **Consecutive chemical activation:** coupled master equations, the output of one feeding the next (PO14 pp. 236–237).

**Steady states (Olzmann)**
- **Final steady state** with physical sinks and no absorbing barrier.
- **Intermediate steady state** with an absorbing barrier, at a user-chosen distance below the threshold (default 10 kT).
- Solvers: banded Cholesky of the symmetrized operator, or BiCGSTAB.
- Every result is checked for its residual and mass balance.

**Observables**
- Yields of products, stabilization and sink (Olzmann, PCCP 4, 3614 (2002), eq. 10).
- Chemical-activation rate coefficients k^ca (González-García, Olzmann, PCCP 12, 12290 (2010), eq. 9).
- Total loss rate coefficients.
- Bimolecular rate coefficients k(R → X) = k∞Φ_X (PR03 eq. 44).

**Eigenvalue analysis (Olzmann's solution method, no absorbing barrier)**
- **k_uni**: the average of k(E) over the thermal eigenvector (GO10, after eq. 12). It is reported beside λ₁, the lowest eigenvalue (GO10 eq. 12), with the sum-rule check λ₁ = k_uni. A warning is given above 1.5%, and the double-precision floor is printed.
- Channel rate coefficients, high-pressure limits, λ₂/k_uni (separation of time scales).
- Association rate coefficients by detailed balance.
- Three eigen-solvers:
  - shifted inverse iteration with the banded Cholesky factor (default, most accurate);
  - Householder/QL (EISPACK tred2/tql2, Olzmann's route);
  - LAPACK DSYEVD.
- **Time-dependent populations** N(t) by eigenvalue expansion (PO14 eqs. 3–4), and the validity windows of the steady-state picture (O02).

**Input**
- Decks in the MESS input format (a subset), so existing decks can be used.
- A MarXus keyword block adds ILT parameters to barriers.
- Exact Eckart tunneling and 1 cm⁻¹ state counting are built from the deck.

**Validation**
- `validation/c2h3_mess_example/` (absorbing barrier) and `validation/c2h3_mess_example_olzmann_eigen/` (eigenvalue analysis) compare H + C₂H₂ ⇌ C₂H₃ at 300–2000 K and 0.1–10 atm with the stored MESS results.
- The eigenvalue route agrees within ±2.5% from 750 to 1750 K. At 300–1000 K the residual offsets are the known exact vs semiclassical Eckart difference.
- Where both routes are valid, they agree with each other to 0.01%.
- All validity checks of the master equation (sum rules, detailed balance, limits, solver agreement) are described in `reports/master_equation_validity_checks.md`.

### Numerical library (`src/numeric/`)
- Banded Cholesky (factor once, solve many) and LDLᵀ with Bunch–Kaufman pivoting.
- BiCGSTAB and GMRES.
- Symmetric eigensolvers: Householder + implicit QL, Jacobi, inverse iteration with deflation and shift.
- LAPACK interface (DSYEVD, system OpenBLAS); tridiagonal solvers; Lanczos Γ function.

---

## Usage example: master equation from a deck

```
cargo run --release --example chemical_activation_from_deck -- deck.inp REACTANT \
    --steady-state all --eigen-solver inverse
```

Options:
- `--steady-state intermediate|final|eigenvalue|both|all`;
- `--barrier-kt X` (absorbing-barrier distance);
- `--eigen-solver inverse|full|lapack`;
- `--sum-rule-tolerance X`.

The output tables begin with comment lines that explain the quantities, with their references. The source file `examples/chemical_activation_from_deck.rs` documents all columns.

---

## Installation and dependencies

MarXus needs **one system library (OpenBLAS)** besides the Rust toolchain. All other dependencies are Rust crates that `cargo` downloads and compiles by itself; nothing else has to be installed or maintained by hand.

### 1. Rust toolchain (required)

Install Rust with rustup (https://rustup.rs):

```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Tested with rustc/cargo 1.98.

### 2. OpenBLAS: BLAS and LAPACK (system library, needed by the default build)

**What uses it.** The default cargo feature `openblas` links the system OpenBLAS. MarXus uses it for the LAPACK symmetric eigensolver DSYEVD, `--eigen-solver lapack` in the master-equation eigenvalue analysis.

**Install it:**

| system | command | status |
|---|---|---|
| Ubuntu / Debian | `sudo apt install libopenblas-dev` | tested (Ubuntu 24.04, libopenblas-dev 0.3.26) |
| Fedora / RHEL | `sudo dnf install openblas-devel` | not tested |
| Arch Linux | `sudo pacman -S openblas` | not tested |
| macOS (Homebrew) | `brew install openblas`, then before building `export RUSTFLAGS="-L $(brew --prefix openblas)/lib"` (Homebrew does not put OpenBLAS on the default library path) | not tested |

**Check that it is found:**

```
ldconfig -p | grep libopenblas        # Linux: should list libopenblas.so
```

**Without OpenBLAS.** If the library is not available, build without it:

```
cargo build --release --no-default-features
```

Everything works except `--eigen-solver lapack`, which then stops with an error naming the missing feature. The default solver (inverse iteration) and the in-house full decomposition (`--eigen-solver full`) need no external library.

**Threads.** OpenBLAS uses all cores by default. To limit it, e.g. to 4 cores:

```
export OPENBLAS_NUM_THREADS=4
```

### 3. Rust crates (automatic)

The current version has no crate dependencies.

The planned higher-precision eigenvalue analysis (`reports/higher_precision_decision.md`) will use pure-Rust crates only:
- `qd`: double-double;
- `faer`: linear algebra;
- `dashu-float`: arbitrary precision.

`cargo build` will download and compile them by itself. They need no system library.

### 4. Python (optional, only for the validation plots)

The plotting scripts in `validation/*/plot_comparison.py` need Python 3 with `numpy` and `matplotlib`:

```
python3 -m venv ~/.venvs/science && source ~/.venvs/science/bin/activate && pip install numpy matplotlib
```

### Build and test

```
cargo build --release -j 4
cargo test -j 4 -- --test-threads=4
```

`-j 4` and `--test-threads=4` limit compilation and tests to 4 cores.

---

---

## To Do (not implemented yet)
- Higher precision for the master equation: double-double assembly and solvers with an arbitrary-precision reference path (planned, `reports/higher_precision_decision.md`).
- Phenomenological rate coefficients of multiwell networks from the chemically significant eigenvalues.
- Treatment of the stepladder model in the eigenvalue analysis when its step spans several grains (independent sub-equations).
- Excited electronic states in the partition functions and state counts.
- A general equilibrium-constant routine (thermochemistry).
- Microcanonical Variational TST (μVTST).
- Canonical Variational TST (CVTST).
- Submerged barrier with a pre-reaction vdW complex: **μ-canonical, J-resolved 2-TST treatment**.

---

## License
This project is licensed under the **GNU General Public License v3.0 (GPL-3.0)**.
