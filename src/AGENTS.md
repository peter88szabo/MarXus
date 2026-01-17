# Repository Guidelines

## Project Structure & Module Organization
This is a Rust crate organized under `src/`.
- Entry points: `src/main.rs` (binary) and `src/lib.rs` (library).
- Core domain modules: `src/rrkm/`, `src/thermal/`, `src/inertia/`, `src/tunneling/`, `src/numeric/`, `src/utils/`, plus `src/molecule.rs`.
- `src/barrierless/sacm/` now has an initial SACM entry point (`sacm_reactant_j_resolved`) for J-resolved states; deeper SACM routines are still being ported.
- Data and reference artifacts live alongside code (e.g., `src/barrierless/sacm/*.inp`, `*.txt`, `*_f90`).
- Example utilities and parsers are in `src/utils/` (e.g., `xyzparser.rs`, `print_rates.rs`).

## Build, Test, and Development Commands
Run these from the repository root (one level above `src/`):
- `cargo build`: compile the library and binary.
- `cargo run -- <args>`: run the main binary with input arguments.
- `cargo test`: run Rust tests (currently minimal).
- `cargo fmt`: format Rust code with rustfmt.
- `cargo clippy`: optional linting for common Rust issues.

## Coding Style & Naming Conventions
- Follow standard Rust style: 4-space indentation, braces on the same line, and trailing commas in multi-line lists.
- Use `snake_case` for functions/modules, `CamelCase` for types, and `SCREAMING_SNAKE_CASE` for constants.
- Keep modules focused on a single physics/chemistry concern; new functionality should live in an existing domain module when possible.

## Testing Guidelines
- Tests are sparse; unit tests live in `src/thermal/thermofuncs.rs`, `src/thermal/old_thermofuncs.rs`, and `src/utils/xyzparser.rs`, plus a parsing check in `src/utils/test_read_parse.rs`.
- Prefer unit tests in the same module file under `#[cfg(test)]`, or add targeted test helpers in `src/utils/`.
- Use descriptive test names such as `parses_xyz_with_units` and run with `cargo test`.

## RRKM Module Overview
- `src/rrkm/sum_and_density.rs` computes vibrational/rotational density and sum of states via Beyer–Swinehart counting and classical rotor formulas.
- `src/rrkm/rrkm_rate.rs` implements the RRKM expression `k(E) = (sigma_ts/sigma_cpx) * W_ts(E - dH0) / rho_cpx(E) / h`, with energy binning controlled by `dE` and `nebin`.
- Callers typically prepare frequency and rotational-constant arrays (in cm⁻¹), compute ZPE/dH0, then pass them into `get_kE`.

## Molecule Data Model
- `src/molecule.rs` defines `MoleculeStruct` with core properties (frequencies, rotational constants, ZPE, energies, symmetry, thermo/tunneling placeholders) and a `MoleculeBuilder` for validated construction.
- `MoleculeBuilder::build()` requires `freq`, either `brot` or `qxyz`, and either `ene` or `dh0`; it computes ZPE and the missing energy term.
- The builder currently uses provided `brot` directly; geometry-based `brot` computation is stubbed (see the placeholder in `build()`).

## Commit & Pull Request Guidelines
- Git history shows short, direct messages (e.g., “xyz parser”, “thermo fixde”); keep commits concise and descriptive.
- PRs should explain the scientific or numerical intent, note any changed inputs (e.g., `*.inp` files), and include sample commands or outputs if behavior changes.

## Configuration & Data Notes
- Input and reference files in `src/barrierless/sacm/` are treated as part of the codebase; document any changes to their format or assumptions.
- If adding new datasets or example inputs, keep them close to the relevant module and describe expected units and coordinate conventions.
