#!/usr/bin/env bash
# MarXus runs for the ZZ-allyl + O2 Gamma Case 2 network (reference: the MESS run of 2025-09-29 in this
# directory, *.inp/*.log/*.out, untouched). At most 4 cores (compilation and OpenBLAS threads).
#   marxus_input/case2_tstlevel_E.inp: the MESS deck with "TSTLevel E" added to the two phase-space-theory
#     cores (B12, B6P7), because the reference MESS run used the E level (its log: "TST level: E"); the
#     MarXus default, as in the current MESS source, is EJ.
#   Outputs (marxus_output/):
#     case2_tstlevel_E_steady_states.out   intermediate (absorbing barrier) and final steady state
#     case2_tstlevel_E_eigenvalue.out      eigenvalue analysis (k_uni, k_inf of every channel)
#     case2_default_EJ_steady_states.out   the original deck as is (PST at the EJ level), for comparison
set -euo pipefail
export CARGO_BUILD_JOBS=4 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 RAYON_NUM_THREADS=4
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
cd "$repo"
cargo build -j 4 --release --example chemical_activation_from_deck
exe=./target/release/examples/chemical_activation_from_deck
deck_e="$here/marxus_input/case2_tstlevel_E.inp"
deck_orig="$here/Gamma-Case2_from_ZZ_allyl+O2_Gamma4-Escape_CCSDT_version-1_2025_09_29_12.7kcal.inp"
mkdir -p "$here/marxus_output"
$exe "$deck_e" R --steady-state both | grep -v "Iterative diagonalization" > "$here/marxus_output/case2_tstlevel_E_steady_states.out"
$exe "$deck_e" R --steady-state eigenvalue 2>/dev/null | grep -v "Iterative diagonalization" > "$here/marxus_output/case2_tstlevel_E_eigenvalue.out"
$exe "$deck_orig" R --steady-state both | grep -v "Iterative diagonalization" > "$here/marxus_output/case2_default_EJ_steady_states.out"
echo "done: $here/marxus_output"
