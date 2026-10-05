#!/usr/bin/env bash
# MarXus runs for the ZZ-allyl + O2 Gamma Case 2 network (reference: the MESS run of 2025-09-29 in this
# directory, *.inp/*.log/*.out, untouched). At most 4 cores (compilation and OpenBLAS threads).
#   marxus_input/case2_tstlevel_E.inp: the MESS deck with "TSTLevel E" added to the two phase-space-theory
#     cores (B12, B6P7), because the reference MESS run used the E level (its log: "TST level: E"); the
#     MarXus default, as in the current MESS source, is EJ.
#   Outputs (marxus_output/):
#     case2_tstlevel_E_steady_states.out   steady-state method: intermediate (absorbing barrier) and final
#                                          steady state; the final one with its thermal rate coefficients from
#                                          the lowest eigenpair of J (GO10 eq. 12: k_uni, k_inf of every channel)
#     case2_default_EJ_steady_states.out   the original deck as is (PST at the EJ level), for comparison
#     case2_tstlevel_E_mess_eckart_*.out   the same with the MESS Eckart tunneling model (--tunneling mess-eckart)
#     case2_tstlevel_E_mess_eckart_cse.out phenomenological rate coefficients from the chemically significant
#                                          eigenvalues (--method cse; Miller, Klippenstein 2006;
#                                          Georgievskii et al. 2013), MESS Eckart model, LAPACK
#     eckart_kappa.csv                     canonical Eckart factors kappa(T) of the tunneling barriers
#                                          (examples/eckart_kappa_from_deck.rs), 100-2000 K: exact and MESS model
set -euo pipefail
export CARGO_BUILD_JOBS=4 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 RAYON_NUM_THREADS=4
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
cd "$repo"
cargo build -j 4 --release --example chemical_activation_from_deck --example eckart_kappa_from_deck
exe=./target/release/examples/chemical_activation_from_deck
deck_e="$here/marxus_input/case2_tstlevel_E.inp"
deck_orig="$here/Gamma-Case2_from_ZZ_allyl+O2_Gamma4-Escape_CCSDT_version-1_2025_09_29_12.7kcal.inp"
mkdir -p "$here/marxus_output"
$exe "$deck_e" R --steady-state both | grep -v "Iterative diagonalization" > "$here/marxus_output/case2_tstlevel_E_steady_states.out"
$exe "$deck_e" R --steady-state both --tunneling mess-eckart | grep -v "Iterative diagonalization" > "$here/marxus_output/case2_tstlevel_E_mess_eckart_steady_states.out"
$exe "$deck_e" R --method cse --tunneling mess-eckart 2>/dev/null | grep -v "Iterative diagonalization" > "$here/marxus_output/case2_tstlevel_E_mess_eckart_cse.out"
$exe "$deck_orig" R --steady-state both | grep -v "Iterative diagonalization" > "$here/marxus_output/case2_default_EJ_steady_states.out"
./target/release/examples/eckart_kappa_from_deck "$deck_e" 2>/dev/null | grep -v "Iterative diagonalization" > "$here/marxus_output/eckart_kappa.csv"
echo "done: $here/marxus_output"
