#!/usr/bin/env bash
# MarXus runs for the ZZ-allyl + O2 Gamma Case 2 network (reference: the MESS run of 2025-09-29 in this
# directory, *.inp/*.log/*.out, untouched). At most 4 cores (compilation, rayon and OpenBLAS threads).
#   marxus_input/case2_tstlevel_E.inp: the MESS deck with "TSTLevel E" added to the two phase-space-theory
#     cores (B12, B6P7), because the reference MESS run used the E level (its log: "TST level: E"); the
#     MarXus default, as in the current MESS source, is EJ.
# The four methods of MarXus (three families), one run each:
#   steady-state-olzmann            final steady state (+ thermal eigenpair, thermal fates of the wells)
#   steady-state-absorbing-barrier  intermediate steady state (absorbing barrier 10 k_BT)
#   cse                             phenomenological rate coefficients (chemically significant eigenvalues)
#   time-integration                direct time integration of a pulse (Rodas4, 1e-12 .. 1e2 s)
# for the exact Eckart tunneling (MarXus default) and the MESS Eckart model (--tunneling mess-eckart, to compare
# with MESS). Outputs (marxus_output/): every run writes its human-readable report (*.out) and its
# machine-readable tables (*.csv and *_tables.csv, --csv), which compare_with_mess.py reads:
#   case2_tstlevel_E_<method>.*               exact Eckart
#   case2_tstlevel_E_mess_eckart_<method>.*   MESS Eckart model
#   case2_default_EJ_absorbing_barrier.*      the original deck as is (PST at the EJ level): capture comparison
#   eckart_kappa.csv                          canonical Eckart factors kappa(T) of the tunneling barriers
#                                             (examples/eckart_kappa_from_deck.rs), 100-2000 K: exact and MESS model
set -euo pipefail
export CARGO_BUILD_JOBS=4 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 RAYON_NUM_THREADS=4
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
cd "$repo"
cargo build -j 4 --release --example chemical_activation_from_deck --example eckart_kappa_from_deck
exe=./target/release/examples/chemical_activation_from_deck
deck_e="$here/marxus_input/case2_tstlevel_E.inp"
deck_orig="$here/Gamma-Case2_from_ZZ_allyl+O2_Gamma4-Escape_CCSDT_version-1_2025_09_29_12.7kcal.inp"
o="$here/marxus_output"
mkdir -p "$o"
run() {  # deck output-stem method [extra options]
    local deck="$1" stem="$2" method="$3"
    shift 3
    echo "running $stem"
    $exe "$deck" R --ncore 4 --method "$method" "$@" --csv "$o/$stem.csv" 2>/dev/null > "$o/$stem.out"
}
for method in steady-state-olzmann steady-state-absorbing-barrier cse time-integration; do
    name=${method#steady-state-}
    name=${name//-/_}
    run "$deck_e" "case2_tstlevel_E_$name" "$method"
    run "$deck_e" "case2_tstlevel_E_mess_eckart_$name" "$method" --tunneling mess-eckart
done
run "$deck_orig" case2_default_EJ_absorbing_barrier steady-state-absorbing-barrier
./target/release/examples/eckart_kappa_from_deck "$deck_e" 2>/dev/null > "$o/eckart_kappa.csv"
echo "done: $o"
