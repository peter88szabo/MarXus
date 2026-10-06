#!/usr/bin/env bash
# Runs MarXus on the three reference decks (input/*.inp, the decks run by both codes) with the reactant P1 (H + C2H2)
# and writes the results into marxus_output/. Run from anywhere; the repository is two levels up.
set -euo pipefail
# At most 4 cores: compilation and runs (OpenBLAS threads included).
export CARGO_BUILD_JOBS=4 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 RAYON_NUM_THREADS=4
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
cd "$repo"
cargo build -j 4 --release --example chemical_activation_from_deck
# The four methods of MarXus (three families), one run each: steady-state-olzmann (final steady state, thermal
# eigenpair), steady-state-absorbing-barrier (intermediate steady state, 10 k_BT), cse, time-integration (pulse).
o="$here/marxus_output"
run() {  # deck output-stem method [extra options]
    local deck="$1" stem="$2" method="$3"
    shift 3
    echo "running $stem"
    ./target/release/examples/chemical_activation_from_deck "$here/input/$deck.inp" P1 --threads 4 --method "$method" "$@" \
        --csv "$o/$stem.csv" 2>/dev/null > "$o/$stem.out"
}
for deck in c2h3_tight c2h3_tight_short c2h3_tight_short_notunneling; do
    for method in steady-state-olzmann steady-state-absorbing-barrier cse time-integration; do
        name=${method#steady-state-}
        run "$deck" "${deck}_${name//-/_}" "$method"
    done
done
# Sensitivity to the absorbing-barrier distance (absorbing-barrier steady state, full deck): 5 and 3 k_BT below
# the lowest threshold instead of the default 10 k_BT.
for kt in 5 3; do
    run c2h3_tight "c2h3_tight_absorbing_barrier_${kt}kT" steady-state-absorbing-barrier --barrier-kt "$kt"
done
