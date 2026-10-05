#!/usr/bin/env bash
# Runs MarXus on the three reference decks (input/*.inp, the decks run by both codes) with the reactant P1 (H + C2H2)
# and writes the results into marxus_output/. Run from anywhere; the repository is two levels up.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
cd "$repo"
cargo build --release --example chemical_activation_from_deck
for deck in c2h3_tight c2h3_tight_short c2h3_tight_short_notunneling; do
    echo "running $deck"
    ./target/release/examples/chemical_activation_from_deck "$here/input/$deck.inp" P1 \
        | grep -v "Iterative diagonalization" > "$here/marxus_output/$deck.out"
done
# Sensitivity to the absorbing-barrier distance (intermediate steady state, full deck): 5 and 3 k_BT below
# the lowest threshold instead of the default 10 k_BT.
for kt in 5 3; do
    echo "running c2h3_tight with the absorbing barrier ${kt} k_BT below the threshold"
    ./target/release/examples/chemical_activation_from_deck "$here/input/c2h3_tight.inp" P1 \
        --steady-state intermediate --barrier-kt "$kt" \
        | grep -v "Iterative diagonalization" > "$here/marxus_output/c2h3_tight_barrier_${kt}kT.out"
done
