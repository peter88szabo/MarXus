#!/usr/bin/env bash
# MarXus runs of the acetyl + O2 decks (input/, written by make_deck.py from the MESMER XML inputs), every deck with
# the four methods, on at most 4 cores. Reactant: R (acetyl + O2).
#   input/acetyl_o2.inp              298 K, 200.72 Torr He; Eckart at TS1 with classical barriers (MESMER example)
#   input/acetyl_o2_zpe_eckart.inp   the same with zero-point Eckart barriers (MESMER Tunnelling Ex1/Ex2)
#   input/acetyl_o2_250K.inp         250 K, 37.48 Torr; no tunneling at TS1 (MESMER reservoirSink)
# Outputs: marxus_output/<deck>_<method>.{out,csv,_tables.csv}
set -euo pipefail
export CARGO_BUILD_JOBS=4 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 RAYON_NUM_THREADS=4
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
cd "$repo"
cargo build -j 4 --release --example chemical_activation_from_deck
exe=./target/release/examples/chemical_activation_from_deck
o="$here/marxus_output"
mkdir -p "$o"
for deck in acetyl_o2 acetyl_o2_zpe_eckart acetyl_o2_250K; do
    for method in steady-state-olzmann steady-state-absorbing-barrier cse time-integration; do
        name=${method#steady-state-}
        name=${name//-/_}
        echo "running ${deck}_$name"
        $exe "$here/input/$deck.inp" R --ncore 4 --method "$method" --csv "$o/${deck}_$name.csv" \
            2> "$o/${deck}_$name.err" > "$o/${deck}_$name.out"
    done
done
echo "done: $o"
