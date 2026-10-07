#!/usr/bin/env bash
# Variant with excited electronic test levels (input/c2h4_ho2_electronic_levels.inp): MESS and MarXus (CSE and
# SteadyStateOlzmann, MESS Eckart model), at most 4 cores.
# Outputs: reference_mess_electronic_levels/, marxus_output/electronic_levels_{cse,olzmann}.*
set -euo pipefail
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 RAYON_NUM_THREADS=4 CARGO_BUILD_JOBS=4
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
mess=/home/peter/Dropbox/Research_Leuven/MESS_kinetics/Source_from_2026/MESS/static/mess
(cd "$here/reference_mess_electronic_levels" && cp "$here/input/c2h4_ho2_electronic_levels.inp" . && "$mess" c2h4_ho2_electronic_levels.inp > mess_stdout.txt 2>&1)
cd "$repo"
cargo build -j 4 --release --example chemical_activation_from_deck
exe=./target/release/examples/chemical_activation_from_deck
o="$here/marxus_output"
for method in cse steady-state-olzmann; do
    name=${method#steady-state-}
    $exe "$here/input/c2h4_ho2_electronic_levels.inp" P1 --ncore 4 --method "$method" --tunneling mess-eckart \
        --csv "$o/electronic_levels_$name.csv" 2> "$o/electronic_levels_$name.err" > "$o/electronic_levels_$name.out"
done
echo "done"
