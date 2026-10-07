#!/usr/bin/env bash
# MarXus runs of input/c2h4_ho2.inp with the four methods, on at most 4 cores. Reactant: P1 (C2H4 + HO2).
# The deck is the MESS deck unchanged: hindered rotors (`Rotor Hindered`), Eckart tunneling, a Dummy product.
#   default                exact Eckart transmission (MarXus default)
#   mess_eckart            the MESS Eckart model (--tunneling mess-eckart), to isolate the rotor treatment
# Outputs: marxus_output/<variant>_<method>.{out,csv,_tables.csv}
set -euo pipefail
export CARGO_BUILD_JOBS=4 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 RAYON_NUM_THREADS=4
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
cd "$repo"
cargo build -j 4 --release --example chemical_activation_from_deck
exe=./target/release/examples/chemical_activation_from_deck
o="$here/marxus_output"
mkdir -p "$o"
for variant in default mess_eckart; do
    extra=()
    [[ $variant == mess_eckart ]] && extra=(--tunneling mess-eckart)
    for method in steady-state-olzmann steady-state-absorbing-barrier cse time-integration; do
        name=${method#steady-state-}
        name=${name//-/_}
        echo "running ${variant}_$name"
        $exe "$here/input/c2h4_ho2.inp" P1 --ncore 4 --method "$method" "${extra[@]}" --csv "$o/${variant}_$name.csv" \
            2> "$o/${variant}_$name.err" > "$o/${variant}_$name.out"
    done
done
echo "done: $o"
