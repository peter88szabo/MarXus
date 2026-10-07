#!/usr/bin/env bash
# Shock heating of C2H3 (input/c2h3_shock_incubation.inp: 300 K population, constant 1000 K, 1 atm bath) at four grain
# widths, EnergyStepOverTemperature 0.4, 0.2, 0.1, 0.05 (grains of about 278, 139, 70 and 35 cm-1). Each run reports the
# incubation time, the vibrational relaxation time, the flux coefficients and the CSE description propagated in time;
# shock_incubation_summary.py collects them into shock_incubation_grain_convergence.csv. Run from anywhere.
set -euo pipefail
# At most 8 cores: compilation and runs.
export CARGO_BUILD_JOBS=8 OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 RAYON_NUM_THREADS=8
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
cd "$repo"
cargo build -j 8 --release --example chemical_activation_from_deck
o="$here/marxus_output"
decks="$(mktemp -d)"
trap 'rm -rf "$decks"' EXIT
for step in 0.4 0.2 0.1 0.05; do
    stem="c2h3_shock_incubation_de${step}"
    sed "s/^EnergyStepOverTemperature .*/EnergyStepOverTemperature           $step/" \
        "$here/input/c2h3_shock_incubation.inp" > "$decks/$stem.inp"
    echo "running $stem"
    ./target/release/examples/chemical_activation_from_deck "$decks/$stem.inp" --ncore 8 \
        --csv "$o/$stem.csv" 2>/dev/null > "$o/$stem.out"
done
python3 "$here/shock_incubation_summary.py" "$o" > "$here/shock_incubation_grain_convergence.csv"
cat "$here/shock_incubation_grain_convergence.csv"
