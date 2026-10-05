#!/usr/bin/env bash
# Olzmann's eigenvalue analysis of J without absorbing barrier for the three reference decks (input/*.inp,
# the decks run by both codes), reactant P1 (H + C2H2). Results in marxus_output/.
#   inverse iteration with the banded Cholesky factor of S + sigma I (default solver): all decks;
#   full decomposition by LAPACK DSYEVD: all decks;
#   full decomposition by the in-house Householder/QL (Olzmann's tred2/tql2 route): the two 1000 K decks
#   only (O(n^3) without blocking; about 2 min per condition for the 2634 grains of the full deck).
# k_uni is the eigenvector average (GO10 after eq. 12); lambda_1, the precision floor and the sum-rule
# deviation are printed beside it, with warnings (no rejection) above the tolerance (default 1.5e-2).
# At most 4 cores: compilation and runs (OpenBLAS threads included).
set -euo pipefail
export CARGO_BUILD_JOBS=4 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 RAYON_NUM_THREADS=4
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
cd "$repo"
cargo build -j 4 --release --example chemical_activation_from_deck
run() {  # deck solver
    echo "running $1 ($2)"
    ./target/release/examples/chemical_activation_from_deck "$here/input/$1.inp" P1 \
        --steady-state eigenvalue --eigen-solver "$2" \
        | grep -v "Iterative diagonalization" > "$here/marxus_output/$1_$2.out"
}
for deck in c2h3_tight c2h3_tight_short c2h3_tight_short_notunneling; do
    run "$deck" inverse
    run "$deck" lapack
done
for deck in c2h3_tight_short c2h3_tight_short_notunneling; do
    run "$deck" full
done
