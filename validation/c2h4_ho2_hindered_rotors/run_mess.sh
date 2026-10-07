#!/usr/bin/env bash
# MESS reference run of input/c2h4_ho2.inp (static MESS binary of the 2026 source), on at most 4 cores.
# Outputs: reference_mess/c2h4_ho2.{out,log}
set -euo pipefail
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mess=/home/peter/Dropbox/Research_Leuven/MESS_kinetics/Source_from_2026/MESS/static/mess
cd "$here/reference_mess"
cp "$here/input/c2h4_ho2.inp" c2h4_ho2.inp
"$mess" c2h4_ho2.inp > mess_stdout.txt 2>&1
echo "done: $here/reference_mess"
