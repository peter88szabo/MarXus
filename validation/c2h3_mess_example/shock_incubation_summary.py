"""Grain-width convergence of the shock-heating observables of C2H3 (run_shock_incubation.sh).

For every run marxus_output/c2h3_shock_incubation_de<step>.csv: the grain width, the incubation time t_inc (its late-time
plateau at k_uni t = 1, Barker, King, J. Chem. Phys. 103, 4953 (1995), eq. 9) in s and in collisions, the vibrational relaxation time
tau_vib at the output time nearest t_inc (as Barker and King compare it), the final steady-state energy E_f, the
late-time flux coefficient r(W1->P1) (= k_uni), and the agreement time t* of the CSE description (relative deviation
within the tolerance) with t* lambda_relax. Usage: shock_incubation_summary.py <marxus_output directory>.
"""
import csv
import glob
import os
import re
import sys


def blocks(path):
    """The titled CSV blocks of a MarXus machine-readable file: {title: (comment lines, header, rows)}. A block starts at
    a line "# prepared time integration..." or "# CSE description..."; further "# " lines before its header are its
    comments; the lines before the first block are skipped."""
    out, title = {}, None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(("# prepared time integration", "# CSE description")):
                title = line[2:]
                out[title] = ([], None, [])
            elif title is None or not line.strip():
                continue
            elif line.startswith("# "):
                out[title][0].append(line[2:])
            elif out[title][1] is None:
                out[title] = (out[title][0], next(csv.reader([line])), [])
            else:
                out[title][2].append([float(x) for x in next(csv.reader([line]))])
    return out


def main(directory):
    print("EnergyStepOverTemperature,grain[cm-1],t_inc[s],t_inc_collisions,tau_vib_at_t_inc[s],E_f[cm-1],"
          "r_late(W1->P1)[1/s],t_star[s],t_star_lambda_relax")
    runs = []
    for path in glob.glob(os.path.join(directory, "c2h3_shock_incubation_de*.csv")):
        match = re.search(r"_de([0-9.]+)\.csv$", path)
        if match is None:  # the _tables.csv files
            continue
        step = float(match.group(1))
        with open(path) as f:
            grain = float(re.search(r"grain ([0-9.]+) cm-1", f.readline()).group(1))
        b = blocks(path)
        _, header, rows = next(v for k, v in b.items() if k.startswith("prepared time integration"))
        col = {name: i for i, name in enumerate(header)}
        # t_inc(t) = t + ln(N/N_ref)/k_inst is a small difference of large numbers at t >> t_inc (its relative error is
        # about (t/t_inc) times that of k_inst): it is read on its plateau at k_uni t = 1.
        k_uni = rows[-1][col["r(W1->P1)"]]
        at = min(rows, key=lambda r: abs(r[col["t[s]"]] * k_uni - 1.0))
        t_inc = at[col["t_inc"]]
        z_rate = rows[-1][col["Z(W1)"]] / rows[-1][col["t[s]"]]
        near = min(rows, key=lambda r: abs(r[col["t[s]"]] - t_inc))
        cse = next((v for k, v in b.items() if k.startswith("CSE description")), None)
        t_star, star_lambda = "", ""
        if cse is not None:
            m = re.search(r"t\* \[s\] ([0-9.e+-]+); t\* lambda_relax ([0-9.e+-]+)", " ".join(cse[0]))
            if m:
                t_star, star_lambda = m.group(1), m.group(2)
        runs.append((step, grain, t_inc, t_inc * z_rate, near[col["tau_vib(W1)"]], near[col["E_f(W1)"]],
                     rows[-1][col["r(W1->P1)"]], t_star, star_lambda))
    for r in sorted(runs, reverse=True):
        print(",".join(f"{x:.6e}" if isinstance(x, float) else str(x) for x in r))


if __name__ == "__main__":
    main(sys.argv[1])
