#!/usr/bin/env python3
"""C2H4 + HO2 hindered-rotor deck: equilibrium constants, species partition functions and kappa, MarXus vs MESS.

Reads
  reference_mess/c2h4_ho2.log, reference_mess_electronic_levels/c2h4_ho2_electronic_levels.log
      MESS logs (LogPrecision 6): "Real equilibrium constants", "partition functions (relative to the ground level,
      units - 1/cm3)", "isomers-to-bimolecular equilibrium coefficients (kappa matrix)"
  marxus_output/mess_eckart_cse.{out,csv}, marxus_output/electronic_levels_cse.out
      MarXus (run_marxus.sh, run_electronic_levels.sh): section PARTITION FUNCTIONS AND EQUILIBRIUM CONSTANTS of the
      report, kappa blocks of the CSE tables
and writes
  equilibrium_constants_comparison.csv   K(W2/P1) = [W2]/([C2H4][HO2]) (cm3), base deck and electronic-level variant,
                                         and the change of K by the electronic levels in each code
  partition_functions_comparison.csv     internal partition functions of W2 and of P1 (with the relative translation,
                                         per cm3): MESS divided by the centre-of-mass translation of each species
  isomer_bimolecular_kappa_comparison.csv   kappa(W2, P1) and kappa(W2, P2) (G13 eq. 34) at every condition

MESS's species partition functions include the centre-of-mass translation of each species per cm3 (P1's fragments
are given separately, P1_0 = HO2, P1_1 = C2H4); it is divided out with the masses of the most abundant isotopes.

Run with the science environment:  source ~/.venvs/science/bin/activate && python3 compare_equilibrium_and_kappa.py
"""
import csv
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
MX = os.path.join(HERE, "marxus_output")
TEMPERATURES = [300.0, 500.0, 700.0, 1000.0, 1500.0, 2000.0]
PRESSURES_BAR = [0.01, 0.1, 1.0, 10.0, 100.0]
H, K_B, AMU = 6.62607015e-34, 1.380649e-23, 1.66053906660e-27
M_H, M_C, M_O = 1.00782503223, 12.0, 15.99491461957
MASS = {"W2": 2 * M_C + 5 * M_H + 2 * M_O, "HO2": M_H + 2 * M_O, "C2H4": 2 * M_C + 4 * M_H}


def translation_per_cm3(mass_amu, t):
    return (2.0 * math.pi * mass_amu * AMU * K_B * t / H**2) ** 1.5 * 1.0e-6


def write_csv(name, rows):
    with open(os.path.join(HERE, name), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def mess_log(path):
    """Real equilibrium constants K(W2/P1) by temperature (in the order of the temperatures), the species partition
    functions {T: {species: Q}} and the kappa matrices in the order of the conditions (T outer, p inner)."""
    lines = open(path).read().split("\n")
    equilibrium, kappas, partition = [], [], {}
    for i, line in enumerate(lines):
        if "Real equilibrium constants" in line:
            equilibrium.append(float(lines[i + 2].split()[2]))
        if "kappa matrix" in line:
            header = lines[i + 1].split()[1:]
            fields = lines[i + 2].split()
            kappas.append(dict(zip(header, map(float, fields[1:]))))
        if line.startswith("partition functions (relative to the ground level"):
            header = lines[i + 1].split()[1:]
            k = i + 2
            while lines[k].split() and lines[k].split()[0].isdigit():
                fields = lines[k].split()
                partition[float(fields[0])] = dict(zip(header, map(float, fields[1:])))
                k += 1
    conditions = [(t, p) for t in TEMPERATURES for p in PRESSURES_BAR]
    return dict(zip(TEMPERATURES, equilibrium)), dict(zip(conditions, kappas)), partition


def marxus_report(path):
    """{T: {"Q": {species: Q}, "K": K(W2/P1)}} from the section PARTITION FUNCTIONS AND EQUILIBRIUM CONSTANTS."""
    text = open(path).read().split("PARTITION FUNCTIONS AND EQUILIBRIUM CONSTANTS")[1].split("_" * 20)[0]
    out = {}
    for part in re.split(r"Temperature = ", text)[1:]:
        t = float(part.split()[0])
        q = {m.group(1): float(m.group(2)) for m in re.finditer(r"^\s+(\S+)\s+(?:well|barrier|bimolecular)\s+\S+\s+(\S+)$", part, re.M)}
        k = float(re.search(r"^W2\s+\S+\s+(\S+)$", part, re.M).group(1))
        out[t] = {"Q": q, "K": k}
    return out


def marxus_kappa(path):
    """{(T, p_bar): {channel: kappa(W2, channel)}} from the CSE tables."""
    out, key, header = {}, None, None
    for line in open(path):
        line = line.rstrip("\n")
        m = re.match(r"# T = (\S+) K, p = (\S+) Torr", line)
        if m:
            key = (float(m.group(1)), round(float(m.group(2)) / 750.0616827041697, 6))
        elif line.startswith("W\\P,"):
            header = line.split(",")[1:]
        elif header and line.startswith("W2,"):
            out[key] = dict(zip(header, map(float, line.split(",")[1:])))
            header = None
    return out


k_mess, kappa_mess, q_mess = mess_log(os.path.join(HERE, "reference_mess", "c2h4_ho2.log"))
k_mess_levels, _, _ = mess_log(os.path.join(HERE, "reference_mess_electronic_levels", "c2h4_ho2_electronic_levels.log"))
mx = marxus_report(os.path.join(MX, "mess_eckart_cse.out"))
mx_levels = marxus_report(os.path.join(MX, "electronic_levels_cse.out"))

rows = []
for t in TEMPERATURES:
    rows.append({
        "T_K": t,
        "K_mess": k_mess[t], "K_marxus": mx[t]["K"], "ratio": mx[t]["K"] / k_mess[t],
        "K_mess_levels": k_mess_levels[t], "K_marxus_levels": mx_levels[t]["K"], "ratio_levels": mx_levels[t]["K"] / k_mess_levels[t],
        "change_by_levels_mess_percent": 100 * (k_mess_levels[t] / k_mess[t] - 1),
        "change_by_levels_marxus_percent": 100 * (mx_levels[t]["K"] / mx[t]["K"] - 1),
    })
write_csv("equilibrium_constants_comparison.csv", rows)

rows = []
for t in TEMPERATURES:
    q_w2 = q_mess[t]["W2"] / translation_per_cm3(MASS["W2"], t)
    q_p1 = q_mess[t]["P1_0"] * q_mess[t]["P1_1"] / translation_per_cm3(MASS["W2"], t)
    rows.append({
        "T_K": t,
        "Q_W2_mess": q_w2, "Q_W2_marxus": mx[t]["Q"]["W2"], "ratio_W2": mx[t]["Q"]["W2"] / q_w2,
        "Q_P1_mess": q_p1, "Q_P1_marxus": mx[t]["Q"]["P1"], "ratio_P1": mx[t]["Q"]["P1"] / q_p1,
        "half_cell_factor_exp(-dE/2kT)": math.exp(-0.5 / (0.69503476 * t)),
    })
write_csv("partition_functions_comparison.csv", rows)

kappa_mx = marxus_kappa(os.path.join(MX, "mess_eckart_cse.csv"))
rows = []
for (t, p), kappa in sorted(kappa_mess.items()):
    x = kappa_mx[(t, p)]
    rows.append({"T_K": t, "p_bar": p, "kappa_W2_P1_mess": kappa["P1"], "kappa_W2_P1_marxus": x["P1"],
                 "kappa_W2_P2_mess": kappa["P2"], "kappa_W2_P2_marxus": x["P2"]})
write_csv("isomer_bimolecular_kappa_comparison.csv", rows)

for r in rows:
    if r["kappa_W2_P1_mess"] or r["kappa_W2_P2_mess"]:
        print(f"T={r['T_K']:.0f} p={r['p_bar']} bar: kappa(W2,P1) MESS {r['kappa_W2_P1_mess']:.4f} MarXus "
              f"{r['kappa_W2_P1_marxus']:.4f}; kappa(W2,P2) MESS {r['kappa_W2_P2_mess']:.4f} MarXus {r['kappa_W2_P2_marxus']:.4f}")
