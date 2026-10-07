"""Comparison of MarXus with MESS for C2H4 + HO2 (P1) -> CH2CH2OOH (W2) -> OH + oxirane (P2), seven hindered rotors.

Inputs:
  reference_mess/c2h4_ho2.out, .log      MESS rate tables and rotor data (OutPrecision, LogPrecision 6)
  marxus_output/<variant>_<method>.out   MarXus reports (rotor lines), variant default (exact Eckart) or
                                         mess_eckart (MESS Eckart model)
  marxus_output/<variant>_<method>_tables.csv
Outputs:
  rotor_comparison.csv            B, ground energy, number of levels and the nine lowest levels of every rotor
  high_pressure_comparison.csv    k_inf of W2 -> P1, W2 -> P2 and the capture P1 -> W2
  rate_comparison.csv             CSE rate coefficients at every (T, p), and the net P1 -> P2 where MESS merges W2
  plots/high_pressure_deviation.png, plots/rate_deviation.png
Run with the science environment: source ~/.venvs/science/bin/activate; python compare_with_mess.py
"""
import csv
import os
import re

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
MESS_OUT = os.path.join(HERE, "reference_mess", "c2h4_ho2.out")
MESS_LOG = os.path.join(HERE, "reference_mess", "c2h4_ho2.log")
MX = os.path.join(HERE, "marxus_output")
BAR_TO_TORR = 750.0616827041697
VARIANTS = {"default": "exact Eckart", "mess_eckart": "MESS Eckart model"}
os.makedirs(os.path.join(HERE, "plots"), exist_ok=True)


def write_csv(name, rows):
    with open(os.path.join(HERE, name), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            writer.writerow({k: (v if isinstance(v, str) else "%.6g" % v) for k, v in r.items()})


def read_mess_out(path):
    """High-pressure tables {T: {from: {to: k}}} and (T, p_bar) tables; '***' is NaN."""
    lines = open(path).read().split("\n")
    high, pressure = {}, {}
    for i, line in enumerate(lines):
        m_tp = re.match(r"Temperature = (\S+) K\s+Pressure = (\S+) bar", line)
        m_t = re.match(r"Temperature = (\S+) K\s*$", line)
        if not (m_tp or (m_t and i + 2 < len(lines) and "High Pressure" in lines[i + 2])):
            continue
        j = i + 1
        while not lines[j].strip().startswith("From"):
            j += 1
        header = lines[j].split()[1:]
        table = {}
        for row in lines[j + 1 : j + 1 + len(header)]:
            fields = row.split()
            table[fields[0]] = {h: (np.nan if v == "***" else float(v)) for h, v in zip(header, fields[1:])}
        if m_tp:
            pressure[(float(m_tp.group(1)), float(m_tp.group(2)))] = table
        else:
            high[float(m_t.group(1))] = table
    return high, pressure


def read_mess_rotors(path):
    """Per rotor in the order of the log: B (1/cm), ground energy (kcal/mol), number of levels, lowest levels."""
    rotors, current = [], None
    for line in open(path):
        if "effective rotational constant[1/cm]" in line:
            current = {"B": float(line.split("=")[1])}
        elif current is not None and "ground energy [kcal/mol]" in line:
            current["ground"] = float(line.split("=")[1])
        elif current is not None and "number of levels" in line:
            current["levels"] = int(line.split("=")[1])
        elif current is not None and "lowest excited states" in line:
            current["lowest"] = [float(x) for x in line.split(":")[1].split()]
            rotors.append(current)
            current = None
    return rotors


def read_marxus_rotors(path):
    pattern = re.compile(
        r"(\S+) rotor (\d+): B (\S+) cm-1, symmetry (\d+), potential (\S+) \.\. (\S+) kcal/mol, ground energy (\S+) "
        r"kcal/mol, (\d+) levels up to (\S+) kcal/mol above the ground; lowest above the ground \(kcal/mol\): (.*)$"
    )
    rotors = []
    for line in open(path):
        m = pattern.search(line)
        if m:
            rotors.append({"species": m.group(1), "rotor": int(m.group(2)), "B": float(m.group(3)),
                           "ground": float(m.group(7)), "levels": int(m.group(8)),
                           "lowest": [float(x) for x in m.group(10).split()]})
    return rotors


def read_tables(path):
    """{title: list of row dicts} of a MarXus tables file; empty fields are NaN."""
    blocks, title, header = {}, None, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("# "):
            title, header = line[2:], None
            blocks[title] = []
        elif not line.strip() or title is None:
            continue
        elif header is None:
            header = line.split(",")
        else:
            blocks[title].append({h: (float(v) if v else np.nan) for h, v in zip(header, line.split(","))})
    return blocks


def block(blocks, prefix):
    for title, rows in blocks.items():
        if title.startswith(prefix):
            return {(r["T[K]"], round(r["P[Torr]"], 6)): r for r in rows}
    raise KeyError(prefix)


mess_high, mess_p = read_mess_out(MESS_OUT)

# 1. Rotors
mess_rotors = read_mess_rotors(MESS_LOG)
marxus_rotors = read_marxus_rotors(os.path.join(MX, "default_cse.out"))
assert len(mess_rotors) == len(marxus_rotors) == 7, (len(mess_rotors), len(marxus_rotors))
rotor_rows = []
for m, x in zip(mess_rotors, marxus_rotors):
    rotor_rows.append({
        "species": x["species"], "rotor": x["rotor"],
        "B_mess_cm1": m["B"], "B_marxus_cm1": x["B"], "B_dev_rel": x["B"] / m["B"] - 1,
        "ground_mess_kcal": m["ground"], "ground_marxus_kcal": x["ground"], "ground_diff_kcal": x["ground"] - m["ground"],
        "levels_mess": m["levels"], "levels_marxus": x["levels"],
        "lowest9_max_abs_diff_kcal": max(abs(a - b) for a, b in zip(m["lowest"], x["lowest"])),
    })
write_csv("rotor_comparison.csv", rotor_rows)

# 2. High-pressure rate coefficients
tables = {v: {m: read_tables(os.path.join(MX, f"{v}_{m}_tables.csv")) for m in ("olzmann", "cse")} for v in VARIANTS}
high_rows = []
for t in sorted(mess_high):
    row = {"T_K": t}
    for v in VARIANTS:
        k_inf = block(tables[v]["olzmann"], "thermal rate coefficients: High-pressure")
        capture = block(tables[v]["cse"], "CSE: Capture, return")
        first = next(r for (tt, _), r in sorted(k_inf.items()) if tt == t)
        cap = next(r for (tt, _), r in sorted(capture.items()) if tt == t)
        for name, mess, marxus in [
            ("W2->P1", mess_high[t]["W2"]["P1"], first["W2->P1"]),
            ("W2->P2", mess_high[t]["W2"]["P2"], first["W2->P2"]),
            ("P1->W2", mess_high[t]["P1"]["W2"], cap["capture"]),
        ]:
            row[f"{name}_mess"] = mess
            row[f"{name}_{v}"] = marxus
            row[f"{name}_{v}_dev_percent"] = 100 * (marxus / mess - 1)
    high_rows.append(row)
write_csv("high_pressure_comparison.csv", high_rows)

# 3. Rate coefficients at every (T, p)
rate_rows = []
for (t, p_bar), mess in sorted(mess_p.items()):
    key = (t, round(p_bar * BAR_TO_TORR, 6))
    mess_kept = not np.isnan(mess["W2"]["P1"])
    for v in VARIANTS:
        w2 = block(tables[v]["cse"], "CSE: Rate coefficients from W2")[key]
        p1p2 = block(tables[v]["cse"], "CSE: Bimolecular-to-bimolecular rate coefficients")[key]
        p1w2 = block(tables[v]["cse"], "CSE: Bimolecular-to-well rate coefficients")[key]
        net = block(tables[v]["cse"], "CSE: Capture, return")[key]
        marxus_kept = not np.isnan(w2["W2->P1"])
        status = {(True, True): "both keep W2", (False, False): "both merge W2",
                  (False, True): "MESS merges W2", (True, False): "MarXus merges W2"}[(mess_kept, marxus_kept)]
        quantities = [("net P1->P2", mess["P1"]["P2"] if not mess_kept else mess["P1"]["P1"], net["net"])]
        if mess_kept and marxus_kept:
            quantities += [("W2->P1", mess["W2"]["P1"], w2["W2->P1"]), ("W2->P2", mess["W2"]["P2"], w2["W2->P2"]),
                           ("P1->W2", mess["P1"]["W2"], p1w2["P1->W2"]), ("P1->P2", mess["P1"]["P2"], p1p2["P1->P2"])]
        for name, m, x in quantities:
            rate_rows.append({"T_K": t, "p_bar": p_bar, "variant": v, "status": status, "quantity": name,
                              "mess": m, "marxus": x, "dev_percent": 100 * (x / m - 1)})
write_csv("rate_comparison.csv", rate_rows)

# 4. Plots
fig, ax = plt.subplots(figsize=(6.5, 4.2))
temps = [r["T_K"] for r in high_rows]
for name, marker in [("W2->P1", "o"), ("W2->P2", "s"), ("P1->W2", "^")]:
    for v, style in [("default", "-"), ("mess_eckart", "--")]:
        ax.plot(temps, [r[f"{name}_{v}_dev_percent"] for r in high_rows], style, marker=marker,
                label=f"{name}, {VARIANTS[v]}")
ax.axhline(0, color="k", lw=0.6)
ax.set_xlabel("T (K)")
ax.set_ylabel("MarXus / MESS − 1 (%)")
ax.set_title("High-pressure rate coefficients, C$_2$H$_4$ + HO$_2$ (7 hindered rotors)")
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "high_pressure_deviation.png"), dpi=150)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
for ax, v in zip(axes, VARIANTS):
    for name, marker in [("net P1->P2", "o"), ("P1->W2", "^"), ("P1->P2", "s"), ("W2->P2", "D"), ("W2->P1", "v")]:
        for t in sorted({r["T_K"] for r in rate_rows}):
            pts = [(r["p_bar"], r["dev_percent"]) for r in rate_rows if r["variant"] == v and r["quantity"] == name
                   and r["T_K"] == t and r["status"] != "MESS merges W2"]
            if pts:
                ax.semilogx(*zip(*pts), marker=marker, ls="-", lw=0.8, label=f"{name}, {t:.0f} K")
    merged = [(r["p_bar"], r["dev_percent"]) for r in rate_rows if r["variant"] == v and r["status"] == "MESS merges W2"]
    ax.semilogx(*zip(*merged), "kx", ms=8, ls="none", label="net P1->P2 where only MESS merges W2 (other quantity)")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("p (bar)")
    ax.set_title(f"MarXus CSE vs MESS ({VARIANTS[v]})")
axes[0].set_ylabel("MarXus / MESS − 1 (%)")
axes[1].legend(fontsize=5.5, ncol=2)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "rate_deviation.png"), dpi=150)

# Summary
print("rotors: max |B dev| %.2e, max |ground diff| %.1e kcal/mol, max |level diff| %.1e kcal/mol, levels equal: %s" % (
    max(abs(r["B_dev_rel"]) for r in rotor_rows), max(abs(r["ground_diff_kcal"]) for r in rotor_rows),
    max(r["lowest9_max_abs_diff_kcal"] for r in rotor_rows),
    all(r["levels_mess"] == r["levels_marxus"] for r in rotor_rows)))
for r in high_rows:
    print("k_inf T=%5.0f K: " % r["T_K"] + "  ".join(
        "%s %+.2f%% / %+.2f%%" % (n, r[f"{n}_default_dev_percent"], r[f"{n}_mess_eckart_dev_percent"])
        for n in ("W2->P1", "W2->P2", "P1->W2")))
for v in VARIANTS:
    for status in ("both keep W2", "MESS merges W2", "both merge W2", "MarXus merges W2"):
        for q in ("net P1->P2", "W2->P1", "W2->P2", "P1->W2", "P1->P2"):
            d = [r["dev_percent"] for r in rate_rows if r["variant"] == v and r["status"] == status and r["quantity"] == q]
            if d:
                print(f"{VARIANTS[v]:>18} | {status:<16} | {q:<10}: {min(d):+7.2f} … {max(d):+7.2f} % ({len(d)})")
