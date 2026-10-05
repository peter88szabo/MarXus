#!/usr/bin/env python3
"""ZZ-allyl + O2, Gamma Case 2: MarXus compared with the reference MESS run of 2025-09-29.

Reads
  Gamma-Case2_..._12.7kcal.out              MESS rate tables (high pressure and per (T, p))
  Gamma-Case2_..._12.7kcal.log              MESS log (tunneling correction factors)
  marxus_output/case2_tstlevel_E_*.out      MarXus (run_marxus.sh)
  marxus_output/case2_default_EJ_steady_states.out   MarXus with the PST cores at the EJ level
and writes
  capture_comparison.csv          k_inf(R -> G2) of the phase-space-theory entrance
  high_pressure_comparison.csv    high-pressure rate coefficients of every channel
  net_yields_comparison.csv       long-time shares of P5, escape (ESC), P1, P7 in the net reaction
  apparent_rates_comparison.csv   apparent bimolecular rate coefficients k(R -> X)
  plots/pes.png, plots/p5_share.png, plots/high_pressure_deviation.png, plots/apparent_rates_760torr.png

Long-time fate from the MESS tables: R forms the wells and end channels with the rate coefficients
k(R -> X) of each (T, p) table; every well w then ends in X with the absorption probability of the
first-order network of the wells (rates of the same table, off-diagonal entries), B = (I - Q)^-1 A
(absorbing Markov chain). The shares are normalized to the net reaction P1 + P5 + P7 + ESC = 1. MarXus:
final steady state (Olzmann) with the escape of G4 as its physical sink; the prompt redissociation to R
is excluded by the same normalization.

Run with the science environment:  source ~/.venvs/science/bin/activate && python3 compare_with_mess.py
"""
import csv
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
STEM = "Gamma-Case2_from_ZZ_allyl+O2_Gamma4-Escape_CCSDT_version-1_2025_09_29_12.7kcal"
WELLS = ["G2", "G3", "G4", "G6"]
ENDS = ["R", "P1", "P5", "P7", "ESC"]
CHANNELS = {"B23": [("G2", "G3"), ("G3", "G2")], "B24": [("G2", "G4"), ("G4", "G2")],
            "B34": [("G3", "G4"), ("G4", "G3")], "B36": [("G3", "G6"), ("G6", "G3")],
            "B4P1": [("G4", "P1")], "B4P5": [("G4", "P5")], "B12": [("G2", "R")], "B6P7": [("G6", "P7")]}
os.makedirs(os.path.join(HERE, "plots"), exist_ok=True)


# ----------------------------------------------------------------------------------------------
# Readers
# ----------------------------------------------------------------------------------------------
def read_mess_out(path):
    """High-pressure tables {T: {from: {to: value}}} and (T, p) tables (last column renamed ESC)."""
    lines = open(path).read().split("\n")
    high, pressure = {}, {}
    for i, line in enumerate(lines):
        m_tp = re.match(r"Temperature = (\S+) K\s+Pressure = (\S+) torr", line)
        m_t = re.match(r"Temperature = (\S+) K\s*$", line)
        if not (m_tp or (m_t and i + 2 < len(lines) and "High Pressure" in lines[i + 2])):
            continue
        j = i + 1
        while not lines[j].strip().startswith("From"):
            j += 1
        header = lines[j].split()[1:]
        if m_tp:
            header[-1] = "ESC"
        table = {}
        k = 1
        while j + k < len(lines) and lines[j + k].split() and lines[j + k].split()[0] in WELLS + ENDS:
            fields = lines[j + k].split()
            table[fields[0]] = {h: (np.nan if v == "***" else float(v)) for h, v in zip(header, fields[1:])}
            k += 1
        if m_tp:
            pressure[(float(m_tp.group(1)), float(m_tp.group(2)))] = table
        else:
            high[float(m_t.group(1))] = table
    return high, pressure


def read_marxus_blocks(path):
    """{block title: list of row dicts} of a MarXus example output."""
    blocks, title, header = {}, None, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("# intermediate steady state") or line.startswith("# final steady state") \
                or line.startswith("# bimolecular rate coefficients") or line.startswith("# eigenvalue analysis"):
            title, header = line[2:], None
            blocks[title] = []
        elif line.startswith("#") or not line.strip() or title is None:
            continue
        elif header is None:
            header = line.split(",")
        else:
            blocks[title].append(dict(zip(header, map(float, line.split(",")))))
    return blocks


def block(blocks, prefix):
    for title, rows in blocks.items():
        if title.startswith(prefix):
            return rows
    raise KeyError(prefix)


def mess_long_time_shares(table):
    k_r = {x: table["R"][x] for x in WELLS + ENDS if x != "R"}
    q, a = np.zeros((4, 4)), np.zeros((4, 5))
    for i, w in enumerate(WELLS):
        rates = {x: table[w][x] for x in WELLS + ENDS if x != w}
        loss = sum(rates.values())
        for j, v in enumerate(WELLS):
            if v != w:
                q[i, j] = rates[v] / loss
        for e, x in enumerate(ENDS):
            a[i, e] = rates[x] / loss
    b = np.linalg.solve(np.eye(4) - q, a)
    final = {x: k_r.get(x, 0.0) for x in ENDS}
    for i, w in enumerate(WELLS):
        for e, x in enumerate(ENDS):
            final[x] += k_r[w] * b[i, e]
    net = sum(final[x] for x in ["P1", "P5", "P7", "ESC"])
    return {x: final[x] / net for x in ["P1", "P5", "P7", "ESC"]}


def write_csv(name, rows):
    with open(os.path.join(HERE, name), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            writer.writerow({k: (v if isinstance(v, str) else "%.6g" % v) for k, v in r.items()})


# ----------------------------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------------------------
mess_high, mess_p = read_mess_out(os.path.join(HERE, STEM + ".out"))
mx_e = read_marxus_blocks(os.path.join(HERE, "marxus_output", "case2_tstlevel_E_steady_states.out"))
mx_ej = read_marxus_blocks(os.path.join(HERE, "marxus_output", "case2_default_EJ_steady_states.out"))
mx_eig = read_marxus_blocks(os.path.join(HERE, "marxus_output", "case2_tstlevel_E_eigenvalue.out"))
temperatures = sorted(mess_high)
pressures = sorted({p for _, p in mess_p})

# 1. Capture (entrance) rate coefficient.
capture = []
k_inf_e = {r["T[K]"]: r["k_inf"] for r in block(mx_e, "bimolecular rate coefficients")}
k_inf_ej = {r["T[K]"]: r["k_inf"] for r in block(mx_ej, "bimolecular rate coefficients")}
for t in temperatures:
    m = mess_high[t]["R"]["G2"]
    capture.append({"T_K": t, "mess": m, "marxus_E": k_inf_e[t], "dev_E_percent": 100 * (k_inf_e[t] / m - 1),
                    "marxus_EJ": k_inf_ej[t], "dev_EJ_percent": 100 * (k_inf_ej[t] / m - 1)})
write_csv("capture_comparison.csv", capture)

# 2. High-pressure rate coefficients of every channel.
eig = {r["T[K]"]: r for r in block(mx_eig, "eigenvalue analysis")}
high = []
for t in temperatures:
    for c, pairs in CHANNELS.items():
        for a, b in pairs:
            m = mess_high[t][a].get(b, np.nan)
            x = eig[t].get(f"k_inf({a}:{c})[1/s]", np.nan)
            high.append({"T_K": t, "channel": c, "from": a, "to": b, "mess": m, "marxus": x,
                         "dev_percent": 100 * (x / m - 1)})
write_csv("high_pressure_comparison.csv", high)

# 3. Long-time shares of the net reaction.
final = {(r["T[K]"], r["P[Torr]"]): r for r in block(mx_e, "final steady state")}
shares = []
for t in temperatures:
    for p in pressures:
        ms = mess_long_time_shares(mess_p[(t, p)])
        f = final[(t, p)]
        xs = {"P1": f["Phi(G4:B4P1)"], "P5": f["Phi(G4:B4P5)"], "P7": f["Phi(G6:B6P7)"], "ESC": f["Phi_sink(G4)"]}
        net = sum(xs.values())
        row = {"T_K": t, "p_torr": p}
        for x in ["P5", "ESC", "P1", "P7"]:
            row[f"mess_{x}"] = ms[x]
            row[f"marxus_{x}"] = xs[x] / net
            row[f"dev_{x}_percent"] = 100 * (xs[x] / net / ms[x] - 1)
        row["marxus_back_to_R"] = f["Phi(G2:B12)"]
        shares.append(row)
write_csv("net_yields_comparison.csv", shares)

# 4. Apparent bimolecular rate coefficients (intermediate steady state) vs the MESS R row.
apparent = []
bimol = {(r["T[K]"], r["P[Torr]"]): r for r in block(mx_e, "bimolecular rate coefficients")}
pairs = [("G2", "k(R->G2)"), ("G3", "k(R->G3)"), ("G4", "k(R->G4)"), ("P5", "k(R->P5 via B4P5)"),
         ("ESC", "k(R->sink of G4)"), ("P1", "k(R->P1 via B4P1)"), ("P7", "k(R->P7 via B6P7)")]
for t in temperatures:
    for p in pressures:
        row = {"T_K": t, "p_torr": p}
        for x, col in pairs:
            row[f"mess_{x}"] = mess_p[(t, p)]["R"][x]
            row[f"marxus_{x}"] = bimol[(t, p)][col]
        apparent.append(row)
write_csv("apparent_rates_comparison.csv", apparent)

# ----------------------------------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------------------------------
# PES diagram (energies relative to R, from the header of the MESS output).
energies = {"R": 0.0, "G2": -16.3, "G3": -18.7, "G4": -19.4, "G6": -19.7, "P1": -5.4, "P5": -27.4, "P7": -14.9}
barriers = {"B12": ("R", "G2", 0.0), "B23": ("G2", "G3", 0.0), "B24": ("G2", "G4", -4.6), "B34": ("G3", "G4", 0.4),
            "B36": ("G3", "G6", -1.4), "B4P1": ("G4", "P1", 1.6), "B4P5": ("G4", "P5", -3.6), "B6P7": ("G6", "P7", -4.7)}
x = {"R": 0, "G2": 2, "G3": 4, "G4": 6, "G6": 8, "P1": 8, "P5": 10, "P7": 10}
xp = {"P1": (6.0, 9.0), "P7": (8.0, 11.0)}
fig, ax = plt.subplots(figsize=(10, 5.5))
labels = {"R": "R\nZZ-allyl + O$_2$", "G2": "G2", "G3": "G3", "G4": "G4\n(escape 2.5e7 s$^{-1}$)", "G6": "G6",
          "P1": "P1 HPALD + HO$_2$", "P5": "P5 IEPOX + OH", "P7": "P7 + HO$_2$"}
pos = {"R": 0, "G2": 2, "G3": 4, "G4": 6, "G6": 8, "P1": 7.5, "P5": 10, "P7": 10.5}
pos["P1"], pos["P7"] = 7.4, 10.6
pos["P5"] = 9.0
for s, e in energies.items():
    ax.hlines(e, pos[s] - 0.35, pos[s] + 0.35, color="k" if s[0] != "P" else "tab:green", lw=3)
    ax.text(pos[s], e - 1.2, f"{labels[s]}\n{e:.1f}", ha="center", va="top", fontsize=8)
for name, (a, b, e) in barriers.items():
    xm = 0.5 * (pos[a] + pos[b])
    color = "tab:blue" if name in ("B12", "B6P7") else "firebrick"
    ax.hlines(e, xm - 0.25, xm + 0.25, color=color, lw=2)
    ax.plot([pos[a] + 0.35, xm - 0.25], [energies[a], e], ":", color="gray", lw=0.8)
    ax.plot([xm + 0.25, pos[b] - 0.35], [e, energies[b]], ":", color="gray", lw=0.8)
    ax.text(xm, e + 0.6, f"{name}\n{e:+.1f}", ha="center", va="bottom", fontsize=7, color=color)
ax.set_ylabel("energy relative to R (kcal/mol)")
ax.set_xticks([])
ax.set_ylim(-32, 6)
ax.set_title("ZZ-allyl + O$_2$, Gamma Case 2: wells, barriers (red: Eckart tunneling; blue: phase-space theory), products")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "pes.png"), dpi=200)
plt.close(fig)

# P5 share of the net reaction.
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
colors = plt.cm.plasma(np.linspace(0.0, 0.8, len(pressures)))
for p, c in zip(pressures, colors):
    sel = [r for r in shares if r["p_torr"] == p]
    axes[0].plot([r["T_K"] for r in sel], [100 * r["mess_P5"] for r in sel], "-o", color=c, mfc="none", label=f"MESS {p:g} Torr")
    axes[0].plot([r["T_K"] for r in sel], [100 * r["marxus_P5"] for r in sel], "s", color=c, label=f"MarXus {p:g} Torr")
    axes[1].plot([r["T_K"] for r in sel], [r["dev_P5_percent"] for r in sel], "-s", color=c, label=f"P5, {p:g} Torr")
    axes[1].plot([r["T_K"] for r in sel], [r["dev_ESC_percent"] for r in sel], "--^", color=c, mfc="none", label=f"escape, {p:g} Torr")
axes[0].set_ylabel("P5 (IEPOX + OH) share of the net reaction (%)")
axes[0].set_title("long-time P5 yield")
axes[1].axhline(0, color="k", lw=0.8)
axes[1].set_ylabel("MarXus / MESS - 1 (%)")
axes[1].set_title("deviation of the shares")
for ax in axes:
    ax.set_xlabel("T (K)")
    ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "p5_share.png"), dpi=200)
plt.close(fig)

# High-pressure deviations per channel.
fig, ax = plt.subplots(figsize=(8.5, 4.8))
for c, pairs_c in CHANNELS.items():
    a, b = pairs_c[0]
    sel = [r for r in high if r["channel"] == c and r["from"] == a]
    ax.plot([r["T_K"] for r in sel], [r["dev_percent"] for r in sel], "-o", label=f"{c} ({a}->{b})")
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("T (K)")
ax.set_ylabel("k$_\\infty$: MarXus / MESS - 1 (%)")
ax.set_title("High-pressure rate coefficients: tunneling barriers differ by the Eckart factor")
ax.legend(fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "high_pressure_deviation.png"), dpi=200)
plt.close(fig)

# Apparent rate coefficients at 760 Torr.
fig, ax = plt.subplots(figsize=(8.5, 5))
sel = [r for r in apparent if r["p_torr"] == 760.0]
for (x, _), c in zip(pairs, plt.cm.tab10(np.arange(len(pairs)))):
    m = [r[f"mess_{x}"] for r in sel]
    ax.plot([r["T_K"] for r in sel], [abs(v) for v in m], "-o", color=c, mfc="none", label=f"MESS R->{x}")
    ax.plot([r["T_K"] for r in sel], [r[f"marxus_{x}"] for r in sel], "s", color=c, label=f"MarXus R->{x}")
ax.set_yscale("log")
ax.set_xlabel("T (K)")
ax.set_ylabel("k(R -> X) (cm$^3$ s$^{-1}$)")
ax.set_title("Apparent rate coefficients at 760 Torr (MarXus: intermediate steady state; |MESS| where negative)")
ax.legend(fontsize=7, ncol=2)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "apparent_rates_760torr.png"), dpi=200)
plt.close(fig)

print("written:", ", ".join(sorted(os.listdir(os.path.join(HERE, "plots")))), "and the four CSV tables")
for r in capture[:1] + capture[-1:]:
    print(f"capture T={r['T_K']:.0f}: E {r['dev_E_percent']:+.2f}%  EJ {r['dev_EJ_percent']:+.2f}%")
for r in shares:
    if r["p_torr"] == 760.0:
        print(f"T={r['T_K']:.0f} 760 Torr: P5 MESS {100*r['mess_P5']:.3f}%  MarXus {100*r['marxus_P5']:.3f}%  ({r['dev_P5_percent']:+.1f}%)")
