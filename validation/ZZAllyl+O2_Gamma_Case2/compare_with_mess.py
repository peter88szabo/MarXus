#!/usr/bin/env python3
"""ZZ-allyl + O2, Gamma Case 2: MarXus compared with the reference MESS run of 2025-09-29.

Reads
  Gamma-Case2_..._12.7kcal.out              MESS rate tables (high pressure and per (T, p))
  Gamma-Case2_..._12.7kcal.log              MESS log (tunneling correction factors)
  marxus_output/case2_tstlevel_E_*.out      MarXus (run_marxus.sh)
  marxus_output/case2_default_EJ_steady_states.out   MarXus with the PST cores at the EJ level
  marxus_output/eckart_kappa.csv            MarXus canonical Eckart factors kappa(T): exact and MESS model
  marxus_output/case2_tstlevel_E_mess_eckart_*.out   MarXus with the MESS Eckart tunneling model, including
                                            the CSE species tables (case2_tstlevel_E_mess_eckart_cse.out)
and writes
  capture_comparison.csv          k_inf(R -> G2) of the phase-space-theory entrance
  high_pressure_comparison.csv    high-pressure rate coefficients of every channel
  net_yields_comparison.csv       long-time shares of P5, escape (ESC), P1, P7 in the net reaction
  apparent_rates_comparison.csv   apparent bimolecular rate coefficients k(R -> X)
  kappa_comparison.csv            tunneling factors kappa(T): MESS log vs MarXus (exact Eckart)
  cse_comparison.csv              every species-to-species rate coefficient: MESS vs the MarXus CSE method
  plots/cse_vs_mess.png           CSE method vs MESS: reactant row at 760 Torr and deviations of all entries
  plots/pes.png                   the network
  plots/iepox_oh_yield.png        P5 (IEPOX + OH): long-time share, deviation, apparent k(R -> P5)
  plots/tunneling_ratio_bars.png  kappa(MarXus)/kappa(MESS) per barrier at 200, 300, 400 K; kappa at 300 K
  plots/high_pressure_deviation.png, plots/apparent_rates_760torr.png

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
                or line.startswith("# bimolecular rate coefficients") \
                or line.startswith("# thermal rate coefficients of the final steady state"):
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


def read_cse(path):
    """MarXus CSE species tables: {(T, p_torr): {from: {to: value}}} (escape(G4) renamed ESC)."""
    out, key, header = {}, None, None
    for line in open(path):
        line = line.rstrip("\n")
        m = re.match(r"# T = (\S+) K, p = (\S+) Torr", line)
        if m:
            key, header = (float(m.group(1)), float(m.group(2))), None
            out[key] = {}
            continue
        if line.startswith("From\\To"):
            header = [("ESC" if h == "escape(G4)" else h) for h in line.split(",")[1:]]
            continue
        if header and line and not line.startswith("#"):
            fields = line.split(",")
            out[key][fields[0]] = dict(zip(header, map(float, fields[1:])))
    return out


def write_csv(name, rows):
    with open(os.path.join(HERE, name), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            writer.writerow({k: (v if isinstance(v, str) else "%.6g" % v) for k, v in r.items()})


def read_mess_kappa(path):
    """Tunneling correction factors of the MESS log: {barrier: {T: kappa}}."""
    lines = open(path).read().split("\n")
    start = next(i for i, l in enumerate(lines) if "tunneling partition function correction factors" in l)
    header = lines[start + 1].split()[1:]          # B23 D B24 D ...
    names = header[0::2]
    kappa = {n: {} for n in names}
    for line in lines[start + 2:]:
        fields = line.split()
        if not fields or not re.match(r"^\d+$", fields[0]):
            break
        t = float(fields[0])
        for k, n in enumerate(names):
            kappa[n][t] = float(fields[1 + 2 * k])
    return kappa


def read_marxus_kappa(path, column="kappa"):
    kappa = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            kappa.setdefault(row["barrier"], {})[float(row["T[K]"])] = float(row[column])
    return kappa


# ----------------------------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------------------------
mess_high, mess_p = read_mess_out(os.path.join(HERE, STEM + ".out"))
mx_e = read_marxus_blocks(os.path.join(HERE, "marxus_output", "case2_tstlevel_E_steady_states.out"))
mx_ej = read_marxus_blocks(os.path.join(HERE, "marxus_output", "case2_default_EJ_steady_states.out"))
mx_me = read_marxus_blocks(os.path.join(HERE, "marxus_output", "case2_tstlevel_E_mess_eckart_steady_states.out"))
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
# Thermal rate coefficients of the final steady state (lowest eigenpair of J, GO10 eq. 12).
eig = {r["T[K]"]: r for r in block(mx_e, "thermal rate coefficients of the final steady state")}
eig_me = {r["T[K]"]: r for r in block(mx_me, "thermal rate coefficients of the final steady state")}
high = []
for t in temperatures:
    for c, pairs in CHANNELS.items():
        for a, b in pairs:
            m = mess_high[t][a].get(b, np.nan)
            x = eig[t].get(f"k_inf({a}:{c})[1/s]", np.nan)
            y = eig_me[t].get(f"k_inf({a}:{c})[1/s]", np.nan)
            high.append({"T_K": t, "channel": c, "from": a, "to": b, "mess": m, "marxus": x,
                         "dev_percent": 100 * (x / m - 1), "marxus_mess_eckart": y,
                         "dev_mess_eckart_percent": 100 * (y / m - 1)})
write_csv("high_pressure_comparison.csv", high)

# 3. Long-time shares of the net reaction.
final = {(r["T[K]"], r["P[Torr]"]): r for r in block(mx_e, "final steady state")}
final_me = {(r["T[K]"], r["P[Torr]"]): r for r in block(mx_me, "final steady state")}


def marxus_shares(f):
    xs = {"P1": f["Phi(G4:B4P1)"], "P5": f["Phi(G4:B4P5)"], "P7": f["Phi(G6:B6P7)"], "ESC": f["Phi_sink(G4)"]}
    net = sum(xs.values())
    return {x: v / net for x, v in xs.items()}


shares = []
for t in temperatures:
    for p in pressures:
        ms = mess_long_time_shares(mess_p[(t, p)])
        xs = marxus_shares(final[(t, p)])
        ys = marxus_shares(final_me[(t, p)])
        row = {"T_K": t, "p_torr": p}
        for x in ["P5", "ESC", "P1", "P7"]:
            row[f"mess_{x}"] = ms[x]
            row[f"marxus_{x}"] = xs[x]
            row[f"dev_{x}_percent"] = 100 * (xs[x] / ms[x] - 1)
            row[f"marxus_mess_eckart_{x}"] = ys[x]
            row[f"dev_mess_eckart_{x}_percent"] = 100 * (ys[x] / ms[x] - 1)
        row["marxus_back_to_R"] = final[(t, p)]["Phi(G2:B12)"]
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

# 5. Tunneling factors.
kappa_mess = read_mess_kappa(os.path.join(HERE, STEM + ".log"))
kappa_mx = read_marxus_kappa(os.path.join(HERE, "marxus_output", "eckart_kappa.csv"))
kappa_mimic = read_marxus_kappa(os.path.join(HERE, "marxus_output", "eckart_kappa.csv"), "kappa_mess")
tunneling_barriers = [b for b in ["B23", "B24", "B34", "B36", "B4P1", "B4P5"] if b in kappa_mx]
kappa_rows = []
for b in tunneling_barriers:
    for t in sorted(kappa_mess[b]):
        if t in kappa_mx[b]:
            kappa_rows.append({"barrier": b, "T_K": t, "kappa_mess": kappa_mess[b][t], "kappa_marxus": kappa_mx[b][t],
                               "ratio": kappa_mx[b][t] / kappa_mess[b][t],
                               "kappa_marxus_mess_model": kappa_mimic[b][t],
                               "ratio_mess_model": kappa_mimic[b][t] / kappa_mess[b][t]})
write_csv("kappa_comparison.csv", kappa_rows)

# 6. CSE species tables vs MESS.
cse = read_cse(os.path.join(HERE, "marxus_output", "case2_tstlevel_E_mess_eckart_cse.out"))
cse_rows = []
for (t, p), table in sorted(cse.items()):
    for a in WELLS + ["R"]:
        for b in WELLS + ENDS:
            if a == b == "R" or b not in table.get(a, {}) or b not in mess_p[(t, p)][a]:
                continue
            m, x = mess_p[(t, p)][a][b], table[a][b]
            cse_rows.append({"T_K": t, "p_torr": p, "from": a, "to": b, "mess": m, "marxus_cse": x,
                             "dev_percent": 100 * (x / m - 1) if m != 0 else np.nan})
write_csv("cse_comparison.csv", cse_rows)

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

# IEPOX + OH (P5): long-time share of the net reaction, its deviation, and the apparent k(R -> P5).
fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
colors = plt.cm.plasma(np.linspace(0.0, 0.8, len(pressures)))
for p, c in zip(pressures, colors):
    sel = [r for r in shares if r["p_torr"] == p]
    ts = [r["T_K"] for r in sel]
    axes[0].plot(ts, [100 * r["mess_P5"] for r in sel], "-o", color=c, mfc="none", label=f"MESS {p:g} Torr")
    axes[0].plot(ts, [100 * r["marxus_P5"] for r in sel], "s", color=c, label=f"MarXus exact Eckart {p:g} Torr")
    axes[0].plot(ts, [100 * r["marxus_mess_eckart_P5"] for r in sel], "^", color=c, label=f"MarXus MESS Eckart {p:g} Torr")
    axes[1].plot(ts, [r["dev_P5_percent"] for r in sel], "-s", color=c, label=f"exact Eckart, {p:g} Torr")
    axes[1].plot(ts, [r["dev_mess_eckart_P5_percent"] for r in sel], "--^", color=c, label=f"MESS Eckart, {p:g} Torr")
    sel_a = [r for r in apparent if r["p_torr"] == p]
    axes[2].plot([r["T_K"] for r in sel_a], [r["mess_P5"] for r in sel_a], "-o", color=c, mfc="none", label=f"MESS {p:g} Torr")
    axes[2].plot([r["T_K"] for r in sel_a], [r["marxus_P5"] for r in sel_a], "s", color=c, label=f"MarXus {p:g} Torr")
axes[0].set_ylabel("IEPOX + OH share of the net reaction (%)")
axes[0].set_title("long-time IEPOX + OH yield\n(MarXus: final steady state; MESS: from its rate tables)", fontsize=10)
axes[1].axhline(0, color="k", lw=0.8)
axes[1].set_ylabel("MarXus / MESS - 1 (%)")
axes[1].set_title("deviation of the IEPOX + OH share from MESS", fontsize=10)
axes[2].set_yscale("log")
axes[2].set_ylabel("k(R -> IEPOX + OH) (cm$^3$ s$^{-1}$)")
axes[2].set_title("apparent rate coefficient\n(MarXus: intermediate steady state)", fontsize=10)
for ax in axes:
    ax.set_xlabel("T (K)")
    ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "iepox_oh_yield.png"), dpi=200)
plt.close(fig)

# Tunneling factors: ratio MarXus/MESS per barrier (grouped bars) and kappa at 300 K.
fig, axes = plt.subplots(1, 3, figsize=(18, 4.8))
group_t = [t for t in (200.0, 300.0, 400.0) if all(t in kappa_mess[b] and t in kappa_mx[b] for b in tunneling_barriers)]
width = 0.8 / len(group_t)
xs = np.arange(len(tunneling_barriers))
for k, (t, c) in enumerate(zip(group_t, ["tab:blue", "tab:orange", "tab:green"])):
    ratios = [kappa_mx[b][t] / kappa_mess[b][t] for b in tunneling_barriers]
    bars = axes[0].bar(xs + (k - (len(group_t) - 1) / 2) * width, [100 * (r - 1) for r in ratios], width, color=c, label=f"{t:.0f} K")
    for bar, r in zip(bars, ratios):
        axes[0].text(bar.get_x() + bar.get_width() / 2, 100 * (r - 1), f"{100 * (r - 1):+.0f}%", ha="center", va="bottom", fontsize=6)
axes[0].set_xticks(xs)
axes[0].set_xticklabels([f"{b}\n{CHANNELS[b][0][0]}->{CHANNELS[b][0][1]}" for b in tunneling_barriers])
axes[0].axhline(0, color="k", lw=0.8)
axes[0].set_ylabel("$\\kappa$(MarXus) / $\\kappa$(MESS) - 1 (%)")
axes[0].set_title("Eckart tunneling factor: MarXus (exact Eckart) vs MESS log", fontsize=10)
axes[0].legend()
k_m = [kappa_mess[b][300.0] for b in tunneling_barriers]
k_x = [kappa_mx[b][300.0] for b in tunneling_barriers]
axes[1].bar(xs - 0.2, k_m, 0.4, color="tab:gray", label="MESS")
axes[1].bar(xs + 0.2, k_x, 0.4, color="tab:red", label="MarXus")
axes[1].set_yscale("log")
axes[1].set_xticks(xs)
axes[1].set_xticklabels(tunneling_barriers)
axes[1].set_ylabel("$\\kappa$(300 K)")
axes[1].set_title("tunneling factor at 300 K", fontsize=10)
axes[1].legend()
for k, (t, c) in enumerate(zip(group_t, ["tab:blue", "tab:orange", "tab:green"])):
    ratios = [kappa_mimic[b][t] / kappa_mess[b][t] for b in tunneling_barriers]
    axes[2].bar(xs + (k - (len(group_t) - 1) / 2) * width, [100 * (r - 1) for r in ratios], width, color=c, label=f"{t:.0f} K")
axes[2].set_xticks(xs)
axes[2].set_xticklabels(tunneling_barriers)
axes[2].axhline(0, color="k", lw=0.8)
axes[2].set_ylabel("$\\kappa$(MarXus MESS model) / $\\kappa$(MESS log) - 1 (%)")
axes[2].set_title("reproduction of the MESS factors by mess_eckart_tunneling", fontsize=10)
axes[2].legend()
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "tunneling_ratio_bars.png"), dpi=200)
plt.close(fig)

# High-pressure deviations per channel, both tunneling models.
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
for ax, key, title in ((axes[0], "dev_percent", "exact Eckart tunneling (MarXus default)"),
                       (axes[1], "dev_mess_eckart_percent", "MESS Eckart tunneling model")):
    for c, pairs_c in CHANNELS.items():
        a, b = pairs_c[0]
        sel = [r for r in high if r["channel"] == c and r["from"] == a]
        ax.plot([r["T_K"] for r in sel], [r[key] for r in sel], "-o", label=f"{c} ({a}->{b})")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("T (K)")
    ax.set_title(title, fontsize=10)
axes[0].set_ylabel("k$_\\infty$: MarXus / MESS - 1 (%)")
axes[1].legend(fontsize=8, ncol=2)
fig.suptitle("High-pressure rate coefficients of every channel")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "high_pressure_deviation.png"), dpi=200)
plt.close(fig)

# CSE method vs MESS.
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
sel = [r for r in cse_rows if r["from"] == "R" and r["p_torr"] == 760.0]
for b, c in zip(["G2", "G3", "G4", "P5", "P1", "P7"], plt.cm.tab10(np.arange(6))):
    pts = [r for r in sel if r["to"] == b]
    axes[0].plot([r["T_K"] for r in pts], [r["mess"] for r in pts], "-o", color=c, mfc="none", label=f"MESS R->{b}")
    axes[0].plot([r["T_K"] for r in pts], [r["marxus_cse"] for r in pts], "s", color=c, label=f"MarXus R->{b}")
axes[0].set_yscale("log")
axes[0].set_xlabel("T (K)")
axes[0].set_ylabel("k(R -> X) (cm$^3$ s$^{-1}$), 760 Torr")
axes[0].set_title("reactant row: MarXus CSE (MESS Eckart, TST level E) vs MESS", fontsize=10)
axes[0].legend(fontsize=7, ncol=2)
# Deviations of the significant entries (|k| above 1e-6 of the largest entry of its row).
significant = []
for r in cse_rows:
    row_max = max(abs(x["mess"]) for x in cse_rows if x["T_K"] == r["T_K"] and x["p_torr"] == r["p_torr"] and x["from"] == r["from"])
    if abs(r["mess"]) > 1e-6 * row_max and not np.isnan(r["dev_percent"]):
        significant.append(r)
pairs_sig = sorted({(r["from"], r["to"]) for r in significant})
for k, (a, b) in enumerate(pairs_sig):
    vals = [r["dev_percent"] for r in significant if (r["from"], r["to"]) == (a, b)]
    axes[1].plot([k] * len(vals), vals, "o", ms=3, color="tab:blue" if a != "R" else "tab:red")
axes[1].set_xticks(range(len(pairs_sig)))
axes[1].set_xticklabels([f"{a}->{b}" for a, b in pairs_sig], rotation=90, fontsize=7)
axes[1].axhline(0, color="k", lw=0.8)
axes[1].set_ylabel("MarXus CSE / MESS - 1 (%)")
axes[1].set_title("all significant entries, 21 conditions (red: from the reactant, cm$^3$ s$^{-1}$)", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "cse_vs_mess.png"), dpi=200)
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

print("written:", ", ".join(sorted(os.listdir(os.path.join(HERE, "plots")))), "and the five CSV tables")
for r in kappa_rows:
    if r["T_K"] == 300.0:
        print(f"kappa 300 K {r['barrier']}: MESS {r['kappa_mess']:.4e} exact {r['kappa_marxus']:.4e} ({r['ratio']:.4f}) "
              f"MESS model {r['kappa_marxus_mess_model']:.4e} ({r['ratio_mess_model']:.5f})")
for r in shares:
    if r["p_torr"] == 760.0:
        print(f"T={r['T_K']:.0f} 760 Torr MESS-Eckart: P5 {100*r['marxus_mess_eckart_P5']:.3f}% ({r['dev_mess_eckart_P5_percent']:+.1f}%)")
for r in capture[:1] + capture[-1:]:
    print(f"capture T={r['T_K']:.0f}: E {r['dev_E_percent']:+.2f}%  EJ {r['dev_EJ_percent']:+.2f}%")
for r in shares:
    if r["p_torr"] == 760.0:
        print(f"T={r['T_K']:.0f} 760 Torr: P5 MESS {100*r['mess_P5']:.3f}%  MarXus {100*r['marxus_P5']:.3f}%  ({r['dev_P5_percent']:+.1f}%)")
