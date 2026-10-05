#!/usr/bin/env python3
"""Comparison of MarXus with the stored MESS results for H + C2H2 <=> C2H3.

Reads
  input/c2h3_tight.inp                       the deck (stationary points for the PES diagram)
  reference_mess_output/*.out                the stored MESS results
  marxus_output/*.out                        the MarXus results (run_marxus.sh)
and writes
  plots/pes.png                              stationary points, Eckart barrier and absorbing barriers
  plots/falloff_P1_W1.png                    k(H + C2H2 -> C2H3) versus pressure
  plots/deviation.png                        MarXus/MESS - 1 for association and dissociation
  plots/high_pressure_limits.png             k_inf of both directions versus 1000/T
  plots/short_decks_1000K.png                1000 K, 1 atm, with and without tunneling
  plots/barrier_distance_sensitivity.png     association deviation for absorbing barriers 10, 5, 3 kT
                                             below the threshold (marxus_output/c2h3_tight_barrier_*kT.out)
  comparison_table.csv                       all compared numbers

MarXus quantities:
  k_inf(P1->W1)      high-pressure association rate coefficient (column k_inf)
  k(P1->W1, T, p)    = k_inf Phi_stab, intermediate steady state (column k(P1->W1))
  k_inf(W1->P1)      final steady state with thermal formation through the only channel, which is
                     equilibrium: k^ca = canonical high-pressure dissociation rate coefficient
  k(W1->P1, T, p)    = k_inf(W1->P1) Phi_stab (detailed balance: k_d(p)/k_a(p) = k_inf,d/k_inf,a)

Run with the science environment:  source ~/.venvs/science/bin/activate && python3 plot_comparison.py
"""
import csv
import math
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
KCAL_PER_CM1 = 2.85914e-3
KB_CM1 = 0.69503476  # cm-1/K
TORR_PER_ATM = 760.0


# ----------------------------------------------------------------------------------------------
# Readers
# ----------------------------------------------------------------------------------------------
def read_mess(path):
    """Species-species tables: high-pressure {T: (k_W1->P1, k_P1->W1)} and {(T, p_atm): (...)}."""
    lines = open(path).read().split("\n")
    high, pressure = {}, {}
    i = 0
    while i < len(lines):
        m_tp = re.match(r"\s*Temperature = (\S+) K\s+Pressure = (\S+) atm", lines[i])
        m_t = re.match(r"\s*Temperature = (\S+) K\s*$", lines[i])
        if m_tp or (m_t and i + 2 < len(lines) and "High Pressure Rate Coefficients" in lines[i + 2]):
            j = i + 1
            while not lines[j].strip().startswith("From\\To"):
                j += 1
            header = lines[j].split()[1:]
            rows = {}
            for k in range(1, len(header) + 1):
                fields = lines[j + k].split()
                rows[fields[0]] = dict(zip(header, fields[1:]))
            value = lambda a, b: float(rows[a][b]) if rows[a][b] != "***" else math.nan  # noqa: E731
            pair = (value("W1", "P1"), value("P1", "W1"))
            if m_tp:
                pressure[(float(m_tp.group(1)), float(m_tp.group(2)))] = pair
            else:
                high[float(m_t.group(1))] = pair
            i = j + len(header)
        i += 1
    return high, pressure


def read_marxus(path):
    """Blocks of the MarXus example output: {block title: list of row dicts}."""
    blocks, title, header = {}, None, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("# intermediate steady state") or line.startswith("# final steady state") \
                or line.startswith("# bimolecular rate coefficients"):
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
    return []


def marxus_quantities(path):
    """k_inf(P1->W1)[T], k(P1->W1)[T,p_atm], Phi_stab[T,p_atm], k_inf(W1->P1)[T] (final steady state)."""
    b = read_marxus(path)
    k_inf_a, k_a = {}, {}
    for r in block(b, "bimolecular rate coefficients"):
        t, p = r["T[K]"], round(r["P[Torr]"] / TORR_PER_ATM, 6)
        k_inf_a[t] = r["k_inf"]
        k_a[(t, p)] = r["k(P1->W1)"]
    phi = {(r["T[K]"], round(r["P[Torr]"] / TORR_PER_ATM, 6)): r["Phi_stab(W1)"] for r in block(b, "intermediate steady state")}
    k_inf_d = {}
    for r in block(b, "final steady state"):
        k_inf_d[r["T[K]"]] = r["k_ca(W1:B1)[1/s]"]  # pressure independent (equilibrium)
    return k_inf_a, k_a, phi, k_inf_d


def read_deck(path):
    """Stationary points of the deck: wells, bimolecular asymptotes, barriers (energies in kcal/mol)."""
    text = [l.split("#")[0].split("!")[0].rstrip() for l in open(path)]
    points, barriers, current, kind = {}, [], None, None
    for line in text:
        tok = line.split()
        if not tok:
            continue
        if tok[0] in ("Well", "Bimolecular"):
            current, kind = tok[1], tok[0]
        elif tok[0] == "Barrier":
            current, kind = tok[1], "Barrier"
            barriers.append({"name": tok[1], "left": tok[2], "right": tok[3]})
        elif tok[0].startswith("ZeroEnergy") or tok[0].startswith("GroundEnergy"):
            value = float(tok[1])
            if "1/cm" in tok[0]:
                value *= KCAL_PER_CM1
            if kind == "Well" and tok[0].startswith("ZeroEnergy"):
                points[current] = ("well", value)
            elif kind == "Bimolecular" and tok[0].startswith("GroundEnergy"):
                points[current] = ("bimolecular", value)
            elif kind == "Barrier" and tok[0].startswith("ZeroEnergy"):
                barriers[-1]["energy"] = value
        elif tok[0].startswith("ImaginaryFrequency") and kind == "Barrier":
            barriers[-1]["imaginary"] = float(tok[1])
    return points, barriers


# ----------------------------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------------------------
mess_high, mess_p = read_mess(os.path.join(HERE, "reference_mess_output", "c2h3_tight.out"))
mx_kinf_a, mx_ka, mx_phi, mx_kinf_d = marxus_quantities(os.path.join(HERE, "marxus_output", "c2h3_tight.out"))
temperatures = sorted({t for t, _ in mess_p})
pressures = sorted({p for _, p in mess_p})
os.makedirs(os.path.join(HERE, "plots"), exist_ok=True)
colors = plt.cm.viridis(np.linspace(0.0, 0.92, len(temperatures)))

# Table of all compared numbers.
rows = []
for t in temperatures:
    for p in pressures:
        mess_d, mess_a = mess_p[(t, p)]
        ka = mx_ka.get((t, p), math.nan)
        kd = mx_kinf_d[t] * mx_phi[(t, p)] if (t in mx_kinf_d and (t, p) in mx_phi) else math.nan
        rows.append({"T_K": t, "p_atm": p, "mess_k_P1_W1": mess_a, "marxus_k_P1_W1": ka,
                     "dev_P1_W1_percent": 100 * (ka / mess_a - 1), "mess_k_W1_P1": mess_d,
                     "marxus_k_W1_P1": kd, "dev_W1_P1_percent": 100 * (kd / mess_d - 1)})
with open(os.path.join(HERE, "comparison_table.csv"), "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    for r in rows:
        writer.writerow({k: ("%.6g" % v) for k, v in r.items()})

# ----------------------------------------------------------------------------------------------
# 1. Potential-energy surface
# ----------------------------------------------------------------------------------------------
points, barriers = read_deck(os.path.join(HERE, "input", "c2h3_tight.inp"))
fig, ax = plt.subplots(figsize=(7.5, 5.5))
x_of = {"W1": 0.0, "B1": 1.0, "P1": 2.0}
labels = {"W1": "C$_2$H$_3$ (W1)", "P1": "C$_2$H$_2$ + H (P1)"}
for name, (kind, e) in points.items():
    x = x_of[name]
    ax.hlines(e, x - 0.25, x + 0.25, color="k", lw=3)
    ax.text(x, e - 2.6, f"{labels.get(name, name)}\n{e:.2f} kcal/mol", ha="center", va="top", fontsize=9)
for b in barriers:
    x, e = x_of[b["name"]], b["energy"]
    ax.hlines(e, x - 0.25, x + 0.25, color="firebrick", lw=3)
    text = f"TS {b['name']}\n{e:.2f} kcal/mol"
    if "imaginary" in b:
        text += f"\nEckart tunneling, {b['imaginary']:.0f}i cm$^{{-1}}$"
    ax.text(x, e + 1.0, text, ha="center", va="bottom", fontsize=9, color="firebrick")
    # Connectors between the stationary points (schematic, not a reaction path).
    for side in (b["left"], b["right"]):
        xs, es = x_of[side], points[side][1]
        if xs < x:
            ax.plot([xs + 0.25, x - 0.25], [es, e], ls=":", color="gray")
        else:
            ax.plot([x + 0.25, xs - 0.25], [e, es], ls=":", color="gray")
# Absorbing barriers of the intermediate steady state: 10 kT below the classical threshold of W1.
well_e = points["W1"][1]
threshold = barriers[0]["energy"]
for t, c in ((1000, "tab:blue"), (1500, "tab:orange"), (2000, "tab:red")):
    eb = threshold - 10 * KB_CM1 * t * KCAL_PER_CM1
    ax.hlines(eb, -0.45, 0.45, color=c, lw=1.2, ls="--")
    note = " (below the well bottom)" if eb < well_e else ""
    ax.text(0.47, eb, f"absorbing barrier, {t} K{note}", va="center", fontsize=8, color=c)
ax.set_xlim(-0.6, 2.6)
ax.set_ylim(min(e for _, e in points.values()) - 6.0, max(b["energy"] for b in barriers) + 9.0)
ax.set_xticks([])
ax.set_ylabel("energy (kcal/mol, zero-point corrected)")
ax.set_title("H + C$_2$H$_2$ $\\rightleftharpoons$ C$_2$H$_3$: stationary points of the deck")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "pes.png"), dpi=200)
plt.close(fig)

# ----------------------------------------------------------------------------------------------
# 2. Fall-off curves of the association
# ----------------------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.5, 5.5))
for t, c in zip(temperatures, colors):
    p = np.array(pressures)
    mess = np.array([mess_p[(t, q)][1] for q in pressures])
    ours = np.array([mx_ka.get((t, q), np.nan) for q in pressures])
    ax.plot(p, mess, "-o", color=c, mfc="none", label=f"{t:.0f} K")
    ax.plot(p, ours, "s", color=c, ms=5)
    ax.plot([25.0], [mess_high[t][1]], "o", color=c, mfc="none")
    if t in mx_kinf_a:
        ax.plot([25.0], [mx_kinf_a[t]], "s", color=c, ms=5)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xticks([0.1, 0.3, 1, 3, 10, 25])
ax.set_xticklabels(["0.1", "0.3", "1", "3", "10", "$\\infty$"])
ax.set_xlabel("pressure (atm)")
ax.set_ylabel("k(H + C$_2$H$_2$ $\\rightarrow$ C$_2$H$_3$) (cm$^3$ s$^{-1}$)")
ax.set_title("Association fall-off: MESS (open circles, lines) vs MarXus (squares)")
ax.legend(fontsize=8, ncol=2, title="T")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "falloff_P1_W1.png"), dpi=200)
plt.close(fig)

# ----------------------------------------------------------------------------------------------
# 3. Deviations
# ----------------------------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharey=False)
pc = plt.cm.plasma(np.linspace(0.0, 0.85, len(pressures)))
for ax, key, title in ((axes[0], "dev_P1_W1_percent", "association H + C$_2$H$_2$ $\\rightarrow$ C$_2$H$_3$"),
                       (axes[1], "dev_W1_P1_percent", "dissociation C$_2$H$_3$ $\\rightarrow$ C$_2$H$_2$ + H")):
    for p, c in zip(pressures, pc):
        d = [(r["T_K"], r[key]) for r in rows if r["p_atm"] == p and not math.isnan(r[key])]
        if d:
            ax.plot(*zip(*d), "-o", color=c, label=f"{p:g} atm")
    ax.axhspan(-5, 5, color="green", alpha=0.08)
    ax.axhline(0, color="k", lw=0.8)
    ax.axvspan(1375, 2100, color="gray", alpha=0.12)
    ax.text(1400, ax.get_ylim()[0] * 0.9 if ax.get_ylim()[0] < 0 else -5, "absorbing barrier inside the\nthermal distribution of W1",
            fontsize=8, va="bottom")
    ax.set_xlabel("T (K)")
    ax.set_ylabel("MarXus / MESS - 1 (%)")
    ax.set_title(title)
    ax.set_xlim(250, 2050)
axes[0].legend(fontsize=8)
fig.suptitle("Deviation of MarXus from MESS (green band: $\\pm$5%)")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "deviation.png"), dpi=200)
plt.close(fig)

# ----------------------------------------------------------------------------------------------
# 4. High-pressure limits
# ----------------------------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1.3]})
t_all = np.array(temperatures)
for col, (index, mx, label) in enumerate(((1, mx_kinf_a, "k$_\\infty$(H + C$_2$H$_2$ $\\rightarrow$ C$_2$H$_3$) (cm$^3$ s$^{-1}$)"),
                                          (0, mx_kinf_d, "k$_\\infty$(C$_2$H$_3$ $\\rightarrow$ C$_2$H$_2$ + H) (s$^{-1}$)"))):
    mess = np.array([mess_high[t][index] for t in temperatures])
    t_mx = np.array([t for t in temperatures if t in mx])
    ours = np.array([mx[t] for t in t_mx])
    axes[0, col].semilogy(1000 / t_all, mess, "-o", mfc="none", label="MESS")
    axes[0, col].semilogy(1000 / t_mx, ours, "s", label="MarXus")
    axes[0, col].set_ylabel(label)
    axes[0, col].legend()
    ratio = ours / np.array([mess_high[t][index] for t in t_mx])
    axes[1, col].plot(1000 / t_mx, 100 * (ratio - 1), "s-", color="tab:orange")
    axes[1, col].axhline(0, color="k", lw=0.8)
    axes[1, col].set_ylabel("MarXus/MESS - 1 (%)")
    axes[1, col].set_xlabel("1000/T (K$^{-1}$)")
axes[0, 1].set_title("dissociation (MarXus: final steady state, where defined)")
axes[0, 0].set_title("association (MarXus: canonical k$_\\infty$ of the entrance)")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "high_pressure_limits.png"), dpi=200)
plt.close(fig)

# ----------------------------------------------------------------------------------------------
# 5. 1000 K, 1 atm, with and without tunneling
# ----------------------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 4.5))
labels, values = [], []
for deck, tag in (("c2h3_tight_short_notunneling", "no tunneling"), ("c2h3_tight_short", "Eckart")):
    mh, mp = read_mess(os.path.join(HERE, "reference_mess_output", deck + ".out"))
    ka_inf, ka, phi, kd_inf = marxus_quantities(os.path.join(HERE, "marxus_output", deck + ".out"))
    for name, mess_value, ours in (
            ("k$_\\infty$ assoc.", mh[1000.0][1], ka_inf[1000.0]),
            ("k(1 atm) assoc.", mp[(1000.0, 1.0)][1], ka[(1000.0, 1.0)]),
            ("k$_\\infty$ dissoc.", mh[1000.0][0], kd_inf[1000.0]),
            ("k(1 atm) dissoc.", mp[(1000.0, 1.0)][0], kd_inf[1000.0] * phi[(1000.0, 1.0)])):
        labels.append(f"{name}\n{tag}")
        values.append(100 * (ours / mess_value - 1))
bars = ax.bar(range(len(values)), values, color=["tab:blue"] * 4 + ["tab:red"] * 4)
for b, v in zip(bars, values):
    ax.text(b.get_x() + b.get_width() / 2, v + (0.08 if v >= 0 else -0.08), f"{v:+.1f}%", ha="center",
            va="bottom" if v >= 0 else "top", fontsize=8)
ax.axhline(0, color="k", lw=0.8)
ax.set_xticks(range(len(values)))
ax.set_xticklabels(labels, fontsize=7)
ax.set_ylabel("MarXus / MESS - 1 (%)")
ax.set_title("1000 K, 1 atm: blue without tunneling, red with Eckart tunneling")
ax.margins(y=0.15)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "short_decks_1000K.png"), dpi=200)
plt.close(fig)

# ----------------------------------------------------------------------------------------------
# 6. Sensitivity to the absorbing-barrier distance (user choice for shallow wells)
# ----------------------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
runs = (("10 kT (default)", "c2h3_tight.out", "tab:blue"),
        ("5 kT", "c2h3_tight_barrier_5kT.out", "tab:orange"),
        ("3 kT", "c2h3_tight_barrier_3kT.out", "tab:green"))
sensitivity = {}
for label, name, color in runs:
    _, ka, _, _ = marxus_quantities(os.path.join(HERE, "marxus_output", name))
    sensitivity[label] = ka
    for p, style in ((0.1, "--"), (10.0, "-")):
        d = [(t, 100 * (ka[(t, p)] / mess_p[(t, p)][1] - 1)) for t in temperatures if (t, p) in ka]
        if d:
            ax.plot(*zip(*d), style, marker="o", color=color, label=f"{label}, {p:g} atm")
ax.axhspan(-5, 5, color="green", alpha=0.08)
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("T (K)")
ax.set_ylabel("k(P1$\\rightarrow$W1): MarXus / MESS - 1 (%)")
ax.set_title("Absorbing-barrier distance below the threshold (intermediate steady state)")
ax.legend(fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "barrier_distance_sensitivity.png"), dpi=200)
plt.close(fig)
with open(os.path.join(HERE, "barrier_distance_sensitivity.csv"), "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["T_K", "p_atm", "mess_k_P1_W1"] + [f"marxus_{l.split()[0]}kT" for l, _, _ in runs]
                    + [f"dev_{l.split()[0]}kT_percent" for l, _, _ in runs])
    for t in temperatures:
        for p in pressures:
            m = mess_p[(t, p)][1]
            ks = [sensitivity[l].get((t, p), math.nan) for l, _, _ in runs]
            writer.writerow([f"{t:g}", f"{p:g}", f"{m:.6g}"] + [f"{k:.6g}" for k in ks]
                            + [f"{100 * (k / m - 1):.2f}" for k in ks])

print("written:", ", ".join(sorted(os.listdir(os.path.join(HERE, "plots")))), "and comparison_table.csv")
