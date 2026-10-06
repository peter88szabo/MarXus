#!/usr/bin/env python3
"""Comparison of MarXus with the stored MESS results for H + C2H2 <=> C2H3.

Reads
  input/c2h3_tight.inp                       the deck (stationary points for the PES diagram)
  reference_mess_output/*.out                the stored MESS results
  marxus_output/<deck>_<method>.csv           the MarXus machine-readable tables of the four methods (run_marxus.sh,
                                               --csv; *.out: the reports; *_tables.csv: every table of a report);
                                               <deck>_olzmann_{lapack,full}.csv: SteadyStateOlzmann with the other
                                               eigen-solvers
  (writes also plots/yields.png: stabilization and prompt redissociation yields, MESS vs MarXus)
and writes
  plots/pes.png                              stationary points, Eckart barrier and absorbing barriers
  plots/falloff_P1_W1.png                    k(H + C2H2 -> C2H3) versus pressure
  plots/deviation.png                        MarXus/MESS - 1 for association and dissociation
  plots/high_pressure_limits.png             k_inf of both directions versus 1000/T
  plots/short_decks_1000K.png                1000 K, 1 atm, with and without tunneling
  plots/barrier_distance_sensitivity.png     association deviation for absorbing barriers 10, 5, 3 kT
                                             below the threshold (marxus_output/c2h3_tight_absorbing_barrier_*kT.csv)
  plots/four_methods_association.png         association: absorbing barrier, Olzmann (detailed balance), CSE vs MESS
  plots/time_evolution_1atm.png              direct time integration of a pulse at 300, 1000, 2000 K
  plots/olzmann_falloff_W1_P1.png            SteadyStateOlzmann k_uni(T, p) of the dissociation versus pressure
  plots/olzmann_deviation.png                k_uni (and lambda_1) and the detailed-balance association vs MESS
  plots/olzmann_vs_absorbing_barrier.png     association: Olzmann versus absorbing barriers 10, 5, 3 kT
  plots/olzmann_sum_rule.png                 |lambda_1 - k_uni|/k_uni (inverse iteration, LAPACK) and lambda_2/k_uni
  plots/olzmann_solvers_1000K.png            1000 K, 1 atm, with and without tunneling, the three eigen-solvers
  comparison_table.csv                       all compared numbers (absorbing barrier, final steady state)
  barrier_distance_sensitivity.csv           association for barrier distances 10, 5, 3 kT
  olzmann_comparison_table.csv               SteadyStateOlzmann: k_uni, lambda_1, sum rule, lambda_2/k_uni, association
  time_integration_decay_vs_k_uni.csv        late-time decay of the time integration against k_uni

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
    return []


def marxus_quantities(*paths):
    """k_inf(P1->W1)[T], k(P1->W1)[T,p_atm], Phi_stab[T,p_atm] (absorbing-barrier run), k_inf(W1->P1)[T] (Olzmann
    run); the two steady-state solvers are separate runs, their machine-readable files are merged."""
    b = {}
    for path in paths:
        b.update(read_marxus(path))
    k_inf_a, k_a = {}, {}
    for r in block(b, "bimolecular rate coefficients of P1 [cm3/s], intermediate"):
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
mx_kinf_a, mx_ka, mx_phi, mx_kinf_d = marxus_quantities(os.path.join(HERE, "marxus_output", "c2h3_tight_absorbing_barrier.csv"),
                                                        os.path.join(HERE, "marxus_output", "c2h3_tight_olzmann.csv"))
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
    ax.plot(p, mess, "-", marker="x", ms=10, mew=2, color=c, label=f"{t:.0f} K")
    ax.plot(p, ours, "o", color=c, ms=5)
    ax.plot([25.0], [mess_high[t][1]], "none", marker="x", ms=10, mew=2, color=c)
    if t in mx_kinf_a:
        ax.plot([25.0], [mx_kinf_a[t]], "o", color=c, ms=5)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xticks([0.1, 0.3, 1, 3, 10, 25])
ax.set_xticklabels(["0.1", "0.3", "1", "3", "10", "$\\infty$"])
ax.set_xlabel("pressure (atm)")
ax.set_ylabel("k(H + C$_2$H$_2$ $\\rightarrow$ C$_2$H$_3$) (cm$^3$ s$^{-1}$)")
ax.set_title("Association fall-off: MESS (x, lines) vs MarXus (circles)")
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
    axes[0, col].semilogy(1000 / t_all, mess, "-", marker="x", ms=10, mew=2, label="MESS")
    axes[0, col].semilogy(1000 / t_mx, ours, "o", ms=5, label="MarXus")
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
    ka_inf, ka, phi, kd_inf = marxus_quantities(os.path.join(HERE, "marxus_output", deck + "_absorbing_barrier.csv"),
                                                os.path.join(HERE, "marxus_output", deck + "_olzmann.csv"))
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
runs = (("10 kT (default)", "c2h3_tight_absorbing_barrier.csv", "tab:blue"),
        ("5 kT", "c2h3_tight_absorbing_barrier_5kT.csv", "tab:orange"),
        ("3 kT", "c2h3_tight_absorbing_barrier_3kT.csv", "tab:green"))
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

# ----------------------------------------------------------------------------------------------
# Yields of the chemically activated C2H3 (intermediate steady state, absorbing barrier 10 kT): stabilization
# and prompt redissociation, % of the formed adducts. MESS: k(P1 -> W1)/k_inf of its own tables.
# ----------------------------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
for t, c in zip(temperatures, colors):
    p = np.array(pressures)
    mess_stab = 100 * np.array([mess_p[(t, q)][1] for q in pressures]) / mess_high[t][1]
    ours = 100 * np.array([mx_phi.get((t, q), np.nan) for q in pressures])
    axes[0].plot(p, mess_stab, "-", marker="x", ms=10, mew=2, color=c, label=f"{t:.0f} K")
    axes[0].plot(p, ours, "o", color=c, ms=5)
    axes[1].plot(p, 100 - mess_stab, "-", marker="x", ms=10, mew=2, color=c, label=f"{t:.0f} K")
    axes[1].plot(p, 100 - ours, "o", color=c, ms=5)
for ax, what in zip(axes, ["stabilization of C$_2$H$_3$", "prompt redissociation to H + C$_2$H$_2$"]):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("pressure (atm)")
    ax.set_ylabel(f"{what} (% of the formed adducts)")
    ax.set_title(f"Yield of {what}: MESS (x, lines; k/k$_\\infty$) vs MarXus (circles)", fontsize=10)
    ax.legend(fontsize=7, ncol=2, title="T")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "yields.png"), dpi=200)
plt.close(fig)

# ----------------------------------------------------------------------------------------------
# The four methods: association k(P1 -> W1) (bimolecular-to-well) at 0.1 and 10 atm, and the time evolution of
# a pulse (direct time integration).
# ----------------------------------------------------------------------------------------------
def read_tables(path):
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


olz = read_marxus(os.path.join(HERE, "marxus_output", "c2h3_tight_olzmann.csv"))
detailed = {(r["T[K]"], round(r["P[Torr]"] / TORR_PER_ATM, 6)): r["k(P1->W1)"] for r in block(olz, "bimolecular rate coefficients of P1 [cm3/s] by detailed balance")}
cse_tab = read_tables(os.path.join(HERE, "marxus_output", "c2h3_tight_cse_tables.csv"))
cse_rows = next(rows for title, rows in cse_tab.items() if title.startswith("CSE: Bimolecular-to-well rate coefficients"))
cse_assoc = {(r["T[K]"], round(r["P[Torr]"] / TORR_PER_ATM, 6)): r["P1->W1"] for r in cse_rows}
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
for ax, p in zip(axes, (0.1, 10.0)):
    ts_ = [t for t in temperatures if (t, p) in mess_p]
    ax.semilogy(ts_, [mess_p[(t, p)][1] for t in ts_], "-", marker="x", ms=10, mew=2, color="k", label="MESS")
    for data, mk, c, label in [(mx_ka, "o", "C0", "steady state, absorbing barrier (10 kT)"),
                               (detailed, "D", "C2", "steady state, Olzmann (k$_{uni}$ K)"),
                               (cse_assoc, "s", "C1", "CSE (G13 eq. 28)")]:
        pts = [(t, data[(t, p)]) for t in ts_ if (t, p) in data and data[(t, p)] > 0]
        if pts:
            ax.semilogy(*zip(*pts), mk, ms=6, mfc="none" if mk in "Ds" else None, color=c, label=label)
    ax.set_xlabel("T (K)")
    ax.set_ylabel("k(H + C$_2$H$_2$ $\\rightarrow$ C$_2$H$_3$) (cm$^3$ s$^{-1}$)")
    ax.set_title(f"Bimolecular-to-well (association) at {p:g} atm: three methods vs MESS", fontsize=11)
    ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "four_methods_association.png"), dpi=200)
plt.close(fig)

evolution, key, header = {}, None, None
for line in open(os.path.join(HERE, "marxus_output", "c2h3_tight_time_integration.csv")):
    line = line.rstrip("\n")
    m = re.match(r"# time evolution: T = (\S+) K, p = (\S+) Torr", line)
    if m:
        key, header = (float(m.group(1)), float(m.group(2))), None
        evolution[key] = []
    elif line.startswith("#"):
        key = None
    elif not line.strip() or key is None:
        continue
    elif header is None:
        header = line.split(",")
    else:
        evolution[key].append(dict(zip(header, map(float, line.split(",")))))
fig, ax = plt.subplots(figsize=(8, 5))
for (t, p_torr), c in zip([(300.0, 760.0), (1000.0, 760.0), (2000.0, 760.0)], ["C0", "C1", "C3"]):
    rows = evolution.get((t, p_torr), [])
    if rows:
        ax.loglog([r["t[s]"] for r in rows], [100 * r["N(W1)"] for r in rows], "-o", ms=3, color=c, label=f"N(C$_2$H$_3$), {t:.0f} K")
        ax.loglog([r["t[s]"] for r in rows], [max(100 * r["W1->P1"], 1e-12) for r in rows], "--", color=c, label=f"H + C$_2$H$_2$ formed back, {t:.0f} K")
ax.set_ylim(1e-6, 200)
ax.set_xlabel("t (s)")
ax.set_ylabel("% of the formed adducts")
ax.set_title("Direct time integration of a pulse of chemically activated C$_2$H$_3$ (1 atm)", fontsize=11)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "time_evolution_1atm.png"), dpi=200)
plt.close(fig)

# Late-time decay of the pulse against the thermal k_uni of SteadyStateOlzmann: once the thermal eigenmode alone is
# left, N(t) ~ exp(-k_uni t). The decay rate is taken from the last two output times with 1e-8 < N < 1e-3 (the
# eigenmode alone, above the integration noise); the time integration uses no eigenvector.
k_uni = {(r["T[K]"], r["P[Torr]"]): r["k_uni[1/s]"] for r in block(olz, "thermal rate coefficients of the final steady state")}
with open(os.path.join(HERE, "time_integration_decay_vs_k_uni.csv"), "w") as f:
    f.write("T_K,p_torr,t1_s,t2_s,N_t1,N_t2,decay_rate_time_integration_s-1,k_uni_olzmann_s-1,ratio_minus_1\n")
    for (t, p_torr), rows in sorted(evolution.items()):
        late = [(r["t[s]"], r["N(W1)"]) for r in rows if 1e-8 < r["N(W1)"] < 1e-3]
        if len(late) < 2 or (t, p_torr) not in k_uni:
            continue
        (t1, n1), (t2, n2) = late[-2], late[-1]
        rate = -np.log(n2 / n1) / (t2 - t1)
        f.write(f"{t:g},{p_torr:g},{t1:.6e},{t2:.6e},{n1:.6e},{n2:.6e},{rate:.6e},{k_uni[(t, p_torr)]:.6e},"
                f"{rate / k_uni[(t, p_torr)] - 1:.3e}\n")

# ----------------------------------------------------------------------------------------------
# SteadyStateOlzmann: thermal rate coefficients from the lowest eigenpair of J (GO10 eq. 12) with the three
# eigen-solvers (marxus_output/<deck>_olzmann.csv: inverse iteration, the default; <deck>_olzmann_lapack.csv:
# LAPACK DSYEVD; <deck>_olzmann_full.csv: Householder/QL, the two 1000 K decks). k_uni(W1->P1, T, p) is the
# eigenvector average, lambda_1 is printed beside it; k(P1->W1) = k_uni k_inf,a/k_inf,d (detailed balance). Where
# the inverse iteration has no result the LAPACK result is used. Writes olzmann_comparison_table.csv and
# plots/olzmann_{falloff_W1_P1,deviation,vs_absorbing_barrier,sum_rule,solvers_1000K}.png.
# ----------------------------------------------------------------------------------------------
TOLERANCE = 1.5e-2  # default sum-rule tolerance of the example (warning threshold)
SOLVERS = (("inverse", "inverse iteration (banded Cholesky)", "tab:blue", "o"),
           ("lapack", "LAPACK DSYEVD (full)", "tab:red", "s"),
           ("full", "Householder/QL (full)", "tab:green", "^"))
SOLVER_STEM = {"inverse": "olzmann", "lapack": "olzmann_lapack", "full": "olzmann_full"}


def read_olzmann_eigen(path):
    """MarXus eigenvalue-analysis output.

    Returns (thermal, association, warned, unavailable):
      thermal[(T, p_atm)]     = row of the thermal table (k_uni, lambda_1, sum_rule_deviation, ...)
      association[(T, p_atm)] = row of the detailed-balance table
      warned                  = set of (T, p_atm) with a sum-rule warning
      unavailable[(T, p_atm)] = message of a condition without result (failed factorization)
    """
    thermal, association, warned, unavailable = {}, {}, set(), {}
    target, header = None, None
    for line in open(path):
        line = line.rstrip("\n")
        # Messages of the thermal block only: the final steady-state table before it has its own
        # "not available" lines (J N = F singular for the deep well without a sink at 300 and 500 K).
        m = re.match(r"# (not available|warning): T = (\S+) K, p = (\S+) Torr: (.*)", line)
        if m and target is thermal:
            key = (float(m.group(2)), round(float(m.group(3)) / TORR_PER_ATM, 6))
            if m.group(1) == "warning":
                warned.add(key)
            else:
                unavailable[key] = m.group(4)
            continue
        if line.startswith("# thermal rate coefficients of the final steady state"):
            target, header = thermal, None
            continue
        if line.startswith("# bimolecular rate coefficients"):
            target, header = association, None
            continue
        if line.startswith("#") or not line.strip() or target is None:
            continue
        if header is None:
            header = line.split(",")
            continue
        row = dict(zip(header, map(float, line.split(","))))
        target[(row["T[K]"], round(row["P[Torr]"] / TORR_PER_ATM, 6))] = row
    return thermal, association, warned, unavailable


def read_csv(path):
    with open(path) as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]



eig_runs = {s: read_olzmann_eigen(os.path.join(HERE, "marxus_output", f"c2h3_tight_{SOLVER_STEM[s]}.csv"))
            for s in ("inverse", "lapack")}
pcolors = plt.cm.plasma(np.linspace(0.0, 0.85, len(pressures)))
eig_previous = {(r["T_K"], r["p_atm"]): r for r in read_csv(os.path.join(HERE, "comparison_table.csv"))}
eig_sensitivity = {(r["T_K"], r["p_atm"]): r for r in read_csv(os.path.join(HERE, "barrier_distance_sensitivity.csv"))}


def eig_best(key):
    """(solver, thermal row, association row): inverse iteration, or LAPACK where the former has no result."""
    for solver in ("inverse", "lapack"):
        thermal, association, _, _ = eig_runs[solver]
        if key in thermal:
            return solver, thermal[key], association.get(key)
    return None, None, None


eig_rows = []
for t in temperatures:
    for p in pressures:
        key = (t, p)
        mess_d, mess_a = mess_p[key]
        solver, row, arow = eig_best(key)
        ku = row["k_uni[1/s]"] if row else math.nan
        lam = row["lambda_1[1/s]"] if row else math.nan
        ka = arow["k(P1->W1)"] if arow else math.nan
        prev = eig_previous.get(key, {}).get("marxus_k_P1_W1", math.nan)
        eig_rows.append({
            "T_K": t, "p_atm": p, "solver": solver or "none",
            "mess_k_W1_P1": mess_d, "marxus_k_uni": ku, "dev_k_uni_percent": 100 * (ku / mess_d - 1),
            "marxus_lambda_1": lam, "dev_lambda_1_percent": 100 * (lam / mess_d - 1),
            "sum_rule_deviation": row["sum_rule_deviation"] if row else math.nan,
            "sum_rule_warning": "yes" if key in eig_runs[solver or "inverse"][2] else "no",
            "lambda_2_over_k_uni": row["lambda_2/k_uni"] if row else math.nan,
            "mess_k_P1_W1": mess_a, "marxus_k_P1_W1_eigen": ka, "dev_P1_W1_eigen_percent": 100 * (ka / mess_a - 1),
            "marxus_k_P1_W1_absorbing_10kT": prev, "dev_P1_W1_absorbing_10kT_percent": 100 * (prev / mess_a - 1),
        })
with open(os.path.join(HERE, "olzmann_comparison_table.csv"), "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(eig_rows[0].keys()))
    writer.writeheader()
    for r in eig_rows:
        writer.writerow({k: (v if isinstance(v, str) else "%.6g" % v) for k, v in r.items()})

# ----------------------------------------------------------------------------------------------
# SteadyStateOlzmann eigen-solvers, Fall-off curves of the dissociation, k_uni
# ----------------------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.8, 5.8))
for t, c in zip(temperatures, colors):
    sel = [r for r in eig_rows if r["T_K"] == t]
    ax.plot(pressures, [mess_p[(t, q)][0] for q in pressures], "-", marker="x", ms=10, mew=2, color=c, label=f"{t:.0f} K")
    ax.plot([r["p_atm"] for r in sel if r["solver"] == "inverse"],
            [r["marxus_k_uni"] for r in sel if r["solver"] == "inverse"], "o", color=c, ms=5)
    ax.plot([r["p_atm"] for r in sel if r["solver"] == "lapack"],
            [r["marxus_k_uni"] for r in sel if r["solver"] == "lapack"], "D", color=c, ms=5)
    ax.plot([25.0], [mess_high[t][0]], "none", marker="x", ms=10, mew=2, color=c)
    k_inf = [eig_best((t, q))[1]["k_inf(W1:B1)[1/s]"] for q in pressures if eig_best((t, q))[1]]
    if k_inf:
        ax.plot([25.0], [k_inf[0]], "o", color=c, ms=5)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xticks([0.1, 0.3, 1, 3, 10, 25])
ax.set_xticklabels(["0.1", "0.3", "1", "3", "10", "$\\infty$"])
ax.set_xlabel("pressure (atm)")
ax.set_ylabel("k(C$_2$H$_3$ $\\rightarrow$ C$_2$H$_2$ + H) (s$^{-1}$)")
ax.set_title("Dissociation fall-off: MESS (x, lines) vs MarXus k$_{uni}$\n"
             "(circles: inverse iteration; diamonds: LAPACK where the Cholesky factor does not exist)", fontsize=10)
ax.legend(fontsize=8, ncol=2, title="T")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "olzmann_falloff_W1_P1.png"), dpi=200)
plt.close(fig)

# ----------------------------------------------------------------------------------------------
# SteadyStateOlzmann eigen-solvers, Deviations from MESS
# ----------------------------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9))
for p, c in zip(pressures, pcolors):
    sel = [r for r in eig_rows if r["p_atm"] == p and r["solver"] != "none"]
    axes[0].plot([r["T_K"] for r in sel], [r["dev_k_uni_percent"] for r in sel], "-o", color=c, label=f"{p:g} atm")
    lam = [(r["T_K"], r["dev_lambda_1_percent"]) for r in sel if r["sum_rule_warning"] == "yes"]
    if lam:
        axes[0].plot(*zip(*lam), "x", color=c, ms=8)
    lap = [(r["T_K"], r["dev_k_uni_percent"]) for r in sel if r["solver"] == "lapack"]
    if lap:
        axes[0].plot(*zip(*lap), "D", color=c, ms=7, mfc="none")
    e = [(r["T_K"], r["dev_P1_W1_eigen_percent"]) for r in sel if not math.isnan(r["dev_P1_W1_eigen_percent"])]
    if e:
        axes[1].plot(*zip(*e), "-o", color=c, label=f"{p:g} atm")
axes[0].plot([], [], "D", color="k", mfc="none", label="LAPACK (no Cholesky factor)")
axes[0].set_ylim(-10, 10)
axes[0].text(310, -9.3, "$\\lambda_1$ at 300 K (sum-rule warning) is off scale: 10$^{9}$-10$^{10}$ %", fontsize=7)
axes[0].set_title("dissociation C$_2$H$_3$ $\\rightarrow$ C$_2$H$_2$ + H, k$_{uni}$")
axes[1].set_title("association by detailed balance, k$_{uni}$ k$_{\\infty,a}$/k$_{\\infty,d}$")
for ax in axes:
    ax.axhspan(-5, 5, color="green", alpha=0.08)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("T (K)")
    ax.set_ylabel("MarXus / MESS - 1 (%)")
    ax.set_xlim(250, 2050)
    ax.legend(fontsize=7)
fig.suptitle("Eigenvalue analysis (no absorbing barrier) vs MESS (green band: $\\pm$5%)")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "olzmann_deviation.png"), dpi=200)
plt.close(fig)

# ----------------------------------------------------------------------------------------------
# SteadyStateOlzmann eigen-solvers, Eigenvalue route versus the absorbing barrier (association)
# ----------------------------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9), sharey=True)
for ax, p in zip(axes, (0.1, 10.0)):
    for label, column, style in (("absorbing barrier 10 kT", "dev_10kT_percent", "tab:gray"),
                                 ("absorbing barrier 5 kT", "dev_5kT_percent", "tab:orange"),
                                 ("absorbing barrier 3 kT", "dev_3kT_percent", "tab:green")):
        d = [(t, eig_sensitivity[(t, p)][column]) for t in temperatures
             if (t, p) in eig_sensitivity and not math.isnan(eig_sensitivity[(t, p)][column])]
        if d:
            ax.plot(*zip(*d), "--o", color=style, mfc="none", label=label)
    e = [(r["T_K"], r["dev_P1_W1_eigen_percent"]) for r in eig_rows
         if r["p_atm"] == p and not math.isnan(r["dev_P1_W1_eigen_percent"])]
    ax.plot(*zip(*e), "-s", color="tab:blue", lw=2, label="eigenvalue analysis (k$_{uni}$, no barrier)")
    ax.axhspan(-5, 5, color="green", alpha=0.08)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_title(f"k(H + C$_2$H$_2$ $\\rightarrow$ C$_2$H$_3$), {p:g} atm")
    ax.set_xlabel("T (K)")
    ax.legend(fontsize=8)
axes[0].set_ylabel("MarXus / MESS - 1 (%)")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "olzmann_vs_absorbing_barrier.png"), dpi=200)
plt.close(fig)

# ----------------------------------------------------------------------------------------------
# SteadyStateOlzmann eigen-solvers, Sum rule and separation of time scales
# ----------------------------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9))
ax = axes[0]
for solver, label, color, marker in SOLVERS[:2]:
    thermal, _, warned, unavailable = eig_runs[solver]
    shift = -12 if solver == "inverse" else 12
    ok = [(t + shift, max(thermal[(t, p)]["sum_rule_deviation"], 1e-16)) for t in temperatures for p in pressures
          if (t, p) in thermal and (t, p) not in warned]
    bad = [(t + shift, thermal[(t, p)]["sum_rule_deviation"]) for t in temperatures for p in pressures
           if (t, p) in thermal and (t, p) in warned]
    fail = [t + shift for t in temperatures for p in pressures if (t, p) in unavailable]
    ax.semilogy(*zip(*ok), marker, color=color, mfc="none", label=label)
    if bad:
        ax.semilogy(*zip(*bad), marker, color=color, label=f"{label}: warning")
    if fail:
        ax.semilogy(fail, [1e14] * len(fail), "x", color=color, ms=8, label=f"{label}: no Cholesky factor")
ax.axhline(TOLERANCE, color="k", ls="--", lw=1)
ax.text(1300, TOLERANCE * 3, "warning threshold 1.5e-2", fontsize=8)
ax.set_xlabel("T (K)  (all pressures)")
ax.set_ylabel("|$\\lambda_1$ - k$_{uni}$| / k$_{uni}$")
ax.set_title("sum rule (GO10 eq. 12)")
ax.legend(fontsize=7)
ax = axes[1]
for p, c in zip(pressures, pcolors):
    d = [(r["T_K"], r["lambda_2_over_k_uni"]) for r in eig_rows
         if r["p_atm"] == p and not math.isnan(r["lambda_2_over_k_uni"])]
    if d:
        ax.semilogy(*zip(*d), "-o", color=c, label=f"{p:g} atm")
ax.set_xlabel("T (K)")
ax.set_ylabel("$\\lambda_2$ / k$_{uni}$")
ax.set_title("separation of the thermal decay from relaxation")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "olzmann_sum_rule.png"), dpi=200)
plt.close(fig)

# ----------------------------------------------------------------------------------------------
# SteadyStateOlzmann eigen-solvers, 1000 K, 1 atm, with and without tunneling, all solvers
# ----------------------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8.5, 4.6))
labels, values, bar_colors = [], [], []
for deck, tag in (("c2h3_tight_short_notunneling", "no tunneling"), ("c2h3_tight_short", "Eckart")):
    mh, mp = read_mess(os.path.join(HERE, "reference_mess_output", deck + ".out"))
    for solver, label, color, _ in SOLVERS:
        thermal, association, _, _ = read_olzmann_eigen(os.path.join(HERE, "marxus_output", f"{deck}_{SOLVER_STEM[solver]}.csv"))
        key = (1000.0, 1.0)
        for name, ours, mess_value in (("k$_{uni}$", thermal[key]["k_uni[1/s]"], mp[key][0]),
                                       ("k assoc.", association[key]["k(P1->W1)"], mp[key][1])):
            labels.append(f"{name}\n{tag}\n{solver}")
            values.append(100 * (ours / mess_value - 1))
            bar_colors.append(color)
bars = ax.bar(range(len(values)), values, color=bar_colors)
for b, v in zip(bars, values):
    ax.text(b.get_x() + b.get_width() / 2, v + (0.03 if v >= 0 else -0.03), f"{v:+.2f}%", ha="center",
            va="bottom" if v >= 0 else "top", fontsize=7)
ax.axhline(0, color="k", lw=0.8)
ax.set_xticks(range(len(values)))
ax.set_xticklabels(labels, fontsize=6)
ax.set_ylabel("MarXus / MESS - 1 (%)")
ax.set_title("1000 K, 1 atm: eigenvalue analysis with the three solvers (blue inverse, red LAPACK, green QL)")
ax.margins(y=0.2)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "olzmann_solvers_1000K.png"), dpi=200)
plt.close(fig)

print("written:", ", ".join(sorted(os.listdir(os.path.join(HERE, "plots")))),
      "and comparison_table.csv, olzmann_comparison_table.csv, time_integration_decay_vs_k_uni.csv")
