#!/usr/bin/env python3
"""Olzmann's eigenvalue analysis in MarXus compared with the stored MESS results for H + C2H2 <=> C2H3.

Reads
  input/c2h3_tight.inp                         the deck (stationary points for the PES diagram)
  reference_mess_output/*.out                  the stored MESS results
  marxus_output/<deck>_<solver>.out            MarXus, --steady-state eigenvalue (run_marxus.sh);
                                               solver: inverse (default), lapack, full
  ../c2h3_mess_example/comparison_table.csv    the earlier absorbing-barrier (intermediate steady state)
  ../c2h3_mess_example/barrier_distance_sensitivity.csv   results, for comparison (read only)
and writes
  plots/pes.png                    stationary points; no absorbing barrier in this route
  plots/falloff_W1_P1.png          k_uni(T, p) versus pressure, with k_inf, against MESS
  plots/deviation.png              MarXus/MESS - 1 of the dissociation (k_uni, and lambda_1 for contrast) and of
                                   the association by detailed balance
  plots/eigen_vs_absorbing_barrier.png   association deviation: eigenvalue route versus absorbing barriers
                                   10, 5, 3 kT below the threshold
  plots/sum_rule.png               relative sum-rule deviation |lambda_1 - k_uni|/k_uni of inverse iteration and
                                   LAPACK, and lambda_2/k_uni
  plots/short_decks_1000K.png      1000 K, 1 atm, with and without tunneling, all solvers
  comparison_table.csv             all compared numbers

MarXus quantities (eigenvalue analysis, no absorbing barrier):
  k_uni(W1->P1, T, p) = sum_j k_j^th: k(E) averaged over the normalized thermal eigenvector of J
                        (Gonzalez-Garcia, Olzmann, PCCP 12, 12290 (2010), text after eq. 12); the reported
                        rate coefficient
  lambda_1            = lowest eigenvalue of J (GO10 eq. 12), equal to k_uni in exact arithmetic; printed for
                        comparison, with a warning when it deviates by more than the tolerance (1.5e-2)
  k(P1->W1, T, p)     = k_uni k_inf,assoc / k_inf,diss (detailed balance)
For every (T, p) the inverse-iteration result is used; where its Cholesky factor does not exist (300 K at
76, 2280, 7600 Torr) the LAPACK result is used, as the output message recommends; such points are marked.

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
PREVIOUS = os.path.join(HERE, "..", "c2h3_mess_example")
KCAL_PER_CM1 = 2.85914e-3
TORR_PER_ATM = 760.0
TOLERANCE = 1.5e-2  # default sum-rule tolerance of the example (warning threshold)
SOLVERS = (("inverse", "inverse iteration (banded Cholesky)", "tab:blue", "o"),
           ("lapack", "LAPACK DSYEVD (full)", "tab:red", "s"),
           ("full", "Householder/QL (full)", "tab:green", "^"))


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


def read_eigen(path):
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
        m = re.match(r"# (not available|warning): T = (\S+) K, p = (\S+) Torr: (.*)", line)
        if m:
            key = (float(m.group(2)), round(float(m.group(3)) / TORR_PER_ATM, 6))
            if m.group(1) == "warning":
                warned.add(key)
            else:
                unavailable[key] = m.group(4)
            continue
        if line.startswith("# eigenvalue analysis"):
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
runs = {s: read_eigen(os.path.join(HERE, "marxus_output", f"c2h3_tight_{s}.out")) for s in ("inverse", "lapack")}
temperatures = sorted({t for t, _ in mess_p})
pressures = sorted({p for _, p in mess_p})
os.makedirs(os.path.join(HERE, "plots"), exist_ok=True)
colors = plt.cm.viridis(np.linspace(0.0, 0.92, len(temperatures)))
pcolors = plt.cm.plasma(np.linspace(0.0, 0.85, len(pressures)))
previous = {(r["T_K"], r["p_atm"]): r for r in read_csv(os.path.join(PREVIOUS, "comparison_table.csv"))}
sensitivity = {(r["T_K"], r["p_atm"]): r for r in read_csv(os.path.join(PREVIOUS, "barrier_distance_sensitivity.csv"))}


def best(key):
    """(solver, thermal row, association row): inverse iteration, or LAPACK where the former has no result."""
    for solver in ("inverse", "lapack"):
        thermal, association, _, _ = runs[solver]
        if key in thermal:
            return solver, thermal[key], association.get(key)
    return None, None, None


rows = []
for t in temperatures:
    for p in pressures:
        key = (t, p)
        mess_d, mess_a = mess_p[key]
        solver, row, arow = best(key)
        k_uni = row["k_uni[1/s]"] if row else math.nan
        lam = row["lambda_1[1/s]"] if row else math.nan
        ka = arow["k(P1->W1)"] if arow else math.nan
        prev = previous.get(key, {}).get("marxus_k_P1_W1", math.nan)
        rows.append({
            "T_K": t, "p_atm": p, "solver": solver or "none",
            "mess_k_W1_P1": mess_d, "marxus_k_uni": k_uni, "dev_k_uni_percent": 100 * (k_uni / mess_d - 1),
            "marxus_lambda_1": lam, "dev_lambda_1_percent": 100 * (lam / mess_d - 1),
            "sum_rule_deviation": row["sum_rule_deviation"] if row else math.nan,
            "sum_rule_warning": "yes" if key in runs[solver or "inverse"][2] else "no",
            "lambda_2_over_k_uni": row["lambda_2/k_uni"] if row else math.nan,
            "mess_k_P1_W1": mess_a, "marxus_k_P1_W1_eigen": ka, "dev_P1_W1_eigen_percent": 100 * (ka / mess_a - 1),
            "marxus_k_P1_W1_absorbing_10kT": prev, "dev_P1_W1_absorbing_10kT_percent": 100 * (prev / mess_a - 1),
        })
with open(os.path.join(HERE, "comparison_table.csv"), "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    for r in rows:
        writer.writerow({k: (v if isinstance(v, str) else "%.6g" % v) for k, v in r.items()})

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
    for side in (b["left"], b["right"]):  # schematic connectors, not a reaction path
        xs, es = x_of[side], points[side][1]
        if xs < x:
            ax.plot([xs + 0.25, x - 0.25], [es, e], ls=":", color="gray")
        else:
            ax.plot([x + 0.25, xs - 0.25], [e, es], ls=":", color="gray")
ax.text(0.0, points["W1"][1] + 4.0, "no absorbing barrier:\nk$_{uni}$ from the\nthermal eigenvector of J", ha="center", fontsize=8,
        color="tab:blue")
ax.set_xlim(-0.6, 2.6)
ax.set_ylim(min(e for _, e in points.values()) - 6.0, max(b["energy"] for b in barriers) + 9.0)
ax.set_xticks([])
ax.set_ylabel("energy (kcal/mol, zero-point corrected)")
ax.set_title("H + C$_2$H$_2$ $\\rightleftharpoons$ C$_2$H$_3$: stationary points of the deck")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "pes.png"), dpi=200)
plt.close(fig)

# ----------------------------------------------------------------------------------------------
# 2. Fall-off curves of the dissociation, k_uni
# ----------------------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.8, 5.8))
for t, c in zip(temperatures, colors):
    sel = [r for r in rows if r["T_K"] == t]
    ax.plot(pressures, [mess_p[(t, q)][0] for q in pressures], "-o", color=c, mfc="none", label=f"{t:.0f} K")
    ax.plot([r["p_atm"] for r in sel if r["solver"] == "inverse"],
            [r["marxus_k_uni"] for r in sel if r["solver"] == "inverse"], "s", color=c, ms=5)
    ax.plot([r["p_atm"] for r in sel if r["solver"] == "lapack"],
            [r["marxus_k_uni"] for r in sel if r["solver"] == "lapack"], "D", color=c, ms=5)
    ax.plot([25.0], [mess_high[t][0]], "o", color=c, mfc="none")
    k_inf = [best((t, q))[1]["k_inf(W1:B1)[1/s]"] for q in pressures if best((t, q))[1]]
    if k_inf:
        ax.plot([25.0], [k_inf[0]], "s", color=c, ms=5)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xticks([0.1, 0.3, 1, 3, 10, 25])
ax.set_xticklabels(["0.1", "0.3", "1", "3", "10", "$\\infty$"])
ax.set_xlabel("pressure (atm)")
ax.set_ylabel("k(C$_2$H$_3$ $\\rightarrow$ C$_2$H$_2$ + H) (s$^{-1}$)")
ax.set_title("Dissociation fall-off: MESS (open circles) vs MarXus k$_{uni}$\n"
             "(squares: inverse iteration; diamonds: LAPACK where the Cholesky factor does not exist)", fontsize=10)
ax.legend(fontsize=8, ncol=2, title="T")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "falloff_W1_P1.png"), dpi=200)
plt.close(fig)

# ----------------------------------------------------------------------------------------------
# 3. Deviations from MESS
# ----------------------------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9))
for p, c in zip(pressures, pcolors):
    sel = [r for r in rows if r["p_atm"] == p and r["solver"] != "none"]
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
fig.savefig(os.path.join(HERE, "plots", "deviation.png"), dpi=200)
plt.close(fig)

# ----------------------------------------------------------------------------------------------
# 4. Eigenvalue route versus the absorbing barrier (association)
# ----------------------------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9), sharey=True)
for ax, p in zip(axes, (0.1, 10.0)):
    for label, column, style in (("absorbing barrier 10 kT", "dev_10kT_percent", "tab:gray"),
                                 ("absorbing barrier 5 kT", "dev_5kT_percent", "tab:orange"),
                                 ("absorbing barrier 3 kT", "dev_3kT_percent", "tab:green")):
        d = [(t, sensitivity[(t, p)][column]) for t in temperatures
             if (t, p) in sensitivity and not math.isnan(sensitivity[(t, p)][column])]
        if d:
            ax.plot(*zip(*d), "--o", color=style, mfc="none", label=label)
    e = [(r["T_K"], r["dev_P1_W1_eigen_percent"]) for r in rows
         if r["p_atm"] == p and not math.isnan(r["dev_P1_W1_eigen_percent"])]
    ax.plot(*zip(*e), "-s", color="tab:blue", lw=2, label="eigenvalue analysis (k$_{uni}$, no barrier)")
    ax.axhspan(-5, 5, color="green", alpha=0.08)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_title(f"k(H + C$_2$H$_2$ $\\rightarrow$ C$_2$H$_3$), {p:g} atm")
    ax.set_xlabel("T (K)")
    ax.legend(fontsize=8)
axes[0].set_ylabel("MarXus / MESS - 1 (%)")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "eigen_vs_absorbing_barrier.png"), dpi=200)
plt.close(fig)

# ----------------------------------------------------------------------------------------------
# 5. Sum rule and separation of time scales
# ----------------------------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9))
ax = axes[0]
for solver, label, color, marker in SOLVERS[:2]:
    thermal, _, warned, unavailable = runs[solver]
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
    d = [(r["T_K"], r["lambda_2_over_k_uni"]) for r in rows
         if r["p_atm"] == p and not math.isnan(r["lambda_2_over_k_uni"])]
    if d:
        ax.semilogy(*zip(*d), "-o", color=c, label=f"{p:g} atm")
ax.set_xlabel("T (K)")
ax.set_ylabel("$\\lambda_2$ / k$_{uni}$")
ax.set_title("separation of the thermal decay from relaxation")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "plots", "sum_rule.png"), dpi=200)
plt.close(fig)

# ----------------------------------------------------------------------------------------------
# 6. 1000 K, 1 atm, with and without tunneling, all solvers
# ----------------------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8.5, 4.6))
labels, values, bar_colors = [], [], []
for deck, tag in (("c2h3_tight_short_notunneling", "no tunneling"), ("c2h3_tight_short", "Eckart")):
    mh, mp = read_mess(os.path.join(HERE, "reference_mess_output", deck + ".out"))
    for solver, label, color, _ in SOLVERS:
        thermal, association, _, _ = read_eigen(os.path.join(HERE, "marxus_output", f"{deck}_{solver}.out"))
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
fig.savefig(os.path.join(HERE, "plots", "short_decks_1000K.png"), dpi=200)
plt.close(fig)

print("written:", ", ".join(sorted(os.listdir(os.path.join(HERE, "plots")))), "and comparison_table.csv")
