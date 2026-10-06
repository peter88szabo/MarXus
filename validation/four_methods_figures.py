#!/usr/bin/env python3
"""The four MarXus methods on both validation systems: against MESS (validation) and against each other
(internal). This script only adds figures; the figures of the other scripts are not touched.

Methods (outputs of the run scripts, the same deck and settings):
  SteadyStateOlzmann            final steady state with its thermal eigenpair (k_uni, thermal fates)
  SteadyStateAbsorbingBarrier   intermediate steady state (prompt yields, stabilization)
  CSE                           phenomenological rate coefficients (Georgievskii et al. 2013)
  TimeIntegration               pulse of chemically activated adducts, populations and yields in time

Writes into each system directory:
  plots/mess_four_methods_*.png        each method against MESS, wherever the method gives the MESS quantity
  plots/internal_four_methods_*.png    the four methods against each other (no MESS)
  four_methods_figures.csv             every plotted comparison: quantity, method, condition, value, reference

Quantities that a method does not give by itself:
  - TimeIntegration, one well: the association is k_inf A, with A the amplitude of the slowest mode of the pulse,
    N(t) = A exp(-r t) at late times, extrapolated to t = 0. For one well A = k(R -> W)/k_inf of CSE (G13 eq. 28)
    exactly: both are (sum_E f1(E)) (sum_E f1(E) k_R(E)) / sum_E k_R(E) f0(E), f1 the slowest eigenvector of J.
    The dissociation is r, where the pulse decays below 1e-3 inside the output window (750-2000 K for C2H3).
  - SteadyStateAbsorbingBarrier: the dissociation of one well is k_inf,d Phi_stab (detailed balance of the pair,
    k_d(p)/k_a(p) = k_inf,d/k_inf,a); its long-time yields are prompt + stabilization x thermal fate, with the
    thermal fates of SteadyStateOlzmann (labelled so).
  - k_inf (the capture rate coefficient) is a property of the deck, the same in every run; the TimeIntegration
    association uses the one of the CSE run.

Run with the science environment:  source ~/.venvs/science/bin/activate && python3 four_methods_figures.py
"""
import csv
import math
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from method_comparison import (FLOOR_RATIO, read_evolution, read_mess, read_mess_high_pressure, read_tables,  # noqa: E402
                               read_thermal, table)

HERE = os.path.dirname(os.path.abspath(__file__))
TORR_PER_ATM = 760.0
NAN = math.nan

STYLE = {
    "olz": dict(name="SteadyStateOlzmann", color="tab:green", marker="D", ms=7),
    "bar": dict(name="SteadyStateAbsorbingBarrier", color="tab:blue", marker="o", ms=5),
    "cse": dict(name="CSE", color="tab:orange", marker="s", ms=6),
    "ti": dict(name="TimeIntegration", color="tab:red", marker="^", ms=8),
}
ORDER = ("olz", "bar", "cse", "ti")


def finite(x):
    return x is not None and isinstance(x, float) and math.isfinite(x)


# ----------------------------------------------------------------------------------------------
# Plot helpers
# ----------------------------------------------------------------------------------------------
def mark(ax, pts, m, label=None, log=False, **kw):
    """Small markers of method m (MarXus); `label` None: the method name."""
    pts = [(x, y) for x, y in pts if finite(float(y)) and (not log or y > 0)]
    if not pts:
        return
    s = STYLE[m]
    ax.plot(*zip(*pts), ls="none", marker=s["marker"], ms=s["ms"], mec=s["color"],
            mfc=s["color"] if m == "bar" else "none", label=s["name"] if label is None else label, **kw)


def mess_line(ax, pts, label="MESS", color="k", log=False):
    """MESS: large x markers on a line."""
    pts = [(x, y) for x, y in pts if finite(float(y)) and (not log or y > 0)]
    if pts:
        ax.plot(*zip(*pts), "-", marker="x", ms=10, mew=2, color=color, label=label)


def at_p(d, p):
    return sorted((t, v) for (t, q), v in d.items() if q == p and finite(v))


def at_t(d, t):
    return sorted((q, v) for (s, q), v in d.items() if s == t and finite(v))


def deviation(a, b):
    """{key: 100 (a/b - 1)} where both are finite."""
    return {k: 100 * (a[k] / b[k] - 1) for k in a if k in b and finite(a[k]) and finite(b[k]) and b[k] != 0}


def abs_deviation(a, b):
    return {k: max(abs(a[k] / b[k] - 1), 1e-17) for k in a if k in b and finite(a[k]) and finite(b[k]) and b[k] != 0}


def band(ax, dev, color, label, marker="o", ls="-"):
    """Line through the mean over the pressures with the min-max band (shaded)."""
    ts = sorted({t for t, _ in dev})
    if not ts:
        return
    lo = [min(v for (s, _), v in dev.items() if s == t) for t in ts]
    hi = [max(v for (s, _), v in dev.items() if s == t) for t in ts]
    ax.plot(ts, [(a + b) / 2 for a, b in zip(lo, hi)], ls=ls, marker=marker, ms=5, mfc="none", color=color, label=label)
    ax.fill_between(ts, lo, hi, color=color, alpha=0.13)


def percent_axes(ax, xlabel="T (K) (band: range over the pressures)"):
    ax.axhline(0, color="k", lw=0.8)
    ax.axhspan(-5, 5, color="green", alpha=0.06)
    ax.set_xlabel(xlabel)


def log_dev_axes(ax):
    ax.set_yscale("log")
    ax.axhspan(1e-17, 1e-6, color="green", alpha=0.07)
    ax.set_xlabel("T (K)")
    ax.set_ylabel("|a/b - 1| per condition (green: below the printed precision)")
    ax.text(0.01, 2.5e-17, "1e-17: equal in all printed digits", transform=ax.get_yaxis_transform(), fontsize=7,
            color="0.3")


def save(fig, d, name):
    fig.tight_layout()
    fig.savefig(os.path.join(d, "plots", name), dpi=200)
    plt.close(fig)


class Record:
    """Every plotted comparison, written to four_methods_figures.csv."""

    def __init__(self):
        self.rows = []

    def add(self, figure, quantity, method, values, reference, reference_label):
        for key in sorted(values):
            ref = reference.get(key, NAN)
            v = values[key]
            self.rows.append({"figure": figure, "quantity": quantity, "method": method, "T_K": key[0],
                              "p_torr": key[1], "value": v, "reference": reference_label, "reference_value": ref,
                              "deviation_percent": 100 * (v / ref - 1) if finite(v) and finite(ref) and ref > 0 else NAN})

    def write(self, d):
        path = os.path.join(d, "four_methods_figures.csv")
        fields = ["figure", "quantity", "method", "T_K", "p_torr", "value", "reference", "reference_value",
                  "deviation_percent"]
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for r in self.rows:
                w.writerow({k: ("%.8g" % r[k]) if isinstance(r[k], float) else r[k] for k in fields})
        return path

    def summary(self):
        """{(figure, quantity, method, reference): (min, max) of the deviation in %}."""
        out = {}
        for r in self.rows:
            if not finite(r["deviation_percent"]):
                continue
            k = (r["figure"], r["quantity"], r["method"], r["reference"])
            lo, hi = out.get(k, (math.inf, -math.inf))
            out[k] = (min(lo, r["deviation_percent"]), max(hi, r["deviation_percent"]))
        return out


def slow_mode(rows, wells):
    """Amplitude A and rate r of the slowest mode of the pulse, N(t) = A exp(-r t) for the total well population,
    from the last two output times with N > 1e-8; `decayed`: N fell below 1e-3 there (r is resolved)."""
    pts = [(r["t[s]"], sum(r[f"N({w})"] for w in wells)) for r in rows]
    late = [(t, n) for t, n in pts if n > 1e-8]
    if len(late) < 2:
        return NAN, NAN, False
    (t1, n1), (t2, n2) = late[-2], late[-1]
    rate = -math.log(n2 / n1) / (t2 - t1)
    return n2 * math.exp(rate * t2), rate, n2 < 1e-3


def absorbing_chain(rates, wells, ends):
    """Thermal fate of each well from well rate coefficients {w: {x: k}}: B = (I - Q)^-1 A (fractions)."""
    n = len(wells)
    q, a = np.zeros((n, n)), np.zeros((n, len(ends)))
    for i, w in enumerate(wells):
        out = {x: rates[w][x] for x in wells + ends if x != w and finite(rates[w].get(x, NAN))}
        loss = sum(out.values())
        for j, v in enumerate(wells):
            if v != w:
                q[i, j] = out.get(v, 0.0) / loss
        for e, x in enumerate(ends):
            a[i, e] = out.get(x, 0.0) / loss
    b = np.linalg.solve(np.eye(n) - q, a)
    return {w: {x: b[i, e] for e, x in enumerate(ends)} for i, w in enumerate(wells)}


# ----------------------------------------------------------------------------------------------
# H + C2H2 <=> C2H3 (one well)
# ----------------------------------------------------------------------------------------------
def c2h3():
    d = os.path.join(HERE, "c2h3_mess_example")
    o = os.path.join(d, "marxus_output")
    stem = "c2h3_tight"
    tabs = {m: read_tables(os.path.join(o, f"{stem}_{m}_tables.csv"))
            for m in ("olzmann", "absorbing_barrier", "cse", "time_integration")}
    evolution = read_evolution(os.path.join(o, f"{stem}_time_integration.csv"))
    thermal = read_thermal(os.path.join(o, f"{stem}_olzmann.csv"))
    mess = read_mess(os.path.join(d, "reference_mess_output", f"{stem}.out"))
    mess_hp = read_mess_high_pressure(os.path.join(d, "reference_mess_output", f"{stem}.out"))
    rec = Record()

    bar_rw = table(tabs["absorbing_barrier"], "intermediate steady state: Bimolecular-to-well rate coefficients")
    bar_cap = table(tabs["absorbing_barrier"], "intermediate steady state: Capture, return and net reaction")
    bar_kca = table(tabs["absorbing_barrier"], "intermediate steady state: Chemical-activation rate coefficients")
    olz_assoc = table(tabs["olzmann"], "thermal rate coefficients: Association by detailed balance")
    olz_th = table(tabs["olzmann"], "thermal rate coefficients: Thermal rate coefficients")
    olz_hp = table(tabs["olzmann"], "thermal rate coefficients: High-pressure rate coefficients k_inf")
    olz_kca = table(tabs["olzmann"], "final steady state: Chemical-activation rate coefficients")
    cse_rw = table(tabs["cse"], "CSE: Bimolecular-to-well rate coefficients")
    cse_w = table(tabs["cse"], "CSE: Rate coefficients from W1")
    cse_cap = table(tabs["cse"], "CSE: Capture, return and net reaction")
    keys = sorted(cse_w)
    temperatures = sorted({t for t, _ in keys})
    pressures = sorted({p for _, p in keys})
    below_floor = {k for k, r in thermal.items() if r["precision_floor[1/s]"] > FLOOR_RATIO * r["k_uni[1/s]"]}

    V = {m: {"assoc": {}, "dissoc": {}, "stab": {}} for m in ORDER}
    amplitude = {}
    for key in keys:
        k_inf = cse_cap[key]["capture"]
        if finite(bar_cap.get(key, {}).get("capture", NAN)) and abs(bar_cap[key]["capture"] / k_inf - 1) > 1e-5:
            print(f"warning: C2H3 {key}: k_inf of the absorbing-barrier and CSE runs differ")
        V["bar"]["assoc"][key] = bar_rw[key]["P1->W1"]
        V["bar"]["stab"][key] = bar_rw[key]["P1->W1"] / k_inf
        V["bar"]["dissoc"][key] = olz_hp[key]["W1->P1"] * V["bar"]["stab"][key]
        V["olz"]["assoc"][key] = olz_assoc[key]["P1->W1"]
        V["olz"]["stab"][key] = olz_assoc[key]["P1->W1"] / k_inf
        V["olz"]["dissoc"][key] = olz_th[key]["k_uni"]
        V["cse"]["assoc"][key] = cse_rw[key]["P1->W1"]
        V["cse"]["stab"][key] = cse_rw[key]["P1->W1"] / k_inf
        V["cse"]["dissoc"][key] = cse_w[key]["W1->P1"]
        if key in evolution:
            a, r, decayed = slow_mode(evolution[key][1], ["W1"])
            amplitude[key] = (a, r)
            V["ti"]["assoc"][key] = k_inf * a
            V["ti"]["stab"][key] = a
            if decayed:
                V["ti"]["dissoc"][key] = r
    M = {"assoc": {}, "dissoc": {}, "stab": {}}
    for key in keys:
        if key in mess:
            M["assoc"][key] = mess[key]["P1"]["W1"]
            M["dissoc"][key] = mess[key]["W1"]["P1"]
            M["stab"][key] = mess[key]["P1"]["W1"] / mess_hp[key[0]]["P1"]["W1"]
    for m in ORDER:
        V[m]["return"] = {k: 1 - v for k, v in V[m]["stab"].items()}
    M["return"] = {k: 1 - v for k, v in M["stab"].items()}

    names = {"assoc": "association k(H + C$_2$H$_2$ $\\rightarrow$ C$_2$H$_3$) (cm$^3$ s$^{-1}$)",
             "dissoc": "dissociation k(C$_2$H$_3$ $\\rightarrow$ C$_2$H$_2$ + H) (s$^{-1}$)"}

    # ---- A1: rates against MESS: absolute at 0.1 and 10 atm, and the deviation over all pressures ----
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    for row, q in enumerate(("assoc", "dissoc")):
        for col, p in enumerate((76.0, 7600.0)):
            ax = axes[row][col]
            mess_line(ax, at_p(M[q], p), log=True)
            for m in ORDER:
                mark(ax, at_p(V[m][q], p), m, log=True)
            ax.set_yscale("log")
            ax.set_xlabel("T (K)")
            ax.set_ylabel(names[q])
            ax.set_title(f"{p / TORR_PER_ATM:g} atm: the four methods (small markers) and MESS (x)", fontsize=10)
            ax.legend(fontsize=7)
        ax = axes[row][2]
        for m in ORDER:
            dev = deviation(V[m][q], M[q])
            band(ax, dev, STYLE[m]["color"], STYLE[m]["name"], STYLE[m]["marker"])
            rec.add("mess_four_methods_rates", q, STYLE[m]["name"], V[m][q], M[q], "MESS")
        percent_axes(ax)
        ax.set_ylabel("MarXus / MESS - 1 (%)")
        ax.set_title(f"{'association' if q == 'assoc' else 'dissociation'}: deviation from MESS, all pressures",
                     fontsize=10)
        ax.legend(fontsize=7)
    fig.suptitle("H + C$_2$H$_2$ $\\rightleftharpoons$ C$_2$H$_3$: rate coefficients of the four methods against MESS. "
                 "TimeIntegration: association k$_\\infty$ A (slowest-mode amplitude), dissociation = late decay "
                 "rate (750-2000 K)", fontsize=10)
    save(fig, d, "mess_four_methods_rates.png")

    # ---- A2: fall-off curves at four temperatures ----
    fig, axes = plt.subplots(2, 4, figsize=(20, 9.5))
    for row, q in enumerate(("assoc", "dissoc")):
        for col, t in enumerate((500.0, 1000.0, 1500.0, 2000.0)):
            ax = axes[row][col]
            pts = [(p / TORR_PER_ATM, v) for p, v in at_t(M[q], t)]
            hp_mess = mess_hp[t]["P1"]["W1"] if q == "assoc" else mess_hp[t]["W1"]["P1"]
            mess_line(ax, pts + [(25.0, hp_mess)], log=True)
            for m in ORDER:
                mark(ax, [(p / TORR_PER_ATM, v) for p, v in at_t(V[m][q], t)], m, log=True)
            hp = cse_cap[(t, 760.0)]["capture"] if q == "assoc" else olz_hp[(t, 760.0)]["W1->P1"]
            ax.plot([25.0], [hp], "*", ms=11, color="purple", label="MarXus k$_\\infty$ (deck)")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xticks([0.1, 0.3, 1, 3, 10, 25])
            ax.set_xticklabels(["0.1", "0.3", "1", "3", "10", "$\\infty$"])
            ax.set_xlabel("pressure (atm)")
            ax.set_title(f"{t:.0f} K", fontsize=10)
            if col == 0:
                ax.set_ylabel(names[q])
                ax.legend(fontsize=7)
    fig.suptitle("Fall-off of both directions: the four methods and MESS (x)", fontsize=11)
    save(fig, d, "mess_four_methods_falloff.png")

    # ---- A3: stabilization and prompt redissociation (chemical activation) yields ----
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    for row, (q, what) in enumerate((("stab", "stabilization of C$_2$H$_3$"),
                                     ("return", "prompt redissociation to H + C$_2$H$_2$"))):
        for col, p in enumerate((76.0, 7600.0)):
            ax = axes[row][col]
            mess_line(ax, [(t, 100 * v) for t, v in at_p(M[q], p)], log=True)
            for m in ORDER:
                mark(ax, [(t, 100 * v) for t, v in at_p(V[m][q], p)], m, log=True)
            ax.set_yscale("log")
            ax.set_xlabel("T (K)")
            ax.set_ylabel(f"{what} (% of the formed adducts)")
            ax.set_title(f"{p / TORR_PER_ATM:g} atm (MESS: k(P1$\\rightarrow$W1)/k$_\\infty$ of its own tables)", fontsize=10)
            ax.legend(fontsize=7)
        ax = axes[row][2]
        for m in ORDER:
            band(ax, deviation(V[m][q], M[q]), STYLE[m]["color"], STYLE[m]["name"], STYLE[m]["marker"])
            rec.add("mess_four_methods_yields", q, STYLE[m]["name"], V[m][q], M[q], "MESS")
        percent_axes(ax)
        ax.set_ylabel("MarXus / MESS - 1 (%)")
        ax.set_title(f"{what}: deviation from MESS", fontsize=10)
        ax.legend(fontsize=7)
    fig.suptitle("Yields of the chemically activated C$_2$H$_3$: the four methods against MESS", fontsize=11)
    save(fig, d, "mess_four_methods_yields.png")

    # ---- A4: deviation per method and pressure (rows: methods) ----
    fig, axes = plt.subplots(4, 2, figsize=(12, 16), sharex=True)
    pc = plt.cm.plasma(np.linspace(0.0, 0.85, len(pressures)))
    for row, m in enumerate(ORDER):
        for col, q in enumerate(("assoc", "dissoc")):
            ax = axes[row][col]
            dev = deviation(V[m][q], M[q])
            for p, c in zip(pressures, pc):
                pts = at_p(dev, p)
                if pts:
                    ax.plot(*zip(*pts), "-" + STYLE[m]["marker"], ms=5, mfc="none", color=c,
                            label=f"{p / TORR_PER_ATM:g} atm")
            percent_axes(ax, "T (K)" if row == 3 else "")
            ax.set_ylabel("MarXus / MESS - 1 (%)")
            ax.set_title(f"{STYLE[m]['name']}: {'association' if q == 'assoc' else 'dissociation'}", fontsize=10)
            if not dev:
                ax.text(0.5, 0.5, "not available", transform=ax.transAxes, ha="center")
            elif row == 0 and col == 0:
                ax.legend(fontsize=7)
    fig.suptitle("Deviation of each method from MESS, per pressure (green band: $\\pm$5%)", fontsize=11)
    save(fig, d, "mess_four_methods_deviation.png")

    # ---- B1: internal: rates of the four methods against SteadyStateOlzmann ----
    fig, axes = plt.subplots(2, 2, figsize=(15, 10.5))
    for col, q in enumerate(("assoc", "dissoc")):
        ax = axes[0][col]
        for m in ("bar", "cse", "ti"):
            band(ax, deviation(V[m][q], V["olz"][q]), STYLE[m]["color"], f"{STYLE[m]['name']} / SteadyStateOlzmann",
                 STYLE[m]["marker"])
            rec.add("internal_four_methods_rates", q, STYLE[m]["name"], V[m][q], V["olz"][q], "SteadyStateOlzmann")
        percent_axes(ax)
        ax.set_ylabel("method / SteadyStateOlzmann - 1 (%)")
        ax.set_title(f"{'association' if q == 'assoc' else 'dissociation'}: each method against SteadyStateOlzmann "
                     f"({'k$_{uni}$ K' if q == 'assoc' else 'k$_{uni}$'})", fontsize=10)
        ax.legend(fontsize=7)
        ax = axes[1][col]
        pairs = [("bar", "olz"), ("cse", "olz"), ("ti", "olz"), ("ti", "cse")]
        for i, (a, b) in enumerate(pairs):
            dev = abs_deviation(V[a][q], V[b][q])
            ok = [(t + 6 * (i - 1.5), v) for (t, p), v in dev.items() if (t, p) not in below_floor]
            floor = [(t + 6 * (i - 1.5), v) for (t, p), v in dev.items() if (t, p) in below_floor]
            label = f"{STYLE[a]['name']} vs {STYLE[b]['name']}"
            if (a, b) == ("ti", "cse"):
                if ok:
                    ax.plot(*zip(*ok), "*", ms=9, color="k", label=label + " (identity)")
            else:
                mark(ax, ok, a, label=label)
            if floor:
                ax.plot(*zip(*floor), "x", ms=5, color="0.55")
            if (a, b) == ("ti", "cse"):  # the other pairs are recorded with the upper panel
                rec.add("internal_four_methods_rates", q, "TimeIntegration", V[a][q], V[b][q], "CSE")
        ax.plot([], [], "x", color="0.55", label="300-500 K: $\\lambda_1$ below the double-precision floor")
        log_dev_axes(ax)
        ax.set_title("the same, |a/b - 1| for every condition (identities: CSE = SteadyStateOlzmann for the "
                     "dissociation;\nTimeIntegration = CSE for the association, and = k$_{uni}$ for the dissociation)",
                     fontsize=9)
        ax.legend(fontsize=7)
    fig.suptitle("H + C$_2$H$_2$ $\\rightleftharpoons$ C$_2$H$_3$: the four methods against each other", fontsize=11)
    save(fig, d, "internal_four_methods_rates.png")

    # ---- B2: internal: yields (stabilization, prompt redissociation) and chemical activation k_ca ----
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.8))
    tc = {500.0: "tab:purple", 1000.0: "tab:blue", 1500.0: "tab:olive", 2000.0: "tab:red"}
    for ax, q, what in ((axes[0], "stab", "stabilization of C$_2$H$_3$"),
                        (axes[1], "return", "prompt redissociation (chemical activation)")):
        for t, c in tc.items():
            for m in ORDER:
                pts = [(p / TORR_PER_ATM, 100 * v) for p, v in at_t(V[m][q], t)]
                s = STYLE[m]
                pts = [(x, y) for x, y in pts if y > 0]
                if pts:
                    ax.plot(*zip(*pts), ls="-" if m == "olz" else "none", lw=0.8, marker=s["marker"], ms=s["ms"],
                            mec=c, mfc=c if m == "bar" else "none", color=c)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("pressure (atm)")
        ax.set_ylabel(f"{what} (% of the formed adducts)")
        ax.set_title(what, fontsize=10)
        for m in ORDER:
            ax.plot([], [], ls="none", marker=STYLE[m]["marker"], mec="k", mfc="none", label=STYLE[m]["name"])
        for t, c in tc.items():
            ax.plot([], [], "-", color=c, label=f"{t:.0f} K")
        ax.legend(fontsize=7, ncol=2)
        for m in ("bar", "cse", "ti"):
            rec.add("internal_four_methods_yields", q, STYLE[m]["name"], V[m][q], V["olz"][q], "SteadyStateOlzmann")
    ax = axes[2]
    kca = {"bar": {k: r["W1->P1"] for k, r in bar_kca.items()}, "olz": {k: r["W1->P1"] for k, r in olz_kca.items()}}
    for m, ls in (("bar", "-"), ("olz", "--")):
        for t, c in tc.items():
            pts = [(p / TORR_PER_ATM, v) for p, v in at_t(kca[m], t) if v > 0]
            if pts:
                ax.plot(*zip(*pts), ls=ls, marker=STYLE[m]["marker"], ms=5, mfc="none", color=c)
    for m, ls in (("bar", "-"), ("olz", "--")):
        ax.plot([], [], ls=ls, marker=STYLE[m]["marker"], color="k", mfc="none",
                label=f"{STYLE[m]['name']} ({'intermediate' if m == 'bar' else 'final'} steady state)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("pressure (atm)")
    ax.set_ylabel("k$_{ca}$(C$_2$H$_3$ $\\rightarrow$ C$_2$H$_2$ + H) (s$^{-1}$)")
    ax.set_title("chemical-activation rate coefficient k$_{ca}$ (GO10 eq. 9): only the steady states define it\n"
                 "(final: the equilibrium, pressure independent)", fontsize=9)
    ax.legend(fontsize=7)
    fig.suptitle("Yields and chemical activation of the four methods (no MESS)", fontsize=11)
    save(fig, d, "internal_four_methods_yields.png")

    # ---- B3: internal: time evolution of the pulse against the phenomenological model of each method ----
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.8))
    for ax, t in zip(axes, (300.0, 1000.0, 2000.0)):
        key = (t, 760.0)
        rows = evolution[key][1]
        ts = np.array([r["t[s]"] for r in rows])
        ax.loglog(ts, [100 * r["N(W1)"] for r in rows], "-", color=STYLE["ti"]["color"], lw=2,
                  label="TimeIntegration: N(C$_2$H$_3$)")
        ax.loglog(ts, [max(100 * r["W1->P1"], 1e-12) for r in rows], ":", color=STYLE["ti"]["color"], lw=2,
                  label="TimeIntegration: H + C$_2$H$_2$ formed back")
        tt = np.logspace(-12, 2, 300)
        models = (("cse", V["cse"]["stab"][key], V["cse"]["dissoc"][key], "--"),
                  ("olz", V["olz"]["stab"][key], V["olz"]["dissoc"][key], "-."))
        for m, a, k, ls in models:
            ax.loglog(tt, np.maximum(100 * a * np.exp(-k * tt), 1e-30), ls, color=STYLE[m]["color"],
                      label=f"{STYLE[m]['name']}: (k(R$\\rightarrow$W)/k$_\\infty$) exp(-k(W$\\rightarrow$P) t)")
        ax.axhline(100 * V["bar"]["stab"][key], color=STYLE["bar"]["color"], lw=1,
                   label="SteadyStateAbsorbingBarrier: $\\Phi_{stab}$")
        ax.set_ylim(1e-6, 200)
        ax.set_xlabel("t (s)")
        ax.set_ylabel("% of the formed adducts")
        ax.set_title(f"{t:.0f} K, 1 atm", fontsize=10)
        ax.legend(fontsize=6.5, loc="lower left")
    fig.suptitle("Pulse of chemically activated C$_2$H$_3$ (TimeIntegration) against the two-state model of the "
                 "rate coefficients of CSE and SteadyStateOlzmann", fontsize=11)
    save(fig, d, "internal_four_methods_time.png")
    return d, rec


# ----------------------------------------------------------------------------------------------
# ZZ-allyl + O2, Case 2 (four wells; MESS Eckart model, as MESS)
# ----------------------------------------------------------------------------------------------
WELLS = ["G2", "G3", "G4", "G6"]
ENDS = ["R", "P1", "P5", "P7", "escape(G4)"]
PRODUCTS = ["P5", "escape(G4)", "P1", "P7"]
EXIT = {"R": "G2->R", "P1": "G4->P1", "P5": "G4->P5", "P7": "G6->P7", "escape(G4)": "escape(G4)"}
LABEL = {"P5": "IEPOX + OH (P5)", "escape(G4)": "escape (G4)", "P1": "P1", "P7": "P7", "R": "back to R",
         "G2": "G2", "G3": "G3", "G4": "G4"}


def read_cse_eigenvalues(path):
    """Chemical eigenvalues of the CSE machine file: {(T, p_torr): [Lambda_1, ...]}."""
    out = {}
    for line in open(path):
        m = re.match(r"# T = (\S+) K, p = (\S+) Torr: chemical eigenvalues \[1/s\]: ([^;]+);", line)
        if m:
            out[(float(m.group(1)), float(m.group(2)))] = [float(x) for x in m.group(3).split(",")]
    return out


def case2():
    d = os.path.join(HERE, "ZZAllyl+O2_Gamma_Case2")
    o = os.path.join(d, "marxus_output")
    stem = "case2_tstlevel_E_mess_eckart"
    tabs = {m: read_tables(os.path.join(o, f"{stem}_{m}_tables.csv"))
            for m in ("olzmann", "absorbing_barrier", "cse", "time_integration")}
    evolution = read_evolution(os.path.join(o, f"{stem}_time_integration.csv"))
    eigen = read_cse_eigenvalues(os.path.join(o, f"{stem}_cse.csv"))
    mess = read_mess(os.path.join(d, "Gamma-Case2_from_ZZ_allyl+O2_Gamma4-Escape_CCSDT_version-1_2025_09_29_12.7kcal.out"))
    rec = Record()

    bar_y = table(tabs["absorbing_barrier"], "intermediate steady state: Yields (% of the formed adducts)")
    bar_net = table(tabs["absorbing_barrier"], "intermediate steady state: Yields without the return to R")
    bar_rp = table(tabs["absorbing_barrier"], "intermediate steady state: Bimolecular-to-bimolecular rate coefficients")
    bar_rw = table(tabs["absorbing_barrier"], "intermediate steady state: Bimolecular-to-well rate coefficients")
    bar_cap = table(tabs["absorbing_barrier"], "intermediate steady state: Capture, return and net reaction")
    olz_y = table(tabs["olzmann"], "final steady state: Yields (% of the formed adducts)")
    olz_net = table(tabs["olzmann"], "final steady state: Yields without the return to R")
    olz_rp = table(tabs["olzmann"], "final steady state: Bimolecular-to-bimolecular rate coefficients, overall")
    olz_th = table(tabs["olzmann"], "thermal rate coefficients: Thermal rate coefficients")
    fates = {w: table(tabs["olzmann"], f"thermal fates: Thermal fate of the molecules thermalized in {w}") for w in WELLS}
    cse_rp = table(tabs["cse"], "CSE: Bimolecular-to-bimolecular rate coefficients")
    cse_rw = table(tabs["cse"], "CSE: Bimolecular-to-well rate coefficients")
    cse_cap = table(tabs["cse"], "CSE: Capture, return and net reaction")
    cse_by = table(tabs["cse"], "CSE: Bimolecular-to-bimolecular yields")
    cse_wy = table(tabs["cse"], "CSE: Bimolecular-to-well yields")
    cse_fate = {w: table(tabs["cse"], f"CSE: Thermal fate of {w}") for w in WELLS}
    cse_from = {w: table(tabs["cse"], f"CSE: Rate coefficients from {w}") for w in WELLS}
    cse_direct = table(tabs["cse"], "CSE: Long-time yields, direct")
    cse_via = table(tabs["cse"], "CSE: Long-time yields, through the wells")
    cse_long = table(tabs["cse"], "CSE: Long-time yields, total")
    ti_net = table(tabs["time_integration"], "time integration: Yields without the return to R")
    ti_rp = table(tabs["time_integration"], "time integration: Bimolecular-to-bimolecular rate coefficients, overall")
    keys = sorted(olz_net)
    pressures = sorted({p for _, p in keys})

    def series():
        return {m: {} for m in ORDER}

    # Long-time shares (% of the eventual net reaction), overall rate coefficients R -> X (cm3/s), prompt branching
    # (% of the net reaction), stabilization rate coefficients R -> W, thermal fates (%).
    share, overall, prompt, chem, thermal_part = ({x: series() for x in PRODUCTS} for _ in range(5))
    stab_rate = {w: series() for w in ["G2", "G3", "G4"]}
    prompt_w = {w: series() for w in ["G2", "G3", "G4"]}
    fate = {(w, x): series() for w in ["G2", "G3", "G4"] for x in ["R", "P5", "escape(G4)"]}
    M_share, M_overall, M_prompt, M_rp, M_fate = ({} for _ in range(5))
    M_rw = {w: {} for w in ["G2", "G3", "G4"]}
    M_prompt_w = {w: {} for w in ["G2", "G3", "G4"]}
    M_fate = {(w, x): {} for w in ["G2", "G3", "G4"] for x in ["R", "P5", "escape(G4)"]}
    M_well = {}
    for key in keys:
        # SteadyStateAbsorbingBarrier with the thermal fates of SteadyStateOlzmann: prompt + stabilization x fate.
        total = {x: bar_y[key][EXIT[x]] + sum(bar_y[key][f"stab({w})"] * fates[w][key][EXIT[x]] / 100 for w in WELLS)
                 for x in PRODUCTS}
        eventual = sum(total.values())
        for x in PRODUCTS:
            share[x]["olz"][key] = olz_net[key][EXIT[x]]
            share[x]["cse"][key] = cse_long[key][f"R->{x}"]
            share[x]["ti"][key] = ti_net[key][EXIT[x]]
            share[x]["bar"][key] = 100 * total[x] / eventual
            overall[x]["olz"][key] = olz_rp[key][f"R->{x}"]
            overall[x]["ti"][key] = ti_rp[key][f"R->{x}"]
            overall[x]["cse"][key] = cse_rp[key][f"R->{x}"] + sum(cse_rw[key][f"R->{w}"] * cse_fate[w][key][f"{w}->{x}"] / 100
                                                                  for w in WELLS)
            overall[x]["bar"][key] = bar_cap[key]["capture"] * total[x] / 100
            prompt[x]["bar"][key] = bar_net[key][EXIT[x]]
            prompt[x]["cse"][key] = cse_by[key][f"R->{x}"]
            chem[x]["bar"][key] = 100 * bar_y[key][EXIT[x]] / eventual
            chem[x]["cse"][key] = cse_direct[key][f"R->{x}"]
            thermal_part[x]["bar"][key] = share[x]["bar"][key] - chem[x]["bar"][key]
            thermal_part[x]["cse"][key] = cse_via[key][f"R->{x}"]
        for w in ["G2", "G3", "G4"]:
            stab_rate[w]["bar"][key] = bar_rw[key][f"R->{w}"]
            stab_rate[w]["cse"][key] = cse_rw[key][f"R->{w}"]
            prompt_w[w]["bar"][key] = bar_net[key][f"stab({w})"]
            prompt_w[w]["cse"][key] = cse_wy[key][f"R->{w}"]
            for x in ["R", "P5", "escape(G4)"]:
                fate[(w, x)]["olz"][key] = fates[w][key][EXIT[x]]
                fate[(w, x)]["cse"][key] = cse_fate[w][key][f"{w}->{x}"]
        # MESS: rate tables, prompt branching, thermal fates (absorbing chain of its well rows), long-time fate.
        if key in mess:
            m = mess[key]
            m_fates = absorbing_chain(m, WELLS, ENDS)
            net = sum(m["R"][x] for x in WELLS + PRODUCTS)
            final = {x: m["R"][x] + sum(m["R"][w] * m_fates[w][x] for w in WELLS) for x in PRODUCTS}
            eventual_m = sum(final.values())
            for x in PRODUCTS:
                M_share.setdefault(x, {})[key] = 100 * final[x] / eventual_m
                M_overall.setdefault(x, {})[key] = final[x]
                M_prompt.setdefault(x, {})[key] = 100 * m["R"][x] / net
                M_rp.setdefault(x, {})[key] = m["R"][x]
            for w in ["G2", "G3", "G4"]:
                M_rw[w][key] = m["R"][w]
                M_prompt_w[w][key] = 100 * m["R"][w] / net
                for x in ["R", "P5", "escape(G4)"]:
                    M_fate[(w, x)][key] = 100 * m_fates[w][x]
            M_well[key] = m

    # ---- A1: rates against MESS (760 Torr absolute; deviations over all pressures) ----
    fig, axes = plt.subplots(2, 3, figsize=(19, 10.5))
    ax = axes[0][0]
    mess_line(ax, at_p(M_rp["P5"], 760.0))
    for m, rows in (("bar", bar_rp), ("cse", cse_rp)):
        mark(ax, at_p({k: r["R->P5"] for k, r in rows.items()}, 760.0), m)
    ax.set_title("bimolecular-to-bimolecular (direct) R $\\rightarrow$ IEPOX + OH, 760 Torr\n"
                 "(SteadyStateOlzmann and TimeIntegration give no direct part)", fontsize=9)
    ax.set_ylabel("k (cm$^3$ s$^{-1}$)")
    ax = axes[0][1]
    for x, c in (("P5", "tab:red"), ("escape(G4)", "tab:purple")):
        mess_line(ax, at_p(M_overall[x], 760.0), label=f"MESS {LABEL[x]}", color=c, log=True)
        for m in ORDER:
            mark(ax, at_p(overall[x][m], 760.0), m, label=f"{STYLE[m]['name']} {LABEL[x]}", log=True)
    ax.set_yscale("log")
    ax.set_title("overall (chemical activation + thermal) R $\\rightarrow$ X, 760 Torr\n"
                 "(MESS: its R row with the thermal fates of its well rows; AbsorbingBarrier: fates of Olzmann)", fontsize=9)
    ax.set_ylabel("k (cm$^3$ s$^{-1}$)")
    ax = axes[0][2]
    for w, c in (("G2", "tab:blue"), ("G3", "tab:orange"), ("G4", "tab:green")):
        mess_line(ax, at_p(M_rw[w], 760.0), label=f"MESS R$\\rightarrow${w}", color=c, log=True)
        for m in ("bar", "cse"):
            mark(ax, at_p(stab_rate[w][m], 760.0), m, label=f"{STYLE[m]['name']} R$\\rightarrow${w}", log=True)
    ax.set_yscale("log")
    ax.set_title("bimolecular-to-well (stabilization), 760 Torr\n(SteadyStateOlzmann, TimeIntegration: no stabilization "
                 "of their own in a multiwell network)", fontsize=9)
    ax.set_ylabel("k (cm$^3$ s$^{-1}$)")
    for ax in axes[0]:
        ax.set_xlabel("T (K)")
        ax.legend(fontsize=6)
    ax = axes[1][0]
    for m, rows in (("bar", bar_rp), ("cse", cse_rp)):
        vals = {k: r["R->P5"] for k, r in rows.items()}
        band(ax, deviation(vals, M_rp["P5"]), STYLE[m]["color"], STYLE[m]["name"], STYLE[m]["marker"])
        rec.add("mess_four_methods_rates", "R->P5 direct", STYLE[m]["name"], vals, M_rp["P5"], "MESS")
    ax.set_title("direct R $\\rightarrow$ IEPOX + OH: deviation from MESS", fontsize=10)
    ax = axes[1][1]
    for x, ls in (("P5", "-"), ("escape(G4)", "--")):
        for m in ORDER:
            band(ax, deviation(overall[x][m], M_overall[x]), STYLE[m]["color"], f"{STYLE[m]['name']} {LABEL[x]}",
                 STYLE[m]["marker"], ls)
            rec.add("mess_four_methods_rates", f"R->{x} overall", STYLE[m]["name"], overall[x][m], M_overall[x], "MESS")
    ax.set_title("overall R $\\rightarrow$ X: deviation from MESS", fontsize=10)
    ax = axes[1][2]
    for w, ls in (("G2", "-"), ("G3", ":"), ("G4", "--")):
        for m in ("bar", "cse"):
            band(ax, deviation(stab_rate[w][m], M_rw[w]), STYLE[m]["color"], f"{STYLE[m]['name']} R$\\rightarrow${w}",
                 STYLE[m]["marker"], ls)
            rec.add("mess_four_methods_rates", f"R->{w}", STYLE[m]["name"], stab_rate[w][m], M_rw[w], "MESS")
    ax.set_title("bimolecular-to-well: deviation from MESS", fontsize=10)
    for ax in axes[1]:
        percent_axes(ax)
        ax.set_ylabel("MarXus / MESS - 1 (%)")
        ax.legend(fontsize=6)
    fig.suptitle("ZZ-allyl + O$_2$ (Case 2, MESS Eckart model): rate coefficients of the four methods against MESS",
                 fontsize=11)
    save(fig, d, "mess_four_methods_rates.png")

    # ---- A2: yields against MESS: long-time shares (all methods) and prompt branching (steady state, CSE) ----
    fig, axes = plt.subplots(2, 2, figsize=(16, 11))
    ax = axes[0][0]
    for x, c in zip(PRODUCTS, ("tab:red", "tab:purple", "tab:green", "tab:gray")):
        mess_line(ax, at_p(M_share[x], 760.0), label=f"MESS {LABEL[x]}", color=c, log=True)
        for m in ORDER:
            mark(ax, at_p(share[x][m], 760.0), m, label=f"{STYLE[m]['name']} {LABEL[x]}" if x == "P5" else "", log=True)
    ax.set_yscale("log")
    ax.set_ylabel("long-time yield (% of the eventual net reaction)")
    ax.set_title("long-time yields at 760 Torr (AbsorbingBarrier: prompt + stabilization x thermal fate of Olzmann)",
                 fontsize=9)
    ax = axes[1][0]
    for x, ls in zip(PRODUCTS, ("-", "--", ":", "-.")):
        for m in ORDER:
            band(ax, deviation(share[x][m], M_share[x]), STYLE[m]["color"], f"{STYLE[m]['name']} {LABEL[x]}",
                 STYLE[m]["marker"], ls)
            rec.add("mess_four_methods_yields", f"long-time {x}", STYLE[m]["name"], share[x][m], M_share[x], "MESS")
    percent_axes(ax)
    ax.set_ylabel("MarXus / MESS - 1 (%)")
    ax.set_title("long-time yields: deviation from MESS", fontsize=10)
    ax = axes[0][1]
    for x, c in (("P5", "tab:red"), ("G2", "tab:blue"), ("G3", "tab:orange"), ("G4", "tab:green")):
        mref = M_prompt[x] if x in PRODUCTS else M_prompt_w[x]
        mess_line(ax, at_p(mref, 760.0), label=f"MESS R$\\rightarrow${x}", color=c, log=True)
        for m in ("bar", "cse"):
            vals = prompt[x][m] if x in PRODUCTS else prompt_w[x][m]
            mark(ax, at_p(vals, 760.0), m, label=f"{STYLE[m]['name']} R$\\rightarrow${x}", log=True)
    ax.set_yscale("log")
    ax.set_ylabel("prompt yield (% of the net reaction of R)")
    ax.set_title("prompt branching of R (chemical activation and stabilization), 760 Torr\n(MESS: its R row; the "
                 "direct R $\\rightarrow$ escape is negative in MESS and CSE and not shown)", fontsize=9)
    ax = axes[1][1]
    for x, ls in (("P5", "-"), ("G2", "--"), ("G3", ":"), ("G4", "-.")):
        mref = M_prompt[x] if x in PRODUCTS else M_prompt_w[x]
        for m in ("bar", "cse"):
            vals = prompt[x][m] if x in PRODUCTS else prompt_w[x][m]
            band(ax, deviation(vals, mref), STYLE[m]["color"], f"{STYLE[m]['name']} R$\\rightarrow${x}",
                 STYLE[m]["marker"], ls)
            rec.add("mess_four_methods_yields", f"prompt R->{x}", STYLE[m]["name"], vals, mref, "MESS")
    percent_axes(ax)
    ax.set_ylabel("MarXus / MESS - 1 (%)")
    ax.set_title("prompt branching: deviation from MESS", fontsize=10)
    for ax in axes.flat:
        ax.set_xlabel("T (K)")
        ax.legend(fontsize=6, ncol=2)
    fig.suptitle("ZZ-allyl + O$_2$ (Case 2, MESS Eckart model): yields of the four methods against MESS", fontsize=11)
    save(fig, d, "mess_four_methods_yields.png")

    # ---- A3: thermal activation against MESS: thermal fates of the wells, well rate coefficients (CSE) ----
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    ax = axes[0]
    colors = {"R": "tab:brown", "P5": "tab:red", "escape(G4)": "tab:purple"}
    for (w, x), mk in zip([(w, x) for w in ["G2", "G3", "G4"] for x in ["R", "P5", "escape(G4)"]],
                          ["-"] * 3 + ["--"] * 3 + [":"] * 3):
        pts = [(t, v) for t, v in at_p(M_fate[(w, x)], 760.0) if v > 1e-6]
        if pts:
            ax.plot(*zip(*pts), mk, marker="x", ms=9, mew=2, color=colors[x], label=f"MESS {w} $\\rightarrow$ {LABEL[x]}")
        for m in ("olz", "cse"):
            mark(ax, [(t, v) for t, v in at_p(fate[(w, x)][m], 760.0) if v > 1e-6], m, label="", log=True)
            rec.add("mess_four_methods_thermal", f"fate {w}->{x}", STYLE[m]["name"], fate[(w, x)][m], M_fate[(w, x)], "MESS")
    for m in ("olz", "cse"):
        ax.plot([], [], ls="none", marker=STYLE[m]["marker"], mec=STYLE[m]["color"], mfc="none", label=STYLE[m]["name"])
    ax.set_yscale("log")
    ax.set_xlabel("T (K)")
    ax.set_ylabel("thermal fate (% of the molecules thermalized in the well)")
    ax.set_title("thermal fates at 760 Torr (solid G2, dashed G3, dotted G4)\nMESS: absorbing chain of its well rows; "
                 "AbsorbingBarrier and TimeIntegration have none", fontsize=9)
    ax.legend(fontsize=6, ncol=2)
    ax = axes[1]
    for x, ls in (("R", "-"), ("P5", "--"), ("escape(G4)", ":")):
        for w in ["G2", "G3", "G4"]:
            for m in ("olz", "cse"):
                dev = deviation(fate[(w, x)][m], M_fate[(w, x)])
                dev = {k: v for k, v in dev.items() if M_fate[(w, x)][k] > 1e-3}
                band(ax, dev, STYLE[m]["color"], f"{STYLE[m]['name']} {w}$\\rightarrow${LABEL[x]}", STYLE[m]["marker"], ls)
    percent_axes(ax)
    ax.set_ylabel("MarXus / MESS - 1 (%)")
    ax.set_title("thermal fates (> 0.001%): deviation from MESS", fontsize=10)
    ax.legend(fontsize=5.5, ncol=2)
    ax = axes[2]
    pairs = [("G4", "P5"), ("G4", "escape(G4)"), ("G4", "P1"), ("G2", "R"), ("G2", "G4"), ("G4", "G2"), ("G3", "G2")]
    for (w, x), c in zip(pairs, plt.cm.tab10(np.linspace(0, 0.9, len(pairs)))):
        vals = {k: cse_from[w][k][f"{w}->{x}"] for k in keys}
        ref = {k: M_well[k][w][x] for k in keys if k in M_well}
        band(ax, deviation(vals, ref), c, f"CSE {w}$\\rightarrow${x}", "s")
        rec.add("mess_four_methods_thermal", f"k({w}->{x})", "CSE", vals, ref, "MESS")
    percent_axes(ax)
    ax.set_ylabel("CSE / MESS - 1 (%)")
    ax.set_title("thermal rate coefficients between the species (only CSE gives them)", fontsize=10)
    ax.legend(fontsize=7)
    fig.suptitle("ZZ-allyl + O$_2$ (Case 2, MESS Eckart model): thermal activation against MESS", fontsize=11)
    save(fig, d, "mess_four_methods_thermal.png")

    # ---- A4: each method against MESS, every quantity it gives (one panel per method) ----
    fig, axes = plt.subplots(1, 4, figsize=(22, 6), sharey=True)
    quantities = {m: [] for m in ORDER}
    for x in PRODUCTS:
        for m in ORDER:
            quantities[m].append((f"long-time {LABEL[x]}", deviation(share[x][m], M_share[x])))
    for m, rows in (("bar", bar_rp), ("cse", cse_rp)):
        quantities[m].append(("direct R$\\rightarrow$P5", deviation({k: r["R->P5"] for k, r in rows.items()}, M_rp["P5"])))
        for w in ["G2", "G3", "G4"]:
            quantities[m].append((f"R$\\rightarrow${w}", deviation(stab_rate[w][m], M_rw[w])))
    for w in ["G2", "G3", "G4"]:
        for x in ["P5", "escape(G4)"]:
            for m in ("olz", "cse"):
                quantities[m].append((f"fate {w}$\\rightarrow${LABEL[x]}", deviation(fate[(w, x)][m], M_fate[(w, x)])))
    for ax, m in zip(axes, ORDER):
        for (label, dev), c in zip(quantities[m], plt.cm.tab20(np.linspace(0, 1, max(len(quantities[m]), 2)))):
            band(ax, dev, c, label, STYLE[m]["marker"])
        percent_axes(ax)
        ax.set_title(STYLE[m]["name"], fontsize=10)
        ax.legend(fontsize=5.5, ncol=1)
    axes[0].set_ylabel("MarXus / MESS - 1 (%)")
    fig.suptitle("ZZ-allyl + O$_2$ (Case 2, MESS Eckart model): each method against MESS, every quantity it gives",
                 fontsize=11)
    save(fig, d, "mess_four_methods_deviation.png")

    # ---- B1: internal: long-time yields of the four methods against SteadyStateOlzmann ----
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    ax = axes[0]
    for x, c in zip(PRODUCTS, ("tab:red", "tab:purple", "tab:green", "tab:gray")):
        pts = at_p(share[x]["olz"], 760.0)
        ax.plot(*zip(*pts), "-", color=c, lw=1, label=f"{LABEL[x]}")
        for m in ORDER:
            mark(ax, at_p(share[x][m], 760.0), m, label=STYLE[m]["name"] if x == "P5" else "", log=True)
    ax.set_yscale("log")
    ax.set_xlabel("T (K)")
    ax.set_ylabel("long-time yield (% of the eventual net reaction)")
    ax.set_title("long-time yields at 760 Torr", fontsize=10)
    ax.legend(fontsize=7, ncol=2)
    ax = axes[1]
    for x, off in zip(PRODUCTS, (-4.5, -1.5, 1.5, 4.5)):
        for m in ("bar", "cse", "ti"):
            dev = abs_deviation(share[x][m], share[x]["olz"])
            mark(ax, [(t + off, v) for (t, p), v in dev.items()], m,
                 label=f"{STYLE[m]['name']} / Olzmann" if x == "P5" else "")
            rec.add("internal_four_methods_yields", f"long-time {x}", STYLE[m]["name"], share[x][m], share[x]["olz"],
                    "SteadyStateOlzmann")
    log_dev_axes(ax)
    ax.set_title("each method against SteadyStateOlzmann, |a/b - 1| (P5, escape, P1, P7 side by side)\nCSE and "
                 "TimeIntegration: identities; AbsorbingBarrier: thermalized-before-reaction assumption", fontsize=9)
    ax.legend(fontsize=7)
    fig.suptitle("ZZ-allyl + O$_2$ (Case 2): long-time yields of the four methods against each other", fontsize=11)
    save(fig, d, "internal_four_methods_yields.png")

    # ---- B2: internal: total formation = chemical activation + thermal, per method ----
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    for ax, x in zip(axes, ("P5", "escape(G4)")):
        for m, part, ls in (("bar", chem, "-"), ("cse", chem, "-")):
            ax.plot(*zip(*at_p(part[x][m], 760.0)), ls, marker=STYLE[m]["marker"], mfc="none", color=STYLE[m]["color"],
                    label=f"{STYLE[m]['name']}: chemical activation (prompt / direct)")
        for m in ("bar", "cse"):
            ax.plot(*zip(*at_p(thermal_part[x][m], 760.0)), ":", marker=STYLE[m]["marker"], mfc="none",
                    color=STYLE[m]["color"], label=f"{STYLE[m]['name']}: thermal (through the stabilized wells)")
        for m in ORDER:
            ax.plot(*zip(*at_p(share[x][m], 760.0)), "--", marker=STYLE[m]["marker"], ms=STYLE[m]["ms"], mfc="none",
                    color=STYLE[m]["color"], lw=0.8, label=f"{STYLE[m]['name']}: total")
            rec.add("internal_four_methods_formation", f"total {x}", STYLE[m]["name"], share[x][m], share[x]["olz"],
                    "SteadyStateOlzmann")
        for m in ("bar", "cse"):
            rec.add("internal_four_methods_formation", f"chemical activation {x}", STYLE[m]["name"], chem[x][m],
                    chem[x]["cse"], "CSE")
            rec.add("internal_four_methods_formation", f"thermal {x}", STYLE[m]["name"], thermal_part[x][m],
                    thermal_part[x]["cse"], "CSE")
        if x == "escape(G4)":
            lo, hi = min(chem[x]["cse"].values()), max(chem[x]["cse"].values())
            ax.text(0.02, 0.40, f"CSE's direct R $\\rightarrow$ escape is negative ({lo:.2f} ... {hi:.2f}%, as MESS's R $\\rightarrow$"
                    " escape entry)\nand is not shown on the log scale; its thermal part is then above the total",
                    transform=ax.transAxes, fontsize=7)
        ax.set_yscale("log")
        ax.set_xlabel("T (K)")
        ax.set_ylabel("% of the eventual net reaction")
        ax.set_title(f"{LABEL[x]} at 760 Torr: chemical activation + thermal = total\n(AbsorbingBarrier: prompt and "
                     "stabilization x thermal fate of Olzmann; CSE: direct and through the wells)", fontsize=9)
        ax.legend(fontsize=6.5)
    fig.suptitle("ZZ-allyl + O$_2$ (Case 2): total formation yield split into chemical activation and thermal reaction",
                 fontsize=11)
    save(fig, d, "internal_four_methods_formation.png")

    # ---- B3: internal: rate coefficients ----
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    ax = axes[0]
    for x, ls in (("P5", "-"), ("P1", ":"), ("P7", "-.")):
        vb = {k: r[f"R->{x}"] for k, r in bar_rp.items()}
        vc = {k: r[f"R->{x}"] for k, r in cse_rp.items()}
        band(ax, deviation(vb, vc), STYLE["bar"]["color"], f"direct R$\\rightarrow${x}", "o", ls)
        rec.add("internal_four_methods_rates", f"direct R->{x}", "SteadyStateAbsorbingBarrier", vb, vc, "CSE")
    for w, ls in (("G2", "-"), ("G3", ":"), ("G4", "--")):
        band(ax, deviation(stab_rate[w]["bar"], stab_rate[w]["cse"]), "tab:cyan", f"R$\\rightarrow${w}", "o", ls)
        rec.add("internal_four_methods_rates", f"R->{w}", "SteadyStateAbsorbingBarrier", stab_rate[w]["bar"],
                stab_rate[w]["cse"], "CSE")
    percent_axes(ax)
    ax.set_ylabel("SteadyStateAbsorbingBarrier / CSE - 1 (%)")
    ax.set_title("prompt and stabilization rate coefficients: absorbing barrier against CSE\n(the only two methods "
                 "that separate them)", fontsize=9)
    ax.legend(fontsize=7)
    ax = axes[1]
    for x, off in (("P5", -2), ("escape(G4)", 2)):
        for m in ("bar", "cse", "ti"):
            dev = abs_deviation(overall[x][m], overall[x]["olz"])
            mark(ax, [(t + off, v) for (t, p), v in dev.items()], m,
                 label=f"{STYLE[m]['name']} / Olzmann" if x == "P5" else "")
            rec.add("internal_four_methods_rates", f"overall R->{x}", STYLE[m]["name"], overall[x][m], overall[x]["olz"],
                    "SteadyStateOlzmann")
    log_dev_axes(ax)
    ax.set_title("overall R $\\rightarrow$ P5 (left) and R $\\rightarrow$ escape (right): each method against "
                 "SteadyStateOlzmann", fontsize=9)
    ax.legend(fontsize=7)
    ax = axes[2]
    lam_cse = {k: v[0] for k, v in eigen.items()}
    k_uni = {k: r["k_uni"] for k, r in olz_th.items()}
    lam_1 = {k: r["lambda_1"] for k, r in olz_th.items()}
    k_ti = {}
    for key, (header, rows) in evolution.items():
        a, r, decayed = slow_mode(rows, WELLS)
        if decayed:
            k_ti[key] = r
    for vals, m, label in ((lam_cse, "cse", "CSE: lowest chemical eigenvalue"),
                           (lam_1, "olz", "SteadyStateOlzmann: $\\lambda_1$"),
                           (k_ti, "ti", "TimeIntegration: late decay rate")):
        dev = abs_deviation(vals, k_uni)
        mark(ax, [(t, v) for (t, p), v in dev.items()], m, label=f"{label} vs k$_{{uni}}$")
        rec.add("internal_four_methods_rates", "thermal decay", label, vals, k_uni, "SteadyStateOlzmann k_uni")
    if not k_ti:
        ax.text(0.03, 0.03, "TimeIntegration: no decay below 1e-3 inside the output window", transform=ax.transAxes,
                fontsize=7)
    log_dev_axes(ax)
    ax.set_title("thermal decay of the network: against k$_{uni}$ of SteadyStateOlzmann", fontsize=10)
    ax.legend(fontsize=7)
    fig.suptitle("ZZ-allyl + O$_2$ (Case 2): rate coefficients of the four methods against each other", fontsize=11)
    save(fig, d, "internal_four_methods_rates.png")

    # ---- B4: internal: thermal fates, SteadyStateOlzmann against CSE ----
    fig, ax = plt.subplots(figsize=(10, 6))
    for (w, x), off in zip([(w, x) for w in ["G2", "G3", "G4"] for x in ["R", "P5", "escape(G4)"]], np.linspace(-4, 4, 9)):
        dev = abs_deviation(fate[(w, x)]["cse"], fate[(w, x)]["olz"])
        dev = {k: v for k, v in dev.items() if fate[(w, x)]["olz"][k] > 1e-3}
        if dev:
            ax.plot([t + off for t, _ in dev], list(dev.values()), "s", mfc="none", ms=5,
                    label=f"{w} $\\rightarrow$ {LABEL[x]}")
        rec.add("internal_four_methods_thermal", f"fate {w}->{x}", "CSE", fate[(w, x)]["cse"], fate[(w, x)]["olz"],
                "SteadyStateOlzmann")
    log_dev_axes(ax)
    ax.set_title("thermal fates of the wells (> 0.001%): CSE (absorbing chain of its rate coefficients) against "
                 "SteadyStateOlzmann\n(J N = thermal source in the well, all conditions)", fontsize=9)
    ax.legend(fontsize=7, ncol=3)
    fig.tight_layout()
    fig.savefig(os.path.join(d, "plots", "internal_four_methods_thermal.png"), dpi=200)
    plt.close(fig)

    # ---- B5: internal: time evolution against the steady-state and CSE yields (300 K, 760 Torr) ----
    key = (300.0, 760.0)
    rows = evolution[key][1]
    ts = [r["t[s]"] for r in rows]
    fig, ax = plt.subplots(figsize=(11, 6.5))
    capture = cse_cap[key]["capture"]
    for x, c in (("P5", "tab:red"), ("escape(G4)", "tab:purple"), ("R", "tab:brown")):
        ax.loglog(ts, [max(100 * r[EXIT[x]], 1e-12) for r in rows], "-", color=c, lw=2,
                  label=f"TimeIntegration: {LABEL[x]}")
        ax.axhline(bar_y[key][EXIT[x]], color=c, ls=":", lw=1.2)
        ax.axhline(olz_y[key][EXIT[x]], color=c, ls="--", lw=1.2)
        if x != "R":
            ax.axhline(100 * cse_rp[key][f"R->{x}"] / capture, color=c, ls="-.", lw=1.0)
    for w, c in (("G2", "tab:blue"), ("G4", "tab:green")):
        ax.loglog(ts, [max(100 * r[f"N({w})"], 1e-12) for r in rows], "-", color=c, lw=1.2, label=f"TimeIntegration: N({w})")
        ax.axhline(bar_y[key][f"stab({w})"], color=c, ls=":", lw=1.2)
    ax.plot([], [], ":", color="k", label="SteadyStateAbsorbingBarrier: prompt yield / stabilization")
    ax.plot([], [], "-.", color="k", label="CSE: direct k(R$\\rightarrow$X)/k$_\\infty$")
    ax.plot([], [], "--", color="k", label="SteadyStateOlzmann: final (long-time) yield")
    ax.set_ylim(1e-6, 200)
    ax.set_xlabel("t (s)")
    ax.set_ylabel("% of the formed adducts")
    ax.set_title("ZZ-allyl + O$_2$, 300 K, 760 Torr: pulse (TimeIntegration) against the yields of the other methods",
                 fontsize=10)
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(d, "plots", "internal_four_methods_time.png"), dpi=200)
    plt.close(fig)
    return d, rec


if __name__ == "__main__":
    for build in (c2h3, case2):
        d, rec = build()
        path = rec.write(d)
        print(f"\n{d}: {path}")
        for (figure, quantity, method, ref), (lo, hi) in rec.summary().items():
            print(f"  {figure:34s} {quantity:24s} {method:42s} vs {ref:28s} {lo:+9.3f} .. {hi:+9.3f} %")
