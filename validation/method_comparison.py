#!/usr/bin/env python3
"""Comparison of the four MarXus methods with each other: identities, agreements and diagnostics.

The four methods (three families) answer different questions, but some of their results must agree exactly
and others only where the time scales separate. This script checks both on the two validation systems, from
the outputs of their run scripts, and writes into each system's directory:

  method_comparison.csv           every compared quantity: system, check, kind, condition, values, deviation
  plots/method_identities.png     exact identities (relative deviation per condition, log scale)
  plots/method_agreement.png      approximate agreements against the eigenvalue separation, and each method vs MESS
  plots/method_diagnostics.png    separation of time scales, sum rule, truncated grains, run time per method

Exact identities (must hold to the printed precision):
  I1  pulse (TimeIntegration at the last time) = final steady state (SteadyStateOlzmann):
      Y_r(inf) = int k_r^T exp(-J t) F dt = k_r^T J^-1 F (README, "Exact identity 1")
  I2  long-time yields from the CSE rate coefficients = final steady state (README, "Exact identity 2";
      Georgievskii et al., J. Phys. Chem. A 117, 12146 (2013), eqs. 21, 25-30)
  I3  one well: CSE k(W -> P) = k_uni of SteadyStateOlzmann (the eigenvector average; Gonzalez-Garcia,
      Olzmann, PCCP 12, 12290 (2010), after eq. 12)
  I4  late-time decay of the pulse = k_uni (the thermal eigenpair; the time integration uses no eigenvector)
  I5  conservation: populations + yields of the pulse = 1 at the last time
  I6  CSE balances: loss balance (G13 eq. 29), capture = R -> wells + R -> products + return (G13 eq. 22)
  I7  sum rule of SteadyStateOlzmann, lambda_1 = k_uni (GO10 eq. 12), within the double-precision floor
Approximate relations (hold where the chemical and relaxation time scales separate):
  A1  prompt R -> P: SteadyStateAbsorbingBarrier k_inf Phi_P against CSE (G13 eq. 21)
  A2  R -> W: SteadyStateAbsorbingBarrier k_inf Phi_stab,W against CSE (G13 eq. 28), and for one well the
      association k_uni K of SteadyStateOlzmann (detailed balance)
  A3  prompt yield + stabilization x thermal fate of each well = final steady state (the stabilized molecules
      are assumed thermalized before they react)
  A4  detailed balance of the CSE pair: k(P -> W)/k(W -> P) against K = k_inf,a/k_inf,d (one well), also for MESS
Each method against MESS, where MESS gives the same quantity.

Run with the science environment:  source ~/.venvs/science/bin/activate && python3 method_comparison.py
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
TORR_PER_ATM = 760.0
NAN = math.nan
ROUNDING_LEVEL = 1e-6  # CSE entries below this fraction of the largest of their row are rounding noise
FLOOR_RATIO = 1e-2  # precision floor eps*max(S_ii) above this fraction of k_uni: lambda_1 not resolved


# ----------------------------------------------------------------------------------------------
# Readers
# ----------------------------------------------------------------------------------------------
def read_tables(path):
    """Tables file of a run (--csv FILE writes FILE_tables.csv): {title: {(T, p_torr): row}}."""
    out, title, header = {}, None, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("# "):
            title, header = line[2:], None
            out[title] = {}
        elif not line.strip() or title is None:
            continue
        elif header is None:
            header = line.split(",")
        else:
            row = {h: (float(v) if v not in ("", "***") else NAN) for h, v in zip(header, line.split(","))}
            out[title][(row["T[K]"], row["P[Torr]"])] = row
    return out


def table(tables, prefix, optional=False):
    """The table whose title starts with `prefix`; {} if `optional` and there is none (e.g. no product channel)."""
    for title, rows in tables.items():
        if title.startswith(prefix):
            return rows
    if optional:
        return {}
    raise KeyError(f"no table '{prefix}'")


def read_evolution(path):
    """Time-evolution blocks of a TimeIntegration machine file: {(T, p_torr): (header, rows)}."""
    out, key, header = {}, None, None
    for line in open(path):
        line = line.rstrip("\n")
        m = re.match(r"# time evolution: T = (\S+) K, p = (\S+) Torr", line)
        if m:
            key, header = (float(m.group(1)), float(m.group(2))), None
            out[key] = (None, [])
        elif line.startswith("#"):
            key = None
        elif not line.strip() or key is None:
            continue
        elif header is None:
            header = line.split(",")
            out[key] = (header, [])
        else:
            out[key][1].append(dict(zip(header, map(float, line.split(",")))))
    return out


def read_thermal(path):
    """Thermal block of a SteadyStateOlzmann machine file: {(T, p_torr): row} (7 significant digits)."""
    out, inside, header = {}, False, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("T[K],P[Torr],k_uni"):
            inside, header = True, line.split(",")
            continue
        if inside:
            if line.startswith("#") or not line.strip():
                break
            row = dict(zip(header, map(float, line.split(","))))
            out[(row["T[K]"], row["P[Torr]"])] = row
    return out


def read_mess(path):
    """MESS species tables: {(T, p_torr): {from: {to: k}}}; '***' -> nan. Pressures in torr or atm."""
    lines = open(path).read().split("\n")
    out, i = {}, 0
    while i < len(lines):
        m = re.match(r"\s*Temperature = (\S+) K\s+Pressure = (\S+) (torr|atm)", lines[i])
        if m:
            p = float(m.group(2)) * (TORR_PER_ATM if m.group(3) == "atm" else 1.0)
            j = i + 1
            while not lines[j].strip().startswith("From\\To"):
                j += 1
            header = []
            for name in lines[j].split()[1:]:  # the escape column of a well carries the well's name
                header.append(f"escape({name})" if name in header else name)
            rows, k = {}, j + 1
            while k < len(lines) and lines[k].split() and lines[k].split()[0] in header:
                fields = lines[k].split()
                rows[fields[0]] = {h: (float(v) if v != "***" else NAN) for h, v in zip(header, fields[1:])}
                k += 1
            out[(float(m.group(1)), p)] = rows
            i = k
        i += 1
    return out


def read_mess_high_pressure(path):
    """MESS high-pressure tables: {T: {from: {to: k}}}."""
    lines = open(path).read().split("\n")
    out, i = {}, 0
    while i < len(lines):
        m = re.match(r"\s*Temperature = (\S+) K\s*$", lines[i])
        if m and i + 2 < len(lines) and "High Pressure Rate Coefficients" in lines[i + 2]:
            j = i + 1
            while not lines[j].strip().startswith("From\\To"):
                j += 1
            header = lines[j].split()[1:]
            rows, k = {}, j + 1
            while k < len(lines) and lines[k].split() and lines[k].split()[0] in header:
                fields = lines[k].split()
                rows[fields[0]] = {h: (float(v) if v != "***" else NAN) for h, v in zip(header, fields[1:])}
                k += 1
            out[float(m.group(1))] = rows
            i = k
        i += 1
    return out


def truncations(out_path):
    """Truncated low-energy grains per T and well, from the 'Collisions:' lines of a report."""
    res = {}
    for line in open(out_path):
        m = re.match(r"\s+T = (\S+) K: (.*)", line)
        if m and "grains (" in line:
            res[float(m.group(1))] = {w: int(g) for w, g in re.findall(r"(\S+) (\d+) grains \(", m.group(2))}
        if line.startswith(" Tunneling:"):
            break
    return res


def run_times(output_dir, stems):
    """Wall time of each run from the modification times of the machine files of the sequential run script:
    the duration of a run is the time since the previous run's file (the first run of the list has none)."""
    times = {s: os.path.getmtime(os.path.join(output_dir, s + ".csv")) for s in stems
             if os.path.exists(os.path.join(output_dir, s + ".csv"))}
    order = sorted(times, key=times.get)
    return {s: times[s] - times[order[k - 1]] for k, s in enumerate(order) if k > 0}


def rel(a, b):
    """Relative deviation a/b - 1 (nan if undefined)."""
    if b == 0 or not (math.isfinite(a) and math.isfinite(b)):
        return NAN
    return a / b - 1.0


# ----------------------------------------------------------------------------------------------
# Collection of the checks
# ----------------------------------------------------------------------------------------------
class Checks:
    def __init__(self, system):
        self.system, self.rows = system, []
        # Conditions where lambda_1 lies below the double-precision floor (floor/k_uni > FLOOR_RATIO): there the
        # identities that involve lambda_1 or the slowest mode are limited by rounding (plotted separately).
        self.below_floor = set()

    def add(self, check, kind, key, quantity, a, b, label_a, label_b, extra=None, deviation=None):
        """`deviation` overrides a/b - 1 where the program prints its own full-precision deviation."""
        self.rows.append({
            "system": self.system, "check": check, "kind": kind, "T_K": key[0], "p_torr": key[1],
            "quantity": quantity, "value_a": a, "value_b": b, "a": label_a, "b": label_b,
            "relative_deviation": rel(a, b) if deviation is None else deviation,
            **({"extra": extra} if extra is not None else {"extra": NAN}),
        })

    def of(self, check):
        return [r for r in self.rows if r["check"] == check]


def late_decay(header, rows, wells):
    """-d ln N/dt of the total well population from the last two output times with 1e-8 < N < 1e-3."""
    pts = [(r["t[s]"], sum(r[f"N({w})"] for w in wells)) for r in rows]
    late = [(t, n) for t, n in pts if 1e-8 < n < 1e-3]
    if len(late) < 2:
        return NAN
    (t1, n1), (t2, n2) = late[-2], late[-1]
    return -math.log(n2 / n1) / (t2 - t1)


def common_checks(c, olz, barrier, cse, ti, evolution, thermal, wells, products, reactant):
    """Identities and relations shared by both systems. `products`: the product columns 'R->X' of the
    bimolecular-to-bimolecular tables."""
    final = table(olz, "final steady state: Yields (% of the formed adducts)")
    final_net = table(olz, "final steady state: Yields without the return to")
    ti_y = table(ti, "time integration: Yields at the last output time")
    ti_left = table(ti, "time integration: Populations left in the wells")
    ti_total = table(ti, "time integration: Total: populations + yields")
    cse_long = table(cse, "CSE: Long-time yields, total")
    cse_diag = table(cse, "CSE: Diagnostics of the CSE solution")
    cse_cap = table(cse, "CSE: Capture, return and net reaction of")
    cse_rw = table(cse, "CSE: Bimolecular-to-well rate coefficients")
    cse_rp = table(cse, "CSE: Bimolecular-to-bimolecular rate coefficients", optional=True)
    bar_rw = table(barrier, "intermediate steady state: Bimolecular-to-well rate coefficients")
    bar_rp = table(barrier, "intermediate steady state: Bimolecular-to-bimolecular rate coefficients", optional=True)
    r_to = reactant + "->"

    for key in sorted(ti_total):
        sep = cse_diag[key]["separation"] if key in cse_diag else NAN
        # I5 conservation of the pulse
        c.add("I5 pulse conservation", "identity", key, "total", ti_total[key]["total"], 100.0, "TimeIntegration", "100%")
        # I1 pulse = final steady state, every exit, where the pulse has ended (populations left < 1e-6 %)
        left = sum(v for k, v in ti_left[key].items() if k.startswith("N("))
        if key in final and left < 1e-6:
            for col in final[key]:
                if col in ("T[K]", "P[Torr]") or not math.isfinite(final[key][col]) or final[key][col] < 1e-4:
                    continue
                c.add("I1 pulse = final steady state", "identity", key, col, ti_y[key][col], final[key][col],
                      "TimeIntegration", "SteadyStateOlzmann")
        # I2 CSE long-time yields = final steady state (% of the net reaction)
        if key in final_net and key in cse_long:
            for col in cse_long[key]:
                if col in ("T[K]", "P[Torr]") or "->" not in col:
                    continue
                product = col.split("->", 1)[1]
                match = [k for k in final_net[key] if k.endswith("->" + product) or k == product]
                if match and math.isfinite(final_net[key][match[0]]) and final_net[key][match[0]] >= 1e-4:
                    c.add("I2 CSE long-time = final steady state", "identity", key, product,
                          cse_long[key][col], final_net[key][match[0]], "CSE", "SteadyStateOlzmann", sep)
        # I4 late-time decay of the pulse = k_uni
        if key in evolution and key in thermal:
            rate = late_decay(*evolution[key], wells)
            if math.isfinite(rate):
                c.add("I4 pulse decay = k_uni", "identity", key, "k_uni", rate, thermal[key]["k_uni[1/s]"],
                      "TimeIntegration", "SteadyStateOlzmann")
        # I6 CSE balances
        if key in cse_diag:
            c.add("I6 CSE loss balance (G13 eq. 29)", "identity", key, "loss balance",
                  1.0 + cse_diag[key]["loss balance"], 1.0, "CSE", "exact", sep)
            net = cse_cap[key]["net"]
            formed = sum(v for k, v in cse_rw[key].items() if k.startswith(r_to)) + \
                sum(v for k, v in cse_rp.get(key, {}).items() if k.startswith(r_to))
            c.add("I6 CSE capture balance (G13 eq. 22)", "identity", key, "net reaction", formed, net, "CSE",
                  "capture - return", sep)
        # I7 sum rule
        if key in thermal:
            th = thermal[key]
            if th["precision_floor[1/s]"] > FLOOR_RATIO * th["k_uni[1/s]"]:
                c.below_floor.add(key)
            c.add("I7 sum rule lambda_1 = k_uni", "identity", key, "lambda_1", th["lambda_1[1/s]"], th["k_uni[1/s]"],
                  "lambda_1", "k_uni", th["precision_floor[1/s]"] / th["k_uni[1/s]"], th["sum_rule_deviation"])
        # A1 prompt R -> P and A2 R -> W: absorbing barrier against CSE
        # CSE entries many orders of magnitude below the largest of their row are rounding noise (and may be
        # negative, as in MESS): only entries above ROUNDING_LEVEL x the row maximum are compared.
        row_max = max([abs(v) for k, v in {**cse_rp.get(key, {}), **cse_rw.get(key, {})}.items() if k.startswith(r_to)] or [0.0])
        if key in bar_rp and key in cse_rp:
            for col in products:
                if col in cse_rp[key] and col in bar_rp[key] and cse_rp[key][col] > ROUNDING_LEVEL * row_max:
                    c.add("A1 prompt R->P: absorbing barrier vs CSE", "approximate", key, col, bar_rp[key][col],
                          cse_rp[key][col], "SteadyStateAbsorbingBarrier", "CSE", sep)
        if key in bar_rw and key in cse_rw:
            total_bar = sum(v for k, v in bar_rw[key].items() if k.startswith(r_to))
            total_cse = sum(v for k, v in cse_rw[key].items() if k.startswith(r_to))
            c.add("A2 total stabilization: absorbing barrier vs CSE", "approximate", key, "sum R->W", total_bar,
                  total_cse, "SteadyStateAbsorbingBarrier", "CSE", sep)
            for col in cse_rw[key]:
                if col.startswith(r_to) and col in bar_rw[key] and cse_rw[key][col] > ROUNDING_LEVEL * row_max:
                    c.add("A2 R->W: absorbing barrier vs CSE", "approximate", key, col, bar_rw[key][col],
                          cse_rw[key][col], "SteadyStateAbsorbingBarrier", "CSE", sep)


# ----------------------------------------------------------------------------------------------
# ZZ-allyl + O2, Case 2 (four wells)
# ----------------------------------------------------------------------------------------------
def case2():
    d = os.path.join(HERE, "ZZAllyl+O2_Gamma_Case2")
    o = os.path.join(d, "marxus_output")
    stem = "case2_tstlevel_E_mess_eckart"
    tabs = {m: read_tables(os.path.join(o, f"{stem}_{m}_tables.csv"))
            for m in ("olzmann", "absorbing_barrier", "cse", "time_integration")}
    evolution = read_evolution(os.path.join(o, f"{stem}_time_integration.csv"))
    thermal = read_thermal(os.path.join(o, f"{stem}_olzmann.csv"))
    wells = ["G2", "G3", "G4", "G6"]
    c = Checks("ZZ-allyl + O2 (Case 2, MESS Eckart model)")
    common_checks(c, tabs["olzmann"], tabs["absorbing_barrier"], tabs["cse"], tabs["time_integration"], evolution,
                  thermal, wells, ["R->P1", "R->P5", "R->P7", "R->escape(G4)"], "R")

    # A3 prompt + stabilization x thermal fate = final steady state (% of the formed adducts)
    prompt = table(tabs["absorbing_barrier"], "intermediate steady state: Yields (% of the formed adducts)")
    final = table(tabs["olzmann"], "final steady state: Yields (% of the formed adducts)")
    fates = {w: table(tabs["olzmann"], f"thermal fates: Thermal fate of the molecules thermalized in {w}") for w in wells}
    sep = table(tabs["cse"], "CSE: Diagnostics of the CSE solution")
    for key in sorted(final):
        for col in final[key]:
            if col in ("T[K]", "P[Torr]") or final[key][col] < 1e-4:
                continue
            total = prompt[key][col] + sum(prompt[key][f"stab({w})"] * fates[w][key][col] / 100.0 for w in wells)
            c.add("A3 prompt + stabilization x fate = final", "approximate", key, col, total, final[key][col],
                  "AbsorbingBarrier + fates", "SteadyStateOlzmann", sep[key]["separation"])

    # Each method against MESS
    mess = read_mess(os.path.join(d, "Gamma-Case2_from_ZZ_allyl+O2_Gamma4-Escape_CCSDT_version-1_2025_09_29_12.7kcal.out"))
    cse_rp = table(tabs["cse"], "CSE: Bimolecular-to-bimolecular rate coefficients")
    cse_rw = table(tabs["cse"], "CSE: Bimolecular-to-well rate coefficients")
    bar_rp = table(tabs["absorbing_barrier"], "intermediate steady state: Bimolecular-to-bimolecular rate coefficients")
    bar_rw = table(tabs["absorbing_barrier"], "intermediate steady state: Bimolecular-to-well rate coefficients")
    olz_rp = table(tabs["olzmann"], "final steady state: Bimolecular-to-bimolecular rate coefficients, overall")
    ti_rp = table(tabs["time_integration"], "time integration: Bimolecular-to-bimolecular rate coefficients, overall")
    for key in sorted(cse_rp):
        m = mess.get(key)
        if not m:
            continue
        for name, rows, label in (("R->P5", cse_rp, "CSE"), ("R->P5", bar_rp, "SteadyStateAbsorbingBarrier"),
                                  ("R->P5", olz_rp, "SteadyStateOlzmann (overall)"),
                                  ("R->P5", ti_rp, "TimeIntegration (overall)"),
                                  ("R->G2", cse_rw, "CSE"), ("R->G2", bar_rw, "SteadyStateAbsorbingBarrier"),
                                  ("R->G4", cse_rw, "CSE"), ("R->G4", bar_rw, "SteadyStateAbsorbingBarrier")):
            if key in rows and name in rows[key]:
                c.add("M vs MESS", "MESS", key, f"{name} {label}", rows[key][name], m["R"][name.split("->")[1]],
                      label, "MESS", sep[key]["separation"])
        for name in ("R->G3", "R->escape(G4)"):
            if cse_rw[key].get(name, cse_rp[key].get(name, NAN)) > ROUNDING_LEVEL * abs(m["R"]["G2"]):
                rows = cse_rw if name in cse_rw[key] else cse_rp
                c.add("M vs MESS", "MESS", key, f"{name} CSE", rows[key][name], m["R"][name.split("->")[1]],
                      "CSE", "MESS", sep[key]["separation"])
    trunc = truncations(os.path.join(o, f"{stem}_olzmann.out"))
    times = run_times(o, [f"{s}_{m}" for m in ("olzmann", "absorbing_barrier", "cse", "time_integration")
                          for s in ("case2_tstlevel_E", stem)])
    return d, c, trunc, times, tabs


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
    thermal_lapack = read_thermal(os.path.join(o, f"{stem}_olzmann_lapack.csv"))
    c = Checks("H + C2H2 = C2H3 (one well)")
    common_checks(c, tabs["olzmann"], tabs["absorbing_barrier"], tabs["cse"], tabs["time_integration"], evolution,
                  thermal, ["W1"], [], "P1")
    cse_w = table(tabs["cse"], "CSE: Rate coefficients from W1")
    cse_rw = table(tabs["cse"], "CSE: Bimolecular-to-well rate coefficients")
    sep = table(tabs["cse"], "CSE: Diagnostics of the CSE solution")
    hp = table(tabs["olzmann"], "thermal rate coefficients: High-pressure rate coefficients k_inf")
    mess = read_mess(os.path.join(d, "reference_mess_output", "c2h3_tight.out"))
    mess_hp = read_mess_high_pressure(os.path.join(d, "reference_mess_output", "c2h3_tight.out"))
    assoc = table(tabs["olzmann"], "thermal rate coefficients: Association by detailed balance")
    assoc_olz = {k: r["P1->W1"] for k, r in assoc.items() if math.isfinite(r["P1->W1"])}
    bar_rw = table(tabs["absorbing_barrier"], "intermediate steady state: Bimolecular-to-well rate coefficients")
    for key in sorted(cse_w):
        s = sep[key]["separation"]
        # I3 one well: CSE k(W1 -> P1) = k_uni
        if key in thermal:
            c.add("I3 one well: CSE k(W->P) = k_uni", "identity", key, "W1->P1", cse_w[key]["W1->P1"],
                  thermal[key]["k_uni[1/s]"], "CSE", "SteadyStateOlzmann", s)
            c.add("E eigen-solvers: k_uni inverse iteration vs LAPACK", "identity", key, "k_uni",
                  thermal[key]["k_uni[1/s]"], thermal_lapack[key]["k_uni[1/s]"], "inverse iteration", "LAPACK", s)
        # A2 association: CSE eq. 28 against SteadyStateOlzmann k_uni K
        if key in assoc_olz:
            c.add("A2 association: CSE vs SteadyStateOlzmann (k_uni K)", "approximate", key, "P1->W1",
                  cse_rw[key]["P1->W1"], assoc_olz[key], "CSE", "SteadyStateOlzmann", s)
            if key in bar_rw and math.isfinite(bar_rw[key].get("P1->W1", NAN)):
                c.add("A2 association: absorbing barrier vs SteadyStateOlzmann", "approximate", key, "P1->W1",
                      bar_rw[key]["P1->W1"], assoc_olz[key], "SteadyStateAbsorbingBarrier", "SteadyStateOlzmann", s)
        # A4 detailed balance of the CSE pair against K = k_inf,a / k_inf,d; also MESS's own pair
        t = key[0]
        if key in assoc_olz and key in thermal:
            k_eq = assoc_olz[key] / thermal[key]["k_uni[1/s]"]  # K = k_inf,a/k_inf,d (k(P1 -> W1) = k_uni K)
            c.add("A4 detailed balance of the pair: CSE", "approximate", key, "k(P1->W1)/k(W1->P1)/K",
                  cse_rw[key]["P1->W1"] / cse_w[key]["W1->P1"], k_eq, "CSE", "K (MarXus)", s)
        if key in mess and t in mess_hp:
            k_eq_mess = mess_hp[t]["P1"]["W1"] / mess_hp[t]["W1"]["P1"]
            c.add("A4 detailed balance of the pair: MESS", "approximate", key, "k(P1->W1)/k(W1->P1)/K",
                  mess[key]["P1"]["W1"] / mess[key]["W1"]["P1"], k_eq_mess, "MESS", "K (MESS)", s)
        # Each method against MESS: association and dissociation
        if key in mess:
            m = mess[key]
            for label, value in (("CSE", cse_rw[key]["P1->W1"]), ("SteadyStateOlzmann (k_uni K)", assoc_olz.get(key, NAN)),
                                 ("SteadyStateAbsorbingBarrier", bar_rw.get(key, {}).get("P1->W1", NAN))):
                c.add("M vs MESS", "MESS", key, f"P1->W1 {label}", value, m["P1"]["W1"], label, "MESS", s)
            for label, value in (("CSE", cse_w[key]["W1->P1"]), ("SteadyStateOlzmann (k_uni)",
                                                                 thermal.get(key, {}).get("k_uni[1/s]", NAN))):
                c.add("M vs MESS", "MESS", key, f"W1->P1 {label}", value, m["W1"]["P1"], label, "MESS", s)
    trunc = truncations(os.path.join(o, f"{stem}_olzmann.out"))
    times = run_times(o, [f"{stem}_{m}" for m in ("olzmann", "absorbing_barrier", "cse", "time_integration")]
                      + [f"{stem}_olzmann_lapack"])
    return d, c, trunc, times, tabs


# ----------------------------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------------------------
def write_csv(d, c):
    path = os.path.join(d, "method_comparison.csv")
    fields = ["system", "check", "kind", "T_K", "p_torr", "quantity", "a", "value_a", "b", "value_b",
              "relative_deviation", "extra", "lambda_1_below_floor"]
    for r in c.rows:
        r["lambda_1_below_floor"] = "yes" if (r["T_K"], r["p_torr"]) in c.below_floor else "no"
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in c.rows:
            w.writerow({k: (("%.8g" % r[k]) if isinstance(r[k], float) else r[k]) for k in fields})
    return path


def summary(c):
    """{check: (count, max |deviation|, condition of the max)} over the finite deviations, without the conditions
    below the double-precision floor for the identities."""
    out = {}
    for r in c.rows:
        dev = abs(r["relative_deviation"])
        if not math.isfinite(dev) or (r["kind"] == "identity" and (r["T_K"], r["p_torr"]) in c.below_floor):
            continue
        n, best, where = out.get(r["check"], (0, -1.0, None))
        if dev > best:
            best, where = dev, (r["T_K"], r["p_torr"], r["quantity"])
        out[r["check"]] = (n + 1, best, where)
    return out


def plot_identities(d, c):
    checks = [k for k in dict.fromkeys(r["check"] for r in c.rows if r["kind"] == "identity")]
    all_devs = [abs(r["relative_deviation"]) for r in c.rows
                if r["kind"] == "identity" and math.isfinite(r["relative_deviation"])]
    x_max = max(1e4, 10 ** math.ceil(math.log10(max(all_devs + [1.0]))) * 1e3)
    fig, ax = plt.subplots(figsize=(12, 0.6 * len(checks) + 2.4))
    for i, check in enumerate(checks):
        rows = [r for r in c.of(check) if math.isfinite(r["relative_deviation"])]
        resolved = [max(abs(r["relative_deviation"]), 1e-17) for r in rows if (r["T_K"], r["p_torr"]) not in c.below_floor]
        floor = [max(abs(r["relative_deviation"]), 1e-17) for r in rows if (r["T_K"], r["p_torr"]) in c.below_floor]
        rng = np.random.default_rng(i)
        line, = ax.semilogx(resolved, i + rng.uniform(-0.25, 0.25, len(resolved)), "o", ms=3.5, alpha=0.6)
        if floor:
            ax.semilogx(floor, i + rng.uniform(-0.25, 0.25, len(floor)), "x", ms=5, color="0.45")
        if resolved:
            text = "equal in all printed digits" if max(resolved) <= 1e-17 else f"max {max(resolved):.1e}"
            ax.text(x_max * 1.5, i, f"{text}  (n = {len(resolved)})", va="center", fontsize=8)
    if c.below_floor:
        ax.plot([], [], "x", color="0.45", label="conditions with lambda_1 below the double-precision floor "
                f"(floor > {FLOOR_RATIO:g} k_uni): rounding-limited, not counted in the max")
        ax.legend(fontsize=7, loc="lower right")
    ax.axvspan(1e-17, 1e-6, color="green", alpha=0.07)
    ax.set_yticks(range(len(checks)))
    ax.set_yticklabels(checks, fontsize=8)
    ax.set_xlim(1e-17, x_max)
    ax.invert_yaxis()
    ax.set_xlabel("|a/b - 1| per condition and quantity (green: below the printed precision, 1e-6)")
    ax.set_title(f"{c.system}: exact identities between the methods", fontsize=10)
    fig.tight_layout(rect=(0, 0, 0.8, 1))
    fig.savefig(os.path.join(d, "plots", "method_identities.png"), dpi=200)
    plt.close(fig)


def plot_agreement(d, c):
    approx = [k for k in dict.fromkeys(r["check"] for r in c.rows if r["kind"] == "approximate")]
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    ax = axes[0]
    for check, marker in zip(approx, "osD^v<>ph"):
        pts = [(r["extra"], abs(r["relative_deviation"])) for r in c.of(check)
               if math.isfinite(r["extra"]) and math.isfinite(r["relative_deviation"])]
        if pts:
            ax.loglog(*zip(*[(x, max(y, 1e-9)) for x, y in pts]), marker, ms=5, mfc="none", label=check)
    ax.set_xlabel("CSE separation of the time scales, $\\Lambda_N/\\Lambda_{N+1}$")
    ax.set_ylabel("|a/b - 1|")
    ax.set_title("approximate relations against the eigenvalue separation", fontsize=10)
    ax.legend(fontsize=7)
    ax = axes[1]
    labels = list(dict.fromkeys(r["quantity"] for r in c.of("M vs MESS")))
    for label, marker in zip(labels, "osD^v<>phx*"):
        pts = sorted((r["T_K"], 100 * r["relative_deviation"]) for r in c.of("M vs MESS")
                     if r["quantity"] == label and math.isfinite(r["relative_deviation"]))
        if pts:
            ts = sorted(set(t for t, _ in pts))
            lo = [min(v for t2, v in pts if t2 == t) for t in ts]
            hi = [max(v for t2, v in pts if t2 == t) for t in ts]
            line, = ax.plot(ts, [(a + b) / 2 for a, b in zip(lo, hi)], "-" + marker, ms=5, mfc="none", label=label)
            ax.fill_between(ts, lo, hi, color=line.get_color(), alpha=0.12)
    ax.axhline(0, color="k", lw=0.8)
    ax.axhspan(-5, 5, color="green", alpha=0.06)
    ax.set_xlabel("T (K) (band: range over the pressures)")
    ax.set_ylabel("MarXus / MESS - 1 (%)")
    ax.set_title("each method against MESS", fontsize=10)
    ax.legend(fontsize=7)
    fig.suptitle(c.system, fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(d, "plots", "method_agreement.png"), dpi=200)
    plt.close(fig)


def plot_diagnostics(d, c, tabs, trunc, times):
    fig, axes = plt.subplots(1, 4, figsize=(19, 4.6))
    diag = table(tabs["cse"], "CSE: Diagnostics of the CSE solution")
    thermal = table(tabs["olzmann"], "thermal rate coefficients: Diagnostics of the thermal eigenpair")
    pressures = sorted({p for _, p in diag})
    for p, color in zip(pressures, plt.cm.plasma(np.linspace(0, 0.85, len(pressures)))):
        pts = sorted((t, diag[(t, q)]["separation"]) for t, q in diag if q == p)
        axes[0].semilogy(*zip(*pts), "-o", ms=4, color=color, label=f"{p:g} Torr")
        pts = sorted((t, thermal[(t, q)]["sum rule"]) for t, q in thermal
                     if q == p and math.isfinite(thermal[(t, q)]["sum rule"]))
        if pts:
            axes[1].semilogy(*zip(*[(t, max(v, 1e-17)) for t, v in pts]), "-o", ms=4, color=color, label=f"{p:g} Torr")
    axes[0].axhline(0.1, color="k", ls="--", lw=1)
    axes[0].set_title("CSE: $\\Lambda_N/\\Lambda_{N+1}$ (warning above 0.1)", fontsize=10)
    axes[1].axhline(1.5e-2, color="k", ls="--", lw=1)
    axes[1].set_title("SteadyStateOlzmann: |$\\lambda_1$ - k$_{uni}$|/k$_{uni}$", fontsize=10)
    for ax in axes[:2]:
        ax.set_xlabel("T (K)")
        ax.legend(fontsize=7)
    ax = axes[2]
    wells = sorted({w for v in trunc.values() for w in v})
    for w, marker in zip(wells, "osD^"):
        pts = sorted((t, v.get(w, 0)) for t, v in trunc.items())
        ax.plot(*zip(*pts), "-" + marker, ms=4, mfc="none", label=w)
    ax.set_xlabel("T (K)")
    ax.set_ylabel("reservoir grains")
    ax.set_title("low-energy reservoir grains (MESMER reservoir state)", fontsize=10)
    ax.legend(fontsize=7)
    ax = axes[3]
    names = sorted(times, key=times.get)
    ax.barh(range(len(names)), [times[n] for n in names], color="tab:blue")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([n.replace("case2_tstlevel_E_", "").replace("c2h3_tight_", "") for n in names], fontsize=7)
    for i, n in enumerate(names):
        ax.text(times[n], i, f" {times[n]:.0f} s", va="center", fontsize=7)
    ax.set_xlabel("wall time on 4 cores (s)")
    ax.set_title("run time per method (all conditions)", fontsize=10)
    fig.suptitle(c.system, fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(d, "plots", "method_diagnostics.png"), dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    for build in (case2, c2h3):
        d, c, trunc, times, tabs = build()
        path = write_csv(d, c)
        plot_identities(d, c)
        plot_agreement(d, c)
        plot_diagnostics(d, c, tabs, trunc, times)
        print(f"\n{c.system}: {path}")
        for check, (n, dev, where) in summary(c).items():
            print(f"  {check:55s} n = {n:4d}  max |dev| = {dev:9.2e}  at {where}")
        print("  run times (s):", ", ".join(f"{k}: {v:.0f}" for k, v in sorted(times.items(), key=lambda x: x[1])))
