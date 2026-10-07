"""Variant with excited electronic test levels (W2: 500 1/cm, g = 2; HO2: 7029 1/cm, g = 2) against MESS.

Compares, with the MESS Eckart model, the high-pressure rate coefficients and the CSE rate coefficients of
  reference_mess_electronic_levels/c2h4_ho2_electronic_levels.out  and  marxus_output/electronic_levels_*,
and the change from the base deck (reference_mess/, marxus_output/mess_eckart_*) in each code.
Writes electronic_levels_comparison.csv and prints a summary (electronic_levels_summary.txt).
"""
import csv, os, re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BAR_TO_TORR = 750.0616827041697


def read_mess_out(path):
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
            f = row.split()
            table[f[0]] = {h: (np.nan if v == "***" else float(v)) for h, v in zip(header, f[1:])}
        if m_tp:
            pressure[(float(m_tp.group(1)), float(m_tp.group(2)))] = table
        else:
            high[float(m_t.group(1))] = table
    return high, pressure


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


def block(blocks, prefix):
    for title, rows in blocks.items():
        if title.startswith(prefix):
            return {(r["T[K]"], round(r["P[Torr]"], 6)): r for r in rows}
    raise KeyError(prefix)


def marxus(stem):
    olz = read_tables(os.path.join(HERE, "marxus_output", f"{stem}_olzmann_tables.csv"))
    cse = read_tables(os.path.join(HERE, "marxus_output", f"{stem}_cse_tables.csv"))
    return {"kinf": block(olz, "thermal rate coefficients: High-pressure"), "cap": block(cse, "CSE: Capture, return"),
            "w2": block(cse, "CSE: Rate coefficients from W2"), "p1w2": block(cse, "CSE: Bimolecular-to-well"),
            "p1p2": block(cse, "CSE: Bimolecular-to-bimolecular rate coefficients")}


mess_variant = read_mess_out(os.path.join(HERE, "reference_mess_electronic_levels", "c2h4_ho2_electronic_levels.out"))
mess_base = read_mess_out(os.path.join(HERE, "reference_mess", "c2h4_ho2.out"))
mx_variant, mx_base = marxus("electronic_levels"), marxus("mess_eckart")

rows = []
for t in sorted(mess_variant[0]):
    first = lambda d: next(r for (tt, _), r in sorted(d.items()) if tt == t)
    for name, mv, mb, xv, xb in [
        ("k_inf W2->P1", mess_variant[0][t]["W2"]["P1"], mess_base[0][t]["W2"]["P1"], first(mx_variant["kinf"])["W2->P1"], first(mx_base["kinf"])["W2->P1"]),
        ("k_inf W2->P2", mess_variant[0][t]["W2"]["P2"], mess_base[0][t]["W2"]["P2"], first(mx_variant["kinf"])["W2->P2"], first(mx_base["kinf"])["W2->P2"]),
        ("k_inf P1->W2", mess_variant[0][t]["P1"]["W2"], mess_base[0][t]["P1"]["W2"], first(mx_variant["cap"])["capture"], first(mx_base["cap"])["capture"]),
    ]:
        rows.append({"T_K": t, "p_bar": np.nan, "quantity": name, "mess": mv, "marxus": xv, "dev_percent": 100 * (xv / mv - 1),
                     "mess_change_percent": 100 * (mv / mb - 1), "marxus_change_percent": 100 * (xv / xb - 1)})
for (t, p), mv in sorted(mess_variant[1].items()):
    key = (t, round(p * BAR_TO_TORR, 6))
    mb = mess_base[1][(t, p)]
    kept = not np.isnan(mv["W2"]["P1"])
    xv_kept = not np.isnan(mx_variant["w2"][key]["W2->P1"])
    if kept != xv_kept:
        rows.append({"T_K": t, "p_bar": p, "quantity": "merging differs", "mess": np.nan, "marxus": np.nan, "dev_percent": np.nan,
                     "mess_change_percent": np.nan, "marxus_change_percent": np.nan})
        continue
    net_m = mv["P1"]["P1"] if kept else mv["P1"]["P2"]
    net_mb = mb["P1"]["P1"] if not np.isnan(mb["W2"]["P1"]) else mb["P1"]["P2"]
    quantities = [("net P1->P2", net_m, net_mb, mx_variant["cap"][key]["net"], mx_base["cap"][key]["net"])]
    if kept:
        quantities += [("W2->P1", mv["W2"]["P1"], mb["W2"]["P1"], mx_variant["w2"][key]["W2->P1"], mx_base["w2"][key]["W2->P1"]),
                       ("W2->P2", mv["W2"]["P2"], mb["W2"]["P2"], mx_variant["w2"][key]["W2->P2"], mx_base["w2"][key]["W2->P2"]),
                       ("P1->W2", mv["P1"]["W2"], mb["P1"]["W2"], mx_variant["p1w2"][key]["P1->W2"], mx_base["p1w2"][key]["P1->W2"])]
    for name, m, m0, x, x0 in quantities:
        rows.append({"T_K": t, "p_bar": p, "quantity": name, "mess": m, "marxus": x, "dev_percent": 100 * (x / m - 1),
                     "mess_change_percent": 100 * (m / m0 - 1), "marxus_change_percent": 100 * (x / x0 - 1)})
with open(os.path.join(HERE, "electronic_levels_comparison.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow({k: (v if isinstance(v, str) else "%.6g" % v) for k, v in r.items()})
for r in rows:
    if r["quantity"] == "merging differs":
        print(f"{r['T_K']:6.0f} {r['p_bar']:7.2f}  merging differs")
    else:
        print(f"{r['T_K']:6.0f} {r['p_bar']:7.2f}  {r['quantity']:<13} MarXus/MESS {r['dev_percent']:+7.3f}%   "
              f"change from base: MESS {r['mess_change_percent']:+8.3f}%  MarXus {r['marxus_change_percent']:+8.3f}%")
