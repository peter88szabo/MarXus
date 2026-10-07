#!/usr/bin/env python3
"""Acetyl + O2: MarXus against the stored MESMER 7.1 results.

Reads
  reference_mesmer/*/mesmer.test, *.test       MESMER outputs: Bartis-Widom rate coefficients (loss, isomerization,
                                               irreversible), canonical rate constants (double-double QA run)
  reference_mesmer/*/*.xml                     MESMER inputs (excess-reactant concentration [O2])
  marxus_output/<deck>_<method>_tables.csv     MarXus runs (run_marxus.sh)
and writes
  comparison_table.csv                         every compared quantity: deck, reference, quantity, MESMER, MarXus, deviation
  plots/deviation_bartis_widom.png             CSE rate coefficients against MESMER's Bartis-Widom ones, all decks
  plots/canonical_rates.png                    high-pressure rate coefficients against MESMER's canonical ones
  plots/rates_298K.png                         the rate coefficients themselves (log scale), MESMER x, MarXus circles
  exact_reference.csv                          partition functions and the canonical rate constants without tunneling
                                               evaluated exactly from the XML data (quantum harmonic oscillators from the
                                               zero-point level, classical rigid rotors, spin multiplicity), against both

MESMER's rate coefficients out of the deficient reactant (acetyl) are pseudo-first-order in the excess reactant O2:
k(acetyl -> X) = k(R -> X) [O2], [O2] = me:excessReactantConc of the XML. "loss" entries are the negative diagonal of
MESMER's Kr matrix: for a well its total loss, for acetyl the net reaction (formation of every well and product).

Run with the science environment:  source ~/.venvs/science/bin/activate && python3 compare_with_mesmer.py
"""
import csv
import math
import os
import re
import sys
import xml.etree.ElementTree as ET

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from method_comparison import read_tables, table  # noqa: E402

REF = os.path.join(HERE, "reference_mesmer")
OUT = os.path.join(HERE, "marxus_output")
def marxus_name(mesmer_name):
    """MarXus species name of a MESMER species (make_deck.py): the deficient reactant is R, a pair of sink products
    X+Y is P_X_Y, wells keep their names."""
    if mesmer_name == "acetyl":
        return "R"
    if "+" in mesmer_name:
        return "P_" + mesmer_name.replace("+", "_")
    return mesmer_name

# (deck, label, MESMER test file, MESMER XML)
CASES = [
    ("acetyl_o2", "298 K, 200.72 Torr, classical Eckart (QA run, double-double)",
     "qa_double_double/mesmer.test", "qa_double_double/Acetyl_O2_association.xml"),
    ("acetyl_o2", "298 K, 200.72 Torr, classical Eckart (example run, double)",
     "example_AcetylO2/mesmer.test", "example_AcetylO2/Acetyl_O2_associationEx.xml"),
    ("acetyl_o2_zpe_eckart", "298 K, 200.72 Torr, zero-point Eckart (Tunnelling Ex1, double)",
     "tunnelling_zpe_eckart/Acetyl_O2_associationEx1.test", "tunnelling_zpe_eckart/Acetyl_O2_associationEx1.xml"),
    ("acetyl_o2_250K", "250 K, 37.48 Torr, no tunneling at TS1 (reservoirSink, double)",
     "reservoir_sink/mesmer.test", "reservoir_sink/reservoirSinkAcetylO2.xml"),
]


def mesmer_rates(path):
    """{'X -> Y': value} and {'X loss': value} of the Bartis-Widom section."""
    out, inside = {}, False
    for line in open(path):
        if line.startswith("First order & pseudo first order rate coefficients"):
            inside = True
            continue
        if inside and line.strip() == "}":
            inside = False
            continue
        if inside and "=" in line:
            key, value = line.split("=")
            out[key.strip()] = float(value)
    return out


def mesmer_canonical(path):
    """Canonical (high-pressure) rate constants of the double-double QA output: {(reaction, kind): value}."""
    out = {}
    for line in open(path):
        m = re.match(r"Canonical (.*?) rate constant (?:of|for the) (?:(?:association|isomerization|irreversible) )?"
                     r"(?:reaction )?(?:reverse of reaction )?(R\d) = (\S+)", line)
        if m:
            out[(m.group(2), m.group(1))] = float(m.group(3))
    return out


def excess_concentration(xml_path):
    root = ET.parse(xml_path).getroot()
    for e in root.iter():
        if e.tag.endswith("excessReactantConc"):
            return float(e.text)
    raise ValueError(xml_path)


def marxus_cse(deck):
    """{'X -> Y': k} in MESMER's names from the MarXus CSE tables: wells in s-1, R rows in cm3/s."""
    t = read_tables(os.path.join(OUT, f"{deck}_cse_tables.csv"))
    out = {}
    for well in ("Int1", "Int2"):
        rows = table(t, f"CSE: Rate coefficients from {well}")
        (key, row), = rows.items()
        for col, v in row.items():
            if "->" in col:
                out[col.replace("->", " -> ")] = v
        out[f"{well} loss"] = row[f"{well} loss"]
    for title in ("CSE: Bimolecular-to-well rate coefficients", "CSE: Bimolecular-to-bimolecular rate coefficients"):
        (key, row), = table(t, title).items()
        for col, v in row.items():
            if col.startswith("R->"):
                out[col.replace("->", " -> ")] = v
    (key, row), = table(t, "CSE: Capture, return and net reaction of R").items()
    out["R net"] = row["net"]
    out["R capture"] = row["capture"]
    return out, key


def compare(deck, label, test, xml, rows):
    ref = mesmer_rates(os.path.join(REF, test))
    o2 = excess_concentration(os.path.join(REF, xml))
    ours, key = marxus_cse(deck)
    for name, value in ref.items():
        if name.endswith(" loss"):
            species = name.split()[0]
            mx = ours.get("R net") * o2 if species == "acetyl" else ours.get(f"{species} loss")
            mesmer = -value
        else:
            a, b = [s.strip() for s in name.split("->")]
            mx = ours.get(f"{marxus_name(a)} -> {marxus_name(b)}")
            if mx is not None and a == "acetyl":
                mx *= o2
            mesmer = value
        rows.append({"deck": deck, "reference": label, "quantity": name.replace("HO_2", "HO2").replace("O_2", "O2"),
                     "T_K": key[0], "p_torr": key[1],
                     "mesmer": mesmer, "marxus": mx if mx is not None else math.nan,
                     "deviation_percent": 100 * (mx / mesmer - 1) if mx is not None and mesmer != 0 else math.nan})


KB_CM1 = 0.69503476
H_SI, KB_SI = 6.62607015e-34, 1.380649e-23


def exact_partition_function(m, temperature):
    """q = q_vib q_rot g_e of the MESMER species model, evaluated in closed form: quantum harmonic oscillators with the
    scaled frequencies counted from the zero-point level, classical rigid rotors (asymmetric: sqrt(pi) (kT)^(3/2) /
    (sigma sqrt(ABC)); linear: kT/(sigma B)), the spin multiplicity."""
    kt = KB_CM1 * temperature
    scale = m.get("frequenciesScaleFactor", 1.0)
    q_vib = 1.0
    for f in m.get("vibFreqs", []):
        q_vib /= 1.0 - math.exp(-f * scale / kt)
    b, sigma = m["rotConsts"], m.get("symmetryNumber", 1.0)
    q_rot = kt / (sigma * b[0]) if len(b) == 1 else math.sqrt(math.pi) * kt ** 1.5 / (sigma * math.sqrt(b[0] * b[1] * b[2]))
    return q_vib * q_rot * m.get("spinMultiplicity", 1.0)


def mesmer_qtot(path):
    """MESMER's "Test rovibronic density of states" tables: {(species, T): qtot}."""
    out, species = {}, None
    for line in open(path):
        m = re.match(r"Test rovibronic density of states for: (\S+)", line)
        if m:
            species = m.group(1)
            continue
        f = line.split()
        if species and len(f) == 4 and f[0].replace(".", "").isdigit():
            out[(species, float(f[0]))] = float(f[1])
        elif line.strip() == "}":
            species = None
    return out


def exact_reference(canonical_rows):
    """Exact partition functions against MESMER's qtot, and exact TST (no tunneling) against both codes."""
    from make_deck import CM1_PER_KJMOL, molecules
    mols = molecules(ET.parse(os.path.join(REF, "qa_double_double/Acetyl_O2_association.xml")).getroot())
    rows = []
    for (species, t), q in sorted(mesmer_qtot(os.path.join(REF, "qa_double_double/mesmer.test")).items()):
        exact = exact_partition_function(mols[species], t)
        rows.append({"quantity": f"q({species})", "T_K": t, "exact": exact, "mesmer": q, "marxus": math.nan,
                     "mesmer_vs_exact_percent": 100 * (q / exact - 1), "marxus_vs_exact_percent": math.nan})
    t = 298.0
    for name, ts, well in (("Int2 -> lactone + OH (k_inf)", "TS3", "Int2"), ("Int1 -> ketene + HO2 (k_inf)", "TS2", "Int1")):
        de = (mols[ts]["zpe_kjmol"] - mols[well]["zpe_kjmol"]) * CM1_PER_KJMOL
        exact = KB_SI * t / H_SI * exact_partition_function(mols[ts], t) / exact_partition_function(mols[well], t) \
            * math.exp(-de / (KB_CM1 * t))
        c = next(r for r in canonical_rows if r["quantity"] == name)
        rows.append({"quantity": name, "T_K": t, "exact": exact, "mesmer": c["mesmer"], "marxus": c["marxus"],
                     "mesmer_vs_exact_percent": 100 * (c["mesmer"] / exact - 1),
                     "marxus_vs_exact_percent": 100 * (c["marxus"] / exact - 1)})
    # R1: the input defines k_inf(T) = A (T/T_inf)^n exp(-E_inf/RT); at T = T_inf = 298 K it is A.
    c = next(r for r in canonical_rows if r["quantity"].startswith("R + O2"))
    rows.append({"quantity": c["quantity"], "T_K": t, "exact": 6.0e-12, "mesmer": c["mesmer"], "marxus": c["marxus"],
                 "mesmer_vs_exact_percent": 100 * (c["mesmer"] / 6.0e-12 - 1),
                 "marxus_vs_exact_percent": 100 * (c["marxus"] / 6.0e-12 - 1)})
    with open(os.path.join(HERE, "exact_reference.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: ("%.7g" % v) if isinstance(v, float) else v for k, v in r.items()})
    fig, ax = plt.subplots(figsize=(10, 5))
    rate_rows = [r for r in rows if not r["quantity"].startswith("q(")]
    x = np.arange(len(rate_rows))
    ax.bar(x - 0.2, [r["mesmer_vs_exact_percent"] for r in rate_rows], 0.4, color="0.4", label="MESMER (100 cm$^{-1}$ grains)")
    ax.bar(x + 0.2, [r["marxus_vs_exact_percent"] for r in rate_rows], 0.4, color="tab:blue", label="MarXus (1 cm$^{-1}$ cells)")
    for i, r in enumerate(rate_rows):
        for dx, v in ((-0.2, r["mesmer_vs_exact_percent"]), (0.2, r["marxus_vs_exact_percent"])):
            ax.text(i + dx, v, f"{v:+.2f}%", ha="center", va="bottom" if v >= 0 else "top", fontsize=8)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([r["quantity"] for r in rate_rows], fontsize=8)
    ax.set_ylabel("canonical rate constant / exact - 1 (%)")
    ax.set_title("Acetyl + O$_2$, 298 K: high-pressure rate constants against their exact values\n(TST without tunneling from "
                 "the XML data; R1: k$_\\infty$(298 K) = A of the input)", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "plots", "canonical_rates_vs_exact.png"), dpi=200)
    plt.close(fig)
    return rows


def main():
    os.makedirs(os.path.join(HERE, "plots"), exist_ok=True)
    rows = []
    for deck, label, test, xml in CASES:
        if os.path.exists(os.path.join(OUT, f"{deck}_cse_tables.csv")):
            compare(deck, label, test, xml, rows)
    # Canonical high-pressure rate constants (double-double QA run) against MarXus k_inf (SteadyStateOlzmann run).
    canonical = mesmer_canonical(os.path.join(REF, "qa_double_double/mesmer.test"))
    t = read_tables(os.path.join(OUT, "acetyl_o2_olzmann_tables.csv"))
    (key, kinf), = table(t, "thermal rate coefficients: High-pressure rate coefficients k_inf").items()
    (_, cap), = table(read_tables(os.path.join(OUT, "acetyl_o2_cse_tables.csv")), "CSE: Capture, return and net reaction of R").items()
    pairs = [(("R1", "bimolecular"), "R + O2 -> Int1 (k_inf, cm3/s)", cap["capture"]),
             (("R1", "first order"), "Int1 -> R (k_inf)", kinf.get("Int1->R")),
             (("R2", "first order forward"), "Int1 -> Int2 (k_inf)", kinf.get("Int1->Int2")),
             (("R2", "first order backward"), "Int2 -> Int1 (k_inf)", kinf.get("Int2->Int1")),
             (("R3", "pseudo first order forward"), "Int2 -> lactone + OH (k_inf)", kinf.get("Int2->P_lactone_OH")),
             (("R4", "pseudo first order forward"), "Int1 -> ketene + HO2 (k_inf)", kinf.get("Int1->P_ketene_HO2"))]
    canonical_rows = []
    for ref_key, name, mx in pairs:
        ref = canonical.get(ref_key)
        if ref is None or mx is None:
            print("missing", ref_key, name, ref, mx)
            continue
        canonical_rows.append({"deck": "acetyl_o2", "reference": "canonical rate constants (QA run, double-double)",
                               "quantity": name, "T_K": key[0], "p_torr": "inf", "mesmer": ref, "marxus": mx,
                               "deviation_percent": 100 * (mx / ref - 1)})
    with open(os.path.join(HERE, "comparison_table.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows + canonical_rows:
            w.writerow({k: ("%.6g" % v) if isinstance(v, float) else v for k, v in r.items()})

    # Deviations of the phenomenological rate coefficients, one row of bars per reference run.
    quantities = list(dict.fromkeys(r["quantity"] for r in rows))
    labels = list(dict.fromkeys(r["reference"] for r in rows))
    fig, ax = plt.subplots(figsize=(15, 6.5))
    width = 0.8 / len(labels)
    for i, label in enumerate(labels):
        vals = [next((r["deviation_percent"] for r in rows if r["reference"] == label and r["quantity"] == q), np.nan)
                for q in quantities]
        ax.bar(np.arange(len(quantities)) + (i - (len(labels) - 1) / 2) * width, vals, width, label=label)
    ax.axhline(0, color="k", lw=0.8)
    ax.axhspan(-5, 5, color="green", alpha=0.07)
    ax.set_xticks(range(len(quantities)))
    ax.set_xticklabels(quantities, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("MarXus (CSE) / MESMER (Bartis-Widom) - 1 (%)")
    ax.set_title("Acetyl + O$_2$: phenomenological rate coefficients against MESMER (green band: $\\pm$5%)", fontsize=11)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "plots", "deviation_bartis_widom.png"), dpi=200)
    plt.close(fig)

    # The rate coefficients themselves at 298 K (QA run).
    qa = [r for r in rows if r["reference"].startswith("298 K, 200.72 Torr, classical Eckart (QA")]
    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(qa))
    ax.semilogy(x, [abs(r["mesmer"]) for r in qa], "x", ms=12, mew=2, color="k", label="MESMER (double-double)")
    ax.semilogy(x, [abs(r["marxus"]) for r in qa], "o", ms=6, mfc="none", color="tab:blue", label="MarXus CSE")
    ax.set_xticks(x)
    ax.set_xticklabels([r["quantity"] for r in qa], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("rate coefficient (s$^{-1}$; acetyl rows pseudo-first-order at [O$_2$])")
    ax.set_title("Acetyl + O$_2$, 298 K, 200.72 Torr He", fontsize=11)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "plots", "rates_298K.png"), dpi=200)
    plt.close(fig)

    # Canonical rate constants.
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(range(len(canonical_rows)), [r["deviation_percent"] for r in canonical_rows], color="tab:blue")
    for i, r in enumerate(canonical_rows):
        ax.text(i, r["deviation_percent"], f"{r['deviation_percent']:+.2f}%", ha="center",
                va="bottom" if r["deviation_percent"] >= 0 else "top", fontsize=8)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(range(len(canonical_rows)))
    ax.set_xticklabels([r["quantity"] for r in canonical_rows], rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("MarXus k$_\\infty$ / MESMER canonical - 1 (%)")
    ax.set_title("Acetyl + O$_2$, 298 K: high-pressure rate coefficients against MESMER's canonical rate constants",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "plots", "canonical_rates.png"), dpi=200)
    plt.close(fig)

    for r in exact_reference(canonical_rows):
        print(f"exact {r['quantity']:30s} T={r['T_K']:5.0f}  exact {r['exact']:12.6e}  MESMER {r['mesmer_vs_exact_percent']:+8.3f}%"
              f"  MarXus {r['marxus_vs_exact_percent']:+8.3f}%")
    for r in rows + canonical_rows:
        print(f"{r['deck']:22s} {r['reference'][:40]:40s} {r['quantity']:26s} MESMER {r['mesmer']:12.5e}  "
              f"MarXus {r['marxus']:12.5e}  {r['deviation_percent']:+8.3f}%")


if __name__ == "__main__":
    main()
