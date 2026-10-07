#!/usr/bin/env python3
"""Write the MarXus input decks (MESS format) of the acetyl + O2 system from the MESMER XML inputs.

Reads (unchanged copies of the MESMER 7.1 distribution, reference_mesmer/):
  example_AcetylO2/Acetyl_O2_associationEx.xml        298 K, 200.72 Torr He; Eckart at TS1 with classical barriers
  tunnelling_zpe_eckart/Acetyl_O2_associationEx1.xml  the same with Eckart barriers from the zero-point levels
  reservoir_sink/reservoirSinkAcetylO2.xml            250 K, 37.48 Torr; no tunneling at TS1
and writes input/acetyl_o2*.inp.

Model of the XML, as MESMER evaluates it:
  - vibrational frequencies x me:frequenciesScaleFactor (0.9854), classical rigid rotors from me:rotConsts with
    me:symmetryNumber, electronic degeneracy = me:spinMultiplicity, energies me:ZPE (kJ/mol, zero-point levels);
  - R1 (acetyl + O2 -> Int1): inverse Laplace transform of k_inf(T) = A (T/T_inf)^n exp(-E_inf/RT);
  - R2 (Int1 -> Int2, TS1): RRKM with Eckart tunneling; the imaginary frequency is not scaled. Barrier heights:
    classical by default, V = E_c(TS) - E_c(well) with E_c = ZPE - 1/2 sum(unscaled frequencies); with me:useZPE
    the zero-point barriers; or the explicit me:BarrierHeights;
  - R3, R4: RRKM to sink products (no product properties are needed: the products are `Dummy` species);
  - exponential down <dE_down> = 130 cm-1 (temperature exponent 0), Lennard-Jones with He (combining rules
    sigma = (s1 + s2)/2, eps = sqrt(e1 e2); Neufeld collision integral);
  - grain 100 cm-1, grid top = highest barrier + 25 kT.
Deck extensions used (MarXus): RotationalConstants[1/cm] and Mass[amu] in place of a geometry, an
InverseLaplaceTransform block in the association barrier, `Dummy` products (as in MESS).

Run with the science environment:  source ~/.venvs/science/bin/activate && python3 make_deck.py
"""
import math
import os
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "reference_mesmer")
KB_CM1 = 0.69503476           # cm-1/K (MarXus constants.rs, KB_CM)
CM1_PER_KJMOL = 83.593543619  # 1 kJ/mol in cm-1 as MESMER converts it: 9864.0381470825287 cm-1 / 118.0 kJ/mol,
                              # the zero-point barrier of TS1 in Tunnelling Ex2 (me:BarrierHeights) over its ZPE difference
K_TO_CM1 = KB_CM1             # Lennard-Jones epsilon: K -> cm-1


def strip(tag):
    return tag.split("}", 1)[1] if "}" in tag else tag


def child(node, name):
    for c in node:
        if strip(c.tag) == name:
            return c
    return None


def children(node, name):
    return [c for c in node if strip(c.tag) == name]


def molecules(root):
    """{id: {property: value}} with ZPE (kJ/mol), rotConsts, symmetryNumber, vibFreqs, scale, MW, spin, epsilon,
    sigma, imFreqs and the energy-transfer <dE_down> (cm-1)."""
    out = {}
    for mol in root.iter():
        if strip(mol.tag) != "molecule" or "id" not in mol.attrib or child(mol, "propertyList") is None:
            continue
        props = {}
        for prop in children(child(mol, "propertyList"), "property"):
            key = prop.attrib.get("dictRef", "").replace("me:", "")
            node = child(prop, "scalar") if child(prop, "scalar") is not None else child(prop, "array")
            if node is None or node.text is None:
                continue
            units = node.attrib.get("units")
            values = [float(x) for x in node.text.split()]
            props[key] = (values, units)
        mol_props = {"id": mol.attrib["id"], "description": mol.attrib.get("description", "")}
        for key, (values, units) in props.items():
            if key == "ZPE":
                assert units in ("kJ/mol", None), (mol.attrib["id"], units)
                mol_props["zpe_kjmol"] = values[0]
            elif key in ("rotConsts", "vibFreqs", "imFreqs"):
                assert units in ("cm-1", None), (mol.attrib["id"], key, units)
                mol_props[key] = values
            elif key in ("symmetryNumber", "frequenciesScaleFactor", "MW", "spinMultiplicity", "epsilon", "sigma"):
                mol_props[key] = values[0]
        etm = None
        for e in mol.iter():
            if strip(e.tag) == "deltaEDown":
                assert e.attrib.get("units", "cm-1") == "cm-1"
                etm = float(e.text)
        if etm is not None:
            mol_props["deltaEDown"] = etm
        for e in mol.iter():
            if strip(e.tag) == "deltaEDownTExponent":
                raise SystemExit("deltaEDownTExponent is not handled")
        out[mol.attrib["id"]] = mol_props
    return out


def reactions(root):
    """Active reactions: id, reactants [(ref, role)], products, transition state, method, ILT, tunneling."""
    out = []
    for r in root.iter():
        if strip(r.tag) != "reaction" or r.attrib.get("active", "true") == "false":
            continue
        rx = {"id": r.attrib["id"], "reactants": [], "products": [], "ts": None, "ilt": None, "tunneling": None}
        for side in ("reactant", "product"):
            for s in children(r, side):
                m = child(s, "molecule")
                rx[side + "s"].append((m.attrib["ref"], m.attrib.get("role")))
        ts = child(r, "transitionState")
        if ts is not None:
            rx["ts"] = child(ts, "molecule").attrib["ref"]
        method = child(r, "MCRCMethod")
        kind = method.attrib.get("name") or method.attrib.get("{http://www.w3.org/2001/XMLSchema-instance}type", "")
        if "ILT" in kind:
            def value(name, parent=r):
                node = child(method, name) if child(method, name) is not None else child(parent, name)
                return node
            pre = value("preExponential")
            ea = value("activationEnergy")
            assert pre.attrib.get("units", "cm3molecule-1s-1") == "cm3molecule-1s-1", pre.attrib
            assert float(ea.text) == 0.0 or ea.attrib.get("units") == "kJ/mol", ea.attrib
            rx["ilt"] = {"A": float(pre.text), "Ea_kjmol": float(ea.text), "T_inf": float(value("TInfinity").text),
                         "n": float(value("nInfinity").text)}
        tun = child(r, "tunneling")
        if tun is not None:
            heights = child(tun, "BarrierHeights")
            mode = "classical"
            if child(tun, "useZPE") is not None:
                mode = "zpe"
            if heights is not None:
                assert heights.attrib.get("units") == "cm-1"
                mode = ("explicit", float(heights.attrib["V0"]), float(heights.attrib["V1"]))
            rx["tunneling"] = mode
        conc = child(r, "excessReactantConc")
        if conc is not None:
            rx["excess_concentration"] = float(conc.text)
        out.append(rx)
    return out


def conditions(root):
    cond = next(e for e in root.iter() if strip(e.tag) == "conditions")
    pairs = [(float(p.attrib["T"]), float(p.attrib["P"]), p.attrib.get("units")) for p in cond.iter()
             if strip(p.tag) == "PTpair"]
    assert all(u == "Torr" for _, _, u in pairs)
    bath = next(e.text.strip() for e in cond.iter() if strip(e.tag) == "bathGas")
    model = next(e for e in root.iter() if strip(e.tag) == "modelParameters")
    grain = next(e for e in model.iter() if strip(e.tag) == "grainSize")
    assert grain.attrib.get("units") == "cm-1"
    above = next(e for e in model.iter() if strip(e.tag) == "energyAboveTheTopHill")
    return pairs, bath, float(grain.text), float(above.text)


def classical_energy_cm1(mol):
    """E_c = ZPE - 1/2 sum of the unscaled frequencies (cm-1)."""
    return mol["zpe_kjmol"] * CM1_PER_KJMOL - 0.5 * sum(mol.get("vibFreqs", []))


def rrho(mol, indent, ref_kjmol, scale_note=True):
    """RRHO block of a molecule: rotational constants, mass, symmetry, scaled frequencies, energy, spin."""
    pad = " " * indent
    scale = mol.get("frequenciesScaleFactor", 1.0)
    freqs = [f * scale for f in mol.get("vibFreqs", [])]
    rot = mol["rotConsts"]
    lines = [f"{pad}RRHO",
             f"{pad}  RotationalConstants[1/cm]      {len(rot)}",
             f"{pad}    " + "  ".join(f"{b:.6g}" for b in rot),
             f"{pad}  Mass[amu]                      {mol['MW']:g}",
             f"{pad}  Core RigidRotor",
             f"{pad}    SymmetryFactor               {mol.get('symmetryNumber', 1.0):g}",
             f"{pad}  End"]
    if freqs:
        note = f"   # x {scale:g}" if scale_note and scale != 1.0 else ""
        lines.append(f"{pad}  Frequencies[1/cm]              {len(freqs)}{note}")
        for i in range(0, len(freqs), 6):
            lines.append(f"{pad}    " + "  ".join(f"{f:.4f}" for f in freqs[i:i + 6]))
    return lines, freqs


def tail(mol, indent, energy_kjmol):
    pad = " " * indent
    return [f"{pad}  ZeroEnergy[kJ/mol]             {energy_kjmol:.4f}",
            f"{pad}  ElectronicLevels[1/cm]         1",
            f"{pad}    0   {mol.get('spinMultiplicity', 1.0):g}",
            f"{pad}End"]


def write_deck(xml_path, deck_name, title):
    root = ET.parse(xml_path).getroot()
    mols = molecules(root)
    rxs = reactions(root)
    pairs, bath, grain_cm1, above_kt = conditions(root)
    source = next(r for r in rxs if r["ilt"] is not None)
    deficient = next(ref for ref, role in source["reactants"] if role == "deficientReactant")
    excess = next(ref for ref, role in source["reactants"] if role == "excessReactant")
    ref_kjmol = mols[deficient]["zpe_kjmol"] + mols[excess]["zpe_kjmol"]
    wells = sorted({ref for r in rxs for ref, role in r["reactants"] + r["products"] if role == "modelled"},
                   key=lambda w: [ref for r in rxs for ref, _ in r["reactants"] + r["products"]].index(w))
    temperatures = sorted({t for t, _, _ in pairs})
    pressures = sorted({p for _, p, _ in pairs})
    kt_min = KB_CM1 * min(temperatures)
    step = grain_cm1 / kt_min * (1 + 1e-6)          # whole 1 cm-1 cells: 100 cells
    well_mol = mols[wells[0]]
    for w in wells[1:]:
        for key in ("epsilon", "sigma", "deltaEDown"):
            assert mols[w][key] == well_mol[key], (w, key)
    he = mols[bath]
    out = [f"# {title}",
           f"# Written by make_deck.py from {os.path.relpath(xml_path, HERE)} (MESMER 7.1 distribution). Energies in",
           f"# kJ/mol relative to {deficient} + {excess} (zero-point levels); frequencies scaled by the XML factor.",
           "TemperatureList[K]                " + "  ".join(f"{t:g}" for t in temperatures),
           "PressureList[torr]                " + "  ".join(f"{p:g}" for p in pressures),
           f"EnergyStepOverTemperature         {step:.8f}      # {grain_cm1:g} cm-1 at {min(temperatures):g} K",
           f"ExcessEnergyOverTemperature       {above_kt:g}",
           "ModelEnergyLimit[kcal/mol]        400",
           "CalculationMethod                 direct",
           "WellCutoff                        10",
           "ChemicalEigenvalueMax             0.2",
           "MarXus",
           "  CollisionIntegral               Neufeld",
           "End",
           "Model",
           "  EnergyRelaxation",
           "    Exponential",
           f"      Factor[1/cm]                {well_mol['deltaEDown']:g}",
           "      Power                       0",
           "      ExponentCutoff              15",
           "    End",
           "  CollisionFrequency",
           "    LennardJones",
           f"      Epsilons[1/cm]              {he['epsilon'] * K_TO_CM1:.4f}  {well_mol['epsilon'] * K_TO_CM1:.4f}   # {he['epsilon']:g} K, {well_mol['epsilon']:g} K",
           f"      Sigmas[angstrom]            {he['sigma']:g}  {well_mol['sigma']:g}",
           f"      Masses[amu]                 {he['MW']:g}  {well_mol['MW']:g}",
           "    End"]
    for w in wells:
        m = mols[w]
        block, _ = rrho(m, 6, ref_kjmol)
        out += [f"  Well     {w}                 # {m['description']}", "    Species"] + block
        out += tail(m, 6, m["zpe_kjmol"] - ref_kjmol)
        out += ["  End"]
    # The bimolecular source and the sink products.
    out += [f"  Bimolecular  R                  # {deficient} + {excess}"]
    for frag in (deficient, excess):
        m = mols[frag]
        block, _ = rrho(m, 6, ref_kjmol)
        out += [f"    Fragment   {frag}"] + block + tail(m, 6, 0.0)
    out += ["    GroundEnergy[kJ/mol]           0.0", "  End"]
    product_names = {}
    for r in rxs:
        if all(role == "sink" for _, role in r["products"]):
            name = "P_" + "_".join(ref for ref, _ in r["products"])
            product_names[r["id"]] = name
            energy = sum(mols[ref]["zpe_kjmol"] for ref, _ in r["products"]) - ref_kjmol
            out += [f"  Bimolecular  {name}       # sink; asymptote at {energy:.2f} kJ/mol (not used)", "    Dummy"]
    # Barriers.
    for r in rxs:
        left = next(ref for ref, role in r["reactants"] if role == "modelled") if r["ilt"] is None else "R"
        right = (next(ref for ref, role in r["products"] if role == "modelled") if r["id"] not in product_names
                 else product_names[r["id"]])
        if r["ilt"] is not None:
            ilt = r["ilt"]
            out += [f"  Barrier  B_{r['id']}  {right}  R        # {deficient} + {excess} -> {right}: inverse Laplace transform",
                    "    RRHO",
                    "      InverseLaplaceTransform",
                    "        Direction                  Association",
                    f"        PreExponential[cm^3/s]     {ilt['A']:g}",
                    f"        TemperatureExponent        {ilt['n']:g}",
                    f"        ReferenceTemperature[K]    {ilt['T_inf']:g}",
                    f"        ActivationEnergy[kJ/mol]   {ilt['Ea_kjmol']:g}",
                    "      End",
                    "      ZeroEnergy[kJ/mol]           0.0        # not used: the ILT threshold is the asymptote + E_inf",
                    "    End"]
            continue
        ts = mols[r["ts"]]
        block, _ = rrho(ts, 4, ref_kjmol)
        lines = [f"  Barrier  B_{r['id']}  {left}  {right}        # {ts['description']}"] + block
        if r["tunneling"] is not None:
            reactant = mols[left]
            product = mols[next(ref for ref, role in r["products"] if role == "modelled")]
            mode = r["tunneling"]
            if mode == "classical":
                v0 = classical_energy_cm1(ts) - classical_energy_cm1(reactant)
                v1 = classical_energy_cm1(ts) - classical_energy_cm1(product)
                note = "classical barriers, E_c = ZPE - 1/2 sum(unscaled frequencies)"
            elif mode == "zpe":
                v0 = (ts["zpe_kjmol"] - reactant["zpe_kjmol"]) * CM1_PER_KJMOL
                v1 = (ts["zpe_kjmol"] - product["zpe_kjmol"]) * CM1_PER_KJMOL
                note = "zero-point barriers"
            else:
                v0, v1 = mode[1], mode[2]
                note = "explicit barrier heights"
            lines += ["      Tunneling    Eckart                 # " + note,
                      f"        ImaginaryFrequency[1/cm]   {ts['imFreqs'][0]:g}     # not scaled",
                      f"        WellDepth[1/cm]            {v0:.4f}",
                      f"        WellDepth[1/cm]            {v1:.4f}",
                      "      End"]
        out += lines + tail(ts, 4, ts["zpe_kjmol"] - ref_kjmol)
    out += ["End", ""]
    os.makedirs(os.path.join(HERE, "input"), exist_ok=True)
    path = os.path.join(HERE, "input", deck_name)
    with open(path, "w") as f:
        f.write("\n".join(out))
    concentration = source.get("excess_concentration")
    print(f"{path}: wells {wells}, T {temperatures}, p {pressures} Torr, [{excess}] = {concentration}")
    return path


if __name__ == "__main__":
    write_deck(os.path.join(REF, "example_AcetylO2", "Acetyl_O2_associationEx.xml"), "acetyl_o2.inp",
               "Acetyl + O2, MESMER example AcetylO2: 298 K, 200.72 Torr He, Eckart at TS1 (classical barriers)")
    write_deck(os.path.join(REF, "tunnelling_zpe_eckart", "Acetyl_O2_associationEx1.xml"), "acetyl_o2_zpe_eckart.inp",
               "Acetyl + O2, MESMER example Tunnelling Ex1: Eckart at TS1 with zero-point barriers")
    write_deck(os.path.join(REF, "reservoir_sink", "reservoirSinkAcetylO2.xml"), "acetyl_o2_250K.inp",
               "Acetyl + O2, MESMER example reservoirSink: 250 K, 37.48 Torr He, no tunneling at TS1")
