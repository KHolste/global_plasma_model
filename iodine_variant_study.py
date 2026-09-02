#!/usr/bin/env python3
"""
iodine_variant_study.py -- Vergleichsstudie der I2-MTCS und Diss.Attachment Varianten.

Nutzt die expliziten Datenvarianten aus chemistry/iodine_lafleur_v1/variants/.
"""
from __future__ import annotations
import sys, copy, json, math
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plasma_chemistry import load_chemistry, ChemistryPackage, ThrusterGeometry, RateModel, E_CH, KB, CONV
from generic_solver import solve_steady_state

BASE = Path(__file__).resolve().parent / "chemistry" / "iodine_lafleur_v1"
GEOM = ThrusterGeometry(
    R=0.06, L=0.10, betai=0.7, betag=0.3, Vgrid=1000.0, sgrid=0.0015,
    frequency=13.56e6, Nw=5, R_ohm=2.0, Rc=0.07, lc=0.10,
)
P_ABS = 400.0
Q0 = 2.1e18


def load_variant_table(path: Path):
    """Lade Ratentabelle aus CSV."""
    Te, K = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("Te"):
                continue
            parts = line.strip().split(",")
            if len(parts) >= 2:
                Te.append(float(parts[0]))
                K.append(float(parts[1]))
    return np.array(Te), np.array(K)


def apply_variant(chem: ChemistryPackage, rxn_id: str, variant_csv: Path) -> ChemistryPackage:
    """Ersetze die Ratentabelle einer Reaktion durch eine Variante."""
    c = copy.deepcopy(chem)
    Te, K = load_variant_table(variant_csv)
    for rxn in c.reactions:
        if rxn.id == rxn_id:
            rxn.rate._table_Te = Te
            rxn.rate._table_K = K
            break
    return c


def run(chem, label=""):
    r = solve_steady_state(chem, GEOM, P_ABS, Q0, max_iter=1200, tol=0.5)
    if not r["converged"]:
        return None
    d = r["densities"]
    nI, nI2 = d.get("I", 0), d.get("I2", 0)
    r["diss_frac"] = nI / (nI + 2*nI2) if (nI + 2*nI2) > 0 else 0
    r["alpha"] = d.get("I-", 0) / r["ne"] if r["ne"] > 0 else 0
    r["label"] = label
    return r


def main():
    print("=" * 90)
    print("  Iodine Variant Comparison Study")
    print(f"  P_abs={P_ABS}W  Q0={Q0:.1e} pps  Geometry: R=6cm L=10cm f=13.56MHz")
    print("=" * 90)

    chem = load_chemistry(BASE / "chemistry.json")
    variants_meta = json.loads((BASE / "variants" / "variants.json").read_text(encoding="utf-8"))

    # Nominal baseline
    nom = run(chem, "NOMINAL")
    if not nom:
        print("Nominal nicht konvergiert!")
        return

    # Collect results
    all_results = [nom]

    # Run all Kelm_I2 variants
    print("\n--- Kelm_I2 Varianten ---")
    for vname, vinfo in variants_meta["Kelm_I2"]["variants"].items():
        vpath = BASE / "variants" / vinfo["file"]
        if not vpath.exists():
            print(f"  {vname}: Datei nicht gefunden")
            continue
        c = apply_variant(chem, "el_I2", vpath)
        r = run(c, f"MTCS:{vname}")
        if r:
            all_results.append(r)
        else:
            print(f"  MTCS:{vname}: nicht konvergiert")

    # Run all Kdissatt_I2 variants
    print("\n--- Kdissatt_I2 Varianten ---")
    for vname, vinfo in variants_meta["Kdissatt_I2"]["variants"].items():
        vpath = BASE / "variants" / vinfo["file"]
        if not vpath.exists():
            print(f"  {vname}: Datei nicht gefunden")
            continue
        c = apply_variant(chem, "dissatt_I2", vpath)
        r = run(c, f"ATT:{vname}")
        if r:
            all_results.append(r)
        else:
            print(f"  ATT:{vname}: nicht konvergiert")

    # Results table
    print("\n" + "=" * 90)
    print(f"{'Variante':>25} {'Te':>6} {'Tg':>5} {'ne':>10} {'n_I':>10} {'n_I2':>10} {'n_I-':>10} {'diss':>6} {'alpha':>7} {'nu_m':>10}")
    print("-" * 90)
    for r in all_results:
        d = r["densities"]
        print(f"{r['label']:>25} {r['Te']:6.2f} {r['Tg']:5.0f} {r['ne']:10.2e} "
              f"{d.get('I',0):10.2e} {d.get('I2',0):10.2e} {d.get('I-',0):10.2e} "
              f"{r['diss_frac']:6.3f} {r['alpha']:7.4f} {r['nu_m']:10.2e}")

    # Delta table
    print("\n" + "=" * 90)
    print("  Relative Aenderung [%] gegenueber Nominal")
    print("=" * 90)
    print(f"{'Variante':>25} {'dTe%':>7} {'dne%':>7} {'ddiss%':>7} {'dalpha%':>9} {'dnu_m%':>8}")
    print("-" * 65)
    for r in all_results[1:]:
        dTe = (r["Te"] - nom["Te"]) / nom["Te"] * 100
        dne = (r["ne"] - nom["ne"]) / nom["ne"] * 100
        dd = (r["diss_frac"] - nom["diss_frac"]) / max(nom["diss_frac"], 1e-6) * 100
        da = (r["alpha"] - nom["alpha"]) / max(nom["alpha"], 1e-6) * 100
        dn = (r["nu_m"] - nom["nu_m"]) / max(nom["nu_m"], 1e-6) * 100
        print(f"{r['label']:>25} {dTe:+7.1f} {dne:+7.1f} {dd:+7.1f} {da:+9.1f} {dn:+8.1f}")

    # Recommendations
    print("\n" + "=" * 90)
    print("  EMPFEHLUNG")
    print("=" * 90)
    print("  Kelm_I2:")
    print("    Nominal: Ambalampitiya 2021 MTCS (neueste BSR, validiert durch Esteves 2022)")
    print("    Unsicherheitsband: [esteves_low, total_elastic] = [0.7x, ~1.5x nominal]")
    print("    Yadav_high (5x): nur fuer extreme Sensitivitaetstests")
    print()
    print("  Kdissatt_I2:")
    print("    Nominal: Ambalampitiya 2021 (BSR, normiert auf Truby)")
    print("    Unsicherheitsband: [hamilton_low, esteves_mid] = [0.3x, 1.5x nominal]")
    print("    Biondi_high (3x): oberes Limit (moegl. HI-kontaminiert)")


if __name__ == "__main__":
    main()
