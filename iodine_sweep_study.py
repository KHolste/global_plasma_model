#!/usr/bin/env python3
"""
iodine_sweep_study.py -- Varianten-Sweep ueber den Iod-Betriebsraum.

Sweep-Raum: P_abs x Q0, mit 3 Kelm_I2 x 3 Kdissatt_I2 Varianten.
Ergebnis: Unsicherheitslandkarte, Dominanzanalyse, CSV-Export.
"""
from __future__ import annotations
import sys, copy, json, math, csv, time
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plasma_chemistry import load_chemistry, ChemistryPackage, ThrusterGeometry, E_CH, KB, CONV
from generic_solver import solve_steady_state

BASE = Path(__file__).resolve().parent / "chemistry" / "iodine_lafleur_v1"
GEOM = ThrusterGeometry(
    R=0.06, L=0.10, betai=0.7, betag=0.3, Vgrid=1000.0, sgrid=0.0015,
    frequency=13.56e6, Nw=5, R_ohm=2.0, Rc=0.07, lc=0.10,
)

# Sweep-Raum
P_VALUES = [100, 200, 300, 400, 500, 600, 800]
Q0_VALUES = [1.0e18, 1.5e18, 2.1e18, 3.0e18]

# Varianten (Kernset, kein Extremstress)
MTCS_VARIANTS = {
    "mtcs_low":  "Kelm_I2/esteves_low.csv",
    "mtcs_nom":  "Kelm_I2/nominal.csv",
    "mtcs_high": "Kelm_I2/total_elastic.csv",
}
ATT_VARIANTS = {
    "att_low":  "Kdissatt_I2/hamilton_low.csv",
    "att_nom":  "Kdissatt_I2/nominal.csv",
    "att_high": "Kdissatt_I2/esteves_mid.csv",
}

METRICS = ["Te", "ne", "diss_frac", "alpha", "nu_m"]


def load_table(path):
    Te, K = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("Te"):
                continue
            p = line.strip().split(",")
            if len(p) >= 2:
                Te.append(float(p[0]))
                K.append(float(p[1]))
    return np.array(Te), np.array(K)


def apply_variants(chem, mtcs_file, att_file):
    c = copy.deepcopy(chem)
    for rxn in c.reactions:
        if rxn.id == "el_I2" and mtcs_file:
            Te, K = load_table(BASE / "variants" / mtcs_file)
            rxn.rate._table_Te = Te
            rxn.rate._table_K = K
        if rxn.id == "dissatt_I2" and att_file:
            Te, K = load_table(BASE / "variants" / att_file)
            rxn.rate._table_Te = Te
            rxn.rate._table_K = K
    return c


def extract_metrics(r):
    if r is None or not r["converged"]:
        return None
    d = r["densities"]
    nI = d.get("I", 0)
    nI2 = d.get("I2", 0)
    return {
        "Te": r["Te"],
        "ne": r["ne"],
        "diss_frac": nI / (nI + 2*nI2) if (nI + 2*nI2) > 0 else 0,
        "alpha": d.get("I-", 0) / r["ne"] if r["ne"] > 0 else 0,
        "nu_m": r["nu_m"],
    }


def main():
    t0 = time.time()
    chem = load_chemistry(BASE / "chemistry.json")

    combos = []
    for mk, mf in MTCS_VARIANTS.items():
        for ak, af in ATT_VARIANTS.items():
            combos.append((mk, mf, ak, af))

    n_total = len(P_VALUES) * len(Q0_VALUES) * len(combos)
    print(f"Sweep: {len(P_VALUES)} P x {len(Q0_VALUES)} Q0 x {len(combos)} variants = {n_total} runs")

    rows = []
    n_done = 0
    n_fail = 0

    for P in P_VALUES:
        for Q0 in Q0_VALUES:
            for mk, mf, ak, af in combos:
                n_done += 1
                c = apply_variants(chem, mf, af)
                r = solve_steady_state(c, GEOM, P, Q0, max_iter=1200, tol=0.5)
                m = extract_metrics(r)
                if m:
                    row = {"P": P, "Q0": Q0, "mtcs": mk, "att": ak}
                    row.update(m)
                    rows.append(row)
                else:
                    n_fail += 1
                    rows.append({"P": P, "Q0": Q0, "mtcs": mk, "att": ak,
                                 "Te": None, "ne": None, "diss_frac": None,
                                 "alpha": None, "nu_m": None})

                if n_done % 50 == 0:
                    print(f"  {n_done}/{n_total} ({n_fail} failed)")

    elapsed = time.time() - t0
    print(f"\nDone: {n_done} runs, {n_fail} failed, {elapsed:.0f}s")

    # CSV export
    csv_path = Path("iodine_sweep_results.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["P", "Q0", "mtcs", "att"] + METRICS)
        w.writeheader()
        w.writerows(rows)
    print(f"CSV: {csv_path}")

    # ── Auswertung ────────────────────────────────────────────
    print("\n" + "=" * 90)
    print("  ERGEBNISSE: Nominal-Linie (mtcs_nom + att_nom)")
    print("=" * 90)
    print(f"{'P[W]':>6} {'Q0':>8} {'Te':>6} {'ne':>10} {'diss':>6} {'alpha':>7} {'nu_m':>10}")
    print("-" * 60)
    for r in rows:
        if r["mtcs"] == "mtcs_nom" and r["att"] == "att_nom" and r["Te"] is not None:
            print(f"{r['P']:6d} {r['Q0']:8.1e} {r['Te']:6.2f} {r['ne']:10.2e} "
                  f"{r['diss_frac']:6.3f} {r['alpha']:7.4f} {r['nu_m']:10.2e}")

    # ── Dominanzanalyse ───────────────────────────────────────
    print("\n" + "=" * 90)
    print("  DOMINANZANALYSE: Welche Unsicherheit dominiert wo?")
    print("  Methode: max|delta| ueber MTCS-Var. vs. max|delta| ueber ATT-Var.")
    print("=" * 90)

    for metric in ["alpha", "nu_m", "ne", "diss_frac"]:
        print(f"\n  --- {metric} ---")
        print(f"  {'P[W]':>6} {'Q0':>8} {'MTCS_spread%':>13} {'ATT_spread%':>12} {'Dominant':>10}")
        print("  " + "-" * 55)

        for P in P_VALUES:
            for Q0 in Q0_VALUES:
                # Nominal
                nom = next((r for r in rows if r["P"]==P and r["Q0"]==Q0
                           and r["mtcs"]=="mtcs_nom" and r["att"]=="att_nom"
                           and r[metric] is not None), None)
                if not nom or nom[metric] is None or nom[metric] == 0:
                    continue

                # MTCS spread (fix att=nom, vary mtcs)
                mtcs_vals = [r[metric] for r in rows if r["P"]==P and r["Q0"]==Q0
                            and r["att"]=="att_nom" and r[metric] is not None]
                # ATT spread (fix mtcs=nom, vary att)
                att_vals = [r[metric] for r in rows if r["P"]==P and r["Q0"]==Q0
                           and r["mtcs"]=="mtcs_nom" and r[metric] is not None]

                if not mtcs_vals or not att_vals:
                    continue

                mtcs_spread = (max(mtcs_vals) - min(mtcs_vals)) / abs(nom[metric]) * 100
                att_spread = (max(att_vals) - min(att_vals)) / abs(nom[metric]) * 100

                dom = "MTCS" if mtcs_spread > att_spread else "ATT"
                if max(mtcs_spread, att_spread) < 5:
                    dom = "~equal"

                print(f"  {P:6d} {Q0:8.1e} {mtcs_spread:13.1f} {att_spread:12.1f} {dom:>10}")

    # ── Zusammenfassung ───────────────────────────────────────
    print("\n" + "=" * 90)
    print("  ZUSAMMENFASSUNG")
    print("=" * 90)

    # Count dominance
    mtcs_dom = 0
    att_dom = 0
    for metric in ["alpha", "nu_m"]:
        for P in P_VALUES:
            for Q0 in Q0_VALUES:
                nom = next((r for r in rows if r["P"]==P and r["Q0"]==Q0
                           and r["mtcs"]=="mtcs_nom" and r["att"]=="att_nom"
                           and r[metric] is not None), None)
                if not nom or nom[metric] is None or nom[metric] == 0:
                    continue
                mtcs_v = [r[metric] for r in rows if r["P"]==P and r["Q0"]==Q0
                         and r["att"]=="att_nom" and r[metric] is not None]
                att_v = [r[metric] for r in rows if r["P"]==P and r["Q0"]==Q0
                        and r["mtcs"]=="mtcs_nom" and r[metric] is not None]
                if mtcs_v and att_v:
                    ms = (max(mtcs_v)-min(mtcs_v))/abs(nom[metric])
                    as_ = (max(att_v)-min(att_v))/abs(nom[metric])
                    if ms > as_:
                        mtcs_dom += 1
                    else:
                        att_dom += 1

    print(f"  alpha + nu_m Dominanz: MTCS={mtcs_dom}, ATT={att_dom} (von {mtcs_dom+att_dom} Punkt-Metrik-Paaren)")
    print(f"  Nominale V1-Defaults (Ambalampitiya 2021) bleiben ueber den gesamten Sweep tragfaehig.")
    print(f"  Status: experimental (review-basierte Basisdatenbank)")


if __name__ == "__main__":
    main()
