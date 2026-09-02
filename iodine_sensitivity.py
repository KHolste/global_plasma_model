#!/usr/bin/env python3
"""
iodine_sensitivity.py -- Erster Iod-Solverlauf und Sensitivitaetsanalyse.

Nutzt chemistry/iodine_lafleur_v1 als Basis-Chemiepaket.
Fuehrt systematische Variationen der 5 kritischsten Prozesse durch.

Betriebsbedingungen orientiert an Grondein 2016 (Table I):
  R=6cm, L=10cm, f=13.56MHz, N=5, P_abs=200-800W, Q0=2.1e18 pps

Referenz: Lafleur et al. 2026, Grondein et al. 2016.
"""
from __future__ import annotations
import sys, copy, math
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plasma_chemistry import (
    load_chemistry, ChemistryPackage, BalanceAssembler, ThrusterGeometry,
    KB, CONV, E_CH, ME, PI
)
from generic_solver import solve_steady_state


# ═══════════════════════════════════════════════════════════════
# Betriebsbedingungen (Grondein 2016, Table I)
# ═══════════════════════════════════════════════════════════════
GEOM = ThrusterGeometry(
    R=0.06, L=0.10, betai=0.7, betag=0.3,
    Vgrid=1000.0, sgrid=0.0015,
    frequency=13.56e6, Nw=5, R_ohm=2.0, Rc=0.07, lc=0.10,
)

P_ABS_W = 400.0   # Absorbierte Leistung [W]
Q0_PPS = 2.1e18   # Teilchenfluss [1/s] (ca. 1 mg/s I2)

CHEM_PATH = Path(__file__).resolve().parent / "chemistry" / "iodine_lafleur_v1" / "chemistry.json"


def run_single(chem: ChemistryPackage, label: str = "nominal",
               P_abs: float = P_ABS_W, Q0: float = Q0_PPS) -> dict | None:
    """Einzelner Solverlauf. Gibt Ergebnisdict oder None bei Misserfolg zurueck."""
    result = solve_steady_state(
        chem, GEOM, P_abs_W=P_abs, Q0_pps=Q0,
        alpha_e_wall=7.0, density_profile_factor=1.0,
        max_iter=800, tol=0.5, verbose=False,
    )
    if not result["converged"]:
        return None

    # Erweiterte Diagnostik
    r = result
    ne = r["ne"]
    dens = r["densities"]
    n_I = dens.get("I", 0)
    n_I2 = dens.get("I2", 0)
    n_Ip = dens.get("I+", 0)
    n_I2p = dens.get("I2+", 0)
    n_Im = dens.get("I-", 0)

    r["label"] = label
    r["diss_frac"] = n_I / (n_I + 2*n_I2) if (n_I + 2*n_I2) > 0 else 0
    r["alpha"] = n_Im / ne if ne > 0 else 0
    r["n_ion_total"] = n_Ip + n_I2p
    r["I_beam_mA"] = E_CH * (n_Ip + n_I2p) * math.sqrt(KB * r["Te"] * CONV / 2.107e-25) * GEOM.Ai * 1000
    return r


def scale_rate(chem: ChemistryPackage, rxn_id: str, factor: float) -> ChemistryPackage:
    """Erzeuge Kopie mit skaliertem Ratenkoeffizienten fuer eine Reaktion."""
    c = copy.deepcopy(chem)
    for rxn in c.reactions:
        if rxn.id == rxn_id:
            if rxn.rate._table_Te is not None:
                rxn.rate._table_K = rxn.rate._table_K * factor
            elif rxn.rate.model.value == "constant":
                rxn.rate.value *= factor
            elif rxn.rate.model.value == "arrhenius":
                rxn.rate.A *= factor
            break
    return c


def main():
    print("=" * 80)
    print("  Iodine V1 Sensitivity Analysis")
    print("  Package: iodine_lafleur_v1 (Lafleur 2026)")
    print(f"  Geometry: R={GEOM.R*100:.0f}cm L={GEOM.L*100:.0f}cm f={GEOM.frequency/1e6:.2f}MHz")
    print(f"  P_abs={P_ABS_W:.0f}W Q0={Q0_PPS:.2e} pps")
    print("=" * 80)

    chem = load_chemistry(CHEM_PATH)

    # ── 1. Nominaler Lauf ─────────────────────────────────────
    print("\n--- Nominal ---")
    nom = run_single(chem, "nominal")
    if nom is None:
        print("FEHLER: Nominaler Lauf nicht konvergiert.")
        # Versuche niedrigere Leistung
        for P in [200, 300, 500, 600]:
            nom = run_single(chem, f"nominal_P{P}", P_abs=P)
            if nom:
                print(f"  Konvergiert bei P_abs={P}W")
                break
        if nom is None:
            print("  Kein konvergierter Nominalfall gefunden. Abbruch.")
            sys.exit(1)

    _print_result(nom)

    # ── 2. Sensitivitaetsanalyse ──────────────────────────────
    print("\n" + "=" * 80)
    print("  SENSITIVITAETSANALYSE")
    print("  Methode: +/- Faktor 2 auf einzelne Prozesse")
    print("=" * 80)

    P_use = float(nom.get("_P_abs", P_ABS_W)) if "_P_abs" in nom else P_ABS_W

    targets = [
        ("el_I",         "K_el,I (atomic elastic, no MTCS)"),
        ("el_I2",        "K_elm,I2 (molecular MTCS)"),
        ("diss_I2",      "K_diss,I2 (dissociation)"),
        ("dissatt_I2",   "K_dissatt,I2 (diss. attachment)"),
        ("exc_I_lumped", "K_exc,I (lumped excitation)"),
    ]

    factors = {"low": 0.5, "nom": 1.0, "high": 2.0}
    all_results = {"nominal": nom}

    for rxn_id, description in targets:
        print(f"\n--- {description} ---")
        for tag, fac in factors.items():
            if tag == "nom":
                continue
            label = f"{rxn_id}_{tag}"
            c_var = scale_rate(chem, rxn_id, fac)
            r = run_single(c_var, label, P_abs=P_use)
            if r:
                all_results[label] = r
                _print_result_compact(r, nom)
            else:
                print(f"  {label}: NICHT KONVERGIERT")
                all_results[label] = None

    # ── 3. Tornado-Tabelle ────────────────────────────────────
    print("\n" + "=" * 80)
    print("  TORNADO-ANALYSE: Relative Aenderung [%] gegenueber Nominal")
    print("=" * 80)

    metrics = [
        ("Te", "Te [eV]"),
        ("Tg", "Tg [K]"),
        ("ne", "ne [m-3]"),
        ("diss_frac", "Diss.frac"),
        ("alpha", "alpha (I-/ne)"),
        ("nu_m", "nu_m [s-1]"),
    ]

    header = f"{'Prozess':>25}"
    for _, mlabel in metrics:
        header += f" | {mlabel:>12}"
    print(header)
    print("-" * len(header))

    # Berechne max(|delta_low|, |delta_high|) fuer jede Kombination
    tornado_data = {}
    for rxn_id, desc in targets:
        low_key = f"{rxn_id}_low"
        high_key = f"{rxn_id}_high"
        r_low = all_results.get(low_key)
        r_high = all_results.get(high_key)

        row_low = f"{'x0.5 ' + rxn_id:>25}"
        row_high = f"{'x2.0 ' + rxn_id:>25}"
        max_impact = 0

        for mkey, mlabel in metrics:
            nom_val = nom.get(mkey, 0)
            if nom_val == 0:
                row_low += f" | {'--':>12}"
                row_high += f" | {'--':>12}"
                continue

            d_low = ((r_low[mkey] - nom_val) / nom_val * 100) if r_low else float('nan')
            d_high = ((r_high[mkey] - nom_val) / nom_val * 100) if r_high else float('nan')

            row_low += f" | {d_low:+12.1f}" if not math.isnan(d_low) else f" | {'FAIL':>12}"
            row_high += f" | {d_high:+12.1f}" if not math.isnan(d_high) else f" | {'FAIL':>12}"

            max_delta = max(abs(d_low) if not math.isnan(d_low) else 0,
                            abs(d_high) if not math.isnan(d_high) else 0)
            max_impact = max(max_impact, max_delta)

        tornado_data[rxn_id] = max_impact
        print(row_low)
        print(row_high)

    # ── 4. Prioritaetsranking ─────────────────────────────────
    print("\n" + "=" * 80)
    print("  PRIORITAETSRANKING: Max. Einfluss auf irgendeine Zielgroesse")
    print("=" * 80)

    sorted_impact = sorted(tornado_data.items(), key=lambda x: -x[1])
    for rank, (rxn_id, impact) in enumerate(sorted_impact, 1):
        desc = next(d for r, d in targets if r == rxn_id)
        print(f"  {rank}. {desc:50s} max |delta| = {impact:.1f}%")

    print("\n" + "=" * 80)
    print("  ZUSAMMENFASSUNG")
    print("=" * 80)
    if sorted_impact:
        top = sorted_impact[0]
        top_desc = next(d for r, d in targets if r == top[0])
        print(f"  Dominanteste Unsicherheit: {top_desc}")
        print(f"  Max. Einfluss: {top[1]:.1f}% Aenderung bei Faktor-2-Variation")
    print(f"  Status: experimental (review-basierte Basisdatenbank, nicht validiert)")


def _print_result(r: dict):
    print(f"  Te={r['Te']:.2f} eV  Tg={r['Tg']:.1f} K  ne={r['ne']:.2e}")
    for sp, n in r["densities"].items():
        print(f"    n_{sp} = {n:.3e}")
    print(f"    diss_frac = {r['diss_frac']:.4f}")
    print(f"    alpha = {r['alpha']:.4f}")
    print(f"    nu_m = {r['nu_m']:.3e}")
    print(f"    merit = {r['merit']:.4e}")


def _print_result_compact(r: dict, nom: dict):
    dTe = (r["Te"] - nom["Te"]) / nom["Te"] * 100 if nom["Te"] > 0 else 0
    dne = (r["ne"] - nom["ne"]) / nom["ne"] * 100 if nom["ne"] > 0 else 0
    print(f"  {r['label']:>25}: Te={r['Te']:.2f} ({dTe:+.1f}%)  "
          f"ne={r['ne']:.2e} ({dne:+.1f}%)  "
          f"diss={r['diss_frac']:.3f}  alpha={r['alpha']:.4f}  "
          f"nu_m={r['nu_m']:.2e}")


if __name__ == "__main__":
    main()
