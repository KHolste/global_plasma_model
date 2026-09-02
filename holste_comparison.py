#!/usr/bin/env python3
"""
holste_comparison.py -- Xe/I2 Vergleich gegen Holste et al. 2018 (RIT-10).

Referenz: Holste et al., Eur. Phys. J. D 72, 9 (2018)
RIT-10: R=5cm, L=6.4cm, 6 turns, V_screen=1500V, f=1.1MHz
I_beam = I_screen - I_accel

Vergleicht modellierte P_RF(Q0) bei festem I_beam mit experimentellen
Performance-Mappings fuer Iod und Xenon.
"""
from __future__ import annotations
import sys, copy, math
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plasma_chemistry import load_chemistry, ThrusterGeometry
from generic_solver import solve_steady_state, solve_for_target_current, _compute_beam_current

CHEM_DIR = Path(__file__).resolve().parent / "chemistry"

# RIT-10 Geometry (Holste 2018, Sec.2.3)
RIT10 = ThrusterGeometry(
    R=0.05, L=0.064, betai=0.6, betag=0.20,
    Vgrid=1500.0, sgrid=0.004,  # 4mm grid spacing
    frequency=1.1e6, Nw=6, R_ohm=1.0, Rc=0.055, lc=0.07,
    eta_opt=0.60,  # Holste 2018: RIT-10, 243 Kanaele
)

# ═══════════════════════════════════════════════════════════════
# Holste 2018 Experimental Reference Data (digitized from Fig.5)
# P_RF [W] vs. mass_flow [mg/s] at 1.1 MHz, fixed I_beam
# ═══════════════════════════════════════════════════════════════
HOLSTE_IODINE = {
    80: [(0.28, 86), (0.35, 72), (0.40, 66), (0.45, 64), (0.50, 66), (0.55, 68), (0.60, 70)],
    85: [(0.35, 78), (0.40, 72), (0.45, 68), (0.50, 68), (0.55, 72), (0.60, 74)],
    90: [(0.40, 78), (0.45, 74), (0.50, 72), (0.55, 74)],
    95: [(0.40, 84), (0.45, 80), (0.50, 78), (0.55, 76)],
}

HOLSTE_XENON = {
    80: [(0.28, 72), (0.35, 66), (0.40, 64), (0.50, 62), (0.60, 60), (0.70, 60), (0.80, 62)],
    85: [(0.35, 74), (0.40, 68), (0.50, 66), (0.60, 66), (0.70, 66), (0.80, 67)],
    90: [(0.40, 77), (0.50, 72), (0.60, 70), (0.70, 70), (0.80, 72)],
    95: [(0.40, 84), (0.50, 78), (0.60, 74), (0.70, 74), (0.80, 75)],
}

SCCM_TO_PPS = 4.477962312e17
M_I2 = 4.2143422e-25  # kg
M_Xe = 2.1801711e-25  # kg


def mg_to_pps(mg_per_s, M_kg):
    """Konvertiere mg/s in particles/s."""
    return mg_per_s * 1e-6 / M_kg


def main():
    print("=" * 90)
    print("  Xe / I2 Performance Comparison: Model vs. Holste 2018 (RIT-10)")
    print("  Geometry: R=5cm, L=6.4cm, f=1.1MHz, V_screen=1500V")
    print("=" * 90)

    # Load chemistry packages
    chem_xe = load_chemistry(CHEM_DIR / "xenon_simple" / "chemistry.json")
    chem_i2 = load_chemistry(CHEM_DIR / "iodine_lafleur_v1" / "chemistry.json")

    # ── 1. Modell: SC bei verschiedenen Leistungen ──────────
    print("\n--- SC-Modus: I_beam vs. P_abs bei mdot=0.5 mg/s ---")
    print(f"{'P_abs':>6} {'I_Xe':>7} {'I_I2':>7}")
    print("-" * 25)

    for P in [10, 20, 50, 100, 200]:
        Q0_xe = mg_to_pps(0.5, M_Xe)
        Q0_i2 = mg_to_pps(0.5, M_I2)
        r_xe = solve_steady_state(chem_xe, RIT10, P, Q0_xe, max_iter=800, tol=0.5)
        r_i2 = solve_steady_state(chem_i2, RIT10, P, Q0_i2, max_iter=800, tol=0.5)
        I_xe = _compute_beam_current(r_xe, RIT10) if r_xe.get("converged") else 0
        I_i2 = _compute_beam_current(r_i2, RIT10) if r_i2.get("converged") else 0
        print(f"{P:6d} {I_xe:7.1f} {I_i2:7.1f}")

    # ── 2. Experimentelle Referenzen ──────────────────────────
    print(f"\n--- Holste 2018 Experimental Reference (1.1 MHz, RIT-10) ---")
    print(f"{'I_beam':>6} {'mdot':>6} {'P_I2':>5} {'P_Xe':>5} {'Ratio':>6}")
    print("-" * 35)

    for I_target in [80, 85, 90, 95]:
        for mdot_i, P_i in HOLSTE_IODINE.get(I_target, []):
            P_xe = None
            for mdot_x, P_x in HOLSTE_XENON.get(I_target, []):
                if abs(mdot_x - mdot_i) < 0.03:
                    P_xe = P_x
                    break
            if P_xe:
                ratio = P_i / P_xe
                print(f"{I_target:6d} {mdot_i:6.2f} {P_i:5.0f} {P_xe:5.0f} {ratio:6.2f}")

    # ── 3. Modell-Limitierung dokumentieren ────────────────────
    print(f"\n--- Modell I_beam Limit bei RIT-10 ---")
    max_I = {}
    for gas, chem, M in [("Xe", chem_xe, M_Xe), ("I2", chem_i2, M_I2)]:
        Q0 = mg_to_pps(0.5, M)
        r = solve_steady_state(chem, RIT10, 200.0, Q0, max_iter=800, tol=0.5)
        I = _compute_beam_current(r, RIT10) if r.get("converged") else 0
        max_I[gas] = I
        print(f"  {gas}: I_beam_max(P=200W, mdot=0.5) = {I:.1f} mA")

    print(f"\n  Experiment: 80-95 mA bei 60-86 W")
    print(f"  Modell:     {max(max_I.values()):.0f} mA bei 200 W (CL + Grid-Optik limitiert)")
    print(f"  -> Faktor {80/max(max(max_I.values()),1):.1f}x zu niedrig")

    # Summary
    print(f"\n{'='*90}")
    print("  HINWEISE")
    print(f"{'='*90}")
    print("  - Holste 2018: P_RF ist RFG-DC-Leistung (inkl. Ohm-Verluste)")
    print("  - Modell: P ist absorbierte Leistung (P_abs < P_RF)")
    print("  - Erwarteter Unterschied: P_model < P_exp (Faktor ~0.5-0.8)")
    print("  - Experiment bei 1.1 MHz, Modell-RF ebenfalls 1.1 MHz")
    print("  - I_beam Def.: Holste: I_screen - I_accel; Modell: Bohm-Extraktion")


if __name__ == "__main__":
    main()
