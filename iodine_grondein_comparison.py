#!/usr/bin/env python3
"""
iodine_grondein_comparison.py -- Referenzvergleich: altes Iodmodell (I.2) vs aktuelles Projekt.

Vergleicht fuer Grondein-kompatible Parameter:
  1. Stationaere Ergebnisse bei festem P_abs
  2. Einfluss der Rekombinationsraten (Grondein vs Lafleur)
  3. Speziesverteilung, Te, Tg, alpha, diss, fIp

Referenz: Iodmodell/I.2 Ausgabe vom 22.07.2019
Geometrie: R=0.06, L=0.1, betai=0.7, betag=0.3, f=13.56MHz
"""
from __future__ import annotations
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from plasma_chemistry import load_chemistry, ThrusterGeometry, KB, CONV, E_CH
from generic_solver import solve_steady_state, _compute_beam_current

# ═══════════════════════════════════════════════════════════════
# Grondein-kompatible Geometrie
# ═══════════════════════════════════════════════════════════════
GRONDEIN_GEOM = ThrusterGeometry(
    R=0.06, L=0.10, betai=0.7, betag=0.3,
    Vgrid=1000.0, sgrid=0.0015,
    frequency=13.56e6, Nw=5, R_ohm=2.0, Rc=0.07, lc=0.10,
    eta_opt=1.0,  # Legacy: kein eta_opt
)

# ═══════════════════════════════════════════════════════════════
# Referenzdaten aus altem Iodmodell I.2 (22.07.2019)
# Betriebspunkt: I-fix 80 mA, verschiedene Q0 in mg/s
# Te ist in Kelvin gespeichert, hier in eV konvertiert
# ═══════════════════════════════════════════════════════════════
M_I2 = 4.2143422e-25  # kg

REF_DATA = [
    # Q0_mg_s, ne,       nIm,      nI2,       nI,        nIp,       nI2p,      Te_K,    Tg_K,  P_RF_W
    (0.30, 7.502e+16, 2.426e+15, 1.109e+18, 2.154e+18, 3.601e+16, 4.143e+16, 60932, 294.1, 72.0),
    (0.40, 8.124e+16, 4.594e+15, 1.717e+18, 3.885e+18, 4.226e+16, 4.357e+16, 50477, 294.7, 52.5),
    (0.50, 8.397e+16, 6.535e+15, 2.188e+18, 5.272e+18, 4.599e+16, 4.452e+16, 46180, 295.1, 47.4),  # ~0.48
    (0.60, 8.697e+16, 9.160e+15, 2.851e+18, 7.389e+18, 5.048e+16, 4.566e+16, 42116, 295.7, 44.2),
    (0.80, 9.010e+16, 1.346e+16, 3.940e+18, 1.089e+19, 5.595e+16, 4.760e+16, 38285, 296.7, 43.2),
    (1.00, 9.180e+16, 1.794e+16, 5.041e+18, 1.436e+19, 6.066e+16, 4.909e+16, 35860, 297.7, 43.6),
]


def mg_to_pps(mg_per_s):
    return mg_per_s * 1e-6 / M_I2


def run_comparison():
    chem = load_chemistry(SCRIPT_DIR / "chemistry" / "iodine_lafleur_v1" / "chemistry.json")

    print("=" * 100)
    print("  Referenzvergleich: Altes Iodmodell I.2 (Grondein 2016) vs Aktuelles Projekt (Lafleur V1)")
    print("  Geometrie: R=60mm, L=100mm, betai=0.7, betag=0.3, f=13.56MHz")
    print("  Modell I.2: Zeitintegration (BDF), Grondein-Raten, K_rec=9.3e-15, K_cex=1.22e-13")
    print("  Aktuell:    Stationaer (LM), Lafleur-Raten, K_rec=5e-14, K_rec_I2p=5e-14")
    print("=" * 100)

    print(f"\n{'Q0':>6} {'Modell':>8} {'Te eV':>7} {'Tg K':>7} {'ne':>10} {'nI-':>10} "
          f"{'nI2':>10} {'nI':>10} {'nI+':>10} {'nI2+':>10} "
          f"{'alpha':>7} {'diss':>6} {'fI+':>6}")
    print("-" * 130)

    for Q0_mg, ne_ref, nIm_ref, nI2_ref, nI_ref, nIp_ref, nI2p_ref, Te_K_ref, Tg_ref, P_ref in REF_DATA:
        Q0_pps = mg_to_pps(Q0_mg)
        Te_eV_ref = Te_K_ref * KB / E_CH

        # Referenzwerte
        alpha_ref = nIm_ref / ne_ref
        diss_ref = nI_ref / (nI_ref + 2 * nI2_ref)
        fIp_ref = nIp_ref / (nIp_ref + nI2p_ref)

        print(f"{Q0_mg:6.2f} {'I.2':>8} {Te_eV_ref:7.2f} {Tg_ref:7.1f} "
              f"{ne_ref:10.3e} {nIm_ref:10.3e} {nI2_ref:10.3e} {nI_ref:10.3e} "
              f"{nIp_ref:10.3e} {nI2p_ref:10.3e} "
              f"{alpha_ref:7.4f} {diss_ref:6.4f} {fIp_ref:6.4f}")

        # Aktuelles Modell bei der RF-Leistung des alten Modells
        # (P_ref ist P_RF, nicht P_abs; im aktuellen Modell ist P = P_abs direkt)
        # Fuer fairen Vergleich: P_abs ~ zeta * P_RF ~ 0.2 * P_RF (typisch)
        # Alternativ: feste Leistung vergleichen
        r = solve_steady_state(chem, GRONDEIN_GEOM, P_ref, Q0_pps, max_iter=800, tol=0.5)

        if r and r["converged"]:
            d = r["densities"]
            ne_cur = r["ne"]
            nI_cur = d.get("I", 0)
            nI2_cur = d.get("I2", 0)
            nIp_cur = d.get("I+", 0)
            nI2p_cur = d.get("I2+", 0)
            nIm_cur = d.get("I-", 0)
            Te_cur = r["Te"]
            Tg_cur = r["Tg"]

            alpha_cur = nIm_cur / ne_cur if ne_cur > 0 else 0
            diss_cur = nI_cur / (nI_cur + 2 * nI2_cur) if (nI_cur + 2 * nI2_cur) > 0 else 0
            fIp_cur = nIp_cur / (nIp_cur + nI2p_cur) if (nIp_cur + nI2p_cur) > 0 else 0

            print(f"{'':>6} {'Lafl.':>8} {Te_cur:7.2f} {Tg_cur:7.1f} "
                  f"{ne_cur:10.3e} {nIm_cur:10.3e} {nI2_cur:10.3e} {nI_cur:10.3e} "
                  f"{nIp_cur:10.3e} {nI2p_cur:10.3e} "
                  f"{alpha_cur:7.4f} {diss_cur:6.4f} {fIp_cur:6.4f}")
        else:
            print(f"{'':>6} {'Lafl.':>8} --- nicht konvergiert ---")

        print()

    # Rekombinationsraten-Audit
    print("\n" + "=" * 80)
    print("  Rekombinationsraten-Audit")
    print("=" * 80)
    print(f"\n{'Prozess':<30} {'Iodmodell I.2':>15} {'Lafleur V1':>15} {'Faktor':>8}")
    print("-" * 70)
    print(f"{'I+ + I- -> 2I':<30} {'9.311e-15':>15} {'5.0e-14':>15} {'5.4x':>8}")
    print(f"{'I2+ + I- -> I + I2':<30} {'1.22e-13':>15} {'5.0e-14':>15} {'2.4x':>8}")
    print()
    print("  Bewertung:")
    print("    Lafleur 2026 (K_rec=5e-14) basiert auf Greaves 2020 Messungen.")
    print("    Grondein 2016 (K_rec=9.3e-15) basiert auf aelteren Yeung-Werten.")
    print("    Der I.2-Wert fuer I2++I- (1.22e-13) ist als Ladungstausch modelliert,")
    print("    nicht als echte Rekombination. Lafleur behandelt beides als Rekombination.")
    print("    Hoehere Rekombinationsrate -> staerkerer I--Abbau -> niedrigeres alpha.")


if __name__ == "__main__":
    run_comparison()
