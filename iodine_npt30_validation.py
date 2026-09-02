#!/usr/bin/env python3
"""
iodine_npt30_validation.py -- Experimentelle Validierung gegen NPT30-I2.

Referenz: Lafleur et al., Plasma Sources Sci. Technol. 31, 114001 (2022)
NPT30-I2: R~1.5cm, L~2.3cm, f~3MHz, 10 turns, mdot=56 ug/s I2

Fig.6: Ion beam composition (I+, I2+, I2+ fractions) vs. P_RF
Fig.7a: Ion beam current (mA) vs. P_RF for mdot~56 ug/s
"""
from __future__ import annotations
import sys, copy, math
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plasma_chemistry import load_chemistry, ThrusterGeometry, E_CH, KB, CONV, PI
from generic_solver import solve_steady_state
from rf_diagnostics import compute_rf_diagnostics

BASE = Path(__file__).resolve().parent / "chemistry" / "iodine_lafleur_v1"

# ═══════════════════════════════════════════════════════════════
# NPT30-I2 Geometrie (Lafleur 2022, Sec.3)
# ═══════════════════════════════════════════════════════════════
NPT30 = ThrusterGeometry(
    R=0.015,        # ~3cm diameter / 2
    L=0.023,        # ~2.3cm length
    betai=0.50,     # "just over 50%" screen transparency
    betag=0.25,     # effective Clausing ~25%
    Vgrid=1000.0,   # acceleration voltage
    sgrid=0.001,    # grid spacing ~1mm
    frequency=3.0e6,  # ~3 MHz
    Nw=10,          # 10-turn antenna
    R_ohm=0.5,      # estimated coil resistance
    Rc=0.018,       # coil slightly larger than chamber
    lc=0.025,       # coil length ~chamber length
)

# Mass flow rate: 56 ug/s I2 -> particles/s
# M_I2 = 253.8 u = 4.214e-25 kg
# 56e-9 kg/s / 4.214e-25 kg = 1.329e17 pps
MDOT_KGS = 56e-9
M_I2_KG = 4.2143422e-25
Q0_NPT30 = MDOT_KGS / M_I2_KG  # ~1.33e17 pps

# ═══════════════════════════════════════════════════════════════
# Experimentelle Referenzdaten (digitalisiert aus Fig.6 und Fig.7a)
# ═══════════════════════════════════════════════════════════════
# Fig.7a: Ion beam current vs. P_RF for iodine, mdot~56 ug/s
# Digitalisiert: P_RF [W], I_beam [mA], uncertainty ~+/-2mA
REF_BEAM_CURRENT = [
    {"P_RF_W": 5,  "I_beam_mA": 3.0,  "unc_mA": 2.0, "source": "Fig.7a"},
    {"P_RF_W": 7,  "I_beam_mA": 8.0,  "unc_mA": 2.0, "source": "Fig.7a"},
    {"P_RF_W": 10, "I_beam_mA": 15.0, "unc_mA": 2.0, "source": "Fig.7a"},
    {"P_RF_W": 13, "I_beam_mA": 20.0, "unc_mA": 3.0, "source": "Fig.7a"},
    {"P_RF_W": 15, "I_beam_mA": 22.0, "unc_mA": 3.0, "source": "Fig.7a"},
    {"P_RF_W": 20, "I_beam_mA": 25.0, "unc_mA": 3.0, "source": "Fig.7a"},
]

# Fig.6: Relative fraction of iodine species in beam, mdot=56 ug/s
# I+ fraction, I2+ fraction vs. P_RF
REF_BEAM_COMPOSITION = [
    {"P_RF_W": 5,  "f_Ip": 0.50, "f_I2p": 0.50, "unc": 0.10, "source": "Fig.6"},
    {"P_RF_W": 7,  "f_Ip": 0.65, "f_I2p": 0.35, "unc": 0.08, "source": "Fig.6"},
    {"P_RF_W": 10, "f_Ip": 0.80, "f_I2p": 0.20, "unc": 0.05, "source": "Fig.6"},
    {"P_RF_W": 13, "f_Ip": 0.87, "f_I2p": 0.13, "unc": 0.05, "source": "Fig.6"},
    {"P_RF_W": 15, "f_Ip": 0.90, "f_I2p": 0.10, "unc": 0.03, "source": "Fig.6"},
    {"P_RF_W": 20, "f_Ip": 0.90, "f_I2p": 0.10, "unc": 0.03, "source": "Fig.6"},
]


def make_lafleur2022_chem(chem):
    """Erstelle Lafleur 2022 Parametersatz."""
    c = copy.deepcopy(chem)
    for rxn in c.reactions:
        if rxn.id == "surfrec_I":
            rxn.surface_gamma = 0.07
        if rxn.id == "rec_Ip_Im":
            rxn.rate.value = 1.47e-14
        if rxn.id == "rec_I2p_Im":
            rxn.rate.value = 5.8e-14
    return c


def run_point(chem, P_abs, Q0=Q0_NPT30):
    """Einzelner Solverlauf, gibt erweiterte Diagnostik zurueck."""
    r = solve_steady_state(chem, NPT30, P_abs, Q0, max_iter=1200, tol=0.5)
    if not r["converged"]:
        return None

    d = r["densities"]
    nIp = d.get("I+", 0)
    nI2p = d.get("I2+", 0)
    nIm = d.get("I-", 0)
    ne = r["ne"]
    nI = d.get("I", 0)
    nI2 = d.get("I2", 0)

    # Beam current: I_beam = e * (Gamma_I+ + Gamma_I2+) * A_i
    # Gamma = h_L * n * u_B (vereinfacht)
    uB_Ip = math.sqrt(KB * r["Te"] * CONV / 2.107e-25)
    uB_I2p = math.sqrt(KB * r["Te"] * CONV / 4.214e-25)
    n_total_neutral = nI + nI2
    sigma_i = 1e-18
    lam = 1.0 / max(n_total_neutral * sigma_i, 1e-10)
    hL = 0.86 * (3 + NPT30.L / (2*lam))**(-0.5)
    Ai = NPT30.Ai
    I_Ip = E_CH * hL * nIp * uB_Ip * Ai * 1000  # mA
    I_I2p = E_CH * hL * nI2p * uB_I2p * Ai * 1000
    I_beam = I_Ip + I_I2p

    # Ion beam fractions
    I_total_ion = I_Ip + I_I2p
    f_Ip = I_Ip / I_total_ion if I_total_ion > 0 else 0
    f_I2p = I_I2p / I_total_ion if I_total_ion > 0 else 0

    rf = compute_rf_diagnostics(ne, r["nu_m"], NPT30)

    return {
        "Te": r["Te"], "ne": ne, "nI": nI, "nI2": nI2,
        "nIp": nIp, "nI2p": nI2p, "nIm": nIm,
        "I_beam_mA": I_beam, "I_Ip_mA": I_Ip, "I_I2p_mA": I_I2p,
        "f_Ip": f_Ip, "f_I2p": f_I2p,
        "diss": nI / (nI + 2*nI2) if (nI + 2*nI2) > 0 else 0,
        "alpha": nIm / ne if ne > 0 else 0,
        "zeta": rf.zeta, "nu_m": r["nu_m"],
    }


def scan_power(chem, label, P_RF_list):
    """Sweep ueber P_RF, berechne P_abs self-consistent via zeta."""
    results = []
    for P_RF in P_RF_list:
        # Iteriere: P_abs = zeta * P_RF
        P_abs = 0.6 * P_RF  # Startschaetzung
        for _ in range(5):
            r = run_point(chem, P_abs)
            if r is None:
                break
            zeta = r["zeta"]
            P_abs_new = max(zeta * P_RF, 0.5)
            if abs(P_abs_new - P_abs) / max(P_abs, 0.1) < 0.05:
                P_abs = P_abs_new
                break
            P_abs = P_abs_new

        if r:
            r["P_RF"] = P_RF
            r["P_abs"] = P_abs
            r["label"] = label
            results.append(r)
        else:
            results.append({"P_RF": P_RF, "label": label, "I_beam_mA": None})
    return results


def main():
    print("=" * 90)
    print("  NPT30-I2 Experimental Validation")
    print("  Reference: Lafleur et al., PSST 31, 114001 (2022)")
    print(f"  Geometry: R={NPT30.R*100:.1f}cm L={NPT30.L*100:.1f}cm f={NPT30.frequency/1e6:.0f}MHz")
    print(f"  Mass flow: {MDOT_KGS*1e6:.0f} ug/s -> Q0={Q0_NPT30:.2e} pps")
    print("=" * 90)

    chem_nom = load_chemistry(BASE / "chemistry.json")
    chem_laf = make_lafleur2022_chem(chem_nom)

    P_RF_list = [5, 7, 10, 13, 15, 20]

    # ── Nominaler Satz ────────────────────────────────────────
    res_nom = scan_power(chem_nom, "Nominal", P_RF_list)

    # ── Lafleur 2022 Satz ─────────────────────────────────────
    res_laf = scan_power(chem_laf, "Lafleur2022", P_RF_list)

    # ── Ergebnistabelle ───────────────────────────────────────
    print("\n--- Model Results ---")
    print(f"{'Set':>12} {'P_RF':>5} {'P_abs':>5} {'I_beam':>7} {'f_I+':>5} {'f_I2+':>6} "
          f"{'Te':>5} {'ne':>10} {'diss':>5} {'alpha':>6} {'zeta':>5}")
    print("-" * 85)
    for r in res_nom + res_laf:
        if r.get("I_beam_mA") is not None:
            print(f"{r['label']:>12} {r['P_RF']:5.0f} {r['P_abs']:5.1f} {r['I_beam_mA']:7.2f} "
                  f"{r['f_Ip']:5.2f} {r['f_I2p']:6.2f} {r['Te']:5.2f} {r['ne']:10.2e} "
                  f"{r['diss']:5.3f} {r['alpha']:6.4f} {r['zeta']:5.3f}")
        else:
            print(f"{r['label']:>12} {r['P_RF']:5.0f}   NC")

    # ── Vergleich: Beam Current (Fig.7a) ──────────────────────
    print("\n" + "=" * 90)
    print("  VERGLEICH: Ion Beam Current vs. Experiment (Fig.7a)")
    print("=" * 90)
    print(f"{'P_RF':>5} {'Exp[mA]':>8} {'Nom[mA]':>8} {'Laf[mA]':>8} {'Nom/Exp':>8} {'Laf/Exp':>8}")
    print("-" * 50)

    scores_nom = []
    scores_laf = []

    for ref in REF_BEAM_CURRENT:
        P = ref["P_RF_W"]
        I_exp = ref["I_beam_mA"]
        r_n = next((r for r in res_nom if r["P_RF"] == P and r.get("I_beam_mA")), None)
        r_l = next((r for r in res_laf if r["P_RF"] == P and r.get("I_beam_mA")), None)

        I_nom = r_n["I_beam_mA"] if r_n else None
        I_laf = r_l["I_beam_mA"] if r_l else None

        rat_n = I_nom / I_exp if I_nom and I_exp > 0 else None
        rat_l = I_laf / I_exp if I_laf and I_exp > 0 else None

        print(f"{P:5d} {I_exp:8.1f} {I_nom or 0:8.1f} {I_laf or 0:8.1f} "
              f"{rat_n:8.2f}" if rat_n else f"  --", end="")
        print(f" {rat_l:8.2f}" if rat_l else "   --")

        if rat_n:
            scores_nom.append(abs(math.log(rat_n)))
        if rat_l:
            scores_laf.append(abs(math.log(rat_l)))

    # ── Vergleich: Beam Composition (Fig.6) ───────────────────
    print("\n" + "=" * 90)
    print("  VERGLEICH: Ion Beam Composition vs. Experiment (Fig.6)")
    print("=" * 90)
    print(f"{'P_RF':>5} {'Exp_fI+':>8} {'Nom_fI+':>8} {'Laf_fI+':>8}")
    print("-" * 35)

    comp_scores_nom = []
    comp_scores_laf = []

    for ref in REF_BEAM_COMPOSITION:
        P = ref["P_RF_W"]
        fIp_exp = ref["f_Ip"]
        r_n = next((r for r in res_nom if r["P_RF"] == P and r.get("f_Ip") is not None), None)
        r_l = next((r for r in res_laf if r["P_RF"] == P and r.get("f_Ip") is not None), None)

        fIp_n = r_n["f_Ip"] if r_n else None
        fIp_l = r_l["f_Ip"] if r_l else None

        print(f"{P:5d} {fIp_exp:8.2f} {fIp_n or 0:8.2f} {fIp_l or 0:8.2f}")

        if fIp_n is not None:
            comp_scores_nom.append(abs(fIp_n - fIp_exp))
        if fIp_l is not None:
            comp_scores_laf.append(abs(fIp_l - fIp_exp))

    # ── Gesamtbewertung ───────────────────────────────────────
    print("\n" + "=" * 90)
    print("  GESAMTBEWERTUNG")
    print("=" * 90)

    avg_log_nom = np.mean(scores_nom) if scores_nom else float('inf')
    avg_log_laf = np.mean(scores_laf) if scores_laf else float('inf')
    avg_comp_nom = np.mean(comp_scores_nom) if comp_scores_nom else float('inf')
    avg_comp_laf = np.mean(comp_scores_laf) if comp_scores_laf else float('inf')

    print(f"  Beam Current (Fig.7a):")
    print(f"    Nominal:     avg |ln(model/exp)| = {avg_log_nom:.3f}")
    print(f"    Lafleur2022: avg |ln(model/exp)| = {avg_log_laf:.3f}")
    better_I = "Lafleur2022" if avg_log_laf < avg_log_nom else "Nominal"
    print(f"    -> {better_I} ist naeher am Experiment")

    print(f"\n  Beam Composition (Fig.6):")
    print(f"    Nominal:     avg |f_I+ error| = {avg_comp_nom:.3f}")
    print(f"    Lafleur2022: avg |f_I+ error| = {avg_comp_laf:.3f}")
    better_f = "Lafleur2022" if avg_comp_laf < avg_comp_nom else "Nominal"
    print(f"    -> {better_f} ist naeher am Experiment")

    print(f"""
  SCHLUSSFOLGERUNG:
    Die NPT30-I2 Geometrie (R=1.5cm) ist deutlich kleiner als die Grondein-
    Geometrie (R=6cm). Oberflaechenrekombination (gamma) hat hier staerkeren
    Einfluss, da A/V groesser ist.

    Beam Current: Beide Saetze sollten den experimentellen Trend qualitativ
    reproduzieren (steigende I_beam mit P_RF). Quantitative Abweichungen
    stammen primaer aus der vereinfachten P_RF -> P_abs Umrechnung.

    Beam Composition: Der I+/I2+ Anteil steigt mit Leistung (mehr Dissoziation).
    Lafleur 2022 (gamma=0.07) kann diesen Trend durch staerkere Wandrekombination
    beeinflussen.

    EMPFEHLUNG:
    - Fuer NPT30-I2 Vergleiche: Lafleur 2022 Satz (gamma=0.07, K_rec angepasst)
    - Fuer allgemeine Studien: Nominaler Satz (Ambalampitiya 2021, gamma=0.02)
    - Beide bleiben als explizite Varianten im Paket verfuegbar
""")


if __name__ == "__main__":
    main()
