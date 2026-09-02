#!/usr/bin/env python3
"""
holste_ifix_diagnosis.py -- Diagnose des I-fix-Modus fuer Holste/RIT-10 mit Iod.

Prueft:
1. I_beam(P) Sweep: Wie steigt der Beamstrom mit der Leistung?
2. SC vs I-fix Konsistenz: Liefern beide Pfade dasselbe I_beam bei gleichem P?
3. Xe vs I2 Vergleich: Qualitativ in Richtung Holste 2018?
4. Extraktionslimits: Was begrenzt den Beamstrom (CL, Bohm, eta_opt)?
5. Parametrische Sensitivitaet: Welche Parameter limitieren am staerksten?

Referenz: Holste et al., Eur. Phys. J. D 72, 9 (2018)
  - RIT-10: I_beam = 80-95 mA bei P_RF = 60-86 W (DC-RFG-Leistung)
"""
from __future__ import annotations
import sys, math
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plasma_chemistry import load_chemistry, ThrusterGeometry, KB, CONV, E_CH, PI, EPS0
from generic_solver import solve_steady_state, _compute_beam_current
from beam_extraction import compute_extraction

CHEM_DIR = Path(__file__).resolve().parent / "chemistry"
SCCM_TO_PPS = 4.477962312e17
M_I2 = 4.2143422e-25  # kg (I2)
M_Xe = 2.1801711e-25  # kg (Xe)
M_I = 2.107e-25        # kg (I)


def holste_geom(eta_opt=0.60):
    """RIT-10 Geometrie nach Holste 2018, Sec. 2.3."""
    return ThrusterGeometry(
        R=0.05, L=0.064, betai=0.6, betag=0.20,
        Vgrid=1500.0, sgrid=0.004,
        frequency=1.1e6, Nw=6, R_ohm=1.0, Rc=0.055, lc=0.07,
        eta_opt=eta_opt,
    )


def mg_to_pps(mg_per_s, M_kg):
    return mg_per_s * 1e-6 / M_kg


def run_sc_point(chem, geom, P_W, Q0_pps):
    """Loese SC-Modus und berechne I_beam plus Extraktionsdetails."""
    r = solve_steady_state(chem, geom, P_W, Q0_pps, max_iter=800, tol=0.5)
    if not r or not r.get("converged"):
        return None

    d = r["densities"]
    I_beam = _compute_beam_current(r, geom)

    # Extraktionsdetails separat berechnen
    ion_dens = {}
    ion_mass = {}
    for sp_id in ("I+", "I2+", "Xe+"):
        if sp_id in d and d[sp_id] > 0:
            ion_dens[sp_id] = d[sp_id]
            if sp_id == "I+":
                ion_mass[sp_id] = M_I
            elif sp_id == "I2+":
                ion_mass[sp_id] = M_I2
            elif sp_id == "Xe+":
                ion_mass[sp_id] = M_Xe

    n_neut = sum(v for k, v in d.items() if k not in ("I+", "I2+", "I-", "Xe+"))

    ex = compute_extraction(
        ion_densities=ion_dens,
        ion_masses=ion_mass,
        Te_eV=r["Te"],
        n_neutral_total=max(n_neut, 1e10),
        sigma_i=1e-18,
        R=geom.R, L=geom.L,
        betai=geom.betai,
        Vgrid=geom.Vgrid,
        sgrid=geom.sgrid,
        eta_opt=geom.eta_opt,
    )

    # Iod-spezifische Groessen
    nI = d.get("I", 0)
    nI2 = d.get("I2", 0)
    nIp = d.get("I+", d.get("Xe+", 0))
    nI2p = d.get("I2+", 0)
    nIm = d.get("I-", 0)
    ne = r["ne"]
    n_ion = nIp + nI2p
    diss = nI / (nI + 2 * nI2) if (nI + 2 * nI2) > 0 else 0
    fIp = nIp / n_ion if n_ion > 0 else 0
    alpha = nIm / ne if ne > 0 else 0

    return {
        "P": P_W, "I_beam": I_beam, "Te": r["Te"], "ne": ne,
        "diss": diss, "fIp": fIp, "alpha": alpha,
        "J_Bohm": ex.J_Bohm_total, "J_CL": ex.J_CL,
        "I_Bohm_limit": ex.I_Bohm_limit_mA, "I_CL_limit": ex.I_CL_limit_mA,
        "limiting": ex.limiting, "eta_opt": ex.eta_opt,
        "nIp": nIp, "nI2p": nI2p, "nIm": nIm,
    }


def section(title):
    print(f"\n{'='*90}")
    print(f"  {title}")
    print(f"{'='*90}")


def main():
    section("HOLSTE RIT-10 I-FIX DIAGNOSE")
    print("  Referenz: Holste et al. 2018, I_beam = 80-95 mA bei P_RF ~ 60-86 W")
    print("  Geom: R=5cm, L=6.4cm, betai=0.6, Vgrid=1500V, sgrid=4mm, f=1.1MHz")

    chem_i2 = load_chemistry(CHEM_DIR / "iodine_lafleur_v1" / "chemistry.json")
    chem_xe = load_chemistry(CHEM_DIR / "xenon_simple" / "chemistry.json")
    geom = holste_geom(eta_opt=0.60)

    mdot = 0.5  # mg/s (typisch Holste)
    Q0_i2 = mg_to_pps(mdot, M_I2)
    Q0_xe = mg_to_pps(mdot, M_Xe)

    # ── 1. I_beam(P) Sweep fuer I2 mit korrektem eta_opt ──────
    section("1. I_beam(P) Sweep: Iod, Holste-Geometrie, eta_opt=0.60")
    print(f"   Massenfluss: {mdot} mg/s, Q0_I2 = {Q0_i2:.3e} pps")
    print()
    print(f"{'P[W]':>7} {'I_beam':>8} {'Te':>6} {'ne':>10} {'diss':>6} "
          f"{'f_I+':>6} {'alpha':>7} {'J_Bohm':>8} {'J_CL':>8} {'Limit':>12} "
          f"{'I_Bohm_lim':>10} {'I_CL_lim':>10}")
    print("-" * 120)

    P_values = [5, 10, 20, 30, 50, 75, 100, 150, 200, 300, 500, 750, 1000]
    results_i2 = []
    for P in P_values:
        pt = run_sc_point(chem_i2, geom, P, Q0_i2)
        if pt:
            results_i2.append(pt)
            print(f"{pt['P']:7.0f} {pt['I_beam']:8.2f} {pt['Te']:6.2f} {pt['ne']:10.3e} "
                  f"{pt['diss']:6.3f} {pt['fIp']:6.3f} {pt['alpha']:7.4f} "
                  f"{pt['J_Bohm']:8.2f} {pt['J_CL']:8.2f} {pt['limiting']:>12} "
                  f"{pt['I_Bohm_limit']:10.2f} {pt['I_CL_limit']:10.2f}")
        else:
            print(f"{P:7.0f}  -- nicht konvergiert --")

    # Pruefen ob 80 mA erreichbar
    max_I = max((pt["I_beam"] for pt in results_i2), default=0)
    P_at_80 = None
    for pt in results_i2:
        if pt["I_beam"] >= 80:
            P_at_80 = pt["P"]
            break
    print()
    if P_at_80:
        print(f"  >>> 80 mA erreichbar ab P ~ {P_at_80} W")
    else:
        print(f"  >>> 80 mA NICHT erreichbar. Maximum: {max_I:.1f} mA")
        if results_i2:
            lim = results_i2[-1]["limiting"]
            print(f"  >>> Limitierender Mechanismus bei P_max: {lim}")

    # ── 2. Vergleich: eta_opt=0.25 (Bug) vs 0.60 (korrekt) ──────
    section("2. Vergleich: eta_opt=0.25 (alter Hardcode) vs. 0.60 (Preset)")
    geom_bug = holste_geom(eta_opt=0.25)
    print(f"{'P[W]':>7} {'I(0.25)':>9} {'I(0.60)':>9} {'Faktor':>8}")
    print("-" * 40)
    for P in [50, 100, 200, 500]:
        pt_bug = run_sc_point(chem_i2, geom_bug, P, Q0_i2)
        pt_fix = run_sc_point(chem_i2, geom, P, Q0_i2)
        I_bug = pt_bug["I_beam"] if pt_bug else 0
        I_fix = pt_fix["I_beam"] if pt_fix else 0
        fac = I_fix / I_bug if I_bug > 0 else float('inf')
        print(f"{P:7.0f} {I_bug:9.2f} {I_fix:9.2f} {fac:8.2f}x")

    # ── 3. Vergleich: Default-Geometrie (Bug) vs Holste-Geometrie ──
    section("3. Vergleich: Default-Geometrie (Bug) vs. Holste-Geometrie")
    geom_default = ThrusterGeometry()  # R=0.02, L=0.04, etc.
    print(f"  Default: R={geom_default.R}, L={geom_default.L}, betai={geom_default.betai}, "
          f"sgrid={geom_default.sgrid}, eta_opt={geom_default.eta_opt}")
    print(f"  Holste:  R={geom.R}, L={geom.L}, betai={geom.betai}, "
          f"sgrid={geom.sgrid}, eta_opt={geom.eta_opt}")
    print()
    print(f"{'P[W]':>7} {'I(default)':>11} {'I(Holste)':>11} {'Faktor':>8}")
    print("-" * 45)
    for P in [50, 100, 200, 500]:
        pt_def = run_sc_point(chem_i2, geom_default, P, Q0_i2)
        pt_hol = run_sc_point(chem_i2, geom, P, Q0_i2)
        I_def = pt_def["I_beam"] if pt_def else 0
        I_hol = pt_hol["I_beam"] if pt_hol else 0
        fac = I_hol / I_def if I_def > 0 else float('inf')
        print(f"{P:7.0f} {I_def:11.2f} {I_hol:11.2f} {fac:8.1f}x")

    # ── 4. Xe vs I2 Vergleich ──────────────────────────────────
    section("4. Xenon vs. Iod bei Holste-Geometrie (eta_opt=0.60)")
    print(f"   mdot = {mdot} mg/s")
    print()
    print(f"{'P[W]':>7} {'I_Xe':>8} {'I_I2':>8} {'Ratio':>7} {'Te_Xe':>6} {'Te_I2':>6}")
    print("-" * 50)
    for P in [20, 50, 100, 200, 500]:
        pt_xe = run_sc_point(chem_xe, geom, P, Q0_xe)
        pt_i2 = run_sc_point(chem_i2, geom, P, Q0_i2)
        I_xe = pt_xe["I_beam"] if pt_xe else 0
        I_i2 = pt_i2["I_beam"] if pt_i2 else 0
        Te_xe = pt_xe["Te"] if pt_xe else 0
        Te_i2 = pt_i2["Te"] if pt_i2 else 0
        ratio = I_i2 / I_xe if I_xe > 0 else float('inf')
        print(f"{P:7.0f} {I_xe:8.2f} {I_i2:8.2f} {ratio:7.2f} {Te_xe:6.2f} {Te_i2:6.2f}")

    # ── 5. CL-Limit-Analyse fuer sgrid ──────────────────────────
    section("5. Child-Langmuir Limit: Sensitivitaet auf sgrid")
    print(f"  Vgrid = 1500 V, M_eff = I+ (127 u)")
    print()
    A_grid = geom.betai * PI * geom.R ** 2
    for sg in [0.0005, 0.001, 0.002, 0.004, 0.006]:
        J_CL = (4.0 / 9.0) * EPS0 * math.sqrt(2 * E_CH / M_I) * 1500.0**1.5 / sg**2
        I_CL = J_CL * A_grid * geom.eta_opt * 1000  # mA
        print(f"  sgrid={sg*1000:.1f}mm: J_CL={J_CL:.1f} A/m², "
              f"I_CL_limit={I_CL:.1f} mA (eta_opt={geom.eta_opt})")

    # ── 6. Parametrische Sensitivitaet ──────────────────────────
    section("6. Parametrische Sensitivitaet: Welcher Parameter limitiert am staerksten?")
    P_test = 100.0  # W
    print(f"  Basis: P={P_test}W, mdot={mdot}mg/s, eta_opt={geom.eta_opt}")
    print()

    base = run_sc_point(chem_i2, geom, P_test, Q0_i2)
    I_base = base["I_beam"] if base else 0
    print(f"  Basis I_beam = {I_base:.2f} mA")
    print()

    variations = [
        ("eta_opt=0.80", {"eta_opt": 0.80}),
        ("eta_opt=1.00", {"eta_opt": 1.00}),
        ("sgrid=1mm",    {"sgrid": 0.001}),
        ("sgrid=0.5mm",  {"sgrid": 0.0005}),
        ("betai=0.80",   {"betai": 0.80}),
        ("Vgrid=2000V",  {"Vgrid": 2000.0}),
    ]
    print(f"  {'Variation':<20} {'I_beam':>8} {'Delta':>8} {'Faktor':>8} {'Limit':>12}")
    print("  " + "-" * 65)
    for name, overrides in variations:
        g = holste_geom(eta_opt=overrides.pop("eta_opt", geom.eta_opt))
        for k, v in overrides.items():
            setattr(g, k, v)
        pt = run_sc_point(chem_i2, g, P_test, Q0_i2)
        I_var = pt["I_beam"] if pt else 0
        delta = I_var - I_base
        fac = I_var / I_base if I_base > 0 else float('inf')
        lim = pt["limiting"] if pt else "?"
        print(f"  {name:<20} {I_var:8.2f} {delta:+8.2f} {fac:8.2f}x {lim:>12}")

    # ── 7. Holste-Referenzdaten vs. Modell ──────────────────────
    section("7. Holste 2018 Referenz vs. Modell")
    HOLSTE_IODINE_80 = [
        (0.28, 86), (0.35, 72), (0.40, 66), (0.45, 64),
        (0.50, 66), (0.55, 68), (0.60, 70),
    ]
    print(f"  I_target = 80 mA (Iod)")
    print(f"  {'mdot[mg/s]':>10} {'P_exp[W]':>9} {'I_model':>9} {'Limit':>12}")
    print("  " + "-" * 50)
    for mdot_pt, P_exp in HOLSTE_IODINE_80:
        Q0 = mg_to_pps(mdot_pt, M_I2)
        pt = run_sc_point(chem_i2, geom, P_exp, Q0)
        I_mod = pt["I_beam"] if pt else 0
        lim = pt["limiting"] if pt else "?"
        print(f"  {mdot_pt:10.2f} {P_exp:9.0f} {I_mod:9.2f} {lim:>12}")

    # ── 8. Zusammenfassung ──────────────────────────────────────
    section("ZUSAMMENFASSUNG")
    print("""
  BUG 1 (KRITISCH, BEHOBEN): _compute_beam_current() hatte eta_opt=0.25 hardcoded.
    Holste-Preset definiert eta_opt=0.60 -> Faktor 2.4x zu niedrig.

  BUG 2 (KRITISCH, BEHOBEN): generic_solver.py main() ignorierte Preset-Geometrie.
    Nutzte ThrusterGeometry()-Defaults (R=2cm) statt Holste (R=5cm).
    Extraktionsflaeche war (2/5)^2 = 16% des Sollwerts.

  BUG 3 (BEHOBEN): P_max war auf 500 W hardcoded statt aus params.txt zu lesen.

  KUMULIERTER EFFEKT: I_beam war ~7% des korrekten Werts.
    Kein Wunder, dass 80 mA bei 500 W unerreichbar waren.

  OFFENE FRAGE: sgrid=4mm im Holste-Preset.
    Die Preset-Notiz sagt "4mm Grid-Loecher" (Apertur), aber sgrid ist der
    Gitter*abstand*. Fuer RIT-10 waere ein Gap von ~0.5-1mm typischer.
    Bei sgrid=4mm ist das CL-Limit stark restriktiv.
    -> Pruefen ob sgrid die Apertur oder den Gap meint.
""")


if __name__ == "__main__":
    main()
