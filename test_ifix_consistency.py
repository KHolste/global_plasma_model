#!/usr/bin/env python3
"""
test_ifix_consistency.py -- Tests fuer I-fix/SC-Konsistenz und eta_opt-Propagation.

Prueft:
  1. eta_opt fliesst von ThrusterGeometry durch _compute_beam_current zu beam_extraction
  2. SC- und I-fix-Modus liefern konsistente I_beam bei gleichem P
  3. Preset-Geometrie wird korrekt aus params.txt gelesen
  4. Holste-Iod-Fall erreicht physikalisch plausible Beamstroeme
  5. Xenon/Iod-Vergleich regressiert nicht

Ausfuehrung:
    python test_ifix_consistency.py
"""
from __future__ import annotations

import sys, math, tempfile, os
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
os.chdir(SCRIPT_DIR)

from plasma_chemistry import load_chemistry, ThrusterGeometry
from generic_solver import (
    solve_steady_state, solve_for_target_current,
    _compute_beam_current, _load_geometry_from_config,
)
from beam_extraction import compute_extraction

CHEM_DIR = SCRIPT_DIR / "chemistry"
SCCM_TO_PPS = 4.477962312e17
M_I2 = 4.2143422e-25
M_Xe = 2.1801711e-25
M_I = 2.107e-25

passed = 0
failed = 0
errors = []


def check(name, condition, detail=""):
    global passed, failed
    if condition:
        passed += 1
        print(f"  PASS: {name}")
    else:
        failed += 1
        msg = f"  FAIL: {name}"
        if detail:
            msg += f" -- {detail}"
        print(msg)
        errors.append(name)


def holste_geom(eta_opt=0.60):
    return ThrusterGeometry(
        R=0.05, L=0.064, betai=0.6, betag=0.20,
        Vgrid=1500.0, sgrid=0.004,
        frequency=1.1e6, Nw=6, R_ohm=1.0, Rc=0.055, lc=0.07,
        eta_opt=eta_opt,
    )


# ═══════════════════════════════════════════════════════════════
# Test 1: eta_opt Propagation
# ═══════════════════════════════════════════════════════════════
def test_eta_opt_propagation():
    print("\n--- Test 1: eta_opt Propagation ---")

    # ThrusterGeometry hat eta_opt Feld
    g1 = ThrusterGeometry()
    check("ThrusterGeometry has eta_opt", hasattr(g1, "eta_opt"))
    check("Default eta_opt is 1.0", g1.eta_opt == 1.0, f"got {g1.eta_opt}")

    g2 = holste_geom(eta_opt=0.60)
    check("Holste eta_opt is 0.60", g2.eta_opt == 0.60, f"got {g2.eta_opt}")

    # _compute_beam_current nutzt geom.eta_opt
    chem_xe = load_chemistry(CHEM_DIR / "xenon_simple" / "chemistry.json")
    Q0 = 0.5 * SCCM_TO_PPS
    r = solve_steady_state(chem_xe, g2, 100.0, Q0, max_iter=800, tol=0.5)
    if r and r.get("converged"):
        I_060 = _compute_beam_current(r, g2)
        g3 = holste_geom(eta_opt=0.25)
        # Solve with same geom except eta_opt for beam calc
        I_025 = _compute_beam_current(r, g3)
        ratio = I_060 / I_025 if I_025 > 0 else 0
        check("eta_opt 0.60 vs 0.25 gives ~2.4x ratio",
              2.0 < ratio < 3.0,
              f"ratio={ratio:.2f}, I(0.60)={I_060:.1f}, I(0.25)={I_025:.1f}")
        check("I_beam(eta_opt=0.60) > 0", I_060 > 0, f"got {I_060:.2f}")
    else:
        check("SC solver converges for Xe/Holste", False, "not converged")


# ═══════════════════════════════════════════════════════════════
# Test 2: SC vs I-fix Konsistenz
# ═══════════════════════════════════════════════════════════════
def test_sc_ifix_consistency():
    print("\n--- Test 2: SC vs I-fix Konsistenz ---")

    chem_xe = load_chemistry(CHEM_DIR / "xenon_simple" / "chemistry.json")
    geom = holste_geom()
    Q0 = 0.5 * SCCM_TO_PPS

    # SC bei P=100W -> bestimme I_beam
    r_sc = solve_steady_state(chem_xe, geom, 100.0, Q0, max_iter=800, tol=0.5)
    if not r_sc or not r_sc.get("converged"):
        check("SC converges at P=100W", False)
        return

    I_sc = _compute_beam_current(r_sc, geom)
    check("SC I_beam > 0 at P=100W", I_sc > 0, f"I_sc={I_sc:.2f}")

    if I_sc < 1.0:
        check("SC I_beam reasonable (>1 mA)", False, f"I_sc={I_sc:.2f}")
        return

    # I-fix: suche P fuer dieses I_beam
    r_ifix, P_found, I_found, status = solve_for_target_current(
        chem_xe, geom, I_sc, Q0, P_min=1.0, P_max=500.0, tol_mA=1.0,
        max_iter=800, tol=0.5,
    )
    check("I-fix converges for I_sc target",
          status == "converged",
          f"status={status}, P={P_found:.1f}W, I={I_found:.2f}mA vs target={I_sc:.2f}mA")

    if status == "converged":
        check("I-fix P within 20% of SC P",
              abs(P_found - 100.0) / 100.0 < 0.20,
              f"P_found={P_found:.1f}W vs 100W")
        check("I-fix I_beam within 1 mA of target",
              abs(I_found - I_sc) < 1.0,
              f"I_found={I_found:.2f} vs target={I_sc:.2f}")


# ═══════════════════════════════════════════════════════════════
# Test 3: params.txt Geometrie-Laden
# ═══════════════════════════════════════════════════════════════
def test_params_txt_loading():
    print("\n--- Test 3: params.txt Geometrie-Laden ---")

    # Schreibe eine temporaere params.txt
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False,
                                     dir=str(SCRIPT_DIR)) as f:
        f.write("# Test config\n")
        f.write("R 0.05\n")
        f.write("L 0.064\n")
        f.write("betai 0.6\n")
        f.write("betag 0.20\n")
        f.write("Vgrid 1500.0\n")
        f.write("sgrid 0.004\n")
        f.write("frequency 1.1e6\n")
        f.write("Nw 6.0\n")
        f.write("R_ohm 1.0\n")
        f.write("Rc 0.055\n")
        f.write("lc 0.07\n")
        f.write("eta_opt 0.60\n")
        tmp_path = f.name

    try:
        geom = _load_geometry_from_config(tmp_path)
        check("R loaded from config", abs(geom.R - 0.05) < 1e-6, f"R={geom.R}")
        check("L loaded from config", abs(geom.L - 0.064) < 1e-6, f"L={geom.L}")
        check("betai loaded from config", abs(geom.betai - 0.6) < 1e-6, f"betai={geom.betai}")
        check("sgrid loaded from config", abs(geom.sgrid - 0.004) < 1e-6, f"sgrid={geom.sgrid}")
        check("eta_opt loaded from config", abs(geom.eta_opt - 0.60) < 1e-6, f"eta_opt={geom.eta_opt}")
        check("frequency loaded from config", abs(geom.frequency - 1.1e6) < 1e3, f"f={geom.frequency}")
    finally:
        os.unlink(tmp_path)

    # Fallback: nicht existierende Datei -> Defaults
    geom_def = _load_geometry_from_config("nonexistent_file_12345.txt")
    check("Missing config returns defaults", geom_def.R == 0.02, f"R={geom_def.R}")
    check("Missing config default eta_opt", geom_def.eta_opt == 1.0, f"eta_opt={geom_def.eta_opt}")


# ═══════════════════════════════════════════════════════════════
# Test 4: Holste I2 liefert plausible Beamstroeme
# ═══════════════════════════════════════════════════════════════
def test_holste_iodine_plausibility():
    print("\n--- Test 4: Holste-Iod Plausibilitaet ---")

    chem_i2 = load_chemistry(CHEM_DIR / "iodine_lafleur_v1" / "chemistry.json")
    geom = holste_geom()

    Q0_i2 = 0.5e-6 / M_I2  # 0.5 mg/s

    # Bei 200W sollte mit korrektem eta_opt ein messbarer Beamstrom herauskommen
    r = solve_steady_state(chem_i2, geom, 200.0, Q0_i2, max_iter=800, tol=0.5)
    if not r or not r.get("converged"):
        check("I2 solver converges at P=200W", False)
        return

    I_beam = _compute_beam_current(r, geom)
    check("I2 I_beam > 0 at P=200W",
          I_beam > 0, f"I_beam={I_beam:.2f}")
    check("I2 I_beam > 5 mA at P=200W (nicht trivial)",
          I_beam > 5, f"I_beam={I_beam:.2f}")

    # I_beam sollte mit P steigen
    I_at_50 = 0
    r50 = solve_steady_state(chem_i2, geom, 50.0, Q0_i2, max_iter=800, tol=0.5)
    if r50 and r50.get("converged"):
        I_at_50 = _compute_beam_current(r50, geom)
    check("I_beam(200W) > I_beam(50W)",
          I_beam > I_at_50,
          f"I(200)={I_beam:.2f} vs I(50)={I_at_50:.2f}")


# ═══════════════════════════════════════════════════════════════
# Test 5: Xe/I2 Nicht-Regression
# ═══════════════════════════════════════════════════════════════
def test_xe_i2_nonregression():
    print("\n--- Test 5: Xe/I2 Vergleich Nicht-Regression ---")

    chem_xe = load_chemistry(CHEM_DIR / "xenon_simple" / "chemistry.json")
    chem_i2 = load_chemistry(CHEM_DIR / "iodine_lafleur_v1" / "chemistry.json")
    geom = holste_geom()

    Q0_xe = 0.5e-6 / M_Xe
    Q0_i2 = 0.5e-6 / M_I2
    P = 100.0

    r_xe = solve_steady_state(chem_xe, geom, P, Q0_xe, max_iter=800, tol=0.5)
    r_i2 = solve_steady_state(chem_i2, geom, P, Q0_i2, max_iter=800, tol=0.5)

    I_xe = _compute_beam_current(r_xe, geom) if r_xe and r_xe.get("converged") else 0
    I_i2 = _compute_beam_current(r_i2, geom) if r_i2 and r_i2.get("converged") else 0

    check("Xe solver converges", r_xe and r_xe.get("converged"))
    check("I2 solver converges", r_i2 and r_i2.get("converged"))
    check("Xe I_beam > 0", I_xe > 0, f"I_xe={I_xe:.2f}")
    check("I2 I_beam > 0", I_i2 > 0, f"I_i2={I_i2:.2f}")

    # Qualitativ: beide sollten messbare Stroeme liefern
    if I_xe > 0 and I_i2 > 0:
        ratio = I_i2 / I_xe
        check("I2/Xe ratio plausible (0.3-3.0)",
              0.3 < ratio < 3.0,
              f"ratio={ratio:.2f}, I_Xe={I_xe:.1f}, I_I2={I_i2:.1f}")


# ═══════════════════════════════════════════════════════════════
# Test 6: beam_extraction respektiert eta_opt
# ═══════════════════════════════════════════════════════════════
def test_extraction_eta_opt():
    print("\n--- Test 6: beam_extraction eta_opt ---")

    ion_dens = {"I+": 1e16, "I2+": 5e15}
    ion_mass = {"I+": M_I, "I2+": M_I2}

    ex1 = compute_extraction(
        ion_dens, ion_mass, Te_eV=4.0, n_neutral_total=1e19, sigma_i=1e-18,
        R=0.05, L=0.064, betai=0.6, Vgrid=1500, sgrid=0.004, eta_opt=0.25,
    )
    ex2 = compute_extraction(
        ion_dens, ion_mass, Te_eV=4.0, n_neutral_total=1e19, sigma_i=1e-18,
        R=0.05, L=0.064, betai=0.6, Vgrid=1500, sgrid=0.004, eta_opt=0.60,
    )

    check("eta_opt stored in result (0.25)", abs(ex1.eta_opt - 0.25) < 1e-6)
    check("eta_opt stored in result (0.60)", abs(ex2.eta_opt - 0.60) < 1e-6)

    ratio = ex2.I_beam_mA / ex1.I_beam_mA if ex1.I_beam_mA > 0 else 0
    check("I_beam scales with eta_opt",
          abs(ratio - 0.60 / 0.25) < 0.01,
          f"ratio={ratio:.3f}, expected={0.60/0.25:.3f}")


def main():
    global passed, failed

    test_eta_opt_propagation()
    test_sc_ifix_consistency()
    test_params_txt_loading()
    test_holste_iodine_plausibility()
    test_xe_i2_nonregression()
    test_extraction_eta_opt()

    print(f"\n{'='*60}")
    print(f"  Ergebnis: {passed} passed, {failed} failed")
    if errors:
        print(f"  Fehlgeschlagen: {', '.join(errors)}")
    print(f"{'='*60}")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
