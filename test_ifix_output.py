#!/usr/bin/env python3
"""
test_ifix_output.py -- Tests fuer I-fix Ergebnisausgabe und GUI-Parsing.

Prueft:
  1. IFIX_RESULT Zeile wird emittiert (Erfolg + Fehler)
  2. Format ist korrekt und parsebar
  3. delta_I-Vorzeichen ist konsistent
  4. Fehlpunkte enthalten gueltige Best-Attempt-Daten
  5. Bestehende RESULT/IODINE_EXT Zeilen regressieren nicht

Ausfuehrung:
    python test_ifix_output.py
"""
from __future__ import annotations
import sys, subprocess, re
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

PYTHON = sys.executable

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


def run_ifix_sweep(chem_json, I_soll, Q0_start, Q0_step, N, P_max=500.0):
    """Fuehre I-fix-Sweep via generic_solver.py aus und gib stdout zurueck."""
    # Erst params.txt schreiben damit Geometrie geladen wird
    params_path = SCRIPT_DIR / "test_ifix_params.txt"
    with open(params_path, "w") as f:
        f.write("R 0.02\nL 0.04\nbetai 0.5\nbetag 0.05145\n")
        f.write("frequency 2.5e6\nNw 6\nR_ohm 0.36\nRc 0.02\nlc 0.04\n")
        f.write("Vgrid 1500\nsgrid 0.001\neta_opt 1.0\n")
        f.write(f"P_RFG_max {P_max}\n")

    cmd = [PYTHON, str(SCRIPT_DIR / "generic_solver.py"),
           str(chem_json), "1", str(I_soll),
           str(Q0_start), str(Q0_step), str(N)]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120,
                            cwd=str(SCRIPT_DIR))
    # Raeumen
    params_path.unlink(missing_ok=True)
    return result.stdout


# ═════════════════════════════════════════════════════════════
# Test 1: IFIX_RESULT wird bei Konvergenz emittiert
# ═════════════════════════════════════════════════════════════
def test_ifix_result_on_success():
    print("\n--- Test 1: IFIX_RESULT bei Konvergenz ---")
    chem = SCRIPT_DIR / "chemistry" / "xenon_simple" / "chemistry.json"
    if not chem.exists():
        check("xenon_simple exists", False, "chemistry.json nicht gefunden")
        return

    # I_soll=10 mA mit Default-Geometrie (R=2cm, eta_opt=1.0)
    out = run_ifix_sweep(chem, I_soll=10, Q0_start=0.4, Q0_step=0.1, N=2)
    lines = out.strip().split("\n")

    ifix_lines = [l for l in lines if l.startswith("IFIX_RESULT")]
    check("IFIX_RESULT emitted", len(ifix_lines) > 0, f"found {len(ifix_lines)}")

    if ifix_lines:
        # Parse erste Zeile
        parts = ifix_lines[0].split()
        check("IFIX_RESULT has 7 fields", len(parts) >= 7, f"parts={len(parts)}")
        if len(parts) >= 7:
            q0 = float(parts[1])
            P = float(parts[2])
            I_target = float(parts[3])
            I_beam = float(parts[4])
            delta_I = float(parts[5])
            status = parts[6]

            check("Q0 > 0", q0 > 0, f"Q0={q0}")
            check("P > 0", P > 0, f"P={P}")
            check("I_target = 10", abs(I_target - 10) < 0.1, f"I_target={I_target}")
            check("I_beam > 0", I_beam > 0, f"I_beam={I_beam}")
            check("delta_I = I_beam - I_target", abs(delta_I - (I_beam - I_target)) < 0.01,
                  f"delta_I={delta_I}, I_beam-I_target={I_beam-I_target}")
            check("status is converged", status == "converged", f"status={status}")

    # RESULT sollte auch noch kommen
    result_lines = [l for l in lines if l.startswith("RESULT")]
    check("RESULT also emitted", len(result_lines) > 0)


# ═════════════════════════════════════════════════════════════
# Test 2: IFIX_RESULT wird bei Fehler (above_P_max) emittiert
# ═════════════════════════════════════════════════════════════
def test_ifix_result_on_failure():
    print("\n--- Test 2: IFIX_RESULT bei Fehler ---")
    chem = SCRIPT_DIR / "chemistry" / "xenon_simple" / "chemistry.json"
    if not chem.exists():
        check("xenon_simple exists", False)
        return

    # Unplausibel hoher Sollstrom -> above_P_max
    out = run_ifix_sweep(chem, I_soll=999, Q0_start=0.4, Q0_step=0.1, N=1, P_max=50)
    lines = out.strip().split("\n")

    ifix_lines = [l for l in lines if l.startswith("IFIX_RESULT")]
    check("IFIX_RESULT emitted on fail", len(ifix_lines) > 0, f"found {len(ifix_lines)}")

    if ifix_lines:
        parts = ifix_lines[0].split()
        if len(parts) >= 7:
            I_target = float(parts[3])
            I_beam = float(parts[4])
            delta_I = float(parts[5])
            status = parts[6]

            check("I_target = 999", abs(I_target - 999) < 1, f"I_target={I_target}")
            check("I_beam < I_target", I_beam < I_target,
                  f"I_beam={I_beam}, I_target={I_target}")
            check("delta_I is negative", delta_I < 0,
                  f"delta_I={delta_I}")
            check("status is not converged", status != "converged", f"status={status}")

    # SOLVER_FAIL sollte auch kommen
    fail_lines = [l for l in lines if l.startswith("SOLVER_FAIL")]
    check("SOLVER_FAIL also emitted", len(fail_lines) > 0)

    # SOLVER_FAIL enthaelt diag=
    if fail_lines:
        check("SOLVER_FAIL has diag=", "diag=" in fail_lines[0])


# ═════════════════════════════════════════════════════════════
# Test 3: Konsistenz zwischen IFIX_RESULT und RESULT
# ═════════════════════════════════════════════════════════════
def test_ifix_result_consistency():
    print("\n--- Test 3: IFIX_RESULT/RESULT Konsistenz ---")
    chem = SCRIPT_DIR / "chemistry" / "xenon_simple" / "chemistry.json"
    if not chem.exists():
        check("xenon_simple exists", False)
        return

    # I_soll=10 mA erreichbar mit Default-Geometrie (eta_opt=1.0)
    out = run_ifix_sweep(chem, I_soll=10, Q0_start=0.4, Q0_step=0.1, N=2)
    lines = out.strip().split("\n")

    ifix_lines = [l for l in lines if l.startswith("IFIX_RESULT")]
    result_lines = [l for l in lines if l.startswith("RESULT")]

    converged_ifix = [l for l in ifix_lines if l.split()[6] == "converged"]
    check("Some points converged", len(converged_ifix) > 0,
          f"converged={len(converged_ifix)}/{len(ifix_lines)}")

    # Fuer konvergierte Punkte: I_beam in IFIX_RESULT ~ I_beam in RESULT
    if converged_ifix and result_lines:
        ifix_I = float(converged_ifix[0].split()[4])
        result_I = float(result_lines[0].split()[5])  # 5. Feld in RESULT
        check("I_beam consistent (IFIX vs RESULT)",
              abs(ifix_I - result_I) < 1.0,
              f"IFIX={ifix_I:.2f}, RESULT={result_I:.2f}")


# ═════════════════════════════════════════════════════════════
# Test 4: Iod-Workflow regressiert nicht
# ═════════════════════════════════════════════════════════════
def test_iodine_ifix():
    print("\n--- Test 4: Iod-I-fix ---")
    chem = SCRIPT_DIR / "chemistry" / "iodine_lafleur_v1" / "chemistry.json"
    if not chem.exists():
        check("iodine_lafleur_v1 exists", False)
        return

    out = run_ifix_sweep(chem, I_soll=8, Q0_start=0.5, Q0_step=0.2, N=2, P_max=100)
    lines = out.strip().split("\n")

    ifix_lines = [l for l in lines if l.startswith("IFIX_RESULT")]
    check("Iodine: IFIX_RESULT emitted", len(ifix_lines) > 0)

    # IODINE_EXT sollte weiterhin kommen
    iod_ext = [l for l in lines if l.startswith("IODINE_EXT")]
    result_lines = [l for l in lines if l.startswith("RESULT")]
    if result_lines:
        check("Iodine: IODINE_EXT emitted", len(iod_ext) > 0,
              f"RESULT={len(result_lines)}, IODINE_EXT={len(iod_ext)}")


# ═════════════════════════════════════════════════════════════
# Test 5: GUI-Parsing-Logik simulieren
# ═════════════════════════════════════════════════════════════
def test_gui_parsing_logic():
    print("\n--- Test 5: GUI-Parsing-Logik ---")

    # Simuliere eine IFIX_RESULT Zeile
    line = "IFIX_RESULT 0.5000 14.32 12.00 12.05 0.05 converged"
    parts = line.split()
    tag = parts[0]

    check("tag is IFIX_RESULT", tag == "IFIX_RESULT")
    check("7 fields", len(parts) >= 7)

    q0 = float(parts[1])
    P = float(parts[2])
    I_target = float(parts[3])
    I_beam = float(parts[4])
    delta_I = float(parts[5])
    status = parts[6]

    check("Q0 parsed", abs(q0 - 0.5) < 0.001)
    check("P parsed", abs(P - 14.32) < 0.01)
    check("I_target parsed", abs(I_target - 12.0) < 0.01)
    check("I_beam parsed", abs(I_beam - 12.05) < 0.01)
    check("delta_I parsed", abs(delta_I - 0.05) < 0.01)
    check("status parsed", status == "converged")

    # Fehlfall
    line2 = "IFIX_RESULT 0.5000 500.00 80.00 40.50 -39.50 above_P_max"
    p2 = line2.split()
    check("fail status parsed", p2[6] == "above_P_max")
    check("fail delta_I negative", float(p2[5]) < 0)


def main():
    global passed, failed

    test_ifix_result_on_success()
    test_ifix_result_on_failure()
    test_ifix_result_consistency()
    test_iodine_ifix()
    test_gui_parsing_logic()

    print(f"\n{'='*60}")
    print(f"  Ergebnis: {passed} passed, {failed} failed")
    if errors:
        print(f"  Fehlgeschlagen: {', '.join(errors)}")
    print(f"{'='*60}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
