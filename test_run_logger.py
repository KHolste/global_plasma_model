#!/usr/bin/env python3
"""
test_run_logger.py -- Tests fuer das strukturierte Python-Run-Log.

Prueft:
  1. RunLogger erzeugt JSONL und Masterlog
  2. JSONL-Format ist gueltig und hat Pflichtfelder
  3. Masterlog hat C++-kompatible Sektionen
  4. Erfolgreiche und fehlgeschlagene Punkte sind korrekt geloggt
  5. Integration: Solver-Lauf erzeugt Log-Dateien mit konsistenten Daten

Ausfuehrung:
    python test_run_logger.py
"""
from __future__ import annotations
import sys, os, json, subprocess, glob
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
os.chdir(SCRIPT_DIR)
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
        print(f"  FAIL: {name} -- {detail}")
        errors.append(name)


# ═════════════════════════════════════════════════════════════
# Test 1: RunLogger Grundfunktion
# ═════════════════════════════════════════════════════════════
def test_logger_basics():
    print("\n--- Test 1: RunLogger Grundfunktion ---")
    from run_logger import RunLogger
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        lg = RunLogger(output_dir=td)
        lg.set_metadata(backend="python", package="test_pkg", solve_mode="SC")
        lg.set_params({"R": 0.02, "L": 0.04, "P": 18.0})

        lg.add_point(idx=0, Q0_sccm=0.4, status="converged", P_W=18.0,
                     I_beam_mA=15.2, Te_eV=3.8, Tg_K=300, ne=1e17, ng=2e19,
                     converged=True, merit=1e-4, iterations=25)
        lg.add_point(idx=1, Q0_sccm=0.5, status="above_P_max", P_W=80.0,
                     I_beam_mA=5.6, I_target_mA=15.0, delta_I_mA=-9.4,
                     diag="Sollstrom nicht erreichbar")
        lg.add_event(0, 0.4, "Test-Event")

        jsonl_path, txt_path = lg.finalize(elapsed=1.23, count_ok=1, count_fail=1)

        # JSONL existiert und ist gueltig
        check("JSONL file exists", Path(jsonl_path).exists())
        check("Masterlog file exists", Path(txt_path).exists())

        # JSONL parsen
        lines = Path(jsonl_path).read_text().strip().split("\n")
        check("JSONL has >= 5 lines", len(lines) >= 5, f"found {len(lines)}")

        records = [json.loads(l) for l in lines]
        types = [r["type"] for r in records]
        check("JSONL has metadata", "metadata" in types)
        check("JSONL has params", "params" in types)
        check("JSONL has points", types.count("point") == 2)
        check("JSONL has event", "event" in types)
        check("JSONL has summary", "summary" in types)

        # Pflichtfelder im Punkt
        pt = next(r for r in records if r["type"] == "point" and r["status"] == "converged")
        for key in ("Q0_sccm", "P_W", "I_beam_mA", "Te_eV", "Tg_K", "ne", "ng", "status"):
            check(f"point has {key}", key in pt, f"keys={list(pt.keys())}")

        # Fehlpunkt hat Diagnose
        fp = next(r for r in records if r["type"] == "point" and r["status"] == "above_P_max")
        check("fail point has diag", fp.get("diag") != "")
        check("fail point has I_target", fp.get("I_target_mA") == 15.0)

        # Summary
        sm = next(r for r in records if r["type"] == "summary")
        check("summary converged=1", sm["converged"] == 1)
        check("summary failed=1", sm["failed"] == 1)

        # Masterlog hat DATA-Zeilen
        txt = Path(txt_path).read_text()
        check("masterlog has DATA_HEADER", "DATA_HEADER|" in txt)
        check("masterlog has DATA|", txt.count("DATA|") >= 2,
              f"found {txt.count('DATA|')} DATA lines")
        check("masterlog has SUMMARY", "SUMMARY" in txt)
        check("masterlog has RESULT TABLE", "RESULT TABLE" in txt)


# ═════════════════════════════════════════════════════════════
# Test 2: Integration — Xenon I-fix erzeugt Log
# ═════════════════════════════════════════════════════════════
def test_xenon_integration():
    print("\n--- Test 2: Xenon I-fix erzeugt Log ---")
    chem = SCRIPT_DIR / "chemistry" / "xenon_simple" / "chemistry.json"
    if not chem.exists():
        check("chemistry exists", False)
        return

    # params.txt schreiben
    with open(SCRIPT_DIR / "params.txt", "w") as f:
        f.write("R 0.02\nL 0.04\nbetai 0.5\nbetag 0.05145\n")
        f.write("frequency 2.5e6\nNw 6\nR_ohm 0.36\nRc 0.02\nlc 0.04\n")
        f.write("Vgrid 1500\nsgrid 0.001\neta_opt 0.30\nP_RFG_max 500\n")

    # Alte Logs aufraeumen
    for old in glob.glob(str(SCRIPT_DIR / "python_run_*.jsonl")):
        os.unlink(old)
    for old in glob.glob(str(SCRIPT_DIR / "simulation_log_*.txt")):
        os.unlink(old)

    cmd = [PYTHON, str(SCRIPT_DIR / "generic_solver.py"),
           str(chem), "1", "3", "0.4", "0.1", "2"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=120,
                       cwd=str(SCRIPT_DIR))

    # Log-Dateien suchen
    jsonl_files = sorted(glob.glob(str(SCRIPT_DIR / "python_run_*.jsonl")))
    txt_files = sorted(glob.glob(str(SCRIPT_DIR / "simulation_log_*.txt")))

    check("JSONL log created", len(jsonl_files) > 0, f"found {len(jsonl_files)}")
    check("Masterlog created", len(txt_files) > 0, f"found {len(txt_files)}")

    if jsonl_files:
        lines = Path(jsonl_files[-1]).read_text().strip().split("\n")
        records = [json.loads(l) for l in lines]
        points = [r for r in records if r["type"] == "point"]
        check("JSONL has 2 points", len(points) == 2, f"found {len(points)}")

        # Mindestens einer sollte Status haben
        statuses = [p["status"] for p in points]
        check("points have status", all(s != "" for s in statuses),
              f"statuses={statuses}")

        # Metadata
        meta = next((r for r in records if r["type"] == "metadata"), {})
        check("metadata has backend=python", meta.get("backend") == "python")
        check("metadata has package", "package" in meta)

        # Params
        params = next((r for r in records if r["type"] == "params"), {})
        check("params has R", "R" in params)
        check("params has eta_opt", "eta_opt" in params)

    # Aufraeumen
    for f in jsonl_files + txt_files:
        os.unlink(f)
    (SCRIPT_DIR / "params.txt").unlink(missing_ok=True)


# ═════════════════════════════════════════════════════════════
# Test 3: Integration — Iod SC erzeugt Log mit Chemie-Feldern
# ═════════════════════════════════════════════════════════════
def test_iodine_integration():
    print("\n--- Test 3: Iod SC erzeugt Log mit Chemie ---")
    chem = SCRIPT_DIR / "chemistry" / "iodine_lafleur_v1" / "chemistry.json"
    if not chem.exists():
        check("iodine chemistry exists", False)
        return

    with open(SCRIPT_DIR / "params.txt", "w") as f:
        f.write("R 0.06\nL 0.10\nbetai 0.7\nbetag 0.3\n")
        f.write("frequency 13.56e6\nNw 5\nR_ohm 2.0\nRc 0.07\nlc 0.10\n")
        f.write("Vgrid 1000\nsgrid 0.0015\neta_opt 0.40\n")

    for old in glob.glob(str(SCRIPT_DIR / "python_run_*.jsonl")):
        os.unlink(old)

    cmd = [PYTHON, str(SCRIPT_DIR / "generic_solver.py"),
           str(chem), "2", "400", "2.0", "0.5", "2"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=120,
                       cwd=str(SCRIPT_DIR))

    jsonl_files = sorted(glob.glob(str(SCRIPT_DIR / "python_run_*.jsonl")))
    check("JSONL log created", len(jsonl_files) > 0)

    if jsonl_files:
        lines = Path(jsonl_files[-1]).read_text().strip().split("\n")
        records = [json.loads(l) for l in lines]
        points = [r for r in records if r["type"] == "point"]
        check("JSONL has points", len(points) > 0)

        if points:
            pt = points[0]
            check("point has diss", pt.get("diss", 0) > 0,
                  f"diss={pt.get('diss')}")
            check("point has fIp", pt.get("fIp", 0) > 0,
                  f"fIp={pt.get('fIp')}")
            check("point has alpha", "alpha" in pt)
            check("point has Te_eV > 0", pt.get("Te_eV", 0) > 0)

    # Aufraeumen
    for f in jsonl_files:
        os.unlink(f)
    for f in glob.glob(str(SCRIPT_DIR / "simulation_log_*.txt")):
        os.unlink(f)
    (SCRIPT_DIR / "params.txt").unlink(missing_ok=True)


# ═════════════════════════════════════════════════════════════
# Test 4: JSONL/stdout Konsistenz
# ═════════════════════════════════════════════════════════════
def test_stdout_jsonl_consistency():
    print("\n--- Test 4: stdout/JSONL Konsistenz ---")
    chem = SCRIPT_DIR / "chemistry" / "xenon_simple" / "chemistry.json"
    if not chem.exists():
        check("chemistry exists", False)
        return

    with open(SCRIPT_DIR / "params.txt", "w") as f:
        f.write("R 0.02\nL 0.04\nbetai 0.5\nbetag 0.05145\n")
        f.write("frequency 2.5e6\nNw 6\nR_ohm 0.36\nRc 0.02\nlc 0.04\n")
        f.write("Vgrid 1500\nsgrid 0.001\neta_opt 0.30\nP_RFG_max 500\n")

    for old in glob.glob(str(SCRIPT_DIR / "python_run_*.jsonl")):
        os.unlink(old)

    cmd = [PYTHON, str(SCRIPT_DIR / "generic_solver.py"),
           str(chem), "1", "3", "0.5", "0.1", "2"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=120,
                       cwd=str(SCRIPT_DIR))

    # stdout RESULT parsen
    stdout_results = []
    for line in r.stdout.strip().split("\n"):
        parts = line.split()
        if parts and parts[0] == "RESULT" and len(parts) >= 7:
            stdout_results.append({
                "ne": float(parts[1]), "I_beam": float(parts[5]), "P": float(parts[6])
            })

    # JSONL parsen
    jsonl_files = sorted(glob.glob(str(SCRIPT_DIR / "python_run_*.jsonl")))
    if jsonl_files:
        records = [json.loads(l) for l in Path(jsonl_files[-1]).read_text().strip().split("\n")]
        log_points = [r for r in records if r["type"] == "point" and r.get("converged")]

        if stdout_results and log_points:
            for i in range(min(len(stdout_results), len(log_points))):
                sr = stdout_results[i]
                lp = log_points[i]
                check(f"consistency[{i}]: I_beam match",
                      abs(sr["I_beam"] - lp["I_beam_mA"]) < 0.1,
                      f"stdout={sr['I_beam']}, log={lp['I_beam_mA']}")
                check(f"consistency[{i}]: P match",
                      abs(sr["P"] - lp["P_W"]) < 0.2,
                      f"stdout={sr['P']}, log={lp['P_W']}")

    # Aufraeumen
    for f in jsonl_files:
        os.unlink(f)
    for f in glob.glob(str(SCRIPT_DIR / "simulation_log_*.txt")):
        os.unlink(f)
    (SCRIPT_DIR / "params.txt").unlink(missing_ok=True)


def main():
    global passed, failed
    test_logger_basics()
    test_xenon_integration()
    test_iodine_integration()
    test_stdout_jsonl_consistency()

    print(f"\n{'='*60}")
    print(f"  Ergebnis: {passed} passed, {failed} failed")
    if errors:
        print(f"  Fehlgeschlagen: {', '.join(errors)}")
    print(f"{'='*60}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
