#!/usr/bin/env python3
"""
test_log_viewer.py -- Tests fuer den vereinheitlichten Log-Viewer.

Prueft:
  1. Formaterkennung (C++ vs Python JSONL)
  2. JSONL-Parser erzeugt korrektes Datenmodell
  3. C++-Parser funktioniert weiterhin
  4. Gemeinsame Spalten sind in beiden Formaten vorhanden
  5. Fehlpunkte haben korrekten Status

Ausfuehrung:
    python test_log_viewer.py
"""
from __future__ import annotations
import sys, json, tempfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from log_viewer import (
    parse_master_log, parse_jsonl_log, detect_log_format,
    load_any_log, get_numeric_columns, get_column_data,
)

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
# Test-Daten erzeugen
# ═════════════════════════════════════════════════════════════

def create_test_jsonl(td):
    """Erzeuge ein Test-JSONL-Log."""
    path = Path(td) / "python_run_test.jsonl"
    records = [
        {"type": "metadata", "backend": "python", "package": "iodine_lafleur_v1",
         "solve_mode": "I-fix", "timestamp_start": "2026-03-30 12:00:00"},
        {"type": "params", "R": 0.05, "L": 0.064, "eta_opt": 0.60, "P_max": 120},
        {"type": "point", "idx": 0, "Q0_sccm": 2.0, "status": "converged",
         "P_W": 71.7, "I_beam_mA": 30.2, "Te_eV": 3.72, "Tg_K": 301.0,
         "ne": 3.6e17, "ng": 2.1e19, "I_target_mA": 30.0, "delta_I_mA": 0.2,
         "diss": 0.13, "fIp": 0.45, "fI2p": 0.55, "alpha": 0.035,
         "converged": True, "merit": 1e-4, "iterations": 25},
        {"type": "point", "idx": 1, "Q0_sccm": 2.5, "status": "above_P_max",
         "P_W": 120.0, "I_beam_mA": 45.1, "Te_eV": 0, "Tg_K": 0,
         "ne": 0, "ng": 0, "I_target_mA": 80.0, "delta_I_mA": -34.9,
         "converged": False, "diag": "Sollstrom 80 mA nicht erreichbar"},
        {"type": "event", "idx": 1, "Q0_sccm": 2.5,
         "message": "above_P_max: I=45.1 bei P=120W"},
        {"type": "summary", "total_points": 2, "converged": 1, "failed": 1,
         "runtime_seconds": 5.67},
    ]
    with open(path, "w") as f:
        for r in records:
            f.write(json.dumps(r) + "\n")
    return str(path)


def create_test_masterlog(td):
    """Erzeuge ein minimales C++ Masterlog."""
    path = Path(td) / "simulation_log_test.txt"
    txt = """==================================================
GLOBAL PLASMA MODEL - SIMULATION LOG
timestamp_start: 2026-03-30 12:00:00
runtime_seconds: 3.45
solver_mode:     stationary (LM), solve_mode=1
==================================================

--------------------------------------------------
SIMULATION PARAMETERS
--------------------------------------------------
Radius                             | 2.000000e-02 | m      | R
Laenge                             | 4.000000e-02 | m      | L

# CL_LIMIT_I_mA=15.5000
# CL_LIMIT_J_A_per_m2=1.234e+02
DATA_HEADER|idx|Q0sccm|status|P_RFG_W|I_mA|Te_eV|Tg_K|n_m3|ng_m3|note
DATA|0|0.4000|CONVERGED|16.01|14.90|4.12|297.0|1.843e+17|2.384e+19|ok
DATA|1|0.4500|CONVERGED|14.40|14.92|3.82|297.9|1.953e+17|3.003e+19|ok
DATA|2|0.5000|NUMERICAL_FAIL|80.00|5.60|0|0|0|0|solver failed

--------------------------------------------------
SUMMARY
--------------------------------------------------
total_points                  | 3
converged                     | 2
numerical_fail                | 1
--------------------------------------------------
"""
    path.write_text(txt)
    return str(path)


# ═════════════════════════════════════════════════════════════
# Test 1: Formaterkennung
# ═════════════════════════════════════════════════════════════
def test_format_detection():
    print("\n--- Test 1: Formaterkennung ---")
    with tempfile.TemporaryDirectory() as td:
        jsonl_path = create_test_jsonl(td)
        cpp_path = create_test_masterlog(td)

        check("detect JSONL", detect_log_format(jsonl_path) == "jsonl")
        check("detect C++", detect_log_format(cpp_path) == "cpp")

        # Suffix-basiert
        check("detect .jsonl suffix", detect_log_format("/tmp/test.jsonl") == "jsonl")


# ═════════════════════════════════════════════════════════════
# Test 2: JSONL-Parser
# ═════════════════════════════════════════════════════════════
def test_jsonl_parser():
    print("\n--- Test 2: JSONL-Parser ---")
    with tempfile.TemporaryDirectory() as td:
        path = create_test_jsonl(td)
        parsed = parse_jsonl_log(path)

        check("has metadata", bool(parsed["metadata"]))
        check("metadata backend=python", parsed["metadata"].get("backend") == "python")
        check("metadata has package", "package" in parsed["metadata"])

        check("has params", bool(parsed["params"]))
        check("params has R", "R" in parsed["params"])
        check("params has eta_opt", "eta_opt" in parsed["params"])

        check("has 2 data points", len(parsed["data"]) == 2)
        check("has columns", len(parsed["columns"]) > 0)
        check("has summary", bool(parsed["summary"]))

        # Spalten pruefen
        cols = parsed["columns"]
        for expected in ("Q0sccm", "P_W", "I_mA", "Te_eV", "Tg_K", "status"):
            check(f"column {expected} present", expected in cols, f"cols={cols}")

        # Punkt 0: konvergiert
        pt0 = parsed["data"][0]
        check("pt0 status=CONVERGED", pt0["status"] == "CONVERGED")
        check("pt0 Q0sccm=2.0", pt0.get("Q0sccm") == 2.0)
        check("pt0 I_mA=30.2", pt0.get("I_mA") == 30.2)
        check("pt0 has diss", pt0.get("diss") == 0.13)
        check("pt0 has I_target_mA", pt0.get("I_target_mA") == 30.0)
        check("pt0 has delta_I_mA", pt0.get("delta_I_mA") == 0.2)

        # Punkt 1: Fehler
        pt1 = parsed["data"][1]
        check("pt1 status=ABOVE_P_MAX", pt1["status"] == "ABOVE_P_MAX")
        check("pt1 has note/diag", pt1.get("note") != "")

        # Events
        check("has events", len(parsed["events"]) > 0)

        # Numerische Spalten
        num_cols = get_numeric_columns(parsed)
        check("numeric columns include Q0sccm", "Q0sccm" in num_cols)
        check("numeric columns include I_mA", "I_mA" in num_cols)


# ═════════════════════════════════════════════════════════════
# Test 3: C++ Parser funktioniert weiterhin
# ═════════════════════════════════════════════════════════════
def test_cpp_parser():
    print("\n--- Test 3: C++ Parser ---")
    with tempfile.TemporaryDirectory() as td:
        path = create_test_masterlog(td)
        parsed = parse_master_log(path)

        check("has metadata", bool(parsed["metadata"]))
        check("has data", len(parsed["data"]) == 3)
        check("has columns", len(parsed["columns"]) > 0)

        # CL-Limit
        check("CL limit parsed", parsed["cl_limit_I_mA"] == 15.5)

        # Status
        check("pt0 converged", parsed["data"][0].get("status") == "CONVERGED")
        check("pt2 failed", parsed["data"][2].get("status") == "NUMERICAL_FAIL")

        # Werte
        check("pt0 Te=4.12", parsed["data"][0].get("Te_eV") == 4.12)

        # Numerische Spalten
        num = get_numeric_columns(parsed)
        check("num has Q0sccm", "Q0sccm" in num)
        check("num has I_mA", "I_mA" in num)


# ═════════════════════════════════════════════════════════════
# Test 4: Auto-Erkennung und load_any_log
# ═════════════════════════════════════════════════════════════
def test_auto_load():
    print("\n--- Test 4: Auto-Erkennung ---")
    with tempfile.TemporaryDirectory() as td:
        jsonl_path = create_test_jsonl(td)
        cpp_path = create_test_masterlog(td)

        p1 = load_any_log(jsonl_path)
        check("JSONL auto-load source=python", p1.get("source") == "python")
        check("JSONL auto-load has data", len(p1["data"]) == 2)

        p2 = load_any_log(cpp_path)
        check("C++ auto-load source=cpp", p2.get("source") == "cpp")
        check("C++ auto-load has data", len(p2["data"]) == 3)


# ═════════════════════════════════════════════════════════════
# Test 5: Gemeinsame Spalten
# ═════════════════════════════════════════════════════════════
def test_common_columns():
    print("\n--- Test 5: Gemeinsame Spalten ---")
    with tempfile.TemporaryDirectory() as td:
        jsonl_path = create_test_jsonl(td)
        cpp_path = create_test_masterlog(td)

        p_py = load_any_log(jsonl_path)
        p_cpp = load_any_log(cpp_path)

        # Beide muessen Q0sccm, I_mA, Te_eV, Tg_K, status haben
        common = ["Q0sccm", "I_mA", "Te_eV", "Tg_K", "status"]
        for col in common:
            check(f"Python has {col}", col in p_py["columns"], f"cols={p_py['columns']}")
            check(f"C++ has {col}", col in p_cpp["columns"], f"cols={p_cpp['columns']}")

        # get_column_data funktioniert fuer beide
        py_Te = get_column_data(p_py, "Te_eV", only_converged=True)
        cpp_Te = get_column_data(p_cpp, "Te_eV", only_converged=True)
        check("Python Te data", len(py_Te) > 0, f"n={len(py_Te)}")
        check("C++ Te data", len(cpp_Te) > 0, f"n={len(cpp_Te)}")

        # Python-spezifische Spalten
        check("Python has diss", "diss" in p_py["columns"])
        check("Python has delta_I_mA", "delta_I_mA" in p_py["columns"])


# ═════════════════════════════════════════════════════════════
# Test 6: Integration mit echtem Solver-Log
# ═════════════════════════════════════════════════════════════
def test_real_solver_log():
    print("\n--- Test 6: Echtes Solver-Log ---")
    import subprocess, glob

    chem = SCRIPT_DIR / "chemistry" / "xenon_simple" / "chemistry.json"
    if not chem.exists():
        check("chemistry exists", False)
        return

    with open(SCRIPT_DIR / "params.txt", "w") as f:
        f.write("R 0.02\nL 0.04\nbetai 0.5\nbetag 0.05145\n")
        f.write("frequency 2.5e6\nNw 6\nR_ohm 0.36\nRc 0.02\nlc 0.04\n")
        f.write("Vgrid 1500\nsgrid 0.001\neta_opt 0.30\nP_RFG_max 500\n")

    for old in glob.glob(str(SCRIPT_DIR / "python_run_*.jsonl")):
        Path(old).unlink()

    cmd = [sys.executable, str(SCRIPT_DIR / "generic_solver.py"),
           str(chem), "2", "18", "0.4", "0.1", "2"]
    subprocess.run(cmd, capture_output=True, text=True, timeout=120,
                   cwd=str(SCRIPT_DIR))

    jsonl_files = sorted(glob.glob(str(SCRIPT_DIR / "python_run_*.jsonl")))
    if not jsonl_files:
        check("JSONL created", False)
        return

    parsed = load_any_log(jsonl_files[-1])
    check("real log source=python", parsed.get("source") == "python")
    check("real log has data", len(parsed["data"]) > 0)

    if parsed["data"]:
        num = get_numeric_columns(parsed)
        check("real log has Te_eV", "Te_eV" in num)
        check("real log has I_mA", "I_mA" in num)

        conv = get_column_data(parsed, "Te_eV", only_converged=True)
        check("real log converged Te > 0", len(conv) > 0 and conv[0] > 0,
              f"Te={conv}")

    # Aufraeumen
    for f in jsonl_files:
        Path(f).unlink()
    for f in glob.glob(str(SCRIPT_DIR / "simulation_log_*.txt")):
        Path(f).unlink()
    (SCRIPT_DIR / "params.txt").unlink(missing_ok=True)


def main():
    global passed, failed
    test_format_detection()
    test_jsonl_parser()
    test_cpp_parser()
    test_auto_load()
    test_common_columns()
    test_real_solver_log()

    print(f"\n{'='*60}")
    print(f"  Ergebnis: {passed} passed, {failed} failed")
    if errors:
        print(f"  Fehlgeschlagen: {', '.join(errors)}")
    print(f"{'='*60}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
