#!/usr/bin/env python3
"""
test_iodine_e2e.py -- End-to-End Validation: Iodine direct vs RF-coupled mode.

Tests the complete chain: RunConfig -> Solver -> stdout -> JSONL log -> consistency.

Test 1: Iodine direct mode (P = P_abs)
Test 2: Iodine RF-coupled mode (P = P_RFG, P_abs derived)
Test 3: Cross-check P_abs, P_RFG, zeta consistency
Test 4: GUI label verification (offscreen)
Test 5: Log metadata consistency
"""
from __future__ import annotations
import sys, os, json, subprocess, glob
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
        print(f"  FAIL: {name} -- {detail}")
        errors.append(name)


def cleanup():
    for pat in ["run_config.json", "params.txt", "python_run_*.jsonl", "simulation_log_*.txt"]:
        for f in glob.glob(str(SCRIPT_DIR / pat)):
            Path(f).unlink(missing_ok=True)


# ═══════════════════════════════════════════════════════════════
# Shared config
# ═══════════════════════════════════════════════════════════════
CHEM = "chemistry/iodine_lafleur_v1/chemistry.json"
GRONDEIN = {
    "R": 0.06, "L": 0.10, "betai": 0.7, "betag": 0.3,
    "Vgrid": 1000, "sgrid": 0.0015, "frequency": 13.56e6,
    "Nw": 5, "R_ohm": 2.0, "Rc": 0.07, "lc": 0.10, "eta_opt": 1.0,
}


def make_config(rf_coupling=False):
    from run_config import RunConfig, GeometryConfig, GridConfig, CoilConfig, OperationConfig, SweepConfig, MetaConfig
    cfg = RunConfig()
    cfg.geometry = GeometryConfig(R=0.06, L=0.10, betai=0.7, betag=0.3)
    cfg.grid = GridConfig(Vgrid=1000, sgrid=0.0015, eta_opt=1.0)
    cfg.coil = CoilConfig(frequency=13.56e6, Nw=5, R_ohm=2.0, Rc=0.07, lc=0.10)
    cfg.operation = OperationConfig(solve_mode=2, P_max=400.0, rf_coupling=rf_coupling)
    cfg.sweep = SweepConfig(Q0_start=2.0, Q0_step=1.0, N=2)
    cfg.meta = MetaConfig(gas="iodine", cs_database="", preset_id="grondein")
    return cfg


def run_solver(cfg):
    """Write config, run solver, return stdout + JSONL log content."""
    cleanup()
    cfg.save_json(SCRIPT_DIR / "run_config.json")
    cfg.to_params_txt(SCRIPT_DIR / "params.txt")
    cmd = [PYTHON, str(SCRIPT_DIR / "generic_solver.py"), str(SCRIPT_DIR / CHEM)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300, cwd=str(SCRIPT_DIR))
    stdout = r.stdout

    # Find JSONL log
    jsonl_files = sorted(glob.glob(str(SCRIPT_DIR / "python_run_*.jsonl")))
    log_records = []
    if jsonl_files:
        for line in Path(jsonl_files[-1]).read_text().strip().split("\n"):
            if line.strip():
                log_records.append(json.loads(line))

    return stdout, log_records


def parse_stdout(stdout):
    """Extract structured data from solver stdout."""
    results = []
    rf_diags = []
    ifix_results = []
    meta_lines = []

    for line in stdout.strip().split("\n"):
        p = line.split()
        if not p:
            continue
        tag = p[0]
        if tag == "RESULT" and len(p) >= 7:
            results.append({
                "ne": float(p[1]), "ng": float(p[2]),
                "Te": float(p[3]), "Tg": float(p[4]),
                "I_beam": float(p[5]), "P": float(p[6]),
            })
        elif tag == "RF_DIAG":
            diag = {}
            for part in p[1:]:
                if "=" in part:
                    k, v = part.split("=", 1)
                    diag[k] = float(v)
            rf_diags.append(diag)
        elif tag == "RF_COUPLING":
            meta_lines.append(line)
        elif tag == "IFIX_RESULT" and len(p) >= 7:
            ifix_results.append({"P": float(p[2]), "status": p[6]})

    return {"results": results, "rf_diags": rf_diags, "meta": meta_lines}


# ═══════════════════════════════════════════════════════════════
# Test 1: Iodine Direct Mode
# ═══════════════════════════════════════════════════════════════
def test_direct_mode():
    print("\n--- Test 1: Iodine Direct Mode (P = P_abs) ---")
    cfg = make_config(rf_coupling=False)
    stdout, log = run_solver(cfg)

    data = parse_stdout(stdout)

    # Solver produces results
    check("direct: has RESULT lines", len(data["results"]) > 0,
          f"found {len(data['results'])}")

    # No RF_DIAG in direct mode
    check("direct: no RF_DIAG", len(data["rf_diags"]) == 0,
          f"found {len(data['rf_diags'])} RF_DIAG lines")

    # No RF_COUPLING announcement
    check("direct: no RF_COUPLING line", len(data["meta"]) == 0)

    if data["results"]:
        r = data["results"][0]
        # P in RESULT is P_abs (= input P = 400W)
        check("direct: P = P_max (400)", abs(r["P"] - 400) < 1,
              f"P={r['P']}")
        check("direct: Te > 0", r["Te"] > 2)
        check("direct: ne > 0", r["ne"] > 1e15)

    # Log metadata
    meta_rec = next((r for r in log if r.get("type") == "metadata"), {})
    check("direct: log has metadata", bool(meta_rec))
    if meta_rec:
        check("direct: log shows config_source",
              "config_source" in meta_rec,
              str(meta_rec.keys()))

    # Log params should NOT mention rf_coupling=True
    params_rec = next((r for r in log if r.get("type") == "params"), {})
    if params_rec:
        rf_flag = params_rec.get("rf_coupling", False)
        check("direct: log rf_coupling=False", rf_flag == False,
              f"rf_coupling={rf_flag}")


# ═══════════════════════════════════════════════════════════════
# Test 2: Iodine RF-Coupled Mode
# ═══════════════════════════════════════════════════════════════
def test_rf_coupled_mode():
    print("\n--- Test 2: Iodine RF-Coupled Mode (P = P_RFG) ---")
    cfg = make_config(rf_coupling=True)
    stdout, log = run_solver(cfg)

    data = parse_stdout(stdout)

    # Solver produces results
    check("rf: has RESULT lines", len(data["results"]) > 0,
          f"found {len(data['results'])}")

    # RF_DIAG present
    check("rf: has RF_DIAG", len(data["rf_diags"]) > 0,
          f"found {len(data['rf_diags'])} RF_DIAG lines")

    # RF_COUPLING announcement present
    check("rf: has RF_COUPLING line", len(data["meta"]) > 0)

    if data["rf_diags"]:
        d = data["rf_diags"][0]
        check("rf: RF_DIAG has P_RFG", "P_RFG" in d, str(d.keys()))
        check("rf: RF_DIAG has P_abs", "P_abs" in d)
        check("rf: RF_DIAG has zeta", "zeta" in d)

        if "P_RFG" in d and "P_abs" in d and "zeta" in d:
            P_RFG = d["P_RFG"]
            P_abs = d["P_abs"]
            zeta = d["zeta"]

            check("rf: P_RFG = 400", abs(P_RFG - 400) < 1, f"P_RFG={P_RFG}")
            check("rf: P_abs < P_RFG", P_abs < P_RFG,
                  f"P_abs={P_abs:.1f} vs P_RFG={P_RFG:.1f}")
            check("rf: 0 < zeta < 1", 0 < zeta < 1, f"zeta={zeta}")
            check("rf: P_abs ~ zeta * P_RFG",
                  abs(P_abs - zeta * P_RFG) / P_RFG < 0.05,
                  f"P_abs={P_abs:.1f}, zeta*P_RFG={zeta*P_RFG:.1f}")

            print(f"  INFO: P_RFG={P_RFG:.1f}W, P_abs={P_abs:.1f}W, zeta={zeta:.4f}")

    # Log metadata shows rf_coupling
    params_rec = next((r for r in log if r.get("type") == "params"), {})
    if params_rec:
        rf_flag = params_rec.get("rf_coupling", False)
        check("rf: log rf_coupling=True", rf_flag == True,
              f"rf_coupling={rf_flag}")


# ═══════════════════════════════════════════════════════════════
# Test 3: Direct vs RF ne comparison
# ═══════════════════════════════════════════════════════════════
def test_direct_vs_rf():
    print("\n--- Test 3: Direct vs RF-coupled comparison ---")
    cfg_d = make_config(rf_coupling=False)
    stdout_d, _ = run_solver(cfg_d)
    data_d = parse_stdout(stdout_d)

    cfg_r = make_config(rf_coupling=True)
    stdout_r, _ = run_solver(cfg_r)
    data_r = parse_stdout(stdout_r)

    if data_d["results"] and data_r["results"]:
        ne_d = data_d["results"][0]["ne"]
        ne_r = data_r["results"][0]["ne"]
        Te_d = data_d["results"][0]["Te"]
        Te_r = data_r["results"][0]["Te"]

        check("comparison: RF ne < direct ne",
              ne_r < ne_d,
              f"ne_rf={ne_r:.2e}, ne_direct={ne_d:.2e}")
        print(f"  INFO: ne ratio = {ne_r/ne_d:.3f}")
        print(f"  INFO: Te direct={Te_d:.2f}eV, Te rf={Te_r:.2f}eV")


# ═══════════════════════════════════════════════════════════════
# Test 4: GUI Label Verification
# ═══════════════════════════════════════════════════════════════
def test_gui_labels():
    print("\n--- Test 4: GUI Label Verification ---")
    os.environ['QT_QPA_PLATFORM'] = 'offscreen'
    from PyQt6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    from gui import SimulatorWindow
    w = SimulatorWindow()

    # Find iodine package
    iod_idx = -1
    for i in range(w.cmb_package.count()):
        if 'iodine' in (w.cmb_package.itemData(i) or '').lower():
            iod_idx = i
            break

    if iod_idx >= 0:
        w.cmb_package.setCurrentIndex(iod_idx)

        # Direct mode (default)
        check("gui: iodine P label = P_abs", w._p_label.text() == "P_abs",
              f"label='{w._p_label.text()}'")

        # RF-coupled mode
        w._chk_rf_coupling.setChecked(True)
        check("gui: rf P label = P_RFG", w._p_label.text() == "P_RFG",
              f"label='{w._p_label.text()}'")

        # Check metric cards exist
        check("gui: P_abs metric card exists", "P_abs" in w.metric_cards)
        check("gui: zeta metric card exists", "zeta" in w.metric_cards)

        # Back to direct
        w._chk_rf_coupling.setChecked(False)
        check("gui: back to P_abs", w._p_label.text() == "P_abs")
    else:
        check("gui: iodine package found", False, "no iodine package")

    # Xenon should always be P_RFG
    for i in range(w.cmb_package.count()):
        if 'biagi' in (w.cmb_package.itemData(i) or '').lower():
            w.cmb_package.setCurrentIndex(i)
            break
    check("gui: xenon P label = P_RFG", w._p_label.text() == "P_RFG",
          f"label='{w._p_label.text()}'")


# ═══════════════════════════════════════════════════════════════
# Test 5: Log metadata consistency
# ═══════════════════════════════════════════════════════════════
def test_log_metadata():
    print("\n--- Test 5: Log Metadata Consistency ---")

    # RF mode run
    cfg = make_config(rf_coupling=True)
    stdout, log = run_solver(cfg)

    meta = next((r for r in log if r.get("type") == "metadata"), {})
    params = next((r for r in log if r.get("type") == "params"), {})
    points = [r for r in log if r.get("type") == "point"]

    check("log: metadata present", bool(meta))
    check("log: params present", bool(params))
    check("log: points present", len(points) > 0)

    if params:
        check("log: rf_coupling in params", "rf_coupling" in params,
              str(list(params.keys())[:10]))

    # Verify RESULT P matches expectations
    data = parse_stdout(stdout)
    if data["results"] and data["rf_diags"]:
        result_P = data["results"][0]["P"]
        rf_P_RFG = data["rf_diags"][0].get("P_RFG", 0)
        # In RF mode, RESULT P should be P_RFG (the input power)
        check("log: RESULT P = P_RFG",
              abs(result_P - rf_P_RFG) < 1,
              f"RESULT P={result_P}, RF_DIAG P_RFG={rf_P_RFG}")


def main():
    global passed, failed
    test_direct_mode()
    test_rf_coupled_mode()
    test_direct_vs_rf()
    test_gui_labels()
    test_log_metadata()
    cleanup()

    print(f"\n{'='*60}")
    print(f"  Result: {passed} passed, {failed} failed")
    if errors:
        print(f"  Failed:")
        for e in errors:
            print(f"    - {e}")
    print(f"{'='*60}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
