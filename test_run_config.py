#!/usr/bin/env python3
"""
test_run_config.py -- Tests fuer die einheitliche Runtime-Konfiguration.

Prueft:
  1. RunConfig Defaults und Struktur
  2. params.txt Roundtrip (write -> read -> compare)
  3. Preset-Integration (Preset -> RunConfig -> params.txt)
  4. Validierung (gueltige und ungueltige Configs)
  5. Cross-Backend Feldkonsistenz (C++ und Python lesen dieselben Felder)
  6. JSON-Serialisierung und Flat-Dict
  7. Konsistenz mit bestehenden Preset-Daten

Ausfuehrung:
    python test_run_config.py
"""
from __future__ import annotations
import sys, json, tempfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from run_config import RunConfig, GeometryConfig, GridConfig, CoilConfig, OperationConfig

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
# Test 1: Defaults
# ═════════════════════════════════════════════════════════════
def test_defaults():
    print("\n--- Test 1: Defaults ---")
    cfg = RunConfig()
    check("R default", cfg.geometry.R == 0.02)
    check("eta_opt default", cfg.grid.eta_opt == 1.0)
    check("solve_mode default", cfg.operation.solve_mode == 1)
    check("I_soll default", cfg.operation.I_soll == 15.0)
    check("P_max default", cfg.operation.P_max == 80.0)
    check("gas default", cfg.meta.gas == "xenon")
    check("rate_model default", cfg.rates.rate_model == 0)


# ═════════════════════════════════════════════════════════════
# Test 2: params.txt Roundtrip
# ═════════════════════════════════════════════════════════════
def test_roundtrip():
    print("\n--- Test 2: params.txt Roundtrip ---")
    cfg = RunConfig()
    cfg.geometry.R = 0.05
    cfg.geometry.L = 0.064
    cfg.grid.eta_opt = 0.60
    cfg.grid.Vgrid = 1500
    cfg.grid.sgrid = 0.004
    cfg.coil.frequency = 1.1e6
    cfg.operation.solve_mode = 1
    cfg.operation.P_max = 120.0
    cfg.operation.I_soll = 80.0
    cfg.sweep.Q0_start = 2.0
    cfg.sweep.Q0_step = 0.5
    cfg.sweep.N = 5
    cfg.rates.rate_model = 0
    cfg.meta.gas = "xenon"
    cfg.meta.preset_id = "holste"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        tmp = f.name
    try:
        cfg.to_params_txt(tmp)

        # Lesen
        cfg2 = RunConfig.from_params_txt(tmp)

        check("R roundtrip", abs(cfg2.geometry.R - 0.05) < 1e-10)
        check("L roundtrip", abs(cfg2.geometry.L - 0.064) < 1e-10)
        check("eta_opt roundtrip", abs(cfg2.grid.eta_opt - 0.60) < 1e-10)
        check("Vgrid roundtrip", abs(cfg2.grid.Vgrid - 1500) < 1e-10)
        check("sgrid roundtrip", abs(cfg2.grid.sgrid - 0.004) < 1e-10)
        check("frequency roundtrip", abs(cfg2.coil.frequency - 1.1e6) < 1e3)
        check("solve_mode roundtrip", cfg2.operation.solve_mode == 1)
        check("P_max roundtrip", abs(cfg2.operation.P_max - 120) < 1e-10)
        check("I_soll roundtrip", abs(cfg2.operation.I_soll - 80) < 1e-10)
        check("Q0_start roundtrip", abs(cfg2.sweep.Q0_start - 2.0) < 1e-10)
        check("Q0_step roundtrip", abs(cfg2.sweep.Q0_step - 0.5) < 1e-10)
        check("N roundtrip", cfg2.sweep.N == 5)
        check("rate_model roundtrip", cfg2.rates.rate_model == 0)
        check("gas roundtrip", cfg2.meta.gas == "xenon")

        # params.txt Format pruefen
        txt = Path(tmp).read_text()
        check("params.txt has R", "R 0.05" in txt)
        check("params.txt has eta_opt", "eta_opt 0.6" in txt)
        check("params.txt has P_RFG_max", "P_RFG_max 120.0" in txt)
        check("params.txt has gas_species", "gas_species xenon" in txt)
    finally:
        Path(tmp).unlink(missing_ok=True)


# ═════════════════════════════════════════════════════════════
# Test 3: Preset-Integration
# ═════════════════════════════════════════════════════════════
def test_preset_integration():
    print("\n--- Test 3: Preset-Integration ---")
    with open(SCRIPT_DIR / "thruster_presets.json", encoding="utf-8") as f:
        data = json.load(f)

    presets = {p["id"]: p for p in data.get("presets", [])}

    # Holste-Preset
    holste = presets["holste"]
    cfg = RunConfig.from_preset(holste, package_gas="iodine",
                                 package_id="py_iodine", cs_database="lafleur_v1")

    check("holste R=0.05", abs(cfg.geometry.R - 0.05) < 1e-10)
    check("holste L=0.064", abs(cfg.geometry.L - 0.064) < 1e-10)
    check("holste eta_opt=0.60", abs(cfg.grid.eta_opt - 0.60) < 1e-10)
    check("holste Vgrid=1500", abs(cfg.grid.Vgrid - 1500) < 1e-10)
    check("holste P_max=120", abs(cfg.operation.P_max - 120) < 1e-10)
    check("holste I_soll=80", abs(cfg.operation.I_soll - 80) < 1e-10)
    check("holste preset_id", cfg.meta.preset_id == "holste")
    check("holste gas=iodine", cfg.meta.gas == "iodine")

    # NPT30-Preset
    npt30 = presets["npt30"]
    cfg2 = RunConfig.from_preset(npt30, package_gas="iodine")
    check("npt30 R=0.015", abs(cfg2.geometry.R - 0.015) < 1e-10)
    check("npt30 eta_opt=0.25", abs(cfg2.grid.eta_opt - 0.25) < 1e-10)

    # Alle Presets muessen valide Config erzeugen
    for pid, preset in presets.items():
        if pid == "custom":
            continue
        c = RunConfig.from_preset(preset)
        errs = c.validate()
        check(f"{pid} preset validates", len(errs) == 0, f"errors: {errs}")


# ═════════════════════════════════════════════════════════════
# Test 4: Validierung
# ═════════════════════════════════════════════════════════════
def test_validation():
    print("\n--- Test 4: Validierung ---")

    # Gueltige Config
    cfg = RunConfig()
    errs = cfg.validate()
    check("default validates", len(errs) == 0, f"errors: {errs}")

    # Ungueltige Geometrie
    cfg2 = RunConfig()
    cfg2.geometry.R = -1
    errs = cfg2.validate()
    check("R<0 detected", any("R" in e for e in errs), f"errors: {errs}")

    # Spule innerhalb Kammer
    cfg3 = RunConfig()
    cfg3.coil.Rc = 0.01  # < R=0.02
    errs = cfg3.validate()
    check("Rc<R detected", any("Rc" in e for e in errs), f"errors: {errs}")

    # I-fix ohne I_soll
    cfg4 = RunConfig()
    cfg4.operation.solve_mode = 1
    cfg4.operation.I_soll = 0
    errs = cfg4.validate()
    check("I-fix without I_soll detected", any("I_soll" in e for e in errs))

    # Ungueltige solve_mode
    cfg5 = RunConfig()
    cfg5.operation.solve_mode = 3
    errs = cfg5.validate()
    check("invalid solve_mode detected", any("solve_mode" in e for e in errs))

    # eta_opt ausserhalb
    cfg6 = RunConfig()
    cfg6.grid.eta_opt = 1.5
    errs = cfg6.validate()
    check("eta_opt>1 detected", any("eta_opt" in e for e in errs))


# ═════════════════════════════════════════════════════════════
# Test 5: Cross-Backend Feldkonsistenz
# ═════════════════════════════════════════════════════════════
def test_cross_backend_fields():
    print("\n--- Test 5: Cross-Backend Feldkonsistenz ---")
    cfg = RunConfig()
    cfg.geometry.R = 0.05
    cfg.grid.eta_opt = 0.60
    cfg.operation.P_max = 120
    cfg.meta.gas = "xenon"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        tmp = f.name
    try:
        cfg.to_params_txt(tmp)
        txt = Path(tmp).read_text()
        lines = {l.split()[0]: l.split()[1] for l in txt.splitlines()
                 if l.strip() and not l.startswith("#") and len(l.split()) >= 2}

        # C++ liest diese Felder aus params.txt:
        cpp_fields = ["R", "L", "betai", "betag", "frequency", "Nw", "R_ohm",
                       "Rc", "lc", "Vgrid", "sgrid", "P_RFG_max", "I_soll",
                       "solve_mode", "rate_model", "density_profile_factor",
                       "gas_species", "cs_database"]
        for f_ in cpp_fields:
            check(f"C++ field '{f_}' in params.txt", f_ in lines, f"keys: {list(lines.keys())}")

        # Python liest diese zusaetzlich:
        py_fields = ["eta_opt"]
        for f_ in py_fields:
            check(f"Python field '{f_}' in params.txt", f_ in lines)

    finally:
        Path(tmp).unlink(missing_ok=True)


# ═════════════════════════════════════════════════════════════
# Test 6: JSON und Flat-Dict
# ═════════════════════════════════════════════════════════════
def test_json_serialization():
    print("\n--- Test 6: JSON/Flat-Dict ---")
    cfg = RunConfig()
    cfg.geometry.R = 0.05
    cfg.grid.eta_opt = 0.60
    cfg.meta.preset_id = "holste"

    j = cfg.to_json()
    check("JSON has geometry", "geometry" in j)
    check("JSON geometry.R", j["geometry"]["R"] == 0.05)
    check("JSON grid.eta_opt", j["grid"]["eta_opt"] == 0.60)
    check("JSON meta.preset_id", j["meta"]["preset_id"] == "holste")

    flat = cfg.to_flat_dict()
    check("flat has R", flat.get("R") == 0.05)
    check("flat has eta_opt", flat.get("eta_opt") == 0.60)
    check("flat has preset_id", flat.get("preset_id") == "holste")

    # JSON roundtrip
    s = json.dumps(j)
    check("JSON serializable", isinstance(s, str) and len(s) > 50)


# ═════════════════════════════════════════════════════════════
# Test 7: Konsistenz mit bestehenden Presets
# ═════════════════════════════════════════════════════════════
def test_preset_params_txt_consistency():
    print("\n--- Test 7: Preset -> params.txt -> from_params_txt ---")
    with open(SCRIPT_DIR / "thruster_presets.json", encoding="utf-8") as f:
        data = json.load(f)

    for preset in data.get("presets", []):
        pid = preset["id"]
        if pid == "custom":
            continue

        cfg = RunConfig.from_preset(preset, package_gas="xenon")

        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            tmp = f.name
        try:
            cfg.to_params_txt(tmp)
            cfg2 = RunConfig.from_params_txt(tmp)

            # Geometrie muss erhalten bleiben
            check(f"{pid}: R preserved",
                  abs(cfg2.geometry.R - cfg.geometry.R) < 1e-10,
                  f"orig={cfg.geometry.R}, read={cfg2.geometry.R}")
            check(f"{pid}: eta_opt preserved",
                  abs(cfg2.grid.eta_opt - cfg.grid.eta_opt) < 1e-10,
                  f"orig={cfg.grid.eta_opt}, read={cfg2.grid.eta_opt}")
            check(f"{pid}: P_max preserved",
                  abs(cfg2.operation.P_max - cfg.operation.P_max) < 1e-10)
        finally:
            Path(tmp).unlink(missing_ok=True)


def test_json_roundtrip():
    print("\n--- Test 8: JSON save/load Roundtrip ---")
    import tempfile
    cfg = RunConfig()
    cfg.geometry.R = 0.05
    cfg.grid.eta_opt = 0.60
    cfg.operation.solve_mode = 1
    cfg.operation.I_soll = 80.0
    cfg.sweep.Q0_start = 2.0
    cfg.sweep.N = 5
    cfg.meta.preset_id = "holste"
    cfg.meta.gas = "iodine"

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        tmp = f.name
    try:
        cfg.save_json(tmp)
        cfg2 = RunConfig.load_json(tmp)
        check("JSON R roundtrip", abs(cfg2.geometry.R - 0.05) < 1e-10)
        check("JSON eta_opt roundtrip", abs(cfg2.grid.eta_opt - 0.60) < 1e-10)
        check("JSON solve_mode roundtrip", cfg2.operation.solve_mode == 1)
        check("JSON I_soll roundtrip", abs(cfg2.operation.I_soll - 80) < 1e-10)
        check("JSON Q0_start roundtrip", abs(cfg2.sweep.Q0_start - 2.0) < 1e-10)
        check("JSON N roundtrip", cfg2.sweep.N == 5)
        check("JSON preset_id roundtrip", cfg2.meta.preset_id == "holste")
        check("JSON gas roundtrip", cfg2.meta.gas == "iodine")
    finally:
        Path(tmp).unlink(missing_ok=True)


def test_python_solver_json_input():
    print("\n--- Test 9: Python-Solver mit run_config.json ---")
    import subprocess, glob
    chem = SCRIPT_DIR / "chemistry" / "xenon_simple" / "chemistry.json"
    if not chem.exists():
        check("chemistry exists", False)
        return

    cfg = RunConfig()
    cfg.geometry.R = 0.02; cfg.geometry.L = 0.04
    cfg.grid.Vgrid = 1500; cfg.grid.sgrid = 0.001; cfg.grid.eta_opt = 0.30
    cfg.coil.frequency = 2.5e6; cfg.coil.Rc = 0.02
    cfg.operation.solve_mode = 2; cfg.operation.P_max = 500
    cfg.sweep.Q0_start = 0.4; cfg.sweep.Q0_step = 0.1; cfg.sweep.N = 2
    cfg.meta.gas = "xenon"
    cfg.save_json(SCRIPT_DIR / "run_config.json")

    cmd = [sys.executable, str(SCRIPT_DIR / "generic_solver.py"), str(chem)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=120,
                       cwd=str(SCRIPT_DIR))

    results = [l for l in r.stdout.strip().split("\n") if l.startswith("RESULT")]
    check("Python JSON: 2 RESULT lines", len(results) == 2, f"got {len(results)}")
    if results:
        parts = results[0].split()
        Te = float(parts[3])
        check("Python JSON: Te > 0", Te > 0, f"Te={Te}")

    # Aufraeumen
    (SCRIPT_DIR / "run_config.json").unlink(missing_ok=True)
    for f in glob.glob(str(SCRIPT_DIR / "python_run_*.jsonl")):
        Path(f).unlink()
    for f in glob.glob(str(SCRIPT_DIR / "simulation_log_*.txt")):
        Path(f).unlink()


def main():
    global passed, failed
    test_defaults()
    test_roundtrip()
    test_preset_integration()
    test_validation()
    test_cross_backend_fields()
    test_json_serialization()
    test_preset_params_txt_consistency()
    test_json_roundtrip()
    test_python_solver_json_input()

    print(f"\n{'='*60}")
    print(f"  Ergebnis: {passed} passed, {failed} failed")
    if errors:
        print(f"  Fehlgeschlagen: {', '.join(errors)}")
    print(f"{'='*60}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
