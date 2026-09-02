"""
test_helpers.py -- Gemeinsame Helferschicht fuer Test- und Referenzlaeufe.

Stellt standardisierte Funktionen bereit, die RunConfig als primaere
Konfigurationsquelle nutzen. Alle Test- und Referenzpfade sollen
diese Helfer verwenden statt eigene params.txt-Schreiblogik.

Architektur:
  make_config()      -- Erzeugt RunConfig aus Preset oder Parametern
  write_config()     -- Schreibt run_config.json + params.txt (Compat)
  run_python_solver() -- Startet Python-Solver mit RunConfig
  run_cpp_solver()   -- Startet C++-Solver mit RunConfig
  parse_results()    -- Parst Solver-stdout
  cleanup()          -- Loescht temporaere Dateien
"""
from __future__ import annotations

import sys
import json
import subprocess
import glob
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PYTHON = sys.executable


def make_config(preset_id: str = "custom", gas: str = "xenon",
                cs_database: str = "biagi", **overrides):
    """Erzeuge RunConfig aus Preset oder freien Parametern.

    Args:
        preset_id: Preset-ID aus thruster_presets.json (oder "custom")
        gas: Gas-Spezies
        cs_database: Cross-Section-Datenbank
        **overrides: Ueberschreibe einzelne Felder (z.B. I_soll=80, P_max=120)

    Returns:
        RunConfig-Objekt, konfigurationsbereit.
    """
    from run_config import RunConfig

    if preset_id != "custom":
        # Aus Preset laden
        with open(SCRIPT_DIR / "thruster_presets.json", encoding="utf-8") as f:
            presets = {p["id"]: p for p in json.load(f).get("presets", [])}
        if preset_id not in presets:
            raise ValueError(f"Preset '{preset_id}' nicht gefunden")
        cfg = RunConfig.from_preset(presets[preset_id], package_gas=gas,
                                     cs_database=cs_database)
    else:
        cfg = RunConfig()
        cfg.meta.gas = gas
        cfg.meta.cs_database = cs_database

    cfg.meta.preset_id = preset_id

    # Overrides anwenden (flache Schluessel -> RunConfig-Felder)
    FIELD_MAP = {
        "R": ("geometry", "R"), "L": ("geometry", "L"),
        "betai": ("geometry", "betai"), "betag": ("geometry", "betag"),
        "Vgrid": ("grid", "Vgrid"), "sgrid": ("grid", "sgrid"),
        "eta_opt": ("grid", "eta_opt"),
        "frequency": ("coil", "frequency"), "Nw": ("coil", "Nw"),
        "R_ohm": ("coil", "R_ohm"), "Rc": ("coil", "Rc"), "lc": ("coil", "lc"),
        "solve_mode": ("operation", "solve_mode"),
        "P_max": ("operation", "P_max"), "I_soll": ("operation", "I_soll"),
        "density_profile_factor": ("operation", "density_profile_factor"),
        "alpha_e_wall": ("operation", "alpha_e_wall"),
        "Q0_start": ("sweep", "Q0_start"), "Q0_step": ("sweep", "Q0_step"),
        "N": ("sweep", "N"),
        "rate_model": ("rates", "rate_model"),
    }
    for key, val in overrides.items():
        if key in FIELD_MAP:
            section, attr = FIELD_MAP[key]
            setattr(getattr(cfg, section), attr, val)

    return cfg


def write_config(cfg, work_dir=None):
    """Schreibe RunConfig als run_config.json + params.txt (Compat).

    Args:
        cfg: RunConfig-Objekt
        work_dir: Arbeitsverzeichnis (Default: SCRIPT_DIR)

    Returns:
        Pfad zur run_config.json
    """
    wd = Path(work_dir or SCRIPT_DIR)
    json_path = wd / "run_config.json"
    txt_path = wd / "params.txt"
    cfg.save_json(json_path)
    cfg.to_params_txt(txt_path)
    return str(json_path)


def run_python_solver(chemistry_json: str, cfg=None, work_dir=None):
    """Starte Python-Solver mit RunConfig.

    Args:
        chemistry_json: Pfad zur chemistry.json (relativ oder absolut)
        cfg: RunConfig (optional, wenn schon geschrieben)
        work_dir: Arbeitsverzeichnis

    Returns:
        (stdout, stderr) Tuple
    """
    wd = str(work_dir or SCRIPT_DIR)
    if cfg:
        write_config(cfg, wd)

    chem_path = str(SCRIPT_DIR / chemistry_json) if not Path(chemistry_json).is_absolute() else chemistry_json

    # Primaerer Pfad: nur chemistry.json, alles aus run_config.json
    cmd = [PYTHON, str(SCRIPT_DIR / "generic_solver.py"), chem_path]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300, cwd=wd)
    return r.stdout, r.stderr


def run_python_solver_legacy(chemistry_json: str, mode: int, param: float,
                              Q0_start: float, Q0_step: float, N: int,
                              work_dir=None):
    """Legacy-Pfad: Python-Solver mit CLI-Argumenten.

    MARKIERT: Dieser Pfad ist Kompatibilitaetsmodus.
    """
    wd = str(work_dir or SCRIPT_DIR)
    chem_path = str(SCRIPT_DIR / chemistry_json) if not Path(chemistry_json).is_absolute() else chemistry_json

    cmd = [PYTHON, str(SCRIPT_DIR / "generic_solver.py"),
           chem_path, str(mode), str(param),
           str(Q0_start), str(Q0_step), str(N)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300, cwd=wd)
    return r.stdout, r.stderr


def parse_results(stdout: str) -> dict:
    """Parse Solver-stdout in strukturiertes Dict."""
    out = {"ifix": [], "results": [], "iodine": [], "fails": []}
    for line in stdout.strip().split("\n"):
        p = line.split()
        if not p:
            continue
        tag = p[0]
        if tag == "IFIX_RESULT" and len(p) >= 7:
            out["ifix"].append({
                "q0": float(p[1]), "P": float(p[2]),
                "I_target": float(p[3]), "I_beam": float(p[4]),
                "delta_I": float(p[5]), "status": p[6],
            })
        elif tag == "RESULT" and len(p) >= 7:
            out["results"].append({
                "ne": float(p[1]), "ng": float(p[2]),
                "Te": float(p[3]), "Tg": float(p[4]),
                "I_beam": float(p[5]), "P": float(p[6]),
            })
        elif tag == "IODINE_EXT" and len(p) >= 10:
            out["iodine"].append({
                "diss": float(p[6]), "fIp": float(p[7]),
                "fI2p": float(p[8]), "alpha": float(p[9]),
            })
        elif tag == "SOLVER_FAIL":
            out["fails"].append(line)
    return out


def cleanup(work_dir=None):
    """Loesche temporaere Konfigurationen und Logs."""
    wd = Path(work_dir or SCRIPT_DIR)
    for pattern in ["run_config.json", "params.txt",
                     "python_run_*.jsonl", "simulation_log_*.txt",
                     "xb_params.txt", "test_ifix_params.txt"]:
        for f in glob.glob(str(wd / pattern)):
            Path(f).unlink(missing_ok=True)
