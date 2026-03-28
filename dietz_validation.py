#!/usr/bin/env python3
"""
dietz_validation.py – Validierungslauf gegen Dietz et al. (2021) RIT-4 Benchmark.

Fuehrt den Solver mit RIT-4-Parametern aus und erzeugt einen CSV-Datensatz
zur spaeteren Gegenuberstellung mit Dietz Fig. 9/10/12.

Ausfuehrung:
    python dietz_validation.py

Voraussetzung:
    Kompiliertes Binary 'chabert' (Linux/WSL) oder 'chabert.exe'
    config_dietz_rit4.txt im selben Verzeichnis
"""
from __future__ import annotations

import os
import sys
import subprocess
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
os.chdir(SCRIPT_DIR)

CONFIG_FILE = "config_dietz_rit4.txt"
CSV_OUTPUT  = "dietz_validation_results.csv"

Q0_VALUES = [0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80]


# ─── Hilfsfunktionen (aus test_solver.py) ─────────────────────

def find_binary() -> str:
    if sys.platform == "win32":
        wsl = subprocess.run(["wsl", "bash", "-c", f'test -f "{wsl_path()}/chabert"'],
                             capture_output=True)
        if wsl.returncode == 0:
            return "wsl"
        exe = SCRIPT_DIR / "chabert.exe"
        if exe.exists():
            return str(exe)
    else:
        native = SCRIPT_DIR / "chabert"
        if native.exists():
            return str(native)
    return ""


def wsl_path() -> str:
    p = str(SCRIPT_DIR).replace("\\", "/")
    if len(p) >= 2 and p[1] == ":":
        return f"/mnt/{p[0].lower()}{p[2:]}"
    return p


def load_base_config() -> dict[str, str]:
    """Lade config_dietz_rit4.txt als Key-Value-Dict."""
    cfg = {}
    config_path = SCRIPT_DIR / CONFIG_FILE
    if not config_path.exists():
        print(f"FEHLER: {CONFIG_FILE} nicht gefunden!")
        sys.exit(1)
    with open(config_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            if len(parts) >= 2:
                key = parts[0]
                val = parts[1].split("#")[0].strip()  # Inline-Kommentare entfernen
                cfg[key] = val
    return cfg


def run_single(base_cfg: dict, q0sccm: float) -> dict:
    """Fuehre einen einzelnen Solver-Lauf fuer einen Q0-Wert durch."""
    # Temporaere Config mit dem spezifischen Q0
    params = dict(base_cfg)
    params["Q0sccm"] = str(q0sccm)
    params["Q0sccm_start"] = str(q0sccm)
    params["Q0sccm_step"] = "0.01"
    params["jjmax"] = "1"

    config_path = SCRIPT_DIR / "_dietz_tmp.txt"
    with open(config_path, "w") as f:
        f.write("# Dietz-Validierung (automatisch generiert)\n")
        for k, v in params.items():
            f.write(f"{k} {v}\n")

    binary = find_binary()
    if not binary:
        return {"status": "NO_BINARY"}

    if binary == "wsl":
        cmd = ["wsl", "bash", "-c",
               f'cd "{wsl_path()}" && ./chabert _dietz_tmp.txt 2>&1']
    else:
        cmd = [binary, str(config_path)]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        config_path.unlink(missing_ok=True)
        return {"status": "TIMEOUT"}

    config_path.unlink(missing_ok=True)

    # Parse stdout
    for line in result.stdout.splitlines():
        parts = line.strip().split()
        if not parts:
            continue

        if parts[0] == "RESULT" and len(parts) >= 7:
            try:
                return {
                    "status": "CONVERGED",
                    "Q0sccm": q0sccm,
                    "n":      float(parts[1]),
                    "ng":     float(parts[2]),
                    "Te":     float(parts[3]),
                    "Tg":     float(parts[4]),
                    "I_mA":   float(parts[5]),
                    "P_RFG":  float(parts[6]),
                }
            except (ValueError, IndexError):
                pass

        if parts[0] in ("NUMERICAL_FAIL", "NO_PHYSICAL_SOLUTION"):
            return {
                "status": parts[0],
                "Q0sccm": q0sccm,
                "reason": " ".join(parts[1:]),
            }

    return {"status": "NO_RESULT", "Q0sccm": q0sccm}


# ─── Hauptprogramm ───────────────────────────────────────────

def main() -> int:
    print("=" * 70)
    print("  Dietz et al. (2021) RIT-4 Validierungslauf")
    print("=" * 70)

    binary = find_binary()
    if not binary:
        print(f"\nFEHLER: Solver-Binary nicht gefunden.")
        print("Bitte zuerst kompilieren.")
        return 1
    print(f"Binary: {binary}")

    base_cfg = load_base_config()
    print(f"Config: {CONFIG_FILE}")
    print(f"Q0-Sweep: {Q0_VALUES} sccm")
    print(f"P_RFG: {base_cfg.get('P_RFG', '?')} W")
    print(f"solve_mode: {base_cfg.get('solve_mode', '?')}")
    print(f"density_profile_factor: {base_cfg.get('density_profile_factor', '1.0')}")

    # Sweep durchfuehren
    results = []
    n_conv = 0
    n_fail = 0

    print(f"\n{'Q0[sccm]':>10} {'Status':>22} {'I_mA':>8} {'Te[eV]':>8} "
          f"{'Tg[K]':>8} {'n':>10} {'ng':>10} {'P_RFG':>8}")
    print("-" * 90)

    for q0 in Q0_VALUES:
        res = run_single(base_cfg, q0)
        results.append(res)

        if res["status"] == "CONVERGED":
            n_conv += 1
            print(f"{q0:10.3f} {'CONVERGED':>22} {res['I_mA']:8.3f} "
                  f"{res['Te']:8.3f} {res['Tg']:8.1f} {res['n']:10.2e} "
                  f"{res['ng']:10.2e} {res['P_RFG']:8.2f}")
        else:
            n_fail += 1
            reason = res.get("reason", res["status"])[:40]
            print(f"{q0:10.3f} {res['status']:>22} {'---':>8} "
                  f"{'---':>8} {'---':>8} {'---':>10} {'---':>10} {'---':>8}  {reason}")

    # CSV schreiben
    csv_path = SCRIPT_DIR / CSV_OUTPUT
    with open(csv_path, "w") as f:
        f.write("Q0sccm,status,I_mA,Te_eV,Tg_K,n_m3,ng_m3,P_RFG_W\n")
        for res in results:
            q0 = res.get("Q0sccm", 0)
            if res["status"] == "CONVERGED":
                f.write(f"{q0},{res['status']},{res['I_mA']},{res['Te']},"
                        f"{res['Tg']},{res['n']},{res['ng']},{res['P_RFG']}\n")
            else:
                f.write(f"{q0},{res['status']},,,,,,\n")

    # Zusammenfassung
    print(f"\n{'=' * 70}")
    print(f"  Ergebnis: {n_conv}/{len(Q0_VALUES)} konvergiert, {n_fail} fehlgeschlagen")
    print(f"  CSV: {CSV_OUTPUT}")
    print(f"{'=' * 70}")

    if n_conv > 0:
        conv = [r for r in results if r["status"] == "CONVERGED"]
        print(f"\n  Bereiche der konvergierten Punkte:")
        print(f"    Q0:   {min(r['Q0sccm'] for r in conv):.2f} - "
              f"{max(r['Q0sccm'] for r in conv):.2f} sccm")
        print(f"    I_mA: {min(r['I_mA'] for r in conv):.2f} - "
              f"{max(r['I_mA'] for r in conv):.2f} mA")
        print(f"    Te:   {min(r['Te'] for r in conv):.2f} - "
              f"{max(r['Te'] for r in conv):.2f} eV")
        print(f"    Tg:   {min(r['Tg'] for r in conv):.1f} - "
              f"{max(r['Tg'] for r in conv):.1f} K")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
