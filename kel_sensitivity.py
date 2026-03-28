#!/usr/bin/env python3
"""
kel_sensitivity.py – Kel-Sensitivitaetsstudie bei Paper-Parametern.
"""
from __future__ import annotations
import subprocess, sys, os
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
os.chdir(SCRIPT_DIR)

PAPER_BASE = {
    "R": 0.06, "L": 0.10, "betai": 0.7, "betag": 0.3,
    "frequency": 13.56e6, "Nw": 5.0, "R_ohm": 2.0,
    "Rc": 0.07, "lc": 0.10, "Vgrid": 1000.0, "sgrid": 0.001,
    "P_RFG_max": 2000.0, "Q0sccm": 29.0,
    "I_soll": 15.0, "solve_mode": 2, "use_paper_kel": 1,
}

KEL_VALUES = [0.25e-13, 0.5e-13, 1.0e-13, 1.5e-13, 2.0e-13]
P_VALUES = [100, 400, 800, 1200]


def wsl_path():
    p = str(SCRIPT_DIR).replace("\\", "/")
    if len(p) >= 2 and p[1] == ":":
        return f"/mnt/{p[0].lower()}{p[2:]}"
    return p


def run_single(params):
    cfg_path = SCRIPT_DIR / "_kel_params.txt"
    with open(cfg_path, "w") as f:
        for k, v in params.items():
            f.write(f"{k} {v}\n")
    cmd = ["wsl", "bash", "-c",
           f'cd "{wsl_path()}" && ./chabert _kel_params.txt 2>&1'] \
        if sys.platform == "win32" else [str(SCRIPT_DIR / "chabert"), str(cfg_path)]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return None
    for line in r.stdout.splitlines():
        parts = line.strip().split()
        if parts and parts[0] == "RESULT" and len(parts) >= 7:
            try:
                return {"n": float(parts[1]), "ng": float(parts[2]),
                        "Te": float(parts[3]), "Tg": float(parts[4]),
                        "I_mA": float(parts[5])}
            except (ValueError, IndexError):
                pass
    cfg_path.unlink(missing_ok=True)
    return None


def main():
    print("=" * 85)
    print("  Kel-Sensitivitaetsstudie (Paper-Geometrie, Q0=29 sccm)")
    print("=" * 85)
    Ai = 0.7 * 3.14159 * 0.06**2

    # Sammle alle Ergebnisse
    data = {}  # (kel, P) -> result
    for kel in KEL_VALUES:
        for p in P_VALUES:
            params = {**PAPER_BASE, "P_RFG": p, "kel_constant": kel,
                      "Q0sccm_start": 29.0, "Q0sccm_step": 1.0, "jjmax": 1}
            r = run_single(params)
            data[(kel, p)] = r

    # ─── Tabelle 1: Te(Kel, P) ─────────────────────────────────
    print(f"\n{'Te [eV]':>20}", end="")
    for p in P_VALUES:
        print(f"  P={p:>4}W", end="")
    print(f"  | Paper")
    print("-" * 75)
    paper_Te = {100: 3.0, 400: 3.5, 800: 4.2, 1200: 5.0}
    for kel in KEL_VALUES:
        print(f"  Kel={kel:.2e}", end="")
        for p in P_VALUES:
            r = data.get((kel, p))
            if r:
                print(f"  {r['Te']:7.2f}", end="")
            else:
                print(f"  {'FAIL':>7}", end="")
        print(f"  | {' / '.join(f'{paper_Te[p]:.1f}' for p in P_VALUES)}")

    # ─── Tabelle 2: Ji(Kel, P) ─────────────────────────────────
    print(f"\n{'Ji [A/m2]':>20}", end="")
    for p in P_VALUES:
        print(f"  P={p:>4}W", end="")
    print(f"  | Paper")
    print("-" * 75)
    paper_Ji = {100: 20, 400: 60, 800: 110, 1200: 150}
    for kel in KEL_VALUES:
        print(f"  Kel={kel:.2e}", end="")
        for p in P_VALUES:
            r = data.get((kel, p))
            if r:
                ji = r['I_mA'] / 1000 / Ai
                print(f"  {ji:7.1f}", end="")
            else:
                print(f"  {'FAIL':>7}", end="")
        print(f"  | {' / '.join(f'{paper_Ji[p]:.0f}' for p in P_VALUES)}")

    # ─── Tabelle 3: Tg(Kel, P) ─────────────────────────────────
    print(f"\n{'Tg [K]':>20}", end="")
    for p in P_VALUES:
        print(f"  P={p:>4}W", end="")
    print(f"  | Paper")
    print("-" * 75)
    paper_Tg = {100: 310, 400: 400, 800: 480, 1200: 500}
    for kel in KEL_VALUES:
        print(f"  Kel={kel:.2e}", end="")
        for p in P_VALUES:
            r = data.get((kel, p))
            if r:
                print(f"  {r['Tg']:7.1f}", end="")
            else:
                print(f"  {'FAIL':>7}", end="")
        print(f"  | {' / '.join(f'{paper_Tg[p]:.0f}' for p in P_VALUES)}")

    # ─── Tabelle 4: Ionisierungsgrad ───────────────────────────
    print(f"\n{'iondeg [%]':>20}", end="")
    for p in P_VALUES:
        print(f"  P={p:>4}W", end="")
    print()
    print("-" * 55)
    for kel in KEL_VALUES:
        print(f"  Kel={kel:.2e}", end="")
        for p in P_VALUES:
            r = data.get((kel, p))
            if r:
                iondeg = r['n'] / max(r['ng'], 1) * 100
                print(f"  {iondeg:7.3f}", end="")
            else:
                print(f"  {'FAIL':>7}", end="")
        print()

    # ─── Zusammenfassung ───────────────────────────────────────
    print(f"\n{'='*85}")
    print(f"  ZUSAMMENFASSUNG")
    print(f"{'='*85}")

    (SCRIPT_DIR / "_kel_params.txt").unlink(missing_ok=True)


if __name__ == "__main__":
    main()
