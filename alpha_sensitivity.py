#!/usr/bin/env python3
"""alpha_sensitivity.py - Sensitivitaet des Wandverlustfaktors alpha_e_wall."""
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
    "I_soll": 15.0, "solve_mode": 2,
    "use_paper_kel": 1, "kel_constant": 1.0e-13,
}

ALPHA_VALUES = [4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
P_VALUES = [100, 400, 800, 1200]

# Paper-Referenzwerte (aus Fig. 2, 3, 4)
PAPER_REF = {
    100:  {"Te": 3.0, "Ji": 20,  "Tg": 310},
    400:  {"Te": 3.5, "Ji": 60,  "Tg": 400},
    800:  {"Te": 4.2, "Ji": 110, "Tg": 480},
    1200: {"Te": 5.0, "Ji": 150, "Tg": 500},
}

def wsl_path():
    p = str(SCRIPT_DIR).replace("\\", "/")
    return f"/mnt/{p[0].lower()}{p[2:]}" if len(p) >= 2 and p[1] == ":" else p

def run_single(params):
    cfg = SCRIPT_DIR / "_alpha_params.txt"
    with open(cfg, "w") as f:
        for k, v in params.items(): f.write(f"{k} {v}\n")
    cmd = (["wsl", "bash", "-c", f'cd "{wsl_path()}" && ./chabert _alpha_params.txt 2>&1']
           if sys.platform == "win32"
           else [str(SCRIPT_DIR / "chabert"), str(cfg)])
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return None
    for line in r.stdout.splitlines():
        p = line.strip().split()
        if p and p[0] == "RESULT" and len(p) >= 7:
            try:
                return {"n": float(p[1]), "ng": float(p[2]),
                        "Te": float(p[3]), "Tg": float(p[4]), "I_mA": float(p[5])}
            except (ValueError, IndexError): pass
    return None

def main():
    Ai = 0.7 * 3.14159 * 0.06**2
    print("=" * 90)
    print("  Wandverlustfaktor-Sensitivitaet (Paper-Geometrie, Q0=29 sccm, Kel=1e-13)")
    print("=" * 90)

    data = {}
    for a in ALPHA_VALUES:
        for p in P_VALUES:
            params = {**PAPER_BASE, "P_RFG": p, "alpha_e_wall": a,
                      "Q0sccm_start": 29.0, "Q0sccm_step": 1.0, "jjmax": 1}
            data[(a, p)] = run_single(params)

    # ─── Te ───────────────────────────────────────────────────
    print(f"\n{'Te [eV]':>14}", end="")
    for p in P_VALUES: print(f"  P={p:>4}W", end="")
    print("  | Paper-Ref")
    print("-" * 80)
    for a in ALPHA_VALUES:
        print(f"  alpha={a:.1f}", end="")
        for p in P_VALUES:
            r = data.get((a, p))
            print(f"  {r['Te']:7.2f}" if r else f"  {'FAIL':>7}", end="")
        ref = ' / '.join(str(PAPER_REF[p]['Te']) for p in P_VALUES)
        print(f"  | {ref}")

    # ─── Ji ───────────────────────────────────────────────────
    print(f"\n{'Ji [A/m2]':>14}", end="")
    for p in P_VALUES: print(f"  P={p:>4}W", end="")
    print("  | Paper-Ref")
    print("-" * 80)
    for a in ALPHA_VALUES:
        print(f"  alpha={a:.1f}", end="")
        for p in P_VALUES:
            r = data.get((a, p))
            if r:
                print(f"  {r['I_mA']/1000/Ai:7.1f}", end="")
            else:
                print(f"  {'FAIL':>7}", end="")
        ref = ' / '.join(str(PAPER_REF[p]['Ji']) for p in P_VALUES)
        print(f"  | {ref}")

    # ─── Tg ───────────────────────────────────────────────────
    print(f"\n{'Tg [K]':>14}", end="")
    for p in P_VALUES: print(f"  P={p:>4}W", end="")
    print("  | Paper-Ref")
    print("-" * 80)
    for a in ALPHA_VALUES:
        print(f"  alpha={a:.1f}", end="")
        for p in P_VALUES:
            r = data.get((a, p))
            print(f"  {r['Tg']:7.1f}" if r else f"  {'FAIL':>7}", end="")
        ref = ' / '.join(str(PAPER_REF[p]['Tg']) for p in P_VALUES)
        print(f"  | {ref}")

    # ─── iondeg ───────────────────────────────────────────────
    print(f"\n{'iondeg [%]':>14}", end="")
    for p in P_VALUES: print(f"  P={p:>4}W", end="")
    print()
    print("-" * 55)
    for a in ALPHA_VALUES:
        print(f"  alpha={a:.1f}", end="")
        for p in P_VALUES:
            r = data.get((a, p))
            if r:
                print(f"  {r['n']/max(r['ng'],1)*100:7.2f}", end="")
            else:
                print(f"  {'FAIL':>7}", end="")
        print()

    # ─── Bester Fit ───────────────────────────────────────────
    print(f"\n{'='*90}")
    print(f"  BESTER FIT: Te-Abweichung zum Paper (RMS ueber 4 Leistungspunkte)")
    print(f"{'='*90}")
    for a in ALPHA_VALUES:
        sq = 0; n = 0
        for p in P_VALUES:
            r = data.get((a, p))
            if r:
                sq += (r["Te"] - PAPER_REF[p]['Te'])**2
                n += 1
        if n: print(f"  alpha={a:.1f}  Te-RMS = {(sq/n)**0.5:.3f} eV")

    (SCRIPT_DIR / "_alpha_params.txt").unlink(missing_ok=True)


if __name__ == "__main__":
    main()
