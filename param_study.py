#!/usr/bin/env python3
"""
param_study.py – Systematische Parameterstudie im selbstkonsistenten Modus.

Kartiert den Arbeitsbereich ueber P_RFG × Q0sccm und erfasst
Konvergenz, Strahlstrom, Temperaturen, Ionisierungsgrad.

Ausfuehrung:
    python param_study.py
"""
from __future__ import annotations
import subprocess, sys, os, re
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
os.chdir(SCRIPT_DIR)

# ─── Konfiguration ────────────────────────────────────────────

P_RFG_VALUES = [5, 10, 15, 20, 30, 50, 75, 100]       # [W]
Q0_START     = 0.15
Q0_STEP      = 0.05
Q0_NPOINTS   = 18                                       # 0.15 .. 1.00 sccm

BASE_PARAMS = {
    "R": 0.02, "L": 0.04, "betai": 0.5, "betag": 0.05145,
    "frequency": 2.5e6, "Nw": 6.0, "R_ohm": 0.36,
    "Rc": 0.02, "lc": 0.04, "Vgrid": 1500.0, "sgrid": 0.001,
    "P_RFG_max": 200.0, "Q0sccm": 0.475,
    "I_soll": 15.0, "solve_mode": 2, "use_paper_kel": 0,
}


def wsl_path() -> str:
    p = str(SCRIPT_DIR).replace("\\", "/")
    if len(p) >= 2 and p[1] == ":":
        return f"/mnt/{p[0].lower()}{p[2:]}"
    return p


def run_sweep(p_rfg: float) -> list[dict]:
    """Fuehre einen Q0-Sweep bei festem P_RFG durch."""
    cfg = {**BASE_PARAMS,
           "P_RFG": p_rfg,
           "Q0sccm_start": Q0_START,
           "Q0sccm_step": Q0_STEP,
           "jjmax": Q0_NPOINTS}

    config_path = SCRIPT_DIR / "_study_params.txt"
    with open(config_path, "w") as f:
        for k, v in cfg.items():
            f.write(f"{k} {v}\n")

    if sys.platform == "win32":
        cmd = ["wsl", "bash", "-c",
               f'cd "{wsl_path()}" && ./chabert _study_params.txt 2>&1']
    else:
        cmd = [str(SCRIPT_DIR / "chabert"), str(config_path)]

    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        stdout = r.stdout
    except subprocess.TimeoutExpired:
        return []

    # Parse RESULT und Fail-Zeilen
    results = []
    q0_idx = 0
    for line in stdout.splitlines():
        parts = line.strip().split()
        if not parts:
            continue

        if parts[0] == "Q0_STEP" and len(parts) >= 2:
            current_q0 = float(parts[1])

        elif parts[0] == "RESULT" and len(parts) >= 7:
            try:
                results.append({
                    "P_RFG": p_rfg,
                    "Q0sccm": current_q0,
                    "status": "CONVERGED",
                    "n": float(parts[1]),
                    "ng": float(parts[2]),
                    "Te": float(parts[3]),
                    "Tg": float(parts[4]),
                    "I_mA": float(parts[5]),
                    "iondeg": float(parts[1]) / max(float(parts[2]), 1) * 100,
                })
            except (ValueError, IndexError):
                pass

        elif parts[0] in ("NUMERICAL_FAIL", "NO_PHYSICAL_SOLUTION"):
            q0_val = float(parts[2]) if len(parts) > 2 else current_q0
            results.append({
                "P_RFG": p_rfg,
                "Q0sccm": q0_val,
                "status": parts[0],
                "n": 0, "ng": 0, "Te": 0, "Tg": 0, "I_mA": 0, "iondeg": 0,
            })

    config_path.unlink(missing_ok=True)
    return results


def main():
    print("=" * 78)
    print("  Parameterstudie: Selbstkonsistenter Modus")
    print(f"  P_RFG: {P_RFG_VALUES} W")
    print(f"  Q0: {Q0_START:.2f} .. {Q0_START + (Q0_NPOINTS-1)*Q0_STEP:.2f} sccm"
          f" ({Q0_NPOINTS} Punkte, Schritt {Q0_STEP})")
    print("=" * 78)

    all_results: list[dict] = []

    for p in P_RFG_VALUES:
        print(f"\n  Sweep P_RFG = {p:6.1f} W ...", end=" ", flush=True)
        res = run_sweep(p)
        n_conv = sum(1 for r in res if r["status"] == "CONVERGED")
        n_fail = len(res) - n_conv
        print(f"{n_conv}/{len(res)} konvergiert" +
              (f", {n_fail} Fails" if n_fail else ""))
        all_results.extend(res)

    # ─── Tabellarische Ausgabe ─────────────────────────────────
    print("\n" + "=" * 78)
    print("  ERGEBNISTABELLE")
    print("=" * 78)
    hdr = f"{'P_RFG':>7} {'Q0sccm':>7} {'Status':>22} {'I_mA':>8} {'Te':>7} {'Tg':>7} {'iondeg':>8} {'n':>10}"
    print(hdr)
    print("-" * len(hdr))

    for r in all_results:
        if r["status"] == "CONVERGED":
            print(f"{r['P_RFG']:7.1f} {r['Q0sccm']:7.3f} {'CONVERGED':>22}"
                  f" {r['I_mA']:8.3f} {r['Te']:7.3f} {r['Tg']:7.1f}"
                  f" {r['iondeg']:8.3f} {r['n']:10.2e}")
        else:
            print(f"{r['P_RFG']:7.1f} {r['Q0sccm']:7.3f} {r['status']:>22}"
                  f" {'---':>8} {'---':>7} {'---':>7} {'---':>8} {'---':>10}")

    # ─── Zusammenfassung ───────────────────────────────────────
    print("\n" + "=" * 78)
    print("  KONVERGENZ-KARTE (C=konvergiert, .=Fail)")
    print("=" * 78)

    q0_values = sorted(set(r["Q0sccm"] for r in all_results))

    # Header
    print(f"{'P\\Q0':>7}", end="")
    for q in q0_values:
        print(f" {q:.2f}", end="")
    print()
    print("-" * (8 + 5 * len(q0_values)))

    for p in P_RFG_VALUES:
        print(f"{p:6.1f}W", end="")
        for q in q0_values:
            match = [r for r in all_results
                     if abs(r["P_RFG"] - p) < 0.1 and abs(r["Q0sccm"] - q) < 0.001]
            if match and match[0]["status"] == "CONVERGED":
                print("    C", end="")
            elif match:
                print("    .", end="")
            else:
                print("    ?", end="")
        print()

    # ─── Stabilitaetsbewertung ─────────────────────────────────
    conv = [r for r in all_results if r["status"] == "CONVERGED"]
    fails = [r for r in all_results if r["status"] != "CONVERGED"]

    print(f"\n{'='*78}")
    print(f"  ZUSAMMENFASSUNG")
    print(f"{'='*78}")
    print(f"  Total Punkte:    {len(all_results)}")
    print(f"  Konvergiert:     {len(conv)}")
    print(f"  Fehlgeschlagen:  {len(fails)}")

    if conv:
        print(f"\n  Konvergierte Bereiche:")
        print(f"    P_RFG:  {min(r['P_RFG'] for r in conv):.0f} .. {max(r['P_RFG'] for r in conv):.0f} W")
        print(f"    Q0:     {min(r['Q0sccm'] for r in conv):.2f} .. {max(r['Q0sccm'] for r in conv):.2f} sccm")
        print(f"    I_mA:   {min(r['I_mA'] for r in conv):.2f} .. {max(r['I_mA'] for r in conv):.2f} mA")
        print(f"    Te:     {min(r['Te'] for r in conv):.2f} .. {max(r['Te'] for r in conv):.2f} eV")
        print(f"    Tg:     {min(r['Tg'] for r in conv):.1f} .. {max(r['Tg'] for r in conv):.1f} K")
        print(f"    iondeg: {min(r['iondeg'] for r in conv):.3f} .. {max(r['iondeg'] for r in conv):.3f} %")

    if fails:
        print(f"\n  Fehlgeschlagene Bereiche:")
        for r in fails:
            print(f"    P={r['P_RFG']:.0f}W Q0={r['Q0sccm']:.2f}sccm: {r['status']}")

    # ─── CSV speichern ─────────────────────────────────────────
    csv_path = SCRIPT_DIR / "param_study_results.csv"
    with open(csv_path, "w") as f:
        f.write("P_RFG_W,Q0sccm,status,I_mA,Te_eV,Tg_K,iondeg_pct,n_m3\n")
        for r in all_results:
            f.write(f"{r['P_RFG']},{r['Q0sccm']},{r['status']},"
                    f"{r['I_mA']},{r['Te']},{r['Tg']},{r['iondeg']},{r['n']}\n")
    print(f"\n  CSV gespeichert: {csv_path.name}")
    print("=" * 78)


if __name__ == "__main__":
    main()
