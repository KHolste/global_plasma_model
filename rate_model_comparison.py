#!/usr/bin/env python3
"""
rate_model_comparison.py – Vollvergleich aller Ratenmodell-Varianten.
"""
from __future__ import annotations
import subprocess, sys, os, csv
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
os.chdir(SCRIPT_DIR)

BASE_PARAMS = {
    "R": 0.0215, "L": 0.041, "betai": 0.4, "betag": 0.0529,
    "frequency": 2.53e6, "Nw": 6.0, "R_ohm": 0.48,
    "Rc": 0.025, "lc": 0.030, "Vgrid": 1200.0, "sgrid": 0.00075,
    "P_RFG": 45.0, "P_RFG_max": 200.0,
    "Q0sccm": 0.5, "I_soll": 30.0, "solve_mode": 2,
    "use_paper_kel": 1, "kel_constant": 1e-13,
    "alpha_e_wall": 7.0, "P_abs_scale": 1.0, "Kex_scale": 1.0,
    "density_profile_factor": 1.0,
}

VARIANTS = [
    ("legacy",         {"ionization_model": 0, "excitation_model": 0, "elastic_model": 0}),
    ("exc_tab",        {"ionization_model": 0, "excitation_model": 1, "elastic_model": 0}),
    ("kiz+kex_tab",    {"ionization_model": 1, "excitation_model": 1, "elastic_model": 0}),
    ("full_tabulated", {"ionization_model": 1, "excitation_model": 1, "elastic_model": 1}),
]

Q0_VALUES = [0.40, 0.50, 0.60, 0.70, 0.80]


def wsl_path():
    p = str(SCRIPT_DIR).replace("\\", "/")
    return f"/mnt/{p[0].lower()}{p[2:]}" if len(p) >= 2 and p[1] == ":" else p


def run_single(params):
    cfg = SCRIPT_DIR / "_rmc_params.txt"
    with open(cfg, "w") as f:
        for k, v in params.items():
            f.write(f"{k} {v}\n")
    cmd = (["wsl", "bash", "-c", f'cd "{wsl_path()}" && ./chabert _rmc_params.txt 2>&1']
           if sys.platform == "win32" else ["./chabert", str(cfg)])
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        cfg.unlink(missing_ok=True)
        return None
    cfg.unlink(missing_ok=True)
    for line in r.stdout.splitlines():
        p = line.strip().split()
        if p and p[0] == "RESULT" and len(p) >= 7:
            try:
                return {"n": float(p[1]), "ng": float(p[2]), "Te": float(p[3]),
                        "Tg": float(p[4]), "I_mA": float(p[5]), "P_RFG": float(p[6]),
                        "status": "CONVERGED"}
            except (ValueError, IndexError):
                pass
        if p and p[0] in ("NUMERICAL_FAIL", "NO_PHYSICAL_SOLUTION"):
            return {"status": p[0]}
    return {"status": "NO_RESULT"}


def main():
    print("=" * 100)
    print("  Ratenmodell-Vollvergleich (RIT-4 Parameter, P_RFG=45W)")
    print("=" * 100)

    all_results = []

    for vname, vopts in VARIANTS:
        print(f"\n--- {vname} ---")
        for q0 in Q0_VALUES:
            params = {**BASE_PARAMS, **vopts, "Q0sccm_start": q0, "Q0sccm_step": 0.01, "jjmax": 1}
            res = run_single(params)
            if res and res.get("status") == "CONVERGED":
                res["variant"] = vname
                res["Q0sccm"] = q0
                res["iondeg"] = res["n"] / max(res["ng"], 1) * 100
                all_results.append(res)
                print(f"  Q0={q0:.2f}: Te={res['Te']:.3f} Tg={res['Tg']:.1f} "
                      f"n={res['n']:.2e} ng={res['ng']:.2e} I={res['I_mA']:.2f} "
                      f"iondeg={res['iondeg']:.2f}%")
            else:
                status = res.get("status", "?") if res else "TIMEOUT"
                all_results.append({"variant": vname, "Q0sccm": q0, "status": status,
                                    "Te": 0, "Tg": 0, "n": 0, "ng": 0, "I_mA": 0, "iondeg": 0})
                print(f"  Q0={q0:.2f}: {status}")

    # CSV speichern
    csv_path = SCRIPT_DIR / "rate_model_comparison.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["variant", "Q0sccm", "status", "Te_eV", "Tg_K", "n_m3", "ng_m3",
                     "I_mA", "iondeg_pct", "P_RFG_W"])
        for r in all_results:
            w.writerow([r["variant"], r["Q0sccm"], r.get("status", ""),
                        r.get("Te", ""), r.get("Tg", ""), r.get("n", ""), r.get("ng", ""),
                        r.get("I_mA", ""), r.get("iondeg", ""), r.get("P_RFG", "")])
    print(f"\nCSV: {csv_path}")

    # Vergleichstabelle
    print(f"\n{'='*100}")
    print(f"  VERGLEICHSTABELLE: Te [eV]")
    print(f"{'='*100}")
    print(f"{'Q0':>6}", end="")
    for vn, _ in VARIANTS:
        print(f" | {vn:>16}", end="")
    print()
    print("-" * 80)
    for q0 in Q0_VALUES:
        print(f"{q0:6.2f}", end="")
        for vn, _ in VARIANTS:
            match = [r for r in all_results if r["variant"] == vn and r["Q0sccm"] == q0 and r.get("Te")]
            if match and match[0]["Te"] > 0:
                print(f" | {match[0]['Te']:16.3f}", end="")
            else:
                print(f" | {'FAIL':>16}", end="")
        print()

    # Delta-Tabelle
    print(f"\n{'='*100}")
    print(f"  DELTA Te vs. Legacy [eV]")
    print(f"{'='*100}")
    print(f"{'Q0':>6}", end="")
    for vn, _ in VARIANTS[1:]:
        print(f" | {vn:>16}", end="")
    print()
    print("-" * 65)
    for q0 in Q0_VALUES:
        leg = [r for r in all_results if r["variant"] == "legacy" and r["Q0sccm"] == q0]
        if not leg or leg[0].get("Te", 0) == 0:
            continue
        te_leg = leg[0]["Te"]
        print(f"{q0:6.2f}", end="")
        for vn, _ in VARIANTS[1:]:
            match = [r for r in all_results if r["variant"] == vn and r["Q0sccm"] == q0]
            if match and match[0].get("Te", 0) > 0:
                delta = match[0]["Te"] - te_leg
                print(f" | {delta:+16.3f}", end="")
            else:
                print(f" | {'FAIL':>16}", end="")
        print()

    # ng-Tabelle
    print(f"\n{'='*100}")
    print(f"  ng [m^-3]")
    print(f"{'='*100}")
    print(f"{'Q0':>6}", end="")
    for vn, _ in VARIANTS:
        print(f" | {vn:>16}", end="")
    print()
    print("-" * 80)
    for q0 in Q0_VALUES:
        print(f"{q0:6.2f}", end="")
        for vn, _ in VARIANTS:
            match = [r for r in all_results if r["variant"] == vn and r["Q0sccm"] == q0]
            if match and match[0].get("ng", 0) > 0:
                print(f" | {match[0]['ng']:16.2e}", end="")
            else:
                print(f" | {'FAIL':>16}", end="")
        print()

    print(f"\n{'='*100}")


if __name__ == "__main__":
    main()
