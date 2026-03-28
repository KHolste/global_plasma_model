#!/usr/bin/env python3
"""
r1_diagnosis.py - Ionisierungsbilanz-Diagnose: Zerlegt r1 in Bausteine
und identifiziert die Ursache der Te-Abweichung zum Chabert-Paper.
"""
import math
import subprocess, sys, os
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
os.chdir(SCRIPT_DIR)

# ─── Physikalische Konstanten (identisch zum C++ Code) ────────
me = 9.10938215e-31; e = 1.602176487e-19; M = 2.1801711e-25
kB = 1.3806504e-23; epsilon0 = 8.854187e-12
pi = 3.141592653589793; sigma_i = 1.0e-18; conv = 11604.5250061657

# Paper-Geometrie
R = 0.06; L = 0.10; betai = 0.7; betag = 0.3
V = pi * R * R * L
A = 2*pi*R*R + 2*pi*R*L
SCCM_TO_PPS = 4.477962312e17

def Kiz(Te_eV):
    """Ionisierungsrate [m^3/s]. Te in eV."""
    TeV = Te_eV  # Te ist bereits in eV in unserem Code
    K1 = 6.73e-15 * math.sqrt(TeV) * (3.97 + 0.643*TeV - 0.0368*TeV**2) * math.exp(-12.127/TeV)
    K2 = 6.73e-15 * math.sqrt(TeV) * (-0.0001031*TeV**2 + 6.386*math.exp(-12.127/TeV))
    return 0.5 * (K1 + K2)

def uB(Te_eV):
    """Bohm-Geschwindigkeit [m/s]. Te in eV."""
    return math.sqrt(kB * Te_eV * conv / M)

def lambda_i(ng):
    return 1.0 / (ng * sigma_i)

def hL(ng):
    lam = lambda_i(ng)
    return 0.86 * (3 + L/(2*lam))**(-0.5)

def hR(ng):
    lam = lambda_i(ng)
    return 0.80 * (4 + R/lam)**(-0.5)

def Aeff(ng):
    return 2*hR(ng)*pi*R*L + 2*hL(ng)*pi*R*R

def r1_balance(Te_eV, ng):
    """r1 = 0 bedeutet: ng*Kiz(Te) = uB(Te)*Aeff/V"""
    lhs = ng * Kiz(Te_eV)             # Ionisierungsproduktion
    rhs = uB(Te_eV) * Aeff(ng) / V    # Wandverluste
    return lhs, rhs, lhs - rhs

def solve_Te_for_ng(ng, Te_lo=0.5, Te_hi=15.0):
    """Finde Te bei dem r1=0 fuer gegebenes ng (Bisection)."""
    for _ in range(200):
        Te_mid = 0.5 * (Te_lo + Te_hi)
        _, _, r = r1_balance(Te_mid, ng)
        if r > 0:
            Te_hi = Te_mid
        else:
            Te_lo = Te_mid
        if Te_hi - Te_lo < 1e-6:
            break
    return 0.5 * (Te_lo + Te_hi)


def main():
    print("=" * 85)
    print("  r1-Ionisierungsbilanz-Diagnose (Paper-Geometrie)")
    print("=" * 85)
    print(f"  R={R*100:.0f}cm  L={L*100:.0f}cm  V={V:.6f}m^3  A={A:.6f}m^2")

    # ─── Solver-Ergebnisse laden ──────────────────────────────
    # Fuehre den Solver mit Paper-Parametern aus
    def run_solver(p_rfg):
        params = {
            "R": R, "L": L, "betai": betai, "betag": betag,
            "frequency": 13.56e6, "Nw": 5.0, "R_ohm": 2.0,
            "Rc": 0.07, "lc": 0.10, "Vgrid": 1000.0, "sgrid": 0.001,
            "P_RFG": p_rfg, "P_RFG_max": 2000.0,
            "Q0sccm": 29.0, "Q0sccm_start": 29.0, "Q0sccm_step": 1.0,
            "jjmax": 1, "I_soll": 15.0, "solve_mode": 2,
            "use_paper_kel": 1, "kel_constant": 1e-13,
        }
        cfg = SCRIPT_DIR / "_r1_params.txt"
        with open(cfg, "w") as f:
            for k, v in params.items(): f.write(f"{k} {v}\n")
        wp = str(SCRIPT_DIR).replace("\\", "/")
        if len(wp) >= 2 and wp[1] == ":":
            wp = f"/mnt/{wp[0].lower()}{wp[2:]}"
        cmd = (["wsl", "bash", "-c", f'cd "{wp}" && ./chabert _r1_params.txt 2>&1']
               if sys.platform == "win32"
               else [str(SCRIPT_DIR / "chabert"), str(cfg)])
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        for line in r.stdout.splitlines():
            p = line.strip().split()
            if p and p[0] == "RESULT" and len(p) >= 7:
                try:
                    return {"n": float(p[1]), "ng": float(p[2]),
                            "Te": float(p[3]), "Tg": float(p[4]), "I_mA": float(p[5])}
                except: pass
        return None

    # ─── Tabelle 1: Modell-Ergebnisse + r1-Zerlegung ─────────
    print(f"\n{'P_RFG':>6} {'Te_mod':>7} {'Te_pap':>7} {'ng':>10} "
          f"{'Kiz':>10} {'uB':>8} {'lam_i':>8} {'hL':>6} {'hR':>6} "
          f"{'Aeff':>8} {'LHS':>10} {'RHS':>10} {'r1':>10}")
    print("-" * 120)

    paper_Te = {100: 3.0, 400: 3.5, 800: 4.2, 1200: 5.0}
    solver_data = {}

    for p in [100, 400, 800, 1200]:
        res = run_solver(p)
        if not res:
            print(f"{p:6} FAIL")
            continue
        solver_data[p] = res
        Te = res["Te"]; ng = res["ng"]
        lhs, rhs_val, r1_val = r1_balance(Te, ng)
        lam = lambda_i(ng)
        print(f"{p:6} {Te:7.3f} {paper_Te[p]:7.1f} {ng:10.2e} "
              f"{Kiz(Te):10.2e} {uB(Te):8.1f} {lam*100:8.3f}cm "
              f"{hL(ng):6.3f} {hR(ng):6.3f} "
              f"{Aeff(ng)*1e4:8.2f}cm2 {lhs:10.2e} {rhs_val:10.2e} {r1_val:10.2e}")

    # ─── Tabelle 2: Was waere Te wenn ng aus Paper? ───────────
    print(f"\n{'='*85}")
    print(f"  Gegenrechnung: Te aus r1=0 fuer verschiedene ng-Werte")
    print(f"{'='*85}")
    print(f"{'ng':>12} {'Te(r1=0)':>10} {'Kiz(Te)':>12} {'uB(Te)':>10} {'Aeff':>10}")
    print("-" * 60)

    # Teste verschiedene ng-Werte um zu sehen wie Te reagiert
    for ng_val in [1e19, 2e19, 3e19, 5e19, 6.4e19, 8e19, 1e20]:
        Te_eq = solve_Te_for_ng(ng_val)
        print(f"{ng_val:12.2e} {Te_eq:10.3f} {Kiz(Te_eq):12.3e} "
              f"{uB(Te_eq):10.1f} {Aeff(ng_val)*1e4:10.2f}cm2")

    # ─── Tabelle 3: Sensitivitaet ────────────────────────────
    print(f"\n{'='*85}")
    print(f"  Sensitivitaet: Te bei ng=6.4e19 (Paper-Wert bei 100W) mit Skalierung")
    print(f"{'='*85}")
    ng_ref = 6.4e19
    Te_ref = solve_Te_for_ng(ng_ref)
    print(f"  Referenz: ng={ng_ref:.1e}, Te(r1=0) = {Te_ref:.3f} eV\n")

    # Aeff-Sensitivitaet
    print(f"  Aeff-Skalierung:")
    for scale in [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]:
        # Temporaer Aeff skalieren
        def Aeff_scaled(ng_v):
            return scale * Aeff(ng_v)
        # r1 mit skaliertem Aeff
        def r1_scaled(Te_eV, ng_v):
            return ng_v * Kiz(Te_eV) - uB(Te_eV) * Aeff_scaled(ng_v) / V
        # Bisection
        lo, hi = 0.5, 15.0
        for _ in range(200):
            mid = 0.5*(lo+hi)
            if r1_scaled(mid, ng_ref) > 0: hi = mid
            else: lo = mid
            if hi-lo < 1e-6: break
        Te_s = 0.5*(lo+hi)
        print(f"    Aeff x {scale:.1f}: Te = {Te_s:.3f} eV (delta = {Te_s-Te_ref:+.3f})")

    # Kiz-Sensitivitaet
    print(f"\n  Kiz-Skalierung:")
    for scale in [0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0]:
        def r1_kscaled(Te_eV, ng_v):
            return ng_v * scale * Kiz(Te_eV) - uB(Te_eV) * Aeff(ng_v) / V
        lo, hi = 0.5, 15.0
        for _ in range(200):
            mid = 0.5*(lo+hi)
            if r1_kscaled(mid, ng_ref) > 0: hi = mid
            else: lo = mid
            if hi-lo < 1e-6: break
        Te_s = 0.5*(lo+hi)
        print(f"    Kiz x {scale:.1f}: Te = {Te_s:.3f} eV (delta = {Te_s-Te_ref:+.3f})")

    # ─── Kernfrage: Was ist der Modell-ng vs Paper-ng? ────────
    print(f"\n{'='*85}")
    print(f"  KERNFRAGE: Liegt das Problem bei ng oder bei der Formel?")
    print(f"{'='*85}")
    print(f"\n  Paper Fig.3 bei P=400W: ng ~ 5e19 m^-3, n ~ 5e17 m^-3")
    print(f"  Paper Fig.4 bei P=400W: Te ~ 3.5 eV")
    print(f"  Modell bei P=400W:      ng = {solver_data.get(400,{}).get('ng',0):.2e}, "
          f"Te = {solver_data.get(400,{}).get('Te',0):.3f}")
    print(f"\n  Te aus r1=0 bei ng=5e19: {solve_Te_for_ng(5e19):.3f} eV")
    print(f"  Te aus r1=0 bei ng=3e19: {solve_Te_for_ng(3e19):.3f} eV")
    print(f"  Te aus r1=0 bei ng=1e19: {solve_Te_for_ng(1e19):.3f} eV")
    print(f"\n  => Um Te=3.5 eV zu erreichen, braeuchte man ng ~ ???")

    # Finde ng bei dem Te=3.5
    for target_Te in [3.0, 3.5, 4.0, 5.0]:
        lo_ng, hi_ng = 1e17, 1e21
        for _ in range(200):
            mid_ng = math.sqrt(lo_ng * hi_ng)
            Te_at = solve_Te_for_ng(mid_ng)
            if Te_at > target_Te: hi_ng = mid_ng
            else: lo_ng = mid_ng
            if hi_ng/lo_ng < 1.001: break
        print(f"  Te={target_Te:.1f} eV erfordert ng = {math.sqrt(lo_ng*hi_ng):.2e} m^-3")

    (SCRIPT_DIR / "_r1_params.txt").unlink(missing_ok=True)
    print()


if __name__ == "__main__":
    main()
