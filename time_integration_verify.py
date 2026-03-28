#!/usr/bin/env python3
"""
time_integration_verify.py - Zeitintegration der Paper-Gleichungen in Python.

Implementiert die exakt selben 4 ODEs (Eq. 4, 5, 11, 12 aus Chabert 2012)
wie der C++-Solver, aber als zeitabhaengige Integration (RK4).
Vergleicht den stationaeren Endzustand mit dem algebraischen Solver.
"""
import math
import subprocess, sys, os
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
os.chdir(SCRIPT_DIR)

# ─── Physikalische Konstanten (identisch zum C++) ─────────────
me = 9.10938215e-31;  e_ch = 1.602176487e-19;  M = 2.1801711e-25
kB = 1.3806504e-23;   epsilon0 = 8.854187e-12
pi = 3.141592653589793; sigma_i = 1e-18;  conv = 11604.5250061657
Eiz = 1.943408035e-18; Eexc = 1.858524725e-18
kappa = 0.0057;  Tg0 = 293.0; mu_0 = 1.256637061e-6
c_light = 299792458.0
SCCM = 4.477962312e17
Kel_val = 1e-13  # Paper-Konstante
alpha_e = 7.0

# ─── Paper-Geometrie ──────────────────────────────────────────
R = 0.06;  L = 0.10;  betai = 0.7;  betag = 0.3
Rc = 0.07; lc = 0.10; Nw = 5.0;    R_ohm = 2.0
Vgrid = 1000.0; sgrid = 0.001; frequency = 13.56e6
Q0 = 29.0 * SCCM

V = pi*R*R*L;  A = 2*pi*R*R + 2*pi*R*L
Ag = betag*pi*R*R;  Ai = betai*pi*R*R
omega = 2*pi*frequency; k0 = omega/c_light
lambda_0 = R/2.405 + L/pi

# ─── Physik-Funktionen (identisch zum C++) ────────────────────
def Kiz(TeV):
    K1 = 6.73e-15*math.sqrt(TeV)*(3.97+0.643*TeV-0.0368*TeV**2)*math.exp(-12.127/TeV)
    K2 = 6.73e-15*math.sqrt(TeV)*(-0.0001031*TeV**2+6.386*math.exp(-12.127/TeV))
    return 0.5*(K1+K2)

def Kex(TeV): return 1.2921e-13*math.exp(-e_ch*11.6/(kB*TeV*conv))
def Kel(TeV): return Kel_val
def uB(TeV): return math.sqrt(kB*TeV*conv/M)
def vg(Tg): return math.sqrt(8*kB*Tg/(pi*M))
def vi(Ti): return math.sqrt(8*kB*Ti/(pi*M))
def lam_i(ng): return 1/(ng*sigma_i)

def hL(ng):
    lam = lam_i(ng)
    return 0.86*(3+L/(2*lam))**(-0.5)

def hR_f(ng):
    lam = lam_i(ng)
    return 0.80*(4+R/lam)**(-0.5)

def Aeff(ng): return 2*hR_f(ng)*pi*R*L + 2*hL(ng)*pi*R*R
def Aeff1(ng): return 2*hR_f(ng)*pi*R*L + (2-betai)*hL(ng)*pi*R*R

def compute_Pabs(n_val, ng_val, Te_val, P_RFG):
    """Berechne P_abs aus dem RF-Kopplungsmodell.
    Ruft den C++-Solver als Blackbox auf fuer die RF-Berechnung,
    oder verwendet die analytische Naeherung fuer stark ueberdichtes Plasma.

    Bei omega_pe >> omega (typisch fuer dieses Modell):
    Die Skin-Tiefe ist viel kleiner als R, und die Kopplung haengt
    primaer von nu_m/omega ab. Fuer die Verifikation verwenden wir
    eine direkte Naeherung: P_abs = zeta * P_RFG, wobei zeta aus
    einer Referenz-Tabelle interpoliert wird.

    Alternativ: einfaches analytisches Modell fuer R_ind.
    """
    wp2 = n_val * e_ch**2 / (me * epsilon0)
    nu_m = Kel(Te_val) * ng_val
    # Fuer stark ueberdichtes Plasma (wp >> omega):
    # R_ind ~ (2*pi*Nw^2)/(eps0*L*omega) * nu_m*omega/(wp2)  (Naeherung)
    # Dies ist die resistive Komponente der Plasma-Impedanz
    if wp2 > omega**2:
        delta = c_light / math.sqrt(wp2)  # Skin-Tiefe (grob)
        # Vereinfachtes Modell: R_ind proportional zu Kollisionsfrequenz
        R_ind_approx = 2*pi*Nw**2/(epsilon0*L*omega) * nu_m*omega / (wp2 + nu_m**2)
        R_ind_approx = max(R_ind_approx, 1e-4)
    else:
        R_ind_approx = 1e-4  # Sehr schwache Kopplung unter Resonanz

    Ic = math.sqrt(2*P_RFG/(R_ind_approx + R_ohm))
    P_abs = 0.5 * R_ind_approx * Ic**2
    return P_abs


def rhs(n, ng, Ug, Ue, P_RFG):
    """Rechte Seite des ODE-Systems. Zustandsvariablen: n, ng, Ug, Ue.
    Ug = 3/2 ng kB Tg (Gasenergiedichte)
    Ue = 3/2 n kB Te conv (Elektronenenergiedichte, Te in eV)
    """
    if n <= 0 or ng <= 0 or Ug <= 0 or Ue <= 0:
        return 0, 0, 0, 0

    Tg = (2/3)*Ug/(ng*kB)
    Te = (2/3)*Ue/(n*kB*conv)

    if Te <= 0 or Tg <= 0:
        return 0, 0, 0, 0

    P_abs = compute_Pabs(n, ng, Te, P_RFG)
    P_vol = P_abs / V

    # Eq. 4: dn/dt
    dn = n*ng*Kiz(Te) - n*uB(Te)*Aeff(ng)/V

    # Eq. 5: dng/dt
    dng = Q0/V + n*uB(Te)*Aeff1(ng)/V - n*ng*Kiz(Te) - 0.25*ng*vg(Tg)*Ag/V

    # Eq. 11: d(3/2 ng kB Tg)/dt
    Pg1 = 3*me/M * kB*(Te*conv - Tg) * n*ng*Kel(Te)
    Pg2 = 0.25*M*uB(Te)**2 * n*ng*sigma_i*vi(Tg)
    Pg3 = kappa*(Tg - Tg0)/lambda_0 * A/V
    dUg = Pg1 + Pg2 - Pg3

    # Eq. 12/13: d(3/2 n kB Te)/dt
    P2 = Eiz * n*ng*Kiz(Te)
    P3 = Eexc * n*ng*Kex(Te)
    P4 = 3*me/M * kB*(Te*conv - Tg) * n*ng*Kel(Te)
    P5 = alpha_e*kB*Te*conv * n*uB(Te)*Aeff(ng)/V
    dUe = P_vol - (P2 + P3 + P4 + P5)

    return dn, dng, dUg, dUe


def rk4_step(n, ng, Ug, Ue, P_RFG, dt):
    """Ein RK4-Schritt."""
    k1 = rhs(n, ng, Ug, Ue, P_RFG)
    k2 = rhs(n+dt/2*k1[0], ng+dt/2*k1[1], Ug+dt/2*k1[2], Ue+dt/2*k1[3], P_RFG)
    k3 = rhs(n+dt/2*k2[0], ng+dt/2*k2[1], Ug+dt/2*k2[2], Ue+dt/2*k2[3], P_RFG)
    k4 = rhs(n+dt*k3[0], ng+dt*k3[1], Ug+dt*k3[2], Ue+dt*k3[3], P_RFG)

    n_new  = n  + dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
    ng_new = ng + dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
    Ug_new = Ug + dt/6*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
    Ue_new = Ue + dt/6*(k1[3]+2*k2[3]+2*k3[3]+k4[3])
    return n_new, ng_new, Ug_new, Ue_new


def integrate_to_steady(P_RFG, n0=1e17, ng0=None, Te0=3.0, Tg0_init=300.0,
                         dt=1e-7, max_steps=5_000_000, tol=1e-6, check_every=2000):
    """Integriere bis zum stationaeren Zustand."""
    if ng0 is None:
        p0 = 4*kB*Tg0*Q0/(vg(Tg0)*Ag)
        ng0 = p0/(kB*Tg0)

    Ug = 1.5*ng0*kB*Tg0_init
    Ue = 1.5*n0*kB*Te0*conv
    n, ng = n0, ng0

    n_prev, ng_prev, Te_prev, Tg_prev = n*0.5, ng*0.5, Te0*0.5, Tg0_init*0.5

    for step in range(max_steps):
        n, ng, Ug, Ue = rk4_step(n, ng, Ug, Ue, P_RFG, dt)

        # Schutz gegen negative Werte
        n = max(n, 1e10); ng = max(ng, 1e14)
        Ug = max(Ug, 1e-10); Ue = max(Ue, 1e-10)

        if step % check_every == 0 and step > 0:
            Tg = (2/3)*Ug/(ng*kB)
            Te = (2/3)*Ue/(n*kB*conv)

            rn  = abs((n-n_prev)/max(n, 1e-30))
            rng = abs((ng-ng_prev)/max(ng, 1e-30))
            rTe = abs((Te-Te_prev)/max(Te, 1e-30))
            rTg = abs((Tg-Tg_prev)/max(Tg, 1e-30))

            if step % (check_every*50) == 0:
                print(f"    step={step:>8}  Te={Te:.3f}eV  Tg={Tg:.1f}K  "
                      f"n={n:.2e}  ng={ng:.2e}  max_rel={max(rn,rng,rTe,rTg):.2e}")

            if max(rn, rng, rTe, rTg) < tol:
                Tg = (2/3)*Ug/(ng*kB)
                Te = (2/3)*Ue/(n*kB*conv)
                return {"converged": True, "steps": step,
                        "n": n, "ng": ng, "Te": Te, "Tg": Tg,
                        "Ug": Ug, "Ue": Ue}

            n_prev, ng_prev, Te_prev, Tg_prev = n, ng, Te, Tg

    Tg = (2/3)*Ug/(ng*kB)
    Te = (2/3)*Ue/(n*kB*conv)
    return {"converged": False, "steps": max_steps,
            "n": n, "ng": ng, "Te": Te, "Tg": Tg, "Ug": Ug, "Ue": Ue}


def run_cpp_solver(P_RFG):
    """Rufe den C++-Solver im selbstkonsistenten Modus auf."""
    params = {
        "R": R, "L": L, "betai": betai, "betag": betag,
        "frequency": frequency, "Nw": Nw, "R_ohm": R_ohm,
        "Rc": Rc, "lc": lc, "Vgrid": Vgrid, "sgrid": sgrid,
        "P_RFG": P_RFG, "P_RFG_max": 2000.0,
        "Q0sccm": 29.0, "Q0sccm_start": 29.0, "Q0sccm_step": 1.0,
        "jjmax": 1, "I_soll": 15.0, "solve_mode": 2,
        "use_paper_kel": 1, "kel_constant": Kel_val, "alpha_e_wall": alpha_e,
    }
    cfg = SCRIPT_DIR / "_verify_params.txt"
    with open(cfg, "w") as f:
        for k, v in params.items(): f.write(f"{k} {v}\n")
    wp = str(SCRIPT_DIR).replace("\\", "/")
    if len(wp) >= 2 and wp[1] == ":":
        wp = f"/mnt/{wp[0].lower()}{wp[2:]}"
    cmd = (["wsl", "bash", "-c", f'cd "{wp}" && ./chabert _verify_params.txt 2>&1']
           if sys.platform == "win32"
           else [str(SCRIPT_DIR / "chabert"), str(cfg)])
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    cfg.unlink(missing_ok=True)
    for line in r.stdout.splitlines():
        p = line.strip().split()
        if p and p[0] == "RESULT" and len(p) >= 7:
            try:
                return {"n": float(p[1]), "ng": float(p[2]),
                        "Te": float(p[3]), "Tg": float(p[4]), "I_mA": float(p[5])}
            except: pass
    return None


def compute_residuals(n, ng, Te, Tg, P_RFG):
    """Berechne rohe Residuen."""
    P_abs = compute_Pabs(n, ng, Te, P_RFG)
    P_vol = P_abs/V
    r1 = n*ng*Kiz(Te) - n*uB(Te)*Aeff(ng)/V
    r2 = Q0/V + n*uB(Te)*Aeff1(ng)/V - n*ng*Kiz(Te) - 0.25*ng*vg(Tg)*Ag/V
    Pg1 = 3*me/M*kB*(Te*conv-Tg)*n*ng*Kel(Te)
    Pg2 = 0.25*M*uB(Te)**2*n*ng*sigma_i*vi(Tg)
    Pg3 = kappa*(Tg-Tg0)/lambda_0*A/V
    r4 = Pg1+Pg2-Pg3
    P2=Eiz*n*ng*Kiz(Te); P3=Eexc*n*ng*Kex(Te)
    P4=3*me/M*kB*(Te*conv-Tg)*n*ng*Kel(Te)
    P5=alpha_e*kB*Te*conv*n*uB(Te)*Aeff(ng)/V
    r3 = P_vol-(P2+P3+P4+P5)
    return r1, r2, r3, r4


def main():
    print("="*85)
    print("  VERIFIKATION: Zeitintegration vs. Stationaerer Solver")
    print("  (Paper-Geometrie, Q0=29 sccm, Kel=1e-13)")
    print("="*85)

    P_values = [100, 400, 800, 1200]

    print(f"\n{'P[W]':>6} | {'Methode':>12} | {'Te[eV]':>7} {'Tg[K]':>7} "
          f"{'n':>10} {'ng':>10} {'Steps':>8}")
    print("-"*80)

    for P_RFG in P_values:
        print(f"\n  === P_RFG = {P_RFG} W ===")

        # Zeitintegration
        print(f"  Zeitintegration laeuft...")
        res_ode = integrate_to_steady(P_RFG)
        tag = "ODE-RK4"
        status = "conv" if res_ode["converged"] else "MAXITER"
        print(f"{P_RFG:6} | {tag:>12} | {res_ode['Te']:7.3f} {res_ode['Tg']:7.1f} "
              f"{res_ode['n']:10.3e} {res_ode['ng']:10.3e} {res_ode['steps']:8} [{status}]")

        # Residuen am ODE-Endzustand
        r1, r2, r3, r4 = compute_residuals(
            res_ode['n'], res_ode['ng'], res_ode['Te'], res_ode['Tg'], P_RFG)
        prod = res_ode['n']*res_ode['ng']*Kiz(res_ode['Te'])
        loss = res_ode['n']*uB(res_ode['Te'])*Aeff(res_ode['ng'])/V
        print(f"       | {'ODE resid':>12} | r1={r1:.2e} r2={r2:.2e} r3={r3:.2e} r4={r4:.2e}")
        print(f"       | {'r1 balance':>12} | prod={prod:.2e} loss={loss:.2e} ratio={prod/max(loss,1e-30):.4f}")

        # C++ stationaerer Solver
        res_cpp = run_cpp_solver(P_RFG)
        if res_cpp:
            tag = "C++ LM"
            print(f"{P_RFG:6} | {tag:>12} | {res_cpp['Te']:7.3f} {res_cpp['Tg']:7.1f} "
                  f"{res_cpp['n']:10.3e} {res_cpp['ng']:10.3e}")
            r1c, r2c, r3c, r4c = compute_residuals(
                res_cpp['n'], res_cpp['ng'], res_cpp['Te'], res_cpp['Tg'], P_RFG)
            print(f"       | {'C++ resid':>12} | r1={r1c:.2e} r2={r2c:.2e} r3={r3c:.2e} r4={r4c:.2e}")

            # Direkter Vergleich
            dTe = res_ode['Te'] - res_cpp['Te']
            dn = (res_ode['n'] - res_cpp['n'])/max(res_cpp['n'], 1e-30)*100
            dng = (res_ode['ng'] - res_cpp['ng'])/max(res_cpp['ng'], 1e-30)*100
            print(f"       | {'Delta':>12} | dTe={dTe:+.3f}eV  dn={dn:+.1f}%  dng={dng:+.1f}%")

    print(f"\n{'='*85}")
    (SCRIPT_DIR / "_verify_params.txt").unlink(missing_ok=True)


if __name__ == "__main__":
    main()
