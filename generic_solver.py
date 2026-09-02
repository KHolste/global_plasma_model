"""
generic_solver.py -- Generischer 0D-Plasma-Solver fuer beliebige Chemiepakete.

Verwendet Levenberg-Marquardt mit Multi-Start in Log-Koordinaten.
Arbeitet mit dem dynamischen Zustandsvektor aus plasma_chemistry.BalanceAssembler.
"""
from __future__ import annotations

import math
import sys
import numpy as np
from pathlib import Path

from plasma_chemistry import (
    ChemistryPackage, BalanceAssembler, ThrusterGeometry,
    load_chemistry, KB, CONV, E_CH, ME, PI, EPS0,
)


def solve_steady_state(
    chem: ChemistryPackage,
    geom: ThrusterGeometry,
    P_abs_W: float,
    Q0_pps: float,
    alpha_e_wall: float = 7.0,
    density_profile_factor: float = 1.0,
    max_iter: int = 200,
    tol: float = 1e-3,
    verbose: bool = False,
) -> dict:
    asm = BalanceAssembler(chem, geom)
    N = asm.state_size
    nd = asm._n_heavy
    P_abs_V = P_abs_W / geom.V

    def make_guess(Te0=4.0, ion_f=0.01, diss_f=0.3):
        x0 = asm.default_state(Q0_pps)
        fd = max((x0[asm._idx[s.id]] for s in chem.feedstock_species), default=1e18)
        for sp in asm._heavy:
            if sp.is_neutral and not sp.is_feedstock:
                x0[asm._idx[sp.id]] = diss_f * fd
            elif sp.is_positive_ion:
                x0[asm._idx[sp.id]] = ion_f * fd
            elif sp.is_negative_ion:
                x0[asm._idx[sp.id]] = 0.01 * ion_f * fd
        x0[asm._Te_idx] = Te0
        x0[asm._Tg_idx] = 310.0
        return x0

    starts = [
        make_guess(4.0, 0.01, 0.3),
        make_guess(3.0, 0.005, 0.1),
        make_guess(5.0, 0.02, 0.5),
        make_guess(2.0, 0.001, 0.05),
    ]

    def to_u(xp):
        u = np.empty(N)
        u[:nd] = np.log(np.clip(xp[:nd], 1e6, None))
        u[nd:] = xp[nd:]
        return u

    def from_u(u):
        xp = np.empty(N)
        xp[:nd] = np.exp(np.clip(u[:nd], math.log(1e6), math.log(1e21)))
        xp[nd:] = u[nd:]
        xp[asm._Te_idx] = np.clip(xp[asm._Te_idx], 0.3, 30.0)
        xp[asm._Tg_idx] = np.clip(xp[asm._Tg_idx], 200.0, 3000.0)
        return xp

    def merit(r):
        return np.sqrt(np.mean(r**2))

    # Adaptive Normierung: berechne Skalen aus erstem Startwert
    _r0 = asm.residual(starts[0], P_abs_V, Q0_pps, alpha_e_wall, density_profile_factor)
    _scale = np.empty(N)
    for i in range(nd):
        _scale[i] = max(abs(_r0[i]), starts[0][i], 1e10)
    _scale[asm._Te_idx] = max(abs(_r0[asm._Te_idx]), abs(P_abs_V), 1e3)
    _scale[asm._Tg_idx] = max(abs(_r0[asm._Tg_idx]), 100.0)

    def F(u):
        xp = from_u(u)
        rr = asm.residual(xp, P_abs_V, Q0_pps, alpha_e_wall, density_profile_factor)
        return rr / _scale

    eps_fd = 1e-4

    def jac(u):
        r0 = F(u)
        J = np.zeros((N, N))
        for j in range(N):
            du = max(abs(u[j]) * eps_fd, 1e-6)
            up = u.copy()
            up[j] += du
            J[:, j] = (F(up) - r0) / du
        return J

    # Multi-start LM
    best_u = to_u(starts[0])
    best_m = 1e30
    converged = False
    total_it = 0

    for si, xs in enumerate(starts):
        u = to_u(xs)
        lam = 1e-2
        per_start = max_iter // len(starts)

        for it in range(per_start):
            r = F(u)
            m = merit(r)

            if verbose and it % 50 == 0:
                xc = from_u(u)
                ne = asm.electron_density(xc)
                print(f"  [{si}] it={it:3d} m={m:.3e} lam={lam:.1e} "
                      f"ne={ne:.2e} Te={xc[asm._Te_idx]:.2f} Tg={xc[asm._Tg_idx]:.0f}")

            if m < best_m:
                best_m = m
                best_u = u.copy()

            if m < tol:
                converged = True
                break

            J = jac(u)
            JtJ = J.T @ J
            Jtr = J.T @ r
            D = np.diag(np.maximum(np.diag(JtJ), 1e-12))

            try:
                delta = np.linalg.solve(JtJ + lam * D, -Jtr)
            except np.linalg.LinAlgError:
                lam *= 10
                continue

            # Step limit
            ms = max(np.max(np.abs(delta)), 1e-30)
            if ms > 1.5:
                delta *= 1.5 / ms

            un = u + delta
            mn = merit(F(un))
            if mn < m:
                u = un
                lam = max(lam * 0.3, 1e-8)
            else:
                lam = min(lam * 4.0, 1e8)

        total_it += it + 1
        if converged:
            break

    x = from_u(best_u)
    ne = asm.electron_density(x)
    return {
        "converged": converged,
        "iterations": total_it,
        "merit": float(merit(F(best_u))),
        "state": x,
        "labels": asm.state_labels,
        "ne": ne,
        "Te": x[asm._Te_idx],
        "Tg": x[asm._Tg_idx],
        "densities": {sp.id: float(x[asm._idx[sp.id]]) for sp in asm._heavy},
        "nu_m": asm.collision_frequency(x),
    }


def solve_steady_state_rf(
    chem: ChemistryPackage,
    geom: ThrusterGeometry,
    P_RFG_W: float,
    Q0_pps: float,
    max_rf_iter: int = 15,
    rf_tol: float = 0.02,
    **solver_kwargs,
) -> dict:
    """Solve with RF coupling: P_RFG -> P_abs via zeta (Bessel-BVP).

    Outer fixed-point iteration:
      1. Estimate P_abs from current zeta
      2. Solve plasma state at P_abs
      3. Compute nu_m from converged state
      4. Compute zeta from rf_diagnostics
      5. Update P_abs = zeta * P_RFG
      6. Repeat until P_abs converges

    This reuses rf_diagnostics.compute_rf_diagnostics() -- no duplicate RF physics.

    Returns dict with extra keys: P_RFG, P_abs, zeta, rf_converged.
    """
    try:
        from rf_diagnostics import compute_rf_diagnostics
    except ImportError:
        # Fallback: no RF coupling available, use P directly
        r = solve_steady_state(chem, geom, P_RFG_W, Q0_pps, **solver_kwargs)
        if r:
            r["P_RFG"] = P_RFG_W
            r["P_abs"] = P_RFG_W
            r["zeta"] = 1.0
            r["rf_converged"] = False
        return r

    asm = BalanceAssembler(chem, geom)

    # Initial guess: zeta ~ 0.3 (typical for ICP)
    zeta = 0.3
    P_abs = zeta * P_RFG_W

    for rf_it in range(max_rf_iter):
        r = solve_steady_state(chem, geom, P_abs, Q0_pps, **solver_kwargs)
        if not r or not r.get("converged"):
            break

        # Compute collision frequency from converged state
        nu_m = asm.collision_frequency(r["state"])
        ne = r["ne"]

        if ne <= 0 or nu_m <= 0:
            break

        # RF diagnostics: compute zeta from plasma state
        rf = compute_rf_diagnostics(ne, nu_m, geom)
        zeta_new = rf.zeta

        if zeta_new <= 0:
            zeta_new = 0.01  # Minimum to avoid zero P_abs

        P_abs_new = zeta_new * P_RFG_W

        # Convergence check
        if abs(P_abs_new - P_abs) / max(P_abs, 1e-10) < rf_tol:
            r["P_RFG"] = P_RFG_W
            r["P_abs"] = P_abs_new
            r["zeta"] = zeta_new
            r["rf_converged"] = True
            r["rf_iterations"] = rf_it + 1
            return r

        # Damped update for stability
        P_abs = 0.5 * P_abs + 0.5 * P_abs_new
        zeta = zeta_new

    # Did not converge RF, return best attempt
    if r:
        r["P_RFG"] = P_RFG_W
        r["P_abs"] = P_abs
        r["zeta"] = zeta
        r["rf_converged"] = False
        r["rf_iterations"] = max_rf_iter
    return r


def solve_for_target_current_rf(
    chem, geom, I_soll_mA, Q0_pps,
    P_RFG_min=5.0, P_RFG_max=2000.0, tol_mA=0.5,
    max_bisect=30, **solver_kwargs,
):
    """RF-coupled I-fix: find P_RFG so that I_beam(P_abs(P_RFG)) = I_soll.

    Outer bisection on P_RFG. For each candidate, solve_steady_state_rf()
    derives P_abs via RF coupling (zeta * P_RFG) and solves the plasma.

    Returns: (result_dict, P_RFG_found, I_found, status)
    result_dict has extra keys: P_RFG, P_abs, zeta, rf_converged.
    """
    def I_at_P_RFG(P_RFG):
        r = solve_steady_state_rf(chem, geom, P_RFG, Q0_pps, **solver_kwargs)
        I = _compute_beam_current(r, geom) if r and r.get("converged") else 0.0
        return I, r

    P_max_user = P_RFG_max
    I_lo, r_lo = I_at_P_RFG(P_RFG_min)
    I_hi, r_hi = I_at_P_RFG(P_RFG_max)

    if I_lo > I_soll_mA:
        return r_lo, P_RFG_min, I_lo, "below_P_min"
    if I_hi < I_soll_mA:
        return r_hi, P_max_user, I_hi, "above_P_max"

    p_lo, p_hi = P_RFG_min, P_RFG_max
    best_r, best_P, best_I = r_lo, P_RFG_min, I_lo

    for it in range(max_bisect):
        P_mid = 0.5 * (p_lo + p_hi)
        I_mid, r_mid = I_at_P_RFG(P_mid)

        if abs(I_mid - I_soll_mA) < abs(best_I - I_soll_mA):
            best_r, best_P, best_I = r_mid, P_mid, I_mid

        if abs(I_mid - I_soll_mA) < tol_mA:
            return r_mid, P_mid, I_mid, "converged"

        if I_mid < I_soll_mA:
            p_lo = P_mid
        else:
            p_hi = P_mid

        if p_hi - p_lo < 0.1:
            return best_r, best_P, best_I, "plateau"

    return best_r, best_P, best_I, "max_iter"


def _compute_beam_current(r, geom):
    """Berechne I_beam [mA] mit vollem Extraktionsmodell.

    Nutzt beam_extraction (CL-Limit, Grid-Optik, geom. Aufteilung)
    falls verfuegbar, sonst Fallback auf einfachen Bohm-Fluss.
    """
    if not r or not r["converged"]:
        return 0.0
    d = r["densities"]

    # Bestimme Ionen und ihre Massen
    ion_dens = {}
    ion_mass = {}
    for sp_id in ("I+", "I2+", "Xe+"):
        if sp_id in d and d[sp_id] > 0:
            ion_dens[sp_id] = d[sp_id]
            if sp_id == "I+":
                ion_mass[sp_id] = 2.107e-25
            elif sp_id == "I2+":
                ion_mass[sp_id] = 4.214e-25
            elif sp_id == "Xe+":
                ion_mass[sp_id] = 2.180e-25

    if not ion_dens:
        return 0.0

    # Neutraldichte fuer Extraktionsmodell
    n_neut = sum(v for k, v in d.items() if k not in ("I+", "I2+", "I-", "Xe+"))

    try:
        from beam_extraction import compute_extraction
        ex = compute_extraction(
            ion_densities=ion_dens,
            ion_masses=ion_mass,
            Te_eV=r["Te"],
            n_neutral_total=max(n_neut, 1e10),
            sigma_i=1e-18,
            R=geom.R, L=geom.L,
            betai=geom.betai,
            Vgrid=geom.Vgrid,
            sgrid=geom.sgrid,
            eta_opt=geom.eta_opt,
        )
        return ex.I_beam_mA
    except ImportError:
        # Fallback: einfacher Bohm-Fluss (nur wenn beam_extraction fehlt)
        n_ion = sum(ion_dens.values())
        M_avg = sum(ion_mass[k]*ion_dens[k] for k in ion_dens) / n_ion
        uB = math.sqrt(KB * r["Te"] * CONV / M_avg)
        return E_CH * n_ion * uB * geom.Ai * 1000


def solve_for_target_current(
    chem, geom, I_soll_mA, Q0_pps,
    P_min=1.0, P_max=500.0, tol_mA=0.5,
    max_bisect=30, **solver_kwargs,
):
    """Finde P_abs so dass I_beam(P) = I_soll [mA].

    Bisection ueber P_abs. Robust gegen Nichtkonvergenz.
    """
    def I_at_P(P):
        r = solve_steady_state(chem, geom, P, Q0_pps, **solver_kwargs)
        return _compute_beam_current(r, geom), r

    # Bracket finden: steigend P -> steigende I
    P_max_user = P_max  # Urspruenglicher Nutzerwert merken
    I_lo, r_lo = I_at_P(P_min)
    I_hi, r_hi = I_at_P(P_max)

    if I_lo > I_soll_mA:
        return r_lo, P_min, I_lo, "below_P_min"

    if I_hi < I_soll_mA:
        # Soll bei P_max nicht erreichbar -> P_max_user zurueckgeben (nicht expandiert!)
        return r_hi, P_max_user, I_hi, "above_P_max"

    # Bisection
    p_lo, p_hi = P_min, P_max
    best_r, best_P, best_I = r_lo, P_min, I_lo

    for it in range(max_bisect):
        P_mid = 0.5 * (p_lo + p_hi)
        I_mid, r_mid = I_at_P(P_mid)

        if abs(I_mid - I_soll_mA) < abs(best_I - I_soll_mA):
            best_r, best_P, best_I = r_mid, P_mid, I_mid

        if abs(I_mid - I_soll_mA) < tol_mA:
            return r_mid, P_mid, I_mid, "converged"

        if I_mid < I_soll_mA:
            p_lo = P_mid
        else:
            p_hi = P_mid

        if p_hi - p_lo < 0.01:
            return best_r, best_P, best_I, "plateau"

    return best_r, best_P, best_I, "max_iter"


def _load_geometry_from_config(config_path="params.txt"):
    """Lese Thruster-Geometrie aus params.txt (von GUI geschrieben).

    Faellt auf Defaults zurueck wenn Datei nicht existiert.
    """
    geom = ThrusterGeometry()
    p = Path(config_path)
    if not p.exists():
        return geom

    params = {}
    for line in p.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            params[parts[0]] = parts[1]

    # Geometrie-Felder aus params.txt uebernehmen
    field_map = {
        "R": "R", "L": "L", "betai": "betai", "betag": "betag",
        "Vgrid": "Vgrid", "sgrid": "sgrid", "frequency": "frequency",
        "Nw": "Nw", "R_ohm": "R_ohm", "Rc": "Rc", "lc": "lc",
        "eta_opt": "eta_opt",
    }
    for cfg_key, attr in field_map.items():
        if cfg_key in params:
            try:
                setattr(geom, attr, float(params[cfg_key]))
            except ValueError:
                pass

    return geom


def _load_run_config():
    """Lade RunConfig: primaer aus run_config.json, Fallback auf params.txt."""
    try:
        from run_config import RunConfig
        rc_path = Path("run_config.json")
        if rc_path.exists():
            return RunConfig.load_json(rc_path), "json"
        pt_path = Path("params.txt")
        if pt_path.exists():
            return RunConfig.from_params_txt(pt_path), "params_txt"
    except ImportError:
        pass
    return None, "none"


def _apply_run_config(rc, geom):
    """Uebertrage RunConfig-Felder auf ThrusterGeometry."""
    geom.R = rc.geometry.R
    geom.L = rc.geometry.L
    geom.betai = rc.geometry.betai
    geom.betag = rc.geometry.betag
    geom.Vgrid = rc.grid.Vgrid
    geom.sgrid = rc.grid.sgrid
    geom.eta_opt = rc.grid.eta_opt
    geom.frequency = rc.coil.frequency
    geom.Nw = rc.coil.Nw
    geom.R_ohm = rc.coil.R_ohm
    geom.Rc = rc.coil.Rc
    geom.lc = rc.coil.lc


def main():
    """CLI-Einstiegspunkt. Unterstuetzt drei Aufrufarten:

    1. Primaer (RunConfig-gesteuert):
      python generic_solver.py <chem.json>
      Liest run_config.json fuer alle Parameter.

    2. Legacy-Sweep (CLI-Argumente):
      python generic_solver.py <chem.json> 2 <P_W> <Q0_start> <Q0_step> <N>
      python generic_solver.py <chem.json> 1 <I_mA> <Q0_start> <Q0_step> <N>
      Geometrie aus run_config.json oder params.txt, Sweep aus CLI.

    3. Einzelpunkt (Debug):
      python generic_solver.py <chem.json> [P_W] [Q0_pps]
    """
    SCCM_TO_PPS = 4.477962312e17

    if len(sys.argv) < 2:
        print("Usage:")
        print("  Primary: python generic_solver.py <chem.json>")
        print("           (reads run_config.json for all parameters)")
        print("  Legacy:  python generic_solver.py <chem.json> 2 <P_W> <Q0_start> <Q0_step> <N>")
        print("  I-fix:   python generic_solver.py <chem.json> 1 <I_mA> <Q0_start> <Q0_step> <N>")
        sys.exit(1)

    cp = Path(sys.argv[1])
    chem = load_chemistry(cp)

    # RunConfig laden (primaer JSON, Fallback params.txt)
    rc, rc_source = _load_run_config()

    # Geometrie: aus RunConfig oder Legacy-Fallback
    if rc:
        geom = ThrusterGeometry()
        _apply_run_config(rc, geom)
        P_max_cfg = rc.operation.P_max
    else:
        geom = _load_geometry_from_config()
        P_max_cfg = 500.0
        p_cfg = Path("params.txt")
        if p_cfg.exists():
            for line in p_cfg.read_text(encoding="utf-8").splitlines():
                parts = line.strip().split()
                if len(parts) >= 2 and parts[0] == "P_RFG_max":
                    try: P_max_cfg = float(parts[1])
                    except ValueError: pass

    # Sweep-Parameter: aus RunConfig (primaer) oder CLI (Legacy)
    if len(sys.argv) >= 7:
        # Legacy-CLI-Pfad (Rueckwaertskompatibilitaet)
        mode = int(sys.argv[2])
        param = float(sys.argv[3])
        Q0_start_sccm = float(sys.argv[4])
        Q0_step_sccm = float(sys.argv[5])
        N = int(float(sys.argv[6]))
    elif rc and len(sys.argv) == 2:
        # Primaerer RunConfig-Pfad: alles aus run_config.json
        mode = rc.operation.solve_mode
        param = rc.operation.I_soll if mode == 1 else rc.operation.P_max
        Q0_start_sccm = rc.sweep.Q0_start
        Q0_step_sccm = rc.sweep.Q0_step
        N = rc.sweep.N
    else:
        # Einzelpunkt-Modus (unveraendert)
        N = 0  # Signal fuer Einzelpunkt

    if N > 0:
        sys.stdout.reconfigure(line_buffering=True) if hasattr(sys.stdout, 'reconfigure') else None

        # Run-Logger initialisieren
        import time as _time
        _t0 = _time.time()
        try:
            from run_logger import RunLogger
            _logger = RunLogger(output_dir=Path.cwd())
            _logger.set_metadata(
                backend="python", package=chem.name,
                chemistry_json=str(cp),
                solve_mode="I-fix" if mode == 1 else "SC",
                param=param,
                config_source=rc_source,
                preset=rc.meta.preset_id if rc else "unknown",
            )
            if rc:
                _logger.set_params(rc.to_flat_dict())
            else:
                _logger.set_params({
                    "R": geom.R, "L": geom.L, "betai": geom.betai, "betag": geom.betag,
                    "Vgrid": geom.Vgrid, "sgrid": geom.sgrid, "eta_opt": geom.eta_opt,
                    "frequency": geom.frequency, "P_max": P_max_cfg,
                    "Q0_start": Q0_start_sccm, "Q0_step": Q0_step_sccm, "N": N,
                })
        except ImportError:
            _logger = None
        _count_ok = 0
        _count_fail = 0

        # RF coupling mode (optional, from RunConfig)
        _rf_coupling = False
        if rc and hasattr(rc.operation, 'rf_coupling'):
            _rf_coupling = rc.operation.rf_coupling
        if _rf_coupling:
            print("RF_COUPLING enabled: P_input = P_RFG, P_abs derived via RF model", flush=True)

        for jj in range(N):
            Q0_sccm = Q0_start_sccm + jj * Q0_step_sccm
            Q0_pps = Q0_sccm * SCCM_TO_PPS

            print(f"Q0_STEP {Q0_sccm:.4f} {jj+1} {N}", flush=True)

            if mode == 1:
                # I-fix: finde P so dass I_beam = I_soll
                I_soll = param
                if _rf_coupling:
                    # RF-coupled: search P_RFG, derive P_abs internally
                    r, P_found, I_found, status = solve_for_target_current_rf(
                        chem, geom, I_soll, Q0_pps,
                        P_RFG_min=1.0, P_RFG_max=P_max_cfg, tol_mA=0.5,
                        max_iter=800, tol=0.5)
                else:
                    # Direct: search P_abs
                    r, P_found, I_found, status = solve_for_target_current(
                        chem, geom, I_soll, Q0_pps,
                        P_min=1.0, P_max=P_max_cfg, tol_mA=0.5,
                        max_iter=800, tol=0.5)
                P_out = P_found
                delta_I = I_found - I_soll

                # Strukturierte I-fix-Ergebniszeile
                print(f"IFIX_RESULT {Q0_sccm:.4f} {P_found:.2f} "
                      f"{I_soll:.2f} {I_found:.2f} {delta_I:.2f} {status}",
                      flush=True)

                # RF diagnostics for converged RF-coupled I-fix
                if _rf_coupling and r and r.get("converged"):
                    P_abs_rf = r.get("P_abs", P_found)
                    zeta_rf = r.get("zeta", 1.0)
                    print(f"RF_DIAG P_RFG={P_found:.1f} P_abs={P_abs_rf:.1f} zeta={zeta_rf:.4f}",
                          flush=True)

                if status == "converged":
                    print(f"CONVERGED 0", flush=True)
                elif status in ("below_P_min", "above_P_max", "plateau", "max_iter"):
                    # Menschenlesbare Diagnose
                    if status == "above_P_max":
                        diag = (f"Sollstrom {I_soll:.0f} mA nicht erreichbar: "
                                f"bei P_max={P_found:.0f} W werden nur {I_found:.1f} mA erreicht")
                    elif status == "below_P_min":
                        diag = (f"Sollstrom {I_soll:.0f} mA unterhalb Minimum: "
                                f"bereits bei P_min werden {I_found:.1f} mA erreicht")
                    elif status == "plateau":
                        diag = (f"I_beam={I_found:.1f} mA bei P={P_found:.0f} W, "
                                f"Plateau erreicht (Ziel: {I_soll:.0f} mA)")
                    else:
                        diag = (f"Max. Iterationen bei P={P_found:.0f} W, "
                                f"I_beam={I_found:.1f} mA (Ziel: {I_soll:.0f} mA)")
                    print(f"SOLVER_FAIL {jj} {Q0_sccm:.4f} {status} "
                          f"I_target={I_soll:.1f}mA I_found={I_found:.2f}mA "
                          f"P_max={P_found:.1f}W diag={diag}",
                          flush=True)
                    # Run-Log: I-fix Fehlpunkt
                    if _logger:
                        _count_fail += 1
                        _logger.add_point(idx=jj, Q0_sccm=Q0_sccm, status=status,
                                          P_W=P_found, I_beam_mA=I_found,
                                          I_target_mA=I_soll, delta_I_mA=delta_I,
                                          diag=diag)
            else:
                # SC: feste Leistung
                P_out = param
                if _rf_coupling:
                    # RF-coupled: P_out is P_RFG, derive P_abs
                    r = solve_steady_state_rf(chem, geom, P_out, Q0_pps, max_iter=800, tol=0.5)
                    if r and r.get("converged"):
                        P_abs_rf = r.get("P_abs", P_out)
                        zeta_rf = r.get("zeta", 1.0)
                        print(f"RF_DIAG P_RFG={P_out:.1f} P_abs={P_abs_rf:.1f} zeta={zeta_rf:.4f}", flush=True)
                else:
                    r = solve_steady_state(chem, geom, P_out, Q0_pps, max_iter=800, tol=0.5)

            if r and r.get("converged", False):
                ne = r["ne"]
                d = r["densities"]
                ng = 0
                for sp in chem.feedstock_species:
                    ng = max(ng, d.get(sp.id, 0))
                if ng == 0:
                    ng = max(d.values()) if d else 0

                I_beam = _compute_beam_current(r, geom)

                print(f"RESULT {ne:.4e} {ng:.4e} {r['Te']:.3f} {r['Tg']:.1f} "
                      f"{I_beam:.2f} {P_out:.1f}", flush=True)

                # Iod-erweiterte Ergebnisse
                nI = d.get("I", 0)
                nI2 = d.get("I2", 0)
                nIp = d.get("I+", d.get("Xe+", 0))
                nI2p = d.get("I2+", 0)
                nIm = d.get("I-", 0)
                n_ion = nIp + nI2p
                diss = nI / (nI + 2*nI2) if (nI + 2*nI2) > 0 else 0
                fIp = nIp / n_ion if n_ion > 0 else 0
                fI2p = nI2p / n_ion if n_ion > 0 else 0
                alpha = nIm / ne if ne > 0 else 0
                print(f"IODINE_EXT {nI:.3e} {nI2:.3e} {nIp:.3e} {nI2p:.3e} "
                      f"{nIm:.3e} {diss:.4f} {fIp:.4f} {fI2p:.4f} {alpha:.4f}",
                      flush=True)

                # Run-Log: Erfolgspunkt
                if _logger:
                    _count_ok += 1
                    lp = dict(idx=jj, Q0_sccm=Q0_sccm, status="converged",
                              P_W=P_out, I_beam_mA=I_beam,
                              Te_eV=r["Te"], Tg_K=r["Tg"],
                              ne=ne, ng=ng, converged=True,
                              merit=r.get("merit", 0), iterations=r.get("iterations", 0),
                              diss=diss, fIp=fIp, fI2p=fI2p, alpha=alpha,
                              nI=nI, nI2=nI2, nIp=nIp, nI2p=nI2p, nIm=nIm)
                    if mode == 1:
                        lp.update(I_target_mA=I_soll, delta_I_mA=I_found - I_soll)
                    _logger.add_point(**lp)

            elif mode != 1:
                print(f"SOLVER_FAIL {jj} {Q0_sccm:.4f} not_converged", flush=True)
                if _logger:
                    _count_fail += 1
                    _logger.add_point(idx=jj, Q0_sccm=Q0_sccm, status="not_converged",
                                      diag="SC solver not converged")

        # Run-Log finalisieren
        if _logger:
            elapsed = _time.time() - _t0
            jsonl_path, txt_path = _logger.finalize(
                elapsed=elapsed, count_ok=_count_ok, count_fail=_count_fail)
            print(f"LOG_FILE {txt_path}", flush=True)

    # Einzelpunkt-Modus
    else:
        P = float(sys.argv[2]) if len(sys.argv) > 2 else 20.0
        Q = float(sys.argv[3]) if len(sys.argv) > 3 else 2e17
        print(f"Chemistry: {chem.name}, Species: {list(chem.species.keys())}")
        r = solve_steady_state(chem, geom, P, Q, verbose=True)
        if r["converged"]:
            ne = r["ne"]
            ng = max(r["densities"].values()) if r["densities"] else 0
            print(f"RESULT {ne:.4e} {ng:.4e} {r['Te']:.3f} {r['Tg']:.1f} 0.00 {P:.1f}")
        print(f"\n{'CONVERGED' if r['converged'] else 'NOT CONVERGED'}")
        print(f"Te={r['Te']:.3f} Tg={r['Tg']:.1f} ne={r['ne']:.3e}")
        for k, v in r["densities"].items():
            print(f"  n_{k} = {v:.3e}")


if __name__ == "__main__":
    main()
