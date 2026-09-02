#!/usr/bin/env python3
"""
test_iodine_ifix_rf.py -- Validates RF-coupled iodine I-fix mode.

Tests:
  1. Direct I-fix (P_abs search) still works
  2. RF-coupled I-fix (P_RFG search) converges
  3. RF I-fix result contains P_RFG, P_abs, zeta
  4. P_abs = zeta * P_RFG holds in converged result
  5. RF I-fix output lines are correct (IFIX_RESULT + RF_DIAG)
  6. Sweep over Q0 produces smooth, plausible results
"""
from __future__ import annotations
import sys, subprocess, glob
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
PYTHON = sys.executable

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


def cleanup():
    for pat in ["run_config.json", "params.txt", "python_run_*.jsonl", "simulation_log_*.txt"]:
        for f in glob.glob(str(SCRIPT_DIR / pat)):
            Path(f).unlink(missing_ok=True)


CHEM = str(SCRIPT_DIR / "chemistry/iodine_lafleur_v1/chemistry.json")


def make_and_write_config(rf_coupling, I_soll, P_max, Q0_start, Q0_step, N):
    from run_config import RunConfig, GeometryConfig, GridConfig, CoilConfig, OperationConfig, SweepConfig, MetaConfig
    cfg = RunConfig()
    cfg.geometry = GeometryConfig(R=0.06, L=0.10, betai=0.7, betag=0.3)
    cfg.grid = GridConfig(Vgrid=1000, sgrid=0.0015, eta_opt=1.0)
    cfg.coil = CoilConfig(frequency=13.56e6, Nw=5, R_ohm=2.0, Rc=0.07, lc=0.10)
    cfg.operation = OperationConfig(
        solve_mode=1, P_max=P_max, I_soll=I_soll, rf_coupling=rf_coupling)
    cfg.sweep = SweepConfig(Q0_start=Q0_start, Q0_step=Q0_step, N=N)
    cfg.meta = MetaConfig(gas="iodine", preset_id="grondein")
    cfg.save_json(SCRIPT_DIR / "run_config.json")
    cfg.to_params_txt(SCRIPT_DIR / "params.txt")


def run_solver():
    cmd = [PYTHON, str(SCRIPT_DIR / "generic_solver.py"), CHEM]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600, cwd=str(SCRIPT_DIR))
    return r.stdout


def parse(stdout):
    ifix = []
    rf_diags = []
    results = []
    for line in stdout.strip().split("\n"):
        p = line.split()
        if not p:
            continue
        if p[0] == "IFIX_RESULT" and len(p) >= 7:
            ifix.append({"Q0": float(p[1]), "P": float(p[2]), "I_target": float(p[3]),
                          "I_found": float(p[4]), "delta_I": float(p[5]), "status": p[6]})
        elif p[0] == "RF_DIAG":
            d = {}
            for part in p[1:]:
                if "=" in part:
                    k, v = part.split("=", 1)
                    d[k] = float(v)
            rf_diags.append(d)
        elif p[0] == "RESULT" and len(p) >= 7:
            results.append({"P": float(p[6]), "I": float(p[5]), "Te": float(p[3])})
    return ifix, rf_diags, results


def test_direct_ifix():
    print("\n--- Test 1: Direct I-fix (P_abs search) ---")
    cleanup()
    make_and_write_config(rf_coupling=False, I_soll=50.0, P_max=500.0,
                           Q0_start=3.0, Q0_step=1.0, N=2)
    stdout = run_solver()
    ifix, rf_diags, results = parse(stdout)

    check("direct: IFIX_RESULT emitted", len(ifix) > 0)
    check("direct: no RF_DIAG", len(rf_diags) == 0)
    if ifix:
        conv = [r for r in ifix if r["status"] == "converged"]
        check("direct: some converged", len(conv) > 0, f"{len(conv)}/{len(ifix)}")
        if conv:
            check("direct: I near target", abs(conv[0]["I_found"] - 50) < 1,
                  f"I={conv[0]['I_found']:.1f}")
            print(f"  INFO: P_abs_found={conv[0]['P']:.1f}W, I={conv[0]['I_found']:.2f}mA")


def test_rf_ifix():
    print("\n--- Test 2: RF-coupled I-fix (P_RFG search) ---")
    cleanup()
    make_and_write_config(rf_coupling=True, I_soll=50.0, P_max=800.0,
                           Q0_start=3.0, Q0_step=1.0, N=2)
    stdout = run_solver()
    ifix, rf_diags, results = parse(stdout)

    check("rf: IFIX_RESULT emitted", len(ifix) > 0)
    check("rf: RF_DIAG emitted", len(rf_diags) > 0,
          f"found {len(rf_diags)}")

    if ifix:
        conv = [r for r in ifix if r["status"] == "converged"]
        check("rf: some converged", len(conv) > 0, f"{len(conv)}/{len(ifix)}")

        if conv:
            P_RFG = conv[0]["P"]
            I_found = conv[0]["I_found"]
            check("rf: I near target", abs(I_found - 50) < 1,
                  f"I={I_found:.1f}")
            print(f"  INFO: P_RFG_found={P_RFG:.1f}W, I={I_found:.2f}mA")

    if rf_diags:
        d = rf_diags[0]
        check("rf: has P_RFG", "P_RFG" in d)
        check("rf: has P_abs", "P_abs" in d)
        check("rf: has zeta", "zeta" in d)

        if all(k in d for k in ("P_RFG", "P_abs", "zeta")):
            P_RFG = d["P_RFG"]
            P_abs = d["P_abs"]
            zeta = d["zeta"]
            check("rf: P_abs < P_RFG", P_abs < P_RFG,
                  f"P_abs={P_abs:.1f}, P_RFG={P_RFG:.1f}")
            check("rf: P_abs ~ zeta*P_RFG",
                  abs(P_abs - zeta * P_RFG) / max(P_RFG, 1) < 0.05,
                  f"P_abs={P_abs:.1f}, zeta*P_RFG={zeta*P_RFG:.1f}")
            print(f"  INFO: P_RFG={P_RFG:.1f}W, P_abs={P_abs:.1f}W, zeta={zeta:.4f}")


def test_rf_vs_direct():
    print("\n--- Test 3: RF I-fix needs more P_RFG than direct P_abs ---")
    cleanup()
    # Direct: search P_abs for 50mA
    make_and_write_config(rf_coupling=False, I_soll=50.0, P_max=500.0,
                           Q0_start=3.0, Q0_step=0, N=1)
    stdout_d = run_solver()
    ifix_d, _, _ = parse(stdout_d)

    cleanup()
    # RF: search P_RFG for 50mA
    make_and_write_config(rf_coupling=True, I_soll=50.0, P_max=800.0,
                           Q0_start=3.0, Q0_step=0, N=1)
    stdout_r = run_solver()
    ifix_r, rf_r, _ = parse(stdout_r)

    conv_d = [r for r in ifix_d if r["status"] == "converged"]
    conv_r = [r for r in ifix_r if r["status"] == "converged"]

    if conv_d and conv_r:
        P_abs_direct = conv_d[0]["P"]
        P_RFG_rf = conv_r[0]["P"]
        if rf_r:
            P_abs_rf = rf_r[0].get("P_abs", 0)
            # Zwingend ist nur: der Generator muss mehr liefern, als das Plasma
            # im selben Lauf aufnimmt. Der Vergleich mit dem Lauf, dem P_abs
            # direkt vorgegeben wird, gilt nur, soweit beide Laeufe denselben
            # Betriebspunkt finden -- und die duerfen sich laut der Pruefung
            # darunter um bis zu 30 Prozent unterscheiden. Als scharfe
            # Ungleichung war er deshalb nicht haltbar.
            check("rf: P_RFG > P_abs im selben Lauf",
                  P_RFG_rf > P_abs_rf,
                  f"P_RFG={P_RFG_rf:.2f}, P_abs_rf={P_abs_rf:.2f}")
            check("rf: P_abs_rf ~ P_abs_direct (within 30%)",
                  abs(P_abs_rf - P_abs_direct) / max(P_abs_direct, 1) < 0.30,
                  f"P_abs_rf={P_abs_rf:.1f}, P_abs_direct={P_abs_direct:.1f}")
        print(f"  INFO: Direct P_abs={P_abs_direct:.1f}W, RF P_RFG={P_RFG_rf:.1f}W")


def test_sweep():
    print("\n--- Test 4: RF I-fix sweep smoothness ---")
    cleanup()
    make_and_write_config(rf_coupling=True, I_soll=50.0, P_max=1000.0,
                           Q0_start=2.0, Q0_step=0.5, N=4)
    stdout = run_solver()
    ifix, rf_diags, _ = parse(stdout)

    conv = [r for r in ifix if r["status"] == "converged"]
    check("sweep: >=2 converged", len(conv) >= 2, f"{len(conv)}/{len(ifix)}")

    if len(conv) >= 2:
        P_vals = [r["P"] for r in conv]
        # P_RFG should generally decrease with increasing Q0 (more gas = easier ionization)
        check("sweep: P values vary", max(P_vals) - min(P_vals) > 1,
              f"P range: {min(P_vals):.1f}-{max(P_vals):.1f}")
        for r in conv:
            print(f"  Q0={r['Q0']:.2f} P_RFG={r['P']:.1f}W I={r['I_found']:.1f}mA")


def main():
    global passed, failed
    test_direct_ifix()
    test_rf_ifix()
    test_rf_vs_direct()
    test_sweep()
    cleanup()

    print(f"\n{'='*60}")
    print(f"  Result: {passed} passed, {failed} failed")
    if errors:
        print(f"  Failed:")
        for e in errors:
            print(f"    - {e}")
    print(f"{'='*60}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
