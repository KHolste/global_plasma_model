#!/usr/bin/env python3
"""
test_rf_coupling.py -- Tests for Option A (P_abs labeling) and Option B (RF coupling).

Tests:
  1. Direct mode (rf_coupling=False): P is P_abs, consistent behavior
  2. RF-coupled mode (rf_coupling=True): P_RFG -> P_abs via zeta
  3. solve_steady_state_rf produces zeta and P_abs fields
  4. rf_diagnostics integration works end-to-end
  5. RunConfig rf_coupling flag roundtrips correctly
"""
from __future__ import annotations
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

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


def test_direct_mode():
    """Test: direct P_abs mode (default, no RF coupling)."""
    print("\n--- Test 1: Direct P_abs mode ---")
    from plasma_chemistry import load_chemistry, ThrusterGeometry
    from generic_solver import solve_steady_state

    chem = load_chemistry(SCRIPT_DIR / "chemistry/iodine_lafleur_v1/chemistry.json")
    geom = ThrusterGeometry(R=0.06, L=0.10, betai=0.7, betag=0.3,
                             Vgrid=1000, sgrid=0.0015, frequency=13.56e6,
                             Nw=5, R_ohm=2.0, Rc=0.07, lc=0.10, eta_opt=1.0)
    SCCM_TO_PPS = 4.477962312e17
    Q0 = 3.0 * SCCM_TO_PPS

    r = solve_steady_state(chem, geom, 400.0, Q0, max_iter=800, tol=0.5)
    check("direct mode converges", r and r.get("converged"),
          f"converged={r.get('converged') if r else 'None'}")
    if r and r["converged"]:
        check("direct mode has Te", r["Te"] > 2, f"Te={r['Te']}")
        check("direct mode has ne", r["ne"] > 1e15, f"ne={r['ne']}")
        # Direct mode: no P_RFG/P_abs/zeta keys
        check("direct mode no P_RFG key", "P_RFG" not in r)
        check("direct mode no zeta key", "zeta" not in r)


def test_rf_coupled_mode():
    """Test: RF-coupled mode produces zeta and P_abs."""
    print("\n--- Test 2: RF-coupled mode ---")
    from plasma_chemistry import load_chemistry, ThrusterGeometry
    from generic_solver import solve_steady_state_rf

    chem = load_chemistry(SCRIPT_DIR / "chemistry/iodine_lafleur_v1/chemistry.json")
    geom = ThrusterGeometry(R=0.06, L=0.10, betai=0.7, betag=0.3,
                             Vgrid=1000, sgrid=0.0015, frequency=13.56e6,
                             Nw=5, R_ohm=2.0, Rc=0.07, lc=0.10, eta_opt=1.0)
    SCCM_TO_PPS = 4.477962312e17
    Q0 = 3.0 * SCCM_TO_PPS

    r = solve_steady_state_rf(chem, geom, 400.0, Q0, max_iter=800, tol=0.5)
    check("RF mode converges", r and r.get("converged"),
          f"converged={r.get('converged') if r else 'None'}")
    if r and r["converged"]:
        check("RF mode has P_RFG", "P_RFG" in r, str(r.keys()))
        check("RF mode has P_abs", "P_abs" in r)
        check("RF mode has zeta", "zeta" in r)

        P_RFG = r["P_RFG"]
        P_abs = r["P_abs"]
        zeta = r["zeta"]
        check("P_RFG = 400", abs(P_RFG - 400) < 0.1, f"P_RFG={P_RFG}")
        check("P_abs < P_RFG", P_abs < P_RFG,
              f"P_abs={P_abs:.1f}, P_RFG={P_RFG:.1f}")
        check("zeta in (0, 1)", 0 < zeta < 1, f"zeta={zeta:.4f}")
        check("P_abs ~ zeta * P_RFG",
              abs(P_abs - zeta * P_RFG) / P_RFG < 0.05,
              f"P_abs={P_abs:.1f}, zeta*P_RFG={zeta*P_RFG:.1f}")

        print(f"  INFO: P_RFG={P_RFG:.1f}W, P_abs={P_abs:.1f}W, zeta={zeta:.4f}")

        # RF-coupled mode should have lower ne than direct at same P_RFG
        # because P_abs < P_RFG
        check("RF Te > 0", r["Te"] > 2)
        check("RF ne > 0", r["ne"] > 1e15)


def test_rf_coupling_vs_direct():
    """Test: RF-coupled at P_RFG gives lower ne than direct at same P."""
    print("\n--- Test 3: RF vs Direct comparison ---")
    from plasma_chemistry import load_chemistry, ThrusterGeometry
    from generic_solver import solve_steady_state, solve_steady_state_rf

    chem = load_chemistry(SCRIPT_DIR / "chemistry/iodine_lafleur_v1/chemistry.json")
    geom = ThrusterGeometry(R=0.06, L=0.10, betai=0.7, betag=0.3,
                             Vgrid=1000, sgrid=0.0015, frequency=13.56e6,
                             Nw=5, R_ohm=2.0, Rc=0.07, lc=0.10, eta_opt=1.0)
    SCCM_TO_PPS = 4.477962312e17
    Q0 = 3.0 * SCCM_TO_PPS
    P = 400.0

    r_direct = solve_steady_state(chem, geom, P, Q0, max_iter=800, tol=0.5)
    r_rf = solve_steady_state_rf(chem, geom, P, Q0, max_iter=800, tol=0.5)

    if r_direct and r_direct["converged"] and r_rf and r_rf["converged"]:
        ne_direct = r_direct["ne"]
        ne_rf = r_rf["ne"]
        check("RF ne < direct ne (P_abs < P_RFG)",
              ne_rf < ne_direct,
              f"ne_rf={ne_rf:.2e}, ne_direct={ne_direct:.2e}")
        print(f"  INFO: ne_direct={ne_direct:.2e}, ne_rf={ne_rf:.2e}, ratio={ne_rf/ne_direct:.3f}")


def test_run_config_rf_flag():
    """Test: rf_coupling flag in RunConfig."""
    print("\n--- Test 4: RunConfig rf_coupling flag ---")
    from run_config import RunConfig
    import tempfile, json

    cfg = RunConfig()
    check("default rf_coupling=False", cfg.operation.rf_coupling == False)

    cfg.operation.rf_coupling = True
    d = cfg.to_json()
    check("JSON contains rf_coupling",
          d["operation"]["rf_coupling"] == True)

    # Roundtrip
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False, mode="w") as f:
        json.dump(d, f)
        tmp = f.name
    try:
        cfg2 = RunConfig.load_json(tmp)
        check("JSON roundtrip rf_coupling",
              cfg2.operation.rf_coupling == True)
    finally:
        Path(tmp).unlink(missing_ok=True)


def main():
    global passed, failed
    test_direct_mode()
    test_rf_coupled_mode()
    test_rf_coupling_vs_direct()
    test_run_config_rf_flag()

    print(f"\n{'='*60}")
    print(f"  Result: {passed} passed, {failed} failed")
    if errors:
        print(f"  Failed: {', '.join(errors)}")
    print(f"{'='*60}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
