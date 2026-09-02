#!/usr/bin/env python3
"""
test_reference_workflows.py -- Datengetriebener Referenz-Workflow-Runner.

Laedt reference_workflows.json, fuehrt jeden Fall aus und prueft die Ergebnisse
gegen die definierten Erwartungen. Toleranzbasiert, nicht bitgenau.

Ausfuehrung:
    python test_reference_workflows.py              # alle Workflows
    python test_reference_workflows.py holste_xe    # nur passende IDs

Ausgabe: PASS/FAIL pro Pruefpunkt, Zusammenfassung am Ende.
"""
from __future__ import annotations
import sys, json, subprocess
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
        print(f"    PASS: {name}")
    else:
        failed += 1
        msg = f"    FAIL: {name}"
        if detail:
            msg += f" -- {detail}"
        print(msg)
        errors.append(name)


# ═════════════════════════════════════════════════════════════
# Konfiguration via RunConfig (primaerer Pfad)
# ═════════════════════════════════════════════════════════════
from test_helpers import make_config, write_config, run_python_solver, cleanup as _cleanup


def load_presets():
    with open(SCRIPT_DIR / "thruster_presets.json", encoding="utf-8") as f:
        data = json.load(f)
    return {p["id"]: p for p in data.get("presets", [])}


def run_solver(chem_json, mode, param, Q0_start, Q0_step, N):
    """Legacy-Wrapper: Wird noch von check_expectations aufgerufen."""
    from test_helpers import run_python_solver_legacy
    stdout, stderr = run_python_solver_legacy(
        chem_json, 1 if mode == "ifix" else 2, param, Q0_start, Q0_step, N)
    return stdout, stderr


# ═════════════════════════════════════════════════════════════
# Output parsen
# ═════════════════════════════════════════════════════════════
def parse_output(stdout):
    ifix_results = []
    results = []
    iodine_exts = []
    solver_fails = []

    for line in stdout.strip().split("\n"):
        parts = line.split()
        if not parts:
            continue
        tag = parts[0]

        if tag == "IFIX_RESULT" and len(parts) >= 7:
            ifix_results.append({
                "q0": float(parts[1]), "P": float(parts[2]),
                "I_target": float(parts[3]), "I_beam": float(parts[4]),
                "delta_I": float(parts[5]), "status": parts[6],
            })
        elif tag == "RESULT" and len(parts) >= 7:
            results.append({
                "ne": float(parts[1]), "ng": float(parts[2]),
                "Te": float(parts[3]), "Tg": float(parts[4]),
                "I_beam": float(parts[5]), "P": float(parts[6]),
            })
        elif tag == "IODINE_EXT" and len(parts) >= 10:
            iodine_exts.append({
                "nI": float(parts[1]), "nI2": float(parts[2]),
                "nIp": float(parts[3]), "nI2p": float(parts[4]),
                "nIm": float(parts[5]), "diss": float(parts[6]),
                "fIp": float(parts[7]), "fI2p": float(parts[8]),
                "alpha": float(parts[9]),
            })
        elif tag == "SOLVER_FAIL":
            solver_fails.append(line)

    return {
        "ifix": ifix_results,
        "results": results,
        "iodine": iodine_exts,
        "fails": solver_fails,
    }


# ═════════════════════════════════════════════════════════════
# Erwartungen pruefen
# ═════════════════════════════════════════════════════════════
def in_range(val, rng):
    return rng[0] <= val <= rng[1]


def check_expectations(wf_id, exp, parsed):
    ifix = parsed["ifix"]
    results = parsed["results"]
    iodine = parsed["iodine"]
    fails = parsed["fails"]

    # Grundlegend: irgendwelche Ergebnisse?
    has_ifix = len(ifix) > 0
    has_results = len(results) > 0

    if "all_converged" in exp:
        if exp["all_converged"]:
            check(f"{wf_id}: output exists", has_results, f"results={len(results)}")
            if has_ifix:
                conv = [r for r in ifix if r["status"] == "converged"]
                check(f"{wf_id}: all converged",
                      len(conv) == len(ifix),
                      f"converged={len(conv)}/{len(ifix)}")
        else:
            if has_ifix:
                conv = [r for r in ifix if r["status"] == "converged"]
                check(f"{wf_id}: not all converged",
                      len(conv) < len(ifix),
                      f"converged={len(conv)}/{len(ifix)}")

    if "some_converged" in exp and exp["some_converged"]:
        if has_ifix:
            conv = [r for r in ifix if r["status"] == "converged"]
            check(f"{wf_id}: some converged",
                  len(conv) > 0,
                  f"converged={len(conv)}/{len(ifix)}")

    # Status-Muster
    if "status_pattern" in exp and has_ifix:
        pattern = exp["status_pattern"]
        matching = [r for r in ifix if r["status"] == pattern]
        check(f"{wf_id}: status={pattern}",
              len(matching) > 0,
              f"found {len(matching)} with status={pattern}")

    # I_beam Bereich
    if "I_beam_range" in exp:
        rng = exp["I_beam_range"]
        # Prüfe gegen alle verfügbaren I_beam-Werte
        all_I = []
        if has_ifix:
            all_I = [r["I_beam"] for r in ifix]
        elif has_results:
            all_I = [r["I_beam"] for r in results]
        if all_I:
            for I_val in all_I:
                if not in_range(I_val, rng):
                    check(f"{wf_id}: I_beam in [{rng[0]},{rng[1]}]",
                          False, f"I_beam={I_val:.2f}")
                    break
            else:
                check(f"{wf_id}: I_beam in [{rng[0]},{rng[1]}]", True)

    # P Bereich (fuer konvergierte Punkte)
    if "P_range" in exp and has_ifix:
        rng = exp["P_range"]
        conv_P = [r["P"] for r in ifix if r["status"] == "converged"]
        if conv_P:
            ok = all(in_range(p, rng) for p in conv_P)
            check(f"{wf_id}: P in [{rng[0]},{rng[1]}]",
                  ok, f"P values: {[f'{p:.1f}' for p in conv_P]}")

    # P am Limit
    if "P_at_limit" in exp and has_ifix:
        limit = exp["P_at_limit"]
        fail_P = [r["P"] for r in ifix if r["status"] != "converged"]
        if fail_P:
            check(f"{wf_id}: fail at P_max={limit}",
                  all(abs(p - limit) < 1.0 for p in fail_P),
                  f"fail P: {[f'{p:.1f}' for p in fail_P]}")

    # Te Bereich
    if "Te_range" in exp and has_results:
        rng = exp["Te_range"]
        ok = all(in_range(r["Te"], rng) for r in results)
        check(f"{wf_id}: Te in [{rng[0]},{rng[1]}]",
              ok, f"Te: {[f'{r['Te']:.2f}' for r in results]}")

    # delta_I
    if "delta_I_below" in exp and has_ifix:
        threshold = exp["delta_I_below"]
        ok = all(r["delta_I"] < threshold for r in ifix)
        check(f"{wf_id}: delta_I < {threshold}",
              ok, f"delta_I: {[f'{r['delta_I']:.1f}' for r in ifix]}")

    # Iod-Erweiterung
    if exp.get("has_iodine_ext") and has_results:
        check(f"{wf_id}: IODINE_EXT present",
              len(iodine) > 0,
              f"IODINE_EXT={len(iodine)}, RESULT={len(results)}")

    # Iod-spezifische Bereiche
    if "diss_range" in exp and iodine:
        rng = exp["diss_range"]
        ok = all(in_range(ie["diss"], rng) for ie in iodine)
        check(f"{wf_id}: diss in [{rng[0]},{rng[1]}]",
              ok, f"diss: {[f'{ie['diss']:.3f}' for ie in iodine]}")

    if "fIp_range" in exp and iodine:
        rng = exp["fIp_range"]
        ok = all(in_range(ie["fIp"], rng) for ie in iodine)
        check(f"{wf_id}: fIp in [{rng[0]},{rng[1]}]",
              ok, f"fIp: {[f'{ie['fIp']:.3f}' for ie in iodine]}")

    # Fail-Meldungen
    if "fail_reason_contains" in exp and fails:
        phrase = exp["fail_reason_contains"]
        found = any(phrase in f for f in fails)
        check(f"{wf_id}: fail contains '{phrase}'",
              found, f"fails: {len(fails)} lines")


# ═════════════════════════════════════════════════════════════
# Hauptprogramm
# ═════════════════════════════════════════════════════════════
def main():
    global passed, failed

    # Filter?
    filter_ids = sys.argv[1:] if len(sys.argv) > 1 else []

    # Lade Workflows und Presets
    with open(SCRIPT_DIR / "reference_workflows.json", encoding="utf-8") as f:
        wf_data = json.load(f)
    workflows = wf_data["workflows"]
    presets = load_presets()

    if filter_ids:
        workflows = [w for w in workflows
                     if any(fid in w["id"] for fid in filter_ids)]
        if not workflows:
            print(f"Keine Workflows gefunden fuer Filter: {filter_ids}")
            return 1

    print(f"Referenz-Workflows: {len(workflows)} Faelle")
    print("=" * 70)

    for wf in workflows:
        wf_id = wf["id"]
        print(f"\n--- {wf['name']} ---")
        print(f"    {wf['description']}")

        # RunConfig aus Preset erzeugen (primaerer Pfad)
        preset = presets.get(wf["preset"])
        if not preset:
            check(f"{wf_id}: preset exists", False, f"preset={wf['preset']}")
            continue

        mode = wf["mode"]
        p = wf["params"]
        gas = "iodine" if "iodine" in wf["chemistry"] else "xenon"

        # SC: P_max = gewuenschte Leistung (Solver nutzt P_max als P_abs im SC-Modus)
        # I-fix: P_max = obere Leistungsgrenze fuer Bisection
        if mode == "ifix":
            cfg = make_config(
                preset_id=wf["preset"], gas=gas, solve_mode=1,
                P_max=p.get("P_max", 500.0), I_soll=p["I_soll"],
                Q0_start=p["Q0_start"], Q0_step=p["Q0_step"], N=p["N"])
        else:
            cfg = make_config(
                preset_id=wf["preset"], gas=gas, solve_mode=2,
                P_max=p["P"],  # SC: P ist die Zielleistung
                Q0_start=p["Q0_start"], Q0_step=p["Q0_step"], N=p["N"])
        write_config(cfg)

        # Chemistry pruefen
        chem_path = SCRIPT_DIR / wf["chemistry"]
        check(f"{wf_id}: chemistry exists", chem_path.exists(), str(chem_path))
        if not chem_path.exists():
            continue

        # Solver via RunConfig ausfuehren
        stdout, stderr = run_python_solver(wf["chemistry"], cfg=None)  # cfg already written

        check(f"{wf_id}: solver ran", len(stdout) > 0,
              f"stdout={len(stdout)}B, stderr={stderr[:100] if stderr else ''}")

        if not stdout:
            continue

        # Parsen
        parsed = parse_output(stdout)
        total_points = len(parsed["ifix"]) + len(parsed["results"])
        check(f"{wf_id}: output points > 0", total_points > 0,
              f"ifix={len(parsed['ifix'])}, result={len(parsed['results'])}")

        # Erwartungen pruefen
        check_expectations(wf_id, wf["expected"], parsed)

    # Aufraeumen
    _cleanup()

    print(f"\n{'='*70}")
    print(f"  Ergebnis: {passed} passed, {failed} failed")
    if errors:
        print(f"  Fehlgeschlagen:")
        for e in errors:
            print(f"    - {e}")
    print(f"{'='*70}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
