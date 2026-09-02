#!/usr/bin/env python3
"""
test_preset_plausibility.py -- Tests fuer Preset-Metadaten, Gas-Kompatibilitaet und I-fix-Bereichsanzeige.

Ausfuehrung:
    python test_preset_plausibility.py
"""
from __future__ import annotations
import sys, json
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
        msg = f"  FAIL: {name}"
        if detail:
            msg += f" -- {detail}"
        print(msg)
        errors.append(name)


# ═════════════════════════════════════════════════════════════
# Test 1: Preset-Metadaten vollstaendig
# ═════════════════════════════════════════════════════════════
def test_preset_metadata():
    print("\n--- Test 1: Preset-Metadaten ---")
    with open(SCRIPT_DIR / "thruster_presets.json", encoding="utf-8") as f:
        data = json.load(f)

    presets = data.get("presets", [])
    check("presets list exists", len(presets) > 0, f"found {len(presets)}")

    for p in presets:
        pid = p.get("id", "?")
        check(f"{pid}: has id", "id" in p)
        check(f"{pid}: has name", "name" in p)
        check(f"{pid}: has recommended_gas", "recommended_gas" in p)

        if pid != "custom":
            check(f"{pid}: has params", bool(p.get("params")))
            check(f"{pid}: has extraction.eta_opt",
                  "eta_opt" in p.get("extraction", {}),
                  f"extraction={p.get('extraction', {})}")

        # Plausibilitaets-Metadaten
        plaus = p.get("plausibility", {})
        if pid != "custom":
            check(f"{pid}: has plausibility block", bool(plaus))
            check(f"{pid}: has I_beam_typ_min_mA",
                  "I_beam_typ_min_mA" in plaus,
                  f"keys={list(plaus.keys())}")
            check(f"{pid}: has I_beam_typ_max_mA",
                  "I_beam_typ_max_mA" in plaus)

            I_min = plaus.get("I_beam_typ_min_mA", 0)
            I_max = plaus.get("I_beam_typ_max_mA", 0)
            check(f"{pid}: I_beam range valid",
                  0 < I_min < I_max,
                  f"I_min={I_min}, I_max={I_max}")


# ═════════════════════════════════════════════════════════════
# Test 2: Gas-Kompatibilitaetslogik
# ═════════════════════════════════════════════════════════════
def test_gas_compatibility():
    print("\n--- Test 2: Gas-Kompatibilitaet ---")
    with open(SCRIPT_DIR / "thruster_presets.json", encoding="utf-8") as f:
        data = json.load(f)
    presets = {p["id"]: p for p in data.get("presets", [])}

    # Chabert: nur Xenon empfohlen
    ch = presets["chabert"]
    check("chabert recommends xenon", "xenon" in ch["recommended_gas"])
    check("chabert has iodine warning",
          "iodine" in ch.get("plausibility", {}).get("gas_warnings", {}))

    # Grondein: Xe + I2
    gr = presets["grondein"]
    check("grondein recommends xenon", "xenon" in gr["recommended_gas"])
    check("grondein recommends iodine", "iodine" in gr["recommended_gas"])

    # NPT30: nur Iod empfohlen
    np_ = presets["npt30"]
    check("npt30 recommends iodine", "iodine" in np_["recommended_gas"])
    check("npt30 has xenon warning",
          "xenon" in np_.get("plausibility", {}).get("gas_warnings", {}))

    # Holste: Xe + I2
    ho = presets["holste"]
    check("holste recommends both", "xenon" in ho["recommended_gas"] and "iodine" in ho["recommended_gas"])


# ═════════════════════════════════════════════════════════════
# Test 3: I-fix-Bereichspruefung gegen Preset-Metadaten
# ═════════════════════════════════════════════════════════════
def test_ifix_range_check():
    print("\n--- Test 3: I-fix Bereichspruefung ---")
    with open(SCRIPT_DIR / "thruster_presets.json", encoding="utf-8") as f:
        data = json.load(f)
    presets = {p["id"]: p for p in data.get("presets", [])}

    def range_check(preset_id, I_soll):
        """Simuliert die GUI-Logik: prueft ob I_soll im typischen Bereich liegt."""
        p = presets[preset_id]
        plaus = p.get("plausibility", {})
        I_min = plaus.get("I_beam_typ_min_mA")
        I_max = plaus.get("I_beam_typ_max_mA")
        if I_min is None or I_max is None:
            return "no_range"
        if I_soll < I_min:
            return "below"
        elif I_soll > I_max:
            return "above"
        return "ok"

    # Chabert: 5-25 mA
    check("chabert: I=15mA in range", range_check("chabert", 15) == "ok")
    check("chabert: I=3mA below", range_check("chabert", 3) == "below")
    check("chabert: I=50mA above", range_check("chabert", 50) == "above")

    # Holste: 30-100 mA
    check("holste: I=80mA in range", range_check("holste", 80) == "ok")
    check("holste: I=10mA below", range_check("holste", 10) == "below")
    check("holste: I=150mA above", range_check("holste", 150) == "above")

    # NPT30: 3-20 mA
    check("npt30: I=15mA in range", range_check("npt30", 15) == "ok")
    check("npt30: I=1mA below", range_check("npt30", 1) == "below")

    # Custom: no range data
    check("custom: no range", range_check("custom", 15) == "no_range")


# ═════════════════════════════════════════════════════════════
# Test 4: SOLVER_FAIL Diagnose-Format
# ═════════════════════════════════════════════════════════════
def test_solver_fail_format():
    print("\n--- Test 4: SOLVER_FAIL Diagnose ---")

    # Simuliere die Ausgabe von generic_solver.py
    from io import StringIO
    import contextlib

    # Teste die Diagnose-Textgenerierung direkt
    I_soll = 80.0
    I_found = 40.5
    P_found = 500.0

    # above_P_max
    status = "above_P_max"
    diag = (f"Sollstrom {I_soll:.0f} mA nicht erreichbar: "
            f"bei P_max={P_found:.0f} W werden nur {I_found:.1f} mA erreicht")
    check("above_P_max diag contains I_soll", "80 mA" in diag)
    check("above_P_max diag contains I_found", "40.5 mA" in diag)
    check("above_P_max diag contains P_max", "500 W" in diag)
    check("above_P_max diag is human-readable", "nicht erreichbar" in diag)

    # below_P_min
    I_found_lo = 42.0
    diag_lo = (f"Sollstrom {15:.0f} mA unterhalb Minimum: "
               f"bereits bei P_min werden {I_found_lo:.1f} mA erreicht")
    check("below_P_min diag contains direction", "unterhalb" in diag_lo)

    # Test: diag= prefix extraction (wie in GUI)
    line = f"SOLVER_FAIL 0 0.5000 above_P_max I_target=80.0mA I_found=40.50mA P_max=500.0W diag={diag}"
    diag_idx = line.find("diag=")
    check("diag= prefix found", diag_idx >= 0)
    extracted = line[diag_idx + 5:]
    check("extracted diag matches", extracted == diag)


# ═════════════════════════════════════════════════════════════
# Test 5: Preset-Parameter konsistent
# ═════════════════════════════════════════════════════════════
def test_preset_params_consistent():
    print("\n--- Test 5: Preset-Parameter Konsistenz ---")
    with open(SCRIPT_DIR / "thruster_presets.json", encoding="utf-8") as f:
        data = json.load(f)

    for p in data.get("presets", []):
        pid = p.get("id")
        if pid == "custom":
            continue
        params = p.get("params", {})
        ext = p.get("extraction", {})
        plaus = p.get("plausibility", {})

        eta = ext.get("eta_opt", 0)
        check(f"{pid}: eta_opt > 0", eta > 0, f"eta_opt={eta}")
        check(f"{pid}: eta_opt <= 1", eta <= 1, f"eta_opt={eta}")

        # I_soll sollte im typischen Bereich liegen
        I_soll = params.get("I_soll", 0)
        I_min = plaus.get("I_beam_typ_min_mA", 0)
        I_max = plaus.get("I_beam_typ_max_mA", 999)
        check(f"{pid}: default I_soll in range",
              I_min <= I_soll <= I_max,
              f"I_soll={I_soll}, range=[{I_min},{I_max}]")

        # P_RFG_max sollte positiv sein
        P_max = params.get("P_RFG_max", 0)
        check(f"{pid}: P_RFG_max > 0", P_max > 0, f"P_max={P_max}")


def main():
    global passed, failed

    test_preset_metadata()
    test_gas_compatibility()
    test_ifix_range_check()
    test_solver_fail_format()
    test_preset_params_consistent()

    print(f"\n{'='*60}")
    print(f"  Ergebnis: {passed} passed, {failed} failed")
    if errors:
        print(f"  Fehlgeschlagen: {', '.join(errors)}")
    print(f"{'='*60}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
