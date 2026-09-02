#!/usr/bin/env python3
"""
test_physics_data_viewer.py -- Tests for the physics data inspection window.

Tests:
  1. Cross section discovery for C++ packages
  2. Rate coefficient discovery for C++ and Python packages
  3. Data loading from CSV files
  4. Arrhenius computation for demo packages
  5. Data availability reporting per package
  6. Window creation (offscreen)
"""
from __future__ import annotations
import sys, os, json
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from physics_data_viewer import (
    discover_cross_sections, discover_rate_coefficients,
    _read_csv_2col, compute_arrhenius_curve,
    _build_tooltip, _build_details,
)

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


def test_cs_discovery():
    print("\n--- Test 1: Cross Section Discovery ---")
    biagi = SCRIPT_DIR / "cross_sections" / "xenon" / "biagi"
    if not biagi.exists():
        check("biagi dir exists", False)
        return

    items = discover_cross_sections(biagi)
    check("biagi has cross sections", len(items) > 0, f"found {len(items)}")

    labels = [i["label"] for i in items]
    check("has elastic", any("Elastic" in l for l in labels), str(labels[:3]))
    check("has ionization", any("Ionization" in l for l in labels))
    check("has excitation states", any("Excitation" in l for l in labels))
    check("excitation count > 10", sum(1 for l in labels if "Excitation" in l) > 10,
          f"found {sum(1 for l in labels if 'Excitation' in l)}")

    # No cross sections for Python packages
    iodine = SCRIPT_DIR / "chemistry" / "iodine_lafleur_v1"
    items_i = discover_cross_sections(iodine)
    check("iodine has no cross sections", len(items_i) == 0)


def test_rate_discovery():
    print("\n--- Test 2: Rate Coefficient Discovery ---")

    # C++ package
    biagi = SCRIPT_DIR / "cross_sections" / "xenon" / "biagi"
    if biagi.exists():
        items = discover_rate_coefficients(biagi, {"backend": "cpp"})
        check("biagi has rate tables", len(items) >= 3, f"found {len(items)}")
        labels = [i["label"] for i in items]
        check("has K_el", any("K_el" in l for l in labels))
        check("has K_iz", any("K_iz" in l for l in labels))
        check("has K_ex", any("K_ex" in l for l in labels))

    # Python package with tabulated rates
    iodine = SCRIPT_DIR / "chemistry" / "iodine_lafleur_v1"
    if iodine.exists():
        items = discover_rate_coefficients(iodine, {"backend": "python"})
        check("iodine has rate coefficients", len(items) >= 8, f"found {len(items)}")

    # Demo package (Arrhenius only)
    xe_simple = SCRIPT_DIR / "chemistry" / "xenon_simple"
    if xe_simple.exists():
        items = discover_rate_coefficients(xe_simple, {"backend": "python"})
        check("xenon_simple has rates", len(items) >= 2, f"found {len(items)}")
        reprs = [i.get("representation", "") for i in items]
        check("xenon_simple uses arrhenius", any("arrhenius" in r for r in reprs), str(reprs))


def test_csv_loading():
    print("\n--- Test 3: CSV Data Loading ---")

    # Cross section
    elastic = SCRIPT_DIR / "cross_sections" / "xenon" / "biagi" / "elastic.csv"
    if elastic.exists():
        result = _read_csv_2col(elastic)
        check("elastic.csv loads", result is not None)
        if result:
            x, y, meta = result
            check("elastic has data", len(x) > 50, f"points={len(x)}")
            check("energy > 0", all(v >= 0 for v in x))
            check("sigma > 0", all(v > 0 for v in y[:10]))

    # Rate table
    kiz = SCRIPT_DIR / "cross_sections" / "xenon" / "biagi" / "kiz_table.csv"
    if kiz.exists():
        result = _read_csv_2col(kiz)
        check("kiz_table.csv loads", result is not None)
        if result:
            x, y, meta = result
            check("kiz has data", len(x) > 100, f"points={len(x)}")

    # Iodine rate
    kiz_i2 = SCRIPT_DIR / "chemistry" / "iodine_lafleur_v1" / "rates" / "Kiz_I2.csv"
    if kiz_i2.exists():
        result = _read_csv_2col(kiz_i2)
        check("Kiz_I2.csv loads", result is not None)
        if result:
            x, y, meta = result
            check("Kiz_I2 has data", len(x) > 100, f"points={len(x)}")

    # Nonexistent file
    result = _read_csv_2col(Path("/nonexistent/file.csv"))
    check("missing file returns None", result is None)


def test_arrhenius():
    print("\n--- Test 4: Arrhenius Computation ---")
    Te, K = compute_arrhenius_curve(1.8e-13, 12.127, "arrhenius")
    check("arrhenius returns 200 points", len(Te) == 200)
    check("K > 0 at high Te", K[-1] > 0, f"K(20eV)={K[-1]:.2e}")
    check("K increases with Te", K[-1] > K[0])

    Te_c, K_c = compute_arrhenius_curve(1e-13, 0, "constant")
    check("constant returns flat", abs(K_c[0] - K_c[-1]) < 1e-20)
    check("constant value correct", abs(K_c[0] - 1e-13) < 1e-20)


def test_data_availability():
    print("\n--- Test 5: Data Availability Summary ---")
    packages = {
        "xenon_biagi": {"path": str(SCRIPT_DIR / "cross_sections/xenon/biagi"), "backend": "cpp"},
        "iodine_lafleur": {"path": str(SCRIPT_DIR / "chemistry/iodine_lafleur_v1"), "backend": "python"},
        "xenon_simple": {"path": str(SCRIPT_DIR / "chemistry/xenon_simple"), "backend": "python"},
    }
    for name, info in packages.items():
        p = Path(info["path"])
        if not p.exists():
            continue
        cs = discover_cross_sections(p) if info["backend"] == "cpp" else []
        rates = discover_rate_coefficients(p, info)
        print(f"  {name}: {len(cs)} cross sections, {len(rates)} rate coefficients")
        check(f"{name} has some data", len(cs) + len(rates) > 0)


def test_tooltips_and_details():
    print("\n--- Test 6: Tooltips and Details ---")
    biagi = SCRIPT_DIR / "cross_sections" / "xenon" / "biagi"
    if biagi.exists():
        items = discover_cross_sections(biagi)
        if items:
            tt = _build_tooltip(items[0])
            check("tooltip not empty", len(tt) > 10, f"len={len(tt)}")
            check("tooltip has Type", "Type:" in tt or "elastic" in tt.lower(), tt[:60])

            det = _build_details(items[0])
            check("details has HTML", "<b>" in det)
            check("details has data type", "cross_section" in det)

    # Rate with reaction string
    iodine = SCRIPT_DIR / "chemistry" / "iodine_lafleur_v1"
    if iodine.exists():
        items = discover_rate_coefficients(iodine, {"backend": "python"})
        if items:
            tt = _build_tooltip(items[0])
            check("iodine tooltip not empty", len(tt) > 5)
            det = _build_details(items[0])
            check("iodine details has HTML", "<b>" in det)

    # Arrhenius item
    arr_item = {"label": "test", "A": 1.8e-13, "E_a_eV": 12.127,
                "model": "arrhenius", "representation": "arrhenius",
                "data_type": "rate_coefficient", "process": "ionization"}
    det = _build_details(arr_item)
    check("arrhenius details has A", "1.8" in det, det[:100])


def test_metadata_richness():
    print("\n--- Test 7: Metadata Richness ---")
    biagi = SCRIPT_DIR / "cross_sections" / "xenon" / "biagi"
    if biagi.exists():
        items = discover_cross_sections(biagi)
        for item in items[:3]:
            check(f"{item['label']}: has process", bool(item.get("process")))
            check(f"{item['label']}: has unit_y", item.get("unit_y") == "m^2")

    iodine = SCRIPT_DIR / "chemistry" / "iodine_lafleur_v1"
    if iodine.exists():
        items = discover_rate_coefficients(iodine, {"backend": "python"})
        for item in items[:3]:
            check(f"{item['label']}: has process", bool(item.get("process")))
            check(f"{item['label']}: has unit_y", item.get("unit_y") == "m^3/s")


def test_window_creation():
    print("\n--- Test 8: Window Creation + Axis Controls ---")
    os.environ['QT_QPA_PLATFORM'] = 'offscreen'
    from PyQt6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    from physics_data_viewer import PhysicsDataWindow

    biagi_info = {
        "display_name": "Xenon (Biagi)", "gas": "xenon", "backend": "cpp",
        "database": "biagi", "path": str(SCRIPT_DIR / "cross_sections/xenon/biagi"),
    }
    w1 = PhysicsDataWindow(biagi_info)
    check("C++ window created", w1 is not None)
    check("C++ has cs list", w1._cs_list.count() > 0)
    check("C++ has rate list", w1._rate_list.count() > 0)
    check("C++ cs list enabled", w1._cs_list.isEnabled())
    check("has cs log-x checkbox", hasattr(w1, "_cs_log_x"))
    check("has cs log-y checkbox", hasattr(w1, "_cs_log_y"))
    check("cs log-y default checked", w1._cs_log_y.isChecked())
    check("cs log-x default unchecked", not w1._cs_log_x.isChecked())
    check("has details panel", hasattr(w1, "_details"))

    # Test axis toggle
    w1._cs_log_x.setChecked(True)
    check("cs log-x toggled", w1._cs_log_x.isChecked())
    w1._cs_log_x.setChecked(False)

    # Tooltips on list items
    if w1._cs_list.count() > 0:
        item = w1._cs_list.item(0)
        tt = item.toolTip()
        check("cs item has tooltip", len(tt) > 5, f"tooltip='{tt[:40]}'")

    iod_info = {
        "display_name": "Iodine (Lafleur V1)", "gas": "iodine", "backend": "python",
        "database": "", "path": str(SCRIPT_DIR / "chemistry/iodine_lafleur_v1"),
    }
    w2 = PhysicsDataWindow(iod_info)
    check("Python window created", w2 is not None)
    check("Python cs list disabled", not w2._cs_list.isEnabled())
    check("Python has rate list", w2._rate_list.count() > 0)
    if w2._rate_list.count() > 0:
        item = w2._rate_list.item(0)
        tt = item.toolTip()
        check("rate item has tooltip", len(tt) > 5, f"tooltip='{tt[:40]}'")


def main():
    global passed, failed
    test_cs_discovery()
    test_rate_discovery()
    test_csv_loading()
    test_arrhenius()
    test_data_availability()
    test_tooltips_and_details()
    test_metadata_richness()
    test_window_creation()

    print(f"\n{'='*60}")
    print(f"  Result: {passed} passed, {failed} failed")
    if errors:
        print(f"  Failed: {', '.join(errors)}")
    print(f"{'='*60}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
