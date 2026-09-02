#!/usr/bin/env python3
"""
test_flow_units.py -- Tests fuer Massenfluss-Einheitenumrechnung.

Prueft:
  1. sccm -> mg/s Umrechnung fuer Xenon und Iod
  2. mg/s -> sccm Roundtrip
  3. Iod-Feedstock ist I2 (nicht I)
  4. Konsistenz mit bekannten Literaturwerten
  5. Integration mit RunConfig

Ausfuehrung:
    python test_flow_units.py
"""
from __future__ import annotations
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from flow_units import (
    sccm_to_mg_per_s, mg_per_s_to_sccm, sccm_to_pps,
    feedstock_mass_kg, FEEDSTOCK_SPECIES, SCCM_TO_PPS,
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


# ═══════════════════════════════════════════════════════════
# Test 1: Grundlegende Umrechnung
# ═══════════════════════════════════════════════════════════
def test_basic_conversion():
    print("\n--- Test 1: Grundlegende Umrechnung ---")

    # Xenon: 1 sccm = SCCM_TO_PPS particles/s * M_Xe = mg/s
    xe_mg = sccm_to_mg_per_s(1.0, "xenon")
    check("Xe 1sccm > 0 mg/s", xe_mg > 0, f"got {xe_mg}")
    # 1 sccm Xe = 4.478e17 * 2.180e-25 * 1e6 ~ 0.0976 mg/s
    check("Xe 1sccm ~ 0.098 mg/s", 0.09 < xe_mg < 0.11,
          f"got {xe_mg:.4f}")

    # Iod: 1 sccm = SCCM_TO_PPS * M_I2 * 1e6
    i2_mg = sccm_to_mg_per_s(1.0, "iodine")
    check("I2 1sccm > 0 mg/s", i2_mg > 0, f"got {i2_mg}")
    # 1 sccm I2 = 4.478e17 * 4.214e-25 * 1e6 ~ 0.1887 mg/s
    check("I2 1sccm ~ 0.189 mg/s", 0.18 < i2_mg < 0.20,
          f"got {i2_mg:.4f}")

    # I2 > Xe weil M_I2 > M_Xe
    check("I2 mg/s > Xe mg/s (gleiche sccm)", i2_mg > xe_mg,
          f"I2={i2_mg:.4f}, Xe={xe_mg:.4f}")


# ═══════════════════════════════════════════════════════════
# Test 2: Roundtrip
# ═══════════════════════════════════════════════════════════
def test_roundtrip():
    print("\n--- Test 2: Roundtrip sccm -> mg/s -> sccm ---")

    for gas in ("xenon", "iodine", "krypton", "argon"):
        for sccm in (0.5, 1.0, 3.0, 10.0):
            mg = sccm_to_mg_per_s(sccm, gas)
            sccm_back = mg_per_s_to_sccm(mg, gas)
            check(f"{gas} {sccm}sccm roundtrip",
                  abs(sccm_back - sccm) / sccm < 1e-10,
                  f"sccm={sccm}, mg={mg:.4f}, back={sccm_back:.6f}")


# ═══════════════════════════════════════════════════════════
# Test 3: Iod-Feedstock ist I2
# ═══════════════════════════════════════════════════════════
def test_iodine_feedstock():
    print("\n--- Test 3: Iod-Feedstock = I2 ---")

    check("Iod feedstock species is I2", FEEDSTOCK_SPECIES["iodine"] == "I2")

    M_I2 = feedstock_mass_kg("iodine")
    M_I = 2.1071711e-25  # Atomares Iod

    check("I2 mass ~ 2 * I mass",
          abs(M_I2 - 2 * M_I) / M_I2 < 0.001,
          f"M_I2={M_I2:.4e}, 2*M_I={2*M_I:.4e}")

    # Stelle sicher dass NICHT atomare Masse verwendet wird
    check("Feedstock mass > 4e-25 (I2, nicht I)",
          M_I2 > 4e-25,
          f"M={M_I2:.4e}")


# ═══════════════════════════════════════════════════════════
# Test 4: Literaturvergleich
# ═══════════════════════════════════════════════════════════
def test_literature():
    print("\n--- Test 4: Literaturvergleich ---")

    # Holste 2018: mdot = 0.5 mg/s Xe → Q0 ~ 5.1 sccm
    sccm_holste = mg_per_s_to_sccm(0.5, "xenon")
    check("Holste Xe 0.5mg/s ~ 5.1 sccm",
          4.5 < sccm_holste < 5.5,
          f"got {sccm_holste:.2f}")

    # NPT30: mdot = 0.056 mg/s I2 → Q0 ~ 0.297 sccm
    sccm_npt30 = mg_per_s_to_sccm(0.056, "iodine")
    check("NPT30 I2 0.056mg/s ~ 0.30 sccm",
          0.25 < sccm_npt30 < 0.35,
          f"got {sccm_npt30:.3f}")

    # Grondein: Q0 = 3 sccm I2 → mdot ~ 0.57 mg/s
    mg_grondein = sccm_to_mg_per_s(3.0, "iodine")
    check("Grondein I2 3sccm ~ 0.57 mg/s",
          0.5 < mg_grondein < 0.65,
          f"got {mg_grondein:.3f}")


# ═══════════════════════════════════════════════════════════
# Test 5: sccm_to_pps gasunabhaengig
# ═══════════════════════════════════════════════════════════
def test_pps_conversion():
    print("\n--- Test 5: sccm -> pps (gasunabhaengig) ---")

    pps = sccm_to_pps(1.0)
    check("1 sccm ~ 4.478e17 pps", abs(pps - SCCM_TO_PPS) < 1e10)

    # pps ist gasunabhaengig
    check("SCCM_TO_PPS value", abs(SCCM_TO_PPS - 4.477962312e17) < 1e10)


def main():
    global passed, failed
    test_basic_conversion()
    test_roundtrip()
    test_iodine_feedstock()
    test_literature()
    test_pps_conversion()

    print(f"\n{'='*60}")
    print(f"  Ergebnis: {passed} passed, {failed} failed")
    if errors:
        print(f"  Fehlgeschlagen: {', '.join(errors)}")
    print(f"{'='*60}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
