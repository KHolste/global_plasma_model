#!/usr/bin/env python3
"""
test_registry.py -- Tests fuer Paketverwaltung und Backend-Routing.

Testet:
1. Discovery findet C++ und Python-Pakete
2. Metadaten korrekt geladen
3. Backend-Routing (C++ fuer Xenon/Biagi, Python fuer Iod)
4. Default-Paket ist Xenon/Biagi
5. Rueckwaertskompatibilitaet (alte Config-Felder mappbar)
6. Fehlerfaelle: fehlende Metadaten, unbekannte Pakete
7. Config-Generierung fuer C++-Backend
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from package_registry import (
    PackageInfo, PackageStatus, BackendType,
    discover_packages, resolve_backend, get_default_package,
    generate_cpp_config,
)

BASE = Path(__file__).resolve().parent

passed = 0
failed = 0


def test(name, condition, detail=""):
    global passed, failed
    if condition:
        passed += 1
        print(f"  PASS: {name}")
    else:
        failed += 1
        print(f"  FAIL: {name} -- {detail}")


# ═══════════════════════════════════════════════════════════════
# 1. Discovery
# ═══════════════════════════════════════════════════════════════
print("\n=== Test 1: Paket-Discovery ===")

packages = discover_packages(BASE)
test("Pakete gefunden", len(packages) >= 4, f"got {len(packages)}")

# Erwartete Pakete
ids = {p.id for p in packages}
test("xenon_biagi vorhanden", "xenon_biagi" in ids, f"ids={ids}")
test("xenon_hayashi vorhanden", "xenon_hayashi" in ids)
test("py_xenon_simple vorhanden", "py_xenon_simple" in ids)
test("py_iodine_grondein vorhanden", "py_iodine_grondein" in ids)

# ═══════════════════════════════════════════════════════════════
# 2. Metadaten
# ═══════════════════════════════════════════════════════════════
print("\n=== Test 2: Paket-Metadaten ===")

xe_biagi = next(p for p in packages if p.id == "xenon_biagi")
test("xe_biagi Gas = xenon", xe_biagi.gas == "xenon")
test("xe_biagi DB = biagi", xe_biagi.database == "biagi")
test("xe_biagi Backend = C++", xe_biagi.backend == BackendType.CPP)
test("xe_biagi Status = production", xe_biagi.status == PackageStatus.PRODUCTION)
test("xe_biagi Pfad existiert", xe_biagi.path.exists())

xe_hay = next(p for p in packages if p.id == "xenon_hayashi")
test("xe_hayashi Gas = xenon", xe_hay.gas == "xenon")
test("xe_hayashi DB = hayashi", xe_hay.database == "hayashi")

py_iod = next(p for p in packages if p.id == "py_iodine_grondein")
test("iod Gas = iodine", py_iod.gas == "iodine")
test("iod Backend = Python", py_iod.backend == BackendType.PYTHON)
test("iod Status = demo", py_iod.status == PackageStatus.DEMO)
test("iod display_name", "Iodine" in py_iod.display_name)

py_xe = next(p for p in packages if p.id == "py_xenon_simple")
test("py_xe Gas = xenon", py_xe.gas == "xenon")
test("py_xe Backend = Python", py_xe.backend == BackendType.PYTHON)

# ═══════════════════════════════════════════════════════════════
# 3. Backend-Routing
# ═══════════════════════════════════════════════════════════════
print("\n=== Test 3: Backend-Routing ===")

test("Xenon/Biagi -> C++", resolve_backend(xe_biagi) == BackendType.CPP)
test("Xenon/Hayashi -> C++", resolve_backend(xe_hay) == BackendType.CPP)
test("Iod/Python -> Python", resolve_backend(py_iod) == BackendType.PYTHON)
test("Xenon/Python -> Python", resolve_backend(py_xe) == BackendType.PYTHON)

# ═══════════════════════════════════════════════════════════════
# 4. Default-Paket
# ═══════════════════════════════════════════════════════════════
print("\n=== Test 4: Default-Paket ===")

default = get_default_package(packages)
test("Default ist xenon_biagi", default is not None and default.id == "xenon_biagi")

# ═══════════════════════════════════════════════════════════════
# 5. Rueckwaertskompatibilitaet
# ═══════════════════════════════════════════════════════════════
print("\n=== Test 5: Rueckwaertskompatibilitaet ===")

# Alte Config-Felder (gas_species, cs_database, rate_model) muessen auf Pakete mappbar sein
old_gas = "xenon"
old_db = "biagi"
matches = [p for p in packages if p.gas == old_gas and p.database == old_db and p.is_cpp]
test("Alte Config (xenon/biagi) findet Paket", len(matches) == 1)
test("Alte Config mapped auf C++ Backend", matches[0].backend == BackendType.CPP if matches else False)

old_gas2 = "xenon"
old_db2 = "hayashi"
matches2 = [p for p in packages if p.gas == old_gas2 and p.database == old_db2]
test("Alte Config (xenon/hayashi) findet Paket", len(matches2) == 1)

# ═══════════════════════════════════════════════════════════════
# 6. Fehlerfaelle
# ═══════════════════════════════════════════════════════════════
print("\n=== Test 6: Fehlerfaelle ===")

# Nicht-existierendes Verzeichnis
empty_pkgs = discover_packages(Path("/nonexistent"))
test("Leeres Verzeichnis -> leere Liste", len(empty_pkgs) == 0)

# PackageInfo mit ungueltigem Status
import tempfile, os
tmp = Path(tempfile.mkdtemp())
bad_chem = tmp / "chemistry" / "broken"
bad_chem.mkdir(parents=True)
(bad_chem / "chemistry.json").write_text('{"not_valid": true}', encoding="utf-8")
pkgs_tmp = discover_packages(tmp)
test("Kaputtes JSON -> kein Absturz", True)  # Darf nicht crashen
(bad_chem / "chemistry.json").unlink()
bad_chem.rmdir()
(tmp / "chemistry").rmdir()
tmp.rmdir()

# ═══════════════════════════════════════════════════════════════
# 7. Config-Generierung
# ═══════════════════════════════════════════════════════════════
print("\n=== Test 7: Config-Generierung ===")

geom = {"R": 0.02, "L": 0.04, "betai": 0.5}
sweep = {"Q0sccm_start": 0.3, "Q0sccm_step": 0.01, "jjmax": 50}
cfg = generate_cpp_config(xe_biagi, geom, sweep, solve_mode=2, rate_model=2)

test("Config enthaelt gas_species", "gas_species xenon" in cfg)
test("Config enthaelt cs_database", "cs_database biagi" in cfg)
test("Config enthaelt rate_model", "rate_model 2" in cfg)
test("Config enthaelt solve_mode", "solve_mode 2" in cfg)
test("Config enthaelt R", "R 0.02" in cfg)
test("Config enthaelt jjmax", "jjmax 50" in cfg)

# ═══════════════════════════════════════════════════════════════
# Zusammenfassung
# ═══════════════════════════════════════════════════════════════
print(f"\n{'='*50}")
print(f"  Tests: {passed + failed} | Passed: {passed} | Failed: {failed}")
print(f"{'='*50}")

sys.exit(1 if failed > 0 else 0)
