"""
package_registry.py -- Zentrale Paketverwaltung und Backend-Routing.

Erkennt verfuegbare Chemiepakete aus zwei Quellen:
1. chemistry/<name>/chemistry.json  (generische Python-Pakete)
2. cross_sections/<gas>/<db>/db_info.json  (C++-kompatible Cross-Section-Pakete)

Stellt eine einheitliche PackageInfo-Schnittstelle bereit, die von der GUI
und dem Backend-Routing verwendet wird.

Architektur:
- PackageInfo: Einheitliche Metadaten fuer alle Pakettypen
- discover_packages(): Automatische Erkennung aller Pakete
- resolve_backend(): Entscheidet C++ oder Python basierend auf Metadaten
- generate_cpp_config(): Erzeugt params.txt fuer C++-Backend
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional


class BackendType(str, Enum):
    CPP = "cpp"           # Bestehender C++-Solver (chabert)
    PYTHON = "python"     # Generischer Python-Solver (generic_solver.py)


class PackageStatus(str, Enum):
    PRODUCTION = "production"    # Validiert, produktionsreif
    EXPERIMENTAL = "experimental"  # Funktional, nicht vollstaendig validiert
    DEMO = "demo"                # Architektur-Demo, nicht fuer Vorhersagen


@dataclass
class PackageInfo:
    """Einheitliche Metadaten fuer ein Chemie-/Cross-Section-Paket."""
    id: str                              # Eindeutiger Schluessel, z.B. "xenon_biagi_legacy"
    display_name: str                    # Anzeigename in der GUI
    gas: str                             # Gasspezies, z.B. "xenon", "iodine"
    database: str = ""                   # Datenbank, z.B. "biagi", "hayashi"
    data_origin: str = ""                # "experimental", "theoretical", "mixed"
    description: str = ""                # Kurzbeschreibung
    backend: BackendType = BackendType.CPP
    status: PackageStatus = PackageStatus.PRODUCTION
    version: str = ""
    path: Path = field(default_factory=Path)  # Pfad zum Paketverzeichnis

    # C++-spezifisch
    rate_models: list[int] = field(default_factory=lambda: [0, 1, 2])
    # rate_model-Werte die dieses Paket unterstuetzt (0=Legacy, 1=Conservative, 2=Full)

    # Python-spezifisch
    chemistry_json: str = ""             # Pfad zur chemistry.json relativ zum Paket

    @property
    def is_production(self) -> bool:
        return self.status == PackageStatus.PRODUCTION

    @property
    def is_cpp(self) -> bool:
        return self.backend == BackendType.CPP

    @property
    def is_python(self) -> bool:
        return self.backend == BackendType.PYTHON

    @property
    def status_label(self) -> str:
        return {
            PackageStatus.PRODUCTION: "",
            PackageStatus.EXPERIMENTAL: " [exp]",
            PackageStatus.DEMO: " [demo]",
        }.get(self.status, "")


def discover_packages(base_dir: Path) -> list[PackageInfo]:
    """Erkenne alle verfuegbaren Chemiepakete.

    Scannt zwei Quellen:
    1. cross_sections/<gas>/<db>/db_info.json → C++-kompatible Pakete
    2. chemistry/<name>/chemistry.json → Generische Python-Pakete
    """
    packages = []

    # ── Quelle 1: C++ Cross-Section-Pakete ────────────────────
    cs_dir = base_dir / "cross_sections"
    if cs_dir.is_dir():
        for gas_dir in sorted(cs_dir.iterdir()):
            if not gas_dir.is_dir() or gas_dir.name in ("tests",):
                continue
            for db_dir in sorted(gas_dir.iterdir()):
                info_path = db_dir / "db_info.json"
                if not db_dir.is_dir() or not info_path.exists():
                    continue
                try:
                    info = _load_json(info_path)
                    pkg = PackageInfo(
                        id=f"{gas_dir.name}_{db_dir.name}",
                        display_name=f"{gas_dir.name.title()} / {info.get('display_name', db_dir.name)}",
                        gas=gas_dir.name,
                        database=info.get("key", db_dir.name),
                        data_origin=info.get("type", ""),
                        description=info.get("comment", ""),
                        backend=BackendType.CPP,
                        status=PackageStatus.PRODUCTION,
                        version=info.get("version", ""),
                        path=db_dir,
                        rate_models=[0, 1, 2],
                    )
                    packages.append(pkg)
                except (json.JSONDecodeError, OSError, KeyError):
                    pass

    # ── Quelle 2: Python-Chemiepakete ─────────────────────────
    chem_dir = base_dir / "chemistry"
    if chem_dir.is_dir():
        for pkg_dir in sorted(chem_dir.iterdir()):
            chem_json = pkg_dir / "chemistry.json"
            if not pkg_dir.is_dir() or not chem_json.exists():
                continue
            try:
                info = _load_json(chem_json)
                name = info.get("name", pkg_dir.name)

                # Status aus Beschreibung ableiten
                desc = info.get("description", "")
                status = PackageStatus.DEMO
                if "demo" in desc.lower() or "architektur" in desc.lower():
                    status = PackageStatus.DEMO
                elif "validated" in desc.lower() or "production" in desc.lower():
                    status = PackageStatus.PRODUCTION
                else:
                    status = PackageStatus.EXPERIMENTAL

                # Gas aus Spezies ableiten
                species = info.get("species", [])
                gas = _infer_gas(name, species)

                pkg = PackageInfo(
                    id=f"py_{pkg_dir.name}",
                    display_name=f"{name} (Python)",
                    gas=gas,
                    data_origin="mixed",
                    description=desc[:120],
                    backend=BackendType.PYTHON,
                    status=status,
                    version="",
                    path=pkg_dir,
                    chemistry_json=str(chem_json),
                )
                packages.append(pkg)
            except (json.JSONDecodeError, OSError, KeyError):
                pass

    return packages


def _load_json(path: Path) -> dict:
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _infer_gas(name: str, species: list) -> str:
    """Versuche den Gasnamen aus Paketnamen oder Spezies abzuleiten."""
    name_lower = name.lower()
    for gas in ("xenon", "krypton", "argon", "iodine", "oxygen", "nitrogen"):
        if gas in name_lower:
            return gas
    # Fallback: aus Spezies-IDs
    sp_ids = {s.get("id", "").lower() for s in species}
    if "xe" in sp_ids or "xe+" in sp_ids:
        return "xenon"
    if "i2" in sp_ids or "i+" in sp_ids:
        return "iodine"
    return "unknown"


def resolve_backend(pkg: PackageInfo) -> BackendType:
    """Bestimme das Backend fuer ein Paket."""
    return pkg.backend


def get_default_package(packages: list[PackageInfo]) -> Optional[PackageInfo]:
    """Finde das Standard-Paket (Xenon/Biagi C++ Production)."""
    for p in packages:
        if p.gas == "xenon" and p.database == "biagi" and p.is_cpp:
            return p
    # Fallback: erstes Produktionspaket
    for p in packages:
        if p.is_production:
            return p
    return packages[0] if packages else None


# ═══════════════════════════════════════════════════════════════
# Backend-Konfiguration
# ═══════════════════════════════════════════════════════════════

def generate_cpp_config(
    pkg: PackageInfo,
    geometry: dict[str, float],
    sweep: dict[str, float],
    solve_mode: int = 1,
    rate_model: int = 0,
) -> str:
    """Erzeuge params.txt-Inhalt fuer C++-Backend."""
    lines = ["# Auto-generated by package_registry"]
    for k, v in geometry.items():
        lines.append(f"{k} {v}")
    for k, v in sweep.items():
        lines.append(f"{k} {v}")
    lines.append(f"method 4")
    lines.append(f"solve_mode {solve_mode}")
    lines.append(f"gas_species {pkg.gas}")
    lines.append(f"cs_database {pkg.database}")
    lines.append(f"rate_model {rate_model}")
    lines.append("use_paper_kel 1")
    return "\n".join(lines) + "\n"
