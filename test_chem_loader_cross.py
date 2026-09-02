#!/usr/bin/env python3
"""
test_chem_loader_cross.py -- Vergleicht den C++-Lader gegen den Python-Lader.

Beide Rechenwege sollen dieselbe Chemie aus derselben Datei lesen. Dieser Test
laesst das uebersetzte Testprogramm jedes vorhandene Chemiepaket ausgeben und
prueft Spezies, Reaktionen und Ratenwerte gegen das, was plasma_chemistry.py
aus derselben Datei macht. Weicht eine Seite ab, ist die gemeinsame Quelle
nicht mehr gemeinsam -- genau das soll hier auffallen.

Aufruf: python test_chem_loader_cross.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from plasma_chemistry import load_chemistry, RateModel  # noqa: E402

#: Stuetzstellen, an denen die Raten verglichen werden. Dieselbe Reihe wie im
#: C++-Testprogramm, in derselben Reihenfolge.
TE_PUNKTE = [0.7, 1.0, 2.0, 3.75, 5.0, 8.0, 12.0, 30.0]

#: Relative Schranke fuer den Ratenvergleich. Beide Seiten interpolieren
#: linear zwischen denselben Stuetzstellen, es bleibt reines Rundungsrauschen.
REL_TOL = 1e-9

passed = 0
failed = 0
fehlernamen: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    global passed, failed
    if cond:
        passed += 1
    else:
        failed += 1
        print(f"  FAIL: {name} -- {detail}")
        fehlernamen.append(name)


def testprogramm() -> Path | None:
    for name in ("test_chem_loader.exe", "test_chem_loader"):
        p = SCRIPT_DIR / name
        if p.exists():
            return p
    return None


def pakete() -> list[Path]:
    """Alle Chemiepakete, alphabetisch."""
    return sorted((SCRIPT_DIR / "chemistry").glob("*/chemistry.json"))


def dump_lesen(prog: Path, paket: Path) -> dict:
    """Ausgabe des C++-Laders einlesen."""
    rel = paket.relative_to(SCRIPT_DIR).as_posix()
    r = subprocess.run([str(prog), "--dump", rel], cwd=SCRIPT_DIR,
                       capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise RuntimeError(f"C++-Lader meldet Fehler fuer {rel}:\n{r.stdout}{r.stderr}")

    out: dict = {"species": [], "reactions": []}
    for zeile in r.stdout.splitlines():
        f = zeile.split("\t")
        if f[0] == "PKG":
            out["name"] = f[1]
            out["wall_temperature_K"] = float(f[2])
            out["sigma_i"] = float(f[3])
        elif f[0] == "SPECIES":
            out["species"].append({
                "id": f[1], "type": f[2], "mass_kg": float(f[3]),
                "charge": int(f[4]), "thermal_conductivity": float(f[5]),
                "is_feedstock": f[6] == "1", "is_beam_extracted": f[7] == "1",
                "wall_products": f[8] if len(f) > 8 else "",
            })
        elif f[0] == "REACTION":
            out["reactions"].append({
                "id": f[1], "type": f[2],
                "reactants": f[3], "products": f[4],
                "energy_eV": float(f[5]),
                "is_electron_impact": f[6] == "1",
                "elastic_heating": f[7] == "1",
                "nu_m": f[8] == "1",
                "surface_gamma": float(f[9]),
                "rate_model": f[10],
                "K": [float(x) for x in f[11:11 + len(TE_PUNKTE)]],
            })
    return out


def stoich_text(d: dict) -> str:
    """Stoechiometrie als sortierte Zeichenkette, wie sie C++ ausgibt."""
    return ",".join(f"{k}:{v}" for k, v in sorted(d.items()))


def nahe(a: float, b: float) -> bool:
    if a == b:
        return True
    skala = max(abs(a), abs(b))
    return abs(a - b) <= REL_TOL * skala


def vergleiche(paket: Path, cpp: dict) -> None:
    """Ein Paket auf beiden Wegen laden und gegenueberstellen."""
    kurz = paket.parent.name
    py = load_chemistry(paket)  # leitet fehlende Wandprodukte selbst ab

    check(f"{kurz}: Name", cpp["name"] == py.name, f"{cpp['name']} vs {py.name}")
    check(f"{kurz}: Wandtemperatur",
          nahe(cpp["wall_temperature_K"], py.wall_temperature_K),
          f"{cpp['wall_temperature_K']} vs {py.wall_temperature_K}")
    check(f"{kurz}: Stossquerschnitt", nahe(cpp["sigma_i"], py.sigma_i),
          f"{cpp['sigma_i']} vs {py.sigma_i}")

    # ── Spezies: Reihenfolge und Inhalt ──────────────────────
    py_sp = py.heavy_species  # ohne Elektron, in Dateireihenfolge
    check(f"{kurz}: Speziesanzahl", len(cpp["species"]) == len(py_sp),
          f"{len(cpp['species'])} vs {len(py_sp)}")
    for c, p in zip(cpp["species"], py_sp):
        check(f"{kurz}/{p.id}: Reihenfolge", c["id"] == p.id, f"{c['id']} vs {p.id}")
        check(f"{kurz}/{p.id}: Typ", c["type"] == p.species_type.value,
              f"{c['type']} vs {p.species_type.value}")
        check(f"{kurz}/{p.id}: Masse", nahe(c["mass_kg"], p.mass_kg),
              f"{c['mass_kg']} vs {p.mass_kg}")
        check(f"{kurz}/{p.id}: Ladung", c["charge"] == p.charge,
              f"{c['charge']} vs {p.charge}")
        check(f"{kurz}/{p.id}: Waermeleitfaehigkeit",
              nahe(c["thermal_conductivity"], p.thermal_conductivity),
              f"{c['thermal_conductivity']} vs {p.thermal_conductivity}")
        check(f"{kurz}/{p.id}: Feedstock", c["is_feedstock"] == p.is_feedstock)
        check(f"{kurz}/{p.id}: Extraktion",
              c["is_beam_extracted"] == p.is_beam_extracted)
        check(f"{kurz}/{p.id}: Wandprodukte",
              c["wall_products"] == stoich_text(p.wall_products),
              f"{c['wall_products']} vs {stoich_text(p.wall_products)}")

    # ── Reaktionen: Struktur und Raten ───────────────────────
    check(f"{kurz}: Reaktionsanzahl", len(cpp["reactions"]) == len(py.reactions),
          f"{len(cpp['reactions'])} vs {len(py.reactions)}")
    for c, p in zip(cpp["reactions"], py.reactions):
        check(f"{kurz}/{p.id}: Kennung", c["id"] == p.id, f"{c['id']} vs {p.id}")
        check(f"{kurz}/{p.id}: Art", c["type"] == p.reaction_type.value,
              f"{c['type']} vs {p.reaction_type.value}")
        check(f"{kurz}/{p.id}: Edukte", c["reactants"] == stoich_text(p.reactants),
              f"{c['reactants']} vs {stoich_text(p.reactants)}")
        check(f"{kurz}/{p.id}: Produkte", c["products"] == stoich_text(p.products),
              f"{c['products']} vs {stoich_text(p.products)}")
        check(f"{kurz}/{p.id}: Energie", nahe(c["energy_eV"], p.energy_eV),
              f"{c['energy_eV']} vs {p.energy_eV}")
        check(f"{kurz}/{p.id}: Elektronenstoss",
              c["is_electron_impact"] == p.is_electron_impact)
        check(f"{kurz}/{p.id}: Gasheizung",
              c["elastic_heating"] == p.contributes_to_elastic_heating)
        check(f"{kurz}/{p.id}: Stossfrequenz", c["nu_m"] == p.contributes_to_nu_m)
        check(f"{kurz}/{p.id}: Oberflaechenkoeffizient",
              nahe(c["surface_gamma"], p.surface_gamma),
              f"{c['surface_gamma']} vs {p.surface_gamma}")
        check(f"{kurz}/{p.id}: Ratenmodell", c["rate_model"] == p.rate.model.value,
              f"{c['rate_model']} vs {p.rate.model.value}")

        for te, k_cpp in zip(TE_PUNKTE, c["K"]):
            k_py = p.rate.evaluate(te)
            check(f"{kurz}/{p.id}: K bei {te} eV", nahe(k_cpp, k_py),
                  f"{k_cpp:.6e} vs {k_py:.6e}")

    tabelliert = sum(1 for r in py.reactions if r.rate.model == RateModel.TABULATED)
    print(f"  {kurz}: {len(py_sp)} Spezies, {len(py.reactions)} Reaktionen, "
          f"{tabelliert} davon tabelliert -- beide Wege gleich"
          if not fehlernamen else f"  {kurz}: geprueft")


def main() -> int:
    prog = testprogramm()
    if prog is None:
        print("test_chem_loader (C++) fehlt. Mit 'python build.py --tests' uebersetzen.",
              file=sys.stderr)
        return 1

    liste = pakete()
    if not liste:
        print("Keine Chemiepakete gefunden.", file=sys.stderr)
        return 1

    print(f"Vergleiche {len(liste)} Chemiepakete auf beiden Wegen:\n")
    for paket in liste:
        try:
            cpp = dump_lesen(prog, paket)
        except Exception as ex:  # noqa: BLE001 -- Grund soll im Bericht stehen
            check(f"{paket.parent.name}: C++-Lader laeuft", False, str(ex))
            continue
        vergleiche(paket, cpp)

    print(f"\n{'=' * 60}")
    print(f"  Ergebnis: {passed} bestanden, {failed} fehlgeschlagen")
    if fehlernamen:
        print(f"  Fehlgeschlagen: {', '.join(fehlernamen[:10])}"
              + (" ..." if len(fehlernamen) > 10 else ""))
    print(f"{'=' * 60}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
