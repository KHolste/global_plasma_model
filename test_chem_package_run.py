#!/usr/bin/env python3
"""
test_chem_package_run.py -- Der Rechenkern rechnet ueber ein Chemiepaket.

Prueft den angeschlossenen Weg von der Konfiguration bis zum Ergebnis:
nennt die Konfiguration ein Paket, laeuft der generische Loeser darueber und
gibt die Dichte jeder Spezies einzeln aus. Nennt sie keines oder ein
fehlerhaftes, bleibt es beim fest verdrahteten Weg -- ein unbrauchbares Paket
darf einen Lauf nicht verhindern.

Aufruf: python test_chem_package_run.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PARAMS = SCRIPT_DIR / "chem_run_params.txt"

BASIS = """R 0.02
L 0.04
betai 0.5
betag 0.05145
frequency 2.5e6
Nw 6.0
R_ohm 0.36
Rc 0.02
lc 0.04
Vgrid 1500.0
sgrid 0.001
P_RFG 18.0
P_RFG_max 80.0
Q0sccm_start 0.40
Q0sccm_step 0.02
jjmax 2
solve_mode 2
gas_species xenon
rate_model 0
"""

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


def kern() -> Path | None:
    for name in ("chabert.exe", "chabert"):
        p = SCRIPT_DIR / name
        if p.exists():
            return p
    return None


def lauf(prog: Path, zusatz: str = "") -> tuple[str, str]:
    PARAMS.write_text(BASIS + zusatz, encoding="utf-8")
    r = subprocess.run([str(prog), PARAMS.name], cwd=SCRIPT_DIR,
                       capture_output=True, text=True, timeout=600)
    return r.stdout, r.stderr


def zeilen(text: str, marke: str) -> list[list[str]]:
    return [z.split() for z in text.splitlines() if z.startswith(marke + " ")]


def test_mit_paket(prog: Path) -> None:
    out, _ = lauf(prog, "chemistry_package xenon_simple\n")

    check("Weg ist der generische", "SOLVER_PATH generic" in out)
    kopf = zeilen(out, "CHEM_PACKAGE")
    check("Paket wird genannt", bool(kopf) and "xenon_simple" in kopf[0][1],
          str(kopf[:1]))
    umfang = zeilen(out, "CHEM_SPECIES")
    check("2 Spezies, 3 Reaktionen",
          bool(umfang) and umfang[0][1] == "2" and umfang[0][3] == "3",
          str(umfang[:1]))

    ergebnisse = zeilen(out, "RESULT")
    check("beide Betriebspunkte gerechnet", len(ergebnisse) == 2,
          f"{len(ergebnisse)} Ergebniszeilen")

    dichten = zeilen(out, "SPECIES_DENSITY")
    check("Dichte je Spezies ausgegeben", len(dichten) == 2 * 2,
          f"{len(dichten)} Zeilen")
    namen = {d[1] for d in dichten}
    check("Xe und Xe+ benannt", namen == {"Xe", "Xe+"}, str(namen))
    check("alle Dichten positiv", all(float(d[2]) > 0 for d in dichten))

    if ergebnisse:
        n, ng, Te, Tg = (float(x) for x in ergebnisse[0][1:5])
        check("Plasmadichte plausibel", 1e15 < n < 1e19, f"{n:.3e}")
        check("Neutraldichte plausibel", 1e17 < ng < 1e21, f"{ng:.3e}")
        check("Elektronentemperatur plausibel", 2.0 < Te < 8.0, f"{Te}")
        check("Gastemperatur plausibel", 280.0 < Tg < 600.0, f"{Tg}")
        # Die zusammengefasste Plasmadichte ist die Summe der positiven Ionen,
        # hier also genau die Xe+-Dichte.
        xe_plus = [float(d[2]) for d in dichten if d[1] == "Xe+"][0]
        check("Plasmadichte gleich Ionendichte", abs(n - xe_plus) <= 1e-3 * n,
              f"{n:.4e} vs {xe_plus:.4e}")

    zusammen = zeilen(out, "SUMMARY")
    check("beide Punkte konvergiert", bool(zusammen) and zusammen[0][1] == "2",
          str(zusammen[:1]))


def test_ohne_paket(prog: Path) -> None:
    out, _ = lauf(prog)
    check("ohne Angabe der alte Weg", "SOLVER_PATH legacy" in out)
    check("ohne Angabe kein Paketkopf", "CHEM_PACKAGE" not in out)
    check("ohne Angabe Ergebnisse da", len(zeilen(out, "RESULT")) == 2)


def test_kaputtes_paket(prog: Path) -> None:
    out, err = lauf(prog, "chemistry_package gibt_es_nicht\n")
    check("fehlerhaftes Paket faellt zurueck", "SOLVER_PATH legacy" in out)
    check("Rueckfall wird gemeldet", "WARNUNG" in err, err[:200])
    check("Lauf trotzdem vollstaendig", len(zeilen(out, "RESULT")) == 2)


def main() -> int:
    prog = kern()
    if prog is None:
        print("chabert fehlt. Mit 'python build.py' uebersetzen.", file=sys.stderr)
        return 1
    try:
        test_mit_paket(prog)
        test_ohne_paket(prog)
        test_kaputtes_paket(prog)
    finally:
        PARAMS.unlink(missing_ok=True)

    print(f"\n{'=' * 60}")
    print(f"  Ergebnis: {passed} bestanden, {failed} fehlgeschlagen")
    if fehlernamen:
        print(f"  Fehlgeschlagen: {', '.join(fehlernamen)}")
    print(f"{'=' * 60}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
