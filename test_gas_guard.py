#!/usr/bin/env python3
"""
test_gas_guard.py -- Keine Xenon-Raten für andere Gase.

Die eingebauten Ratenanpassungen sind Xenon-Polynome. Für Krypton und Argon
liegen die Stoffkonstanten vor, aber keine Querschnittsdaten; bis 2026-09-02
fiel der Rechenkern still auf die Xenon-Anpassungen zurück und lieferte
Zahlen, die wie Argon aussahen und Xenon waren.

Der Lauf wird jetzt abgebrochen und begründet. Wer es trotzdem will, sagt es
ausdrücklich über `allow_foreign_rate_fits`.

Aufruf: python test_gas_guard.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PARAMS = SCRIPT_DIR / "gas_guard_params.txt"

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
jjmax 1
solve_mode 2
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


def lauf(prog: Path, zusatz: str) -> tuple[int, str, str]:
    PARAMS.write_text(BASIS + zusatz, encoding="utf-8")
    r = subprocess.run([str(prog), PARAMS.name], cwd=SCRIPT_DIR,
                       capture_output=True, text=True, timeout=600)
    return r.returncode, r.stdout, r.stderr


def test_argon_wird_abgelehnt(prog: Path) -> None:
    rc, out, err = lauf(prog, "gas_species argon\nrate_model 0\n")
    check("Argon mit eingebauten Anpassungen bricht ab", rc != 0, f"Exit {rc}")
    check("Grund wird genannt", "gelten nur fuer Xenon" in err, err[:200])
    check("Ionisation benannt", "Ionisation" in err)
    check("Anregung benannt", "Anregung" in err)
    check("elastischer Stoss benannt", "Elastischer Stoss" in err)
    check("Abhilfe wird genannt",
          "allow_foreign_rate_fits" in err and "chemistry_package" in err, err[-200:])
    check("kein Ergebnis ausgegeben", "RESULT " not in out)


def test_argon_mit_ausdruecklicher_erlaubnis(prog: Path) -> None:
    rc, out, err = lauf(prog, "gas_species argon\nrate_model 0\n"
                              "allow_foreign_rate_fits 1\n")
    check("mit Erlaubnis laeuft es", rc == 0, f"Exit {rc}")
    check("Warnung erscheint", "WARNUNG" in out, out[:200])
    check("Lauf wird als solcher gekennzeichnet",
          "FOREIGN_RATE_FITS argon" in out)


def test_eigener_elastischer_wert(prog: Path) -> None:
    """Ein ausdruecklich gesetzter Wert ist eine Entscheidung, keine Uebernahme."""
    rc, out, err = lauf(prog, "gas_species argon\nrate_model 0\nkel_constant 5e-14\n")
    check("Ionisation blockiert weiterhin", rc != 0, f"Exit {rc}")
    check("elastischer Stoss nicht mehr beanstandet",
          "Elastischer Stoss" not in err, err[:300])


def test_xenon_unbehelligt(prog: Path) -> None:
    for zusatz, name in (("gas_species xenon\nrate_model 0\n", "Legacy"),
                         ("gas_species xenon\nrate_model 2\n", "tabelliert")):
        rc, out, err = lauf(prog, zusatz)
        check(f"Xenon {name} laeuft", rc == 0, f"Exit {rc} -- {err[:150]}")
        check(f"Xenon {name} liefert ein Ergebnis", "RESULT " in out)
        check(f"Xenon {name} ohne Beanstandung", "gelten nur fuer Xenon" not in err)


def test_paket_ist_ausgenommen(prog: Path) -> None:
    """Ein Chemiepaket bringt eigene Raten mit und ist nicht betroffen."""
    rc, out, err = lauf(prog, "gas_species xenon\nrate_model 0\n"
                              "chemistry_package xenon_biagi\n")
    check("Lauf ueber ein Paket laeuft", rc == 0, f"Exit {rc}")
    check("Paketweg genommen", "SOLVER_PATH generic" in out)


def main() -> int:
    prog = kern()
    if prog is None:
        print("chabert fehlt. Mit 'python build.py' uebersetzen.", file=sys.stderr)
        return 1
    try:
        test_argon_wird_abgelehnt(prog)
        test_argon_mit_ausdruecklicher_erlaubnis(prog)
        test_eigener_elastischer_wert(prog)
        test_xenon_unbehelligt(prog)
        test_paket_ist_ausgenommen(prog)
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
