#!/usr/bin/env python3
"""
test_power_range.py -- Leistungsbereich und Treibstoffgrenze.

Bis 2026-09-02 brach die Rechnung im selbstkonsistenten Betrieb oberhalb etwa
sechzig bis siebzig Watt ab. Das war kein physikalischer Befund: Die Lösung
existiert, sie liegt nur zu weit von jedem kalten Startwert entfernt. Der
Löser fährt die Leistung deshalb als Rückfall von unten hoch. Dieser Test
hält den wiedergewonnenen Bereich fest.

Zusätzlich zwei Aussagen, die immer gelten müssen:
  - Der Strahlstrom kann nicht größer sein als der zugeführte Teilchenstrom
    mal Elementarladung. Mehr als jedes Atom kann man nicht extrahieren.
  - Ein unerreichbarer Zielstrom ist eine physikalische Aussage und muss als
    solche gemeldet werden, nicht als Rechenfehler.

Aufruf: python test_power_range.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PARAMS = SCRIPT_DIR / "power_range_params.txt"

E_CH = 1.602176487e-19
SCCM_TO_PPS = 4.477962e17
Q0_SCCM = 0.40
Q0_SCHRITT = 0.02

BASIS = f"""R 0.02
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
P_RFG_max 500.0
Q0sccm_start {Q0_SCCM}
Q0sccm_step {Q0_SCHRITT}
jjmax 3
gas_species xenon
rate_model 2
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


def lauf(prog: Path, zusatz: str) -> str:
    PARAMS.write_text(BASIS + zusatz, encoding="utf-8")
    r = subprocess.run([str(prog), PARAMS.name], cwd=SCRIPT_DIR,
                       capture_output=True, text=True, timeout=900)
    return r.stdout


def zeilen(text: str, marke: str) -> list[list[str]]:
    return [z.split() for z in text.splitlines() if z.startswith(marke + " ")]


def grenzstrom_mA(punkt: int) -> float:
    """Hoechster moeglicher Strahlstrom: jedes zugefuehrte Atom einfach geladen.

    Der Durchfluss waechst mit jedem Punkt des Durchlaufs, also auch die Grenze.
    """
    q0 = Q0_SCCM + punkt * Q0_SCHRITT
    return E_CH * q0 * SCCM_TO_PPS * 1000.0


def test_leistungsbereich(prog: Path) -> None:
    """Selbstkonsistenter Betrieb ueber den ganzen sinnvollen Bereich."""
    print(f"  Treibstoffgrenze: {grenzstrom_mA(0):.2f} mA bei {Q0_SCCM} sccm, "
          f"{grenzstrom_mA(2):.2f} mA bei {Q0_SCCM + 2*Q0_SCHRITT:.2f} sccm")

    for P in (18.0, 40.0, 60.0, 80.0, 120.0):
        out = lauf(prog, f"solve_mode 2\nP_RFG {P}\n")
        zus = zeilen(out, "SUMMARY")
        ok = int(zus[0][1]) if zus else -1
        check(f"{P:.0f} W: alle drei Punkte konvergiert", ok == 3,
              f"{ok} von 3 -- {zus[:1]}")

        for i, e in enumerate(zeilen(out, "RESULT")):
            I_mA = float(e[5])
            grenze = grenzstrom_mA(i)
            check(f"{P:.0f} W, Punkt {i+1}: Strahlstrom unter der Treibstoffgrenze",
                  I_mA <= grenze * 1.001,
                  f"{I_mA:.2f} mA gegen {grenze:.2f} mA")


def test_unerreichbarer_strom(prog: Path) -> None:
    """Ein Zielstrom oberhalb der Treibstoffgrenze ist Physik, kein Rechenfehler."""
    ziel = grenzstrom_mA(0) * 2.0
    out = lauf(prog, f"solve_mode 1\nP_RFG 18.0\nI_soll {ziel:.1f}\n")

    ergebnisse = zeilen(out, "IFIX_RESULT")
    check("drei Punkte bewertet", len(ergebnisse) == 3, str(len(ergebnisse)))
    check("als unerreichbar gemeldet",
          all(e[6] == "above_P_max" for e in ergebnisse),
          str([e[6] for e in ergebnisse]))
    check("nicht als Rechenfehler gemeldet",
          all(e[6] != "numerical_fail" for e in ergebnisse))

    zus = zeilen(out, "SUMMARY_DETAIL")
    if zus:
        text = " ".join(zus[0])
        check("als fehlende physikalische Loesung gezaehlt",
              "no_physical_solution=3" in text, text)

    # Der beste gefundene Strom muss unterhalb der Treibstoffgrenze liegen
    for i, e in enumerate(ergebnisse):
        I_best = float(e[4])
        grenze = grenzstrom_mA(i)
        check(f"Punkt {i+1}: bester Strom unter der Treibstoffgrenze",
              I_best <= grenze * 1.001,
              f"{I_best:.2f} mA gegen {grenze:.2f} mA")


def test_erreichbarer_strom(prog: Path) -> None:
    """Ein Zielstrom innerhalb der Grenze muss getroffen werden."""
    out = lauf(prog, "solve_mode 1\nP_RFG 18.0\nI_soll 15.0\n")
    ergebnisse = zeilen(out, "IFIX_RESULT")
    check("Zielstrom erreicht", all(e[6] == "converged" for e in ergebnisse),
          str([e[6] for e in ergebnisse]))
    for e in ergebnisse:
        abweichung = abs(float(e[5]))
        check("Strom trifft die Vorgabe", abweichung < 0.1, f"{abweichung:.3f} mA")


def main() -> int:
    prog = kern()
    if prog is None:
        print("chabert fehlt. Mit 'python build.py' uebersetzen.", file=sys.stderr)
        return 1
    try:
        test_leistungsbereich(prog)
        test_unerreichbarer_strom(prog)
        test_erreichbarer_strom(prog)
    finally:
        PARAMS.unlink(missing_ok=True)

    print(f"\n{'=' * 60}")
    print(f"  Ergebnis: {passed} bestanden, {failed} fehlgeschlagen")
    if fehlernamen:
        print(f"  Fehlgeschlagen: {', '.join(fehlernamen[:8])}")
    print(f"{'=' * 60}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
