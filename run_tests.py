"""
run_tests.py -- Faehrt den gesamten Testbestand und sagt, was laeuft.

Die vorhandenen Tests sind eigenstaendige Programme, die ihr Ergebnis ueber
den Exit-Code melden: null heisst bestanden. Dieser Starter sucht sie, fuehrt
sie einzeln aus und fasst zusammen. Er schreibt die Tests nicht um.

Neben den Python-Tests laufen die uebersetzten C++-Testprogramme mit, sofern
sie vorliegen; erzeugt werden sie mit ``python build.py --tests``.

Aufruf:
    python run_tests.py                  # alles
    python run_tests.py --list           # nur auflisten, nichts ausfuehren
    python run_tests.py --only chemistry # nur Tests, deren Name das enthaelt
    python run_tests.py --build          # vorher den Rechenkern uebersetzen
    python run_tests.py --verbose        # Ausgabe auch bestandener Tests zeigen
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

import build as cpp_build

PROJECT_DIR = Path(__file__).resolve().parent

#: Vorgabe je Test. Die Laeufe mit dem C++-Kern brauchen deutlich laenger als
#: die reinen Rechentests, deshalb grosszuegig.
DEFAULT_TIMEOUT_S = 600

#: Traegt zwar das Namensmuster, ist aber die gemeinsame Helferschicht der
#: Tests und selbst keiner.
KEINE_TESTS = {"test_helpers"}

BESTANDEN = "bestanden"
FEHLER = "fehlgeschlagen"
ZEITUEBERSCHREITUNG = "Zeit ueberschritten"
#: Der Test ist nicht an der Sache gescheitert, sondern an einem fehlenden
#: Baustein im aufrufenden Interpreter. Das ist eine Aussage ueber die
#: Umgebung, nicht ueber das Modell, und wird getrennt ausgewiesen.
UMGEBUNG = "Umgebung unvollstaendig"


class Ergebnis:
    def __init__(self, name: str, status: str, dauer_s: float, ausgabe: str,
                 fehlendes_paket: str = ""):
        self.name = name
        self.status = status
        self.dauer_s = dauer_s
        self.ausgabe = ausgabe
        self.fehlendes_paket = fehlendes_paket

    @property
    def ok(self) -> bool:
        return self.status == BESTANDEN


def _fehlendes_paket(ausgabe: str) -> str:
    """Name des Pakets, an dem ein Lauf gescheitert ist, sonst leer."""
    marke = "ModuleNotFoundError: No module named "
    pos = ausgabe.rfind(marke)
    if pos < 0:
        return ""
    return ausgabe[pos + len(marke):].split("\n", 1)[0].strip().strip("'\"")


def finde_python_tests(project_dir: Path = PROJECT_DIR) -> list[Path]:
    """Alle eigenstaendigen Python-Tests, alphabetisch."""
    return sorted(p for p in project_dir.glob("test_*.py")
                  if p.is_file() and p.stem not in KEINE_TESTS)


def finde_cpp_tests(project_dir: Path = PROJECT_DIR) -> list[Path]:
    """Uebersetzte C++-Testprogramme, sofern vorhanden."""
    gefunden = []
    for name in cpp_build.CPP_TESTS:
        for kandidat in (project_dir / name, project_dir / f"{name}.exe"):
            if kandidat.exists():
                gefunden.append(kandidat)
                break
    return gefunden


def fuehre_aus(pfad: Path, timeout_s: int, project_dir: Path = PROJECT_DIR) -> Ergebnis:
    """Fuehrt einen Test aus und bewertet ihn ueber seinen Exit-Code."""
    befehl = [sys.executable, str(pfad)] if pfad.suffix == ".py" else [str(pfad)]
    start = time.monotonic()
    try:
        r = subprocess.run(befehl, cwd=project_dir, capture_output=True,
                           text=True, timeout=timeout_s)
    except subprocess.TimeoutExpired:
        return Ergebnis(pfad.stem, ZEITUEBERSCHREITUNG, time.monotonic() - start,
                        f"Nach {timeout_s} s abgebrochen.")
    dauer = time.monotonic() - start
    ausgabe = (r.stdout or "") + (r.stderr or "")
    if r.returncode == 0:
        return Ergebnis(pfad.stem, BESTANDEN, dauer, ausgabe)
    paket = _fehlendes_paket(ausgabe)
    if paket:
        return Ergebnis(pfad.stem, UMGEBUNG, dauer, ausgabe, paket)
    return Ergebnis(pfad.stem, FEHLER, dauer, ausgabe)


def _zeile(erg: Ergebnis) -> str:
    zeichen = {BESTANDEN: "ok  ", FEHLER: "FEHL",
               ZEITUEBERSCHREITUNG: "ZEIT", UMGEBUNG: "UMGB"}[erg.status]
    hinweis = f"  ({erg.fehlendes_paket} fehlt)" if erg.fehlendes_paket else ""
    return f"  {zeichen}  {erg.name:<34} {erg.dauer_s:6.1f} s{hinweis}"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Faehrt den Testbestand.")
    ap.add_argument("--list", action="store_true", help="nur auflisten")
    ap.add_argument("--only", metavar="TEXT",
                    help="nur Tests, deren Name TEXT enthaelt")
    ap.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT_S,
                    metavar="S", help=f"Zeitgrenze je Test in Sekunden (Vorgabe {DEFAULT_TIMEOUT_S})")
    ap.add_argument("--build", action="store_true",
                    help="vor dem Lauf den Rechenkern samt Testprogrammen uebersetzen")
    ap.add_argument("--verbose", action="store_true",
                    help="Ausgabe auch bestandener Tests zeigen")
    args = ap.parse_args(argv)

    if args.build:
        print("Uebersetze den Rechenkern ...")
        rc = cpp_build.build(with_tests=True, verbose=False)
        if rc != 0:
            print("Bau fehlgeschlagen, Testlauf abgebrochen.", file=sys.stderr)
            return rc
        print()

    tests = finde_python_tests() + finde_cpp_tests()
    if args.only:
        tests = [t for t in tests if args.only.lower() in t.stem.lower()]
    if not tests:
        print("Keine Tests gefunden.", file=sys.stderr)
        return 1

    if args.list:
        for t in tests:
            print(f"  {t.stem}")
        print(f"\n{len(tests)} Tests.")
        return 0

    kern = PROJECT_DIR / cpp_build.binary_name()
    if not kern.exists():
        print(f"Hinweis: {kern.name} fehlt -- Tests gegen den C++-Kern werden "
              f"scheitern. Mit --build uebersetzen.\n")

    print(f"Interpreter: {sys.executable} ({sys.version.split()[0]})")
    print(f"{len(tests)} Tests:\n")
    ergebnisse = []
    for t in tests:
        erg = fuehre_aus(t, args.timeout)
        ergebnisse.append(erg)
        print(_zeile(erg), flush=True)
        # Bei fehlendem Baustein sagt die Ausgabe nichts ueber das Modell --
        # der Hinweis in der Zeile genuegt.
        zeigen = args.verbose or erg.status in (FEHLER, ZEITUEBERSCHREITUNG)
        if erg.ausgabe.strip() and zeigen:
            for zeile in erg.ausgabe.rstrip().splitlines():
                print(f"        {zeile}")

    bestanden = [e for e in ergebnisse if e.ok]
    echte_fehler = [e for e in ergebnisse if e.status in (FEHLER, ZEITUEBERSCHREITUNG)]
    umgebung = [e for e in ergebnisse if e.status == UMGEBUNG]
    gesamt_s = sum(e.dauer_s for e in ergebnisse)

    print(f"\n{len(bestanden)} von {len(ergebnisse)} bestanden, {gesamt_s:.0f} s gesamt.")
    if echte_fehler:
        print("Nicht bestanden: " + ", ".join(e.name for e in echte_fehler))
    if umgebung:
        fehlt = sorted({e.fehlendes_paket for e in umgebung})
        print(f"Nicht gelaufen, weil im Interpreter fehlt: {', '.join(fehlt)} "
              f"({', '.join(e.name for e in umgebung)})")
    return 0 if not (echte_fehler or umgebung) else 1


if __name__ == "__main__":
    sys.exit(main())
