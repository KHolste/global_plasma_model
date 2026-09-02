"""
build.py -- Uebersetzung des C++-Rechenkerns. Eine Definition fuer alle Aufrufer.

Der Kern besteht aus sieben Uebersetzungseinheiten plus dem Einstiegspunkt.
Welche das sind und mit welchen Schaltern sie uebersetzt werden, steht
ausschliesslich hier. Die Oberflaeche baut sich daraus ihre Prozesszeile,
die Kommandozeile ruft die Schritte direkt auf.

Aufruf von der Kommandozeile:
    python build.py            # uebersetzt den Rechenkern
    python build.py --clean    # entfernt Objektdateien und Programm zuvor
    python build.py --tests    # uebersetzt zusaetzlich die C++-Testprogramme
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent

#: Uebersetzungseinheiten des Rechenkerns, ohne Endung. Reihenfolge ist die
#: Uebersetzungsreihenfolge; sie spielt fuer das Ergebnis keine Rolle, macht
#: aber die Ausgabe nachvollziehbar.
CPP_MODULES = ["bessel_wrapper", "sim_config", "rates", "physics", "solver",
               "sim_logging", "chem_loader"]

#: Einstiegspunkt, wird gegen die Objektdateien gebunden.
CPP_SOURCE = "main.cpp"

#: Uebersetzungsschalter. -march=native bindet das Ergebnis an die bauende
#: Maschine; das ist gewollt, weil das Programm hier laeuft, wo es entsteht.
CXXFLAGS = ["-O3", "-march=native", "-std=c++17"]

#: Eigenstaendige C++-Testprogramme. Jedes bindet gegen dieselben Objekte.
CPP_TESTS = ["test_chem_system", "test_generic_lm", "test_chem_loader",
             "test_extraction"]


def object_suffix(platform: str = sys.platform) -> str:
    """Endung der Objektdateien. Windows-Werkzeuge erwarten .obj."""
    return ".obj" if platform == "win32" else ".o"


def binary_name(platform: str = sys.platform) -> str:
    """Name des erzeugten Programms."""
    return "chabert.exe" if platform == "win32" else "chabert"


def compile_steps(cc: str, obj_ext: str, binary: str,
                  with_tests: bool = False) -> list[list[str]]:
    """Alle Aufrufe, die den Rechenkern erzeugen, als Argumentlisten.

    Jeder Eintrag ist ein vollstaendiger Compileraufruf. Aufrufer koennen sie
    direkt ausfuehren oder zu einer Shell-Zeile verketten -- die Schritte
    selbst sind an beiden Stellen dieselben.
    """
    objs = [f"{m}{obj_ext}" for m in CPP_MODULES]
    steps = [[cc, *CXXFLAGS, "-c", f"{m}.cpp", "-o", obj]
             for m, obj in zip(CPP_MODULES, objs)]
    steps.append([cc, *CXXFLAGS, CPP_SOURCE, *objs, "-o", binary])
    if with_tests:
        for t in CPP_TESTS:
            steps.append([cc, *CXXFLAGS, f"{t}.cpp", *objs, "-o", t])
    return steps


def build_command(cwd: str, cc: str, use_wsl: bool,
                  platform: str = sys.platform) -> tuple[str, list[str]]:
    """Prozess und Argumente fuer einen Bau in einem Rutsch.

    Fuer Aufrufer, die den Bau nicht selbst Schritt fuer Schritt fahren, etwa
    die Oberflaeche mit ihrem eigenen Prozessobjekt. ``cwd`` ist das
    Arbeitsverzeichnis in der Schreibweise, die der aufgerufenen Shell passt --
    bei WSL also der uebersetzte Pfad.
    """
    if use_wsl:
        steps = compile_steps("g++", ".o", "chabert")
        chain = " && ".join(" ".join(s) for s in steps)
        return "wsl", ["bash", "-c", f'cd "{cwd}" && {chain} 2>&1']
    shell = "cmd" if platform == "win32" else "bash"
    flag = "/c" if platform == "win32" else "-c"
    steps = compile_steps(f'"{cc}"', object_suffix(platform), binary_name(platform))
    chain = " && ".join(" ".join(s) for s in steps)
    return shell, [flag, chain]


def clean(project_dir: Path = PROJECT_DIR, platform: str = sys.platform) -> list[str]:
    """Entfernt Objektdateien, Programm und Testprogramme. Gibt zurueck, was weg ist."""
    removed = []
    targets = [f"{m}{object_suffix(platform)}" for m in CPP_MODULES]
    targets.append(binary_name(platform))
    targets.extend(CPP_TESTS)
    targets.extend(f"{t}.exe" for t in CPP_TESTS)
    for name in targets:
        p = project_dir / name
        if p.exists():
            p.unlink()
            removed.append(name)
    return removed


def build(project_dir: Path = PROJECT_DIR, with_tests: bool = False,
          verbose: bool = True) -> int:
    """Uebersetzt den Rechenkern Schritt fuer Schritt. Rueckgabe ist der Exit-Code."""
    cc = shutil.which("g++")
    if not cc:
        print("g++ nicht gefunden. Bitte einen C++-Uebersetzer bereitstellen.",
              file=sys.stderr)
        return 127

    steps = compile_steps(cc, object_suffix(), binary_name(), with_tests)
    for i, step in enumerate(steps, 1):
        if verbose:
            ziel = step[-1]
            print(f"[{i}/{len(steps)}] {ziel}")
        r = subprocess.run(step, cwd=project_dir, capture_output=True, text=True)
        if r.stdout.strip():
            print(r.stdout.rstrip())
        if r.stderr.strip():
            print(r.stderr.rstrip(), file=sys.stderr)
        if r.returncode != 0:
            print(f"Abbruch bei Schritt {i}: {' '.join(step)}", file=sys.stderr)
            return r.returncode
    if verbose:
        print(f"Fertig: {project_dir / binary_name()}")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Uebersetzt den C++-Rechenkern.")
    ap.add_argument("--clean", action="store_true",
                    help="Objektdateien und Programm zuvor entfernen")
    ap.add_argument("--tests", action="store_true",
                    help="zusaetzlich die C++-Testprogramme uebersetzen")
    ap.add_argument("--quiet", action="store_true", help="nur Fehler ausgeben")
    args = ap.parse_args(argv)

    if args.clean:
        removed = clean()
        if not args.quiet:
            print("Entfernt: " + (", ".join(removed) if removed else "nichts"))
    return build(with_tests=args.tests, verbose=not args.quiet)


if __name__ == "__main__":
    sys.exit(main())
