#!/usr/bin/env python3
"""
test_solver.py – Minimale Testsuite fuer den Global Xenon Model Solver.

Fuehrt den kompilierten Solver als Blackbox aus und prueft:
  1. Referenzlauf mit Standardparametern
  2. Konvergenzverhalten
  3. Saubere Behandlung unphysikalischer Parameter
  4. Regression gegen gespeicherte Referenzwerte

Ausfuehrung:
    python test_solver.py

Voraussetzung:
    Kompiliertes Binary 'chabert' (Linux/WSL) oder 'chabert.exe' (Windows)
    im selben Verzeichnis wie dieses Skript.
"""
from __future__ import annotations

import os
import sys
import subprocess
import re
import math
from pathlib import Path
from dataclasses import dataclass, field

SCRIPT_DIR = Path(__file__).resolve().parent
os.chdir(SCRIPT_DIR)

# ─── Hilfsfunktionen ─────────────────────────────────────────

def find_binary() -> str:
    """Finde das Solver-Binary (nativ oder WSL)."""
    if sys.platform == "win32":
        # Versuche WSL
        wsl = subprocess.run(["wsl", "bash", "-c", f'test -f "{wsl_path()}/chabert"'],
                             capture_output=True)
        if wsl.returncode == 0:
            return "wsl"
        exe = SCRIPT_DIR / "chabert.exe"
        if exe.exists():
            return str(exe)
    else:
        native = SCRIPT_DIR / "chabert"
        if native.exists():
            return str(native)
    return ""


def wsl_path() -> str:
    """Konvertiere Windows-Pfad zu WSL-Pfad."""
    p = str(SCRIPT_DIR).replace("\\", "/")
    if len(p) >= 2 and p[1] == ":":
        return f"/mnt/{p[0].lower()}{p[2:]}"
    return p


def run_solver(params: dict[str, float | int | str]) -> tuple[int, str]:
    """Schreibe params.txt, starte Solver, gib (returncode, stdout) zurueck."""
    config_path = SCRIPT_DIR / "test_params.txt"
    with open(config_path, "w") as f:
        f.write("# Testparameter (automatisch generiert)\n")
        for k, v in params.items():
            f.write(f"{k} {v}\n")

    binary = find_binary()
    if not binary:
        return -1, "ERROR: Solver binary not found"

    if binary == "wsl":
        cmd = ["wsl", "bash", "-c",
               f'cd "{wsl_path()}" && ./chabert test_params.txt 2>&1']
    else:
        cmd = [binary, str(config_path)]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        return result.returncode, result.stdout
    except subprocess.TimeoutExpired:
        return -2, "ERROR: Solver timed out (120s)"
    finally:
        config_path.unlink(missing_ok=True)


def parse_results(stdout: str) -> list[dict]:
    """Parse RESULT-Zeilen aus stdout."""
    results = []
    for line in stdout.splitlines():
        parts = line.strip().split()
        if len(parts) >= 7 and parts[0] == "RESULT":
            try:
                results.append({
                    "n":     float(parts[1]),
                    "ng":    float(parts[2]),
                    "Te":    float(parts[3]),
                    "Tg":    float(parts[4]),
                    "I_mA":  float(parts[5]),
                    "P_RFG": float(parts[6]),
                })
            except (ValueError, IndexError):
                pass
    return results


def parse_summary(stdout: str) -> dict:
    """Parse SUMMARY_DETAIL-Zeile."""
    for line in stdout.splitlines():
        if line.startswith("SUMMARY_DETAIL"):
            d = {}
            for token in line.split():
                if "=" in token:
                    k, v = token.split("=", 1)
                    try:
                        d[k] = int(v)
                    except ValueError:
                        d[k] = v
            return d
    return {}


def count_tag(stdout: str, tag: str) -> int:
    """Zaehle Zeilen die mit einem bestimmten Tag beginnen."""
    return sum(1 for line in stdout.splitlines() if line.strip().startswith(tag))


# ─── Test-Infrastruktur ──────────────────────────────────────

@dataclass
class TestResult:
    name: str
    passed: bool
    message: str = ""

test_results: list[TestResult] = []


def report(name: str, passed: bool, msg: str = ""):
    """Registriere ein Testergebnis."""
    test_results.append(TestResult(name, passed, msg))
    status = "\033[32mPASS\033[0m" if passed else "\033[31mFAIL\033[0m"
    print(f"  [{status}] {name}" + (f" — {msg}" if msg else ""))


def assert_in_range(name: str, value: float, lo: float, hi: float, unit: str = ""):
    """Pruefe ob ein Wert im erwarteten Bereich liegt."""
    ok = lo <= value <= hi
    msg = f"{value:.4g} {'in' if ok else 'NOT in'} [{lo}, {hi}]"
    if unit:
        msg += f" {unit}"
    report(name, ok, msg)
    return ok


# ─── Standard-Testparameter ──────────────────────────────────

DEFAULT_PARAMS = {
    "R": 0.02,
    "L": 0.04,
    "betai": 0.5,
    "betag": 0.05145,
    "frequency": 2.5e6,
    "Nw": 6.0,
    "R_ohm": 0.36,
    "Rc": 0.02,
    "lc": 0.04,
    "Vgrid": 1500.0,
    "sgrid": 0.001,
    "P_RFG": 18.0,
    "P_RFG_max": 80.0,
    "Q0sccm": 0.475,
    "Q0sccm_start": 0.30,
    "Q0sccm_step": 0.05,
    "jjmax": 5,
    "I_soll": 15.0,
    "solve_mode": 1,
}


# ─── Test 1: Referenzlauf ─────────────────────────────────────

def test_reference_run():
    """Solver laeuft durch und liefert plausible Ergebnisse."""
    print("\n--- Test 1: Referenzlauf ---")

    rc, stdout = run_solver(DEFAULT_PARAMS)
    report("Solver startet und beendet", rc == 0,
           f"returncode={rc}" if rc != 0 else "")

    if rc != 0:
        report("Referenzlauf abgebrochen", False, stdout[:200])
        return

    summary = parse_summary(stdout)
    total = summary.get("total", 0)
    report("SUMMARY vorhanden", total > 0, f"total={total}")

    results = parse_results(stdout)
    report("Mindestens 1 RESULT", len(results) >= 1, f"{len(results)} Punkte")

    if results:
        r = results[0]
        assert_in_range("Te plausibel", r["Te"], 0.5, 15.0, "eV")
        assert_in_range("Tg plausibel", r["Tg"], 250, 2000, "K")
        assert_in_range("n plausibel", r["n"], 1e14, 1e20, "m^-3")
        assert_in_range("ng plausibel", r["ng"], 1e16, 1e22, "m^-3")
        assert_in_range("I_mA plausibel", r["I_mA"], 0.1, 100, "mA")
        assert_in_range("P_RFG plausibel", r["P_RFG"], 1, 5000, "W")


# ─── Test 2: Konvergenz ──────────────────────────────────────

def test_convergence():
    """Pruefe dass bei Standardparametern Konvergenz erreicht wird."""
    print("\n--- Test 2: Konvergenz ---")

    params = {**DEFAULT_PARAMS, "jjmax": 3, "Q0sccm_start": 0.40}
    rc, stdout = run_solver(params)

    if rc != 0:
        report("Solver laeuft", False, f"returncode={rc}")
        return

    n_converged = count_tag(stdout, "CONVERGED")
    report("Mindestens 1 Punkt konvergiert", n_converged >= 1,
           f"{n_converged} CONVERGED-Tags")

    summary = parse_summary(stdout)
    conv = summary.get("converged", 0)
    report("Summary zaehlt konvergierte Punkte", conv >= 1,
           f"converged={conv}")


# ─── Test 3: Unphysikalische Parameter ───────────────────────

def test_unphysical_params():
    """Solver stuerzt nicht ab bei unphysikalischen Parametern."""
    print("\n--- Test 3: Unphysikalische Parameter ---")

    # Extrem niedriger Massendurchsatz → kein physikalisches Gleichgewicht erwartet
    params = {**DEFAULT_PARAMS,
              "Q0sccm_start": 0.001,
              "Q0sccm_step": 0.001,
              "jjmax": 2,
              "I_soll": 100.0}  # Unerreichbar hoher Zielstrom

    rc, stdout = run_solver(params)
    report("Solver stuerzt nicht ab", rc == 0, f"returncode={rc}")

    n_nophys = count_tag(stdout, "NO_PHYSICAL_SOLUTION")
    n_numfail = count_tag(stdout, "NUMERICAL_FAIL")
    report("Saubere Fehlermeldung statt Crash",
           n_nophys + n_numfail >= 1,
           f"NO_PHYSICAL={n_nophys}, NUMERICAL_FAIL={n_numfail}")

    summary = parse_summary(stdout)
    report("SUMMARY trotzdem vorhanden",
           "total" in summary, str(summary))


# ─── Test 4: Regression ──────────────────────────────────────

# Referenzwerte fuer den Regressions-Check (bei Aenderung bewusst aktualisieren).
# Format: (Q0sccm_start, erwartetes I_mA, erwartetes P_RFG, Toleranz_relativ)
# Die Werte muessen beim ersten erfolgreichen Lauf hier eingetragen werden.
REGRESSION_REF = None  # wird dynamisch beim ersten Lauf gesetzt


def test_regression():
    """Pruefe numerische Stabilitaet gegen Referenzwerte."""
    print("\n--- Test 4: Regressions-Check ---")

    params = {**DEFAULT_PARAMS, "Q0sccm_start": 0.40, "Q0sccm_step": 0.1, "jjmax": 1}
    rc, stdout = run_solver(params)

    if rc != 0:
        report("Solver laeuft", False)
        return

    results = parse_results(stdout)
    if not results:
        report("RESULT vorhanden", False, "Keine RESULT-Zeilen")
        return

    r = results[0]
    ref_file = SCRIPT_DIR / "test_reference.txt"

    if ref_file.exists():
        # Vergleiche mit gespeicherten Referenzwerten
        ref = {}
        for line in ref_file.read_text().splitlines():
            if "=" in line and not line.startswith("#"):
                k, v = line.split("=", 1)
                ref[k.strip()] = float(v.strip())

        tol = 0.05  # 5% relative Toleranz
        for key in ["Te", "I_mA", "P_RFG", "n"]:
            if key in ref and key in r:
                actual = r[key]
                expected = ref[key]
                if expected == 0:
                    continue
                rel_err = abs(actual - expected) / abs(expected)
                ok = rel_err < tol
                report(f"Regression {key}",
                       ok,
                       f"erwartet={expected:.4g} ist={actual:.4g} rel_err={rel_err:.2%}")
    else:
        # Erster Lauf: Referenzwerte speichern
        with open(ref_file, "w") as f:
            f.write("# Referenzwerte fuer Regressions-Test\n")
            f.write(f"# Erzeugt von test_solver.py\n")
            f.write(f"# Parameter: Q0sccm_start=0.40, jjmax=1, solve_mode=1\n")
            for key in ["Te", "Tg", "n", "ng", "I_mA", "P_RFG"]:
                f.write(f"{key} = {r[key]}\n")
        report("Referenz erstellt", True,
               f"Referenzwerte in {ref_file.name} gespeichert (naechster Lauf vergleicht)")


# ─── Test 5: Selbstkonsistenter Modus ────────────────────────

def test_self_consistent_mode():
    """Pruefe solve_mode=2 (P_RFG fest, I ergibt sich)."""
    print("\n--- Test 5: Selbstkonsistenter Modus ---")

    params = {**DEFAULT_PARAMS,
              "solve_mode": 2,
              "P_RFG": 50.0,
              "jjmax": 3}

    rc, stdout = run_solver(params)
    report("Solver startet im Modus 2", rc == 0, f"returncode={rc}")

    if rc != 0:
        return

    results = parse_results(stdout)
    report("Ergebnisse vorhanden", len(results) >= 1, f"{len(results)} Punkte")

    if results:
        # I_mA muss sich ergeben (nicht I_soll sein)
        i_values = [r["I_mA"] for r in results]
        all_same_as_target = all(abs(i - 15.0) < 0.1 for i in i_values)
        report("I_mA ist NICHT fest auf I_soll",
               not all_same_as_target or len(results) < 2,
               f"I_mA-Werte: {[f'{i:.2f}' for i in i_values]}")


# ─── Hauptprogramm ───────────────────────────────────────────

def main() -> int:
    print("=" * 60)
    print("  Global Xenon Model — Solver Test Suite")
    print("=" * 60)

    binary = find_binary()
    if not binary:
        print("\nFEHLER: Solver-Binary nicht gefunden.")
        print("Bitte zuerst kompilieren:")
        print("  g++ -O2 -std=c++17 chabert_modified.cpp -o chabert")
        return 1
    print(f"\nBinary: {binary}")

    test_reference_run()
    test_convergence()
    test_unphysical_params()
    test_regression()
    test_self_consistent_mode()

    # Zusammenfassung
    n_pass = sum(1 for t in test_results if t.passed)
    n_fail = sum(1 for t in test_results if not t.passed)
    n_total = len(test_results)

    print("\n" + "=" * 60)
    print(f"  Ergebnis: {n_pass}/{n_total} bestanden, {n_fail} fehlgeschlagen")
    print("=" * 60)

    if n_fail > 0:
        print("\nFehlgeschlagene Tests:")
        for t in test_results:
            if not t.passed:
                print(f"  - {t.name}: {t.message}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
