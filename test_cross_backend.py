#!/usr/bin/env python3
"""
test_cross_backend.py -- Cross-Backend-Regression: C++ vs Python fuer Xenon.

Fuehrt identische Konfigurationen auf beiden Backends aus und vergleicht
die Ergebnisse toleranzbasiert. Prueft Drift, nicht Identitaet.

Systematische Unterschiede (dokumentiert, kein Testfail):
  - Abbruchschranke: C++ 1e-4 auf dem groessten skalierten Residuum, Python
    eine anders normierte RMS-Groesse mit sehr viel lockererer Schranke
  - P_abs: C++ berechnet via Bessel-BVP (P_abs < P_RFG), Python nimmt P direkt
  - I_beam: C++ nutzt einfachen Bohm-Fluss, Python nutzt Bohm+CL+eta_opt
  - ne: Koppelt an P_abs, daher systematisch verschoben

Fair vergleichbar:
  - Te: Beide loesen dieselbe Energiebilanz (Chabert Eq.13)
  - Tg: Beide loesen dieselbe Thermalbilanz (Chabert Eq.11)
  - Konvergenzstatus: Sollte uebereinstimmen
  - Qualitatives Verhalten: Gleiche Trends ueber Q0-Sweep

Ausfuehrung:
    python test_cross_backend.py                    # alle Faelle
    python test_cross_backend.py chabert_sc_18      # einzelner Fall

Voraussetzung:
    Kompiliertes C++ Binary 'chabert' (WSL oder nativ) im Projektverzeichnis.
"""
from __future__ import annotations
import sys, os, subprocess, json, shutil
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
os.chdir(SCRIPT_DIR)
PYTHON = sys.executable

passed = 0
failed = 0
skipped = 0
errors = []


def check(name, condition, detail=""):
    global passed, failed
    if condition:
        passed += 1
        print(f"    PASS: {name}")
    else:
        failed += 1
        print(f"    FAIL: {name} -- {detail}")
        errors.append(name)


def note(name, detail):
    """Dokumentierter Unterschied, kein Fail."""
    print(f"    NOTE: {name} -- {detail}")


# ═════════════════════════════════════════════════════════════
# C++-Binary finden
# ═════════════════════════════════════════════════════════════
def _wsl_path():
    """Konvertiere Windows-Pfad zu WSL-Pfad."""
    p = str(SCRIPT_DIR).replace("\\", "/")
    if len(p) >= 2 and p[1] == ":":
        return f"/mnt/{p[0].lower()}{p[2:]}"
    return p

def find_cpp_binary():
    # 1. Windows-natives Binary (.exe)
    exe = SCRIPT_DIR / "chabert.exe"
    if exe.exists():
        return "native", str(exe)

    # 2. Linux-natives Binary (direkt aufrufbar, nicht unter Windows)
    if sys.platform != "win32":
        native = SCRIPT_DIR / "chabert"
        if native.exists() and os.access(str(native), os.X_OK):
            return "native", str(native)

    # 3. WSL (Windows): Linux-Binary ueber WSL ausfuehren
    for wsl_cmd in ["wsl", "/c/Windows/System32/wsl.exe"]:
        try:
            wp = _wsl_path()
            r = subprocess.run([wsl_cmd, "bash", "-c", f'test -f "{wp}/chabert"'],
                               capture_output=True, timeout=10)
            if r.returncode == 0:
                return "wsl", wp
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    return None, None


# ═════════════════════════════════════════════════════════════
# Solver ausfuehren (via RunConfig)
# ═════════════════════════════════════════════════════════════
from test_helpers import make_config, write_config, cleanup as _cleanup

def write_run_config(cfg):
    """Schreibe RunConfig als JSON + params.txt (fuer C++ Compat)."""
    write_config(cfg)
    # Zusaetzlich xb_params.txt fuer C++ (Legacy-Compat)
    cfg.to_params_txt(SCRIPT_DIR / "xb_params.txt")


def run_cpp(binary_type, binary_path):
    # C++ liest run_config.json (primaer) oder xb_params.txt (Fallback)
    if binary_type == "wsl":
        for wsl_cmd in ["wsl", "/c/Windows/System32/wsl.exe"]:
            try:
                cmd = [wsl_cmd, "bash", "-c",
                       f'cd "{binary_path}" && ./chabert xb_params.txt 2>/dev/null']
                r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
                if r.returncode == 0 or r.stdout:
                    return r.stdout
            except FileNotFoundError:
                continue
        return ""
    else:
        cmd = [binary_path, str(SCRIPT_DIR / "xb_params.txt")]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        return r.stdout


def run_python(chem_json, mode, param, Q0_start, Q0_step, N):
    # Python liest primaer run_config.json (schon geschrieben)
    cmd = [PYTHON, str(SCRIPT_DIR / "generic_solver.py"),
           str(SCRIPT_DIR / chem_json)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=120,
                       cwd=str(SCRIPT_DIR))
    return r.stdout


def parse_results(stdout):
    """Extrahiere RESULT-Zeilen als Liste von Dicts."""
    results = []
    for line in stdout.strip().split("\n"):
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "RESULT" and len(parts) >= 7:
            results.append({
                "ne": float(parts[1]),
                "ng": float(parts[2]),
                "Te": float(parts[3]),
                "Tg": float(parts[4]),
                "I_mA": float(parts[5]),
                "P": float(parts[6]),
            })
    return results


# ═════════════════════════════════════════════════════════════
# Testfaelle
# ═════════════════════════════════════════════════════════════

# Gemeinsame Chabert-Geometrie (beide Backends identisch)
CHABERT_GEOM = {
    "R": 0.02, "L": 0.04, "betai": 0.5, "betag": 0.05145,
    "frequency": 2.5e6, "Nw": 6, "R_ohm": 0.36, "Rc": 0.02, "lc": 0.04,
    "Vgrid": 1500, "sgrid": 0.001,
    "alpha_e_wall": 7.0, "density_profile_factor": 1.0,
}

CASES = [
    {
        "id": "chabert_sc_18",
        "name": "Chabert SC 18W — Paper-Referenz",
        "why": "Grundlegender SC-Fall bei Standardgeometrie. Ratenmodelle identisch.",
        "geom": {**CHABERT_GEOM, "P_RFG": 18, "Q0sccm": 0.475},
        "cpp_extra": {"solve_mode": 2, "rate_model": 0, "Q0sccm_start": 0.40,
                      "Q0sccm_step": 0.05, "jjmax": 3, "debug_level": 0},
        "py_mode": "sc", "py_param": 18.0,
        "py_sweep": (0.40, 0.05, 3),
        "tolerances": {
            "Te_rel": 0.30,    # 30% — P_abs-Kopplung (C++ via RF-BVP vs Py direkt)
            "Tg_rel": 0.05,    # 5% — schwach gekoppelt
        },
    },
    {
        "id": "chabert_sc_30",
        "name": "Chabert SC 30W — Hoehere Leistung",
        "why": "Prueft Skalierungsverhalten bei erhoehter Leistung.",
        "geom": {**CHABERT_GEOM, "P_RFG": 30, "Q0sccm": 0.475},
        "cpp_extra": {"solve_mode": 2, "rate_model": 0, "Q0sccm_start": 0.40,
                      "Q0sccm_step": 0.05, "jjmax": 3, "debug_level": 0},
        "py_mode": "sc", "py_param": 30.0,
        "py_sweep": (0.40, 0.05, 3),
        "tolerances": {
            "Te_rel": 0.30, "Tg_rel": 0.05,
            "ne_rel": 0.40, "ng_rel": 0.25,
        },
    },
    {
        "id": "chabert_sc_sweep",
        "name": "Chabert SC 18W — Breiterer Q0-Sweep",
        "why": "Prueft Stabiliaet ueber breiteren Massenflussbereich.",
        "geom": {**CHABERT_GEOM, "P_RFG": 18, "Q0sccm": 0.475},
        "cpp_extra": {"solve_mode": 2, "rate_model": 0, "Q0sccm_start": 0.30,
                      "Q0sccm_step": 0.1, "jjmax": 5, "debug_level": 0},
        "py_mode": "sc", "py_param": 18.0,
        "py_sweep": (0.30, 0.1, 5),
        "tolerances": {
            "Te_rel": 0.30, "Tg_rel": 0.05,
            "ne_rel": 0.40, "ng_rel": 0.25,
        },
    },
]

CHEM_XE = "chemistry/xenon_simple/chemistry.json"


# ═════════════════════════════════════════════════════════════
# Vergleichslogik
# ═════════════════════════════════════════════════════════════
def rel_diff(a, b):
    if abs(a) + abs(b) < 1e-30:
        return 0.0
    return abs(a - b) / max(abs(a), abs(b))


def compare_results(case_id, cpp_res, py_res, tol):
    n_points = min(len(cpp_res), len(py_res))
    check(f"{case_id}: same point count",
          len(cpp_res) == len(py_res),
          f"C++={len(cpp_res)}, Py={len(py_res)}")

    if n_points == 0:
        check(f"{case_id}: at least one point", False, "no results from either backend")
        return

    for i in range(n_points):
        c = cpp_res[i]
        p = py_res[i]
        prefix = f"{case_id}[{i}]"

        # ── Fair vergleichbar: Te, Tg ──────────────────────
        # Beide loesen dieselbe Energiebilanz, aber C++ hat niedrigeres
        # P_abs (RF-Kopplung), was Te beeinflusst. 30% Toleranz fuer SC.
        rd = rel_diff(c["Te"], p["Te"])
        check(f"{prefix}: Te within {tol['Te_rel']*100:.0f}%",
              rd <= tol["Te_rel"],
              f"C++={c['Te']:.3f} Py={p['Te']:.3f} diff={rd*100:.1f}%")

        rd = rel_diff(c["Tg"], p["Tg"])
        check(f"{prefix}: Tg within {tol['Tg_rel']*100:.0f}%",
              rd <= tol["Tg_rel"],
              f"C++={c['Tg']:.1f} Py={p['Tg']:.1f} diff={rd*100:.1f}%")

        # ── Strukturelle Offsets: ne, ng, I_beam ───────────
        # C++ bekommt weniger P_abs (Bessel-BVP mit Ohm-Verlusten),
        # Python bekommt P direkt. Daher ne(C++) < ne(Py) systematisch.
        # Das ist ein ARCHITEKTURUNTERSCHIED, kein Bug.
        # Wir tracken das Ratio um Drift zu erkennen.
        ne_ratio = c["ne"] / p["ne"] if p["ne"] > 0 else 0
        ng_ratio = c["ng"] / p["ng"] if p["ng"] > 0 else 0

        # ne-Ratio: C++ hat weniger P_abs (RF-BVP), Ratio variiert mit Q0
        check(f"{prefix}: ne ratio stable",
              0.10 < ne_ratio < 0.70,
              f"C++/Py={ne_ratio:.3f} (C++={c['ne']:.2e}, Py={p['ne']:.2e})")

        # ng-Ratio: variiert leicht mit Betriebspunkt. Obergrenze am 2026-09-02
        # von 0.85 auf 0.90 angehoben: die C++-Seite rechnet seit der Umstellung
        # der Abbruchschranke von 1e-2 auf 1e-4 genauer, wodurch sich das
        # Verhaeltnis an den oberen Sweep-Punkten auf 0.867 verschoben hat. Die
        # Python-Seite laeuft weiterhin mit ihrer eigenen, deutlich lockereren
        # Schranke; das Verhaeltnis bleibt ein Driftwaechter, kein Massstab.
        check(f"{prefix}: ng ratio stable",
              0.25 < ng_ratio < 0.90,
              f"C++/Py={ng_ratio:.3f} (C++={c['ng']:.2e}, Py={p['ng']:.2e})")

        # I_beam: C++ Bohm-direkt (hoch), Python Bohm+CL+eta_opt (niedrig)
        I_ratio = c["I_mA"] / p["I_mA"] if p["I_mA"] > 0 else 0
        note(f"{prefix}: I_beam ratio",
             f"C++/Py={I_ratio:.2f} (C++={c['I_mA']:.2f}mA Py={p['I_mA']:.2f}mA, "
             f"Restdifferenz: P_abs-Offset durch RF-BVP)")


# ═════════════════════════════════════════════════════════════
# Hauptprogramm
# ═════════════════════════════════════════════════════════════
def main():
    global passed, failed, skipped

    filter_ids = sys.argv[1:] if len(sys.argv) > 1 else []

    # C++ Binary suchen
    bin_type, bin_path = find_cpp_binary()
    if not bin_type:
        print("SKIP: C++ Binary 'chabert' nicht gefunden. Cross-Backend-Tests uebersprungen.")
        print("      Zum Kompilieren: GUI -> Kompilieren, oder manuell:")
        print("      g++ -O3 -std=c++17 -c bessel_wrapper.cpp sim_config.cpp rates.cpp "
              "physics.cpp solver.cpp sim_logging.cpp")
        print("      g++ -O3 -std=c++17 main.cpp *.o -o chabert")
        return 0  # Skip, kein Fail

    # Chemistry pruefen
    chem_path = SCRIPT_DIR / CHEM_XE
    if not chem_path.exists():
        print(f"SKIP: {CHEM_XE} nicht gefunden.")
        return 0

    cases = CASES
    if filter_ids:
        cases = [c for c in cases if any(f in c["id"] for f in filter_ids)]

    print(f"Cross-Backend Regression: C++ vs Python")
    print(f"  Binary: {bin_type} ({bin_path})")
    print(f"  Faelle: {len(cases)}")
    print()
    print("  Systematische Unterschiede (dokumentiert):")
    print("    P_abs: C++ berechnet via RF-BVP, Python nimmt P direkt")
    print("    I_beam: C++ Bohm-Fluss, Python Bohm+CL+eta_opt")
    print("    ne: Gekoppelt an P_abs, daher verschoben")
    print("  Fair vergleichbar: Te, Tg (selbe Bilanzgleichungen)")
    print("=" * 70)

    for case in cases:
        cid = case["id"]
        print(f"\n--- {case['name']} ---")
        print(f"    {case['why']}")

        # RunConfig erzeugen und schreiben
        g = case["geom"]
        ex = case["cpp_extra"]
        cfg = make_config(
            R=g["R"], L=g["L"], betai=g["betai"], betag=g["betag"],
            Vgrid=g["Vgrid"], sgrid=g["sgrid"], eta_opt=0.25,
            frequency=g["frequency"], Nw=g.get("Nw", 6), R_ohm=g.get("R_ohm", 0.36),
            Rc=g.get("Rc", 0.02), lc=g.get("lc", 0.04),
            solve_mode=ex.get("solve_mode", 2),
            rate_model=ex.get("rate_model", 0),
            Q0_start=ex.get("Q0sccm_start", 0.4),
            Q0_step=ex.get("Q0sccm_step", 0.05),
            N=ex.get("jjmax", 3),
        )
        write_run_config(cfg)

        # C++ ausfuehren
        cpp_out = run_cpp(bin_type, bin_path)
        cpp_res = parse_results(cpp_out)
        check(f"{cid}: C++ converged", len(cpp_res) > 0,
              f"results={len(cpp_res)}")

        # Python ausfuehren
        q0s, q0d, N = case["py_sweep"]
        py_out = run_python(CHEM_XE, case["py_mode"], case["py_param"],
                            q0s, q0d, N)
        py_res = parse_results(py_out)
        check(f"{cid}: Python converged", len(py_res) > 0,
              f"results={len(py_res)}")

        # Vergleichen
        if cpp_res and py_res:
            compare_results(cid, cpp_res, py_res, case["tolerances"])

        # Konvergenzstatus
        if cpp_res and py_res:
            check(f"{cid}: both backends converged",
                  len(cpp_res) > 0 and len(py_res) > 0)

    # Aufraeumen
    _cleanup()
    (SCRIPT_DIR / "xb_params.txt").unlink(missing_ok=True)

    print(f"\n{'='*70}")
    print(f"  Ergebnis: {passed} passed, {failed} failed")
    if errors:
        print(f"  Fehlgeschlagen:")
        for e in errors:
            print(f"    - {e}")
    else:
        print("  Beide Backends konsistent (innerhalb dokumentierter Toleranzen)")
    print(f"{'='*70}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
