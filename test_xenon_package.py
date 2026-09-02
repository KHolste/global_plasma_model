#!/usr/bin/env python3
"""
test_xenon_package.py -- Das Xenon-Paket aus den Biagi-Daten gegen den
fest verdrahteten tabellierten Weg.

Das Paket führt die fünfzig Anregungsprozesse einzeln, der fest verdrahtete
Weg fasst sie zu einem Term zusammen. Beides muss auf dasselbe hinauslaufen:
die Summe der Einzelraten ist die Sammelrate, die Summe der Einzelverluste
ist der Sammelverlust, und beide Rechenwege müssen denselben Betriebspunkt
finden.

Geprüft wird bei enger Abbruchschranke. Bei der voreingestellten Schranke von
einem Prozent unterscheiden sich die beiden Wege um mehrere Prozent, weil die
Lösung dort noch gar nicht festliegt -- das ist eine Aussage über die Schranke
und nicht über die Chemie.

Aufruf: python test_xenon_package.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PAKET = SCRIPT_DIR / "chemistry" / "xenon_biagi"
PARAMS = SCRIPT_DIR / "xenon_package_params.txt"
KEX_TABELLE = SCRIPT_DIR / "cross_sections" / "xenon" / "biagi" / "kex_table.csv"
E_CH = 1.602176487e-19

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
rate_model 2
newton_tol 1e-4
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


def lies_tabelle(pfad: Path) -> tuple[list[float], list[float]]:
    Te, K = [], []
    for zeile in pfad.read_text(encoding="utf-8").splitlines():
        if not zeile or zeile.startswith("#") or zeile.startswith("Te"):
            continue
        teile = zeile.split(",")
        Te.append(float(teile[0]))
        K.append(float(teile[1]))
    return Te, K


def test_aufbau() -> dict:
    """Struktur des Pakets."""
    datei = PAKET / "chemistry.json"
    check("Paket vorhanden", datei.exists(), str(datei))
    if not datei.exists():
        return {}
    pkg = json.loads(datei.read_text(encoding="utf-8"))

    arten = [r["type"] for r in pkg["reactions"]]
    check("ein elastischer Prozess", arten.count("elastic") == 1)
    check("ein ionisierender Prozess", arten.count("ionization") == 1)
    check("50 Anregungsprozesse", arten.count("excitation") == 50,
          str(arten.count("excitation")))
    check("alle Raten tabelliert",
          all(r["rate"]["model"] == "tabulated" for r in pkg["reactions"]
              if r["type"] != "surface_recombination"))
    check("jede Anregung mit eigener Schwelle",
          len({r["energy_eV"] for r in pkg["reactions"] if r["type"] == "excitation"}) > 40)
    check("jede Anregung mit eigener Tabelle",
          len({r["rate"]["file"] for r in pkg["reactions"]}) == len(pkg["reactions"]))
    for r in pkg["reactions"]:
        pfad = PAKET / r["rate"]["file"]
        if not pfad.exists():
            check(f"Tabelle {r['id']} vorhanden", False, str(pfad))
            return pkg
    check("alle Ratentabellen vorhanden", True)
    return pkg


def test_summe(pkg: dict) -> None:
    """Summe der Einzelprozesse gegen die vorhandene Sammeltabelle."""
    if not pkg or not KEX_TABELLE.exists():
        check("Sammeltabelle vorhanden", False, str(KEX_TABELLE))
        return

    Te_ref, K_ref, P_ref = [], [], []
    for zeile in KEX_TABELLE.read_text(encoding="utf-8").splitlines():
        if not zeile or zeile.startswith("#") or zeile.startswith("Te"):
            continue
        a, b, c = zeile.split(",")
        Te_ref.append(float(a))
        K_ref.append(float(b))
        P_ref.append(float(c))

    summe_K = [0.0] * len(Te_ref)
    summe_P = [0.0] * len(Te_ref)
    for r in pkg["reactions"]:
        if r["type"] != "excitation":
            continue
        Te, K = lies_tabelle(PAKET / r["rate"]["file"])
        check(f"{r['id']}: gleiches Te-Gitter", Te == Te_ref, f"{len(Te)} Punkte")
        if Te != Te_ref:
            return
        for i, k in enumerate(K):
            summe_K[i] += k
            summe_P[i] += k * r["energy_eV"] * E_CH

    d_K = max(abs(a - b) / b for a, b in zip(summe_K, K_ref) if b > 0)
    d_P = max(abs(a - b) / b for a, b in zip(summe_P, P_ref) if b > 0)
    check("Summe der Anregungsraten trifft die Sammeltabelle", d_K < 1e-5, f"{d_K:.2e}")
    check("Summe der Energieverluste trifft die Sammeltabelle", d_P < 1e-5, f"{d_P:.2e}")
    print(f"  groesste Abweichung: Raten {d_K:.1e}, Energieverluste {d_P:.1e}")


def lauf(prog: Path, zusatz: str) -> list[list[float]]:
    PARAMS.write_text(BASIS + zusatz, encoding="utf-8")
    r = subprocess.run([str(prog), PARAMS.name], cwd=SCRIPT_DIR,
                       capture_output=True, text=True, timeout=600)
    return [[float(x) for x in z.split()[1:7]]
            for z in r.stdout.splitlines() if z.startswith("RESULT ")]


def test_gleiches_ergebnis(prog: Path) -> None:
    """Beide Wege muessen denselben Betriebspunkt finden."""
    fest = lauf(prog, "")
    paket = lauf(prog, "chemistry_package xenon_biagi\n")

    check("fester Weg liefert zwei Punkte", len(fest) == 2, str(len(fest)))
    check("Paketweg liefert zwei Punkte", len(paket) == 2, str(len(paket)))
    if len(fest) != len(paket) or not fest:
        return

    namen = ["Plasmadichte", "Neutraldichte", "Te", "Tg", "Strahlstrom"]
    for k, (a, b) in enumerate(zip(fest, paket)):
        for i, name in enumerate(namen):
            rel = abs(a[i] - b[i]) / max(abs(a[i]), 1e-30)
            check(f"Punkt {k+1}, {name} stimmt ueberein (<0.5%)", rel < 5e-3,
                  f"{a[i]:.6g} vs {b[i]:.6g} ({rel*100:.3f}%)")


def main() -> int:
    prog = kern()
    if prog is None:
        print("chabert fehlt. Mit 'python build.py' uebersetzen.", file=sys.stderr)
        return 1
    try:
        pkg = test_aufbau()
        test_summe(pkg)
        test_gleiches_ergebnis(prog)
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
