#!/usr/bin/env python3
"""
noble_gas_comparison.py -- Vergleichsrechnungen fuer Xenon, Krypton und Argon.

Alle drei Edelgase laufen ueber denselben generischen Weg, dieselbe Geometrie
und denselben Betriebspunkt; unterschiedlich ist allein das Chemiepaket. Damit
ist der Vergleich einer der Stoffdaten und nicht einer der Rechenwege.

Zwei Auswertungen:

  1. Gasvergleich       -- Xe, Kr und Ar aus der Biagi-Datenbank ueber einen
                           Durchflussbereich bei fester Generatorleistung.
  2. Datenbankvergleich -- je Gas alle vorhandenen vollstaendigen Datenbanken
                           an einem Punkt. Zeigt, wieviel allein die Wahl der
                           Querschnittsquelle ausmacht.

Der Durchfluss wird von unten nach oben durchfahren. Der Kern faehrt sich
seine Durchfluss-Leiter bei Bedarf selbst, deshalb tragen auch die kalten
Startpunkte der leichten Gase.

Aufruf:
    python noble_gas_comparison.py
    python noble_gas_comparison.py --csv ergebnis.csv
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PARAMS = SCRIPT_DIR / "vergleich_params.txt"

#: Geometrie und Spule des Standardtriebwerks, fuer alle Laeufe gleich.
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
P_RFG_max 200.0
solve_mode 2
rate_model 0
"""

P_RFG = 18.0          # Generatorleistung [W]
Q_UNTEN = 0.25        # unterster Durchfluss [sccm]
Q_SCHRITT = 0.05
Q_ANZAHL = 12
Q_PUNKT = 0.50        # Punkt fuer den Datenbankvergleich [sccm]

#: Gas -> Datenbanken mit vollstaendigem Satz und eindeutiger Zuordnung.
#: Die erste ist die, mit der der Gasvergleich gefahren wird.
DATENBANKEN = {
    "xenon":   ["biagi", "hayashi"],
    "krypton": ["biagi", "biagi_v7_1", "morgan", "siglo"],
    "argon":   ["biagi", "biagi_v7_1", "bsr", "hayashi", "iaa", "ist_lisbon",
                "morgan", "puech"],
}

SYMBOL = {"xenon": "Xe", "krypton": "Kr", "argon": "Ar"}


class Punkt:
    """Ein konvergierter Betriebspunkt."""
    __slots__ = ("Q0", "P", "n", "ng", "Te", "Tg", "I", "F", "eta_m", "xi", "icp")

    def __init__(self, Q0, P, n, ng, Te, Tg, I):
        self.Q0, self.P = Q0, P
        self.n, self.ng, self.Te, self.Tg, self.I = n, ng, Te, Tg, I
        self.F = self.eta_m = self.xi = self.icp = float("nan")


def kern() -> Path:
    for name in ("chabert.exe", "chabert"):
        p = SCRIPT_DIR / name
        if p.exists():
            return p
    print("chabert fehlt. Mit 'python build.py' uebersetzen.", file=sys.stderr)
    sys.exit(1)


def lauf(prog: Path, gas: str, db: str) -> list[Punkt]:
    """Einen Durchflusslauf rechnen und die konvergierten Punkte zurueckgeben."""
    text = (BASIS
            + f"P_RFG {P_RFG}\n"
            + f"Q0sccm_start {Q_UNTEN}\n"
            + f"Q0sccm_step {Q_SCHRITT}\n"
            + f"jjmax {Q_ANZAHL}\n"
            + f"gas_species {gas}\n"
            + f"chemistry_package {gas}_{db}\n")
    PARAMS.write_text(text, encoding="utf-8")
    r = subprocess.run([str(prog), PARAMS.name], cwd=SCRIPT_DIR,
                       capture_output=True, text=True, timeout=1800)

    # Ein fehlendes oder fehlerhaftes Paket laesst den Kern auf die fest
    # verdrahtete Xenon-Physik zurueckfallen. Das Ergebnis saehe dann nach
    # einem Lauf dieses Gases aus und waere keiner. Hier ist das ein Abbruch.
    if "SOLVER_PATH generic" not in r.stdout:
        print(f"FEHLER: {gas}/{db} lief nicht ueber das Chemiepaket, sondern "
              f"ueber den fest verdrahteten Weg.\n{r.stderr.strip()[:400]}",
              file=sys.stderr)
        sys.exit(1)

    punkte: list[Punkt] = []
    for zeile in r.stdout.splitlines():
        f = zeile.split()
        if zeile.startswith("RESULT ") and len(f) >= 7:
            punkte.append(Punkt(float("nan"), float(f[6]), float(f[1]), float(f[2]),
                                float(f[3]), float(f[4]), float(f[5])))
        elif zeile.startswith("RESULT_EXT ") and punkte and len(f) >= 13:
            p = punkte[-1]
            p.Q0 = float(f[1])
            p.F = float(f[8])       # Gesamtschub [mN]
            p.icp = float(f[9])     # Einkoppelwirkungsgrad
            p.xi = float(f[11])     # Schub je Leistung [mN/kW]
            p.eta_m = float(f[12])  # Massenwirkungsgrad
    return punkte


def bei(punkte: list[Punkt], Q0: float) -> Punkt | None:
    for p in punkte:
        if abs(p.Q0 - Q0) < 1e-6:
            return p
    return None


def anregungszahl(gas: str, db: str) -> int:
    pfad = SCRIPT_DIR / "chemistry" / f"{gas}_{db}" / "chemistry.json"
    if not pfad.exists():
        return 0
    paket = json.loads(pfad.read_text(encoding="utf-8"))
    return sum(1 for r in paket["reactions"] if r["type"] == "excitation")


def tabelle_gasvergleich(daten: dict[str, list[Punkt]]) -> list[str]:
    z = [f"Gasvergleich, Biagi-Querschnitte, P_RFG = {P_RFG:.0f} W", ""]
    spalte = f"{'Te':>6} {'n_e':>10} {'I':>7} {'F':>6} {'eta_m':>6}"
    einheit = f"{'eV':>6} {'1/m^3':>10} {'mA':>7} {'mN':>6} {'-':>6}"
    kopf = f"{'Q0':>6} | " + " | ".join(f"{SYMBOL[g]:^38}" for g in daten)
    z.append(kopf)
    z.append(f"{'sccm':>6} | " + " | ".join(spalte for _ in daten))
    z.append(f"{'':>6} | " + " | ".join(einheit for _ in daten))
    z.append("-" * len(kopf))

    alle_q = sorted({round(p.Q0, 4) for ps in daten.values() for p in ps}, reverse=True)
    leer = f"{'--':>6} {'--':>10} {'--':>7} {'--':>6} {'--':>6}"
    for q in alle_q:
        teile = []
        for ps in daten.values():
            p = bei(ps, q)
            teile.append(leer if p is None else
                         f"{p.Te:6.2f} {p.n:10.3e} {p.I:7.2f} {p.F:6.3f} {p.eta_m:6.3f}")
        z.append(f"{q:6.2f} | " + " | ".join(teile))
    return z


def spanne_prozent(werte: list[float]) -> float:
    """Groesster Wert gegen kleinsten, in Prozent."""
    klein, gross = min(werte), max(werte)
    return (gross / klein - 1.0) * 100.0 if klein > 0 else float("nan")


def tabelle_datenbanken(daten: dict[str, dict[str, Punkt | None]]) -> list[str]:
    z = [f"Datenbankvergleich bei Q0 = {Q_PUNKT:.2f} sccm, P_RFG = {P_RFG:.0f} W", ""]
    z.append(f"{'Gas':<9}{'Datenbank':<12}{'Anreg':>6}{'Te':>8}{'n_e':>12}"
             f"{'n_g':>12}{'I':>9}{'F':>8}{'eta_m':>8}{'xi':>9}")
    z.append(f"{'':<9}{'':<12}{'':>6}{'eV':>8}{'1/m^3':>12}"
             f"{'1/m^3':>12}{'mA':>9}{'mN':>8}{'-':>8}{'mN/kW':>9}")
    z.append("-" * 93)
    for gas, dbs in daten.items():
        werte = [p for p in dbs.values() if p is not None]
        for i, (db, p) in enumerate(dbs.items()):
            name = gas if i == 0 else ""
            n_exc = anregungszahl(gas, db)
            if p is None:
                z.append(f"{name:<9}{db:<12}{n_exc:>6}{'nicht konvergiert':>56}")
            else:
                z.append(f"{name:<9}{db:<12}{n_exc:>6}{p.Te:8.2f}{p.n:12.3e}"
                         f"{p.ng:12.3e}{p.I:9.2f}{p.F:8.3f}{p.eta_m:8.3f}{p.xi:9.2f}")
        if len(werte) > 1:
            z.append(f"{'':<9}{'Spanne':<12}{'':>6}"
                     f"{spanne_prozent([p.Te for p in werte]):7.1f}%"
                     f"{spanne_prozent([p.n for p in werte]):11.1f}%"
                     f"{spanne_prozent([p.ng for p in werte]):11.1f}%"
                     f"{spanne_prozent([p.I for p in werte]):8.1f}%"
                     f"{spanne_prozent([p.F for p in werte]):7.1f}%"
                     f"{spanne_prozent([p.eta_m for p in werte]):7.1f}%"
                     f"{spanne_prozent([p.xi for p in werte]):8.1f}%")
        z.append("")
    return z


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Vergleichsrechnungen fuer die drei Edelgase.")
    ap.add_argument("--csv", default=None, help="Ergebnisse zusaetzlich als CSV")
    args = ap.parse_args(argv)

    prog = kern()
    alle: dict[tuple[str, str], list[Punkt]] = {}
    try:
        for gas, dbs in DATENBANKEN.items():
            for db in dbs:
                print(f"rechne {gas} / {db} ...", file=sys.stderr)
                alle[(gas, db)] = lauf(prog, gas, db)
    finally:
        PARAMS.unlink(missing_ok=True)

    gasdaten = {gas: alle[(gas, dbs[0])] for gas, dbs in DATENBANKEN.items()}
    dbdaten = {gas: {db: bei(alle[(gas, db)], Q_PUNKT) for db in dbs}
               for gas, dbs in DATENBANKEN.items()}

    print("\n".join(tabelle_gasvergleich(gasdaten) + [""]
                    + tabelle_datenbanken(dbdaten)))

    if args.csv:
        with open(args.csv, "w", encoding="utf-8", newline="\n") as f:
            f.write("gas,datenbank,Q0_sccm,P_RFG_W,Te_eV,n_e_m3,n_g_m3,Tg_K,"
                    "I_mA,F_mN,eta_mass,xi_mN_kW,icp_eff\n")
            for (gas, db), ps in alle.items():
                for p in ps:
                    f.write(f"{gas},{db},{p.Q0},{p.P},{p.Te},{p.n},{p.ng},"
                            f"{p.Tg},{p.I},{p.F},{p.eta_m},{p.xi},{p.icp}\n")
        print(f"\nCSV geschrieben: {args.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
