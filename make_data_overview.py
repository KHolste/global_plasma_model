#!/usr/bin/env python3
"""
make_data_overview.py -- Übersicht über alle vorhandenen Wirkungsquerschnitte
und Ratenkoeffizienten, nach Gassorte geordnet.

Erzeugt unter data_overview/ je Gassorte einen Ordner mit

    cross_sections/<datenbank>/<prozess>.csv    Energie [eV], Querschnitt [m^2]
    rate_coefficients/<quelle>/<prozess>.csv    Te [eV], K [m^3/s]
    figures/...                                 je eine Abbildung pro Prozess
                                                plus Übersichtsbilder

Die CSV-Dateien haben eine Kopfzeile mit der Quelle und danach zwei Spalten.
Analytisch gegebene Ratenkoeffizienten werden auf einem festen Te-Gitter
ausgewertet, damit sie im selben Format vorliegen.

Aufruf:
    python make_data_overview.py
    python make_data_overview.py --nur-csv     (ohne Abbildungen)
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
ZIEL = SCRIPT_DIR / "data_overview"

#: Te-Gitter für analytisch gegebene Ratenkoeffizienten.
TE_GITTER = [0.5 + 0.05 * i for i in range(391)]


# ═══════════════════════════════════════════════════════════════
# Lesen
# ═══════════════════════════════════════════════════════════════

def lies_zwei_spalten(pfad: Path) -> tuple[list[float], list[float]]:
    x, y = [], []
    for zeile in pfad.read_text(encoding="utf-8").splitlines():
        z = zeile.strip()
        if not z or z.startswith("#") or z[0].isalpha():
            continue
        teile = z.split(",")
        if len(teile) < 2:
            continue
        try:
            x.append(float(teile[0]))
            y.append(float(teile[1]))
        except ValueError:
            continue
    return x, y


def schreibe_csv(pfad: Path, kopf: str, spalten: tuple[str, str],
                 x: list[float], y: list[float]) -> None:
    pfad.parent.mkdir(parents=True, exist_ok=True)
    with open(pfad, "w", encoding="utf-8", newline="\n") as f:
        f.write(f"# {kopf}\n")
        f.write(f"{spalten[0]},{spalten[1]}\n")
        for a, b in zip(x, y):
            f.write(f"{a:.6g},{b:.6e}\n")


def dateiname(text: str) -> str:
    """Aus einer Prozessbezeichnung einen brauchbaren Dateinamen machen."""
    erlaubt = []
    for c in text:
        if c.isalnum() or c in "._-+":
            erlaubt.append(c)
        elif c in " /\\":
            erlaubt.append("_")
    name = "".join(erlaubt).strip("_")
    while "__" in name:
        name = name.replace("__", "_")
    return name or "prozess"


# ═══════════════════════════════════════════════════════════════
# Abbildungen
# ═══════════════════════════════════════════════════════════════

def bild_einzeln(pfad: Path, x, y, titel: str, untertitel: str,
                 xlabel: str, ylabel: str, logx: bool) -> None:
    import matplotlib.pyplot as plt

    pfad.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.plot(x, y, lw=1.4)
    if logx:
        ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(titel, fontsize=10, pad=20)
    if untertitel:
        ax.text(0.0, 1.01, untertitel, transform=ax.transAxes,
                fontsize=7.5, color="0.35", va="bottom")
    ax.grid(True, which="both", lw=0.3, alpha=0.5)
    fig.tight_layout()
    fig.savefig(pfad, dpi=130)
    plt.close(fig)


def bild_uebersicht(pfad: Path, kurven: list[tuple[str, list, list]],
                    titel: str, untertitel: str,
                    xlabel: str, ylabel: str, logx: bool) -> None:
    import matplotlib.pyplot as plt

    if not kurven:
        return
    pfad.parent.mkdir(parents=True, exist_ok=True)
    breite = 9.0 if len(kurven) > 12 else 7.0
    fig, ax = plt.subplots(figsize=(breite, 5.2))
    for name, x, y in kurven:
        ax.plot(x, y, lw=1.0, label=name)
    if logx:
        ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(titel, fontsize=11, pad=22)
    if untertitel:
        ax.text(0.0, 1.01, untertitel, transform=ax.transAxes,
                fontsize=8, color="0.35", va="bottom")
    ax.grid(True, which="both", lw=0.3, alpha=0.5)
    spalten = 1 if len(kurven) <= 14 else (2 if len(kurven) <= 34 else 3)
    ax.legend(fontsize=5.5 if len(kurven) > 20 else 7, ncol=spalten,
              loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
    fig.tight_layout()
    fig.savefig(pfad, dpi=130, bbox_inches="tight")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════
# Wirkungsquerschnitte
# ═══════════════════════════════════════════════════════════════

def quelle_der_datenbank(db_dir: Path) -> str:
    info = db_dir / "db_info.json"
    if info.exists():
        d = json.loads(info.read_text(encoding="utf-8"))
        teile = [d.get("display_name", db_dir.name), d.get("source", "")]
        if d.get("retrieved"):
            teile.append(f"abgerufen {d['retrieved']}")
        return ", ".join(t for t in teile if t)
    return db_dir.name


def querschnitte(gas: str, gas_dir: Path, bilder: bool) -> int:
    anzahl = 0
    for db_dir in sorted(p for p in gas_dir.iterdir() if p.is_dir()):
        meta = db_dir / "metadata.json"
        if not meta.exists():
            continue
        quelle = quelle_der_datenbank(db_dir)
        prozesse = json.loads(meta.read_text(encoding="utf-8"))["processes"]

        kurven = []
        for p in prozesse:
            if not p.get("file"):
                continue          # Prozesstyp, den der Wandler nicht ablegt
            quelldatei = db_dir / p["file"]
            if not quelldatei.exists():
                continue
            x, y = lies_zwei_spalten(quelldatei)
            if len(x) < 2:
                continue

            art = p["type"].lower()
            zustand = p.get("state")
            schwelle = p.get("threshold_eV")
            bezeichnung = p.get("process", art)
            kurz = f"{art}_{dateiname(zustand)}" if zustand else art
            if schwelle:
                bezeichnung += f", Schwelle {schwelle} eV"

            schreibe_csv(ZIEL / gas / "cross_sections" / db_dir.name / f"{dateiname(kurz)}.csv",
                         f"{bezeichnung} | {quelle}",
                         ("energy_eV", "cross_section_m2"), x, y)
            anzahl += 1

            if bilder:
                bild_einzeln(
                    ZIEL / gas / "figures" / "cross_sections" / db_dir.name / f"{dateiname(kurz)}.png",
                    x, y, bezeichnung, quelle,
                    "Elektronenenergie [eV]", "Wirkungsquerschnitt [m$^2$]", logx=True)
            kurven.append((zustand or art, x, y))

        if bilder and kurven:
            bild_uebersicht(
                ZIEL / gas / "figures" / f"uebersicht_querschnitte_{db_dir.name}.png",
                kurven, f"{gas.capitalize()}: alle Wirkungsquerschnitte ({len(kurven)} Prozesse)",
                quelle, "Elektronenenergie [eV]", "Wirkungsquerschnitt [m$^2$]", logx=True)
    return anzahl


# ═══════════════════════════════════════════════════════════════
# Ratenkoeffizienten
# ═══════════════════════════════════════════════════════════════

def werte_rate_aus(rate: dict, paket_dir: Path) -> tuple[list[float], list[float], str] | None:
    """Ratenkoeffizient auf ein Te-Gitter bringen, unabhaengig vom Modell."""
    modell = rate.get("model")
    if modell == "tabulated":
        pfad = paket_dir / rate.get("file", "")
        if not pfad.exists():
            return None
        x, y = lies_zwei_spalten(pfad)
        return (x, y, "tabelliert") if len(x) >= 2 else None
    if modell == "arrhenius":
        A, Ea = float(rate["A"]), float(rate["E_a_eV"])
        y = [A * math.exp(-Ea / Te) for Te in TE_GITTER]
        return TE_GITTER, y, f"Arrhenius, A={A:.4g} m^3/s, E_a={Ea:.4g} eV"
    if modell == "constant":
        v = float(rate["value"])
        return TE_GITTER, [v] * len(TE_GITTER), f"konstant, {v:.4g} m^3/s"
    if modell == "polynomial":
        c = [float(k) for k in rate["coeffs"]]
        y = [sum(k * Te**i for i, k in enumerate(c)) for Te in TE_GITTER]
        return TE_GITTER, y, "Polynom"
    return None


def raten_aus_paket(gas: str, paket_dir: Path, bilder: bool) -> int:
    pkg = json.loads((paket_dir / "chemistry.json").read_text(encoding="utf-8"))
    quelle = pkg.get("name", paket_dir.name)
    anzahl = 0
    kurven = []

    for rxn in pkg["reactions"]:
        ergebnis = werte_rate_aus(rxn["rate"], paket_dir)
        if ergebnis is None:
            continue
        x, y, art = ergebnis
        if not any(v > 0 for v in y):
            continue
        bezeichnung = rxn.get("name", rxn["id"])
        energie = rxn.get("energy_eV")
        if energie:
            bezeichnung += f", Schwelle {energie} eV"

        schreibe_csv(ZIEL / gas / "rate_coefficients" / paket_dir.name / f"{dateiname(rxn['id'])}.csv",
                     f"{bezeichnung} | {quelle} | {art}",
                     ("Te_eV", "K_m3s"), x, y)
        anzahl += 1

        if bilder:
            bild_einzeln(
                ZIEL / gas / "figures" / "rate_coefficients" / paket_dir.name / f"{dateiname(rxn['id'])}.png",
                x, y, bezeichnung, f"{quelle} ({art})",
                "Elektronentemperatur [eV]", "Ratenkoeffizient [m$^3$/s]", logx=False)
        kurven.append((rxn["id"], x, y))

    if bilder and kurven:
        bild_uebersicht(
            ZIEL / gas / "figures" / f"uebersicht_raten_{paket_dir.name}.png",
            kurven, f"{gas.capitalize()}: alle Ratenkoeffizienten ({len(kurven)} Reaktionen)",
            quelle, "Elektronentemperatur [eV]", "Ratenkoeffizient [m$^3$/s]", logx=False)
    return anzahl


def raten_aus_varianten(gas: str, paket_dir: Path, bilder: bool) -> int:
    """Alternative Ratensaetze, die als Varianten danebenliegen."""
    var_dir = paket_dir / "variants"
    if not var_dir.is_dir():
        return 0
    anzahl = 0
    for prozess_dir in sorted(p for p in var_dir.iterdir() if p.is_dir()):
        kurven = []
        for datei in sorted(prozess_dir.glob("*.csv")):
            x, y = lies_zwei_spalten(datei)
            if len(x) < 2:
                continue
            name = f"{prozess_dir.name}_{datei.stem}"
            schreibe_csv(
                ZIEL / gas / "rate_coefficients" / f"{paket_dir.name}_varianten" / f"{dateiname(name)}.csv",
                f"{prozess_dir.name}, Variante {datei.stem} | {paket_dir.name}",
                ("Te_eV", "K_m3s"), x, y)
            anzahl += 1
            if bilder:
                bild_einzeln(
                    ZIEL / gas / "figures" / "rate_coefficients" /
                    f"{paket_dir.name}_varianten" / f"{dateiname(name)}.png",
                    x, y, f"{prozess_dir.name}, Variante {datei.stem}", paket_dir.name,
                    "Elektronentemperatur [eV]", "Ratenkoeffizient [m$^3$/s]", logx=False)
            kurven.append((datei.stem, x, y))
        if bilder and len(kurven) > 1:
            bild_uebersicht(
                ZIEL / gas / "figures" / f"varianten_{dateiname(prozess_dir.name)}.png",
                kurven, f"{prozess_dir.name}: Varianten im Vergleich", paket_dir.name,
                "Elektronentemperatur [eV]", "Ratenkoeffizient [m$^3$/s]", logx=False)
    return anzahl


def raten_aus_tabellen(gas: str, gas_dir: Path, bilder: bool) -> int:
    """Die vorab gerechneten Sammeltabellen neben den Querschnitten."""
    benennung = {
        "kiz_table.csv": ("Ionisation, integriert ueber die Querschnitte", 1),
        "kel_table.csv": ("Elastischer Stoss, integriert ueber die Querschnitte", 1),
        "kex_table.csv": ("Anregung, Summe aller Prozesse", 1),
    }
    anzahl = 0
    for db_dir in sorted(p for p in gas_dir.iterdir() if p.is_dir()):
        quelle = quelle_der_datenbank(db_dir)
        kurven = []
        for datei, (bezeichnung, spalte) in benennung.items():
            pfad = db_dir / datei
            if not pfad.exists():
                continue
            x, y = lies_zwei_spalten(pfad)
            if len(x) < 2:
                continue
            kurz = datei.replace("_table.csv", "")
            schreibe_csv(ZIEL / gas / "rate_coefficients" / f"{db_dir.name}_sammeltabellen" / f"{kurz}.csv",
                         f"{bezeichnung} | {quelle}", ("Te_eV", "K_m3s"), x, y)
            anzahl += 1
            if bilder:
                bild_einzeln(
                    ZIEL / gas / "figures" / "rate_coefficients" /
                    f"{db_dir.name}_sammeltabellen" / f"{kurz}.png",
                    x, y, bezeichnung, quelle,
                    "Elektronentemperatur [eV]", "Ratenkoeffizient [m$^3$/s]", logx=False)
            kurven.append((kurz, x, y))
        if bilder and kurven:
            bild_uebersicht(
                ZIEL / gas / "figures" / f"uebersicht_sammeltabellen_{db_dir.name}.png",
                kurven, f"{gas.capitalize()}: zusammengefasste Ratenkoeffizienten",
                quelle, "Elektronentemperatur [eV]", "Ratenkoeffizient [m$^3$/s]", logx=False)
    return anzahl


# ═══════════════════════════════════════════════════════════════
# Hauptlauf
# ═══════════════════════════════════════════════════════════════

def gas_eines_pakets(paket_dir: Path) -> str:
    meta = paket_dir / "metadata.json"
    if meta.exists():
        d = json.loads(meta.read_text(encoding="utf-8"))
        if d.get("gas"):
            return str(d["gas"])
    name = paket_dir.name.lower()
    for gas in ("xenon", "iodine", "argon", "krypton"):
        if name.startswith(gas) or gas in name:
            return gas
    return name


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Datenübersicht erzeugen.")
    ap.add_argument("--nur-csv", action="store_true", help="ohne Abbildungen")
    args = ap.parse_args(argv)
    bilder = not args.nur_csv

    if bilder:
        try:
            import matplotlib
            matplotlib.use("Agg")
        except ImportError:
            print("matplotlib fehlt -- erzeuge nur die CSV-Dateien.", file=sys.stderr)
            bilder = False

    if ZIEL.exists():
        shutil.rmtree(ZIEL)
    ZIEL.mkdir(parents=True)

    bestand: dict[str, dict[str, int]] = {}

    # Gassorten aus beiden Quellen sammeln
    cs_wurzel = SCRIPT_DIR / "cross_sections"
    gase = {p.name for p in cs_wurzel.iterdir()
            if p.is_dir() and p.name != "tests"} if cs_wurzel.is_dir() else set()
    pakete = [p for p in (SCRIPT_DIR / "chemistry").iterdir()
              if p.is_dir() and (p / "chemistry.json").exists()]
    for p in pakete:
        gase.add(gas_eines_pakets(p))

    for gas in sorted(gase):
        eintrag = {"querschnitte": 0, "raten": 0}
        gas_dir = cs_wurzel / gas
        if gas_dir.is_dir():
            eintrag["querschnitte"] += querschnitte(gas, gas_dir, bilder)
            eintrag["raten"] += raten_aus_tabellen(gas, gas_dir, bilder)
        for paket in sorted(pakete):
            if gas_eines_pakets(paket) != gas:
                continue
            eintrag["raten"] += raten_aus_paket(gas, paket, bilder)
            eintrag["raten"] += raten_aus_varianten(gas, paket, bilder)
        bestand[gas] = eintrag

        ordner = ZIEL / gas
        ordner.mkdir(parents=True, exist_ok=True)
        if eintrag["querschnitte"] == 0 and eintrag["raten"] == 0:
            (ordner / "README.md").write_text(
                f"# {gas.capitalize()}\n\n"
                "Für dieses Gas liegen im Projekt weder Wirkungsquerschnitte noch\n"
                "Ratenkoeffizienten vor. Im Rechenkern sind lediglich die\n"
                "Stoffkonstanten hinterlegt (Masse, Ionisierungs- und\n"
                "Anregungsenergie, Wärmeleitfähigkeit).\n",
                encoding="utf-8", newline="\n")
        print(f"  {gas:<10} {eintrag['querschnitte']:>3} Querschnitte, "
              f"{eintrag['raten']:>3} Ratenkoeffizienten")

    zeilen = ["# Datenübersicht", "",
              "Erzeugt mit `python make_data_overview.py`. Alles hier ist abgeleitet;",
              "die Quellen liegen unter `cross_sections/` und `chemistry/`.", "",
              "| Gas | Wirkungsquerschnitte | Ratenkoeffizienten |",
              "|-----|---------------------|--------------------|"]
    for gas, e in sorted(bestand.items()):
        zeilen.append(f"| {gas} | {e['querschnitte']} | {e['raten']} |")
    zeilen += ["",
               "Jede CSV-Datei hat eine Kopfzeile mit Prozess und Quelle und danach",
               "zwei Spalten. Bei den Querschnitten sind das Energie in eV und",
               "Wirkungsquerschnitt in m^2, bei den Raten Elektronentemperatur in eV",
               "und Ratenkoeffizient in m^3/s. Analytisch gegebene Raten sind auf dem",
               "Gitter von 0.5 bis 20 eV ausgewertet.", "",
               "Die Abbildungen unter `figures/` werden nicht mitversioniert.", ""]
    (ZIEL / "README.md").write_text("\n".join(zeilen), encoding="utf-8", newline="\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
