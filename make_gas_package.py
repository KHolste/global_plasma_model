#!/usr/bin/env python3
"""
make_gas_package.py -- Erzeugt aus Querschnittsdaten ein Chemiepaket.

Der fest verdrahtete Weg fasst die Anregung zu einem einzigen Term zusammen.
Die Querschnittsdaten enthalten aber viele einzelne Anregungsprozesse mit je
eigener Schwelle. Dieses Skript integriert jeden Prozess einzeln über eine
Maxwell-Verteilung und schreibt daraus ein Chemiepaket, in dem jeder Prozess
eine eigene Reaktion mit eigener Ratentabelle ist.

Erzeugt:
    chemistry/<gas>_<datenbank>/chemistry.json
    chemistry/<gas>_<datenbank>/rates/*.csv

Zur Kontrolle wird die Summe der Einzelprozesse gegen die bereits vorhandene
zusammengefasste Tabelle gestellt: sowohl die Summe der Ratenkoeffizienten als
auch die Summe der Energieverluste müssen übereinstimmen.

Aufruf:
    python make_gas_package.py
    python make_gas_package.py --gas argon --db hayashi
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import cs_auswahl  # noqa: E402
from rate_coefficients import CrossSectionData, compute_rate_coefficient  # noqa: E402

E_CH = 1.602176487e-19
SIGMA_I = 1e-18

#: Stoffwerte je Gas: Symbol, Masse [kg], Waermeleitfaehigkeit [W/(m K)].
#: Dieselben Werte wie im Rechenkern.
GASE = {
    "xenon":   ("Xe", 2.1801711e-25, 0.0057),
    "krypton": ("Kr", 1.3914984e-25, 0.0094),
    "argon":   ("Ar", 6.6335209e-26, 0.0177),
}

#: Dasselbe Gitter wie die vorhandenen zusammengefassten Tabellen.
TE_GRID = np.arange(0.5, 20.05, 0.05)


def rate_table(cs: CrossSectionData) -> np.ndarray:
    return np.array([compute_rate_coefficient(cs, float(Te)) for Te in TE_GRID])


def schreibe_tabelle(pfad: Path, kopf: str, werte: np.ndarray) -> None:
    pfad.parent.mkdir(parents=True, exist_ok=True)
    with open(pfad, "w", encoding="utf-8", newline="\n") as f:
        f.write(f"# {kopf}\n")
        f.write("Te_eV,K_m3s\n")
        for Te, K in zip(TE_GRID, werte):
            f.write(f"{Te:.3f},{K:.6e}\n")


def kennung(symbol: str, state: str) -> str:
    """Reaktionskennung aus dem Zustandsnamen, z.B. 1S5 -> exc_Xe_1s5."""
    return f"exc_{symbol}_" + state.lower().replace(" ", "").replace("/", "_")


def baue(gas: str, cs_dir: Path, out_dir: Path) -> int:
    symbol, masse, kappa = GASE[gas]
    ion = f"{symbol}+"
    meta = json.loads((cs_dir / "metadata.json").read_text(encoding="utf-8"))
    prozesse = meta["processes"]

    # Welcher Satz gilt, wenn die Datenbank mehrere mitbringt, steht in
    # auswahl.json neben den Daten. Ohne Entscheidung wird nichts gebaut --
    # das Paket wuerde sonst zwei Rechnungen desselben Gases vermengen.
    if not cs_auswahl.ist_entschieden(cs_dir):
        print(f"FEHLER: {cs_dir} bringt mehrere Saetze mit "
              f"({', '.join(cs_auswahl.zweitsaetze(cs_dir))}).",
              file=sys.stderr)
        print(f"Welcher gilt, gehoert nach {cs_dir / 'auswahl.json'}.",
              file=sys.stderr)
        return 1

    el_datei = cs_auswahl.elastic_datei(cs_dir).name
    iz_datei = cs_auswahl.ionization_datei(cs_dir).name
    exc_gewaehlt = {p.as_posix().split("/")[-1]
                    for p in cs_auswahl.excitation_dateien(cs_dir)}

    elastisch = next((p for p in prozesse if p["file"] == el_datei), None)
    ionisation = next((p for p in prozesse if p["file"] == iz_datei), None)
    anregungen = [p for p in prozesse if p["type"] == "EXCITATION"
                  and Path(p["file"]).name in exc_gewaehlt]
    if not elastisch or not ionisation:
        print("FEHLER: elastischer oder ionisierender Prozess fehlt.", file=sys.stderr)
        return 1

    print(f"Quelle: {cs_dir}")
    print(f"  1 elastisch ({el_datei}), 1 ionisierend ({iz_datei}), "
          f"{len(anregungen)} Anregungen")
    print(f"  Te-Gitter: {TE_GRID[0]:.2f} bis {TE_GRID[-1]:.2f} eV, {len(TE_GRID)} Punkte")

    rates_dir = out_dir / "rates"
    reaktionen: list[dict] = []

    # ── Elastischer Stoss ────────────────────────────────────────
    cs_el = CrossSectionData.from_csv(cs_dir / elastisch["file"])
    if cs_el is None:
        print("FEHLER: elastische Querschnitte nicht lesbar.", file=sys.stderr)
        return 1
    K_el = rate_table(cs_el)
    schreibe_tabelle(rates_dir / f"Kel_{symbol}.csv",
                     f"Elastischer Impulsuebertrag, {cs_dir.name}", K_el)
    reaktionen.append({
        "id": f"el_{symbol}",
        "name": f"Elastisch: e + {symbol} -> e + {symbol}",
        "type": "elastic",
        "reactants": {"e": 1, symbol: 1},
        "products": {"e": 1, symbol: 1},
        "rate": {"model": "tabulated", "file": f"rates/Kel_{symbol}.csv"},
        "energy_eV": 0.0,
        "is_electron_impact": True,
        "elastic_heating": True,
        "nu_m": True,
    })

    # ── Ionisation ───────────────────────────────────────────────
    cs_iz = CrossSectionData.from_csv(cs_dir / ionisation["file"])
    if cs_iz is None:
        print("FEHLER: Ionisationsquerschnitte nicht lesbar.", file=sys.stderr)
        return 1
    K_iz = rate_table(cs_iz)
    schreibe_tabelle(rates_dir / f"Kiz_{symbol}.csv",
                     f"Einfachionisation, {cs_dir.name}", K_iz)
    E_iz = float(cs_iz.threshold_eV or ionisation["threshold_eV"])
    reaktionen.append({
        "id": f"iz_{symbol}",
        "name": f"Ionisation: e + {symbol} -> 2e + {ion}",
        "type": "ionization",
        "reactants": {"e": 1, symbol: 1},
        "products": {"e": 2, ion: 1},
        "rate": {"model": "tabulated", "file": f"rates/Kiz_{symbol}.csv"},
        "energy_eV": E_iz,
        "is_electron_impact": True,
    })

    # ── Anregungen, jede fuer sich ───────────────────────────────
    summe_K = np.zeros_like(TE_GRID)
    summe_P = np.zeros_like(TE_GRID)
    uebersprungen = []
    for p in anregungen:
        cs = CrossSectionData.from_csv(cs_dir / p["file"])
        if cs is None or len(cs.energy_eV) < 2:
            uebersprungen.append(p.get("state", p["file"]))
            continue
        K = rate_table(cs)
        if not np.any(K > 0):
            uebersprungen.append(p.get("state", p["file"]))
            continue
        schwelle = float(cs.threshold_eV or p["threshold_eV"])
        state = p.get("state") or Path(p["file"]).stem
        rid = kennung(symbol, state)
        datei = f"rates/K_{rid}.csv"
        schreibe_tabelle(out_dir / datei,
                         f"Anregung {symbol}({state}), Schwelle {schwelle} eV, {cs_dir.name}", K)
        reaktionen.append({
            "id": rid,
            "name": f"Anregung: e + {symbol} -> e + {symbol}({state}) bei {schwelle} eV",
            "type": "excitation",
            "reactants": {"e": 1, symbol: 1},
            "products": {"e": 1, symbol: 1},
            "rate": {"model": "tabulated", "file": datei},
            "energy_eV": schwelle,
            "is_electron_impact": True,
        })
        summe_K += K
        summe_P += K * schwelle * E_CH

    print(f"  {len(reaktionen) - 2} Anregungen uebernommen"
          + (f", {len(uebersprungen)} ohne Daten uebersprungen" if uebersprungen else ""))

    # ── Chemiepaket schreiben ────────────────────────────────────
    paket = {
        "name": f"{gas.capitalize()} ({cs_dir.name})",
        "description": (
            f"{gas.capitalize()} aus {meta.get('source', 'LXCat')}. Elastischer Stoss, Einfachionisation "
            f"und {len(reaktionen) - 2} einzelne Anregungsprozesse mit je eigener Schwelle "
            f"und eigener Ratentabelle, integriert ueber eine Maxwell-Verteilung."
        ),
        "wall_temperature_K": 293.0,
        "sigma_i": SIGMA_I,
        "species": [
            {"id": "e", "name": "Elektron", "type": "electron",
             "mass_kg": 9.10938215e-31, "charge": -1},
            {"id": symbol, "name": gas.capitalize(), "type": "neutral_atom",
             "mass_kg": masse, "charge": 0, "is_feedstock": True,
             "thermal_conductivity": kappa},
            {"id": ion, "name": ion, "type": "positive_ion",
             "mass_kg": masse, "charge": 1, "is_beam_extracted": True,
             "wall_products": {symbol: 1}},
        ],
        "reactions": reaktionen,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    ziel = out_dir / "chemistry.json"
    with open(ziel, "w", encoding="utf-8", newline="\n") as f:
        json.dump(paket, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print(f"  geschrieben: {ziel}")

    # ── Kontrolle gegen die zusammengefasste Tabelle ─────────────
    kex = cs_dir / "kex_table.csv"
    if not kex.exists():
        print("  Hinweis: kex_table.csv fehlt, keine Gegenprobe moeglich.")
        return 0

    zeilen = [z for z in kex.read_text(encoding="utf-8").splitlines()
              if z and not z.startswith("#") and not z.startswith("Te")]
    tab = np.array([[float(x) for x in z.split(",")] for z in zeilen])
    Te_t, K_t, P_t = tab[:, 0], tab[:, 1], tab[:, 2]
    K_ref = np.interp(TE_GRID, Te_t, K_t)
    P_ref = np.interp(TE_GRID, Te_t, P_t)
    maske = K_ref > 0
    d_K = np.max(np.abs(summe_K[maske] - K_ref[maske]) / K_ref[maske])
    d_P = np.max(np.abs(summe_P[maske] - P_ref[maske]) / P_ref[maske])
    print(f"  Gegenprobe gegen kex_table.csv:")
    print(f"    groesste Abweichung der Ratensumme:        {d_K:.2e}")
    print(f"    groesste Abweichung der Energieverluste:   {d_P:.2e}")
    if max(d_K, d_P) > 1e-6:
        print("    WARNUNG: Einzelprozesse und Sammelt abelle weichen ab.", file=sys.stderr)
        return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Chemiepaket aus Querschnittsdaten erzeugen.")
    ap.add_argument("--gas", default="xenon", help="Gasverzeichnis unter cross_sections/")
    ap.add_argument("--db", default="biagi", help="Datenbankverzeichnis")
    ap.add_argument("--out", default=None, help="Zielverzeichnis des Pakets")
    args = ap.parse_args(argv)

    cs_dir = SCRIPT_DIR / "cross_sections" / args.gas / args.db
    if not cs_dir.is_dir():
        print(f"FEHLER: {cs_dir} nicht gefunden.", file=sys.stderr)
        return 1
    out_dir = Path(args.out) if args.out else SCRIPT_DIR / "chemistry" / f"{args.gas}_{args.db}"
    if args.gas not in GASE:
        print(f"FEHLER: fuer {args.gas} sind keine Stoffwerte hinterlegt.", file=sys.stderr)
        return 1
    return baue(args.gas, cs_dir, out_dir)


if __name__ == "__main__":
    sys.exit(main())
