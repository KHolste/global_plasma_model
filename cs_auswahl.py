"""
cs_auswahl.py -- Welcher Querschnittssatz einer Datenbank gilt.

Manche Datenbanken liefern zu einem Gas mehr als einen in sich geschlossenen
Satz oder mehrere Varianten desselben Prozesses. BSR fuehrt fuer Argon zwei
vollstaendige Rechnungen nebeneinander, IAA vier elastische Querschnitte mit
verschiedener Bedeutung. Wer sie unbesehen aufaddiert, zaehlt doppelt.

Welcher Satz gilt, ist eine physikalische Entscheidung und keine Umwandlung.
Sie steht deshalb in einer eigenen Datei `auswahl.json` neben den Daten und
wird von Hand gepflegt. Der Wandler schreibt `metadata.json` und `db_info.json`
bei jedem Lauf neu und wuerde die Entscheidung mitnehmen; `auswahl.json`
ueberlebt.

Aufbau:

    {
      "begruendung": "warum dieser Satz und nicht der andere",
      "elastic":     "elastic_2.csv",
      "ionization":  "ionization_1.csv",
      "excitation":  ["excitation/....csv", ...]
    }

Jeder Schluessel darf fehlen. Fehlt er, gilt die Vorbelegung: `elastic.csv`,
`ionization.csv` und alle Dateien unter `excitation/`.
"""
from __future__ import annotations

import json
from pathlib import Path

DATEINAME = "auswahl.json"


def lade(base_dir: Path) -> dict:
    """Die Auswahl einer Datenbank, leer wenn keine getroffen wurde."""
    pfad = Path(base_dir) / DATEINAME
    if not pfad.exists():
        return {}
    try:
        return json.loads(pfad.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}


def elastic_datei(base_dir: Path) -> Path:
    return Path(base_dir) / lade(base_dir).get("elastic", "elastic.csv")


def ionization_datei(base_dir: Path) -> Path:
    return Path(base_dir) / lade(base_dir).get("ionization", "ionization.csv")


def excitation_dateien(base_dir: Path) -> list[Path]:
    """Die gewaehlten Anregungsdateien, sonst alle vorhandenen."""
    base_dir = Path(base_dir)
    gewaehlt = lade(base_dir).get("excitation")
    if gewaehlt is not None:
        return [base_dir / name for name in gewaehlt]
    exc_dir = base_dir / "excitation"
    return sorted(exc_dir.glob("excitation_*.csv")) if exc_dir.is_dir() else []


def zweitsaetze(base_dir: Path) -> list[str]:
    """Namen der Zweit- und Folgesaetze einer Datenbank.

    Der Wandler haengt an jeden weiteren elastischen oder ionisierenden Block
    eine laufende Nummer. Deren Vorhandensein ist das Zeichen, dass die
    Datenbank mehr als einen Satz mitbringt.
    """
    return sorted(f.name for f in Path(base_dir).glob("*.csv")
                  if f.name.startswith(("elastic_", "ionization_"))
                  and f.stem.rsplit("_", 1)[-1].isdigit())


def ist_entschieden(base_dir: Path) -> bool:
    """Ist klar, welcher Satz gilt?

    Entweder bringt die Datenbank nur einen mit, oder es liegt eine Auswahl
    daneben, die alle drei Prozessarten benennt.
    """
    if not zweitsaetze(base_dir):
        return True
    a = lade(base_dir)
    return all(k in a for k in ("elastic", "ionization", "excitation"))
