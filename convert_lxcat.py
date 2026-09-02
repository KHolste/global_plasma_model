#!/usr/bin/env python3
"""
convert_lxcat.py – Konvertiert LXCat-Querschnittsdaten in das interne CSV-Format.

Liest eine LXCat-Blockdatei und erzeugt pro Prozess eine CSV-Datei
im cross_sections/<gas>/ Verzeichnis.

Ausfuehrung:
    python convert_lxcat.py cross_sections/tests/xenon_all.txt xenon
"""
from __future__ import annotations
import sys, json, re
from pathlib import Path
from dataclasses import dataclass, field


@dataclass
class LXCatBlock:
    """Ein einzelner Stossprozess aus einer LXCat-Datei."""
    block_type: str         # ELASTIC, EXCITATION, IONIZATION
    target_line: str        # z.B. "Xe -> Xe(1S5)(8.315eV)"
    threshold_or_mM: float  # Schwellenenergie [eV] oder m/M
    process_desc: str = ""
    comment: str = ""
    updated: str = ""
    species: str = ""
    param: str = ""
    database: str = ""       # z.B. "Hayashi database"
    permlink: str = ""       # z.B. "www.lxcat.net/Hayashi"
    db_comment: str = ""     # Quellenangabe der Datenbank
    energies: list[float] = field(default_factory=list)
    cross_sections: list[float] = field(default_factory=list)


def parse_lxcat(filepath: str | Path) -> list[LXCatBlock]:
    """Parse eine LXCat-Datei und extrahiere alle Bloecke."""
    blocks = []
    lines = Path(filepath).read_text(encoding="utf-8").splitlines()

    # Eine LXCat-Datei kann mehrere Datenbanken enthalten. Jeder Block gehoert
    # zu der Datenbank, deren Kopf zuletzt kam.
    db_name = db_permlink = db_comment = ""
    db_comment_offen = False

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("DATABASE:"):
            db_name = line[9:].strip()
            db_permlink = db_comment = ""
            db_comment_offen = False
        elif line.startswith("PERMLINK:"):
            db_permlink = line[9:].strip()
        elif line.startswith("COMMENT:") and not db_comment:
            db_comment = line[8:].strip()
            db_comment_offen = True
        elif db_comment_offen:
            if line and not line.startswith(("*", "x", "-")) and ":" not in line[:12]:
                db_comment += " " + line
            else:
                db_comment_offen = False

        # Block-Start erkennen
        if line in ("ELASTIC", "EXCITATION", "IONIZATION", "EFFECTIVE", "ATTACHMENT"):
            db_comment_offen = False
            block = LXCatBlock(block_type=line, target_line="", threshold_or_mM=0.0,
                               database=db_name, permlink=db_permlink,
                               db_comment=db_comment)

            # Zeile 2: Target
            i += 1
            if i < len(lines):
                block.target_line = lines[i].strip()

            # Zeile 3: Schwellenenergie oder m/M
            i += 1
            if i < len(lines):
                try:
                    block.threshold_or_mM = float(lines[i].strip())
                except ValueError:
                    pass

            # Kommentarzeilen bis zur Datentabelle
            i += 1
            while i < len(lines):
                cl = lines[i].strip()
                if cl.startswith("-----"):
                    break
                if cl.startswith("SPECIES:"):
                    block.species = cl[8:].strip()
                elif cl.startswith("PROCESS:"):
                    block.process_desc = cl[8:].strip()
                elif cl.startswith("COMMENT:"):
                    block.comment = cl[8:].strip()
                elif cl.startswith("UPDATED:"):
                    block.updated = cl[8:].strip()
                elif cl.startswith("PARAM.:"):
                    block.param = cl[7:].strip()
                i += 1

            # Datentabelle lesen
            i += 1  # Ueberspringe "-----"
            while i < len(lines):
                dl = lines[i].strip()
                if dl.startswith("-----"):
                    break
                parts = dl.split()
                if len(parts) >= 2:
                    try:
                        block.energies.append(float(parts[0]))
                        block.cross_sections.append(float(parts[1]))
                    except ValueError:
                        pass
                i += 1

            blocks.append(block)
        i += 1

    return blocks


def sanitize_filename(name: str) -> str:
    """Erzeuge einen sauberen Dateinamen aus einem Prozessnamen."""
    # Entferne Sonderzeichen, behalte Buchstaben/Zahlen/Unterstriche
    name = name.replace("->", "_to_").replace("<->", "_rev_")
    name = re.sub(r"[^a-zA-Z0-9_]", "_", name)
    name = re.sub(r"_+", "_", name).strip("_")
    return name.lower()


def extract_state_name(target_line: str) -> str:
    """Extrahiere den Zustandsnamen aus der Target-Zeile."""
    parts = re.findall(r"\(([^)]+)\)", target_line)
    for p in parts:
        if "eV" not in p:
            return p
    return ""


# Zaehler fuer eindeutige Dateinamen bei Kollisionen
_filename_counter: dict[str, int] = {}

def unique_filename(base: str) -> str:
    """Erzeuge einen eindeutigen Dateinamen falls Kollision."""
    if base not in _filename_counter:
        _filename_counter[base] = 0
        return base
    _filename_counter[base] += 1
    return f"{base}_{_filename_counter[base]}"


def write_csv(block: LXCatBlock, filepath: Path):
    """Schreibe einen Block als CSV."""
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(f"# {block.process_desc}\n")
        f.write(f"# Source: LXCat, {block.database or 'unbekannte Datenbank'}\n")
        f.write(f"# Target: {block.target_line}\n")
        if block.block_type in ("EXCITATION", "IONIZATION"):
            f.write(f"# Threshold: {block.threshold_or_mM} eV\n")
        elif block.block_type == "ELASTIC":
            f.write(f"# m/M: {block.threshold_or_mM}\n")
        if block.comment:
            f.write(f"# Comment: {block.comment}\n")
        if block.updated:
            f.write(f"# Updated: {block.updated}\n")
        f.write(f"# Data points: {len(block.energies)}\n")
        f.write("energy_eV,cross_section_m2\n")
        for e, cs in zip(block.energies, block.cross_sections):
            f.write(f"{e},{cs}\n")


def datenbank_kurz(name: str) -> str:
    """Aus "Hayashi database" wird "hayashi"."""
    wort = name.replace("database", "").strip().split()
    return sanitize_filename(wort[0].lower()) if wort else "unbekannt"


def erzeugt_am(lxcat_file: Path) -> str:
    """Das Erzeugungsdatum aus dem Dateikopf, sonst leer."""
    for zeile in lxcat_file.read_text(encoding="utf-8").splitlines()[:5]:
        if zeile.startswith("Generated on"):
            rest = zeile.replace("Generated on", "").strip()
            return rest.split(".")[0].strip()
    return ""


def schreibe_datenbank(blocks, out_dir: Path, gas_name: str,
                       lxcat_file: Path, datum: str) -> dict:
    """Alle Bloecke einer Datenbank in ein Verzeichnis schreiben."""
    global _filename_counter
    _filename_counter = {}

    out_dir.mkdir(parents=True, exist_ok=True)
    exc_dir = out_dir / "excitation"
    exc_dir.mkdir(exist_ok=True)

    db_name = blocks[0].database or out_dir.name
    metadata = {
        "source": f"LXCat, {db_name}",
        "origin_file": str(lxcat_file),
        "retrieved": datum,
        "gas": gas_name,
        "processes": [],
    }

    zahl = {"ELASTIC": 0, "EXCITATION": 0, "IONIZATION": 0, "EFFECTIVE": 0, "sonst": 0}

    for block in blocks:
        entry = {
            "type": block.block_type,
            "target": block.target_line,
            "process": block.process_desc,
            "threshold_eV": block.threshold_or_mM if block.block_type not in ("ELASTIC", "EFFECTIVE") else None,
            "mM": block.threshold_or_mM if block.block_type in ("ELASTIC", "EFFECTIVE") else None,
            "data_points": len(block.energies),
        }

        if block.block_type == "ELASTIC":
            write_csv(block, out_dir / "elastic.csv")
            entry["file"] = "elastic.csv"
            zahl["ELASTIC"] += 1

        elif block.block_type == "EFFECTIVE":
            # Der effektive Querschnitt ist der gesamte Impulsuebertrag, also
            # elastisch plus inelastisch. Er wird abgelegt, damit er sichtbar
            # ist, aber nicht als elastischer Querschnitt verwendet.
            write_csv(block, out_dir / "effective.csv")
            entry["file"] = "effective.csv"
            entry["note"] = ("Gesamter Impulsuebertrag, nicht der elastische. "
                             "Nicht als elastischer Querschnitt verwendbar.")
            zahl["EFFECTIVE"] += 1

        elif block.block_type == "EXCITATION":
            state = extract_state_name(block.target_line) or f"exc_{zahl['EXCITATION']}"
            fname = f"{unique_filename('excitation_' + sanitize_filename(state))}.csv"
            write_csv(block, exc_dir / fname)
            entry["file"] = f"excitation/{fname}"
            entry["state"] = state
            zahl["EXCITATION"] += 1

        elif block.block_type == "IONIZATION":
            name = "ionization.csv" if zahl["IONIZATION"] == 0 else f"ionization_{zahl['IONIZATION']}.csv"
            write_csv(block, out_dir / name)
            entry["file"] = name
            zahl["IONIZATION"] += 1

        else:
            entry["file"] = None
            entry["note"] = "nicht unterstuetzter Typ"
            zahl["sonst"] += 1

        metadata["processes"].append(entry)

    with open(out_dir / "metadata.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)

    db_info = {
        "key": out_dir.name,
        "display_name": db_name,
        "source": f"{db_name}, via LXCat",
        "reference": blocks[0].permlink or "www.lxcat.net",
        "comment": blocks[0].db_comment,
        "retrieved": datum,
        "gas": gas_name,
        "processes": {
            "elastic": zahl["ELASTIC"] > 0,
            "effective": zahl["EFFECTIVE"] > 0,
            "ionization": zahl["IONIZATION"] > 0,
            "excitation": zahl["EXCITATION"] > 0,
            "excitation_count": zahl["EXCITATION"],
        },
    }
    with open(out_dir / "db_info.json", "w", encoding="utf-8") as f:
        json.dump(db_info, f, indent=2, ensure_ascii=False)

    print(f"  {db_name:<22} -> {out_dir}")
    print(f"      elastisch {zahl['ELASTIC']}, effektiv {zahl['EFFECTIVE']}, "
          f"Anregung {zahl['EXCITATION']}, Ionisation {zahl['IONIZATION']}"
          + (f", uebersprungen {zahl['sonst']}" if zahl["sonst"] else ""))
    return zahl


def main():
    if len(sys.argv) < 3:
        print("Aufruf: python convert_lxcat.py <lxcat_datei> <zielordner>")
        print("Enthaelt die Datei mehrere Datenbanken, entsteht je Datenbank")
        print("ein Unterordner:")
        print("  python convert_lxcat.py argon.txt cross_sections/argon")
        print("  -> cross_sections/argon/hayashi/, cross_sections/argon/phelps/")
        sys.exit(1)

    lxcat_file = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    if not lxcat_file.exists():
        print(f"FEHLER: {lxcat_file} nicht gefunden!")
        sys.exit(1)

    blocks = parse_lxcat(lxcat_file)
    datum = erzeugt_am(lxcat_file)
    gas_name = out_dir.name if out_dir.parent.name == "cross_sections" else out_dir.parent.name

    nach_db = {}
    for b in blocks:
        nach_db.setdefault(b.database or "", []).append(b)

    print(f"{lxcat_file}: {len(blocks)} Bloecke aus {len(nach_db)} Datenbank(en)"
          + (f", erzeugt am {datum}" if datum else ""))

    for db_name, db_blocks in nach_db.items():
        ziel = out_dir if len(nach_db) == 1 else out_dir / datenbank_kurz(db_name)
        schreibe_datenbank(db_blocks, ziel, gas_name, lxcat_file, datum)


if __name__ == "__main__":
    main()
