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
    energies: list[float] = field(default_factory=list)
    cross_sections: list[float] = field(default_factory=list)


def parse_lxcat(filepath: str | Path) -> list[LXCatBlock]:
    """Parse eine LXCat-Datei und extrahiere alle Bloecke."""
    blocks = []
    lines = Path(filepath).read_text(encoding="utf-8").splitlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # Block-Start erkennen
        if line in ("ELASTIC", "EXCITATION", "IONIZATION", "EFFECTIVE", "ATTACHMENT"):
            block = LXCatBlock(block_type=line, target_line="", threshold_or_mM=0.0)

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
        f.write(f"# Source: LXCat / Biagi database (MagBoltz V8.97)\n")
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


def main():
    if len(sys.argv) < 3:
        print("Usage: python convert_lxcat.py <lxcat_file> <gas_name>")
        print("Example: python convert_lxcat.py cross_sections/tests/xenon_all.txt xenon")
        sys.exit(1)

    lxcat_file = Path(sys.argv[1])
    gas_name = sys.argv[2].lower()

    if not lxcat_file.exists():
        print(f"FEHLER: {lxcat_file} nicht gefunden!")
        sys.exit(1)

    print(f"Parsing: {lxcat_file}")
    blocks = parse_lxcat(lxcat_file)
    print(f"Gefunden: {len(blocks)} Bloecke")

    # Zielverzeichnis
    out_dir = Path("cross_sections") / gas_name
    out_dir.mkdir(parents=True, exist_ok=True)
    exc_dir = out_dir / "excitation"
    exc_dir.mkdir(exist_ok=True)

    metadata = {
        "source": "LXCat / Biagi database (MagBoltz V8.97)",
        "origin_file": str(lxcat_file),
        "retrieved": "2026-03-28",
        "gas": gas_name,
        "processes": [],
    }

    n_elastic = 0
    n_excitation = 0
    n_ionization = 0
    n_skipped = 0

    for block in blocks:
        entry = {
            "type": block.block_type,
            "target": block.target_line,
            "process": block.process_desc,
            "threshold_eV": block.threshold_or_mM if block.block_type != "ELASTIC" else None,
            "mM": block.threshold_or_mM if block.block_type == "ELASTIC" else None,
            "data_points": len(block.energies),
        }

        if block.block_type == "ELASTIC":
            filepath = out_dir / "elastic.csv"
            write_csv(block, filepath)
            entry["file"] = "elastic.csv"
            n_elastic += 1
            print(f"  ELASTIC: {len(block.energies)} Punkte -> {filepath}")

        elif block.block_type == "EXCITATION":
            state = extract_state_name(block.target_line)
            if not state:
                state = f"exc_{n_excitation}"
            base = f"excitation_{sanitize_filename(state)}"
            fname = f"{unique_filename(base)}.csv"
            filepath = exc_dir / fname
            write_csv(block, filepath)
            entry["file"] = f"excitation/{fname}"
            entry["state"] = state
            n_excitation += 1
            print(f"  EXCITATION ({state}): {len(block.energies)} Punkte, "
                  f"E_thr={block.threshold_or_mM} eV -> {filepath}")

        elif block.block_type == "IONIZATION":
            filepath = out_dir / "ionization.csv"
            write_csv(block, filepath)
            entry["file"] = "ionization.csv"
            n_ionization += 1
            print(f"  IONIZATION: {len(block.energies)} Punkte, "
                  f"E_thr={block.threshold_or_mM} eV -> {filepath}")

        else:
            n_skipped += 1
            entry["file"] = None
            entry["note"] = "skipped (unsupported type)"
            print(f"  {block.block_type}: UEBERSPRUNGEN")

        metadata["processes"].append(entry)

    # Metadaten speichern
    meta_path = out_dir / "metadata.json"
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)
    print(f"\nMetadaten: {meta_path}")

    # Zusammenfassung
    print(f"\n{'='*60}")
    print(f"  Zusammenfassung: {gas_name}")
    print(f"{'='*60}")
    print(f"  Elastisch:    {n_elastic}")
    print(f"  Anregung:     {n_excitation}")
    print(f"  Ionisation:   {n_ionization}")
    print(f"  Uebersprungen: {n_skipped}")
    print(f"  Gesamt:       {len(blocks)}")
    print(f"  Zielordner:   {out_dir}")


if __name__ == "__main__":
    main()
