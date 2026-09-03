#!/usr/bin/env python3
"""
precompute_rates.py -- Berechnet tabellierte Ratenkoeffizienten aus
Wirkungsquerschnitten und speichert sie als Lookup-Tabelle fuer den C++-Solver.

Erzeugt:
  <base_dir>/kiz_table.csv
  <base_dir>/kel_table.csv
  <base_dir>/kex_table.csv

Ausfuehrung:
    python precompute_rates.py                                    # Default: xenon/biagi
    python precompute_rates.py cross_sections/xenon/biagi
    python precompute_rates.py cross_sections/xenon/hayashi
"""
from __future__ import annotations
import sys
import numpy as np
from pathlib import Path
import cs_auswahl
from rate_coefficients import CrossSectionData, compute_rate_coefficient

E_CH = 1.602176487e-19  # J/eV


def precompute(base_dir: Path, trotzdem: bool = False):
    print(f"Basis: {base_dir}")

    # Bringt die Datenbank mehrere Saetze mit, wuerden die Sammeltabellen sie
    # aufaddieren und damit doppelt zaehlen. Welcher gilt, steht in
    # auswahl.json; ohne sie wird hier nichts gerechnet.
    if not cs_auswahl.ist_entschieden(base_dir) and not trotzdem:
        print("  ABBRUCH: die Datenbank bringt mehr als einen elastischen oder")
        print("  ionisierenden Satz mit:")
        for name in cs_auswahl.zweitsaetze(base_dir):
            print(f"    {name}")
        print("  Die Sammeltabellen wuerden die Saetze aufaddieren und damit")
        print("  doppelt zaehlen. Welcher Satz gilt, gehoert nach auswahl.json")
        print("  neben die Daten -- oder mit --trotzdem rechnen.")
        return False

    if cs_auswahl.lade(base_dir):
        print(f"  Auswahl: {cs_auswahl.elastic_datei(base_dir).name}, "
              f"{cs_auswahl.ionization_datei(base_dir).name}, "
              f"{len(cs_auswahl.excitation_dateien(base_dir))} Anregungen")

    Te_grid = np.arange(0.5, 20.05, 0.05)
    print(f"Te-Gitter: {Te_grid[0]:.2f} - {Te_grid[-1]:.2f} eV ({len(Te_grid)} Punkte)")

    # Kiz
    cs_iz = CrossSectionData.from_csv(cs_auswahl.ionization_datei(base_dir))
    if cs_iz and len(cs_iz.energy_eV) >= 2:
        out = base_dir / "kiz_table.csv"
        with open(out, "w") as f:
            f.write(f"# Ionisation rate ({base_dir.name})\nTe_eV,Kiz_m3s\n")
            for Te in Te_grid:
                f.write(f"{Te:.3f},{compute_rate_coefficient(cs_iz, float(Te)):.6e}\n")
        print(f"  kiz_table.csv: {len(Te_grid)} Zeilen")
    else:
        print(f"  WARNUNG: ionization.csv nicht gefunden oder leer")

    # Kel
    cs_el = CrossSectionData.from_csv(cs_auswahl.elastic_datei(base_dir))
    if cs_el and len(cs_el.energy_eV) >= 2:
        out = base_dir / "kel_table.csv"
        with open(out, "w") as f:
            f.write(f"# Elastic rate ({base_dir.name})\nTe_eV,Kel_m3s\n")
            for Te in Te_grid:
                f.write(f"{Te:.3f},{compute_rate_coefficient(cs_el, float(Te)):.6e}\n")
        print(f"  kel_table.csv: {len(Te_grid)} Zeilen")
    else:
        print(f"  WARNUNG: elastic.csv nicht gefunden oder leer")

    # Kex (Summe aller Anregungsprozesse)
    exc_files = cs_auswahl.excitation_dateien(base_dir)
    if exc_files:
        procs = []
        for fp in exc_files:
            cs = CrossSectionData.from_csv(fp)
            if cs and len(cs.energy_eV) >= 2:
                procs.append((fp.stem, cs))
        print(f"  {len(procs)} Anregungsprozesse geladen")

        out = base_dir / "kex_table.csv"
        with open(out, "w") as f:
            f.write(f"# Excitation rates ({base_dir.name}, {len(procs)} Prozesse)\n")
            f.write("Te_eV,Kex_total_m3s,Pexc_coeff_Jm3s\n")
            for Te in Te_grid:
                Kex_tot = 0.0
                Pexc = 0.0
                for _, cs in procs:
                    K = compute_rate_coefficient(cs, float(Te))
                    if K > 0:
                        Kex_tot += K
                        Pexc += K * cs.threshold_eV * E_CH
                f.write(f"{Te:.3f},{Kex_tot:.6e},{Pexc:.6e}\n")
        print(f"  kex_table.csv: {len(Te_grid)} Zeilen")
    else:
        print(f"  WARNUNG: keine Anregungsquerschnitte gefunden")

    print("Fertig.")
    return True


def main():
    argumente = [a for a in sys.argv[1:] if a != "--trotzdem"]
    trotzdem = "--trotzdem" in sys.argv
    base = Path(argumente[0]) if argumente else Path("cross_sections/xenon/biagi")

    if not base.is_dir():
        print(f"FEHLER: {base} nicht gefunden!")
        sys.exit(1)

    if not precompute(base, trotzdem):
        sys.exit(2)


if __name__ == "__main__":
    main()
