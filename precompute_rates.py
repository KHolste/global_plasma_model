#!/usr/bin/env python3
"""
precompute_rates.py – Berechnet tabellierte Ratenkoeffizienten aus
Wirkungsquerschnitten und speichert sie als Lookup-Tabelle fuer den C++-Solver.

Erzeugt:
  cross_sections/xenon/kex_table.csv
    Spalten: Te_eV, Kex_total_m3s, Pexc_coeff_Wm3
    Pexc_coeff = Summe_i [Kex_i * Eexc_i] (in W*m^3 wenn mit n*ng multipliziert)

Ausfuehrung:
    python precompute_rates.py
"""
from __future__ import annotations
import numpy as np
from pathlib import Path
from rate_coefficients import CrossSectionData, compute_rate_coefficient, kex_legacy

E_CH = 1.602176487e-19  # J/eV

def main():
    exc_dir = Path("cross_sections/xenon/excitation")
    out_path = Path("cross_sections/xenon/kex_table.csv")

    # Alle Anregungsprozesse laden
    exc_files = sorted(exc_dir.glob("excitation_*.csv"))
    processes = []
    for f in exc_files:
        cs = CrossSectionData.from_csv(f)
        if cs and len(cs.energy_eV) >= 2:
            processes.append((f.stem, cs))

    print(f"Geladene Anregungsprozesse: {len(processes)}")

    # Te-Gitter: 0.5 bis 20 eV in 0.05 eV Schritten
    Te_grid = np.arange(0.5, 20.05, 0.05)

    print(f"Te-Gitter: {Te_grid[0]:.2f} bis {Te_grid[-1]:.2f} eV ({len(Te_grid)} Punkte)")
    print("Berechne Ratenkoeffizienten...")

    rows = []
    for Te in Te_grid:
        Kex_total = 0.0
        Pexc_coeff = 0.0  # Summe_i [Kex_i * Eexc_i_J]

        for name, cs in processes:
            K_i = compute_rate_coefficient(cs, float(Te))
            if K_i > 0:
                Kex_total += K_i
                Pexc_coeff += K_i * cs.threshold_eV * E_CH  # Eexc in Joule

        rows.append((float(Te), Kex_total, Pexc_coeff))

    # CSV schreiben
    with open(out_path, "w") as f:
        f.write("# Vorberechnete Xenon-Anregungsraten (Biagi/LXCat, 48 Prozesse)\n")
        f.write("# Kex_total = Summe_i Kex_i(Te)\n")
        f.write("# Pexc_coeff = Summe_i [Kex_i(Te) * Eexc_i]  (multipliziere mit n*ng fuer W/m^3)\n")
        f.write("# Te in eV, Kex in m^3/s, Pexc_coeff in J*m^3/s\n")
        f.write("Te_eV,Kex_total_m3s,Pexc_coeff_Jm3s\n")
        for Te, Kex, Pexc in rows:
            f.write(f"{Te:.3f},{Kex:.6e},{Pexc:.6e}\n")

    print(f"Geschrieben: {out_path} ({len(rows)} Zeilen)")

    # Verifikation
    print(f"\nVerifikation gegen Legacy-Fit:")
    print(f"{'Te[eV]':>8} {'Kex_tab':>12} {'Kex_leg':>12} {'Ratio':>8} {'Pexc_coeff':>14}")
    print("-" * 60)
    for Te, Kex, Pexc in rows:
        if Te in [2.0, 3.0, 4.0, 5.0, 8.0, 10.0]:
            kl = kex_legacy(Te)
            print(f"{Te:8.1f} {Kex:12.4e} {kl:12.4e} {Kex/kl:8.3f} {Pexc:14.4e}")


if __name__ == "__main__":
    main()
