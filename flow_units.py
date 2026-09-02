"""
flow_units.py -- Massenfluss-Einheitenumrechnung fuer Xenon und Iod.

Source of Truth: Q0 in sccm (gasunabhaengig, Volumenstrom bei STP).
Darstellung: sccm oder mg/s (gasabhaengig ueber Feedstock-Masse).

Feedstock-Basis:
  Xenon: Xe (131.293 u)
  Iod:   I2 (253.809 u) -- molekular, NICHT atomar

Umrechnung:
  mg/s = sccm * SCCM_TO_PPS * M_feedstock_kg * 1e6
  sccm = mg_per_s / (SCCM_TO_PPS * M_feedstock_kg * 1e6)
"""
from __future__ import annotations

SCCM_TO_PPS = 4.477962312e17  # sccm -> particles/s (gasunabhaengig, STP)

# Feedstock-Massen [kg] -- Basis fuer mg/s-Umrechnung
FEEDSTOCK_MASS = {
    "xenon":  2.1801711e-25,   # Xe: 131.293 u
    "iodine": 4.2143422e-25,   # I2: 253.809 u (MOLEKULAR, nicht atomar!)
    "krypton": 1.3914984e-25,  # Kr: 83.798 u
    "argon":  6.6335209e-26,   # Ar: 39.948 u
}

# Feedstock-Spezies-IDs
FEEDSTOCK_SPECIES = {
    "xenon": "Xe",
    "iodine": "I2",   # Molekulares Iod (Einspeisespezies)
    "krypton": "Kr",
    "argon": "Ar",
}


def feedstock_mass_kg(gas: str) -> float:
    """Feedstock-Teilchenmasse in kg. Fuer Iod: I2 (nicht I!)."""
    return FEEDSTOCK_MASS.get(gas, FEEDSTOCK_MASS["xenon"])


def sccm_to_mg_per_s(sccm: float, gas: str) -> float:
    """Konvertiere sccm -> mg/s fuer gegebenes Gas."""
    M = feedstock_mass_kg(gas)
    return sccm * SCCM_TO_PPS * M * 1e6  # kg -> mg = 1e6


def mg_per_s_to_sccm(mg_per_s: float, gas: str) -> float:
    """Konvertiere mg/s -> sccm fuer gegebenes Gas."""
    M = feedstock_mass_kg(gas)
    return mg_per_s / (SCCM_TO_PPS * M * 1e6)


def sccm_to_pps(sccm: float) -> float:
    """sccm -> particles/s (gasunabhaengig)."""
    return sccm * SCCM_TO_PPS


def pps_to_sccm(pps: float) -> float:
    """particles/s -> sccm (gasunabhaengig)."""
    return pps / SCCM_TO_PPS
