#!/usr/bin/env python3
"""
rate_coefficients.py – Berechnung von Ratenkoeffizienten aus tabellierten
Wirkungsquerschnitten fuer das Global Xenon Plasma Model.

Laedt Querschnittsdaten aus CSV-Dateien und berechnet daraus temperatur-
abhaengige Ratenkoeffizienten ueber Maxwell-Boltzmann-Integration.

Unterstuetzt:
  - Ionisation (Kiz)
  - Anregung (Kex)
  - Elastische Stoesse (Kel)

Physik:
  Ratenkoeffizient K(Te) = <sigma * v> = Integral ueber Maxwell-Verteilung:
  K(Te) = integral_0^inf  sigma(E) * v(E) * f_MB(E, Te) dE

  mit:
    v(E) = sqrt(2E/me)                   Elektronengeschwindigkeit
    f_MB(E, Te) = 2/sqrt(pi) * (1/kTe)^(3/2) * sqrt(E) * exp(-E/kTe)
                                          Maxwell-Boltzmann Energieverteilung
    kTe = kB * Te [J]                    Thermische Energie

Einheiten:
  - Energie E:            eV (in CSV) -> J (intern)
  - Querschnitt sigma:    m^2
  - Temperatur Te:        eV (Eingabe) -> K (intern via conv=11604.5)
  - Ratenkoeffizient K:   m^3/s
"""
from __future__ import annotations

import math
import os
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

# Physikalische Konstanten
ME = 9.10938215e-31      # Elektronenmasse [kg]
E_CH = 1.602176487e-19   # Elementarladung [C = J/eV]
KB = 1.3806504e-23       # Boltzmann [J/K]
CONV = 11604.5250061657  # eV -> K


@dataclass
class CrossSectionData:
    """Tabellierte Wirkungsquerschnitte fuer einen Stossprozess."""
    energy_eV: np.ndarray      # Energiepunkte [eV]
    cross_section_m2: np.ndarray  # Querschnitte [m^2]
    process_name: str = ""
    threshold_eV: float = 0.0
    source: str = ""

    @staticmethod
    def from_csv(path: str | Path) -> Optional["CrossSectionData"]:
        """Lade Querschnittsdaten aus einer CSV-Datei."""
        path = Path(path)
        if not path.exists():
            return None

        energies = []
        cross_sections = []
        process_name = ""
        threshold = 0.0
        source = ""

        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    # Metadaten aus Kommentaren
                    if "Threshold:" in line:
                        try:
                            threshold = float(line.split("Threshold:")[1].strip().split()[0])
                        except (ValueError, IndexError):
                            pass
                    if "Process:" in line:
                        process_name = line.split("Process:")[1].strip()
                    if "Source:" in line:
                        source = line.split("Source:")[1].strip()
                    continue
                if line.startswith("energy_eV"):
                    continue  # Header

                parts = line.split(",")
                if len(parts) >= 2:
                    try:
                        energies.append(float(parts[0]))
                        cross_sections.append(float(parts[1]))
                    except ValueError:
                        continue

        if not energies:
            return None

        return CrossSectionData(
            energy_eV=np.array(energies),
            cross_section_m2=np.array(cross_sections),
            process_name=process_name,
            threshold_eV=threshold,
            source=source,
        )


def compute_rate_coefficient(cs: CrossSectionData, Te_eV: float) -> float:
    """Berechne den Ratenkoeffizienten K(Te) durch Maxwell-Boltzmann-Integration.

    K(Te) = integral sigma(E) * v(E) * f_MB(E, Te) dE

    Args:
        cs: Tabellierte Querschnittsdaten
        Te_eV: Elektronentemperatur in eV

    Returns:
        Ratenkoeffizient in m^3/s
    """
    if Te_eV <= 0 or len(cs.energy_eV) < 2:
        return 0.0

    kTe_J = KB * Te_eV * CONV  # kB * Te [J], Te in eV -> K -> J

    E_eV = cs.energy_eV
    sigma = cs.cross_section_m2
    E_J = E_eV * E_CH  # eV -> J

    # Maxwell-Boltzmann Integrand:
    # f(E) = sigma(E) * v(E) * f_MB(E)
    # v(E) = sqrt(2E/me)
    # f_MB(E) = 2/sqrt(pi) * (1/kTe)^(3/2) * sqrt(E) * exp(-E/kTe)
    #
    # Zusammen: sigma * sqrt(2E/me) * 2/sqrt(pi) * (1/kTe)^1.5 * sqrt(E) * exp(-E/kTe)
    # = sigma * 2/sqrt(pi) * (1/kTe)^1.5 * E * sqrt(2/me) * exp(-E/kTe)

    # Vermeidung von Overflow: exp(-E/kTe) kann fuer grosse E sehr klein werden
    with np.errstate(over="ignore", under="ignore"):
        integrand = (
            sigma
            * 2.0 / math.sqrt(math.pi)
            * kTe_J ** (-1.5)
            * E_J  # sqrt(E) * sqrt(E) = E
            * math.sqrt(2.0 / ME)
            * np.exp(-E_J / kTe_J)
            * E_CH  # dE in eV -> dE in J (Jacobian)
        )

    # Ersetze NaN/Inf durch 0
    integrand = np.where(np.isfinite(integrand), integrand, 0.0)

    # Trapez-Integration
    K = np.trapezoid(integrand, E_eV)
    return float(K) if math.isfinite(K) else 0.0


@dataclass
class GasRateModel:
    """Komplettes Ratenmodell fuer ein Gas, bestehend aus tabellierten Querschnitten."""
    gas_name: str = "unknown"
    ionization: Optional[CrossSectionData] = None
    excitation: Optional[CrossSectionData] = None
    elastic: Optional[CrossSectionData] = None

    # Gecachte Werte (LUT)
    _kiz_cache: dict = field(default_factory=dict)
    _kex_cache: dict = field(default_factory=dict)
    _kel_cache: dict = field(default_factory=dict)

    @staticmethod
    def load(gas_dir: str | Path) -> Optional["GasRateModel"]:
        """Lade alle verfuegbaren Querschnittsdaten fuer ein Gas."""
        gas_dir = Path(gas_dir)
        if not gas_dir.is_dir():
            return None

        model = GasRateModel(gas_name=gas_dir.name)
        model.ionization = CrossSectionData.from_csv(gas_dir / "ionization.csv")
        model.excitation = CrossSectionData.from_csv(gas_dir / "excitation.csv")
        model.elastic = CrossSectionData.from_csv(gas_dir / "elastic.csv")
        return model

    def Kiz(self, Te_eV: float) -> Optional[float]:
        """Ionisierungsratenkoeffizient [m^3/s]."""
        if self.ionization is None:
            return None
        # Cache mit 0.01 eV Aufloesung
        key = round(Te_eV, 2)
        if key not in self._kiz_cache:
            self._kiz_cache[key] = compute_rate_coefficient(self.ionization, Te_eV)
        return self._kiz_cache[key]

    def Kex(self, Te_eV: float) -> Optional[float]:
        """Anregungsratenkoeffizient [m^3/s]."""
        if self.excitation is None or len(self.excitation.energy_eV) < 2:
            return None
        key = round(Te_eV, 2)
        if key not in self._kex_cache:
            self._kex_cache[key] = compute_rate_coefficient(self.excitation, Te_eV)
        return self._kex_cache[key]

    def Kel(self, Te_eV: float) -> Optional[float]:
        """Elastischer Stossratenkoeffizient [m^3/s]."""
        if self.elastic is None:
            return None
        key = round(Te_eV, 2)
        if key not in self._kel_cache:
            self._kel_cache[key] = compute_rate_coefficient(self.elastic, Te_eV)
        return self._kel_cache[key]

    def clear_cache(self):
        self._kiz_cache.clear()
        self._kex_cache.clear()
        self._kel_cache.clear()


# ─── Legacy-Fits (Chabert 2012) ───────────────────────────────

def kiz_legacy(Te_eV: float) -> float:
    """Chabert Polynomfit fuer Kiz [m^3/s]."""
    K1 = 6.73e-15 * math.sqrt(Te_eV) * (3.97 + 0.643*Te_eV - 0.0368*Te_eV**2) * math.exp(-12.127/Te_eV)
    K2 = 6.73e-15 * math.sqrt(Te_eV) * (-0.0001031*Te_eV**2 + 6.386*math.exp(-12.127/Te_eV))
    return 0.5 * (K1 + K2)

def kex_legacy(Te_eV: float) -> float:
    """Chabert Arrhenius-Fit fuer Kex [m^3/s]."""
    return 1.2921e-13 * math.exp(-E_CH * 11.6 / (KB * Te_eV * CONV))

def kel_legacy_constant() -> float:
    """Paper-Konstante fuer Kel [m^3/s]."""
    return 1e-13


# ─── Verifikation ─────────────────────────────────────────────

def verify_against_legacy(gas_model: GasRateModel):
    """Vergleiche tabellierte Raten mit Legacy-Fits."""
    print(f"Verifikation: {gas_model.gas_name}")
    print(f"{'Te[eV]':>8} | {'Kiz_tab':>11} {'Kiz_leg':>11} {'Ratio':>7} | "
          f"{'Kel_tab':>11} {'Kel_leg':>11} {'Ratio':>7}")
    print("-" * 75)

    for Te in [2.0, 3.0, 4.0, 5.0, 8.0, 10.0]:
        kiz_t = gas_model.Kiz(Te)
        kiz_l = kiz_legacy(Te)
        kel_t = gas_model.Kel(Te)
        kel_l = kel_legacy_constant()

        r_kiz = kiz_t / kiz_l if kiz_t and kiz_l else 0
        r_kel = kel_t / kel_l if kel_t and kel_l else 0

        print(f"{Te:8.1f} | {kiz_t or 0:11.3e} {kiz_l:11.3e} {r_kiz:7.3f} | "
              f"{kel_t or 0:11.3e} {kel_l:11.3e} {r_kel:7.2f}")


if __name__ == "__main__":
    import sys
    base = Path(__file__).resolve().parent / "cross_sections" / "xenon"
    model = GasRateModel.load(base)
    if model:
        print(f"Gas: {model.gas_name}")
        print(f"Ionisation: {len(model.ionization.energy_eV) if model.ionization else 0} Punkte")
        print(f"Elastisch:  {len(model.elastic.energy_eV) if model.elastic else 0} Punkte")
        print(f"Anregung:   {'vorhanden' if model.excitation and len(model.excitation.energy_eV) > 0 else 'nicht verfuegbar (Arrhenius)'}")
        print()
        verify_against_legacy(model)
    else:
        print("FEHLER: Konnte Xenon-Daten nicht laden!")
        sys.exit(1)
