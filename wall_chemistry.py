"""
wall_chemistry.py -- Erweiterte Wandchemie fuer Iod-ICP-Globalmodelle.

Modelliert die Konkurrenz zwischen gasphasiger atomarer Bildung und
oberflaecheninduzierter Rueckrekombination. Getrennte Behandlung
verschiedener Wandflaechentypen.

Wandtypen:
  mantel:    Zylindrische Mantelflaeche (Oxidkeramik, Edelstahl, ...)
  back_wall: Geschlossene Rueckwand
  grid:      Grid-Flaeche (metallisch, teilweise offen)

Physik:
  Oberflaechenrekombinationsrate:
    R_wall = gamma * (1/4) * n_I * v_th * A_eff / V

  Dabei ist A_eff flaechengewichtet:
    A_eff = sum(gamma_i * A_i) / gamma_ref

  Auf der Grid-Flaeche: Ein Anteil beta_g der Flaeche ist offen.
  I-Atome die durch offene Loecher fliegen, rekombinieren NICHT.
  Nur I auf der geschlossenen Grid-Flaeche (1-beta_g)*A_grid rekombiniert.

Referenzen:
  Lafleur 2022: gamma = 0.07 (gefittet, global)
  Grondein 2016: gamma = 0.02 (konservativ)
  Lafleur 2026 Review Sec.4: gamma = 0.005 - 0.44 (materialabhaengig)
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

PI = math.pi
KB = 1.3806504e-23


@dataclass
class WallChemistry:
    """Flaechenspezifische Wandrekombination."""

    # Rekombinationskoeffizienten [-]
    gamma_mantel: float = 0.02       # Zylindrischer Mantel
    gamma_back: float = 0.02         # Geschlossene Rueckwand
    gamma_grid_closed: float = 0.02  # Geschlossener Teil der Grid-Flaeche

    # Geometrie (aus ThrusterGeometry berechnet)
    R: float = 0.015
    L: float = 0.023
    betag: float = 0.25  # Gastransparenz des Grids

    @property
    def A_mantel(self) -> float:
        return 2 * PI * self.R * self.L

    @property
    def A_back(self) -> float:
        return PI * self.R**2

    @property
    def A_grid_total(self) -> float:
        return PI * self.R**2

    @property
    def A_grid_open(self) -> float:
        return self.betag * self.A_grid_total

    @property
    def A_grid_closed(self) -> float:
        return (1 - self.betag) * self.A_grid_total

    @property
    def V(self) -> float:
        return PI * self.R**2 * self.L

    @property
    def A_V_ratio(self) -> float:
        """Flaechen-zu-Volumen-Verhaeltnis [m^-1]."""
        A_total = self.A_mantel + self.A_back + self.A_grid_total
        return A_total / self.V

    def effective_gamma_A(self) -> float:
        """Effektive gamma * A fuer Gesamtrekombination [m^2]."""
        return (self.gamma_mantel * self.A_mantel +
                self.gamma_back * self.A_back +
                self.gamma_grid_closed * self.A_grid_closed)

    def wall_loss_rate(self, n_I: float, Tg: float, M_I: float = 2.107e-25) -> float:
        """Atomare Wandverlustrate [m^-3 s^-1].

        Returns: dn_I/dt Verlust durch Wandrekombination
        """
        v_th = math.sqrt(8 * KB * Tg / (PI * M_I))
        gamma_A_eff = self.effective_gamma_A()
        return 0.25 * v_th * n_I * gamma_A_eff / self.V

    def grid_outflow_rate(self, n_I: float, Tg: float, M_I: float = 2.107e-25) -> float:
        """Atomarer Verlust durch offene Grid-Loecher [m^-3 s^-1].

        I-Atome die durch das Grid fliegen, gehen verloren (keine Rekombination).
        """
        v_th = math.sqrt(8 * KB * Tg / (PI * M_I))
        return 0.25 * v_th * n_I * self.A_grid_open / self.V


@dataclass
class WallDiagnostics:
    """Diagnostik fuer die Konkurrenz Gasphase vs. Wand."""

    # Raten [m^-3 s^-1]
    R_diss_direct: float = 0.0      # Direkte Dissoziation
    R_diss_prediss: float = 0.0     # Praedissoziation
    R_diss_dissatt: float = 0.0     # Diss. Attachment (produziert auch I)
    R_diss_dissiz: float = 0.0      # Diss. Ionisation (produziert auch I)
    R_diss_total: float = 0.0       # Gesamt atomare Produktion

    R_wall_rec: float = 0.0         # Wandrekombination I -> 0.5 I2
    R_grid_outflow: float = 0.0     # Verlust durch Grid (kein Rekombi)
    R_wall_total: float = 0.0       # Gesamter atomarer Wandverlust

    # Verhaeltnisse
    ratio_wall_to_production: float = 0.0  # R_wall / R_diss_total
    regime: str = ""                        # "gas_dominated" oder "wall_dominated"

    A_V: float = 0.0                # A/V [m^-1]
    gamma_eff: float = 0.0          # Effektives gamma [-]

    def summary(self) -> str:
        return (f"R_diss={self.R_diss_total:.2e}, R_wall={self.R_wall_total:.2e}, "
                f"ratio={self.ratio_wall_to_production:.2f}, {self.regime}, "
                f"A/V={self.A_V:.0f}")


def diagnose_wall_competition(
    state: dict,
    chem_rates: dict[str, float],
    wall: WallChemistry,
    Tg: float = 300.0,
) -> WallDiagnostics:
    """Diagnostiziere die Konkurrenz zwischen Gasphase und Wandrekombination.

    Args:
        state: {"I": n_I, "I2": n_I2, "ne": ne, ...}
        chem_rates: {"diss_I2": R, "exc_I2_prediss": R, "dissatt_I2": R, "dissiz_I2": R}
        wall: WallChemistry Konfiguration
        Tg: Gastemperatur [K]

    Returns:
        WallDiagnostics
    """
    diag = WallDiagnostics()

    # Gasphasen-Produktion von I
    diag.R_diss_direct = chem_rates.get("diss_I2", 0)
    diag.R_diss_prediss = chem_rates.get("exc_I2_prediss", 0)
    diag.R_diss_dissatt = chem_rates.get("dissatt_I2", 0)
    diag.R_diss_dissiz = chem_rates.get("dissiz_I2", 0)
    diag.R_diss_total = (diag.R_diss_direct + diag.R_diss_prediss +
                          diag.R_diss_dissatt + diag.R_diss_dissiz)

    # Wandverluste
    n_I = state.get("I", 0)
    diag.R_wall_rec = wall.wall_loss_rate(n_I, Tg)
    diag.R_grid_outflow = wall.grid_outflow_rate(n_I, Tg)
    diag.R_wall_total = diag.R_wall_rec + diag.R_grid_outflow

    # Konkurrenz-Verhaeltnis
    if diag.R_diss_total > 0:
        diag.ratio_wall_to_production = diag.R_wall_total / diag.R_diss_total
    else:
        diag.ratio_wall_to_production = float('inf')

    diag.regime = "wall_dominated" if diag.ratio_wall_to_production > 1.0 else "gas_dominated"
    diag.A_V = wall.A_V_ratio
    diag.gamma_eff = wall.effective_gamma_A() / (wall.A_mantel + wall.A_back + wall.A_grid_total) if (wall.A_mantel + wall.A_back + wall.A_grid_total) > 0 else 0

    return diag
