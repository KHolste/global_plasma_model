"""
grid_composition.py -- Gridnahes Kompositionsmodell fuer Iod-ICP-Thruster.

Berechnet die extrahierte Ionenkomposition nicht aus dem globalen
Volumenmittel, sondern aus einer physikalisch korrigierten gridnahen
Zusammensetzung. Beruecksichtigt:

1. Differentielle Diffusion: Leichtere I-Atome diffundieren schneller
   zum Grid als schwerere I2-Molekuele. Am Grid ist daher n_I/n_I2
   hoeher als im Volumenmittel.

2. Lokale Ionisation: Die am Grid ankommenden Neutralen werden lokal
   ionisiert. Da n_I am Grid angereichert ist, steigt auch I+/I2+.

Physik:
  Die Diffusionsgeschwindigkeit skaliert mit sqrt(T/M).
  Fuer I vs I2: v_diff,I / v_diff,I2 = sqrt(M_I2/M_I) = sqrt(2) ~ 1.41

  Das Verhältnis n_I/n_I2 am Grid ist daher um einen Faktor f_enrich
  hoeher als im Volumenmittel:

    f_enrich = sqrt(M_I2 / M_I)     (Diffusionsfaktor)

  Fuer kleine Geometrien (L vergleichbar mit mittlerer freier Weglaenge)
  ist der Effekt staerker. Fuer grosse Geometrien (L >> lambda) geht
  er gegen 1 (gut durchmischt).

  Korrekturfaktor fuer die Geometrieabhaengigkeit:
    f_geom = 1 + (f_diff - 1) * min(lambda_eff / L, 1.0)

  wobei lambda_eff die effektive mittlere freie Weglaenge ist.

Referenz:
  Lieberman & Lichtenberg 2005, Ch.5 (Diffusion in discharges)
  Chabert & Braithwaite 2011, Ch.3 (Species transport)

Annahmen:
  - Einheitliche Gastemperatur (0D Tg)
  - Nur I und I2 als relevante Neutrale
  - Keine separate radiale Diffusionsbehandlung
  - Ionenzusammensetzung folgt aus lokaler Neutralzusammensetzung
"""
from __future__ import annotations

import math
from dataclasses import dataclass

KB = 1.3806504e-23
PI = math.pi


@dataclass
class GridComposition:
    """Gridnahe Zusammensetzung im Vergleich zum Volumenmittel."""

    # Globale (Volumen) Werte
    n_I_bulk: float = 0.0
    n_I2_bulk: float = 0.0
    n_Ip_bulk: float = 0.0
    n_I2p_bulk: float = 0.0
    diss_bulk: float = 0.0       # n_I / (n_I + 2*n_I2)
    f_Ip_bulk: float = 0.0       # n_I+ / (n_I+ + n_I2+)

    # Gridnahe Werte (korrigiert)
    n_I_grid: float = 0.0
    n_I2_grid: float = 0.0
    n_Ip_grid: float = 0.0
    n_I2p_grid: float = 0.0
    diss_grid: float = 0.0
    f_Ip_grid: float = 0.0

    # Korrekturfaktoren
    f_enrich: float = 1.0        # n_I/n_I2 Anreicherungsfaktor
    f_diff: float = 1.0          # Reiner Diffusionsfaktor sqrt(M_I2/M_I)
    f_geom: float = 1.0          # Geometrie-Korrektur

    def summary(self) -> str:
        return (f"f_I+: bulk={self.f_Ip_bulk:.3f} grid={self.f_Ip_grid:.3f} "
                f"(+{(self.f_Ip_grid/max(self.f_Ip_bulk,1e-10)-1)*100:.1f}%), "
                f"f_enrich={self.f_enrich:.3f}")


def compute_grid_composition(
    n_I: float,
    n_I2: float,
    n_Ip: float,
    n_I2p: float,
    M_I: float = 2.107e-25,
    M_I2: float = 4.214e-25,
    Te_eV: float = 4.0,
    Tg: float = 300.0,
    L: float = 0.023,
    sigma_II2: float = 5.5e-19,   # I-I2 Stossquerschnitt [m^2]
    n_total_neutral: float = 0.0,  # Fuer lambda-Berechnung
) -> GridComposition:
    """Berechne gridnahe Zusammensetzung aus Volumenmittel.

    Args:
        n_I, n_I2: Volumenmittel-Dichten der Neutralen [m^-3]
        n_Ip, n_I2p: Volumenmittel-Dichten der Ionen [m^-3]
        M_I, M_I2: Massen [kg]
        Te_eV: Elektronentemperatur (fuer Ionisationsraten)
        Tg: Gastemperatur [K]
        L: Kammerlaenge [m]
        sigma_II2: I-I2 Stossquerschnitt [m^2]
        n_total_neutral: Gesamte Neutraldichte [m^-3]

    Returns:
        GridComposition mit Bulk- und Grid-Werten
    """
    gc = GridComposition()

    # Bulk-Werte
    gc.n_I_bulk = n_I
    gc.n_I2_bulk = n_I2
    gc.n_Ip_bulk = n_Ip
    gc.n_I2p_bulk = n_I2p

    n_total = n_I + 2 * n_I2
    gc.diss_bulk = n_I / n_total if n_total > 0 else 0
    n_ion_total = n_Ip + n_I2p
    gc.f_Ip_bulk = n_Ip / n_ion_total if n_ion_total > 0 else 0

    # ── Diffusionsbasierte Anreicherung ───────────────────────
    # Reiner Massenfaktor: leichtere Spezies diffundiert schneller
    gc.f_diff = math.sqrt(M_I2 / M_I)  # = sqrt(2) ~ 1.414 fuer Iod

    # Geometrie-Korrekturfaktor: staerker bei kleinen L
    # lambda_eff = mittlere freie Weglaenge fuer I in I2
    if n_total_neutral > 0 and n_total_neutral > 1e10:
        lambda_eff = 1.0 / (n_total_neutral * sigma_II2)
    else:
        lambda_eff = L  # Freimolekularer Grenzfall

    # Knudsen-Zahl Kn = lambda/L
    Kn = lambda_eff / L if L > 0 else 1.0

    # Bei Kn >> 1 (freier Flug): voller Effekt
    # Bei Kn << 1 (viskos): Diffusion mittelt aus -> weniger Effekt
    # Uebergangsfunktion: f_geom = 1 + (f_diff - 1) * min(Kn, 1.0)
    gc.f_geom = 1.0 + (gc.f_diff - 1.0) * min(Kn, 1.0)

    gc.f_enrich = gc.f_geom

    # ── Gridnahe Neutralzusammensetzung ───────────────────────
    # (n_I/n_I2)_grid = (n_I/n_I2)_bulk * f_enrich
    # Normiert: n_I_grid + n_I2_grid = n_I + n_I2 (Masse erhalten)
    if n_I2 > 0:
        ratio_bulk = n_I / n_I2
        ratio_grid = ratio_bulk * gc.f_enrich
        # Aus ratio_grid und n_total = n_I + n_I2:
        n_neutral_total = n_I + n_I2
        gc.n_I2_grid = n_neutral_total / (1 + ratio_grid)
        gc.n_I_grid = n_neutral_total - gc.n_I2_grid
    else:
        gc.n_I_grid = n_I
        gc.n_I2_grid = 0

    gc.diss_grid = gc.n_I_grid / (gc.n_I_grid + 2*gc.n_I2_grid) if (gc.n_I_grid + 2*gc.n_I2_grid) > 0 else 0

    # ── Gridnahe Ionenkomposition ─────────────────────────────
    # Die Ionenkomposition am Grid folgt aus dem gridnahen Neutral-
    # verhaeltnis. Da die Ionisationsraten K_iz,I und K_iz,I2 fuer
    # beide Spezies aehnlich sind (gleiche Groessenordnung), skaliert
    # das I+/I2+ Verhaeltnis naeherungsweise mit n_I/n_I2:
    #
    #   (I+/I2+)_grid = (I+/I2+)_bulk * (n_I/n_I2)_grid / (n_I/n_I2)_bulk
    #                  = (I+/I2+)_bulk * f_enrich
    #
    # Das ist physikalisch: mehr atomare Neutrale am Grid -> mehr I+ Produktion.

    if n_I2 > 0 and n_I2p > 0:
        ion_ratio_bulk = n_Ip / n_I2p
        ion_ratio_grid = ion_ratio_bulk * gc.f_enrich

        # Normiert: n_Ip_grid + n_I2p_grid = n_Ip + n_I2p
        gc.n_I2p_grid = n_ion_total / (1 + ion_ratio_grid)
        gc.n_Ip_grid = n_ion_total - gc.n_I2p_grid
    else:
        gc.n_Ip_grid = n_Ip
        gc.n_I2p_grid = n_I2p

    gc.f_Ip_grid = gc.n_Ip_grid / n_ion_total if n_ion_total > 0 else gc.f_Ip_bulk

    return gc
