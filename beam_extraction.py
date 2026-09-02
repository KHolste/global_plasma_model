"""
beam_extraction.py -- Verbessertes Beam-Extraktions- und Wandverlustmodell.

Trennt den Ionenverlustpfad in:
1. Bohm-Fluss zur Plasma-Sheath-Grenze (physikalisch)
2. Child-Langmuir-Limitierung der Extraktion (space-charge)
3. Grid-Optik-Effizienz (device-spezifisch)
4. Wandneutralisation (Bilanz: was nicht extrahiert wird)

Referenzen:
- Chabert et al. 2012, Eq.4-5 (hL, hR, Aeff)
- Grondein et al. 2016, Sec.II.A (Child-Langmuir)
- Lafleur et al. 2022, Sec.2-3 (NPT30-I2 Grid-Optik)
- Goebel & Katz 2008, Ch.5 (Grid-Extraktion)

Annahmen (Version 1):
- Elektropositive Sheath (alpha < 1)
- Gleiche hL, hR fuer alle positiven Ionen
- Gemischter CL-Strom: gewichteter Mittelwert ueber Ionenmassen
- Grid-Optik als skalarer Effizienzfaktor (device-spezifisch)
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

# Konstanten
E_CH = 1.602176487e-19
EPS0 = 8.854187817e-12
KB = 1.3806504e-23
CONV = 11604.5250061657
PI = math.pi


@dataclass
class ExtractionResult:
    """Ergebnis der Beam-Extraktionsberechnung."""
    # Pro Ionenart
    ion_species: list[str] = field(default_factory=list)
    Gamma_Bohm: dict[str, float] = field(default_factory=dict)  # Bohm-Fluss [m^-2 s^-1]
    J_Bohm: dict[str, float] = field(default_factory=dict)      # Stromdichte [A/m^2]
    uB: dict[str, float] = field(default_factory=dict)           # Bohm-Geschwindigkeit [m/s]

    # Aggregiert
    J_Bohm_total: float = 0.0     # Totale Bohm-Stromdichte [A/m^2]
    J_CL: float = 0.0             # Child-Langmuir Limit [A/m^2]
    J_extracted: float = 0.0      # Tatsaechlich extrahiert [A/m^2]
    I_beam_mA: float = 0.0        # Beam-Strom [mA]
    I_CL_limit_mA: float = 0.0    # CL-Limit [mA]
    I_Bohm_limit_mA: float = 0.0  # Plasma-Limit [mA]

    # Verlustaufteilung
    frac_extracted: float = 0.0   # Anteil extrahiert [-]
    frac_wall: float = 0.0        # Anteil Wandneutralisation [-]
    limiting: str = ""             # "plasma" oder "space_charge"

    # Grid-Optik
    eta_opt: float = 0.0          # Grid-Optik-Effizienz [-]

    # Pro-Spezies Beam
    I_beam_species: dict[str, float] = field(default_factory=dict)  # [mA]
    f_beam_species: dict[str, float] = field(default_factory=dict)  # Fraktionen [-]


def compute_extraction(
    ion_densities: dict[str, float],    # {species_id: density [m^-3]}
    ion_masses: dict[str, float],       # {species_id: mass [kg]}
    Te_eV: float,
    n_neutral_total: float,
    sigma_i: float,                     # Ion-Neutral QS [m^2]
    R: float,                           # Kammerradius [m]
    L: float,                           # Kammerlaenge [m]
    betai: float,                       # Ionentransparenz [-]
    Vgrid: float,                       # Gitterspannung [V]
    sgrid: float,                       # Gitterabstand [m]
    eta_opt: float = 0.25,              # Grid-Optik-Effizienz [-]
) -> ExtractionResult:
    """Berechne Beam-Extraktion mit CL-Limitierung und Grid-Optik.

    Args:
        ion_densities: Dichten der positiven Ionen
        ion_masses: Massen der positiven Ionen
        Te_eV: Elektronentemperatur [eV]
        n_neutral_total: Gesamte Neutraldichte [m^-3]
        sigma_i: Ion-Neutral-Stossquerschnitt [m^2]
        R, L: Kammergeometrie [m]
        betai: Ionengittertransparenz
        Vgrid: Gitterspannung [V]
        sgrid: Gitterabstand [m]
        eta_opt: Grid-Optik-Effizienz (0-1)

    Returns:
        ExtractionResult mit allen Diagnosegroessen
    """
    res = ExtractionResult()
    res.eta_opt = eta_opt
    res.ion_species = list(ion_densities.keys())

    if not ion_densities or Te_eV <= 0:
        return res

    # ── 1. Bohm-Fluss und geometrische Aufteilung ───────────
    # hL, hR nach Chabert/Lieberman
    lam = 1.0 / max(n_neutral_total * sigma_i, 1e-10)
    hL = 0.86 * (3 + L / (2 * lam))**(-0.5)
    hR = 0.80 * (4 + R / lam)**(-0.5)

    # Verlustflaechen (geometrisch getrennt)
    A_end = PI * R**2             # Eine Endflaeche
    A_mantel = 2 * PI * R * L    # Mantelflaeche

    # Effektive Verlustfluesse (hL/hR-gewichtet)
    # Grid-Seite: eine Endflaeche mit Gitter
    Aeff_grid_side = hL * A_end
    # Geschlossene Endwand: andere Endflaeche
    Aeff_back_wall = hL * A_end
    # Mantelwand
    Aeff_mantel = hR * A_mantel
    # Gesamte Ionenverlustflaeche
    Aeff_total = Aeff_grid_side + Aeff_back_wall + Aeff_mantel

    # Fraktionaler Anteil, der zur Grid-Seite fliesst
    f_grid = Aeff_grid_side / Aeff_total if Aeff_total > 0 else 0.5

    # Offene Gitterflaeche
    A_grid_open = betai * A_end

    for sp_id, n_ion in ion_densities.items():
        M = ion_masses.get(sp_id, 2.107e-25)
        uB_sp = math.sqrt(KB * Te_eV * CONV / M)
        # Gesamter Bohm-Fluss (alle Waende)
        Gamma_total = n_ion * uB_sp  # Flussdichte [m^-2 s^-1]
        # Davon: nur Grid-seitiger Anteil ist beamverfuegbar
        Gamma_grid = f_grid * hL * n_ion * uB_sp
        J_grid = E_CH * Gamma_grid   # Effektive beamverfuegbare Stromdichte [A/m^2]

        res.uB[sp_id] = uB_sp
        res.Gamma_Bohm[sp_id] = Gamma_grid  # Nur Grid-seitiger Anteil
        res.J_Bohm[sp_id] = J_grid
        res.J_Bohm_total += J_grid

    # ── 2. Child-Langmuir Limit ───────────────────────────────
    # J_CL = (4/9) * eps0 * sqrt(2e/M_eff) * Vgrid^1.5 / sgrid^2
    # Fuer gemischte Ionen: effektive Masse gewichtet nach Stromdichteanteil
    if res.J_Bohm_total > 0:
        M_eff = 0.0
        for sp_id in res.ion_species:
            J_sp = res.J_Bohm.get(sp_id, 0)
            M_sp = ion_masses.get(sp_id, 2.107e-25)
            M_eff += (J_sp / res.J_Bohm_total) * M_sp
    else:
        M_eff = 2.107e-25  # Default: atomares I

    if Vgrid > 0 and sgrid > 0:
        res.J_CL = (4.0/9.0) * EPS0 * math.sqrt(2*E_CH/M_eff) * Vgrid**1.5 / sgrid**2
    else:
        res.J_CL = 1e10  # Kein Limit

    # ── 3. Tatsaechliche Extraktion ───────────────────────────
    # J_extracted = min(J_Bohm, J_CL) * eta_opt
    J_available = min(res.J_Bohm_total, res.J_CL)
    res.J_extracted = J_available * eta_opt

    # Beam-Strom
    res.I_beam_mA = res.J_extracted * A_grid_open * 1000  # mA
    res.I_CL_limit_mA = res.J_CL * A_grid_open * eta_opt * 1000
    res.I_Bohm_limit_mA = res.J_Bohm_total * A_grid_open * eta_opt * 1000

    # Limitierender Mechanismus
    if res.J_Bohm_total <= res.J_CL:
        res.limiting = "plasma"
    else:
        res.limiting = "space_charge"

    # ── 4. Pro-Spezies Beam-Aufteilung ────────────────────────
    # Proportional zum Bohm-Stromanteil
    for sp_id in res.ion_species:
        if res.J_Bohm_total > 0:
            frac = res.J_Bohm[sp_id] / res.J_Bohm_total
        else:
            frac = 0.0
        res.I_beam_species[sp_id] = frac * res.I_beam_mA
        res.f_beam_species[sp_id] = frac

    # ── 5. Verlustaufteilung ──────────────────────────────────
    # Gesamter Ionenverlust (alle Waende) fuer Bilanz
    n_ion_total = sum(ion_densities.values())
    M_avg = sum(ion_masses[s]*ion_densities[s] for s in ion_densities) / max(n_ion_total, 1) if n_ion_total > 0 else 2.107e-25
    uB_avg = math.sqrt(KB * Te_eV * CONV / M_avg)
    I_total_all_walls = E_CH * n_ion_total * uB_avg * Aeff_total * 1000  # mA (gesamter Ionenverlust)
    if I_total_all_walls > 0:
        res.frac_extracted = res.I_beam_mA / I_total_all_walls
        res.frac_wall = 1.0 - res.frac_extracted
    else:
        res.frac_extracted = 0.0
        res.frac_wall = 1.0

    return res
