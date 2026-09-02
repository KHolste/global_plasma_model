// beam_extraction_cpp.hpp -- Beam-Extraktionsmodell fuer den C++-Kern.
//
// Implementiert dasselbe physikalische Modell wie beam_extraction.py:
//   1. Bohm-Fluss mit geometrischer Aufteilung (hL, hR)
//   2. Child-Langmuir Space-Charge-Limit
//   3. Grid-Optik-Effizienz (eta_opt)
//   4. Limitierender Mechanismus (plasma vs space_charge)
//
// Ersetzt den bisherigen einfachen Bohm-Fluss (hL * n * uB * Ai).
// Der alte Pfad ignorierte CL-Limit und eta_opt — was zu systematischen
// Faktor-3-9x-Abweichungen zum Python-Pfad fuehrte.
//
// Referenzen:
//   Chabert et al. 2012, Eq.4-5 (hL, hR, Aeff)
//   Grondein et al. 2016, Sec.II.A (Child-Langmuir)
//   Goebel & Katz 2008, Ch.5 (Grid-Extraktion)
#ifndef BEAM_EXTRACTION_CPP_HPP
#define BEAM_EXTRACTION_CPP_HPP

#include "sim_context.hpp"
#include <cmath>
#include <string>
#include <algorithm>

namespace PhysConst {
    // Hier nur als Reminder — echte Definitionen in phys_const.hpp
}

// ═════════════════════════════════════════════════════════════
// Extraktionsergebnis
// ═════════════════════════════════════════════════════════════
struct ExtractionResult {
    double J_Bohm_total = 0;   // Totale Bohm-Stromdichte am Grid [A/m^2]
    double J_CL = 0;           // Child-Langmuir Limit [A/m^2]
    double J_extracted = 0;    // Tatsaechlich extrahiert [A/m^2]
    double I_beam_mA = 0;      // Beam-Strom [mA]
    double I_CL_limit_mA = 0;  // CL-Limit als Strom [mA]
    double I_Bohm_limit_mA = 0; // Plasma-Limit als Strom [mA]
    double eta_opt = 0;        // Verwendete Grid-Optik-Effizienz [-]
    double f_grid = 0;         // Fraktionaler Grid-seitiger Anteil [-]
    double frac_extracted = 0; // Anteil extrahiert vs. Gesamtverlust [-]
    std::string limiting;      // "plasma" oder "space_charge"
};

// ═════════════════════════════════════════════════════════════
// Beam-Extraktion berechnen
//
// Identisch zur Python-Implementierung in beam_extraction.py:
//   J_Bohm = f_grid * hL * n_ion * uB (grid-seitiger Bohm-Fluss)
//   J_CL = (4/9) * eps0 * sqrt(2e/M) * V^1.5 / d^2
//   J_extracted = min(J_Bohm, J_CL) * eta_opt
//   I_beam = J_extracted * A_grid_open * 1000
// ═════════════════════════════════════════════════════════════

inline ExtractionResult compute_extraction_cpp(
    double n_ion,           // Ionendichte [m^-3]
    double M_ion,           // Ionenmasse [kg]
    double Te_eV,           // Elektronentemperatur [eV]
    double n_neutral,       // Neutraldichte [m^-3]
    double sigma_i,         // Ion-Neutral-Stossquerschnitt [m^2]
    const Thruster& t)      // Thruster-Geometrie (inkl. eta_opt)
{
    using namespace PhysConst;
    ExtractionResult res;
    res.eta_opt = t.eta_opt;

    if (n_ion <= 0 || Te_eV <= 0) return res;

    // ── 1. Bohm-Fluss mit geometrischer Aufteilung ──────────
    double lam = 1.0 / std::max(n_neutral * sigma_i, 1e-10);
    double hL = 0.86 * std::pow(3.0 + t.L / (2.0 * lam), -0.5);
    double hR = 0.80 * std::pow(4.0 + t.R / lam, -0.5);

    double A_end = pi * t.R * t.R;           // Eine Endflaeche [m^2]
    double A_mantel = 2.0 * pi * t.R * t.L;  // Mantelflaeche [m^2]

    double Aeff_grid = hL * A_end;            // Grid-seitige Verlustflaeche
    double Aeff_back = hL * A_end;            // Rueckwand-Verlustflaeche
    double Aeff_wall = hR * A_mantel;         // Mantel-Verlustflaeche
    double Aeff_total = Aeff_grid + Aeff_back + Aeff_wall;

    res.f_grid = (Aeff_total > 0) ? (Aeff_grid / Aeff_total) : 0.5;

    // Bohm-Geschwindigkeit
    double uB_val = std::sqrt(kB * Te_eV * conv / M_ion);

    // Bohm-Stromdichte am Grid
    // Verwende den vollen hL-Fluss zur Grid-Endwand (wie Legacy und Chabert 2012).
    // Die f_grid-Aufteilung (Grid vs Rueckwand vs Mantel) ist fuer die Verlustbilanz
    // relevant, aber der Beam-Strom wird physikalisch durch den direkten Bohm-Fluss
    // zur Grid-Seite bestimmt, nicht durch den fraktionalen Gesamtverlustanteil.
    double Gamma_grid = hL * n_ion * uB_val;
    res.J_Bohm_total = e * Gamma_grid;

    // ── 2. Child-Langmuir Limit ──────────────────────────────
    if (t.Vgrid > 0 && t.sgrid > 0) {
        res.J_CL = (4.0/9.0) * epsilon0 * std::sqrt(2.0*e/M_ion)
                    * std::pow(t.Vgrid, 1.5) / (t.sgrid * t.sgrid);
    } else {
        res.J_CL = 1e10;  // Kein Limit
    }

    // ── 3. Tatsaechliche Extraktion ──────────────────────────
    double J_available = std::min(res.J_Bohm_total, res.J_CL);
    res.J_extracted = J_available * t.eta_opt;

    // Offene Gitterflaeche
    double A_grid_open = t.betai * A_end;

    res.I_beam_mA = res.J_extracted * A_grid_open * 1000.0;
    res.I_CL_limit_mA = res.J_CL * A_grid_open * t.eta_opt * 1000.0;
    res.I_Bohm_limit_mA = res.J_Bohm_total * A_grid_open * t.eta_opt * 1000.0;

    // Limitierender Mechanismus
    res.limiting = (res.J_Bohm_total <= res.J_CL) ? "plasma" : "space_charge";

    // Verlustaufteilung
    double I_total_walls = e * n_ion * uB_val * Aeff_total * 1000.0;
    if (I_total_walls > 0)
        res.frac_extracted = res.I_beam_mA / I_total_walls;

    return res;
}

// ═════════════════════════════════════════════════════════════
// Hauptfunktion: Beam-Strom aus SimContext + PlasmaState
//
// Ersetzt die alte beam_current_mA() mit vollem Extraktionsmodell.
// Rueckwaertskompatibilitaet: Wenn eta_opt=1.0 und sgrid klein,
// naehert sich das Ergebnis dem alten Bohm-Fluss an.
// ═════════════════════════════════════════════════════════════

inline ExtractionResult compute_beam_extraction(const SimContext& ctx, const PlasmaState& s) {
    return compute_extraction_cpp(
        s.n,                // Ionendichte (= Elektronendichte fuer Xenon)
        ctx.gas.M,          // Ionenmasse
        s.Te,               // Elektronentemperatur [eV]
        s.ng,               // Neutraldichte
        ctx.gas.sigma_i,    // Stossquerschnitt
        ctx.thruster);      // Geometrie inkl. eta_opt
}

// Kurzform: nur I_beam in mA (kompatibel mit alter Schnittstelle)
inline double beam_current_mA_full(const SimContext& ctx, const PlasmaState& s) {
    return compute_beam_extraction(ctx, s).I_beam_mA;
}

#endif // BEAM_EXTRACTION_CPP_HPP
