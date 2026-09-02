// beam_extraction_cpp.hpp -- Extraktion, Schub und Wirkungsgrade aus einer Quelle.
//
// Frueher standen hier zwei getrennte Rechnungen: der Strahlstrom kam aus dem
// vollen Extraktionsmodell mit Raumladungsgrenze, Schub und Massenwirkungsgrad
// dagegen aus dem reinen Bohm-Fluss. Bei begrenzter Extraktion widersprachen
// sich die Ausgaben. Jetzt kommt alles aus einer Rechnung, und diese summiert
// ueber die Ionensorten und fuehrt die Ladungszahl mit.
//
// Modell:
//   1. Bohm-Fluss je Ionensorte zur Gitterseite, mit deren eigener Masse und
//      Ladungszahl:  u_B,i = sqrt(Z_i * k * Te / M_i),  Gamma_i = hL * n_i * u_B,i
//   2. Stromdichte je Sorte:  J_i = Z_i * e * Gamma_i
//   3. Raumladungsgrenze je Sorte nach Child-Langmuir mit deren Z und M. Fuer
//      ein Gemisch ist die Grenze erreicht, wenn sum_i J_i / J_CL,i = 1; ist
//      die Summe groesser, wird der gesamte Strahl mit demselben Faktor
//      gedrosselt, die Zusammensetzung bleibt erhalten.
//   4. Gitteroptik eta_opt wirkt auf alles gleichermassen.
//   5. Schub aus dem tatsaechlich extrahierten Teilchenstrom je Sorte mal
//      deren Masse und Austrittsgeschwindigkeit v_i = sqrt(2 Z_i e U / M_i).
//   6. Massenwirkungsgrad als Verhaeltnis von ausgetragenem zu zugefuehrtem
//      Massenstrom -- nicht als Teilchenverhaeltnis, das bei molekularem
//      Treibstoff seine Bedeutung verliert.
//
// Referenzen:
//   Chabert et al. 2012, Eq.4-5 (hL, hR, Aeff)
//   Grondein et al. 2016, Sec.II.A (Child-Langmuir)
//   Goebel & Katz 2008, Ch.5 (Grid-Extraktion, Mehrfachladung)
#ifndef BEAM_EXTRACTION_CPP_HPP
#define BEAM_EXTRACTION_CPP_HPP

#include "sim_context.hpp"
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

// ═════════════════════════════════════════════════════════════
// Eingabe: der Zustand, soweit die Extraktion ihn braucht
// ═════════════════════════════════════════════════════════════

struct BeamIon {
    std::string id;
    double n = 0;        // Dichte [m^-3]
    double M = 0;        // Masse [kg]
    int    Z = 1;        // Ladungszahl
};

struct BeamNeutral {
    std::string id;
    double n = 0;        // Dichte [m^-3]
    double M = 0;        // Masse [kg]
};

struct ExtractionInput {
    std::vector<BeamIon> ions;
    std::vector<BeamNeutral> neutrals;
    double Te_eV = 0;
    double Tg_K = 0;
    double sigma_i = 1e-18;      // Ion-Neutral-Stossquerschnitt [m^2]
    double mdot_in_kg_s = 0;     // zugefuehrter Massenstrom [kg/s]
};

// ═════════════════════════════════════════════════════════════
// Ergebnis
// ═════════════════════════════════════════════════════════════

struct IonBeamShare {
    std::string id;
    int    Z = 1;
    double u_B = 0;          // Bohm-Geschwindigkeit [m/s]
    double J_bohm = 0;       // Stromdichte am Gitter, ungedrosselt [A/m^2]
    double J_extracted = 0;  // tatsaechlich extrahiert [A/m^2]
    double I_mA = 0;         // Strom dieser Sorte [mA]
    double thrust_N = 0;     // Schubbeitrag [N]
    double mdot_kg_s = 0;    // ausgetragener Massenstrom [kg/s]
    double v_exhaust = 0;    // Austrittsgeschwindigkeit [m/s]
};

struct ExtractionResult {
    std::vector<IonBeamShare> ions;

    double J_Bohm_total = 0;    // Summe der ungedrosselten Stromdichten [A/m^2]
    double J_extracted = 0;     // Summe der extrahierten Stromdichten [A/m^2]
    double I_beam_mA = 0;       // Strahlstrom [mA]
    double I_CL_limit_mA = 0;   // Strom an der Raumladungsgrenze [mA]
    double I_Bohm_limit_mA = 0; // Strom, den das Plasma anbietet [mA]
    double space_charge_ratio = 0;  // sum_i J_i / J_CL,i
    double throttle = 1.0;      // Drosselfaktor aus der Raumladungsgrenze
    double eta_opt = 1.0;

    double thrust_ions_N = 0;
    double thrust_neutrals_N = 0;
    double thrust_total_N = 0;

    double mdot_out_kg_s = 0;   // ausgetragene Ionenmasse [kg/s]
    double eta_mass = 0;        // ausgetragener zu zugefuehrtem Massenstrom
    double P_beam_W = 0;        // Strahlleistung, fuer den Kopplungswirkungsgrad

    double hL = 0, hR = 0;
    std::string limiting;       // "plasma" oder "space_charge"
};

// ═════════════════════════════════════════════════════════════
// Rechnung
// ═════════════════════════════════════════════════════════════

inline ExtractionResult compute_extraction(const ExtractionInput& in, const Thruster& t) {
    using namespace PhysConst;
    ExtractionResult res;
    res.eta_opt = t.eta_opt;
    res.limiting = "plasma";
    if (in.Te_eV <= 0) return res;

    // ── Randschichtfaktoren aus der mittleren freien Weglaenge ──
    double n_neut_total = 0;
    for (const auto& nu : in.neutrals) n_neut_total += nu.n;
    const double lam = 1.0 / std::max(n_neut_total * in.sigma_i, 1e-10);
    res.hL = 0.86 * std::pow(3.0 + t.L / (2.0 * lam), -0.5);
    res.hR = 0.80 * std::pow(4.0 + t.R / lam, -0.5);

    const double A_end = pi * t.R * t.R;
    const double A_open = t.betai * A_end;     // offene Gitterflaeche

    // ── Bohm-Fluss und Raumladungsgrenze je Ionensorte ──────────
    struct Zwischen { double J = 0, J_CL = 0; };
    std::vector<Zwischen> zw(in.ions.size());

    for (size_t k = 0; k < in.ions.size(); ++k) {
        const BeamIon& ion = in.ions[k];
        IonBeamShare sh;
        sh.id = ion.id;
        sh.Z = ion.Z;
        if (ion.n <= 0 || ion.M <= 0 || ion.Z <= 0) { res.ions.push_back(sh); continue; }

        // Bohm-Geschwindigkeit mit Ladungszahl: die Vorschicht beschleunigt
        // ein Z-fach geladenes Ion durch das Z-fache Potentialgefaelle.
        sh.u_B = std::sqrt(ion.Z * kB * in.Te_eV * conv / ion.M);
        zw[k].J = ion.Z * e * res.hL * ion.n * sh.u_B;
        sh.J_bohm = zw[k].J;

        zw[k].J_CL = (t.Vgrid > 0 && t.sgrid > 0)
            ? (4.0/9.0) * epsilon0 * std::sqrt(2.0 * ion.Z * e / ion.M)
              * std::pow(t.Vgrid, 1.5) / (t.sgrid * t.sgrid)
            : 1e30;

        sh.v_exhaust = std::sqrt(2.0 * ion.Z * e * t.Vgrid / ion.M);
        res.ions.push_back(sh);
        res.J_Bohm_total += zw[k].J;
        if (zw[k].J_CL > 0) res.space_charge_ratio += zw[k].J / zw[k].J_CL;
    }

    // ── Drosselung, wenn die Raumladung nicht mehr traegt ────────
    res.throttle = (res.space_charge_ratio > 1.0) ? 1.0 / res.space_charge_ratio : 1.0;
    if (res.space_charge_ratio > 1.0) res.limiting = "space_charge";

    // ── Stroeme, Schub und Massenstrom ──────────────────────────
    for (size_t k = 0; k < in.ions.size(); ++k) {
        IonBeamShare& sh = res.ions[k];
        const BeamIon& ion = in.ions[k];
        if (sh.J_bohm <= 0) continue;

        sh.J_extracted = sh.J_bohm * res.throttle * t.eta_opt;
        sh.I_mA = sh.J_extracted * A_open * 1000.0;

        // Teilchenstrom aus der Stromdichte, mit der Ladungszahl zurueck
        const double N_dot = sh.J_extracted * A_open / (ion.Z * e);   // 1/s
        sh.mdot_kg_s = N_dot * ion.M;
        sh.thrust_N = sh.mdot_kg_s * sh.v_exhaust;

        res.J_extracted += sh.J_extracted;
        res.I_beam_mA += sh.I_mA;
        res.thrust_ions_N += sh.thrust_N;
        res.mdot_out_kg_s += sh.mdot_kg_s;
        res.P_beam_W += 0.5 * sh.mdot_kg_s * sh.v_exhaust * sh.v_exhaust;
    }

    res.I_Bohm_limit_mA = res.J_Bohm_total * A_open * t.eta_opt * 1000.0;
    res.I_CL_limit_mA = (res.space_charge_ratio > 0)
        ? res.J_Bohm_total / res.space_charge_ratio * A_open * t.eta_opt * 1000.0
        : 0.0;

    // ── Schub der ausstroemenden Neutralteilchen ────────────────
    for (const auto& nu : in.neutrals) {
        if (nu.n <= 0 || nu.M <= 0 || in.Tg_K <= 0) continue;
        const double v_mean = std::sqrt(8 * kB * in.Tg_K / (pi * nu.M));
        const double Gamma = 0.25 * nu.n * v_mean;
        const double mdot = Gamma * t.Ag * nu.M;
        res.thrust_neutrals_N += mdot * v_mean;
        res.P_beam_W += 0.5 * mdot * v_mean * v_mean;
    }
    res.thrust_total_N = res.thrust_ions_N + res.thrust_neutrals_N;

    // ── Massenwirkungsgrad als Massenverhaeltnis ────────────────
    res.eta_mass = (in.mdot_in_kg_s > 0) ? res.mdot_out_kg_s / in.mdot_in_kg_s : 0.0;

    return res;
}

// ═════════════════════════════════════════════════════════════
// Anschluss des fest verdrahteten Wegs: eine Ionensorte, Ladung eins
// ═════════════════════════════════════════════════════════════

inline ExtractionInput extraction_input_legacy(const SimContext& ctx, const PlasmaState& s) {
    ExtractionInput in;
    in.ions.push_back({ctx.gas.species, s.n, ctx.gas.M, 1});
    in.neutrals.push_back({ctx.gas.species, s.ng, ctx.gas.M});
    in.Te_eV = s.Te;
    in.Tg_K = s.Tg;
    in.sigma_i = ctx.gas.sigma_i;
    in.mdot_in_kg_s = ctx.thruster.Q0 * ctx.gas.M;
    return in;
}

inline ExtractionResult compute_beam_extraction(const SimContext& ctx, const PlasmaState& s) {
    return compute_extraction(extraction_input_legacy(ctx, s), ctx.thruster);
}

inline double beam_current_mA_full(const SimContext& ctx, const PlasmaState& s) {
    return compute_beam_extraction(ctx, s).I_beam_mA;
}

#endif // BEAM_EXTRACTION_CPP_HPP
