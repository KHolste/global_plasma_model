// chem_system.hpp -- Generisches Chemie-System fuer den C++-Kern.
//
// Ersetzt den hart kodierten Xenon-Spezialfall durch ein abstrahiertes
// Spezies-/Reaktionsmodell mit generischem Residualaufbau.
//
// Architektur (analog zu Python BalanceAssembler):
//   ChemSpecies     -- Spezies mit Typ, Masse, Ladung, Eigenschaften
//   ChemReaction    -- Reaktion mit Stoechiometrie, Rate, Energieverlust
//   ChemSystem      -- Sammlung von Spezies + Reaktionen
//   GenericResidual -- Assembliert Residualvektor aus ChemSystem
//
// Zustandsvektor: [n_species_0, ..., n_species_N, Te, Tg]
// Elektronendichte aus Quasineutralitaet.
//
// Der Xenon-Fall wird als erster Nutzer dieser Struktur aufgesetzt.
#ifndef CHEM_SYSTEM_HPP
#define CHEM_SYSTEM_HPP

#include "phys_const.hpp"
#include "sim_context.hpp"
#include "rates.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <map>

// ═════════════════════════════════════════════════════════════
// Spezies-Typen
// ═════════════════════════════════════════════════════════════
enum class SpeciesType {
    NEUTRAL_ATOM,
    NEUTRAL_MOLECULE,
    POSITIVE_ION,
    NEGATIVE_ION,
};

struct ChemSpecies {
    std::string id;           // z.B. "Xe", "Xe+", "I", "I2", "I+"
    SpeciesType type;
    double mass_kg;           // Teilchenmasse [kg]
    int charge = 0;           // Ladungszahl (0, +1, -1, ...)
    double thermal_cond = 0;  // Waermeleitfaehigkeit [W/(m*K)]
    bool is_feedstock = false;
    double sigma_i = 1e-18;   // Ion-Neutral Stossquerschnitt [m^2]
    // Nachgereichte Felder stehen am Ende, damit die bestehenden
    // Klammerinitialisierungen unveraendert gueltig bleiben.
    std::string name;         // Anzeigename, sonst gleich der id
    bool is_beam_extracted = false;  // Wird ueber das Gitter extrahiert?

    bool is_neutral() const { return type == SpeciesType::NEUTRAL_ATOM || type == SpeciesType::NEUTRAL_MOLECULE; }
    bool is_positive_ion() const { return type == SpeciesType::POSITIVE_ION; }
    bool is_negative_ion() const { return type == SpeciesType::NEGATIVE_ION; }
};

// ═════════════════════════════════════════════════════════════
// Ratenkoeffizient (verschiedene Modelle)
// ═════════════════════════════════════════════════════════════
enum class RateType { CONSTANT, ARRHENIUS, TABULATED, POLYNOMIAL,
                      CTX_KIZ, CTX_KEL, CTX_KEX };

struct RateCoeff {
    RateType type = RateType::CONSTANT;
    double value = 0;         // Fuer CONSTANT
    double A = 0, E_a_eV = 0; // Fuer ARRHENIUS: K = A * exp(-E_a / Te)

    // Fuer TABULATED: Dateiname relativ zum Paketverzeichnis und die
    // geladene Stuetzstellenreihe. Ausserhalb des Tabellenbereichs wird auf
    // den Randwert festgehalten -- dasselbe Verhalten wie auf der Python-Seite.
    std::string table_file;
    std::vector<double> table_Te, table_K;

    // Fuer POLYNOMIAL: K = sum_i coeffs[i] * Te^i
    std::vector<double> poly_coeffs;

    // CTX_KIZ/KEL/KEX: nutzen die Ratenfunktionen aus rates.hpp und brauchen
    // deshalb den SimContext.

    double interp_table(double Te_eV) const {
        if (table_Te.empty()) return 0.0;
        if (Te_eV <= table_Te.front()) return table_K.front();
        if (Te_eV >= table_Te.back())  return table_K.back();
        size_t lo = 0, hi = table_Te.size() - 1;
        while (hi - lo > 1) {
            size_t m = (lo + hi) / 2;
            if (table_Te[m] <= Te_eV) lo = m; else hi = m;
        }
        double d = table_Te[hi] - table_Te[lo];
        if (d <= 0) return table_K[lo];
        double f = (Te_eV - table_Te[lo]) / d;
        return table_K[lo] * (1.0 - f) + table_K[hi] * f;
    }

    double eval_poly(double Te_eV) const {
        double acc = 0, pw = 1;
        for (double c : poly_coeffs) { acc += c * pw; pw *= Te_eV; }
        return acc;
    }

    double evaluate(double Te_eV) const {
        switch (type) {
            case RateType::CONSTANT:   return value;
            case RateType::ARRHENIUS:  return (Te_eV > 0) ? A * std::exp(-E_a_eV / Te_eV) : 0;
            case RateType::TABULATED:  return interp_table(Te_eV);
            case RateType::POLYNOMIAL: return eval_poly(Te_eV);
            default: return value;
        }
    }

    double evaluate(const SimContext& ctx, double Te_eV) const {
        switch (type) {
            case RateType::CTX_KIZ: return Kiz(ctx, Te_eV);
            case RateType::CTX_KEL: return Kel(ctx, Te_eV);
            case RateType::CTX_KEX: return Kex(ctx, Te_eV);
            default: return evaluate(Te_eV);
        }
    }
};

// ═════════════════════════════════════════════════════════════
// Reaktion
// ═════════════════════════════════════════════════════════════
struct ChemReaction {
    std::string id;              // Kennung aus dem Paket, z.B. "iz_Xe"
    std::string name;
    std::string type;            // Prozessart als Klartext, z.B. "ionization"
    std::map<std::string, int> reactants;  // {species_id: stoechiometrie}
    std::map<std::string, int> products;
    RateCoeff rate;
    double energy_eV = 0;        // Energieverlust pro Ereignis [eV]
    bool is_electron_impact = false;
    bool contributes_elastic = false;
    bool contributes_nu_m = false;  // Beitrag zur Stossfrequenz
    double surface_gamma = 0;    // Oberflaechenrekombinationskoeffizient

    // Netto-Stoechiometrie: products - reactants
    std::map<std::string, int> net_stoichiometry() const {
        std::map<std::string, int> net;
        for (auto& [id, s] : reactants) net[id] -= s;
        for (auto& [id, s] : products) net[id] += s;
        return net;
    }
};

// ═════════════════════════════════════════════════════════════
// Chemie-System (Spezies + Reaktionen)
// ═════════════════════════════════════════════════════════════
struct ChemSystem {
    std::string name;
    std::string description;
    std::string source_path;          // Woher das Paket geladen wurde, sonst leer
    double wall_temperature_K = 293.0;
    double sigma_i = 1e-18;           // Voreinstellung fuer alle Spezies
    std::vector<ChemSpecies> species;
    std::vector<ChemReaction> reactions;

    // Index-Lookup
    int species_index(const std::string& id) const {
        for (int i = 0; i < (int)species.size(); ++i)
            if (species[i].id == id) return i;
        return -1;
    }

    // Hilfsfunktionen
    std::vector<const ChemSpecies*> neutrals() const {
        std::vector<const ChemSpecies*> r;
        for (auto& s : species) if (s.is_neutral()) r.push_back(&s);
        return r;
    }
    std::vector<const ChemSpecies*> positive_ions() const {
        std::vector<const ChemSpecies*> r;
        for (auto& s : species) if (s.is_positive_ion()) r.push_back(&s);
        return r;
    }
    const ChemSpecies* feedstock() const {
        for (auto& s : species) if (s.is_feedstock) return &s;
        return nullptr;
    }

    int state_size() const { return (int)species.size() + 2; }  // N_species + Te + Tg
    int Te_idx() const { return (int)species.size(); }
    int Tg_idx() const { return (int)species.size() + 1; }
};

// ═════════════════════════════════════════════════════════════
// Stossfrequenz aus dem Reaktionsnetz
//
// Summiert ueber alle Reaktionen, die als Beitrag zur Stossfrequenz
// gekennzeichnet sind: nu_m = sum_r K_r(Te) * n(Stosspartner). Ist keine
// Reaktion gekennzeichnet, kommt -1 zurueck; der Aufrufer faellt dann auf
// den bisherigen Weg mit dem einen elastischen Ratenkoeffizienten zurueck.
// ═════════════════════════════════════════════════════════════

inline double chem_nu_m(const ChemSystem& sys, const std::vector<double>& state,
                        const SimContext& ctx, double Te) {
    bool markiert = false;
    double nu = 0;
    for (const auto& rxn : sys.reactions) {
        if (!rxn.contributes_nu_m) continue;
        markiert = true;
        double K = rxn.rate.evaluate(ctx, Te);
        if (K <= 0) continue;
        for (const auto& kv : rxn.reactants) {
            if (kv.first == "e") continue;
            int idx = sys.species_index(kv.first);
            if (idx >= 0) nu += K * state[idx];
        }
    }
    return markiert ? nu : -1.0;
}

// ═════════════════════════════════════════════════════════════
// Generischer Residual-Assembler
// ═════════════════════════════════════════════════════════════

inline double electron_density(const ChemSystem& sys, const std::vector<double>& state) {
    double ne = 0;
    for (int i = 0; i < (int)sys.species.size(); ++i) {
        if (sys.species[i].is_positive_ion())
            ne += std::abs(sys.species[i].charge) * state[i];
        else if (sys.species[i].is_negative_ion())
            ne -= std::abs(sys.species[i].charge) * state[i];
    }
    return std::max(ne, 0.0);
}

inline std::vector<double> assemble_residual(
    const ChemSystem& sys,
    const std::vector<double>& state,
    double P_abs_V,           // Absorbierte Leistungsdichte [W/m^3]
    double Q0_pps,            // Teilchenfluss [s^-1]
    const Thruster& geom,
    const SimContext& ctx,    // Fuer kontextabhaengige Raten (Kiz/Kel/Kex)
    double alpha_e_wall = 7.0,
    double density_profile_factor = 1.0)
{
    using namespace PhysConst;

    int N = sys.state_size();
    int Te_i = sys.Te_idx();
    int Tg_i = sys.Tg_idx();
    double Te = state[Te_i];
    double Tg = state[Tg_i];
    double V = geom.V;
    double ne = electron_density(sys, state);

    std::vector<double> resid(N, 0.0);

    // Gesamte Neutraldichte (fuer mittlere freie Weglaenge)
    double n_neut_total = 0;
    double sigma_avg = 1e-18;
    for (int i = 0; i < (int)sys.species.size(); ++i) {
        if (sys.species[i].is_neutral()) {
            n_neut_total += state[i];
            sigma_avg = sys.species[i].sigma_i;
        }
    }

    // ── 1. Volumenreaktionen ────────────────────────────────
    for (auto& rxn : sys.reactions) {
        if (rxn.surface_gamma > 0) continue;  // Oberflaechenreaktionen separat

        double K = rxn.rate.evaluate(ctx, Te);
        if (K <= 0) continue;

        // Reaktionsrate: K * Produkt(n_reactant)
        double rate = K;
        for (auto& [sp_id, stoech] : rxn.reactants) {
            if (sp_id == "e") {
                rate *= std::pow(ne, stoech);
            } else {
                int idx = sys.species_index(sp_id);
                if (idx >= 0) rate *= std::pow(state[idx], stoech);
            }
        }

        // Netto-Stoechiometrie -> Speziesbilanzen
        auto net = rxn.net_stoichiometry();
        for (auto& [sp_id, delta] : net) {
            if (sp_id == "e") continue;
            int idx = sys.species_index(sp_id);
            if (idx >= 0) resid[idx] += delta * rate;
        }

        // Energieverlust -> Elektronenenergiebilanz
        if (rxn.energy_eV > 0 && rxn.is_electron_impact) {
            double E_J = rxn.energy_eV * e;
            double n_eff = ne * density_profile_factor;
            double rate_e = K;
            for (auto& [sp_id, stoech] : rxn.reactants) {
                if (sp_id == "e") rate_e *= std::pow(n_eff, stoech);
                else {
                    int idx = sys.species_index(sp_id);
                    if (idx >= 0) rate_e *= std::pow(state[idx], stoech);
                }
            }
            resid[Te_i] -= E_J * rate_e;
        }

        // Elastische Heizung -> Gastemperatur
        if (rxn.contributes_elastic) {
            for (auto& [sp_id, stoech] : rxn.reactants) {
                if (sp_id == "e") continue;
                int idx = sys.species_index(sp_id);
                if (idx >= 0) {
                    double M_sp = sys.species[idx].mass_kg;
                    double Pg = 3.0 * me / M_sp * kB * (Te*conv - Tg) * ne * state[idx] * K;
                    resid[Tg_i] += Pg;
                }
            }
        }
    }

    // ── 2. Oberflaechenreaktionen ───────────────────────────
    for (auto& rxn : sys.reactions) {
        if (rxn.surface_gamma <= 0) continue;
        for (auto& [sp_id, stoech] : rxn.reactants) {
            if (sp_id == "e") continue;
            int idx = sys.species_index(sp_id);
            if (idx < 0) continue;
            double v_th = std::sqrt(8*kB*Tg / (pi*sys.species[idx].mass_kg));
            double surf_rate = rxn.surface_gamma * 0.25 * v_th * state[idx] * geom.A / V;
            resid[idx] -= surf_rate;
            for (auto& [pid, ps] : rxn.products) {
                int pi_ = sys.species_index(pid);
                if (pi_ >= 0) resid[pi_] += ps * surf_rate;
            }
        }
    }

    // ── 3. Feedstock-Zufuhr ─────────────────────────────────
    for (auto& sp : sys.species) {
        if (sp.is_feedstock) {
            int idx = sys.species_index(sp.id);
            if (idx >= 0) resid[idx] += Q0_pps / V;
        }
    }

    // ── 4. Ionenwandverluste (Bohm) ─────────────────────────
    double lam = 1.0 / std::max(n_neut_total * sigma_avg, 1e-10);
    double hL = 0.86 * std::pow(3.0 + geom.L/(2*lam), -0.5);
    double hR = 0.80 * std::pow(4.0 + geom.R/lam, -0.5);
    double Aeff = 2*hR*pi*geom.R*geom.L + 2*hL*pi*geom.R*geom.R;
    // Aeff1: fuer Neutralrueckstrom (ohne Grid-Anteil)
    double Aeff1 = 2*hR*pi*geom.R*geom.L + (2 - geom.betai)*hL*pi*geom.R*geom.R;

    for (int i = 0; i < (int)sys.species.size(); ++i) {
        auto& sp = sys.species[i];
        if (sp.is_positive_ion()) {
            double uBi = std::sqrt(kB*Te*conv / sp.mass_kg);
            double wall_loss = state[i] * uBi * Aeff / V;
            resid[i] -= wall_loss;

            // Neutralrueckstrom: Ionen die an der Wand neutralisiert werden
            // tragen ueber Aeff1 zur Neutralbilanz bei (Chabert r2 Term 2)
            for (auto& n_sp : sys.species) {
                if (n_sp.is_neutral() && std::abs(n_sp.mass_kg - sp.mass_kg)/sp.mass_kg < 0.01) {
                    int ni = sys.species_index(n_sp.id);
                    if (ni >= 0) resid[ni] += state[i] * uBi * Aeff1 / V;
                    break;
                }
            }

            // Ion-Neutral Stoss-Heizung (Chabert Pg2)
            double vi_val = std::sqrt(8*kB*Tg / (pi*sp.mass_kg));
            double Pg2 = 0.25 * sp.mass_kg * uBi * uBi * state[i] * n_neut_total * sigma_avg * vi_val;
            resid[Tg_i] += Pg2;

            // Elektronen-Wandverlust-Energie
            resid[Te_i] -= alpha_e_wall * kB * Te * conv * wall_loss;
        }
        else if (sp.is_neutral()) {
            // Neutralverlust durch Gitter
            double v_mean = std::sqrt(8*kB*Tg / (pi*sp.mass_kg));
            double outflow = 0.25 * state[i] * v_mean * geom.Ag / V;
            resid[i] -= outflow;
        }
    }

    // ── 5. RF-Leistungseintrag ──────────────────────────────
    resid[Te_i] += P_abs_V;

    // ── 6. Gaswaermeleitung ─────────────────────────────────
    double kappa_eff = 0;
    if (n_neut_total > 0) {
        for (auto& sp : sys.species) {
            if (sp.is_neutral() && sp.thermal_cond > 0) {
                int idx = sys.species_index(sp.id);
                if (idx >= 0) kappa_eff += sp.thermal_cond * (state[idx] / n_neut_total);
            }
        }
    }
    if (kappa_eff > 0) {
        double Pg_cond = kappa_eff * (Tg - Tg0) / geom.lambda_0 * geom.A / V;
        resid[Tg_i] -= Pg_cond;
    }

    return resid;
}

// ═════════════════════════════════════════════════════════════
// Xenon-Spezialfall: ChemSystem aus SimContext aufbauen
//
// Erster Nutzer der generischen Struktur. Nutzt die bestehenden
// Ratenkoeffizienten aus rates.hpp (Kiz, Kel, Kex) als Callbacks.
// ═════════════════════════════════════════════════════════════

inline ChemSystem build_xenon_system(const SimContext& ctx) {
    ChemSystem sys;
    sys.name = "Xenon (generic)";

    // Spezies
    sys.species.push_back({"Xe", SpeciesType::NEUTRAL_ATOM, ctx.gas.M, 0,
                           ctx.gas.kappa, true, ctx.gas.sigma_i});
    sys.species.push_back({"Xe+", SpeciesType::POSITIVE_ION, ctx.gas.M, +1,
                           0, false, ctx.gas.sigma_i});

    // Ionisation: e + Xe -> 2e + Xe+
    ChemReaction ioniz;
    ioniz.name = "ionization";
    ioniz.reactants = {{"e", 1}, {"Xe", 1}};
    ioniz.products = {{"e", 2}, {"Xe+", 1}};
    ioniz.energy_eV = ctx.gas.Eiz / PhysConst::e;  // J -> eV
    ioniz.is_electron_impact = true;
    ioniz.rate.type = RateType::CTX_KIZ;
    sys.reactions.push_back(ioniz);

    // Anregung: e + Xe -> e + Xe* (Energieverlust, keine Speziesaenderung)
    ChemReaction excit;
    excit.name = "excitation";
    excit.reactants = {{"e", 1}, {"Xe", 1}};
    excit.products = {{"e", 1}, {"Xe", 1}};
    excit.energy_eV = ctx.gas.Eexc / PhysConst::e;
    excit.is_electron_impact = true;
    excit.rate.type = RateType::CTX_KEX;
    sys.reactions.push_back(excit);

    // Elastisch: e + Xe -> e + Xe (Gasheizung)
    ChemReaction elastic;
    elastic.name = "elastic";
    elastic.reactants = {{"e", 1}, {"Xe", 1}};
    elastic.products = {{"e", 1}, {"Xe", 1}};
    elastic.energy_eV = 0;
    elastic.contributes_elastic = true;
    elastic.contributes_nu_m = true;
    elastic.rate.type = RateType::CTX_KEL;
    sys.reactions.push_back(elastic);

    return sys;
}

#endif // CHEM_SYSTEM_HPP
