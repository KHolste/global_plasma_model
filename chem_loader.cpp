// chem_loader.cpp -- Einlesen der Chemiepakete. Einzige Uebersetzungseinheit,
// die die JSON-Bibliothek einbindet.
#include "chem_loader.hpp"

#include "third_party/json.hpp"

#include <fstream>
#include <sstream>
#include <set>
#include <map>
#include <cmath>

using json = nlohmann::json;

// ═════════════════════════════════════════════════════════════
// Hilfsfunktionen
// ═════════════════════════════════════════════════════════════

namespace {

std::string parent_dir(const std::string& path) {
    size_t p = path.find_last_of("/\\");
    return (p == std::string::npos) ? std::string(".") : path.substr(0, p);
}

std::string join_path(const std::string& dir, const std::string& rel) {
    if (rel.empty()) return dir;
    if (dir.empty() || dir == ".") return rel;
    return dir + "/" + rel;
}

// Wert lesen, wenn vorhanden, sonst die Vorgabe. Ein null-Eintrag zaehlt als
// nicht vorhanden -- die Pakete verwenden ihn stellenweise so.
double get_double(const json& j, const char* key, double fallback) {
    auto it = j.find(key);
    if (it == j.end() || it->is_null()) return fallback;
    return it->get<double>();
}

bool get_bool(const json& j, const char* key, bool fallback) {
    auto it = j.find(key);
    if (it == j.end() || it->is_null()) return fallback;
    return it->get<bool>();
}

std::string get_string(const json& j, const char* key, const std::string& fallback) {
    auto it = j.find(key);
    if (it == j.end() || it->is_null()) return fallback;
    return it->get<std::string>();
}

std::map<std::string, int> get_stoich(const json& j, const char* key) {
    std::map<std::string, int> out;
    auto it = j.find(key);
    if (it == j.end() || !it->is_object()) return out;
    for (auto e = it->begin(); e != it->end(); ++e)
        out[e.key()] = e.value().get<int>();
    return out;
}

// Ratentabelle lesen: zwei Spalten, Komma getrennt, Kopf- und Kommentarzeilen
// werden uebergangen. Dasselbe Format, das die Python-Seite liest.
bool load_rate_table(const std::string& path,
                     std::vector<double>& Te, std::vector<double>& K) {
    std::ifstream f(path);
    if (!f) return false;
    Te.clear(); K.clear();
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == 'T') continue;
        std::istringstream ss(line);
        std::string a, b;
        if (!std::getline(ss, a, ',')) continue;
        if (!std::getline(ss, b, ',')) continue;
        try {
            double te = std::stod(a), k = std::stod(b);
            Te.push_back(te); K.push_back(k);
        } catch (const std::exception&) { continue; }
    }
    return !Te.empty();
}

bool parse_species_type(const std::string& s, SpeciesType& out) {
    if (s == "neutral_atom")     { out = SpeciesType::NEUTRAL_ATOM;     return true; }
    if (s == "neutral_molecule") { out = SpeciesType::NEUTRAL_MOLECULE; return true; }
    if (s == "positive_ion")     { out = SpeciesType::POSITIVE_ION;     return true; }
    if (s == "negative_ion")     { out = SpeciesType::NEGATIVE_ION;     return true; }
    return false;
}

}  // namespace

// ═════════════════════════════════════════════════════════════
// Fehlerausgabe
// ═════════════════════════════════════════════════════════════

std::string ChemLoadResult::error_text() const {
    std::ostringstream os;
    for (size_t i = 0; i < errors.size(); ++i) {
        if (i) os << "\n";
        os << "  - " << errors[i];
    }
    return os.str();
}

// ═════════════════════════════════════════════════════════════
// Pfadaufloesung
// ═════════════════════════════════════════════════════════════

std::string resolve_chem_package(const std::string& name_or_path) {
    if (name_or_path.size() > 5 &&
        name_or_path.compare(name_or_path.size() - 5, 5, ".json") == 0)
        return name_or_path;
    return "chemistry/" + name_or_path + "/chemistry.json";
}

// ═════════════════════════════════════════════════════════════
// Laden
// ═════════════════════════════════════════════════════════════

ChemLoadResult load_chem_package(const std::string& json_path) {
    ChemLoadResult res;
    ChemSystem& sys = res.system;
    sys.source_path = json_path;

    std::ifstream f(json_path);
    if (!f) {
        res.errors.push_back("Datei nicht lesbar: " + json_path);
        return res;
    }

    json doc;
    try {
        f >> doc;
    } catch (const std::exception& ex) {
        res.errors.push_back(std::string("JSON nicht lesbar: ") + ex.what());
        return res;
    }
    if (!doc.is_object()) {
        res.errors.push_back("Wurzelelement ist kein Objekt");
        return res;
    }

    const std::string base = parent_dir(json_path);

    sys.name = get_string(doc, "name", "");
    sys.description = get_string(doc, "description", "");
    sys.wall_temperature_K = get_double(doc, "wall_temperature_K", 293.0);
    sys.sigma_i = get_double(doc, "sigma_i", 1e-18);
    if (sys.name.empty()) res.errors.push_back("Feld 'name' fehlt");

    // ── Spezies ──────────────────────────────────────────────
    auto sp_it = doc.find("species");
    if (sp_it == doc.end() || !sp_it->is_array()) {
        res.errors.push_back("Feld 'species' fehlt oder ist keine Liste");
        return res;
    }
    std::set<std::string> known_ids;
    for (const auto& sd : *sp_it) {
        std::string id = get_string(sd, "id", "");
        std::string type_str = get_string(sd, "type", "");
        if (id.empty()) { res.errors.push_back("Spezies ohne 'id'"); continue; }

        if (type_str == "electron") {
            // Nicht als Zustandsvariable fuehren -- die Elektronendichte folgt
            // aus der Quasineutralitaet.
            ++res.electron_entries;
            known_ids.insert(id);
            continue;
        }

        SpeciesType st;
        if (!parse_species_type(type_str, st)) {
            res.errors.push_back("Spezies '" + id + "': unbekannter Typ '" + type_str + "'");
            continue;
        }
        if (sd.find("mass_kg") == sd.end()) {
            res.errors.push_back("Spezies '" + id + "': 'mass_kg' fehlt");
            continue;
        }

        ChemSpecies sp;
        sp.id = id;
        sp.type = st;
        sp.mass_kg = get_double(sd, "mass_kg", 0.0);
        sp.charge = (int)get_double(sd, "charge", 0.0);
        sp.thermal_cond = get_double(sd, "thermal_conductivity", 0.0);
        sp.is_feedstock = get_bool(sd, "is_feedstock", false);
        sp.sigma_i = get_double(sd, "sigma_i", sys.sigma_i);
        sp.name = get_string(sd, "name", id);
        sp.is_beam_extracted = get_bool(sd, "is_beam_extracted", false);
        sp.wall_products = get_stoich(sd, "wall_products");

        if (sp.mass_kg <= 0)
            res.errors.push_back("Spezies '" + id + "': Masse ist nicht positiv");
        if (st == SpeciesType::POSITIVE_ION && sp.charge <= 0)
            res.errors.push_back("Spezies '" + id + "': positives Ion mit Ladung " +
                                 std::to_string(sp.charge));
        if (st == SpeciesType::NEGATIVE_ION && sp.charge >= 0)
            res.errors.push_back("Spezies '" + id + "': negatives Ion mit Ladung " +
                                 std::to_string(sp.charge));

        known_ids.insert(id);
        sys.species.push_back(sp);
    }

    // ── Wandprodukte der Ionen ───────────────────────────────
    // Was ein Ion an der Wand hinterlaesst, gehoert in das Paket. Fehlt die
    // Angabe, wird sie aus der Masse abgeleitet -- das trifft den einatomigen
    // Fall, aber nicht den molekularen, und wird deshalb gezaehlt und
    // gemeldet. Die Massenbilanz wird in jedem Fall geprueft.
    for (auto& sp : sys.species) {
        if (!sp.is_positive_ion()) continue;

        if (sp.wall_products.empty()) {
            for (const auto& n : sys.species) {
                if (!n.is_neutral()) continue;
                if (std::fabs(n.mass_kg - sp.mass_kg) <= 0.01 * sp.mass_kg) {
                    sp.wall_products[n.id] = 1;
                    ++res.wall_products_derived;
                    break;
                }
            }
            if (sp.wall_products.empty()) {
                res.errors.push_back("Spezies '" + sp.id + "': 'wall_products' fehlt "
                                     "und laesst sich nicht aus der Masse ableiten");
                continue;
            }
        }

        double masse = 0;
        bool brauchbar = true;
        for (const auto& wp : sp.wall_products) {
            int idx = sys.species_index(wp.first);
            if (idx < 0) {
                res.errors.push_back("Spezies '" + sp.id + "': Wandprodukt '" +
                                     wp.first + "' ist nicht definiert");
                brauchbar = false;
                break;
            }
            if (!sys.species[idx].is_neutral()) {
                res.errors.push_back("Spezies '" + sp.id + "': Wandprodukt '" +
                                     wp.first + "' ist nicht neutral");
                brauchbar = false;
                break;
            }
            if (wp.second <= 0) {
                res.errors.push_back("Spezies '" + sp.id + "': Wandprodukt '" +
                                     wp.first + "' mit Anzahl " +
                                     std::to_string(wp.second));
                brauchbar = false;
                break;
            }
            masse += wp.second * sys.species[idx].mass_kg;
        }
        if (brauchbar && std::fabs(masse - sp.mass_kg) > 0.01 * sp.mass_kg)
            res.errors.push_back("Spezies '" + sp.id + "': die Wandprodukte wiegen "
                                 "nicht so viel wie das Ion");
    }

    // ── Reaktionen ───────────────────────────────────────────
    auto rx_it = doc.find("reactions");
    if (rx_it == doc.end() || !rx_it->is_array()) {
        res.errors.push_back("Feld 'reactions' fehlt oder ist keine Liste");
        return res;
    }
    for (const auto& rd : *rx_it) {
        ChemReaction rxn;
        rxn.id = get_string(rd, "id", "");
        rxn.name = get_string(rd, "name", rxn.id);
        rxn.type = get_string(rd, "type", "");
        if (rxn.id.empty()) { res.errors.push_back("Reaktion ohne 'id'"); continue; }

        rxn.reactants = get_stoich(rd, "reactants");
        rxn.products = get_stoich(rd, "products");
        rxn.energy_eV = get_double(rd, "energy_eV", 0.0);
        rxn.is_electron_impact = get_bool(rd, "is_electron_impact", true);
        rxn.contributes_elastic = get_bool(rd, "elastic_heating", false);
        rxn.contributes_nu_m = get_bool(rd, "nu_m", false);
        rxn.surface_gamma = get_double(rd, "surface_gamma", 0.0);

        // Rate
        auto rate_it = rd.find("rate");
        if (rate_it == rd.end() || !rate_it->is_object()) {
            res.errors.push_back("Reaktion '" + rxn.id + "': 'rate' fehlt");
            continue;
        }
        const json& jr = *rate_it;
        std::string model = get_string(jr, "model", "");
        if (model == "constant") {
            rxn.rate.type = RateType::CONSTANT;
            rxn.rate.value = get_double(jr, "value", 0.0);
        } else if (model == "arrhenius") {
            rxn.rate.type = RateType::ARRHENIUS;
            rxn.rate.A = get_double(jr, "A", 0.0);
            rxn.rate.E_a_eV = get_double(jr, "E_a_eV", 0.0);
        } else if (model == "polynomial") {
            rxn.rate.type = RateType::POLYNOMIAL;
            auto c_it = jr.find("coeffs");
            if (c_it == jr.end() || !c_it->is_array()) {
                res.errors.push_back("Reaktion '" + rxn.id + "': 'coeffs' fehlt");
                continue;
            }
            for (const auto& c : *c_it) rxn.rate.poly_coeffs.push_back(c.get<double>());
        } else if (model == "tabulated") {
            rxn.rate.type = RateType::TABULATED;
            rxn.rate.table_file = get_string(jr, "file", "");
            if (rxn.rate.table_file.empty()) {
                res.errors.push_back("Reaktion '" + rxn.id + "': Tabellendatei fehlt");
                continue;
            }
            const std::string tp = join_path(base, rxn.rate.table_file);
            if (load_rate_table(tp, rxn.rate.table_Te, rxn.rate.table_K))
                ++res.tables_loaded;
            else
                res.errors.push_back("Reaktion '" + rxn.id + "': Tabelle nicht lesbar: " + tp);
        } else {
            res.errors.push_back("Reaktion '" + rxn.id + "': unbekanntes Ratenmodell '" + model + "'");
            continue;
        }

        sys.reactions.push_back(rxn);
    }

    // ── Pruefung, gleiche Regeln wie auf der Python-Seite ─────
    if (res.electron_entries != 1)
        res.errors.push_back("Genau 1 Elektron-Spezies erwartet, " +
                             std::to_string(res.electron_entries) + " gefunden");

    int n_neutral = 0, n_pos = 0, n_feed = 0;
    for (const auto& sp : sys.species) {
        if (sp.is_neutral()) ++n_neutral;
        if (sp.is_positive_ion()) ++n_pos;
        if (sp.is_feedstock) ++n_feed;
    }
    if (n_neutral == 0) res.errors.push_back("Mindestens eine neutrale Spezies erforderlich");
    if (n_pos == 0)     res.errors.push_back("Mindestens ein positives Ion erforderlich");
    if (n_feed == 0)    res.errors.push_back("Mindestens eine Feedstock-Spezies erforderlich");

    for (const auto& rxn : sys.reactions) {
        for (const auto* side : {&rxn.reactants, &rxn.products}) {
            for (const auto& kv : *side) {
                if (kv.first == "e") continue;
                if (!known_ids.count(kv.first))
                    res.errors.push_back("Reaktion '" + rxn.id + "': Spezies '" +
                                         kv.first + "' nicht definiert");
            }
        }
    }

    res.ok = res.errors.empty();
    return res;
}
