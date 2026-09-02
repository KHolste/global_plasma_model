// sim_config.cpp -- Config-Laden und Anwendung auf SimContext.
//
// Unterstuetzt zwei Formate:
//   1. run_config.json (primaer, strukturiertes JSON)
//   2. params.txt (Legacy, flaches Key-Value)
//
// Auto-Detection: Wenn die Eingabedatei mit '{' beginnt, wird JSON geparst.
#include "sim_config.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <set>
#include <algorithm>
#include <cstdlib>

// ═════════════════════════════════════════════════════════════
// Minimaler JSON-Value-Extraktor (kein vollstaendiger Parser)
//
// Extrahiert Werte aus dem flachen RunConfig-JSON-Format:
//   {"geometry": {"R": 0.05, "L": 0.064}, "grid": {"Vgrid": 1500}, ...}
// Werte werden als Key-Value-Paare in ConfigData.numeric/strings abgelegt,
// wobei der Sektionsprefix wegfaellt (R statt geometry.R).
// ═════════════════════════════════════════════════════════════

static ConfigData parseJsonConfig(const std::string& text) {
    ConfigData cd;
    // Felder, die als String interpretiert werden
    static const std::set<std::string> STR_KEYS = {
        "backend", "package_id", "gas", "cs_database", "preset_id"
    };
    // Mapping: JSON-Feld -> params.txt-aequivalenter Schluessel
    // (nur wo noetig; die meisten Felder werden direkt uebernommen)
    static const std::map<std::string, std::string> REMAP = {
        {"P_max", "P_RFG_max"},
        {"Q0_start", "Q0sccm_start"},
        {"Q0_step", "Q0sccm_step"},
        {"gas", "gas_species"},
    };

    // Einfache Extraktion: suche nach "key": value Paare
    // Funktioniert fuer das bekannte RunConfig-Format (2 Ebenen, keine Arrays)
    size_t pos = 0;
    while (pos < text.size()) {
        // Finde naechsten String-Schluessel
        size_t q1 = text.find('"', pos);
        if (q1 == std::string::npos) break;
        size_t q2 = text.find('"', q1 + 1);
        if (q2 == std::string::npos) break;
        std::string key = text.substr(q1 + 1, q2 - q1 - 1);
        pos = q2 + 1;

        // Finde ':'
        size_t colon = text.find(':', pos);
        if (colon == std::string::npos) break;
        pos = colon + 1;

        // Ueberspringe Whitespace
        while (pos < text.size() && (text[pos] == ' ' || text[pos] == '\n' || text[pos] == '\r' || text[pos] == '\t'))
            pos++;

        if (pos >= text.size()) break;

        // Wert: { (Sektion) oder Skalar
        if (text[pos] == '{') {
            // Sektion: ueberspringe, Felder werden auf naechster Ebene gefunden
            pos++;
            continue;
        }

        // Ueberspringe Sektions-Schluessel (geometry, grid, coil, operation, sweep, rates, meta)
        if (key == "geometry" || key == "grid" || key == "coil" || key == "operation" ||
            key == "sweep" || key == "rates" || key == "meta") {
            continue;
        }

        // Skalarwert extrahieren
        if (text[pos] == '"') {
            // String-Wert
            size_t sq1 = pos;
            size_t sq2 = text.find('"', sq1 + 1);
            if (sq2 != std::string::npos) {
                std::string val = text.substr(sq1 + 1, sq2 - sq1 - 1);
                std::string mapped = REMAP.count(key) ? REMAP.at(key) : key;
                cd.strings[mapped] = val;
                pos = sq2 + 1;
            }
        } else {
            // Numerischer Wert (oder bool/null — ignorieren wir)
            size_t end = text.find_first_of(",}\n\r", pos);
            if (end == std::string::npos) end = text.size();
            std::string valstr = text.substr(pos, end - pos);
            // Trim
            while (!valstr.empty() && (valstr.back() == ' ' || valstr.back() == '\t'))
                valstr.pop_back();

            if (valstr == "true" || valstr == "false" || valstr == "null") {
                pos = end;
                continue;
            }

            std::string mapped = REMAP.count(key) ? REMAP.at(key) : key;

            // Versuche als double
            try {
                size_t idx;
                double d = std::stod(valstr, &idx);
                if (idx > 0) {
                    if (STR_KEYS.count(key)) {
                        cd.strings[mapped] = valstr;
                    } else {
                        cd.numeric[mapped] = d;
                    }
                }
            } catch (...) {
                // Kein numerischer Wert
            }
            pos = end;
        }
    }
    return cd;
}

// ═════════════════════════════════════════════════════════════
// Legacy params.txt Parser
// ═════════════════════════════════════════════════════════════

static ConfigData parseParamsTxt(const std::string& text) {
    ConfigData cd;
    static const std::set<std::string> STR = {"gas_species", "cs_database", "chemistry_package"};
    std::istringstream stream(text);
    std::string line;
    while (std::getline(stream, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string key;
        if (!(ss >> key)) continue;
        if (STR.count(key)) {
            std::string v; if (ss >> v) cd.strings[key] = v;
        } else {
            double v; if (ss >> v) cd.numeric[key] = v;
        }
    }
    return cd;
}

// ═════════════════════════════════════════════════════════════
// Haupt-Ladefunktion mit Auto-Detection
// ═════════════════════════════════════════════════════════════

ConfigData loadConfig(const std::string& filename) {
    // 1. Wenn explizit eine .json-Datei angegeben: direkt JSON lesen
    bool is_json = (filename.size() > 5 && filename.substr(filename.size() - 5) == ".json");

    // HINWEIS: Kein automatisches Laden von run_config.json mehr.
    // Stale JSON-Dateien koennen sonst den C++-Solver korrumpieren.
    // JSON wird nur geladen wenn explizit als Argument uebergeben.

    // Datei oeffnen
    std::ifstream f(filename);
    if (!f) {
        std::cerr << "[Config] Datei '" << filename << "' nicht gefunden." << std::endl;
        return ConfigData{};
    }

    std::string text((std::istreambuf_iterator<char>(f)),
                      std::istreambuf_iterator<char>());

    // Auto-Detect: erstes Non-Whitespace-Zeichen
    size_t start = text.find_first_not_of(" \t\n\r");
    if (start != std::string::npos && text[start] == '{') {
        return parseJsonConfig(text);
    }
    return parseParamsTxt(text);
}

// ═════════════════════════════════════════════════════════════
// Anwendung auf SimContext (unveraenderte Semantik)
// ═════════════════════════════════════════════════════════════

void applyConfig(SimContext& ctx, const ConfigData& cd) {
    const auto& c = cd.numeric;
    auto& g = ctx.gas;
    auto& t = ctx.thruster;
    auto& r = ctx.rates;
    auto& s = ctx.solver;

    // Gas
    if (cd.strings.count("gas_species")) {
        std::string gs = cd.strings.at("gas_species");
        if (!g.set_species(gs))
            std::cerr << "FEHLER: Unbekanntes Gas '" << gs << "'" << std::endl;
    }
    if (cd.strings.count("cs_database")) ctx.cs_database = cd.strings.at("cs_database");
    if (cd.strings.count("chemistry_package")) ctx.chem_package = cd.strings.at("chemistry_package");
    if (cd.strings.count("preset_id")) ctx.meta_preset_id = cd.strings.at("preset_id");

    // Thruster
    if (c.count("R"))         t.R = c.at("R");
    if (c.count("L"))         t.L = c.at("L");
    if (c.count("betai"))     t.betai = c.at("betai");
    if (c.count("betag"))     t.betag = c.at("betag");
    if (c.count("frequency")) t.frequency = c.at("frequency");
    if (c.count("Nw"))        t.Nw = c.at("Nw");
    if (c.count("R_ohm"))     t.R_ohm = c.at("R_ohm");
    if (c.count("Rc"))        t.Rc = c.at("Rc");
    if (c.count("lc"))        t.lc = c.at("lc");
    if (c.count("Vgrid"))     t.Vgrid = c.at("Vgrid");
    if (c.count("sgrid"))     t.sgrid = c.at("sgrid");
    if (c.count("eta_opt"))   t.eta_opt = c.at("eta_opt");  // NEU: C++ liest jetzt eta_opt
    if (c.count("P_RFG"))     t.P_RFG = c.at("P_RFG");
    if (c.count("Q0sccm"))    t.Q0sccm = c.at("Q0sccm");
    ctx.recompute();

    // Solver
    if (c.count("solve_mode"))  s.solve_mode = std::max(1, std::min(2, (int)c.at("solve_mode")));
    if (c.count("I_soll"))      s.I_soll = c.at("I_soll");
    if (c.count("Q0sccm_start")) s.Q0sccm_start = c.at("Q0sccm_start");
    if (c.count("Q0sccm_step"))  s.Q0sccm_step = c.at("Q0sccm_step");
    if (c.count("jjmax"))       s.jjmax = (int)c.at("jjmax");
    if (c.count("N"))           s.jjmax = (int)c.at("N");  // RunConfig-Name
    if (c.count("P_RFG_max"))   s.P_RFG_max = c.at("P_RFG_max");
    if (c.count("P_abs_scale")) s.P_abs_scale = c.at("P_abs_scale");
    if (c.count("density_profile_factor")) s.density_profile_factor = c.at("density_profile_factor");
    if (c.count("alpha_e_wall")) s.alpha_e_wall = c.at("alpha_e_wall");
    if (c.count("wall_energy_model")) s.wall_energy_model = (int)c.at("wall_energy_model");
    if (c.count("newton_tol"))  s.newton_tol = c.at("newton_tol");
    if (c.count("newton_max_iter")) s.newton_max_iter = (int)c.at("newton_max_iter");
    if (c.count("newton_max_log_step")) s.newton_max_log_step = c.at("newton_max_log_step");
    if (c.count("newton_fd_eps")) s.newton_fd_eps = c.at("newton_fd_eps");
    if (c.count("power_tol_mA")) s.power_tol_mA = c.at("power_tol_mA");
    if (c.count("power_max_iter")) s.power_max_iter = (int)c.at("power_max_iter");
    if (c.count("power_min"))   s.power_min = c.at("power_min");
    if (c.count("n_min"))  s.n_min  = c.at("n_min");
    if (c.count("n_max"))  s.n_max  = c.at("n_max");
    if (c.count("ng_min")) s.ng_min = c.at("ng_min");
    if (c.count("ng_max")) s.ng_max = c.at("ng_max");
    if (c.count("Te_min")) s.Te_min = c.at("Te_min");
    if (c.count("Te_max")) s.Te_max = c.at("Te_max");
    if (c.count("Tg_min")) s.Tg_min = c.at("Tg_min");
    if (c.count("Tg_max")) s.Tg_max = c.at("Tg_max");

    // Ratenmodell
    if (c.count("rate_model")) {
        ctx.rate_model = (int)c.at("rate_model");
        if (ctx.rate_model == 1)      { r.ionization_model=1; r.excitation_model=1; r.elastic_model=0; }
        else if (ctx.rate_model == 2) { r.ionization_model=1; r.excitation_model=1; r.elastic_model=1; }
        else                          { r.ionization_model=0; r.excitation_model=0; r.elastic_model=0; }
    }
    if (c.count("ionization_model")) r.ionization_model = (int)c.at("ionization_model");
    if (c.count("elastic_model"))    r.elastic_model = (int)c.at("elastic_model");
    if (c.count("excitation_model")) r.excitation_model = (int)c.at("excitation_model");
    if (c.count("kel_constant"))   { r.kel_constant = c.at("kel_constant");
                                     r.kel_constant_explicit = true; }
    if (c.count("allow_foreign_rate_fits"))
        ctx.allow_foreign_rate_fits = c.at("allow_foreign_rate_fits") != 0.0;
    if (c.count("Kex_scale"))        r.Kex_scale = c.at("Kex_scale");
    if (c.count("use_paper_kel") && (int)c.at("use_paper_kel") != 0)
        r.elastic_model = 0;

    // Debug
    if (c.count("debug_level")) ctx.debug_level = std::max(0, std::min(3, (int)c.at("debug_level")));
}
