// test_chem_loader.cpp -- Test fuer das Einlesen der Chemiepakete.
//
// Prueft, dass der C++-Kern dieselben Pakete liest wie die Python-Seite:
// Spezies mit Typ, Masse und Ladung, Reaktionen mit Stoechiometrie, sowie
// die vier Ratenmodelle einschliesslich der tabellierten Raten aus CSV.
// Geprueft wird ausserdem, dass fehlerhafte Pakete abgelehnt statt still
// halb geladen werden.
//
// Zweiter Betriebsmodus: mit "--dump <paket.json>" gibt das Programm den
// geladenen Inhalt maschinenlesbar aus. Darauf setzt der Vergleich gegen den
// Python-Lader auf (test_chem_loader_cross.py).
//
// Build:
//   g++ -O3 -std=c++17 test_chem_loader.cpp sim_config.o rates.o physics.o \
//       solver.o sim_logging.o bessel_wrapper.o chem_loader.o -o test_chem_loader
#include "chem_loader.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <string>

using namespace std;

static int passed = 0, failed = 0;

static void check(const string& name, bool cond, const string& detail = "") {
    if (cond) { ++passed; cout << "  PASS: " << name << endl; }
    else { ++failed; cout << "  FAIL: " << name << " -- " << detail << endl; }
}

static const ChemSpecies* find_species(const ChemSystem& sys, const string& id) {
    int i = sys.species_index(id);
    return (i >= 0) ? &sys.species[i] : nullptr;
}

static const ChemReaction* find_reaction(const ChemSystem& sys, const string& id) {
    for (const auto& r : sys.reactions) if (r.id == id) return &r;
    return nullptr;
}

static string type_name(SpeciesType t) {
    switch (t) {
        case SpeciesType::NEUTRAL_ATOM:     return "neutral_atom";
        case SpeciesType::NEUTRAL_MOLECULE: return "neutral_molecule";
        case SpeciesType::POSITIVE_ION:     return "positive_ion";
        case SpeciesType::NEGATIVE_ION:     return "negative_ion";
    }
    return "?";
}

static string rate_name(RateType t) {
    switch (t) {
        case RateType::CONSTANT:   return "constant";
        case RateType::ARRHENIUS:  return "arrhenius";
        case RateType::TABULATED:  return "tabulated";
        case RateType::POLYNOMIAL: return "polynomial";
        default: return "context";
    }
}

static string stoich_text(const map<string, int>& m) {
    string s;
    for (const auto& kv : m) {
        if (!s.empty()) s += ",";
        s += kv.first + ":" + to_string(kv.second);
    }
    return s;
}

// ═══ Maschinenlesbare Ausgabe fuer den Vergleich mit Python ═══
static int dump(const string& path) {
    ChemLoadResult r = load_chem_package(path);
    if (!r.ok) {
        cerr << "Paket nicht ladbar:\n" << r.error_text() << endl;
        return 1;
    }
    const ChemSystem& sys = r.system;
    cout << setprecision(12) << scientific;
    cout << "PKG\t" << sys.name << "\t" << sys.wall_temperature_K << "\t" << sys.sigma_i << "\n";
    for (const auto& sp : sys.species)
        cout << "SPECIES\t" << sp.id << "\t" << type_name(sp.type) << "\t"
             << sp.mass_kg << "\t" << sp.charge << "\t" << sp.thermal_cond << "\t"
             << (sp.is_feedstock ? 1 : 0) << "\t" << (sp.is_beam_extracted ? 1 : 0) << "\t"
             << stoich_text(sp.wall_products) << "\n";
    for (const auto& rx : sys.reactions) {
        cout << "REACTION\t" << rx.id << "\t" << rx.type << "\t"
             << stoich_text(rx.reactants) << "\t" << stoich_text(rx.products) << "\t"
             << rx.energy_eV << "\t" << (rx.is_electron_impact ? 1 : 0) << "\t"
             << (rx.contributes_elastic ? 1 : 0) << "\t" << (rx.contributes_nu_m ? 1 : 0) << "\t"
             << rx.surface_gamma << "\t" << rate_name(rx.rate.type);
        for (double Te : {0.7, 1.0, 2.0, 3.75, 5.0, 8.0, 12.0, 30.0})
            cout << "\t" << rx.rate.evaluate(Te);
        cout << "\n";
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc >= 3 && string(argv[1]) == "--dump") return dump(argv[2]);

    // ── Test 1: Xenon-Paket ──────────────────────────────────
    cout << "\n--- Test 1: xenon_simple laden ---" << endl;
    ChemLoadResult xe = load_chem_package(resolve_chem_package("xenon_simple"));
    check("Paket geladen", xe.ok, xe.error_text());
    if (!xe.ok) { cout << "\nAbbruch: Grundpaket fehlt." << endl; return 1; }

    const ChemSystem& xs = xe.system;
    check("genau ein Elektroneneintrag", xe.electron_entries == 1, to_string(xe.electron_entries));
    check("Elektron nicht als Zustandsvariable", xs.species_index("e") < 0);
    check("2 schwere Spezies", xs.species.size() == 2, to_string(xs.species.size()));
    check("3 Reaktionen", xs.reactions.size() == 3, to_string(xs.reactions.size()));
    check("Zustandslaenge 4", xs.state_size() == 4, to_string(xs.state_size()));

    const ChemSpecies* xe_n = find_species(xs, "Xe");
    const ChemSpecies* xe_i = find_species(xs, "Xe+");
    check("Xe vorhanden", xe_n != nullptr);
    check("Xe+ vorhanden", xe_i != nullptr);
    if (xe_n && xe_i) {
        check("Xe ist neutrales Atom", xe_n->type == SpeciesType::NEUTRAL_ATOM);
        check("Xe ist Feedstock", xe_n->is_feedstock);
        check("Xe-Masse", fabs(xe_n->mass_kg - 2.1801711e-25) < 1e-31, to_string(xe_n->mass_kg));
        check("Xe-Waermeleitfaehigkeit", fabs(xe_n->thermal_cond - 0.0057) < 1e-9,
              to_string(xe_n->thermal_cond));
        check("Xe+ ist positives Ion", xe_i->type == SpeciesType::POSITIVE_ION);
        check("Xe+ Ladung 1", xe_i->charge == 1, to_string(xe_i->charge));
        check("Xe+ wird extrahiert", xe_i->is_beam_extracted);
        check("Reihenfolge wie in der Datei",
              xs.species_index("Xe") == 0 && xs.species_index("Xe+") == 1);
    }
    check("Stossquerschnitt aus dem Paket", fabs(xs.sigma_i - 1e-18) < 1e-24, to_string(xs.sigma_i));

    // ── Test 2: Reaktionen und Arrhenius-Rate ────────────────
    cout << "\n--- Test 2: Reaktionen und Raten ---" << endl;
    const ChemReaction* iz = find_reaction(xs, "iz_Xe");
    const ChemReaction* el = find_reaction(xs, "el_Xe");
    check("Ionisation vorhanden", iz != nullptr);
    check("elastische Reaktion vorhanden", el != nullptr);
    if (iz) {
        check("Ionisation ist Elektronenstoss", iz->is_electron_impact);
        check("Ionisation Edukte e+Xe", stoich_text(iz->reactants) == "Xe:1,e:1",
              stoich_text(iz->reactants));
        check("Ionisation Produkte 2e+Xe+", stoich_text(iz->products) == "Xe+:1,e:2",
              stoich_text(iz->products));
        check("Ionisationsenergie 12.127 eV", fabs(iz->energy_eV - 12.127) < 1e-9,
              to_string(iz->energy_eV));
        check("Arrhenius-Modell", iz->rate.type == RateType::ARRHENIUS);
        double Te = 3.75, erwartet = 1.8e-13 * exp(-12.127 / Te);
        double ist = iz->rate.evaluate(Te);
        check("Arrhenius-Wert", fabs(ist - erwartet) < 1e-9 * fabs(erwartet),
              to_string(ist) + " statt " + to_string(erwartet));
        // Netto-Stoechiometrie: ein Xe verschwindet, ein Xe+ entsteht
        auto net = iz->net_stoichiometry();
        check("Netto-Bilanz Xe", net["Xe"] == -1, to_string(net["Xe"]));
        check("Netto-Bilanz Xe+", net["Xe+"] == 1, to_string(net["Xe+"]));
    }
    if (el) {
        check("elastisch heizt das Gas", el->contributes_elastic);
        check("elastisch traegt zur Stossfrequenz bei", el->contributes_nu_m);
        check("konstante Rate", el->rate.type == RateType::CONSTANT);
        check("Ratenwert 1e-13", fabs(el->rate.evaluate(3.75) - 1e-13) < 1e-20,
              to_string(el->rate.evaluate(3.75)));
    }

    // ── Test 3: Iod-Paket mit Molekuel, negativem Ion, Tabellen ──
    cout << "\n--- Test 3: iodine_lafleur_v1 laden ---" << endl;
    ChemLoadResult io = load_chem_package(resolve_chem_package("iodine_lafleur_v1"));
    check("Iod-Paket geladen", io.ok, io.error_text());
    if (io.ok) {
        const ChemSystem& is = io.system;
        check("5 schwere Spezies", is.species.size() == 5, to_string(is.species.size()));
        check("13 Reaktionen", is.reactions.size() == 13, to_string(is.reactions.size()));
        const ChemSpecies* i2 = find_species(is, "I2");
        const ChemSpecies* im = find_species(is, "I-");
        check("I2 ist Molekuel", i2 && i2->type == SpeciesType::NEUTRAL_MOLECULE);
        check("I2 ist Feedstock", i2 && i2->is_feedstock);
        check("I- ist negatives Ion", im && im->type == SpeciesType::NEGATIVE_ION);
        check("I- Ladung -1", im && im->charge == -1);
        check("Tabellen geladen", io.tables_loaded == 10, to_string(io.tables_loaded));

        const ChemReaction* iz_i2 = find_reaction(is, "iz_I2");
        check("iz_I2 tabelliert", iz_i2 && iz_i2->rate.type == RateType::TABULATED);
        if (iz_i2) {
            check("Tabelle nicht leer", !iz_i2->rate.table_Te.empty(),
                  to_string(iz_i2->rate.table_Te.size()));
            // Monotone Stuetzstellen und positive Raten im Betriebsbereich
            bool monoton = true;
            for (size_t k = 1; k < iz_i2->rate.table_Te.size(); ++k)
                if (iz_i2->rate.table_Te[k] <= iz_i2->rate.table_Te[k-1]) monoton = false;
            check("Stuetzstellen monoton", monoton);
            check("Rate bei 5 eV positiv", iz_i2->rate.evaluate(5.0) > 0,
                  to_string(iz_i2->rate.evaluate(5.0)));
            // Ausserhalb des Bereichs wird auf den Randwert festgehalten
            double unten = iz_i2->rate.evaluate(1e-3);
            check("unterhalb der Tabelle Randwert",
                  fabs(unten - iz_i2->rate.table_K.front()) < 1e-30 * (1 + fabs(unten)));
            double oben = iz_i2->rate.evaluate(1e6);
            check("oberhalb der Tabelle Randwert",
                  fabs(oben - iz_i2->rate.table_K.back()) < 1e-30 * (1 + fabs(oben)));
            // Stuetzstelle selbst wird exakt getroffen
            size_t mid = iz_i2->rate.table_Te.size() / 2;
            double an_stelle = iz_i2->rate.evaluate(iz_i2->rate.table_Te[mid]);
            check("Stuetzstelle exakt",
                  fabs(an_stelle - iz_i2->rate.table_K[mid]) <= 1e-12 * fabs(iz_i2->rate.table_K[mid]),
                  to_string(an_stelle));
        }
        const ChemReaction* surf = find_reaction(is, "surfrec_I");
        check("Oberflaechenreaktion vorhanden", surf != nullptr);
        check("Oberflaechenkoeffizient 0.02",
              surf && fabs(surf->surface_gamma - 0.02) < 1e-12);
    }

    // ── Test 4: fehlerhafte Pakete werden abgelehnt ──────────
    cout << "\n--- Test 4: Fehlerbehandlung ---" << endl;
    ChemLoadResult fehlt = load_chem_package("chemistry/gibt_es_nicht/chemistry.json");
    check("fehlende Datei gemeldet", !fehlt.ok && !fehlt.errors.empty());

    const char* kaputt = "test_chem_loader_tmp.json";
    {
        ofstream f(kaputt);
        f << R"({"name":"kaputt","species":[)"
             R"({"id":"e","name":"Elektron","type":"electron","mass_kg":9.10938215e-31,"charge":-1},)"
             R"({"id":"X","name":"X","type":"neutral_atom","mass_kg":1e-25,"is_feedstock":true}],)"
             R"("reactions":[{"id":"r1","type":"ionization","reactants":{"e":1,"X":1},)"
             R"("products":{"e":2,"X+":1},"rate":{"model":"arrhenius","A":1e-13,"E_a_eV":10.0}}]})";
    }
    ChemLoadResult kr = load_chem_package(kaputt);
    check("unvollstaendiges Paket abgelehnt", !kr.ok);
    bool nennt_ion = false, nennt_spezies = false;
    for (const auto& e : kr.errors) {
        if (e.find("positives Ion erforderlich") != string::npos) nennt_ion = true;
        if (e.find("'X+' nicht definiert") != string::npos) nennt_spezies = true;
    }
    check("fehlendes Ion benannt", nennt_ion, kr.error_text());
    check("undefinierte Spezies benannt", nennt_spezies, kr.error_text());
    remove(kaputt);

    const char* syntax = "test_chem_loader_tmp2.json";
    { ofstream f(syntax); f << "{ das ist kein JSON"; }
    ChemLoadResult sr = load_chem_package(syntax);
    check("Syntaxfehler gemeldet", !sr.ok && !sr.errors.empty());
    remove(syntax);

    // ── Test 5: geladenes Paket ist benutzbar ────────────────
    // Das Laden allein nuetzt nichts, wenn der Bilanzassembler mit dem
    // Ergebnis nicht rechnen kann. Deshalb hier eine Auswertung mit dem
    // geladenen Xenon-Paket.
    cout << "\n--- Test 5: geladenes Paket im Bilanzassembler ---" << endl;
    {
        SimContext ctx;
        ctx.recompute();
        const ChemSystem& s = xe.system;
        vector<double> state(s.state_size());
        state[s.species_index("Xe")]  = 1e19;
        state[s.species_index("Xe+")] = 1e17;
        state[s.Te_idx()] = 3.75;
        state[s.Tg_idx()] = 400.0;

        check("Elektronendichte gleich Ionendichte",
              fabs(electron_density(s, state) - 1e17) < 1e3,
              to_string(electron_density(s, state)));

        auto r = assemble_residual(s, state, 1e6, ctx.thruster.Q0, ctx.thruster, ctx);
        check("Residuum hat Zustandslaenge", (int)r.size() == s.state_size(),
              to_string(r.size()));
        bool endlich = true;
        for (double v : r) if (!isfinite(v)) endlich = false;
        check("Residuum endlich", endlich);
        bool wirkt = false;
        for (double v : r) if (v != 0.0) wirkt = true;
        check("Residuum nicht durchweg null", wirkt);
    }

    // ── Test 6: mehrfach geladene Ionen ──────────────────────
    // Der eigentliche Zweck der Uebung: ein Paket mit Xe2+ muss die Ladung
    // zwei tragen und in der Quasineutralitaet doppelt zaehlen.
    cout << "\n--- Test 6: mehrfach geladenes Ion ---" << endl;
    const char* zweifach = "test_chem_loader_tmp3.json";
    {
        ofstream f(zweifach);
        f << R"({"name":"Xe mit Xe2+","sigma_i":1e-18,"species":[)"
             R"({"id":"e","type":"electron","mass_kg":9.10938215e-31,"charge":-1},)"
             R"({"id":"Xe","type":"neutral_atom","mass_kg":2.1801711e-25,"charge":0,)"
             R"("is_feedstock":true,"thermal_conductivity":0.0057},)"
             R"({"id":"Xe+","type":"positive_ion","mass_kg":2.1801711e-25,"charge":1,)"
             R"("is_beam_extracted":true},)"
             R"({"id":"Xe2+","type":"positive_ion","mass_kg":2.1801711e-25,"charge":2,)"
             R"("is_beam_extracted":true}],)"
             R"("reactions":[)"
             R"({"id":"iz_Xe","type":"ionization","reactants":{"e":1,"Xe":1},)"
             R"("products":{"e":2,"Xe+":1},"energy_eV":12.127,)"
             R"("rate":{"model":"arrhenius","A":1.8e-13,"E_a_eV":12.127}},)"
             R"({"id":"iz_Xe2","type":"ionization","reactants":{"e":1,"Xe+":1},)"
             R"("products":{"e":2,"Xe2+":1},"energy_eV":21.21,)"
             R"("rate":{"model":"arrhenius","A":9.0e-14,"E_a_eV":21.21}}]})";
    }
    ChemLoadResult zw = load_chem_package(zweifach);
    check("Paket mit Xe2+ geladen", zw.ok, zw.error_text());
    if (zw.ok) {
        const ChemSystem& zs = zw.system;
        const ChemSpecies* xe2 = find_species(zs, "Xe2+");
        check("Xe2+ vorhanden", xe2 != nullptr);
        check("Ladung 2", xe2 && xe2->charge == 2, xe2 ? to_string(xe2->charge) : "-");
        check("3 schwere Spezies", zs.species.size() == 3, to_string(zs.species.size()));
        check("Zustandslaenge 5", zs.state_size() == 5, to_string(zs.state_size()));

        vector<double> st(zs.state_size());
        st[zs.species_index("Xe")]   = 1e19;
        st[zs.species_index("Xe+")]  = 1e17;
        st[zs.species_index("Xe2+")] = 2e16;
        st[zs.Te_idx()] = 4.0;
        st[zs.Tg_idx()] = 400.0;
        // n_e = 1*n(Xe+) + 2*n(Xe2+)
        double ne_soll = 1e17 + 2 * 2e16;
        check("Quasineutralitaet zaehlt die Ladung doppelt",
              fabs(electron_density(zs, st) - ne_soll) < 1e3,
              to_string(electron_density(zs, st)) + " statt " + to_string(ne_soll));

        // Zweite Ionisationsstufe verbraucht Xe+ und erzeugt Xe2+
        const ChemReaction* r2 = find_reaction(zs, "iz_Xe2");
        check("zweite Stufe vorhanden", r2 != nullptr);
        if (r2) {
            auto net = r2->net_stoichiometry();
            check("Xe+ wird verbraucht", net["Xe+"] == -1, to_string(net["Xe+"]));
            check("Xe2+ entsteht", net["Xe2+"] == 1, to_string(net["Xe2+"]));
        }
    }
    remove(zweifach);

    // ── Test 7: Wandprodukte ─────────────────────────────────
    // Was ein Ion an der Wand hinterlaesst, gehoert in das Paket und darf
    // nicht aus der Masse geraten werden. Fuer einatomige Ionen faellt beides
    // zusammen, fuer molekulare nicht.
    cout << "\n--- Test 7: Wandprodukte ---" << endl;
    {
        const ChemSpecies* xe_ion2 = find_species(xe.system, "Xe+");
        check("Xe+ hinterlaesst ein Xe",
              xe_ion2 && stoich_text(xe_ion2->wall_products) == "Xe:1",
              xe_ion2 ? stoich_text(xe_ion2->wall_products) : "-");
        check("nichts abgeleitet, weil erklaert", xe.wall_products_derived == 0,
              to_string(xe.wall_products_derived));

        // Molekuelion, das an der Wand in zwei Atome zerfaellt
        const char* mol = "test_chem_loader_tmp4.json";
        {
            ofstream f(mol);
            f << R"({"name":"Iod mit dissoziativer Neutralisation","sigma_i":1e-18,"species":[)"
                 R"({"id":"e","type":"electron","mass_kg":9.10938215e-31,"charge":-1},)"
                 R"({"id":"I2","type":"neutral_molecule","mass_kg":4.2143422e-25,"charge":0,)"
                 R"("is_feedstock":true,"thermal_conductivity":0.0039},)"
                 R"({"id":"I","type":"neutral_atom","mass_kg":2.1071711e-25,"charge":0,)"
                 R"("thermal_conductivity":0.0039},)"
                 R"({"id":"I2+","type":"positive_ion","mass_kg":4.2143422e-25,"charge":1,)"
                 R"("is_beam_extracted":true,"wall_products":{"I":2}}],)"
                 R"("reactions":[)"
                 R"({"id":"iz_I2","type":"ionization","reactants":{"e":1,"I2":1},)"
                 R"("products":{"e":2,"I2+":1},"energy_eV":9.31,)"
                 R"("rate":{"model":"arrhenius","A":1e-13,"E_a_eV":9.31}}]})";
        }
        ChemLoadResult mr = load_chem_package(mol);
        check("Paket mit zwei Wandprodukten geladen", mr.ok, mr.error_text());
        if (mr.ok) {
            const ChemSpecies* i2p = find_species(mr.system, "I2+");
            check("I2+ hinterlaesst zwei I",
                  i2p && stoich_text(i2p->wall_products) == "I:2",
                  i2p ? stoich_text(i2p->wall_products) : "-");
            check("nichts abgeleitet", mr.wall_products_derived == 0);

            // Der Rueckstrom muss doppelt in die Atombilanz gehen
            SimContext ctx2;
            ctx2.recompute();
            const ChemSystem& ms = mr.system;
            vector<double> st(ms.state_size());
            st[ms.species_index("I2")]  = 1e19;
            st[ms.species_index("I")]   = 1e18;
            st[ms.species_index("I2+")] = 1e17;
            st[ms.Te_idx()] = 4.0;
            st[ms.Tg_idx()] = 400.0;
            auto r = assemble_residual(ms, st, 1e6, ctx2.thruster.Q0, ctx2.thruster, ctx2);
            // Ohne Reaktionen auf I bleibt nur Rueckstrom minus Ausstroemen
            double uB = sqrt(1 * PhysConst::kB * 4.0 * PhysConst::conv / 4.2143422e-25);
            double lam = 1.0 / (1.1e19 * 1e-18);
            double hL = 0.86 * pow(3.0 + ctx2.thruster.L/(2*lam), -0.5);
            double hR = 0.80 * pow(4.0 + ctx2.thruster.R/lam, -0.5);
            double Aeff1 = 2*hR*PhysConst::pi*ctx2.thruster.R*ctx2.thruster.L
                         + (2 - ctx2.thruster.betai)*hL*PhysConst::pi*ctx2.thruster.R*ctx2.thruster.R;
            double rueck = 2.0 * 1e17 * uB * Aeff1 / ctx2.thruster.V;
            double v_mean = sqrt(8*PhysConst::kB*400.0/(PhysConst::pi*2.1071711e-25));
            double raus = 0.25 * 1e18 * v_mean * ctx2.thruster.Ag / ctx2.thruster.V;
            double erwartet = rueck - raus;
            double ist = r[ms.species_index("I")];
            check("Rueckstrom geht doppelt in die Atombilanz",
                  fabs(ist - erwartet) <= 1e-9 * fabs(erwartet),
                  to_string(ist) + " statt " + to_string(erwartet));
        }
        remove(mol);

        // Fehlerfaelle
        const char* schlecht = "test_chem_loader_tmp5.json";
        auto schreibe = [&](const char* wp) {
            ofstream f(schlecht);
            f << R"({"name":"kaputt","sigma_i":1e-18,"species":[)"
                 R"({"id":"e","type":"electron","mass_kg":9.10938215e-31,"charge":-1},)"
                 R"({"id":"Xe","type":"neutral_atom","mass_kg":2.1801711e-25,"charge":0,)"
                 R"("is_feedstock":true},)"
                 R"({"id":"Xe+","type":"positive_ion","mass_kg":2.1801711e-25,"charge":1,)"
                 R"("wall_products":)" << wp << R"(}],)"
                 R"("reactions":[{"id":"iz","type":"ionization","reactants":{"e":1,"Xe":1},)"
                 R"("products":{"e":2,"Xe+":1},"energy_eV":12.1,)"
                 R"("rate":{"model":"constant","value":1e-15}}]})";
        };

        schreibe(R"({"Kr":1})");
        ChemLoadResult e1 = load_chem_package(schlecht);
        check("unbekanntes Wandprodukt abgelehnt", !e1.ok,
              e1.ok ? "durchgelassen" : "");
        check("unbekanntes Wandprodukt benannt",
              e1.error_text().find("nicht definiert") != string::npos, e1.error_text());

        schreibe(R"({"Xe+":1})");
        ChemLoadResult e2 = load_chem_package(schlecht);
        check("nichtneutrales Wandprodukt abgelehnt", !e2.ok);
        check("nichtneutrales Wandprodukt benannt",
              e2.error_text().find("nicht neutral") != string::npos, e2.error_text());

        schreibe(R"({"Xe":2})");
        ChemLoadResult e3 = load_chem_package(schlecht);
        check("Massenbilanz an der Wand geprueft", !e3.ok);
        check("Massenfehler benannt",
              e3.error_text().find("wiegen nicht") != string::npos, e3.error_text());
        remove(schlecht);

        // Fehlt die Angabe, wird sie abgeleitet und gezaehlt
        const char* ohne = "test_chem_loader_tmp6.json";
        {
            ofstream f(ohne);
            f << R"({"name":"ohne Angabe","sigma_i":1e-18,"species":[)"
                 R"({"id":"e","type":"electron","mass_kg":9.10938215e-31,"charge":-1},)"
                 R"({"id":"Xe","type":"neutral_atom","mass_kg":2.1801711e-25,"charge":0,)"
                 R"("is_feedstock":true},)"
                 R"({"id":"Xe+","type":"positive_ion","mass_kg":2.1801711e-25,"charge":1}],)"
                 R"("reactions":[{"id":"iz","type":"ionization","reactants":{"e":1,"Xe":1},)"
                 R"("products":{"e":2,"Xe+":1},"energy_eV":12.1,)"
                 R"("rate":{"model":"constant","value":1e-15}}]})";
        }
        ChemLoadResult ar = load_chem_package(ohne);
        check("ohne Angabe abgeleitet", ar.ok, ar.error_text());
        check("Ableitung wird gezaehlt", ar.wall_products_derived == 1,
              to_string(ar.wall_products_derived));
        const ChemSpecies* abg = ar.ok ? find_species(ar.system, "Xe+") : nullptr;
        check("Ableitung trifft das Neutralteilchen",
              abg && stoich_text(abg->wall_products) == "Xe:1",
              abg ? stoich_text(abg->wall_products) : "-");
        remove(ohne);
    }

    // ── Test 8: Kurzname zu Pfad ─────────────────────────────
    cout << "\n--- Test 8: Pfadaufloesung ---" << endl;
    check("Kurzname aufgeloest",
          resolve_chem_package("xenon_simple") == "chemistry/xenon_simple/chemistry.json",
          resolve_chem_package("xenon_simple"));
    check("vollstaendiger Pfad bleibt",
          resolve_chem_package("a/b/c.json") == "a/b/c.json");

    cout << "\n========================================" << endl;
    cout << "  Bestanden: " << passed << "   Fehlgeschlagen: " << failed << endl;
    cout << "========================================" << endl;
    return failed == 0 ? 0 : 1;
}
