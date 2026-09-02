// test_chem_system.cpp -- Test fuer den generischen Chemie-Kern.
//
// Vergleicht den generischen Residualassembler (assemble_residual) mit
// dem bisherigen Xenon-Hardcode (residual_raw) fuer identische Zustaende.
// Prueft Spezies-/Reaktionsstruktur und numerische Konsistenz.
//
// Build:
//   g++ -O3 -std=c++17 test_chem_system.cpp sim_config.o rates.o physics.o \
//       solver.o sim_logging.o bessel_wrapper.o -o test_chem_system
// Run:
//   ./test_chem_system
#include "chem_system.hpp"
#include "physics.hpp"
#include "sim_config.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

using namespace std;

int passed = 0, failed = 0;

void check(const string& name, bool cond, const string& detail = "") {
    if (cond) { ++passed; cout << "  PASS: " << name << endl; }
    else { ++failed; cout << "  FAIL: " << name << " -- " << detail << endl; }
}

int main() {
    // ── 1. SimContext mit Chabert-Defaults ────────────────────
    SimContext ctx;
    ctx.recompute();

    // ── Test 1: ChemSystem Aufbau ────────────────────────────
    cout << "\n--- Test 1: ChemSystem Aufbau ---" << endl;
    ChemSystem sys = build_xenon_system(ctx);

    check("system name", sys.name == "Xenon (generic)");
    check("2 species", sys.species.size() == 2, to_string(sys.species.size()));
    check("3 reactions", sys.reactions.size() == 3, to_string(sys.reactions.size()));
    check("Xe is neutral", sys.species[0].is_neutral());
    check("Xe is feedstock", sys.species[0].is_feedstock);
    check("Xe+ is positive ion", sys.species[1].is_positive_ion());
    check("state_size = 4", sys.state_size() == 4);  // Xe, Xe+, Te, Tg
    check("Te_idx = 2", sys.Te_idx() == 2);
    check("Tg_idx = 3", sys.Tg_idx() == 3);
    check("species_index Xe = 0", sys.species_index("Xe") == 0);
    check("species_index Xe+ = 1", sys.species_index("Xe+") == 1);

    // ── Test 2: Ratenkoeffizienten ───────────────────────────
    cout << "\n--- Test 2: Ratenkoeffizienten ---" << endl;
    double Te = 3.75;
    double K_iz = sys.reactions[0].rate.evaluate(ctx, Te);
    double K_ex = sys.reactions[1].rate.evaluate(ctx, Te);
    double K_el = sys.reactions[2].rate.evaluate(ctx, Te);

    check("Kiz > 0", K_iz > 0, to_string(K_iz));
    check("Kex > 0", K_ex > 0, to_string(K_ex));
    check("Kel > 0", K_el > 0, to_string(K_el));

    // Vergleich mit Legacy-Funktionen
    double K_iz_legacy = Kiz(ctx, Te);
    double K_ex_legacy = Kex(ctx, Te);
    double K_el_legacy = Kel(ctx, Te);

    check("Kiz matches legacy", abs(K_iz - K_iz_legacy) / K_iz_legacy < 1e-10,
          to_string(K_iz) + " vs " + to_string(K_iz_legacy));
    check("Kex matches legacy", abs(K_ex - K_ex_legacy) / K_ex_legacy < 1e-10);
    check("Kel matches legacy", abs(K_el - K_el_legacy) / K_el_legacy < 1e-10);

    // ── Test 3: Generischer Residual vs Hardcode ─────────────
    cout << "\n--- Test 3: Residual-Vergleich ---" << endl;

    // Typischer konvergierter Zustand
    double n = 2e17, ng = 3e19;
    double Tg_val = 300.0;
    PlasmaState ps{n, ng, Te, Tg_val};

    // Generischer Residual
    vector<double> state = {ng, n, Te, Tg_val};  // Xe, Xe+, Te, Tg
    double P_RFG = 18.0;
    RFState rf = compute_rf(ctx, n, ng, Te, P_RFG);

    double P_abs_V = 0;
    if (rf.valid) P_abs_V = rf.P_abs * ctx.solver.P_abs_scale / ctx.thruster.V;

    auto gen_r = assemble_residual(sys, state, P_abs_V, ctx.thruster.Q0,
                                    ctx.thruster, ctx,
                                    ctx.solver.alpha_e_wall,
                                    ctx.solver.density_profile_factor);

    // Legacy Residual
    auto leg_r = residual_raw(ctx, ps, P_RFG, nullptr);

    cout << "  Generic residual: ["
         << gen_r[0] << ", " << gen_r[1] << ", " << gen_r[2] << ", " << gen_r[3] << "]" << endl;
    cout << "  Legacy residual:  ["
         << leg_r[0] << ", " << leg_r[1] << ", " << leg_r[2] << ", " << leg_r[3] << "]" << endl;

    // Vergleich (Mappung: generic[0]=Xe=r2(neutral), generic[1]=Xe+=r1(ion))
    // Legacy: r1=Ionenbilanz, r2=Neutralbilanz, r3=Eenergie, r4=Gasenergie
    // Generic: [Xe-Bilanz, Xe+-Bilanz, Te-Bilanz, Tg-Bilanz]
    //   Xe-Bilanz  = Legacy r2 (Neutralgasbilanz)
    //   Xe+-Bilanz = Legacy r1 (Ionenbilanz)
    //   Te-Bilanz  = Legacy r3 (Elektronenenergie)
    //   Tg-Bilanz  = Legacy r4 (Gasenergie)

    double tol = 0.05;  // 5% relative Toleranz (unterschiedliche Detailformulierung)

    // Ionenbilanz: generic Xe+ vs legacy r1
    double diff_ion = abs(gen_r[1] - leg_r[0]) / (abs(leg_r[0]) + 1e-30);
    check("Ion balance matches (<5%)", diff_ion < tol,
          "gen=" + to_string(gen_r[1]) + " leg=" + to_string(leg_r[0]) + " diff=" + to_string(diff_ion*100) + "%");

    // Neutralbilanz: generic Xe vs legacy r2
    double diff_neut = abs(gen_r[0] - leg_r[1]) / (abs(leg_r[1]) + 1e-30);
    check("Neutral balance matches (<5%)", diff_neut < tol,
          "gen=" + to_string(gen_r[0]) + " leg=" + to_string(leg_r[1]) + " diff=" + to_string(diff_neut*100) + "%");

    // Elektronenenergie: generic Te vs legacy r3
    double diff_ee = abs(gen_r[2] - leg_r[2]) / (abs(leg_r[2]) + 1e-30);
    check("Electron energy matches (<5%)", diff_ee < tol,
          "gen=" + to_string(gen_r[2]) + " leg=" + to_string(leg_r[2]) + " diff=" + to_string(diff_ee*100) + "%");

    // Gasenergie: generic Tg vs legacy r4
    double diff_ge = abs(gen_r[3] - leg_r[3]) / (abs(leg_r[3]) + 1e-30);
    check("Gas energy matches (<5%)", diff_ge < tol,
          "gen=" + to_string(gen_r[3]) + " leg=" + to_string(leg_r[3]) + " diff=" + to_string(diff_ge*100) + "%");

    // ── Test 4: Quasineutralitaet ────────────────────────────
    cout << "\n--- Test 4: Quasineutralitaet ---" << endl;
    double ne = electron_density(sys, state);
    check("ne = n(Xe+)", abs(ne - n) / n < 1e-10, to_string(ne) + " vs " + to_string(n));

    // Mit negativen Ionen (hypothetisch)
    ChemSystem sys2 = sys;
    sys2.species.push_back({"X-", SpeciesType::NEGATIVE_ION, ctx.gas.M, -1});
    vector<double> state2 = {ng, n, 1e15, Te, Tg_val};  // Xe, Xe+, X-, Te, Tg
    double ne2 = electron_density(sys2, state2);
    check("ne with negative ions", abs(ne2 - (n - 1e15)) < 1e5,
          to_string(ne2) + " vs " + to_string(n - 1e15));

    // ── Zusammenfassung ──────────────────────────────────────
    cout << "\n" << string(60, '=') << endl;
    cout << "  Ergebnis: " << passed << " passed, " << failed << " failed" << endl;
    cout << string(60, '=') << endl;

    return failed > 0 ? 1 : 0;
}
