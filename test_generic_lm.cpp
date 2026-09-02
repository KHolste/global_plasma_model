// test_generic_lm.cpp -- Test: Generischer LM-Solver fuer Xenon.
//
// Prueft:
//   1. Generischer Solver konvergiert fuer Xenon
//   2. Ergebnis ist konsistent mit Legacy-Solver
//   3. StateLayout korrekt
//   4. PlasmaState-Extraktion funktioniert
//
// Build:
//   g++ -O3 -std=c++17 test_generic_lm.cpp sim_config.o rates.o physics.o \
//       solver.o sim_logging.o bessel_wrapper.o -o test_generic_lm
#include "generic_lm.hpp"
#include "solver.hpp"
#include "sim_config.hpp"
#include <iostream>
#include <cmath>
#include <string>

using namespace std;

int passed = 0, failed = 0;

void check(const string& name, bool cond, const string& detail = "") {
    if (cond) { ++passed; cout << "  PASS: " << name << endl; }
    else { ++failed; cout << "  FAIL: " << name << " -- " << detail << endl; }
}

int main() {
    SimContext ctx;
    ctx.recompute();

    // ── Test 1: StateLayout ──────────────────────────────────
    cout << "\n--- Test 1: StateLayout ---" << endl;
    ChemSystem sys = build_xenon_system(ctx);
    StateLayout lay(sys);

    check("N=4", lay.N == 4);
    check("n_species=2", lay.n_species == 2);
    check("Te_idx=2", lay.Te_idx == 2);
    check("Tg_idx=3", lay.Tg_idx == 3);

    // ── Test 2: Default-Startzustand ─────────────────────────
    cout << "\n--- Test 2: Default-Startzustand ---" << endl;
    auto x0 = gen_default_state(sys, ctx.thruster.Q0, ctx.thruster.V);
    check("state size=4", (int)x0.size() == 4);
    check("Xe > 0", x0[0] > 0, to_string(x0[0]));
    check("Xe+ > 0", x0[1] > 0, to_string(x0[1]));
    check("Te=3.75", abs(x0[2] - 3.75) < 0.01);
    check("Tg=300", abs(x0[3] - 300) < 1);

    // ── Test 3: Generischer LM-Solver konvergiert ────────────
    cout << "\n--- Test 3: Generischer LM-Solver ---" << endl;
    double P_RFG = 18.0;

    // Multi-Start mit verschiedenen Anfangswerten
    auto x1 = x0;
    x1[0] *= 0.5; x1[1] *= 2;  // Andere Dichteverhaeltnisse
    auto x2 = x0;
    x2[0] *= 2; x2[1] *= 0.1; x2[2] = 5.0;  // Hoeheres Te
    auto x3 = x0;
    x3[0] *= 0.1; x3[1] *= 10; x3[2] = 2.5;  // Niedrigeres Te

    GenSolveResult gr = gen_solve_multistart(sys, ctx, P_RFG, {x0, x1, x2, x3});
    check("converged", gr.converged, "reason=" + gr.reason);
    check("resid < 0.1", gr.resid_norm < 0.1, to_string(gr.resid_norm));
    check("iterations > 0", gr.iterations > 0, to_string(gr.iterations));
    check("RF valid", gr.rf.valid);

    if (gr.converged) {
        cout << "  State: Xe=" << gr.state[0] << " Xe+=" << gr.state[1]
             << " Te=" << gr.state[2] << " Tg=" << gr.state[3] << endl;

        check("Xe > 1e16", gr.state[0] > 1e16);
        check("Xe+ > 1e15", gr.state[1] > 1e15);
        check("Te in [2,8]", gr.state[2] > 2 && gr.state[2] < 8,
              to_string(gr.state[2]));
        check("Tg in [250,400]", gr.state[3] > 250 && gr.state[3] < 400,
              to_string(gr.state[3]));
    }

    // ── Test 4: Vergleich mit Legacy-Solver ──────────────────
    cout << "\n--- Test 4: Legacy-Vergleich ---" << endl;
    PlasmaState guess = safe_defaults(ctx);
    SolveResult lr = solve_lm(ctx, P_RFG, guess);
    check("legacy converged", lr.converged, "reason=" + lr.reason);

    if (gr.converged && lr.converged) {
        // Generisch: state[1]=Xe+ (n), state[0]=Xe (ng)
        // Legacy: state.n (Xe+), state.ng (Xe)
        PlasmaState gps = gr.to_plasma_state(sys);

        double rel_n = abs(gps.n - lr.state.n) / lr.state.n;
        double rel_ng = abs(gps.ng - lr.state.ng) / lr.state.ng;
        double rel_Te = abs(gps.Te - lr.state.Te) / lr.state.Te;
        double rel_Tg = abs(gps.Tg - lr.state.Tg) / lr.state.Tg;

        cout << "  Generic: n=" << gps.n << " ng=" << gps.ng
             << " Te=" << gps.Te << " Tg=" << gps.Tg << endl;
        cout << "  Legacy:  n=" << lr.state.n << " ng=" << lr.state.ng
             << " Te=" << lr.state.Te << " Tg=" << lr.state.Tg << endl;
        cout << "  Diff:    n=" << rel_n*100 << "% ng=" << rel_ng*100
             << "% Te=" << rel_Te*100 << "% Tg=" << rel_Tg*100 << "%" << endl;

        // Toleranz: <10% (verschiedene Skalierungen, verschiedene Startpunkte)
        check("n within 10%", rel_n < 0.10, to_string(rel_n*100) + "%");
        check("ng within 15%", rel_ng < 0.15, to_string(rel_ng*100) + "%");
        check("Te within 10%", rel_Te < 0.10, to_string(rel_Te*100) + "%");
        check("Tg within 5%", rel_Tg < 0.05, to_string(rel_Tg*100) + "%");

        // Beam-Strom
        double I_gen = gen_beam_current_mA(sys, ctx, gr.state);
        double I_leg = beam_current_mA(ctx, lr.state);
        double rel_I = abs(I_gen - I_leg) / std::max(I_leg, 0.01);
        cout << "  I_beam: gen=" << I_gen << " leg=" << I_leg << " diff=" << rel_I*100 << "%" << endl;
        check("I_beam within 15%", rel_I < 0.15, to_string(rel_I*100) + "%");
    }

    // ── Test 5: Mehrere Betriebspunkte ───────────────────────
    cout << "\n--- Test 5: Mehrere Q0-Punkte ---" << endl;
    int conv_count = 0;
    for (double Q0sccm : {0.3, 0.4, 0.5, 0.6, 0.7}) {
        ctx.thruster.Q0sccm = Q0sccm;
        ctx.thruster.Q0 = Q0sccm * PhysConst::SCCM_TO_PPS;
        auto s0 = gen_default_state(sys, ctx.thruster.Q0, ctx.thruster.V);
        auto s1 = s0; s1[0] *= 0.5; s1[1] *= 2;
        auto s2 = s0; s2[0] *= 2; s2[1] *= 0.1; s2[2] = 5.0;
        GenSolveResult r = gen_solve_multistart(sys, ctx, 18.0, {s0, s1, s2});
        if (r.converged) conv_count++;
        cout << "  Q0=" << Q0sccm << " conv=" << r.converged
             << " iter=" << r.iterations << " resid=" << r.resid_norm << endl;
    }
    check(">=3 of 5 converged", conv_count >= 3, to_string(conv_count) + "/5");

    // ── Zusammenfassung ──────────────────────────────────────
    cout << "\n" << string(60, '=') << endl;
    cout << "  Ergebnis: " << passed << " passed, " << failed << " failed" << endl;
    cout << string(60, '=') << endl;
    return failed > 0 ? 1 : 0;
}
