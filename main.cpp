// main.cpp -- Einstiegspunkt. Kein globaler Zustand, alles ueber SimContext.
//
// Build (6 Kompilationseinheiten):
//   g++ -O3 -std=c++17 -c bessel_wrapper.cpp -o bessel_wrapper.o
//   g++ -O3 -std=c++17 -c sim_config.cpp -o sim_config.o
//   g++ -O3 -std=c++17 -c rates.cpp -o rates.o
//   g++ -O3 -std=c++17 -c physics.cpp -o physics.o
//   g++ -O3 -std=c++17 -c solver.cpp -o solver.o
//   g++ -O3 -std=c++17 -c sim_logging.cpp -o sim_logging.o
//   g++ -O3 -std=c++17 main.cpp *.o -o chabert
#include "solver.hpp"
#include "sim_logging.hpp"
#include "sim_config.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <vector>

using namespace std;

int main(int argc, char** argv) {
    string configFile = "params.txt";
    if (argc >= 2) configFile = argv[1];

    // Einziger SimContext -- kein globaler Zustand
    SimContext ctx;
    auto cd = loadConfig(configFile);
    applyConfig(ctx, cd);

    auto start_time = chrono::high_resolution_clock::now();
    int count_ok = 0, count_nophys = 0, count_numfail = 0;

    vector<SimLogRow> log_rows;
    vector<SimLogEvent> log_events;

    string red="\033[31m", green="\033[32m", reset="\033[0m";

    try {
        auto& t = ctx.thruster;
        auto& sp = ctx.solver;
        auto& g = ctx.gas;
        auto& r = ctx.rates;

        cout << green
             << "####################################################\n"
             << "#          Global Plasma Model - " << g.species << "\n"
             << "#   Solver: Newton (stationaer)                    #\n"
             << "####################################################\n"
             << reset << endl;

        cout << "GAS_SPECIES " << g.species << endl;
        cout << "GAS_MASS " << scientific << setprecision(6) << g.M << " kg" << endl;
        cout << "SOLVE_MODE " << sp.solve_mode << " "
             << (sp.solve_mode == 2 ? "selbstkonsistent" : "fester_strahlstrom") << endl;

        const char* rmn = (ctx.rate_model==2) ? "Full tabulated (Biagi/LXCat)"
                        : (ctx.rate_model==1) ? "Conservative tabulated (Kiz+Kex tab, Kel legacy)"
                        :                       "Legacy (paper-compatible)";
        cout << "RATE_MODEL " << ctx.rate_model << " " << rmn << endl;

        string cs_base = "cross_sections/" + g.species + "/" + ctx.cs_database + "/";
        cout << "CS_DATABASE " << ctx.cs_database << " (" << cs_base << ")" << endl;

        // Ratentabellen laden
        if (r.elastic_model == 1) {
            if (load_kel_table(r, cs_base + "kel_table.csv"))
                cout << "ELASTIC_MODEL tabulated (" << r.kel.size() << " Te-Punkte)" << endl;
            else { cerr << "WARNUNG: kel_table.csv nicht gefunden, Legacy-Modus" << endl; r.elastic_model = 0; }
        }
        if (r.elastic_model == 0) cout << "ELASTIC_MODEL legacy (Kel=" << r.kel_constant << ")" << endl;

        if (r.ionization_model == 1) {
            if (load_kiz_table(r, cs_base + "kiz_table.csv"))
                cout << "IONIZATION_MODEL tabulated (" << r.kiz.size() << " Te-Punkte)" << endl;
            else { cerr << "WARNUNG: kiz_table.csv nicht gefunden, Legacy-Modus" << endl; r.ionization_model = 0; }
        }
        if (r.ionization_model == 0) cout << "IONIZATION_MODEL legacy (Chabert Polynomfit)" << endl;

        if (r.excitation_model == 1) {
            if (load_kex_table(r, cs_base + "kex_table.csv"))
                cout << "EXCITATION_MODEL tabulated (" << r.kex.size() << " Te-Punkte)" << endl;
            else { cerr << "WARNUNG: kex_table.csv nicht gefunden, Legacy-Modus" << endl; r.excitation_model = 0; }
        }
        if (r.excitation_model == 0) cout << "EXCITATION_MODEL legacy (Chabert Arrhenius)" << endl;

        ofstream datei("output_kh.txt");
        if (!datei) { cerr << "Fehler: Ausgabedatei!" << endl; return 1; }
        datei << "Method, Q0sccm, Te, Tg, n, ng, iondeg, P_RFG, P_abs, I_extr_mA, "
              << "collision_freq, R_induktiv, I_coil, epsilon_p_real, epsilon_p_imag, "
              << "u_Bohm, J_i, zeta, gamma, xi, eta, plasmafrequenz, frequency_MHz, "
              << "thrust_ions_mN, thrust_atoms_mN, thrust_total_mN, icp_power_efficiency, P_RF_W, "
              << "n_eff, density_profile_factor\n";

        if (t.Rc < t.R) {
            cout << red << "FEHLER: Rc=" << t.Rc << " < R=" << t.R << " - Spule innerhalb Kammer!" << reset << endl;
            return 1;
        }
        cout << "CL_LIMIT " << scientific << setprecision(6) << t.J_CL << " " << t.J_CL*t.Ai*1000 << endl;

        double P_RFG_start = t.P_RFG;
        PlasmaState prev{}, last_good{};
        bool have_prev = false, have_good = false;
        double last_good_q0 = 0;

        // ═══ Q0-Sweep ══════════════════════════════════════
        for (int jj = 0; jj < sp.jjmax; ++jj) {
            t.Q0sccm = sp.Q0sccm_start + jj * sp.Q0sccm_step;
            t.Q0 = t.Q0sccm * PhysConst::SCCM_TO_PPS;
            cout << "Q0_STEP " << fixed << setprecision(4) << t.Q0sccm << " " << (jj+1) << " " << sp.jjmax << endl;

            PlasmaState guess;
            if (have_good && fabs(t.Q0sccm - last_good_q0) <= 20*sp.Q0sccm_step) guess = last_good;
            else if (have_prev) guess = prev;
            else guess = safe_defaults(ctx);
            if (!state_finite_positive(guess)) guess = safe_defaults(ctx);

            PowerSolveResult ps;
            if (sp.solve_mode == 2)
                ps = solve_at_fixed_power(ctx, t.P_RFG, guess);
            else
                ps = solve_for_target_current(ctx, guess);

            // IFIX_RESULT: strukturierte Ausgabe pro I-fix-Punkt (analog zu Python)
            if (sp.solve_mode == 1) {
                double I_found = isfinite(ps.I_mA) ? ps.I_mA : 0;
                double P_found = isfinite(ps.P_RFG_sol) ? ps.P_RFG_sol : 0;
                double delta_I = I_found - sp.I_soll;
                string status_str = ps.converged ? "converged" :
                    (ps.fail_type == SolveFailType::NO_PHYSICAL_SOLUTION ? "above_P_max" : "numerical_fail");
                cout << "IFIX_RESULT " << fixed << setprecision(4) << t.Q0sccm << " "
                     << setprecision(2) << P_found << " " << sp.I_soll << " "
                     << I_found << " " << delta_I << " " << status_str << endl;
            }

            if (!power_result_valid(ctx, ps)) {
                string ft = (ps.fail_type == SolveFailType::NO_PHYSICAL_SOLUTION)
                             ? "NO_PHYSICAL_SOLUTION" : "NUMERICAL_FAIL";
                if (ps.fail_type == SolveFailType::NO_PHYSICAL_SOLUTION) {
                    count_nophys++;
                    cout << "NO_PHYSICAL_SOLUTION " << jj << " " << fixed << setprecision(4)
                         << t.Q0sccm << " " << ps.reason
                         << " I_best=" << (isfinite(ps.I_mA)?ps.I_mA:-1)
                         << " P_max_tried=" << (isfinite(ps.P_trial_last)?ps.P_trial_last:-1) << endl;
                } else {
                    count_numfail++;
                    cout << "NUMERICAL_FAIL " << jj << " " << fixed << setprecision(4)
                         << t.Q0sccm << " " << ps.reason << endl;
                }
                log_rows.push_back({jj, t.Q0sccm, ft, ft,
                    isfinite(ps.P_trial_last)?ps.P_trial_last:0, isfinite(ps.I_mA)?ps.I_mA:0,
                    0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,0, 0,0,0, 0,0,0,0, ps.reason, false});
                log_events.push_back({jj, t.Q0sccm, ft + ": " + ps.reason});
                continue;
            }
            count_ok++;
            prev = ps.state; have_prev = true;
            if (isfinite(ps.inner_resid_norm) && ps.inner_resid_norm < sp.newton_tol) {
                last_good = ps.state; have_good = true; last_good_q0 = t.Q0sccm;
            }
            t.P_RFG = ps.P_RFG_sol;

            double n = ps.state.n, ng = ps.state.ng, Te = ps.state.Te, Tg = ps.state.Tg;
            DerivedQuantities dq = compute_derived(ctx, n, ng, Te, Tg, ps.rf.R_ind, ps.rf.I_coil, ps.rf.P_abs);

            if (isfinite(ps.inner_resid_norm) && ps.inner_resid_norm >= sp.newton_tol)
                cout << "SOFT_ACCEPT " << fixed << setprecision(4) << t.Q0sccm << " " << scientific << ps.inner_resid_norm << endl;

            cout << "RESULT " << scientific << setprecision(4) << n << " " << ng << " "
                 << fixed << setprecision(3) << Te << " " << Tg << " " << dq.I_extr_mA << " " << t.P_RFG << endl;

            cout << "RESULT_EXT " << fixed << setprecision(6) << t.Q0sccm << " " << t.P_RFG << " "
                 << scientific << setprecision(4) << n << " " << ng << " " << n << " "
                 << fixed << setprecision(4)
                 << dq.T_i_N*1e3 << " " << dq.T_n_N*1e3 << " " << dq.T_total_N*1e3 << " "
                 << dq.icp_eff << " " << dq.gamma_eff << " " << dq.xi_mN_kW << " " << dq.eta_mass << endl;

            log_rows.push_back({jj, t.Q0sccm, "CONVERGED", "NONE",
                t.P_RFG, dq.I_extr_mA, Te, Tg, n, ng, ps.inner_resid_norm,
                dq.iondeg, ps.rf.P_abs, dq.cf, ps.rf.R_ind, ps.rf.I_coil,
                dq.eps_p_real, dq.eps_p_imag, dq.u_Bohm, dq.J_i, dq.pf, dq.P_RF,
                dq.T_i_N*1e3, dq.T_n_N*1e3, dq.T_total_N*1e3,
                dq.icp_eff, dq.gamma_eff, dq.xi_mN_kW, dq.eta_mass,
                "ok", true});

            emit_csv_row(datei, "Stationary", t.Q0sccm, n, ng, Te, Tg,
                         t.P_RFG, ps.rf.P_abs, ps.rf.R_ind, ps.rf.I_coil, ctx, dq);

            t.P_RFG = P_RFG_start;
        }
    } catch (const exception& ex) { cerr << "Fehler: " << ex.what() << endl; }

    auto end_time = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration<double>(end_time - start_time).count();

    cout << "\nSUMMARY " << count_ok << " " << count_nophys << " " << count_numfail << endl;
    cout << "SUMMARY_DETAIL converged=" << count_ok
         << " no_physical_solution=" << count_nophys
         << " numerical_fail=" << count_numfail
         << " total=" << ctx.solver.jjmax << endl;
    cout << "\nDer Code hat " << fixed << setprecision(3) << elapsed << " Sekunden gebraucht." << endl;

    write_masterlog(ctx, configFile, elapsed, count_ok, count_nophys, count_numfail, log_rows, log_events);

    return 0;
}
