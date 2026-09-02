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
#include "chem_loader.hpp"
#include "generic_lm.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <vector>

using namespace std;

// ═════════════════════════════════════════════════════════════
// Generischer Rechenweg ueber ein geladenes Chemiepaket
// ═════════════════════════════════════════════════════════════

// Zusammengefasster Zustand fuer die abgeleiteten Groessen. Diese rechnen
// weiterhin mit einer Ionensorte und einer Masse; die Verallgemeinerung auf
// Ladungszahl und mehrere Sorten steht noch aus.
static PlasmaState lumped_state(const ChemSystem& sys, const vector<double>& st) {
    PlasmaState ps{};
    for (int i = 0; i < (int)sys.species.size(); ++i) {
        if (sys.species[i].is_positive_ion()) ps.n += st[i];
        else if (sys.species[i].is_neutral()) ps.ng += st[i];
    }
    ps.Te = st[sys.Te_idx()];
    ps.Tg = st[sys.Tg_idx()];
    return ps;
}

// Stoffwerte aus dem Paket in den Kontext uebernehmen, soweit die abgeleiteten
// Groessen sie brauchen: Masse aus dem extrahierten Ion, Waermeleitfaehigkeit
// aus der zugefuehrten Sorte, Stossquerschnitt aus dem Paketkopf.
static void adopt_package_properties(SimContext& ctx, const ChemSystem& sys) {
    const ChemSpecies* ion = nullptr;
    for (const auto& sp : sys.species) {
        if (!sp.is_positive_ion()) continue;
        if (!ion) ion = &sp;
        if (sp.is_beam_extracted) { ion = &sp; break; }
    }
    const ChemSpecies* feed = sys.feedstock();
    if (ion) ctx.gas.M = ion->mass_kg;
    if (feed && feed->thermal_cond > 0) ctx.gas.kappa = feed->thermal_cond;
    ctx.gas.sigma_i = sys.sigma_i;
    ctx.recompute();
}

static void run_generic_sweep(SimContext& ctx, const ChemSystem& sys, ofstream& datei,
                              vector<SimLogRow>& log_rows, vector<SimLogEvent>& log_events,
                              int& count_ok, int& count_nophys, int& count_numfail) {
    auto& t = ctx.thruster;
    auto& sp = ctx.solver;
    const double P_RFG_start = t.P_RFG;

    vector<double> prev, last_good;
    bool have_prev = false, have_good = false;
    double last_good_q0 = 0;

    for (int jj = 0; jj < sp.jjmax; ++jj) {
        t.Q0sccm = sp.Q0sccm_start + jj * sp.Q0sccm_step;
        t.Q0 = t.Q0sccm * PhysConst::SCCM_TO_PPS;
        cout << "Q0_STEP " << fixed << setprecision(4) << t.Q0sccm << " "
             << (jj + 1) << " " << sp.jjmax << endl;

        vector<double> guess;
        if (have_good && fabs(t.Q0sccm - last_good_q0) <= 20 * sp.Q0sccm_step) guess = last_good;
        else if (have_prev) guess = prev;
        else guess = gen_default_state(sys, t.Q0, t.V);

        GenPowerResult ps = (sp.solve_mode == 2)
            ? gen_solve_at_fixed_power(sys, ctx, t.P_RFG, guess)
            : gen_solve_for_target_current(sys, ctx, guess);

        if (sp.solve_mode == 1) {
            double I_found = isfinite(ps.I_mA) ? ps.I_mA : 0;
            double P_found = isfinite(ps.P_RFG_sol) ? ps.P_RFG_sol : 0;
            string status_str = ps.converged ? "converged" :
                (ps.fail_type == SolveFailType::NO_PHYSICAL_SOLUTION ? "above_P_max" : "numerical_fail");
            cout << "IFIX_RESULT " << fixed << setprecision(4) << t.Q0sccm << " "
                 << setprecision(2) << P_found << " " << sp.I_soll << " "
                 << I_found << " " << (I_found - sp.I_soll) << " " << status_str << endl;
        }

        bool brauchbar = ps.converged && (int)ps.state.size() == sys.state_size()
                         && ps.rf.valid && isfinite(ps.P_RFG_sol) && ps.P_RFG_sol > 0;
        if (!brauchbar) {
            string ft = (ps.fail_type == SolveFailType::NO_PHYSICAL_SOLUTION)
                        ? "NO_PHYSICAL_SOLUTION" : "NUMERICAL_FAIL";
            if (ps.fail_type == SolveFailType::NO_PHYSICAL_SOLUTION) {
                count_nophys++;
                cout << "NO_PHYSICAL_SOLUTION " << jj << " " << fixed << setprecision(4)
                     << t.Q0sccm << " " << ps.reason
                     << " I_best=" << (isfinite(ps.I_mA) ? ps.I_mA : -1)
                     << " P_max_tried=" << (isfinite(ps.P_trial_last) ? ps.P_trial_last : -1) << endl;
            } else {
                count_numfail++;
                cout << "NUMERICAL_FAIL " << jj << " " << fixed << setprecision(4)
                     << t.Q0sccm << " " << ps.reason << endl;
            }
            log_rows.push_back({jj, t.Q0sccm, ft, ft,
                isfinite(ps.P_trial_last) ? ps.P_trial_last : 0, isfinite(ps.I_mA) ? ps.I_mA : 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                ps.reason, false});
            log_events.push_back({jj, t.Q0sccm, ft + ": " + ps.reason});
            continue;
        }

        count_ok++;
        prev = ps.state; have_prev = true;
        if (isfinite(ps.inner_resid_norm) && ps.inner_resid_norm < sp.newton_tol) {
            last_good = ps.state; have_good = true; last_good_q0 = t.Q0sccm;
        }
        t.P_RFG = ps.P_RFG_sol;

        PlasmaState st = lumped_state(sys, ps.state);
        ExtractionResult ex = gen_extraction(sys, ctx, ps.state);
        DerivedQuantities dq = compute_derived_from(ctx, ex, st.n, st.ng, st.Te, st.Tg,
                                                     ps.rf.R_ind, ps.rf.I_coil, ps.rf.P_abs);

        if (isfinite(ps.inner_resid_norm) && ps.inner_resid_norm >= sp.newton_tol)
            cout << "SOFT_ACCEPT " << fixed << setprecision(4) << t.Q0sccm << " "
                 << scientific << ps.inner_resid_norm << endl;

        if (sp.solve_mode == 2) {
            cout << "PID_DONE " << fixed << setprecision(4) << ps.I_mA << " 0.0000 "
                 << ps.P_RFG_sol << " " << st.Te << " " << st.Tg << endl;
            cout << "CONVERGED " << ps.iterations << endl;
        }

        cout << "RESULT " << scientific << setprecision(4) << st.n << " " << st.ng << " "
             << fixed << setprecision(3) << st.Te << " " << st.Tg << " "
             << dq.I_extr_mA << " " << t.P_RFG << endl;

        // Aufgeschluesselt, was der zusammengefasste Zustand verdeckt
        for (int i = 0; i < (int)sys.species.size(); ++i)
            cout << "SPECIES_DENSITY " << sys.species[i].id << " "
                 << scientific << setprecision(4) << ps.state[i] << endl;
        for (const auto& sh : ex.ions)
            cout << "BEAM_SHARE " << sh.id << " Z " << sh.Z
                 << " I_mA " << fixed << setprecision(4) << sh.I_mA
                 << " Schub_mN " << sh.thrust_N * 1e3
                 << " v_aus " << sh.v_exhaust << endl;
        if (ex.limiting == "space_charge")
            cout << "BEAM_LIMIT space_charge Drosselung " << fixed << setprecision(4)
                 << ex.throttle << endl;

        cout << "RESULT_EXT " << fixed << setprecision(6) << t.Q0sccm << " " << t.P_RFG << " "
             << scientific << setprecision(4) << st.n << " " << st.ng << " " << st.n << " "
             << fixed << setprecision(4)
             << dq.T_i_N * 1e3 << " " << dq.T_n_N * 1e3 << " " << dq.T_total_N * 1e3 << " "
             << dq.icp_eff << " " << dq.gamma_eff << " " << dq.xi_mN_kW << " " << dq.eta_mass << endl;

        log_rows.push_back({jj, t.Q0sccm, "CONVERGED", "NONE",
            t.P_RFG, dq.I_extr_mA, st.Te, st.Tg, st.n, st.ng, ps.inner_resid_norm,
            dq.iondeg, ps.rf.P_abs, dq.cf, ps.rf.R_ind, ps.rf.I_coil,
            dq.eps_p_real, dq.eps_p_imag, dq.u_Bohm, dq.J_i, dq.pf, dq.P_RF,
            dq.T_i_N * 1e3, dq.T_n_N * 1e3, dq.T_total_N * 1e3,
            dq.icp_eff, dq.gamma_eff, dq.xi_mN_kW, dq.eta_mass,
            "ok", true});

        emit_csv_row(datei, "Generic", t.Q0sccm, st.n, st.ng, st.Te, st.Tg,
                     t.P_RFG, ps.rf.P_abs, ps.rf.R_ind, ps.rf.I_coil, ctx, dq);

        t.P_RFG = P_RFG_start;
    }
}

int main(int argc, char** argv) {
    string configFile = "params.txt";
    if (argc >= 2) configFile = argv[1];

    // Einziger SimContext -- kein globaler Zustand
    SimContext ctx;
    auto cd = loadConfig(configFile);
    applyConfig(ctx, cd);

    // Chemiepaket, falls die Konfiguration eines nennt. Scheitert das Laden,
    // laeuft der bisherige fest verdrahtete Weg weiter: ein fehlerhaftes Paket
    // darf einen Lauf nicht verhindern, muss aber auffallen.
    ChemSystem chem_sys;
    bool chem_geladen = false;
    if (!ctx.chem_package.empty()) {
        const string chem_pfad = resolve_chem_package(ctx.chem_package);
        ChemLoadResult cr = load_chem_package(chem_pfad);
        if (cr.ok) {
            chem_sys = cr.system;
            chem_geladen = true;
            adopt_package_properties(ctx, chem_sys);
        } else {
            cerr << "WARNUNG: Chemiepaket " << chem_pfad << " nicht geladen:" << endl
                 << cr.error_text() << endl
                 << "Falle auf die fest verdrahtete Physik zurueck." << endl;
        }
    }

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

        if (chem_geladen) {
            cout << "CHEM_PACKAGE " << chem_sys.source_path << " (" << chem_sys.name << ")" << endl;
            cout << "CHEM_SPECIES " << chem_sys.species.size()
                 << " CHEM_REACTIONS " << chem_sys.reactions.size() << endl;
            for (const auto& csp : chem_sys.species)
                cout << "CHEM_SPECIES_ITEM " << csp.id << " Ladung " << csp.charge
                     << " Masse " << scientific << setprecision(6) << csp.mass_kg << endl;
            cout << "SOLVER_PATH generic" << endl;
        } else {
            cout << "SOLVER_PATH legacy" << endl;
        }

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

        // Ratenanpassungen gegen die Gasauswahl absichern. Ein geladenes
        // Chemiepaket bringt seine eigenen Raten mit und ist davon nicht
        // betroffen.
        if (!chem_geladen) {
            auto probleme = rate_fit_gas_problems(ctx);
            if (!probleme.empty()) {
                ostream& aus = ctx.allow_foreign_rate_fits ? cout : cerr;
                aus << (ctx.allow_foreign_rate_fits ? "WARNUNG" : red + "FEHLER")
                    << ": fuer " << g.species << " liegen keine tabellierten Raten vor, "
                    << "und die eingebauten Anpassungen gelten nur fuer Xenon."
                    << (ctx.allow_foreign_rate_fits ? "" : reset) << endl;
                for (const auto& p : probleme) aus << "  - " << p << endl;
                if (!ctx.allow_foreign_rate_fits) {
                    cerr << "Abhilfe: Querschnittsdaten unter cross_sections/" << g.species
                         << "/ hinterlegen, ein Chemiepaket ueber chemistry_package waehlen, "
                         << "oder mit allow_foreign_rate_fits 1 bewusst mit Xenon-Raten rechnen."
                         << endl;
                    return 2;
                }
                cout << "FOREIGN_RATE_FITS " << g.species << endl;
            }
        }

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
        if (chem_geladen) {
            run_generic_sweep(ctx, chem_sys, datei, log_rows, log_events,
                              count_ok, count_nophys, count_numfail);
        } else
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
