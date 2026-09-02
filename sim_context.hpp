// sim_context.hpp -- Zentraler Simulationskontext. Traegt allen mutierbaren Zustand.
//
// SimContext ersetzt alle bisherigen globalen Variablen:
//   PhysConst::M/Eiz/...  -> ctx.gas.*
//   g_thruster / Const::R  -> ctx.thruster.*
//   Const::solve_mode/...  -> ctx.solver.*
//   g_kiz_table / ...      -> ctx.rates.*
//   debug_level            -> ctx.debug_level
//
// Wird in main() erzeugt und per Referenz an alle Funktionen uebergeben.
#ifndef SIM_CONTEXT_HPP
#define SIM_CONTEXT_HPP

#include "phys_const.hpp"
#include <vector>
#include <string>
#include <map>
#include <cmath>

// ═════════════════════════════════════════════════════════════
// Gasspezifische Eigenschaften
// ═════════════════════════════════════════════════════════════
struct GasProps {
    double M        = 2.1801711e-25;   // Atommasse [kg]
    double Eiz      = 1.943408035e-18; // Ionisierungsenergie [J]
    double Eexc     = 1.858524725e-18; // Anregungsenergie [J]
    double sigma_i  = 1.0e-18;        // Stossquerschnitt [m^2]
    double kappa    = 0.0057;          // Waermeleitfaehigkeit [W/(m*K)]
    std::string species = "xenon";

    // Die fest eingebauten Ratenanpassungen sind Xenon-Polynome aus Chabert
    // 2012. Fuer die anderen Gase gibt es sie nicht -- das haelt legacy_fits
    // fest, damit niemand still mit Xenon-Raten fuer Argon rechnet.
    bool has_legacy_fits = true;

    struct Entry { double M_kg, Eiz_J, Eexc_J, kappa_WmK; bool legacy_fits; };
    static inline const std::map<std::string, Entry> DB = {
        {"xenon",   {2.1801711e-25, 1.943408035e-18, 1.858524725e-18, 0.0057,  true}},
        {"krypton", {1.3914984e-25, 2.24009e-18,     1.60218e-18,     0.0094, false}},
        {"argon",   {6.6335209e-26, 2.52473e-18,     1.85853e-18,     0.0177, false}},
    };

    bool set_species(const std::string& name) {
        auto it = DB.find(name);
        if (it == DB.end()) return false;
        species = name;
        M = it->second.M_kg; Eiz = it->second.Eiz_J;
        Eexc = it->second.Eexc_J; kappa = it->second.kappa_WmK;
        has_legacy_fits = it->second.legacy_fits;
        return true;
    }
};

// ═════════════════════════════════════════════════════════════
// Triebwerksgeometrie + abgeleitete Groessen
// ═════════════════════════════════════════════════════════════
struct Thruster {
    // Primaer
    double R = 0.02, L = 0.04;
    double betai = 0.5, betag = 0.05145;
    double frequency = 2.5e6;
    double Nw = 6.0, R_ohm = 0.36, Rc = 0.02, lc = 0.04;
    double Vgrid = 1500.0, sgrid = 0.001;
    double eta_opt = 1.0;             // Grid-Optik-Effizienz [-] (1.0=Legacy, <1 bei explizitem Preset)
    double P_RFG = 18.0, Q0sccm = 0.475;
    // Abgeleitet
    double omega = 0, Q0 = 0, lambda_0 = 0, L_coil = 0;
    double A = 0, Ag = 0, Ai = 0, V = 0, k_0 = 0, J_CL = 0;

    void recompute(double M_gas) {
        using namespace PhysConst;
        omega    = 2.0 * pi * frequency;
        Q0       = Q0sccm * SCCM_TO_PPS;
        lambda_0 = R / 2.405 + L / pi;
        L_coil   = mu_0 * pi * Rc * Rc * Nw * Nw / lc;
        A        = 2.0 * pi * R * R + 2.0 * pi * R * L;
        Ag       = betag * pi * R * R;
        Ai       = betai * pi * R * R;
        V        = pi * R * R * L;
        k_0      = omega / c;
        J_CL     = (4.0/9.0) * epsilon0 * std::sqrt(2.0*e/M_gas) * std::pow(Vgrid, 1.5) / (sgrid*sgrid);
    }
};

// ═════════════════════════════════════════════════════════════
// Ratentabellen
// ═════════════════════════════════════════════════════════════
struct KizEntry { double Te_eV, Kiz; };
struct KelEntry { double Te_eV, Kel; };
struct KexEntry { double Te_eV, Kex_total, Pexc_coeff; };

struct RateConfig {
    int ionization_model = 0;   // 0=legacy, 1=tabulated
    int elastic_model    = 0;
    int excitation_model = 0;
    double Kex_scale     = 1.0;
    double kel_constant  = 1.0e-13;
    // Wurde der Wert bewusst gesetzt? Dann ist er eine Entscheidung des
    // Anwenders und keine stillschweigend uebernommene Xenon-Zahl.
    bool   kel_constant_explicit = false;

    std::vector<KizEntry> kiz;
    std::vector<KelEntry> kel;
    std::vector<KexEntry> kex;
};

// ═════════════════════════════════════════════════════════════
// Solver-Parameter
// ═════════════════════════════════════════════════════════════
struct SolverParams {
    int    solve_mode   = 1;
    double I_soll       = 15.0;
    double Q0sccm_start = 0.27;
    double Q0sccm_step  = 0.01;
    int    jjmax        = 73;
    double P_RFG_max    = 80.0;
    double P_abs_scale  = 1.0;
    double density_profile_factor = 1.0;
    // Energie, die ein Elektron-Ion-Paar an die Wand traegt, in Einheiten von
    // Te. Modell 1 rechnet sie aus dem Randschichtpotential, Modell 0 nimmt
    // die feste Zahl alpha_e_wall wie bisher.
    int    wall_energy_model = 1;
    double alpha_e_wall = 7.0;

    // Abbruchschranke auf dem groessten skalierten Residuum. Bei 1e-2 liegt
    // der Betriebspunkt noch nicht fest: die Ionisationsrate haengt
    // exponentiell von Te ab, ein Prozent Restfehler in der Leistungsbilanz
    // verschiebt die Dichte um mehrere Prozent.
    double newton_tol        = 1e-4;
    int    newton_max_iter   = 45;
    double newton_max_log_step = 0.35;
    double newton_fd_eps     = 1e-5;

    double power_tol_mA  = 0.05;
    int    power_max_iter = 35;
    double power_min     = 1.0;

    // Zustandsgrenzen. Sie sind Fangnetze fuer den Loeser, keine Physik: sie
    // sollen ein davonlaufendes Verfahren abfangen, aber niemals das Ergebnis
    // bestimmen. Sitzt eine konvergierte Loesung dicht an einer Grenze, wird
    // das gemeldet. Die Gastemperatur brauchte mehr Luft, weil der Weg zur
    // Loesung bei molekularen Gasen voruebergehend durch sehr heisse
    // Zwischenzustaende fuehrt, obwohl die Loesung selbst kuehl ist.
    double n_min  = 1e12, n_max  = 1e20;
    double ng_min = 1e16, ng_max = 1e22;
    double Te_min = 0.3,  Te_max = 20.0;
    double Tg_min = 200.0, Tg_max = 10000.0;

    int    ptc_max_iter     = 80;
    double ptc_start_gain   = 0.20;
    double ptc_min_gain     = 1e-4;
    double ptc_switch_merit = 5e-3;
    double ptc_accept_ratio = 0.98;

    double soft_abs_resid   = 5.0;
    double soft_rel_improve = 0.80;
};

// ═════════════════════════════════════════════════════════════
// Plasmazustand und RF-Zustand (Datenstrukturen)
// ═════════════════════════════════════════════════════════════
struct PlasmaState { double n, ng, Te, Tg; };

struct RFState {
    double P_abs = 0, R_ind = 0, I_coil = 0;
    bool valid = false;
};

struct DerivedQuantities {
    double I_extr_mA = 0, iondeg = 0, J_i = 0;
    double cf = 0, u_Bohm = 0, pf = 0;
    double eps_p_real = 0, eps_p_imag = 0;
    double zeta = 0, icp_eff = 0, P_RF = 0;
    double T_i_N = 0, T_n_N = 0, T_total_N = 0;
    double gamma_eff = 0, xi_mN_kW = 0, eta_mass = 0;
    // Extraktionsdiagnostik (neu)
    double I_CL_limit_mA = 0;     // Child-Langmuir Limit [mA]
    double I_Bohm_limit_mA = 0;   // Plasma-seitiges Limit [mA]
    double eta_opt_used = 0;       // Tatsaechlich verwendetes eta_opt
    std::string beam_limiting;     // "plasma" oder "space_charge"
};

// ═════════════════════════════════════════════════════════════
// SimContext -- gesamter Simulationszustand
// ═════════════════════════════════════════════════════════════
struct SimContext {
    GasProps     gas;
    Thruster     thruster;
    RateConfig   rates;
    SolverParams solver;

    std::string  cs_database     = "biagi";
    // Bewusst mit den Xenon-Anpassungen fuer ein anderes Gas rechnen.
    bool         allow_foreign_rate_fits = false;
    // Chemiepaket fuer den generischen Rechenweg. Leer heisst: fest
    // verdrahtete Xenon-Physik wie bisher.
    std::string  chem_package    = "";
    std::string  meta_preset_id  = "custom";
    int          rate_model   = 0;
    int          debug_level  = 2;

    // Abgeleitete Groessen nach Config-Aenderung aktualisieren
    void recompute() { thruster.recompute(gas.M); }
};

// ═════════════════════════════════════════════════════════════
// Solver-Ergebnistypen
// ═════════════════════════════════════════════════════════════
struct SolveResult {
    bool converged = false;
    bool soft_ok   = false;
    PlasmaState state{};
    RFState rf{};
    int iterations = 0;
    double resid_norm = 1e30;
    std::string reason;
};

enum class SolveFailType { NONE = 0, NO_PHYSICAL_SOLUTION, NUMERICAL_FAIL };

struct PowerSolveResult {
    bool converged = false;
    bool hit_limit = false;
    SolveFailType fail_type = SolveFailType::NONE;
    PlasmaState state{};
    RFState rf{};
    double P_RFG_sol = 0, P_trial_last = 0;
    double I_mA = 0, err_mA = 0;
    double inner_resid_norm = 1e30;
    int iterations = 0;
    std::string reason;
};

// ═════════════════════════════════════════════════════════════
// Log-Strukturen
// ═════════════════════════════════════════════════════════════
struct SimLogRow {
    int idx; double Q0sccm;
    std::string status, fail_type_str;
    double P_sol, I_mA, Te, Tg, n, ng, resid;
    double iondeg, P_abs, collision_freq, R_induktiv, I_coil;
    double eps_p_real, eps_p_imag, u_Bohm, J_i, plasmafrequenz, P_RF;
    double thrust_ions_mN, thrust_atoms_mN, thrust_total_mN;
    double icp_eff, gamma_eff, xi_mN_kW, eta_mass;
    std::string note; bool has_data;
};

struct SimLogEvent { int idx; double Q0sccm; std::string message; };

#endif // SIM_CONTEXT_HPP
