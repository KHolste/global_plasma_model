#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <stdexcept>
#include <chrono>
#include <cmath>
#include <complex>
#include <climits>
#include <cstdlib>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <map>
#include <set>
#include <string>
#include <limits>
#include "bessel_wrapper.hpp"  // Schlanker Wrapper statt voller 11k-Zeilen Bessel-Library

using namespace std::literals::complex_literals;
using namespace std;

// ============================================================
// Strukturierte Simulationslogdatei
// ============================================================
static std::string make_timestamp() {
    auto now = std::chrono::system_clock::now();
    auto t = std::chrono::system_clock::to_time_t(now);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y%m%d_%H%M%S", std::localtime(&t));
    return std::string(buf);
}
static std::string make_timestamp_readable() {
    auto now = std::chrono::system_clock::now();
    auto t = std::chrono::system_clock::to_time_t(now);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
    return std::string(buf);
}

struct SimLogRow {
    int idx;
    double Q0sccm;
    std::string status;       // CONVERGED / NO_PHYSICAL_SOLUTION / NUMERICAL_FAIL
    std::string fail_type;    // NONE / NO_PHYSICAL_SOLUTION / NUMERICAL_FAIL
    // Primaere Zustandsgroessen
    double P_sol, I_mA, Te, Tg, n, ng, resid;
    // Abgeleitete Groessen (alle nur bei has_data=true gueltig)
    double iondeg, P_abs, collision_freq, R_induktiv, I_coil;
    double eps_p_real, eps_p_imag, u_Bohm, J_i, plasmafrequenz, P_RF;
    // Performance-Groessen
    double thrust_ions_mN, thrust_atoms_mN, thrust_total_mN;
    double icp_eff, gamma_eff, xi_mN_kW, eta_mass;
    std::string note;
    bool has_data;
};

struct SimLogEvent {
    int idx;
    double Q0sccm;
    std::string message;
};

static std::vector<SimLogRow> g_log_rows;
static std::vector<SimLogEvent> g_log_events;

static void simlog_add_event(int idx, double q0, const std::string& msg) {
    g_log_events.push_back({idx, q0, msg});
}

enum class DebugLevel { OFF=0, ERROR=1, INFO=2, DETAIL=3 };
static int debug_level = 2;
static std::string debug_log_path = "chabert_debug.log";

static inline const char* dbg_level_name(int lvl) {
    switch (lvl) {
        case 1: return "ERROR";
        case 2: return "INFO";
        case 3: return "DETAIL";
        default: return "OFF";
    }
}

static void debug_emit(int lvl, const std::string& tag, const std::string& msg, bool also_stderr=false) {
    if (debug_level < lvl) return;
    std::ostringstream os;
    os << "DBG " << tag << " " << msg;
    const std::string line = os.str();
    std::cout << line << std::endl;
    std::cout.flush();
    if (also_stderr || lvl <= 1) {
        std::cerr << line << std::endl;
        std::cerr.flush();
    }
    std::ofstream log(debug_log_path, std::ios::app);
    if (log) log << line << '\n';
}

#define DBG(lvl, expr) \
    do { if (debug_level >= (lvl)) { \
        std::ostringstream _dbg_os; _dbg_os << expr; \
        debug_emit((lvl), "TRACE", _dbg_os.str()); \
    } } while(0)

static std::string state_summary(double n, double ng, double Te, double Tg, double Ug, double Ue) {
    std::ostringstream os;
    os << std::scientific << std::setprecision(6)
       << "n=" << n << " ng=" << ng << " Ug=" << Ug << " Ue=" << Ue << " "
       << std::fixed << std::setprecision(6)
       << "Te=" << Te << " Tg=" << Tg;
    return os.str();
}

// ============================================================
// Config-Einlese-Funktion
// ============================================================
struct ConfigData {
    std::map<std::string, double> numeric;
    std::map<std::string, std::string> strings;
};

ConfigData loadConfig(const std::string& filename) {
    ConfigData cd;
    std::ifstream f(filename);
    if (!f) {
        std::cerr << "[Config] Datei '" << filename << "' nicht gefunden – Standardwerte werden verwendet." << std::endl;
        return cd;
    }
    // String-Parameter, die nicht als double geparsed werden
    static const std::set<std::string> STRING_KEYS = {"gas_species"};
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string key, val_str;
        if (!(ss >> key)) continue;
        if (STRING_KEYS.count(key)) {
            if (ss >> val_str) cd.strings[key] = val_str;
        } else {
            double val;
            if (ss >> val) cd.numeric[key] = val;
        }
    }
    return cd;
}
#define CFG_LOAD(cfg, name) if ((cfg).count(#name)) (name) = (cfg).at(#name)

// ============================================================
// Physikalische Konstanten und Triebwerksparameter
// ============================================================
// ============================================================
// Gruppe A: Physikalische Naturkonstanten (unveraenderlich)
// ============================================================
namespace PhysConst {
    constexpr double me       = 9.10938215e-31;    // Elektronenmasse [kg]
    constexpr double e        = 1.602176487e-19;    // Elementarladung [C]
    constexpr double kB       = 1.3806504e-23;      // Boltzmann-Konstante [J/K]
    constexpr double epsilon0 = 8.854187e-12;       // Vakuumpermittivitaet [F/m]
    constexpr double pi       = 3.141592653589793;
    constexpr double c        = 299792458.0;        // Lichtgeschwindigkeit [m/s]
    constexpr double mu_0     = 1.256637061e-6;     // Vakuumpermeabilitaet [H/m]

    constexpr double SCCM_TO_PPS = 4.477962312e17;  // sccm → particles/s

    constexpr double Tg0    = 293.00;               // Wandtemperatur [K]
    constexpr double conv   = 11604.5250061657;     // eV → K Umrechnung

    // Gasspezifische Konstanten – Xenon-Defaults (werden bei anderem gas_species ueberschrieben)
    double M        = 2.1801711e-25;      // Atommasse [kg] (Xe: 131.293 u)
    double Eiz      = 1.943408035e-18;    // Ionisierungsenergie [J] (Xe: 12.127 eV)
    double Eexc     = 1.858524725e-18;    // Anregungsenergie [J] (Xe: 11.6 eV)
    double sigma_i  = 1.0e-18;            // Stossquerschnitt [m²]
    double sigmae   = 1.0e-18;            // Elastischer Stossquerschnitt [m²]
    double kappa    = 0.0057;             // Waermeleitfaehigkeit [W/(m·K)]

    // Gasspezifische Lookup: (M_kg, Eiz_J, Eexc_J, kappa_W_mK)
    // Eexc: niedrigste dominante Anregungsenergie fuer Legacy-Arrhenius-Fit
    struct GasProperties {
        double M_kg;
        double Eiz_J;      // Ionisierungsenergie [J]
        double Eexc_J;     // Anregungsenergie [J] (Legacy)
        double kappa_WmK;  // Waermeleitfaehigkeit [W/(m·K)]
    };

    // Datenquelle: NIST, Lieberman & Lichtenberg
    static const std::map<std::string, GasProperties> GAS_DB = {
        {"xenon",   {2.1801711e-25, 1.943408035e-18, 1.858524725e-18, 0.0057}},  // 131.293 u, 12.127 eV, 11.6 eV
        {"krypton", {1.3914984e-25, 2.24009e-18,     1.60218e-18,     0.0094}},  // 83.798 u, 13.9996 eV, 10.0 eV
        {"argon",   {6.6335209e-26, 2.52473e-18,     1.85853e-18,     0.0177}},  // 39.948 u, 15.7596 eV, 11.6 eV
    };

    void set_gas(const std::string& gas_name) {
        auto it = GAS_DB.find(gas_name);
        if (it == GAS_DB.end()) {
            std::cerr << "FEHLER: Unbekanntes Gas '" << gas_name << "', verwende Xenon-Defaults" << std::endl;
            return;
        }
        const auto& gp = it->second;
        M       = gp.M_kg;
        Eiz     = gp.Eiz_J;
        Eexc    = gp.Eexc_J;
        kappa   = gp.kappa_WmK;
    }
}

// ============================================================
// Gruppe B+C: Triebwerksparameter + abgeleitete Groessen
// ============================================================
struct ThrusterParams {
    // --- Primaere Eingabeparameter (Gruppe B) ---
    double R        = 0.02;       // Kammerradius [m]
    double L        = 0.04;       // Kammerlaenge [m]
    double betai    = 0.5;        // Ionen-Gittertransparenz [-]
    double betag    = 0.05145;    // Gas-Gittertransparenz [-]
    double frequency= 2.5e6;     // RF-Anregefrequenz [Hz]
    double Nw       = 6.00;       // Spulenwindungen [-]
    double R_ohm    = 0.36;       // Ohmscher Spulenwiderstand [Ohm]
    double Rc       = 0.02;       // Spulenradius [m]
    double lc       = 0.04;       // Spulenlaenge [m]
    double Vgrid    = 1500.00;    // Gitterspannung [V]
    double sgrid    = 0.001;      // Gitterabstand [m]
    double P_RFG    = 18.00;      // RF-Generatorleistung [W]
    double Q0sccm   = 0.475;      // Massenfluss [sccm]

    // --- Abgeleitete Groessen (Gruppe C, aus B berechnet) ---
    double omega    = 0.0;        // Kreisfrequenz [rad/s]
    double Q0       = 0.0;        // Teilchenfluss [s⁻¹]
    double lambda_0 = 0.0;        // Thermische Diffusionslaenge [m]
    double L_coil   = 0.0;        // Spuleninduktivitaet [H]
    double A        = 0.0;        // Gesamte Wandflaeche [m²]
    double Ag       = 0.0;        // Gasaustrittsflaeche [m²]
    double Ai       = 0.0;        // Ionenaustrittsflaeche [m²]
    double V        = 0.0;        // Kammervolumen [m³]
    double k_0      = 0.0;        // Vakuum-Wellenzahl [1/m]
    double J_CL     = 0.0;        // Child-Langmuir-Stromdichte [A/m²]

    // Berechne alle abgeleiteten Groessen aus den Primaerparametern.
    // Muss nach jeder Aenderung eines Primaerparameters aufgerufen werden.
    void recompute_derived() {
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
        J_CL     = (4.0/9.0) * epsilon0 * sqrt(2.0*e/M) * pow(Vgrid, 1.5) / (sgrid*sgrid);
    }

    // Lade Primaerparameter aus Config und berechne Abgeleitete
    void load_from_config(const std::map<std::string, double>& cfg) {
        CFG_LOAD(cfg, R);        CFG_LOAD(cfg, L);        CFG_LOAD(cfg, betai);
        CFG_LOAD(cfg, betag);    CFG_LOAD(cfg, frequency); CFG_LOAD(cfg, Nw);
        CFG_LOAD(cfg, R_ohm);    CFG_LOAD(cfg, Rc);       CFG_LOAD(cfg, lc);
        CFG_LOAD(cfg, Vgrid);    CFG_LOAD(cfg, sgrid);    CFG_LOAD(cfg, P_RFG);
        CFG_LOAD(cfg, Q0sccm);
        recompute_derived();
    }
};

// Globale Instanz (Rueckwaertskompatibilitaet: wird ueber Const-Namespace gespiegelt)
static ThrusterParams g_thruster = [](){
    ThrusterParams tp;
    tp.recompute_derived();  // Abgeleitete Groessen aus Standardwerten berechnen
    return tp;
}();

// ============================================================
// Gruppe D: Solver-/Sweep-Parameter + Rueckwaertskompatibilitaet
// ============================================================
namespace Const {
    // Physikalische Konstanten (Aliase fuer Rueckwaertskompatibilitaet)
    using namespace PhysConst;

    // Triebwerksparameter: Referenzen auf die globale ThrusterParams-Instanz.
    // Bestehender Code kann weiterhin 'R', 'L', 'V' etc. ohne Qualifizierung nutzen.
    // Langfristig sollten Funktionen const ThrusterParams& entgegennehmen.
    double& R         = g_thruster.R;
    double& L         = g_thruster.L;
    double& betai     = g_thruster.betai;
    double& betag     = g_thruster.betag;
    double& frequency = g_thruster.frequency;
    double& omega     = g_thruster.omega;
    double& Nw        = g_thruster.Nw;
    double& R_ohm     = g_thruster.R_ohm;
    double& Rc        = g_thruster.Rc;
    double& lc        = g_thruster.lc;
    double& Vgrid     = g_thruster.Vgrid;
    double& sgrid     = g_thruster.sgrid;
    double& P_RFG     = g_thruster.P_RFG;
    double& Q0sccm    = g_thruster.Q0sccm;
    double& Q0        = g_thruster.Q0;
    double& lambda_0  = g_thruster.lambda_0;
    double& L_coil    = g_thruster.L_coil;
    double& A         = g_thruster.A;
    double& Ag        = g_thruster.Ag;
    double& Ai        = g_thruster.Ai;
    double& V         = g_thruster.V;
    double& k_0       = g_thruster.k_0;
    double& J_CL      = g_thruster.J_CL;

    // Solver-Parameter (bleiben vorerst direkt im Namespace)
    int    solve_mode   = 1;
    double I_soll       = 15.00;
    double Q0sccm_start = 0.27;
    double Q0sccm_step  = 0.01;
    int    jjmax        = 73;
    int    use_paper_kel = 0;
    double kel_constant  = 1.0e-13;  // Konstanter Kel-Wert [m^3/s] wenn use_paper_kel=1

    // Gasspezies fuer Cross-Section-Pfadaufloesung
    // Bestimmt: cross_sections/<gas_species>/{kel,kiz,kex}_table.csv
    // und die physikalischen Konstanten (M, Eiz, Eexc, kappa)
    std::string gas_species = "xenon";

    // Ratenmodell-Preset (bequeme Gesamtsteuerung):
    // 0 = legacy (alle Raten aus Chabert-Fits, paper-kompatibel)
    // 1 = conservative_tabulated (Kiz+Kex tabelliert, Kel legacy)
    // 2 = full_tabulated (alle Raten aus Biagi/LXCat-Tabellen)
    // Wird in applyConfig auf die Einzelschalter abgebildet.
    // Einzelschalter koennen danach noch ueberschrieben werden.
    int    rate_model = 0;

    // Einzelschalter (werden von rate_model gesetzt, aber auch einzeln konfigurierbar)
    int    ionization_model = 0;  // 0=legacy, 1=tabulated
    int    elastic_model    = 0;  // 0=legacy, 1=tabulated
    int    excitation_model = 0;  // 0=legacy, 1=tabulated

    // Diagnostische Skalierungsfaktoren (1.0 = unveraendert).
    double P_abs_scale = 1.0;
    double Kex_scale   = 1.0;   // 0.0 = Anregung deaktiviert (wie Xenon_Code.py)

    // Dichte-Profilkorrekturfaktor (Dietz et al. 2021, S.25):
    // Das 0D-Modell ueberschaetzt die mittlere Elektronendichte, weil es ein
    // homogenes Profil annimmt. Die effektive mittlere Dichte ist n_eff = factor * n.
    // factor = 1.0: unveraendertes 0D-Verhalten (Standard)
    // factor = 0.82: Korrektur nach Dietz (3D-Mittelwert / Zentralwert)
    // Wirkt NUR auf Volumenverlustterme (Ionisation, Anregung, elastische Stoesse),
    // NICHT auf Wandfluss (Bohm-Fluss) oder Randbedingungen.
    double density_profile_factor = 1.0;

    // Elektronen-Wandverlustfaktor in P5 (Eq. 13 im Paper):
    //   P5 = alpha_e * kB * Te * n * uB * Aeff / V
    // Der Faktor alpha_e fasst kinetische Energie (5/2 kBTe) + Sheath-Potenzial
    // (typisch ~2 kBTe bei Xe) + Ionenenergie zusammen.
    // Lieberman & Lichtenberg: alpha_e ≈ 5.2 (ohne Ionen) bis 7.2 (mit Ionen)
    // Chabert 2012: alpha_e = 7 (Eq. 13)
    double alpha_e_wall = 7.0;

    double P_RFG_max    = 80.00;
    int    newton_max_iter   = 45;
    double newton_tol        = 1e-2;
    double power_tol_mA      = 0.05;
    int    power_max_iter    = 35;
    double power_min         = 1.0;

    double n_min   = 1e12, n_max   = 1e20;
    double ng_min  = 1e16, ng_max  = 1e22;
    double Te_min  = 0.3,  Te_max  = 20.0;
    double Tg_min  = 200.0, Tg_max = 2500.0;
    double newton_max_log_step = 0.35;
    double newton_fd_eps = 1e-5;

    int    ptc_max_iter       = 80;
    double ptc_start_gain     = 0.20;
    double ptc_min_gain       = 1e-4;
    double ptc_switch_merit   = 5e-3;
    double ptc_accept_ratio   = 0.98;

    double stationary_soft_abs_resid = 5.0;
    double stationary_soft_rel_improve = 0.80;

    void applyConfig(const ConfigData& cd) {
        const auto& cfg = cd.numeric;

        // gas_species (String-Parameter) → physikalische Konstanten setzen
        if (cd.strings.count("gas_species")) {
            gas_species = cd.strings.at("gas_species");
            PhysConst::set_gas(gas_species);
            g_thruster.recompute_derived();  // M hat sich ggf. geaendert (J_CL)
        }

        // Triebwerksparameter → ThrusterParams-Struct
        g_thruster.load_from_config(cfg);

        // Solver-Parameter (bleiben im Namespace)
        CFG_LOAD(cfg, P_RFG_max);
        CFG_LOAD(cfg, I_soll);
        CFG_LOAD(cfg, Q0sccm_start); CFG_LOAD(cfg, Q0sccm_step);
        CFG_LOAD(cfg, power_tol_mA); CFG_LOAD(cfg, power_min);
        CFG_LOAD(cfg, n_min); CFG_LOAD(cfg, n_max); CFG_LOAD(cfg, ng_min); CFG_LOAD(cfg, ng_max);
        CFG_LOAD(cfg, Te_min); CFG_LOAD(cfg, Te_max); CFG_LOAD(cfg, Tg_min); CFG_LOAD(cfg, Tg_max);
        CFG_LOAD(cfg, newton_max_log_step); CFG_LOAD(cfg, newton_fd_eps);
        if (cfg.count("jjmax"))  jjmax  = (int)cfg.at("jjmax");
        if (cfg.count("solve_mode")) solve_mode = std::max(1, std::min(2, (int)cfg.at("solve_mode")));
        if (cfg.count("use_paper_kel")) use_paper_kel = ((int)cfg.at("use_paper_kel") != 0) ? 1 : 0;
        CFG_LOAD(cfg, kel_constant);
        CFG_LOAD(cfg, alpha_e_wall);
        CFG_LOAD(cfg, P_abs_scale);
        // rate_model Preset: setzt alle Einzelschalter auf einmal
        if (cfg.count("rate_model")) {
            rate_model = (int)cfg.at("rate_model");
            if (rate_model == 1) { ionization_model = 1; excitation_model = 1; elastic_model = 0; }
            else if (rate_model == 2) { ionization_model = 1; excitation_model = 1; elastic_model = 1; }
            else { ionization_model = 0; excitation_model = 0; elastic_model = 0; }
        }
        // Einzelschalter koennen das Preset danach ueberschreiben
        if (cfg.count("ionization_model")) ionization_model = (int)cfg.at("ionization_model");
        if (cfg.count("elastic_model")) elastic_model = (int)cfg.at("elastic_model");
        if (cfg.count("excitation_model")) excitation_model = (int)cfg.at("excitation_model");
        CFG_LOAD(cfg, Kex_scale);
        CFG_LOAD(cfg, density_profile_factor);
        if (cfg.count("newton_max_iter")) newton_max_iter = (int)cfg.at("newton_max_iter");
        if (cfg.count("power_max_iter"))  power_max_iter  = (int)cfg.at("power_max_iter");
        if (cfg.count("debug_level"))     debug_level = std::max(0, std::min(3, (int)cfg.at("debug_level")));
    }
}
using namespace Const;

// ============================================================
// Tabellierte Stossraten (Biagi/LXCat)
// ============================================================

// Kiz-Tabelle
struct KizTableEntry { double Te_eV, Kiz; };
static std::vector<KizTableEntry> g_kiz_table;

static bool load_kiz_table(const std::string& path) {
    std::ifstream f(path);
    if (!f) return false;
    g_kiz_table.clear();
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == 'T') continue;
        std::istringstream ss(line);
        KizTableEntry entry{};
        char comma;
        if (ss >> entry.Te_eV >> comma >> entry.Kiz) {
            g_kiz_table.push_back(entry);
        }
    }
    return !g_kiz_table.empty();
}

static double interpolate_kiz_table(double Te_eV) {
    if (g_kiz_table.empty() || Te_eV <= 0.0) return 0.0;
    if (Te_eV <= g_kiz_table.front().Te_eV) return g_kiz_table.front().Kiz;
    if (Te_eV >= g_kiz_table.back().Te_eV)  return g_kiz_table.back().Kiz;
    size_t lo = 0, hi = g_kiz_table.size() - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (g_kiz_table[mid].Te_eV <= Te_eV) lo = mid; else hi = mid;
    }
    double t = (Te_eV - g_kiz_table[lo].Te_eV) / (g_kiz_table[hi].Te_eV - g_kiz_table[lo].Te_eV);
    return g_kiz_table[lo].Kiz * (1.0 - t) + g_kiz_table[hi].Kiz * t;
}

// Kel-Tabelle
struct KelTableEntry { double Te_eV, Kel; };
static std::vector<KelTableEntry> g_kel_table;

static bool load_kel_table(const std::string& path) {
    std::ifstream f(path);
    if (!f) return false;
    g_kel_table.clear();
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == 'T') continue;
        std::istringstream ss(line);
        KelTableEntry entry{};
        char comma;
        if (ss >> entry.Te_eV >> comma >> entry.Kel) {
            g_kel_table.push_back(entry);
        }
    }
    return !g_kel_table.empty();
}

static double interpolate_kel_table(double Te_eV) {
    if (g_kel_table.empty() || Te_eV <= 0.0) return 0.0;
    if (Te_eV <= g_kel_table.front().Te_eV) return g_kel_table.front().Kel;
    if (Te_eV >= g_kel_table.back().Te_eV)  return g_kel_table.back().Kel;
    size_t lo = 0, hi = g_kel_table.size() - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (g_kel_table[mid].Te_eV <= Te_eV) lo = mid; else hi = mid;
    }
    double t = (Te_eV - g_kel_table[lo].Te_eV) / (g_kel_table[hi].Te_eV - g_kel_table[lo].Te_eV);
    return g_kel_table[lo].Kel * (1.0 - t) + g_kel_table[hi].Kel * t;
}

// Kex-Tabelle
struct KexTableEntry { double Te_eV, Kex_total, Pexc_coeff; };
static std::vector<KexTableEntry> g_kex_table;

static bool load_kex_table(const std::string& path) {
    std::ifstream f(path);
    if (!f) return false;
    g_kex_table.clear();
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == 'T') continue;
        std::istringstream ss(line);
        KexTableEntry entry{};
        char comma;
        if (ss >> entry.Te_eV >> comma >> entry.Kex_total >> comma >> entry.Pexc_coeff) {
            g_kex_table.push_back(entry);
        }
    }
    return !g_kex_table.empty();
}

// Lineare Interpolation in der Kex-Tabelle
static double interpolate_kex_table(double Te_eV, bool use_pexc_coeff) {
    if (g_kex_table.empty() || Te_eV <= 0.0) return 0.0;
    if (Te_eV <= g_kex_table.front().Te_eV)
        return use_pexc_coeff ? g_kex_table.front().Pexc_coeff : g_kex_table.front().Kex_total;
    if (Te_eV >= g_kex_table.back().Te_eV)
        return use_pexc_coeff ? g_kex_table.back().Pexc_coeff : g_kex_table.back().Kex_total;
    // Binaere Suche
    size_t lo = 0, hi = g_kex_table.size() - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (g_kex_table[mid].Te_eV <= Te_eV) lo = mid; else hi = mid;
    }
    double t = (Te_eV - g_kex_table[lo].Te_eV) / (g_kex_table[hi].Te_eV - g_kex_table[lo].Te_eV);
    if (use_pexc_coeff)
        return g_kex_table[lo].Pexc_coeff * (1.0 - t) + g_kex_table[hi].Pexc_coeff * t;
    return g_kex_table[lo].Kex_total * (1.0 - t) + g_kex_table[hi].Kex_total * t;
}

// ============================================================
// Hilfsfunktionen
// ============================================================
double vg(double Tg)       { return sqrt(8 * kB * Tg / (pi * M)); }
double vi(double Ti)       { return sqrt(8 * kB * Ti / (pi * M)); }
double Gamma_g(double ng, double vg_val) { return 0.25 * ng * vg_val; }

double Kex(double Te) { return Kex_scale * 1.2921e-13 * exp(-e * 11.6 / (kB * Te * conv)); }

double Kiz(double Te) {
    if (ionization_model == 1 && !g_kiz_table.empty()) {
        return interpolate_kiz_table(Te);  // Te ist in eV
    }
    // Legacy: Chabert Polynomfit
    double TeV = kB * Te * conv / e;
    double K1 = 6.73e-15 * sqrt(TeV) * (3.97 + 0.643*TeV - 0.0368*TeV*TeV) * exp(-12.127/TeV);
    double K2 = 6.73e-15 * sqrt(TeV) * (-0.0001031*TeV*TeV + 6.386*exp(-12.127/TeV));
    return 0.5 * (K1 + K2);
}

// Elastischer Elektron-Neutral-Stosskoeffizient [m³/s]
//
// elastic_model == 1 (tabulated):
//   Biagi/LXCat Momentum-Transfer-Querschnitt, Maxwell-Boltzmann-integriert.
//   Kel(3eV) ≈ 2.3e-13, Kel(5eV) ≈ 2.7e-13 m³/s fuer Xenon.
//   Physikalisch am genauesten (Audit 2026-03-28 bestaetigt).
//
// elastic_model == 0 (legacy):
//   Kel = kel_constant (default 1e-13 m³/s), wie in Chabert et al.,
//   Phys. Plasmas 19, 073512 (2012), S.2.
//   Unterschaetzt Kel im Betriebsbereich um Faktor 2-3 (Audit-Ergebnis).
//
// Kel wirkt auf:
//   - Kollisionsfrequenz nu_m = Kel * ng  (RF-Kopplung, eps_p)
//   - Gasheizung (Pg1 in Eq.11)
//   - Elektronenkuehlung (P4 in Eq.13)
double Kel(double Te) {
    double val;
    if (elastic_model == 1 && !g_kel_table.empty()) {
        val = interpolate_kel_table(Te);  // Te ist in eV
    } else {
        // Legacy: konstanter Wert (Paper-kompatibel)
        val = kel_constant;
    }
    // Sicherheitspruefung: Kel muss positiv sein
    if (val <= 0.0) {
        cerr << "FEHLER: Kel(" << Te << ") = " << val << " <= 0! "
             << "elastic_model=" << elastic_model
             << ", kel_constant=" << kel_constant << endl;
        val = 1e-15;  // Minimaler Sicherheitswert, keine stille 0
    }
    return val;
}

double uB(double Te)         { return sqrt(kB * Te * conv / M); }
double lambda_i(double ng)   { return 1.0 / (ng * sigma_i); }
double plasma_freq(double n) { return sqrt(n * e * e / (me * epsilon0)); }
double coll_freq(double ng, double Te) { return Kel(Te) * ng; }

double Aeff(double lambda) {
    double hL = 0.86 * pow(3 + L/(2*lambda), -0.5);
    double hR = 0.80 * pow(4 + R/lambda,     -0.5);
    return 2*hR*pi*R*L + 2*hL*pi*R*R;
}
double Aeff1(double lambda) {
    double hL = 0.86 * pow(3 + L/(2*lambda), -0.5);
    double hR = 0.80 * pow(4 + R/lambda,     -0.5);
    return 2*hR*pi*R*L + (2 - betai)*hL*pi*R*R;
}
double Gamma_i_func(double lambda, double Te, double n) {
    double hL = 0.86 * pow(3 + L/(2*lambda), -0.5);
    return hL * n * uB(Te);
}
double Thrust(double Gamma_i_val) {
    return Gamma_i_val * M * Ai * sqrt(2 * e * Vgrid / M);
}
complex<double> my_calc_eps_p(double n, double ng, double Te) {
    double a = plasma_freq(n), b = omega, cf = coll_freq(ng, Te);
    return complex<double>(1 - a*a/(b*b+cf*cf), -a*a*cf/(b*(b*b+cf*cf)));
}
template <typename T>
int signum(T val) { return (T(0) < val) - (val < T(0)); }


struct StationaryRFState {
    double P_abs = std::numeric_limits<double>::quiet_NaN();
    double R_ind = std::numeric_limits<double>::quiet_NaN();
    double I_coil = std::numeric_limits<double>::quiet_NaN();
    bool valid = false;
};

struct StationaryPlasmaState {
    double n;
    double ng;
    double Te;
    double Tg;
};

static StationaryRFState stationary_compute_rf(double n, double ng, double Te, double P_RFG_local) {
    StationaryRFState out{};
    if (!(std::isfinite(n) && std::isfinite(ng) && std::isfinite(Te) && std::isfinite(P_RFG_local)) ||
        n <= 0.0 || ng <= 0.0 || Te <= 0.0 || P_RFG_local <= 0.0) {
        return out;
    }

    double a = plasma_freq(n), b = omega, cf = coll_freq(ng, Te);
    double denom = b*b + cf*cf;
    if (!std::isfinite(a) || !std::isfinite(cf) || denom <= 0.0 || !std::isfinite(denom)) return out;

    double it = a*a*cf / (b*denom);
    double rt = 1.0 - a*a / denom;

    double rabs   = sqrt(it*it + rt*rt);
    double arg_rt = max(0.0, (rabs + rt) / 2.0);
    double arg_it = max(0.0, (rabs - rt) / 2.0);
    double o_rt   = sqrt(arg_rt);
    double o_it   = signum(it) * sqrt(arg_it);

    complex<double> k_1(o_rt*k_0, -o_it*k_0);
    complex<double> kR = R * k_1;
    complex<double> eps_p(rt, -it);
    complex<double> j0 = bessel_cyl_j(0, kR);
    complex<double> j1 = bessel_cyl_j(1, kR);
    complex<double> denom_c = eps_p * j0;

    // Numerische Stabilisierung: |eps_p * J0| kann nahe Null werden bei
    // omega_p ≈ omega (Plasmareso nanz). Statt NaN/Inf zu erzeugen, wird
    // der Nenner auf einen Minimalwert begrenzt. Das haelt R_ind endlich
    // und P_abs stetig differenzierbar — entscheidend fuer den FD-Jacobian.
    if (!std::isfinite(denom_c.real()) || !std::isfinite(denom_c.imag())) return out;
    double denom_c_abs = std::abs(denom_c);
    if (denom_c_abs < 1e-30) {
        // Richtung beibehalten, Betrag auf Minimum setzen
        denom_c = denom_c / std::max(denom_c_abs, 1e-300) * 1e-30;
    }

    complex<double> result = (1i * kR * j1) / denom_c;
    double R_ind_raw = 2*pi*Nw*Nw / (epsilon0*L*omega) * result.real();

    // R_ind clampen statt abbrechen: Bei omega_p ≈ omega kann R_ind
    // negativ werden (Vorzeichenwechsel von result.real()). Ein harter
    // Abbruch erzeugt NaN im FD-Jacobian, weil die perturbierte Seite
    // (n+h vs n-h) den Sprung trifft. Stattdessen: stetige Begrenzung.
    // R_ind_min ≈ 1e-4 Ohm: physikalisch vernachlaessigbar gegenueber
    // R_ohm = 0.36 Ohm, aber numerisch stabil.
    const double R_ind_min = 1e-4;
    double R_ind = std::max(R_ind_min, R_ind_raw);
    if (!std::isfinite(R_ind)) return out;

    double Ic = sqrt(2.0 * P_RFG_local / (R_ind + R_ohm));
    double P_abs = 0.5 * R_ind * Ic * Ic;
    if (!std::isfinite(Ic) || !std::isfinite(P_abs)) return out;

    out.P_abs = P_abs;
    out.R_ind = R_ind;
    out.I_coil = Ic;
    out.valid = true;
    return out;
}

static bool stationary_state_finite_positive(const StationaryPlasmaState& s) {
    return std::isfinite(s.n) && std::isfinite(s.ng) && std::isfinite(s.Te) && std::isfinite(s.Tg)
        && s.n > 0.0 && s.ng > 0.0 && s.Te > 0.0 && s.Tg > 0.0;
}

// Maximaler Ionisierungsgrad n/ng fuer Zustandsgrenzen-Pruefung.
// Numerische Schutzgrenze — verhindert, dass der LM-Solver in
// Regionen mit n > ng gelangt (physikalisch unplausibel fuer
// schwach ionisiertes Plasma). Wert 0.95 statt 0.5, weil:
// - Typischer Betrieb: iondeg < 10% (Paper: ~7% am CL-Limit)
// - Hohe Leistung / niedriger Fluss: bis ~30-50% denkbar
// - 0.5 verwarf den Solver-Suchraum unnoetig bei Hochleistung
// - 0.95 erlaubt den Solver, durch Hochionisierungs-Regionen
//   zu iterieren, verhindert aber n ≈ ng (numerische Probleme
//   in der Neutralgasbilanz, da Nenner → 0)
constexpr double iondeg_max = 0.95;

static bool stationary_state_in_bounds(const StationaryPlasmaState& s) {
    if (!stationary_state_finite_positive(s)) return false;
    if (s.n  < n_min  || s.n  > n_max)  return false;
    if (s.ng < ng_min || s.ng > ng_max) return false;
    if (s.Te < Te_min || s.Te > Te_max) return false;
    if (s.Tg < Tg_min || s.Tg > Tg_max) return false;
    double iondeg = s.n / std::max(1.0, s.ng);
    if (!std::isfinite(iondeg) || iondeg < 0.0 || iondeg > iondeg_max) return false;
    return true;
}

static StationaryPlasmaState stationary_safe_defaults_for_q(double q0_particles) {
    StationaryPlasmaState s;
    double p_null = 4 * kB * Tg0 * q0_particles / (vg(Tg0) * Ag);
    s.ng = max(1e16, p_null / (kB * Tg0));
    s.n  = 1.0e17;
    s.Te = 3.75;
    s.Tg = 300.0;
    return s;
}

// ============================================================
// Abgeleitete Groessen — einheitlich fuer alle Solver-Pfade
//
// Wird von compute_derived() aus dem konvergierten Plasmazustand
// (n, ng, Te, Tg) und den RF-Groessen (R_ind, I_coil, P_abs)
// berechnet. Wird vom stationaeren Solver und der Ausgabe verwendet.
//
// Referenz: Chabert et al., Phys. Plasmas 19, 073512 (2012)
// ============================================================
struct DerivedQuantities {
    // --- Teilchenstroeme und Dichten ---
    double I_extr_mA;    // [mA]         Strahlstrom: Ai * e * Gamma_i * 1000
    double iondeg;       // [%]          Ionisierungsgrad: n / ng * 100
    double J_i;          // [A/m²]       Ionenstromdichte: e * Gamma_i

    // --- Plasmaparameter ---
    double cf;           // [s⁻¹]        Elektron-Neutral-Kollisionsfrequenz: Kel(Te) * ng
    double u_Bohm;       // [m/s]        Bohm-Geschwindigkeit: sqrt(kB * Te / M)
    double pf;           // [rad/s]      Plasmafrequenz: sqrt(n * e² / (me * eps0))
    double eps_p_real;   // [-]          Re(epsilon_p) = 1 - omega_pe² / (omega² + nu_m²)
    double eps_p_imag;   // [-]          Im(epsilon_p) = -omega_pe² * nu_m / (omega * (omega² + nu_m²))

    // --- RF-Kopplung ---
    double zeta;         // [-]          ICP-Kopplungseffizienz: R_ind / (R_ind + R_coil)  [Eq.22]
    double icp_eff;      // [-]          = zeta (Alias fuer Ausgabe-Kompatibilitaet)
    double P_RF;         // [W]          Gesamt-RF-Leistung: 0.5 * (R_ind + R_coil) * I_coil²  [Eq.26]

    // --- Schub ---
    double T_i_N;        // [N]          Ionenschub: Gamma_i * M * v_beam * Ai  [Eq.17]
    double T_n_N;        // [N]          Neutralteilchenschub: Gamma_g * M * vg * Ag  [Eq.20]
    double T_total_N;    // [N]          Gesamtschub: T_i + T_n

    // --- Effizienzen ---
    double gamma_eff;    // [-]          Schubleistungseffizienz: (P_i + P_n) / (P_i + P_n + P_RF)  [Eq.23]
    double xi_mN_kW;     // [mN/kW]      Schubeffizienz: 1000 * T_total / P_RF  [Eq.24]
    double eta_mass;     // [-]          Massennutzungsgrad: Gamma_i * Ai / Q0  [Eq.25]
};

// Berechne alle abgeleiteten Groessen aus dem konvergierten Zustand.
// Input:  n, ng [m⁻³], Te [eV], Tg [K], R_ind [Ohm], I_coil [A], P_abs [W]
// Output: DerivedQuantities struct (siehe Felddokumentation oben)
//
// HINWEIS: Te wird intern ueber conv (=11604.5 K/eV) nach Kelvin umgerechnet,
// wo noetig (uB, Energieterme). In der Ausgabe bleibt Te in eV.
static DerivedQuantities compute_derived(double n, double ng, double Te, double Tg,
                                         double R_ind, double I_coil_val, double P_abs_val) {
    DerivedQuantities d{};

    // Ionenfluss: Gamma_i = hL * n * uB  [m⁻² s⁻¹]
    double lam  = lambda_i(ng);          // mittlere freie Weglaenge [m]
    double Gi_f = Gamma_i_func(lam, Te, n);

    // Strahlstrom: I = Ai * e * Gamma_i  [A], ausgegeben in [mA]
    d.I_extr_mA   = Ai * e * Gi_f * 1000.0;

    // Ionisierungsgrad: n / ng * 100  [%]
    d.iondeg      = n / std::max(ng, 1.0) * 100.0;

    // Kollisionsfrequenz: nu_m = Kel(Te) * ng  [s⁻¹]
    d.cf          = coll_freq(ng, Te);

    // Bohm-Geschwindigkeit: uB = sqrt(kB * Te / M)  [m/s]
    d.u_Bohm      = uB(Te);

    // Ionenstromdichte: J_i = e * Gamma_i  [A/m²]
    d.J_i         = e * Gi_f;

    // Plasmafrequenz: omega_pe = sqrt(n * e² / (me * eps0))  [rad/s]
    d.pf          = plasma_freq(n);

    // Komplexe Plasmapermittivitaet: eps_p = 1 - omega_pe² / (omega * (omega - i*nu_m))
    complex<double> ep = my_calc_eps_p(n, ng, Te);
    d.eps_p_real  = real(ep);
    d.eps_p_imag  = imag(ep);

    // ICP-Kopplungseffizienz: zeta = R_ind / (R_ind + R_coil)  [Eq.22]
    d.zeta        = (R_ind > 1e-10) ? R_ind / (R_ind + R_ohm) : 0.0;
    d.icp_eff     = d.zeta;

    // Neutralteilchenfluss: Gamma_g = 0.25 * ng * vg  [m⁻² s⁻¹]
    double vg_f   = vg(Tg);
    double Gn_f   = 0.25 * ng * vg_f;

    // Ionenschub: T_i = Gamma_i * M * v_beam * Ai  [N]  [Eq.17]
    d.T_i_N       = Thrust(Gi_f);

    // Neutralteilchenschub: T_n = Gamma_g * M * vg * Ag  [N]  [Eq.20]
    d.T_n_N       = Gn_f * M * vg_f * Ag;
    d.T_total_N   = d.T_i_N + d.T_n_N;

    // Gesamt-RF-Leistung: P_RF = 0.5 * (R_ind + R_coil) * I_coil²  [W]  [Eq.26]
    d.P_RF        = 0.5 * (R_ind + R_ohm) * I_coil_val * I_coil_val;

    // Schubleistungseffizienz: gamma = (P_i + P_n) / (P_i + P_n + P_RF)  [Eq.23]
    double v_extr = sqrt(2.0 * e * Vgrid / M);                // Extraktionsgeschw. [m/s]
    double pow_i  = 0.5 * M * v_extr * v_extr * Gi_f * Ai;   // Ionenschubleistung [W]  [Eq.18]
    double pow_n  = 0.5 * M * vg_f * vg_f * Gn_f * Ag;       // Neutralschubleistung [W] [Eq.21]
    d.gamma_eff   = (d.P_RF > 1e-10) ? (pow_i + pow_n) / (pow_i + pow_n + d.P_RF) : 0.0;

    // Schubeffizienz: xi = (T_i + T_n) / P_RF  [mN/kW]  [Eq.24]
    d.xi_mN_kW    = (d.P_RF > 1e-10) ? 1000.0 * d.T_total_N / d.P_RF : 0.0;

    // Massennutzungsgrad: eta = Gamma_i * Ai / Q0  [-]  [Eq.25]
    d.eta_mass    = (Q0 > 0.0) ? Gi_f * Ai / Q0 : 0.0;

    return d;
}

// Schreibt eine CSV-Datenzeile in die Ausgabedatei (fuer alle Solver-Pfade)
static void emit_csv_row(std::ofstream& datei, const std::string& method_str,
                         double Q0sccm_val, double n, double ng, double Te, double Tg,
                         double P_RFG_val, double P_abs_val,
                         double R_ind_val, double I_coil_val,
                         const DerivedQuantities& d) {
    datei << method_str  << ", " << Q0sccm_val << ", " << Te << ", " << Tg << ", "
          << std::scientific << n << ", " << ng << ", "
          << std::fixed << d.iondeg << ", " << P_RFG_val << ", "
          << P_abs_val << ", " << d.I_extr_mA << ", "
          << d.cf << ", " << R_ind_val << ", " << I_coil_val << ", "
          << d.eps_p_real << ", " << d.eps_p_imag << ", "
          << d.u_Bohm << ", " << d.J_i << ", "
          << d.zeta << ", " << d.gamma_eff << ", "
          << d.xi_mN_kW << ", " << d.eta_mass << ", "
          << d.pf << ", " << frequency/1e6 << ", "
          << d.T_i_N*1e3 << ", " << d.T_n_N*1e3 << ", " << d.T_total_N*1e3 << ", "
          << d.icp_eff << ", " << d.P_RF << ", "
          << density_profile_factor * n << ", " << density_profile_factor << "\n";
}

static std::array<double,4> stationary_residual_raw(const StationaryPlasmaState& s,
                                                    double P_RFG_local,
                                                    StationaryRFState* rf_out = nullptr) {
    StationaryRFState rf = stationary_compute_rf(s.n, s.ng, s.Te, P_RFG_local);
    if (rf_out) *rf_out = rf;
    if (!rf.valid) {
        double nan = std::numeric_limits<double>::quiet_NaN();
        return {nan,nan,nan,nan};
    }

    double P_vol = rf.P_abs * P_abs_scale / V;  // P_abs_scale: diagnostischer Faktor

    // n_eff: effektive mittlere Dichte fuer Volumenverlustterme
    // Bei density_profile_factor = 1.0 ist n_eff = n (unveraendertes 0D-Verhalten).
    double n_eff = density_profile_factor * s.n;

    // r1: Ionenbilanz — Produktion (Volumenterm) verwendet n_eff,
    //     Wandverlust (Bohm-Fluss) verwendet n (Randdichte bleibt unveraendert)
    double r1 = n_eff*s.ng*Kiz(s.Te) - s.n*uB(s.Te)*Aeff(lambda_i(s.ng))/V;

    // r2: Neutralgasbilanz — unveraendert (kein Dichtekorrektur hier)
    double r2 = Q0/V + s.n*uB(s.Te)*Aeff1(lambda_i(s.ng))/V - s.n*s.ng*Kiz(s.Te) - Gamma_g(s.ng,vg(s.Tg))*Ag/V;

    // r4: Gasenergiebilanz — unveraendert
    double Pg1 = 3.0*me/M * kB*(s.Te*conv-s.Tg) * s.n*s.ng*Kel(s.Te);
    double Pg2 = 0.25*M*uB(s.Te)*uB(s.Te) * s.n*s.ng*sigma_i*vi(s.Tg);
    double Pg3 = kappa*(s.Tg-Tg0)/lambda_0 * A/V;
    double r4 = Pg1 + Pg2 - Pg3;

    // r3: Elektronenenergiebilanz — Volumenverlustterme verwenden n_eff,
    //     Wandverlustterm P5 verwendet n (Bohm-Fluss = Randdichte)
    double P2 = Eiz  * n_eff*s.ng*Kiz(s.Te);
    // P3: Anregungsverluste — Legacy (Chabert) oder tabelliert (Biagi)
    double P3;
    if (excitation_model == 1 && !g_kex_table.empty()) {
        // Tabelliert: Pexc_coeff = Summe_i [Kex_i * Eexc_i] (prozessaufgeloest)
        P3 = interpolate_kex_table(s.Te, true) * n_eff * s.ng;
    } else {
        // Legacy: Eexc * Kex(Te) mit globaler Schwellenenergie
        P3 = Eexc * n_eff*s.ng*Kex(s.Te);
    }
    double P4 = 3.0*me/M * kB*(s.Te*conv-s.Tg) * n_eff*s.ng*Kel(s.Te);
    double P5 = alpha_e_wall*kB*s.Te*conv * s.n*uB(s.Te)*Aeff(lambda_i(s.ng))/V;
    double r3 = P_vol - (P2 + P3 + P4 + P5);

    return {r1, r2, r3, r4};
}

static std::array<double,4> stationary_residual_scaled(const StationaryPlasmaState& s,
                                                       double P_RFG_local,
                                                       StationaryRFState* rf_out = nullptr) {
    if (!stationary_state_in_bounds(s)) {
        double nan = std::numeric_limits<double>::quiet_NaN();
        return {nan,nan,nan,nan};
    }
    auto raw = stationary_residual_raw(s, P_RFG_local, rf_out);
    if (!std::isfinite(raw[0]) || !std::isfinite(raw[1]) || !std::isfinite(raw[2]) || !std::isfinite(raw[3])) {
        double nan = std::numeric_limits<double>::quiet_NaN();
        return {nan,nan,nan,nan};
    }
    StationaryRFState rf = rf_out ? *rf_out : stationary_compute_rf(s.n, s.ng, s.Te, P_RFG_local);
    if (!rf.valid || !std::isfinite(rf.P_abs) || rf.P_abs < 0.0) {
        double nan = std::numeric_limits<double>::quiet_NaN();
        return {nan,nan,nan,nan};
    }
    // Physikalisch motivierte Skalierung: jedes Residuum wird durch die
    // Groessenordnung seines dominanten Terms geteilt, damit alle vier
    // Gleichungen im Solver gleich stark wiegen.
    double scale1 = std::max(1e-20, std::fabs(s.n * s.ng * Kiz(s.Te)));
    double scale2 = std::max(1e-20, Q0 / V);
    double scale3 = std::max(1e-6,  rf.P_abs / V);
    double scale4 = std::max(1e-6,  std::fabs(kappa * (s.Tg - Tg0) / lambda_0 * A / V)
                           + std::fabs(3.0 * me / M * kB * (s.Te * conv - s.Tg) * s.n * s.ng * Kel(s.Te)));
    return {raw[0]/scale1, raw[1]/scale2, raw[2]/scale3, raw[3]/scale4};
}

static bool stationary_solve_linear_4x4(double A_[4][4], double b[4], double x[4]) {
    double M_[4][5];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) M_[i][j] = A_[i][j];
        M_[i][4] = b[i];
    }

    for (int col = 0; col < 4; ++col) {
        int piv = col;
        for (int r = col+1; r < 4; ++r)
            if (fabs(M_[r][col]) > fabs(M_[piv][col])) piv = r;
        if (fabs(M_[piv][col]) < 1e-20 || !std::isfinite(M_[piv][col])) return false;
        if (piv != col) for (int c = col; c < 5; ++c) std::swap(M_[piv][c], M_[col][c]);

        double div = M_[col][col];
        for (int c = col; c < 5; ++c) M_[col][c] /= div;
        for (int r = 0; r < 4; ++r) {
            if (r == col) continue;
            double fac = M_[r][col];
            for (int c = col; c < 5; ++c) M_[r][c] -= fac * M_[col][c];
        }
    }

    for (int i = 0; i < 4; ++i) x[i] = M_[i][4];
    return true;
}

static double stationary_merit(const std::array<double,4>& r, const StationaryPlasmaState& s) {
    if (!stationary_state_in_bounds(s)) return std::numeric_limits<double>::infinity();
    double m = 0.0;
    for (double v : r) {
        if (!std::isfinite(v)) return std::numeric_limits<double>::infinity();
        m = std::max(m, std::fabs(v));
    }
    return m;
}

struct StationarySolveResult {
    bool converged = false;
    bool soft_ok = false;   // Soft-Accept: brauchbarer Zustand, aber Energiebilanz nicht erfuellt
    StationaryPlasmaState state{};
    StationaryRFState rf{};
    int iterations = 0;
    double resid_norm = std::numeric_limits<double>::infinity();
    std::string reason;
};

static bool stationary_soft_accept(const StationarySolveResult& r, double initial_merit) {
    if (!stationary_state_finite_positive(r.state)) return false;
    if (!stationary_state_in_bounds(r.state)) return false;
    if (!r.rf.valid) return false;
    if (!std::isfinite(r.resid_norm)) return false;
    if (!std::isfinite(initial_merit) || initial_merit <= 0.0) return false;

    double rel_target = stationary_soft_rel_improve * initial_merit;
    double thresh = std::min(stationary_soft_abs_resid, rel_target);
    bool enough_improvement = r.resid_norm <= rel_target;
    bool low_enough_abs = r.resid_norm <= stationary_soft_abs_resid;
    bool enough_iterations = r.iterations >= 2;
    bool accept = enough_improvement && low_enough_abs && enough_iterations;

    std::ostringstream os;
    os << "resid=" << std::scientific << r.resid_norm
       << " initial=" << initial_merit
       << " abs_thr=" << stationary_soft_abs_resid
       << " rel_thr=" << rel_target
       << " iter=" << r.iterations
       << " accept=" << accept;
    debug_emit(3, accept ? "SOFT_ACCEPT_OK" : "SOFT_ACCEPT_REJECT", os.str());
    return accept;
}


// ============================================================
// Levenberg-Marquardt Solver fuer das 4D Plasmagleichgewicht
//
// Ersetzt Newton + Line-Search. Gruende:
// 1. Newton scheitert bei schlecht konditioniertem Jacobian
//    (n kuerzt sich aus der Ionisierungsbilanz → J-Spalte ≈ 0)
// 2. Diskrete Line-Search (7 alpha-Werte) findet bei flachen
//    Merit-Landschaften keinen Abstieg
// 3. LM interpoliert automatisch zwischen Steepest Descent
//    (grosses λ, robust) und Newton (kleines λ, schnell)
//
// Algorithmus:
//   Loesung von (J^T J + λ D) dx = -J^T F
//   D = diag(J^T J) fuer Skaleninvarianz (Marquardt-Variante)
//   λ-Update ueber Gain Ratio ρ = actual / predicted reduction
//   Akzeptanz: jeder Schritt mit Kostenreduktion
// ============================================================

static StationarySolveResult stationary_solve_lm(double P_RFG_local,
                                                  const StationaryPlasmaState& initial) {
    StationarySolveResult out;
    out.state = initial;

    if (!stationary_state_in_bounds(initial)) {
        out.reason = "invalid initial state";
        return out;
    }

    std::array<double,4> x = {log(initial.n), log(initial.ng), log(initial.Te), log(initial.Tg)};

    StationaryPlasmaState s{exp(x[0]), exp(x[1]), exp(x[2]), exp(x[3])};
    StationaryRFState rf;
    auto F = stationary_residual_scaled(s, P_RFG_local, &rf);

    if (!rf.valid || !std::isfinite(F[0]) || !std::isfinite(F[1]) ||
        !std::isfinite(F[2]) || !std::isfinite(F[3])) {
        out.reason = "invalid initial residual";
        out.resid_norm = std::numeric_limits<double>::infinity();
        return out;
    }

    double cost = 0.0;
    for (double v : F) cost += v * v;
    cost *= 0.5;

    // LM-Steuerparameter
    double lambda = 1e-2;
    const double lambda_min = 1e-12;
    const double lambda_max = 1e8;
    const int max_iter = 80;
    const int max_lm_tries = 15;
    int stagnation_count = 0;
    const int max_stagnation = 8;
    double prev_cost = cost;

    for (int iter = 0; iter < max_iter; ++iter) {
        double merit = stationary_merit(F, s);

        // Konvergenzpruefung: merit < newton_tol bedeutet alle 4 skalierten
        // Residuen (inkl. r3/Energiebilanz) sind unter der Schwelle.
        if (merit < newton_tol) {
            out.converged = true;
            out.state = s;
            out.rf = rf;
            out.iterations = iter;
            out.resid_norm = merit;
            out.reason = "ok";
            debug_emit(2, "LM_OK", "iter=" + std::to_string(iter) +
                " merit=" + std::to_string(merit) + " lambda=" + std::to_string(lambda));
            return out;
        }

        // Fortschritt (kompatibel mit GUI-Logging)
        std::cout << "NEWTON_IT " << iter << " "
                  << std::fixed << std::setprecision(4) << P_RFG_local << " "
                  << std::scientific << std::setprecision(6)
                  << merit << " " << s.n << " " << s.ng << " "
                  << std::fixed << std::setprecision(6)
                  << s.Te << " " << s.Tg << std::endl;

        if (!std::isfinite(merit)) {
            out.reason = "invalid residual";
            out.state = s; out.rf = rf;
            out.iterations = iter; out.resid_norm = merit;
            return out;
        }

        // --- Jacobian (zentrale Finite Differenzen) ---
        double J[4][4];
        bool jac_ok = true;
        for (int j = 0; j < 4 && jac_ok; ++j) {
            std::array<double,4> xp = x, xm = x;
            double h = newton_fd_eps * std::max(1.0, std::fabs(x[j]));
            xp[j] += h;
            xm[j] -= h;
            StationaryPlasmaState sp{exp(xp[0]), exp(xp[1]), exp(xp[2]), exp(xp[3])};
            StationaryPlasmaState sm{exp(xm[0]), exp(xm[1]), exp(xm[2]), exp(xm[3])};
            auto rp = stationary_residual_scaled(sp, P_RFG_local, nullptr);
            auto rm = stationary_residual_scaled(sm, P_RFG_local, nullptr);
            for (int i = 0; i < 4; ++i) {
                if (!std::isfinite(rp[i]) || !std::isfinite(rm[i])) { jac_ok = false; break; }
                J[i][j] = (rp[i] - rm[i]) / (2.0 * h);
            }
        }
        if (!jac_ok) {
            out.reason = "jacobian nan";
            out.state = s; out.rf = rf;
            out.iterations = iter; out.resid_norm = merit;
            return out;
        }

        // --- J^T J und J^T F aufbauen ---
        double JtJ[4][4] = {};
        double JtF[4] = {};
        double diag_max = 0.0;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                double sum = 0.0;
                for (int k = 0; k < 4; ++k) sum += J[k][i] * J[k][j];
                JtJ[i][j] = sum;
            }
            double sf = 0.0;
            for (int k = 0; k < 4; ++k) sf += J[k][i] * F[k];
            JtF[i] = sf;
            diag_max = std::max(diag_max, JtJ[i][i]);
        }

        // Gradient-Norm-Check: lokales Minimum von ||F||² aber F≠0
        double grad_norm_sq = 0.0;
        for (int i = 0; i < 4; ++i) grad_norm_sq += JtF[i] * JtF[i];
        if (grad_norm_sq < 1e-28 * std::max(1.0, cost * cost)) {
            out.reason = "local minimum (grad~0)";
            out.state = s; out.rf = rf;
            out.iterations = iter; out.resid_norm = merit;
            DBG(1, "LM_LOCAL_MIN iter=" << iter << " merit=" << merit
                << " grad=" << sqrt(grad_norm_sq));
            return out;
        }

        // --- LM-Schritt: (J^T J + λ D) dx = -J^T F ---
        // Versuche mit steigendem λ bis Kostenreduktion erreicht
        bool step_accepted = false;
        for (int lm_try = 0; lm_try < max_lm_tries; ++lm_try) {
            double A[4][4];
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) A[i][j] = JtJ[i][j];
                // Marquardt-Daempfung: skaliert mit Diagonale fuer Invarianz
                A[i][i] += lambda * std::max(JtJ[i][i], 1e-10 * std::max(diag_max, 1e-20));
            }

            double neg_JtF[4] = {-JtF[0], -JtF[1], -JtF[2], -JtF[3]};
            double dx[4];
            if (!stationary_solve_linear_4x4(A, neg_JtF, dx)) {
                lambda = std::min(lambda * 4.0, lambda_max);
                if (lambda >= lambda_max) break;
                continue;
            }

            // Schritt im Log-Raum begrenzen
            for (int i = 0; i < 4; ++i)
                dx[i] = std::max(-newton_max_log_step, std::min(newton_max_log_step, dx[i]));

            // Trial-Punkt auswerten
            std::array<double,4> x_trial;
            for (int i = 0; i < 4; ++i) x_trial[i] = x[i] + dx[i];

            StationaryPlasmaState st{exp(x_trial[0]), exp(x_trial[1]),
                                     exp(x_trial[2]), exp(x_trial[3])};
            if (!stationary_state_in_bounds(st)) {
                lambda = std::min(lambda * 4.0, lambda_max);
                if (lambda >= lambda_max) break;
                continue;
            }

            StationaryRFState rf_trial;
            auto F_trial = stationary_residual_scaled(st, P_RFG_local, &rf_trial);

            double cost_trial = 0.0;
            bool F_ok = rf_trial.valid;
            for (double v : F_trial) {
                if (!std::isfinite(v)) { F_ok = false; break; }
                cost_trial += v * v;
            }
            cost_trial *= 0.5;

            if (!F_ok || !std::isfinite(cost_trial)) {
                lambda = std::min(lambda * 4.0, lambda_max);
                if (lambda >= lambda_max) break;
                continue;
            }

            // Gain Ratio: ρ = tatsaechliche / vorhergesagte Reduktion
            // pred = 0.5 dx^T J^T J dx + λ dx^T D dx  (immer >= 0)
            double actual = cost - cost_trial;
            double pred;
            {
                double jtj_part = 0.0, damp_part = 0.0;
                for (int i = 0; i < 4; ++i) {
                    double jj_dx = 0.0;
                    for (int j2 = 0; j2 < 4; ++j2) jj_dx += JtJ[i][j2] * dx[j2];
                    jtj_part += dx[i] * jj_dx;
                    damp_part += lambda * std::max(JtJ[i][i], 1e-10 * std::max(diag_max, 1e-20))
                                 * dx[i] * dx[i];
                }
                pred = 0.5 * jtj_part + damp_part;
                if (pred <= 0.0) pred = 1e-30;
            }

            double rho = actual / pred;

            DBG(3, "LM_TRY iter=" << iter << " try=" << lm_try
                << " lam=" << lambda << " cost=" << cost
                << " trial=" << cost_trial << " rho=" << rho);

            if (actual > 0.0) {
                // Schritt akzeptiert: Kosten gesunken
                x = x_trial;
                s = st;
                rf = rf_trial;
                F = F_trial;
                cost = cost_trial;
                step_accepted = true;

                // λ-Update basierend auf Gain Ratio
                if (rho > 0.75)      lambda = std::max(lambda / 3.0, lambda_min);
                else if (rho < 0.25) lambda = std::min(lambda * 2.0, lambda_max);
                break;
            } else {
                lambda = std::min(lambda * 4.0, lambda_max);
                if (lambda >= lambda_max) break;
            }
        }

        if (!step_accepted) {
            out.reason = "lm step rejected";
            out.state = s; out.rf = rf;
            out.iterations = iter; out.resid_norm = merit;
            DBG(1, "LM_REJECT iter=" << iter << " merit=" << merit << " lam=" << lambda);
            return out;
        }

        // Stagnations-Erkennung
        double rel_reduction = (prev_cost - cost) / std::max(prev_cost, 1e-30);
        if (rel_reduction < 1e-10) {
            stagnation_count++;
            if (stagnation_count >= max_stagnation) {
                double final_merit = stationary_merit(F, s);
                out.reason = "stagnation";
                out.state = s; out.rf = rf;
                out.iterations = iter; out.resid_norm = final_merit;
                DBG(1, "LM_STAGNATION iter=" << iter << " merit=" << final_merit);
                return out;
            }
        } else {
            stagnation_count = 0;
        }
        prev_cost = cost;
    }

    // max_iter erreicht
    double final_merit = stationary_merit(F, s);
    out.converged = final_merit < newton_tol;
    out.state = s;
    out.rf = rf;
    out.iterations = max_iter;
    out.resid_norm = final_merit;
    out.reason = out.converged ? "ok" : "lm max iter";
    return out;
}

static double stationary_beam_current_mA(const StationaryPlasmaState& s) {
    double lam = lambda_i(s.ng);
    return Ai * e * Gamma_i_func(lam, s.Te, s.n) * 1000.0;
}


static StationarySolveResult stationary_solve_ptc_then_newton(double P_RFG_local,
                                                              const StationaryPlasmaState& initial) {
    StationarySolveResult best;
    best.state = initial;

    if (!stationary_state_in_bounds(initial)) {
        best.reason = "invalid initial state";
        return best;
    }

    std::array<double,4> x0 = {log(initial.n), log(initial.ng), log(initial.Te), log(initial.Tg)};
    StationaryPlasmaState s0 = initial;
    StationaryRFState rf0;
    auto r0 = stationary_residual_scaled(s0, P_RFG_local, &rf0);
    double merit_initial = stationary_merit(r0, s0);

    // 1) Erst direkter LM-Versuch (schneller Pfad)
    best = stationary_solve_lm(P_RFG_local, initial);
    if (best.converged) return best;
    if (stationary_soft_accept(best, merit_initial)) {
        best.converged = false;
        best.soft_ok = true;
        best.reason = "soft-ok-direct";
        return best;
    }

    // 2) Pseudo-transient continuation in Log-Variablen
    std::array<double,4> x = x0;
    StationaryPlasmaState s = initial;
    StationaryRFState rf = rf0;  // BUG-FIX: war vorher uninitialisiert → PTC lief nie
    auto r = r0;
    double merit = merit_initial;

    if (!std::isfinite(merit) || !rf.valid) {
        return best;
    }

    double gain = ptc_start_gain;
    StationarySolveResult best_local = best;
    if (std::isfinite(merit) && merit < best_local.resid_norm) {
        best_local.state = s;
        best_local.rf = rf;
        best_local.resid_norm = merit;
        best_local.iterations = 0;
        best_local.reason = "ptc-start";
    }

    for (int it = 0; it < ptc_max_iter; ++it) {
        bool accepted = false;
        std::array<double,4> x_trial = x;

        for (double alpha : {1.0, 0.5, 0.25, 0.125, 0.0625}) {
            for (int j = 0; j < 4; ++j) {
                double step = -alpha * gain * r[j];
                step = std::max(-0.12, std::min(0.12, step));
                x_trial[j] = x[j] + step;
            }

            StationaryPlasmaState st{exp(x_trial[0]), exp(x_trial[1]), exp(x_trial[2]), exp(x_trial[3])};
            if (!stationary_state_in_bounds(st)) continue;

            StationaryRFState rf_trial;
            auto r_trial = stationary_residual_scaled(st, P_RFG_local, &rf_trial);
            double merit_trial = stationary_merit(r_trial, st);

            if (!std::isfinite(merit_trial) || !rf_trial.valid) continue;

            if (merit_trial < ptc_accept_ratio * merit) {
                x = x_trial;
                s = st;
                r = r_trial;
                rf = rf_trial;
                merit = merit_trial;
                accepted = true;

                if (merit < best_local.resid_norm) {
                    best_local.state = s;
                    best_local.rf = rf;
                    best_local.resid_norm = merit;
                    best_local.iterations = it + 1;
                    best_local.reason = "ptc-progress";
                }
                break;
            }
        }

        if (!accepted) {
            gain *= 0.5;
            if (gain < ptc_min_gain) break;
            continue;
        }

        if (merit < ptc_switch_merit) {
            StationarySolveResult polished = stationary_solve_lm(P_RFG_local, s);
            if (polished.converged) return polished;
            if (stationary_soft_accept(polished, merit_initial)) {
                polished.converged = false;
                polished.soft_ok = true;
                polished.reason = "soft-ok-polished";
                return polished;
            }
            if (polished.resid_norm < best_local.resid_norm) best_local = polished;
        }
    }

    // 3) Letzter Newton-Finish vom besten gefundenen Zustand
    StationarySolveResult final_try = stationary_solve_lm(P_RFG_local, best_local.state);
    if (final_try.converged) return final_try;
    if (stationary_soft_accept(final_try, merit_initial)) {
        final_try.converged = false;
        final_try.soft_ok = true;
        final_try.reason = "soft-ok-final";
        return final_try;
    }
    if (stationary_soft_accept(best_local, merit_initial)) {
        best_local.converged = false;
        best_local.soft_ok = true;
        if (best_local.reason.empty() || best_local.reason == "all starts failed")
            best_local.reason = "soft-ok-best";
        return best_local;
    }
    if (final_try.resid_norm < best_local.resid_norm) return final_try;
    return best_local;
}

static StationarySolveResult stationary_solve_power_robust(double P_RFG_local,
                                                           const std::vector<StationaryPlasmaState>& starts) {
    StationarySolveResult best;
    best.reason = "all starts failed";
    int start_idx = 0;
    for (const auto& st : starts) {
        if (!stationary_state_in_bounds(st)) {
            debug_emit(3, "START_SKIP", "P_RFG=" + std::to_string(P_RFG_local) + " idx=" + std::to_string(start_idx) + " reason=out_of_bounds");
            ++start_idx;
            continue;
        }
        debug_emit(3, "START_TRY", "P_RFG=" + std::to_string(P_RFG_local) + " idx=" + std::to_string(start_idx) + " " + state_summary(st.n, st.ng, st.Te, st.Tg, 0.0, 0.0));
        StationarySolveResult cur = stationary_solve_ptc_then_newton(P_RFG_local, st);
        std::ostringstream os;
        os << "P_RFG=" << std::fixed << std::setprecision(6) << P_RFG_local
           << " idx=" << start_idx
           << " converged=" << cur.converged
           << " resid=" << std::scientific << cur.resid_norm
           << " I_mA=" << std::fixed << std::setprecision(6) << (stationary_state_finite_positive(cur.state) ? stationary_beam_current_mA(cur.state) : std::numeric_limits<double>::quiet_NaN())
           << " reason=" << cur.reason;
        debug_emit(cur.converged ? 2 : 3, cur.converged ? "START_OK" : "START_FAIL", os.str());
        if (cur.converged) return cur;
        if (cur.resid_norm < best.resid_norm) best = cur;
        ++start_idx;
    }
    return best;
}

enum class SolveFailType {
    NONE = 0,               // Kein Fehler (konvergiert)
    NO_PHYSICAL_SOLUTION,   // Physikalisch keine Loesung im P_RFG-Bereich
    NUMERICAL_FAIL          // Numerisches Versagen (Loesung koennte existieren)
};

struct StationaryPowerSolveResult {
    bool converged = false;
    bool hit_limit = false;
    SolveFailType fail_type = SolveFailType::NONE;
    StationaryPlasmaState state{};
    StationaryRFState rf{};
    double P_RFG_sol = std::numeric_limits<double>::quiet_NaN();
    double P_trial_last = std::numeric_limits<double>::quiet_NaN();
    double I_mA = std::numeric_limits<double>::quiet_NaN();
    double err_mA = std::numeric_limits<double>::quiet_NaN();
    double inner_resid_norm = std::numeric_limits<double>::infinity();
    int iterations = 0;
    std::string reason;
};

static bool stationary_power_result_valid(const StationaryPowerSolveResult& ps) {
    return ps.converged &&
           stationary_state_finite_positive(ps.state) &&
           stationary_state_in_bounds(ps.state) &&
           ps.rf.valid &&
           std::isfinite(ps.P_RFG_sol) &&
           ps.P_RFG_sol > 0.0 &&
           std::isfinite(ps.I_mA);
}

// ============================================================
// Selbstkonsistenter Modus (solve_mode == 2):
// P_RFG ist fest vorgegeben. Der 4D-LM-Solver loest das
// gekoppelte Gleichungssystem (r1-r4) bei diesem P_RFG.
// Der Strahlstrom I_mA ergibt sich als Ausgabe, nicht als Ziel.
// ============================================================

static StationaryPowerSolveResult stationary_solve_at_fixed_power(
    double P_RFG_fixed, const StationaryPlasmaState& guess) {

    StationaryPowerSolveResult out;

    debug_emit(2, "SELFCONSISTENT_BEGIN",
        "Q0sccm=" + std::to_string(Q0 / SCCM_TO_PPS) +
        " P_RFG=" + std::to_string(P_RFG_fixed));

    // Versuche mehrere Startwerte
    StationaryPlasmaState safe = stationary_safe_defaults_for_q(Q0);
    std::vector<StationaryPlasmaState> starts = {guess, safe};
    if (stationary_state_in_bounds(guess)) {
        starts.push_back(StationaryPlasmaState{
            std::sqrt(std::max(1.0, guess.n * safe.n)),
            std::sqrt(std::max(1.0, guess.ng * safe.ng)),
            0.5 * (guess.Te + safe.Te),
            0.5 * (guess.Tg + safe.Tg)
        });
    }

    StationarySolveResult best;
    best.reason = "all starts failed";
    for (const auto& st : starts) {
        if (!stationary_state_in_bounds(st)) continue;
        // PTC + LM: robuster als direkter LM-Aufruf bei schwierigen Startwerten
        StationarySolveResult cur = stationary_solve_ptc_then_newton(P_RFG_fixed, st);
        if (cur.converged) { best = cur; break; }
        if (cur.resid_norm < best.resid_norm) best = cur;
    }

    out.iterations = best.iterations;
    out.inner_resid_norm = best.resid_norm;
    out.P_RFG_sol = P_RFG_fixed;
    out.P_trial_last = P_RFG_fixed;

    if (best.converged && stationary_state_finite_positive(best.state)) {
        out.converged = true;
        out.state = best.state;
        out.rf = best.rf;
        out.I_mA = stationary_beam_current_mA(best.state);
        out.err_mA = 0.0;  // Kein Zielstrom im selbstkonsistenten Modus
        out.reason = "ok";

        std::cout << "PID_DONE " << std::fixed << std::setprecision(4)
                  << out.I_mA << " 0.0000 " << P_RFG_fixed << " "
                  << best.state.Te << " " << best.state.Tg << std::endl;
        std::cout << "CONVERGED " << best.iterations << std::endl;

        debug_emit(2, "SELFCONSISTENT_OK",
            "P_RFG=" + std::to_string(P_RFG_fixed) +
            " I_mA=" + std::to_string(out.I_mA) +
            " resid=" + std::to_string(best.resid_norm));
    } else {
        out.converged = false;
        out.fail_type = SolveFailType::NUMERICAL_FAIL;
        out.state = best.state;
        out.rf = best.rf;
        out.I_mA = stationary_state_finite_positive(best.state)
                   ? stationary_beam_current_mA(best.state)
                   : std::numeric_limits<double>::quiet_NaN();
        out.err_mA = std::numeric_limits<double>::quiet_NaN();
        out.reason = "self-consistent solve failed: " + best.reason;
        debug_emit(1, "SELFCONSISTENT_FAIL",
            "P_RFG=" + std::to_string(P_RFG_fixed) +
            " reason=" + best.reason, true);
    }
    return out;
}

// ============================================================
// Modus 1: Fester Strahlstrom — 4D-Solver mit Power-Bisection
//
// Der LM-Solver loest alle 4 Gleichungen fuer gegebenes P_RFG.
// Die aeussere Bisection variiert P_RFG bis I_mA = I_soll.
// ============================================================

static StationaryPowerSolveResult stationary_solve_for_target_current(const StationaryPlasmaState& guess) {
    StationaryPowerSolveResult out;

    debug_emit(2, "TARGET_CURRENT_BEGIN",
        "Q0sccm=" + std::to_string(Q0 / SCCM_TO_PPS) +
        " I_soll=" + std::to_string(I_soll) +
        " P_RFG_seed=" + std::to_string(P_RFG));

    // --- Innerer 4D-Solve fuer ein gegebenes P_RFG ---
    // Gibt konvergierten Zustand + I_mA zurueck, oder Fail.
    auto solve_at_power = [&](double p_try, const StationaryPlasmaState& start) -> StationarySolveResult {
        // Versuche mehrere Startwerte
        StationaryPlasmaState safe = stationary_safe_defaults_for_q(Q0);
        std::vector<StationaryPlasmaState> starts = {start, safe};
        // Geometrisches Mittel als dritter Startwert
        if (stationary_state_in_bounds(start)) {
            starts.push_back(StationaryPlasmaState{
                std::sqrt(std::max(1.0, start.n * safe.n)),
                std::sqrt(std::max(1.0, start.ng * safe.ng)),
                0.5 * (start.Te + safe.Te),
                0.5 * (start.Tg + safe.Tg)
            });
        }

        StationarySolveResult best;
        best.reason = "all starts failed";
        for (const auto& st : starts) {
            if (!stationary_state_in_bounds(st)) continue;
            StationarySolveResult cur = stationary_solve_lm(p_try, st);
            if (cur.converged) return cur;
            if (cur.resid_norm < best.resid_norm) best = cur;
        }
        return best;
    };

    // --- Bracket-Suche: I_mA(P_lo) < I_soll < I_mA(P_hi) ---
    //
    // Strategie: Von unten aufbauen statt von oben herab.
    // Grund: Bei tabellierten Ratenmodellen (v.a. full_tabulated) hat der
    // innere Solver eine deutlich niedrigere Konvergenzgrenze in P_RFG
    // (~60-70W bei Xe), weil der hoehere Kel zu steilerer Ionisierung fuehrt
    // und der Zustandsraum bei hoher Leistung kollabiiert (ng->0, iondeg->max).
    // Wenn p_hi direkt auf 100W gesetzt wird und dort der Solver scheitert,
    // findet der alte Algorithmus nie ein Bracket und expandiert sinnlos.
    //
    // Neuer Ansatz:
    // 1. Starte mit p_lo=2W, evaluiere
    // 2. Verdopple p_hi schrittweise (4, 8, 16, 32, 64, ...) bis
    //    entweder ein Bracket gefunden wird (f_lo*f_hi < 0)
    //    oder p_hi nicht mehr konvergiert
    // 3. Wenn p_hi nicht konvergiert: halbiere rueckwaerts bis konvergent
    // 4. Fallback auf P_RFG_max als absolute Obergrenze

    StationaryPlasmaState warm = guess;
    if (!stationary_state_in_bounds(warm))
        warm = stationary_safe_defaults_for_q(Q0);

    double p_lo = 5.0;
    StationarySolveResult s_lo = solve_at_power(p_lo, warm);
    double I_lo = (s_lo.converged && stationary_state_finite_positive(s_lo.state))
                  ? stationary_beam_current_mA(s_lo.state) : 0.0;
    double f_lo = I_lo - I_soll;
    if (s_lo.converged) warm = s_lo.state;

    // Aufwaerts-Suche: verdopple p_hi bis Bracket oder Konvergenzgrenze
    double p_hi = 10.0;
    StationarySolveResult s_hi;
    double I_hi = 0.0, f_hi = -I_soll;
    double p_last_good = p_lo;
    StationaryPlasmaState state_last_good = warm;
    double I_last_good = I_lo;

    const double P_expand_max = std::max(P_RFG_max, 200.0);
    int expand_count = 0;
    bool bracket_found = false;

    while (p_hi <= P_expand_max && expand_count < 20) {
        s_hi = solve_at_power(p_hi, state_last_good);
        bool hi_ok = s_hi.converged && stationary_state_finite_positive(s_hi.state);
        I_hi = hi_ok ? stationary_beam_current_mA(s_hi.state) : 0.0;
        f_hi = I_hi - I_soll;

        std::cout << "POWER_BRACKET " << std::fixed << std::setprecision(4)
                  << p_lo << " " << p_hi << " " << f_lo << " " << f_hi
                  << " hi_conv=" << hi_ok << std::endl;

        if (hi_ok) {
            // Solver konvergiert bei p_hi
            p_last_good = p_hi;
            state_last_good = s_hi.state;
            I_last_good = I_hi;

            if (f_lo * f_hi < 0.0) {
                // Bracket gefunden
                bracket_found = true;
                break;
            }
            // Noch kein Bracket → weiter verdoppeln
            p_hi = std::min(p_hi * 2.0, P_expand_max + 1.0);
        } else {
            // Solver scheitert bei p_hi → Konvergenzgrenze erreicht.
            // Der Bracket muss zwischen p_last_good und p_hi liegen
            // (oder die Physik hat keine Loesung bei I_soll).
            // Halbiere rueckwaerts um den hoechsten konvergenten Punkt zu finden.
            double p_probe = 0.5 * (p_last_good + p_hi);
            for (int refine = 0; refine < 8; ++refine) {
                StationarySolveResult s_probe = solve_at_power(p_probe, state_last_good);
                bool probe_ok = s_probe.converged && stationary_state_finite_positive(s_probe.state);
                double I_probe = probe_ok ? stationary_beam_current_mA(s_probe.state) : 0.0;

                if (probe_ok) {
                    p_last_good = p_probe;
                    state_last_good = s_probe.state;
                    I_last_good = I_probe;
                    p_probe = 0.5 * (p_probe + p_hi);  // Suche weiter oben
                } else {
                    p_hi = p_probe;
                    p_probe = 0.5 * (p_last_good + p_probe);  // Suche weiter unten
                }
            }
            // Setze p_hi auf den hoechsten konvergenten Punkt
            p_hi = p_last_good;
            s_hi = solve_at_power(p_hi, state_last_good);
            I_hi = (s_hi.converged && stationary_state_finite_positive(s_hi.state))
                   ? stationary_beam_current_mA(s_hi.state) : I_last_good;
            f_hi = I_hi - I_soll;

            if (f_lo * f_hi < 0.0) bracket_found = true;
            break;
        }
        expand_count++;
    }

    debug_emit(2, "POWER_BRACKET_INIT",
        "p_lo=" + std::to_string(p_lo) + " I_lo=" + std::to_string(I_lo) +
        " p_hi=" + std::to_string(p_hi) + " I_hi=" + std::to_string(I_hi) +
        " lo_conv=" + std::to_string(s_lo.converged) +
        " hi_conv=" + std::to_string((int)bracket_found) +
        " p_last_good=" + std::to_string(p_last_good));

    // Pruefen: Bracket gefunden?
    if (!bracket_found) {
        out.hit_limit = true;
        out.P_trial_last = p_hi;
        StationarySolveResult& best_s = (std::fabs(f_lo) < std::fabs(f_hi)) ? s_lo : s_hi;
        if (stationary_state_finite_positive(best_s.state)) {
            out.state = best_s.state;
            out.rf = best_s.rf;
            out.I_mA = stationary_beam_current_mA(best_s.state);
            out.err_mA = I_soll - out.I_mA;
        }
        // Klassifikation: Wenn beide Bracket-Enden konvergiert sind,
        // dann hat die Physik kein Gleichgewicht bei I_soll → NO_PHYSICAL.
        // Wenn mindestens ein Ende numerisch gescheitert ist → NUMERICAL.
        bool both_converged = s_lo.converged && s_hi.converged;
        if (both_converged) {
            out.fail_type = SolveFailType::NO_PHYSICAL_SOLUTION;
            out.reason = "I(P) does not cross I_soll in [" + std::to_string(p_lo) + ", "
                         + std::to_string(p_hi) + "] W (I_lo=" + std::to_string(I_lo)
                         + " I_hi=" + std::to_string(I_hi) + ")";
        } else {
            out.fail_type = SolveFailType::NUMERICAL_FAIL;
            out.reason = "bracket incomplete (lo_conv=" + std::to_string(s_lo.converged)
                         + " hi_conv=" + std::to_string(s_hi.converged)
                         + " I_lo=" + std::to_string(I_lo)
                         + " I_hi=" + std::to_string(I_hi) + ")";
        }
        debug_emit(1, "TARGET_CURRENT_FAIL", out.reason, true);
        return out;
    }

    // --- Bisection + Regula Falsi auf P_RFG ---
    //
    // Abbruchkriterien (Prioritaet):
    // 1. |I_mA - I_soll| < power_tol_mA              → strict convergence
    // 2. |I_mA - I_soll| < 10*power_tol_mA UND
    //    Intervallbreite < 0.1 W                      → near-hit acceptance
    // 3. Plateau: |p_hi - p_lo| < 0.01 W             → P-Intervall erschoepft
    // 4. iter >= max_iter                             → Timeout mit Restdiagnose
    //
    // Bei 2-4: bestes Ergebnis akzeptieren, wenn |error| < 1 mA.
    // Sonst: sauberer FAIL mit Diagnose.

    const int bisect_max_iter = 60;
    const double near_hit_tol = 10.0 * power_tol_mA;   // 0.5 mA
    const double P_interval_min = 0.01;                 // 0.01 W

    StationaryPlasmaState last_good = s_lo.converged ? s_lo.state : warm;

    // Bestes Ergebnis ueber alle Iterationen merken
    double best_abs_error = std::numeric_limits<double>::infinity();
    StationarySolveResult best_s;
    double best_p = 0.0, best_I = 0.0;
    int inner_fail_count = 0;

    for (int iter = 0; iter < bisect_max_iter; ++iter) {
        // Secant-Schritt mit Bisection-Fallback
        double p_mid;
        double denom_sec = f_hi - f_lo;
        if (std::fabs(denom_sec) > 1e-12) {
            p_mid = p_lo - f_lo * (p_hi - p_lo) / denom_sec;
            if (!(p_mid > p_lo && p_mid < p_hi) || !std::isfinite(p_mid))
                p_mid = 0.5 * (p_lo + p_hi);
        } else {
            p_mid = 0.5 * (p_lo + p_hi);
        }

        std::cout << "PID_START " << iter << " "
                  << std::fixed << std::setprecision(4) << p_mid << std::endl;

        StationarySolveResult s_mid = solve_at_power(p_mid, last_good);
        if (!s_mid.converged) {
            // Fallback: Bisection-Mitte mit safe guess
            p_mid = 0.5 * (p_lo + p_hi);
            s_mid = solve_at_power(p_mid, stationary_safe_defaults_for_q(Q0));
            if (!s_mid.converged) {
                inner_fail_count++;
                // Nicht sofort aufgeben — Intervall halbieren und weiter
                if (inner_fail_count >= 5) {
                    out.fail_type = SolveFailType::NUMERICAL_FAIL;
                    out.reason = "inner solve failed " + std::to_string(inner_fail_count)
                                 + "x at P=" + std::to_string(p_mid)
                                 + " (" + s_mid.reason + ")";
                    out.P_trial_last = p_mid;
                    // Bestes bisheriges Ergebnis trotzdem eintragen
                    if (best_abs_error < 1e30) {
                        out.state = best_s.state; out.rf = best_s.rf;
                        out.inner_resid_norm = best_s.resid_norm;
                        out.P_RFG_sol = best_p; out.I_mA = best_I;
                        out.err_mA = I_soll - best_I;
                    }
                    return out;
                }
                // Halbiere das Intervall blind und versuche naechste Iteration
                p_hi = 0.5 * (p_lo + p_hi);
                f_hi = f_lo; // pessimistisch, wird beim naechsten Eval korrigiert
                continue;
            }
        }

        double I_mid = stationary_beam_current_mA(s_mid.state);
        double error = I_soll - I_mid;
        double f_mid = I_mid - I_soll;
        double abs_error = std::fabs(error);

        std::cout << "PID_DONE " << std::fixed << std::setprecision(4)
                  << I_mid << " " << error << " " << p_mid << " "
                  << s_mid.state.Te << " " << s_mid.state.Tg << std::endl;

        last_good = s_mid.state;

        // Bestes Ergebnis tracken
        if (abs_error < best_abs_error) {
            best_abs_error = abs_error;
            best_s = s_mid;
            best_p = p_mid;
            best_I = I_mid;
        }

        out.iterations = iter;
        out.state = s_mid.state;
        out.rf = s_mid.rf;
        out.inner_resid_norm = s_mid.resid_norm;
        out.P_RFG_sol = p_mid;
        out.P_trial_last = p_mid;
        out.I_mA = I_mid;
        out.err_mA = error;

        // --- Abbruchkriterium 1: strikte Konvergenz ---
        if (abs_error < power_tol_mA) {
            out.converged = true;
            out.reason = "ok";
            debug_emit(2, "TARGET_CURRENT_OK",
                "iter=" + std::to_string(iter) +
                " P_sol=" + std::to_string(p_mid) +
                " I_mA=" + std::to_string(I_mid));
            std::cout << "CONVERGED " << iter << std::endl;
            return out;
        }

        double P_width = p_hi - p_lo;

        // --- Abbruchkriterium 2: near-hit bei engem Intervall ---
        if (abs_error < near_hit_tol && P_width < 0.1) {
            // Verwende bestes Ergebnis
            out.converged = true;
            out.state = best_s.state; out.rf = best_s.rf;
            out.inner_resid_norm = best_s.resid_norm;
            out.P_RFG_sol = best_p; out.I_mA = best_I;
            out.err_mA = I_soll - best_I;
            out.reason = "near-hit (err=" + std::to_string(best_abs_error)
                         + " mA, dP=" + std::to_string(P_width) + " W)";
            debug_emit(2, "TARGET_CURRENT_NEARHIT",
                "iter=" + std::to_string(iter) +
                " P_sol=" + std::to_string(best_p) +
                " I_mA=" + std::to_string(best_I) +
                " err=" + std::to_string(best_abs_error));
            std::cout << "CONVERGED " << iter << std::endl;
            return out;
        }

        // --- Abbruchkriterium 3: P-Intervall erschoepft (Plateau) ---
        if (P_width < P_interval_min) {
            if (best_abs_error < 1.0) {
                // Nahe genug — akzeptieren
                out.converged = true;
                out.state = best_s.state; out.rf = best_s.rf;
                out.inner_resid_norm = best_s.resid_norm;
                out.P_RFG_sol = best_p; out.I_mA = best_I;
                out.err_mA = I_soll - best_I;
                out.reason = "plateau-accept (err=" + std::to_string(best_abs_error)
                             + " mA, dP=" + std::to_string(P_width) + " W)";
                std::cout << "CONVERGED " << iter << std::endl;
            } else {
                out.fail_type = SolveFailType::NO_PHYSICAL_SOLUTION;
                out.reason = "P-plateau without I_soll match (best_err="
                             + std::to_string(best_abs_error) + " mA)";
            }
            return out;
        }

        // Bracket-Update
        if (f_lo * f_mid <= 0.0) { p_hi = p_mid; f_hi = f_mid; }
        else                      { p_lo = p_mid; f_lo = f_mid; }
    }

    // --- Abbruchkriterium 4: max iter erreicht ---
    // Wenn bestes Ergebnis nahe genug, trotzdem akzeptieren
    if (best_abs_error < near_hit_tol) {
        out.converged = true;
        out.state = best_s.state; out.rf = best_s.rf;
        out.inner_resid_norm = best_s.resid_norm;
        out.P_RFG_sol = best_p; out.I_mA = best_I;
        out.err_mA = I_soll - best_I;
        out.reason = "max-iter-accept (err=" + std::to_string(best_abs_error) + " mA)";
        std::cout << "CONVERGED " << bisect_max_iter << std::endl;
        return out;
    }

    // Echtes Scheitern: logge Restdiagnose
    out.fail_type = SolveFailType::NUMERICAL_FAIL;
    {
        std::ostringstream os;
        os << "power bisection max iter"
           << " p_lo=" << std::fixed << std::setprecision(4) << p_lo
           << " p_hi=" << p_hi
           << " dP=" << (p_hi - p_lo)
           << " best_I=" << best_I
           << " best_err=" << best_abs_error
           << " I_lo=" << (p_lo > 0 ? std::to_string(I_lo) : "?")
           << " I_hi=" << (p_hi > 0 ? std::to_string(I_hi) : "?");
        out.reason = os.str();
    }
    out.P_trial_last = out.P_RFG_sol;
    debug_emit(1, "BISECTION_MAXITER", out.reason, true);
    return out;
}



// ============================================================
// main
// ============================================================
int main(int argc, char** argv) {
    std::string configFile = "params.txt";
    if (argc >= 2) configFile = argv[1];

    auto cfg = loadConfig(configFile);
    Const::applyConfig(cfg);

    auto start_time = std::chrono::high_resolution_clock::now();

    int count_converged = 0;
    int count_no_physical = 0;
    int count_numerical_fail = 0;

    std::string red   = "\033[31m";
    std::string green = "\033[32m";
    std::string cyan  = "\033[36m";
    std::string reset = "\033[0m";
    // HINWEIS: system("reset") entfernt – stoert subprocess-Ausgaben

    try {
        cout << green
             << "####################################################\n"
             << "#          Global Plasma Model – " << gas_species << "\n"
             << "#   Solver: Newton (stationaer)                    #\n"
             << "####################################################\n"
             << reset << endl;

        cout << "GAS_SPECIES " << gas_species << endl;
        cout << "GAS_MASS " << std::scientific << std::setprecision(6) << M << " kg" << endl;
        cout << "SOLVE_MODE " << solve_mode << " "
             << (solve_mode == 2 ? "selbstkonsistent" : "fester_strahlstrom") << endl;

        // Rate model preset name
        const char* rate_model_name = (rate_model == 2) ? "Full tabulated (Biagi/LXCat)"
                                    : (rate_model == 1) ? "Conservative tabulated (Kiz+Kex tab, Kel legacy)"
                                    :                     "Legacy (paper-compatible)";
        cout << "RATE_MODEL " << rate_model << " " << rate_model_name << endl;

        // Cross-section Basispfad fuer das gewaehlte Gas
        std::string cs_base = "cross_sections/" + gas_species + "/";

        // Kel-Tabelle laden falls tabellierter Modus aktiv
        if (elastic_model == 1) {
            std::string kel_path = cs_base + "kel_table.csv";
            if (load_kel_table(kel_path)) {
                cout << "ELASTIC_MODEL tabulated (" << g_kel_table.size() << " Te-Punkte)" << endl;
            } else {
                cerr << "WARNUNG: " << kel_path << " nicht gefunden, verwende Legacy-Modus" << endl;
                elastic_model = 0;
            }
        }
        if (elastic_model == 0) {
            cout << "ELASTIC_MODEL legacy (constant Kel=" << kel_constant << " m^3/s)" << endl;
        }

        // Kiz-Tabelle laden falls tabellierter Modus aktiv
        if (ionization_model == 1) {
            std::string kiz_path = cs_base + "kiz_table.csv";
            if (load_kiz_table(kiz_path)) {
                cout << "IONIZATION_MODEL tabulated (" << g_kiz_table.size() << " Te-Punkte)" << endl;
            } else {
                cerr << "WARNUNG: " << kiz_path << " nicht gefunden, verwende Legacy-Modus" << endl;
                ionization_model = 0;
            }
        }
        if (ionization_model == 0) {
            cout << "IONIZATION_MODEL legacy (Chabert Polynomfit)" << endl;
        }

        // Kex-Tabelle laden falls tabellierter Modus aktiv
        if (excitation_model == 1) {
            std::string kex_path = cs_base + "kex_table.csv";
            if (load_kex_table(kex_path)) {
                cout << "EXCITATION_MODEL tabulated (" << g_kex_table.size() << " Te-Punkte)" << endl;
            } else {
                cerr << "WARNUNG: " << kex_path << " nicht gefunden, verwende Legacy-Modus" << endl;
                excitation_model = 0;
            }
        }
        if (excitation_model == 0) {
            cout << "EXCITATION_MODEL legacy (Chabert Arrhenius)" << endl;
        }

        g_log_rows.clear();
        g_log_events.clear();

        std::ofstream datei("output_kh.txt");
        if (!datei) { cerr << "Fehler: Ausgabedatei nicht oeffenbar!" << endl; return 1; }

        datei << "Method, Q0sccm, Te, Tg, n, ng, iondeg, P_RFG, P_abs, I_extr_mA, "
              << "collision_freq, R_induktiv, I_coil, epsilon_p_real, epsilon_p_imag, "
              << "u_Bohm, J_i, zeta, gamma, xi, eta, plasmafrequenz, frequency_MHz, "
              << "thrust_ions_mN, thrust_atoms_mN, thrust_total_mN, icp_power_efficiency, P_RF_W, "
              << "n_eff, density_profile_factor\n";

        // Zustandsvariablen fuer den Q0-Sweep
        double Te = 3.75, Tg = 300.0, n = 1.0e17, ng = 1.0e19;
        double P_abs = 0.0, R_induktiv = 0.0, I_coil = 0.0;
        double P_RFG_start = P_RFG;

        StationaryPlasmaState stationary_prev = stationary_safe_defaults_for_q(Q0sccm_start * SCCM_TO_PPS);
        bool stationary_have_prev = false;
        StationaryPlasmaState stationary_last_good{};
        bool stationary_have_good = false;
        double stationary_last_good_q0 = 0.0;

        // Zaehler sind im aeusseren main()-Scope deklariert (vor try)

        // BUG-FIX: crash_of_radius-Logik war im Original invertiert.
        // Rc = Spulenradius, R = Kammerradius.
        // Spule muss AUSSERHALB liegen (Rc >= R).
        if (Rc < R) {
            cout << red << "FEHLER: Spulenradius Rc=" << Rc << " < Kammerradius R=" << R
                 << " – Spule liegt innerhalb der Kammer! Programm stoppt." << reset << endl;
            return 1;
        }

        // Child-Langmuir-Limit: Stromdichte und zugehoerige RF-Leistung
        // J_CL ist bereits berechnet in Const::applyConfig().
        // Die CL-Leistung ist der Punkt, an dem J_i = J_CL, also Gamma_i = J_CL/e.
        // Bei CL-Limit: I_beam = J_CL * Ai → P_RF_CL haengt vom Zustand ab.
        // Wir geben J_CL als Referenzwert aus; P_RF_CL wird aus konvergierten
        // Punkten interpoliert (wenn J_i den Wert ueberschreitet).
        cout << "CL_LIMIT " << scientific << setprecision(6)
             << J_CL << " " << J_CL * Ai * 1000.0 << endl;
        // CL_LIMIT <J_CL_A_per_m2> <I_CL_mA>

        // Aeussere Schleife: Q0sccm variieren
        for (int jj = 0; jj < jjmax; ++jj) {
            Q0sccm = Q0sccm_start + jj * Q0sccm_step;
            Q0     = Q0sccm * SCCM_TO_PPS;
            cout << "Q0_STEP " << fixed << setprecision(4) << Q0sccm
                 << " " << (jj+1) << " " << jjmax << endl;
            debug_emit(2, "Q0_STEP", "jj=" + std::to_string(jj) + " Q0sccm=" + std::to_string(Q0sccm) + " Q0=" + std::to_string(Q0));

            {
                StationaryPlasmaState guess;
                if (stationary_have_good && std::fabs(Q0sccm - stationary_last_good_q0) <= 20.0 * Q0sccm_step) {
                    guess = stationary_last_good;
                    debug_emit(2, "GUESS_SOURCE", "last_good (Q0_good=" + std::to_string(stationary_last_good_q0) + ")");
                } else if (stationary_have_prev) {
                    guess = stationary_prev;
                    debug_emit(2, "GUESS_SOURCE", "prev" + std::string(stationary_have_good ? " (last_good too far)" : " (no last_good yet)"));
                } else {
                    guess = stationary_safe_defaults_for_q(Q0);
                    debug_emit(2, "GUESS_SOURCE", "safe_defaults");
                }
                if (!stationary_state_finite_positive(guess)) guess = stationary_safe_defaults_for_q(Q0);
                debug_emit(2, "STATIONARY_GUESS", state_summary(guess.n, guess.ng, guess.Te, guess.Tg, 0.0, 0.0));

                StationaryPowerSolveResult ps;
                if (solve_mode == 2) {
                    // Selbstkonsistenter Modus: P_RFG fest, I ergibt sich
                    ps = stationary_solve_at_fixed_power(P_RFG, guess);
                } else {
                    // Standardmodus: I_soll fest, P_RFG wird gesucht
                    ps = stationary_solve_for_target_current(guess);
                }
                if (!stationary_power_result_valid(ps)) {
                    // Differenzierte Fehlerausgabe
                    if (ps.fail_type == SolveFailType::NO_PHYSICAL_SOLUTION) {
                        count_no_physical++;
                        cout << "NO_PHYSICAL_SOLUTION " << jj << " " << fixed << setprecision(4)
                             << Q0sccm << " " << ps.reason
                             << " I_best=" << (std::isfinite(ps.I_mA) ? ps.I_mA : -1.0)
                             << " P_max_tried=" << (std::isfinite(ps.P_trial_last) ? ps.P_trial_last : -1.0)
                             << endl;
                    } else {
                        count_numerical_fail++;
                        cout << "NUMERICAL_FAIL " << jj << " " << fixed << setprecision(4)
                             << Q0sccm << " " << ps.reason
                             << " P_try=" << (std::isfinite(ps.P_trial_last) ? ps.P_trial_last : -1.0)
                             << " I_mA=" << (std::isfinite(ps.I_mA) ? ps.I_mA : -1.0)
                             << endl;
                    }
                    debug_emit(1, ps.fail_type == SolveFailType::NO_PHYSICAL_SOLUTION
                                  ? "NO_PHYSICAL_SOLUTION" : "NUMERICAL_FAIL",
                               "jj=" + std::to_string(jj) + " Q0sccm=" + std::to_string(Q0sccm)
                               + " reason=" + ps.reason, true);
                    // Fehlgeschlagener Punkt wird NICHT als Startwert gespeichert
                    // (stationary_prev und stationary_last_good bleiben unveraendert)
                    {
                        std::string ft = (ps.fail_type == SolveFailType::NO_PHYSICAL_SOLUTION)
                                         ? "NO_PHYSICAL_SOLUTION" : "NUMERICAL_FAIL";
                        g_log_rows.push_back({jj, Q0sccm, ft, ft,
                            std::isfinite(ps.P_trial_last) ? ps.P_trial_last : 0.0,
                            std::isfinite(ps.I_mA) ? ps.I_mA : 0.0,
                            0, 0, 0, 0, 0,  // Te, Tg, n, ng, resid
                            0, 0, 0, 0, 0,  // iondeg, P_abs, cf, R_ind, I_coil
                            0, 0, 0, 0, 0, 0,  // eps_r, eps_i, uB, Ji, pf, P_RF
                            0, 0, 0,  // thrust_ions, atoms, total
                            0, 0, 0, 0,  // icp_eff, gamma, xi, eta
                            ps.reason, false});
                        simlog_add_event(jj, Q0sccm, ft + ": " + ps.reason);
                    }
                    continue;
                }
                count_converged++;

                stationary_prev = ps.state;
                stationary_have_prev = true;
                if (std::isfinite(ps.inner_resid_norm) && ps.inner_resid_norm < newton_tol) {
                    stationary_last_good = ps.state;
                    stationary_have_good = true;
                    stationary_last_good_q0 = Q0sccm;
                }
                P_RFG = ps.P_RFG_sol;
                P_abs = ps.rf.P_abs;
                R_induktiv = ps.rf.R_ind;
                I_coil = ps.rf.I_coil;
                n = ps.state.n;
                ng = ps.state.ng;
                Te = ps.state.Te;
                Tg = ps.state.Tg;

                DerivedQuantities dq = compute_derived(n, ng, Te, Tg, R_induktiv, I_coil, P_abs);

                if (std::isfinite(ps.inner_resid_norm) && ps.inner_resid_norm >= newton_tol) {
                    cout << "SOFT_ACCEPT " << fixed << setprecision(4)
                         << Q0sccm << " " << scientific << ps.inner_resid_norm << endl;
                }

                cout << "RESULT " << scientific << setprecision(4)
                     << n << " " << ng << " " << fixed << setprecision(3)
                     << Te << " " << Tg << " " << dq.I_extr_mA << " " << P_RFG << endl;

                cout << "RESULT_EXT " << fixed << setprecision(6)
                     << Q0sccm << " " << P_RFG << " "
                     << scientific << setprecision(4)
                     << n << " " << ng << " " << n << " "
                     << fixed << setprecision(4)
                     << dq.T_i_N*1e3 << " " << dq.T_n_N*1e3 << " " << dq.T_total_N*1e3 << " "
                     << dq.icp_eff << " " << dq.gamma_eff << " " << dq.xi_mN_kW << " " << dq.eta_mass
                     << endl;

                g_log_rows.push_back({jj, Q0sccm, "CONVERGED", "NONE",
                    P_RFG, dq.I_extr_mA, Te, Tg, n, ng, ps.inner_resid_norm,
                    dq.iondeg, P_abs, dq.cf, R_induktiv, I_coil,
                    dq.eps_p_real, dq.eps_p_imag, dq.u_Bohm, dq.J_i, dq.pf, dq.P_RF,
                    dq.T_i_N*1e3, dq.T_n_N*1e3, dq.T_total_N*1e3,
                    dq.icp_eff, dq.gamma_eff, dq.xi_mN_kW, dq.eta_mass,
                    "ok", true});

                emit_csv_row(datei, "Stationary", Q0sccm, n, ng, Te, Tg,
                             P_RFG, P_abs, R_induktiv, I_coil, dq);

                P_RFG = P_RFG_start;
            }
        } // for (jj)

    } catch (const exception& ex) {
        cerr << "Fehler: " << ex.what() << endl;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end_time - start_time).count();
    std::string end_ts = make_timestamp_readable();

    cout << "\nSUMMARY " << count_converged << " " << count_no_physical << " " << count_numerical_fail << endl;
    cout << "SUMMARY_DETAIL converged=" << count_converged
         << " no_physical_solution=" << count_no_physical
         << " numerical_fail=" << count_numerical_fail
         << " total=" << jjmax << endl;
    cout << "\nDer Code hat " << fixed << setprecision(3)
         << elapsed << " Sekunden gebraucht." << endl;

    // ========== Strukturierte Logdatei schreiben ==========
    {
        std::string logname = "simulation_log_" + make_timestamp() + ".txt";
        std::ofstream lf(logname);
        if (lf) {
            cout << "LOG_FILE " << logname << endl;

            // BLOCK 1 — Header
            lf << "==================================================\n"
               << "GLOBAL PLASMA MODEL — SIMULATION LOG\n"
               << "gas_species:     " << gas_species << "\n"
               << "timestamp_start: " << make_timestamp_readable() << "\n"
               << "timestamp_end:   " << end_ts << "\n"
               << "runtime_seconds: " << std::fixed << std::setprecision(3) << elapsed << "\n"
               << "solver_mode:     stationary (LM), solve_mode=" << solve_mode << "\n"
               << "rate_model:      " << rate_model << " ("
               << ((rate_model == 2) ? "Full tabulated" : (rate_model == 1) ? "Conservative tabulated" : "Legacy")
               << ")\n"
               << "config_file:     " << (argc >= 2 ? argv[1] : "params.txt") << "\n"
               << "==================================================\n\n";

            // BLOCK 2 — Parameter
            lf << "--------------------------------------------------\n"
               << "SIMULATION PARAMETERS\n"
               << "--------------------------------------------------\n";
            auto pp = [&](const char* label, double val, const char* unit, const char* key) {
                lf << std::left << std::setw(35) << label
                   << "| " << std::right << std::setw(14) << std::scientific << std::setprecision(6) << val
                   << " | " << std::setw(6) << unit
                   << " | " << key << "\n";
            };
            pp("Radius Entladungsgefaess",       R,         "m",    "R");
            pp("Laenge Entladungsgefaess",        L,         "m",    "L");
            pp("Transparenz Ionen",              betai,     "-",    "betai");
            pp("Transparenz Neutralgas",         betag,     "-",    "betag");
            pp("Anregefrequenz",                 frequency, "Hz",   "frequency");
            pp("Spulenwindungen",                Nw,        "-",    "Nw");
            pp("Ohmscher Widerstand Spule",      R_ohm,     "Ohm",  "R_ohm");
            pp("Spulenradius",                   Rc,        "m",    "Rc");
            pp("Spulenlaenge",                    lc,        "m",    "lc");
            pp("Screengitter-Spannung",          Vgrid,     "V",    "Vgrid");
            pp("Gitterabstand",                  sgrid,     "m",    "sgrid");
            pp("DC-Leistung RFG (Startwert)",    P_RFG,     "W",    "P_RFG");
            pp("Max RF-Leistung",                P_RFG_max, "W",    "P_RFG_max");
            pp("Massenfluss (Startwert)",         Q0sccm,   "sccm", "Q0sccm");
            lf << "--------------------------------------------------\n";
            pp("Q0sccm Start",                   Q0sccm_start, "sccm", "Q0sccm_start");
            pp("Q0sccm Schritt",                 Q0sccm_step,  "sccm", "Q0sccm_step");
            pp("Anzahl Schritte",                (double)jjmax, "-",   "jjmax");
            pp("Ziel-Strahlstrom",               I_soll,    "mA",   "I_soll");
            lf << "--------------------------------------------------\n";
            pp("Volumen V",                      V,         "m^3",  "V");
            pp("Wandflaeche A",                   A,         "m^2",  "A");
            pp("omega",                          omega,     "rad/s","omega");
            pp("k0",                             k_0,       "1/m",  "k_0");
            pp("Ag (Gasaustritt)",               Ag,        "m^2",  "Ag");
            pp("Ai (Ionenaustritt)",             Ai,        "m^2",  "Ai");
            pp("Lambda_0 (therm. Diffusion)",    lambda_0,  "m",    "lambda_0");
            pp("L_coil",                         L_coil,    "H",    "L_coil");
            pp("J_CL (Child-Langmuir)",          J_CL,      "A/m^2","J_CL");
            pp("newton_tol",                     newton_tol,"rel",  "newton_tol");
            lf << std::left << std::setw(35) << "gas_species"
               << "| " << gas_species << "\n";
            lf << std::left << std::setw(35) << "Kel-Modell"
               << "| " << (elastic_model == 1 ? "tabulated (Biagi/LXCat)" : "constant (legacy)")
               << ", kel_constant=" << kel_constant
               << "\n";
            lf << "\n";

            // BLOCK 3 — Ergebnistabelle
            lf << "--------------------------------------------------\n"
               << "RESULT TABLE\n"
               << "--------------------------------------------------\n";
            // Header
            lf << std::left
               << std::setw(5)  << "idx" << "| "
               << std::setw(9)  << "Q0_sccm" << "| "
               << std::setw(22) << "status" << "| "
               << std::setw(14) << "P_sol_W" << "| "
               << std::setw(10) << "I_mA" << "| "
               << std::setw(8)  << "Te_eV" << "| "
               << std::setw(8)  << "Tg_K" << "| "
               << std::setw(12) << "n_m-3" << "| "
               << std::setw(12) << "ng_m-3" << "| "
               << std::setw(10) << "resid" << "| "
               << "note" << "\n";
            lf << std::string(130, '-') << "\n";

            for (const auto& row : g_log_rows) {
                lf << std::left << std::setw(5) << row.idx << "| ";
                lf << std::fixed << std::setprecision(4) << std::setw(9) << row.Q0sccm << "| ";
                lf << std::left << std::setw(22) << row.status << "| ";
                if (row.has_data) {
                    lf << std::right << std::fixed << std::setprecision(2) << std::setw(14) << row.P_sol << "| "
                       << std::setprecision(3) << std::setw(10) << row.I_mA << "| "
                       << std::setprecision(3) << std::setw(8) << row.Te << "| "
                       << std::setprecision(1) << std::setw(8) << row.Tg << "| "
                       << std::scientific << std::setprecision(2) << std::setw(12) << row.n << "| "
                       << std::setw(12) << row.ng << "| "
                       << std::setprecision(1) << std::setw(10) << row.resid << "| ";
                } else {
                    // Fehlschlag: Teildaten oder -
                    lf << std::right << std::setw(14) << "-" << "| "
                       << std::setw(10) << (row.I_mA > 0 ? std::to_string(row.I_mA).substr(0,8) : "-") << "| "
                       << std::setw(8) << "-" << "| "
                       << std::setw(8) << "-" << "| "
                       << std::setw(12) << "-" << "| "
                       << std::setw(12) << "-" << "| "
                       << std::setw(10) << "-" << "| ";
                }
                lf << std::left << row.note << "\n";
            }
            lf << "\n";

            // BLOCK 4 — Ereignisdetails
            lf << "--------------------------------------------------\n"
               << "EVENT DETAILS\n"
               << "--------------------------------------------------\n";
            int last_idx = -1;
            for (const auto& ev : g_log_events) {
                if (ev.idx != last_idx) {
                    lf << "\n[Q0 idx=" << ev.idx << " | Q0_sccm="
                       << std::fixed << std::setprecision(4) << ev.Q0sccm << "]\n";
                    last_idx = ev.idx;
                }
                lf << "  - " << ev.message << "\n";
            }
            if (g_log_events.empty()) lf << "(keine Ereignisse)\n";
            lf << "\n";

            // BLOCK 5 — Summary
            lf << "--------------------------------------------------\n"
               << "SUMMARY\n"
               << "--------------------------------------------------\n";
            lf << std::left << std::setw(30) << "total_points"         << "| " << jjmax << "\n"
               << std::setw(30) << "converged"              << "| " << count_converged << "\n"
               << std::setw(30) << "no_physical_solution"   << "| " << count_no_physical << "\n"
               << std::setw(30) << "numerical_fail"         << "| " << count_numerical_fail << "\n"
               << std::setw(30) << "runtime_seconds"        << "| " << std::fixed << std::setprecision(3) << elapsed << "\n";

            // Min/Max der erfolgreichen Punkte
            double p_min=1e30, p_max=-1e30, te_min_v=1e30, te_max_v=-1e30, i_min=1e30, i_max=-1e30;
            double q0_first_ok = -1, q0_last_ok = -1;
            for (const auto& row : g_log_rows) {
                if (!row.has_data) continue;
                if (q0_first_ok < 0) q0_first_ok = row.Q0sccm;
                q0_last_ok = row.Q0sccm;
                p_min  = std::min(p_min,  row.P_sol);
                p_max  = std::max(p_max,  row.P_sol);
                te_min_v = std::min(te_min_v, row.Te);
                te_max_v = std::max(te_max_v, row.Te);
                i_min  = std::min(i_min,  row.I_mA);
                i_max  = std::max(i_max,  row.I_mA);
            }
            if (count_converged > 0) {
                lf << std::setw(30) << "q0_first_success_sccm"  << "| " << std::fixed << std::setprecision(4) << q0_first_ok << "\n"
                   << std::setw(30) << "q0_last_success_sccm"   << "| " << q0_last_ok << "\n"
                   << std::setw(30) << "P_solution_range_W"     << "| " << std::setprecision(2) << p_min << " .. " << p_max << "\n"
                   << std::setw(30) << "Te_range_eV"            << "| " << std::setprecision(3) << te_min_v << " .. " << te_max_v << "\n"
                   << std::setw(30) << "I_mA_range"             << "| " << std::setprecision(3) << i_min << " .. " << i_max << "\n";
            }
            lf << "--------------------------------------------------\n\n";

            // BLOCK 6 — Maschinenlesbare Datentabelle (CSV mit | Separator)
            // Parser-Hinweis: Zeilen die mit "DATA|" beginnen sind Datensaetze.
            // Zeile mit "DATA_HEADER|" ist der Spaltenkopf.
            lf << "--------------------------------------------------\n"
               << "MACHINE-READABLE DATA TABLE\n"
               << "--------------------------------------------------\n";
            lf << "# CL_LIMIT_J_A_per_m2=" << std::scientific << std::setprecision(6) << J_CL << "\n";
            lf << "# CL_LIMIT_I_mA=" << std::fixed << std::setprecision(4) << J_CL * Ai * 1000.0 << "\n";
            lf << "DATA_HEADER|idx|Q0sccm|status|P_RFG_W|P_abs_W|P_RF_W|I_mA|Te_eV|Tg_K"
               << "|n_m3|ng_m3|iondeg_pct|collision_freq|R_induktiv_Ohm|I_coil_A"
               << "|eps_p_real|eps_p_imag|u_Bohm_ms|J_i_A_m2|plasmafrequenz_rad_s"
               << "|thrust_ions_mN|thrust_atoms_mN|thrust_total_mN"
               << "|icp_power_efficiency|gamma_thrust_eff|xi_mN_per_kW|eta_mass_util"
               << "|n_eff|density_profile_factor"
               << "|resid|note\n";

            for (const auto& row : g_log_rows) {
                lf << "DATA|" << row.idx << "|"
                   << std::fixed << std::setprecision(4) << row.Q0sccm << "|"
                   << row.status << "|";
                if (row.has_data) {
                    lf << std::fixed << std::setprecision(4) << row.P_sol << "|"
                       << row.P_abs << "|" << row.P_RF << "|"
                       << std::setprecision(4) << row.I_mA << "|"
                       << std::setprecision(4) << row.Te << "|"
                       << std::setprecision(2) << row.Tg << "|"
                       << std::scientific << std::setprecision(6)
                       << row.n << "|" << row.ng << "|"
                       << std::fixed << std::setprecision(4)
                       << row.iondeg << "|"
                       << std::scientific << std::setprecision(6)
                       << row.collision_freq << "|"
                       << row.R_induktiv << "|"
                       << row.I_coil << "|"
                       << row.eps_p_real << "|"
                       << row.eps_p_imag << "|"
                       << row.u_Bohm << "|"
                       << row.J_i << "|"
                       << row.plasmafrequenz << "|"
                       << std::fixed << std::setprecision(6)
                       << row.thrust_ions_mN << "|"
                       << row.thrust_atoms_mN << "|"
                       << row.thrust_total_mN << "|"
                       << std::setprecision(6)
                       << row.icp_eff << "|"
                       << row.gamma_eff << "|"
                       << row.xi_mN_kW << "|"
                       << row.eta_mass << "|"
                       << density_profile_factor * row.n << "|"
                       << density_profile_factor << "|"
                       << std::scientific << row.resid;
                } else {
                    // 28 leere Felder fuer fehlgeschlagene Punkte (26 + n_eff + dpf)
                    for (int f = 0; f < 28; ++f) lf << "|";
                }
                lf << "|" << row.note << "\n";
            }
            lf << "--------------------------------------------------\n";
        }
    }

    return 0;
}
