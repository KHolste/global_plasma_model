#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <stdexcept>
#include <chrono>
#include <thread>
#include <math.h>
#include <random>
#include <numeric>
#include <complex>
#include <climits>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <map>
#include <string>
#include <limits>
#include "bessel-library.hpp"

using namespace std::literals::complex_literals;
using namespace std;

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

static bool bad_value(double x) {
    return !std::isfinite(x);
}

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
std::map<std::string, double> loadConfig(const std::string& filename) {
    std::map<std::string, double> cfg;
    std::ifstream f(filename);
    if (!f) {
        std::cerr << "[Config] Datei '" << filename << "' nicht gefunden – Standardwerte werden verwendet." << std::endl;
        return cfg;
    }
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string key; double val;
        if (ss >> key >> val) cfg[key] = val;
    }
    return cfg;
}
#define CFG_LOAD(cfg, name) if ((cfg).count(#name)) (name) = (cfg).at(#name)

// ============================================================
// Physikalische Konstanten und Triebwerksparameter
// ============================================================
namespace Const {
    double me       = 9.10938215e-31;
    double e        = 1.602176487e-19;
    double M        = 2.1801711e-25;
    double kB       = 1.3806504e-23;
    double epsilon0 = 8.854187e-12;
    double pi       = 3.141592653589793;
    double c        = 299792458.0;
    double mu_0     = 1.256637061e-6;

    double Tg0    = 293.00;
    double Eiz    = 1.943408035e-18;
    double Eexc   = 1.858524725e-18;
    double sigma_i= 1.0e-18;
    double sigmae = 1.0e-18;
    double kappa  = 0.0057;
    double conv   = 11604.5250061657;

    // Triebwerksparameter (Standardwerte KK)
    double R        = 0.02;
    double L        = 0.04;
    double betai    = 0.5;
    double betag    = 0.05145;
    double frequency= 2.5e6;
    double omega    = 2 * pi * frequency;
    double Nw       = 6.00;
    double R_ohm    = 0.36;
    double Rc       = 0.02;
    double lc       = 0.04;
    double Vgrid    = 1500.00;
    double sgrid    = 0.001;
    double P_RFG    = 18.00;
    double Q0sccm   = 0.475;
    double Q0       = Q0sccm * 4.477962312e17;
    double lambda_0 = R / 2.405 + L / pi;
    double L_coil   = mu_0 * pi * Rc * Rc * Nw * Nw / lc;

    double A    = 2 * pi * R * R + 2 * pi * R * L;
    double Ag   = betag * pi * R * R;
    double Ai   = betai * pi * R * R;
    double V    = pi * R * R * L;
    double k_0  = omega / c;
    double J_CL = (4.0/9.0) * epsilon0 * sqrt(2*e/M) * pow(Vgrid, 1.5) / (sgrid*sgrid);

    // Solver-Parameter
    int    method       = 1;      // 1=Euler, 2=RK4, 3=RK45, 4=Newton (stationaer)
    double I_soll       = 15.00;  // Ziel-Strahlstrom (mA)
    double Q0sccm_start = 0.27;
    double Q0sccm_step  = 0.01;
    int    jjmax        = 73;

    // Stationaerer Newton-Solver
    double P_RFG_max    = 80.00;
    int    newton_max_iter   = 45;
    double newton_tol        = 1e-2;   // Relative Toleranz (mit physik. Skalierung: ~1% Fehler)
    double power_tol_mA      = 0.05;
    int    power_max_iter    = 35;
    double power_min         = 1.0;

    // physikalische Schranken / Daempfung
    double n_min   = 1e12, n_max   = 1e20;
    double ng_min  = 1e16, ng_max  = 1e22;
    double Te_min  = 0.3,  Te_max  = 20.0;
    double Tg_min  = 200.0, Tg_max = 2500.0;
    double newton_max_log_step = 0.35;
    double newton_fd_eps = 1e-5;

    // Pseudo-transient continuation (robuster Warmstart vor Newton)
    int    ptc_max_iter       = 80;
    double ptc_start_gain     = 0.20;
    double ptc_min_gain       = 1e-4;
    double ptc_switch_merit   = 5e-3;
    double ptc_accept_ratio   = 0.98;

    // "brauchbar genug"-Akzeptanz fuer stationaere Startloesungen:
    // verhindert, dass der aeussere Power-Solver wegen eines zu strengen
    // inneren Newton-Kriteriums sofort mit I_mA=-1 abbricht.
    double stationary_soft_abs_resid = 5.0;
    double stationary_soft_rel_improve = 0.80;

    void applyConfig(const std::map<std::string, double>& cfg) {
        CFG_LOAD(cfg, R);       CFG_LOAD(cfg, L);        CFG_LOAD(cfg, betai);
        CFG_LOAD(cfg, betag);   CFG_LOAD(cfg, frequency);CFG_LOAD(cfg, Nw);
        CFG_LOAD(cfg, R_ohm);   CFG_LOAD(cfg, Rc);       CFG_LOAD(cfg, lc);
        CFG_LOAD(cfg, Vgrid);   CFG_LOAD(cfg, sgrid);    CFG_LOAD(cfg, P_RFG);
        CFG_LOAD(cfg, P_RFG_max);
        CFG_LOAD(cfg, Q0sccm);  CFG_LOAD(cfg, I_soll);
        CFG_LOAD(cfg, Q0sccm_start); CFG_LOAD(cfg, Q0sccm_step);
        CFG_LOAD(cfg, power_tol_mA); CFG_LOAD(cfg, power_min);
        CFG_LOAD(cfg, n_min); CFG_LOAD(cfg, n_max); CFG_LOAD(cfg, ng_min); CFG_LOAD(cfg, ng_max);
        CFG_LOAD(cfg, Te_min); CFG_LOAD(cfg, Te_max); CFG_LOAD(cfg, Tg_min); CFG_LOAD(cfg, Tg_max);
        CFG_LOAD(cfg, newton_max_log_step); CFG_LOAD(cfg, newton_fd_eps);
        if (cfg.count("jjmax"))  jjmax  = (int)cfg.at("jjmax");
        if (cfg.count("method")) method = (int)cfg.at("method");
        if (cfg.count("newton_max_iter")) newton_max_iter = (int)cfg.at("newton_max_iter");
        if (cfg.count("power_max_iter"))  power_max_iter  = (int)cfg.at("power_max_iter");
        if (cfg.count("debug_level"))     debug_level = std::max(0, std::min(3, (int)cfg.at("debug_level")));

        Q0       = Q0sccm * 4.477962312e17;
        omega    = 2 * pi * frequency;
        lambda_0 = R / 2.405 + L / pi;
        L_coil   = mu_0 * pi * Rc * Rc * Nw * Nw / lc;
        A        = 2 * pi * R * R + 2 * pi * R * L;
        Ag       = betag * pi * R * R;
        Ai       = betai * pi * R * R;
        V        = pi * R * R * L;
        k_0      = omega / c;
        J_CL     = (4.0/9.0) * epsilon0 * sqrt(2*e/M) * pow(Vgrid, 1.5) / (sgrid*sgrid);
    }
}
using namespace Const;

// ============================================================
// Hilfsfunktionen
// ============================================================
double vg(double Tg)       { return sqrt(8 * kB * Tg / (pi * M)); }
double vi(double Ti)       { return sqrt(8 * kB * Ti / (pi * M)); }
double Gamma_g(double ng, double vg_val) { return 0.25 * ng * vg_val; }

double Kex(double Te) { return 1.2921e-13 * exp(-e * 11.6 / (kB * Te * conv)); }

double Kiz(double Te) {
    double TeV = kB * Te * conv / e;
    double K1 = 6.73e-15 * sqrt(TeV) * (3.97 + 0.643*TeV - 0.0368*TeV*TeV) * exp(-12.127/TeV);
    double K2 = 6.73e-15 * sqrt(TeV) * (-0.0001031*TeV*TeV + 6.386*exp(-12.127/TeV));
    return 0.5 * (K1 + K2);
}

double Kel(double Te) {
    return -1.45239e-13 + 2.92063e-13*Te - 7.59346e-14*Te*Te
           + 9.78729e-15*pow(Te,3) - 6.3466e-16*pow(Te,4) + 1.64868e-19*pow(Te,5);
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

// ============================================================
// RF-Leistungskopplung
// ============================================================
void do_the_RF_magic(double& n, double& ng, double& Te,
                     double* I_coil_ptr, double* P_abs_ptr, double* R_induktiv_ptr) {
    double a = plasma_freq(n), b = omega, cf = coll_freq(ng, Te);
    double it = a*a*cf / (b*(b*b+cf*cf));
    double rt = 1 - a*a / (b*b+cf*cf);
    double o_rt = sqrt((sqrt(it*it+rt*rt)+rt)/2.0);
    double o_it = signum(it) * sqrt((sqrt(it*it+rt*rt)-rt)/2.0);

    complex<double> k_1(o_rt*k_0, -o_it*k_0);
    complex<double> kR = R * k_1;
    complex<double> eps_p(rt, -it);
    complex<double> result = (1i * kR * bessel::cyl_j(1,kR)) / (eps_p * bessel::cyl_j(0,kR));

    double R_ind = 2*pi*Nw*Nw / (epsilon0*L*omega) * result.real();
    double Ic    = sqrt(2*P_RFG / (R_ind + R_ohm));
    *R_induktiv_ptr = R_ind;
    *I_coil_ptr     = Ic;
    *P_abs_ptr      = 0.5 * R_ind * Ic * Ic;
}


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
    complex<double> j0 = bessel::cyl_j(0, kR);
    complex<double> j1 = bessel::cyl_j(1, kR);
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

static bool stationary_state_in_bounds(const StationaryPlasmaState& s) {
    if (!stationary_state_finite_positive(s)) return false;
    if (s.n  < n_min  || s.n  > n_max)  return false;
    if (s.ng < ng_min || s.ng > ng_max) return false;
    if (s.Te < Te_min || s.Te > Te_max) return false;
    if (s.Tg < Tg_min || s.Tg > Tg_max) return false;
    double iondeg = s.n / std::max(1.0, s.ng);
    if (!std::isfinite(iondeg) || iondeg < 0.0 || iondeg > 0.5) return false;
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

static std::array<double,4> stationary_residual_raw(const StationaryPlasmaState& s,
                                                    double P_RFG_local,
                                                    StationaryRFState* rf_out = nullptr) {
    StationaryRFState rf = stationary_compute_rf(s.n, s.ng, s.Te, P_RFG_local);
    if (rf_out) *rf_out = rf;
    if (!rf.valid) {
        double nan = std::numeric_limits<double>::quiet_NaN();
        return {nan,nan,nan,nan};
    }

    double P_vol = rf.P_abs / V;

    double r1 = s.n*s.ng*Kiz(s.Te) - s.n*uB(s.Te)*Aeff(lambda_i(s.ng))/V;
    double r2 = Q0/V + s.n*uB(s.Te)*Aeff1(lambda_i(s.ng))/V - s.n*s.ng*Kiz(s.Te) - Gamma_g(s.ng,vg(s.Tg))*Ag/V;

    double Pg1 = 3.0*me/M * kB*(s.Te*conv-s.Tg) * s.n*s.ng*Kel(s.Te);
    double Pg2 = 0.25*M*uB(s.Te)*uB(s.Te) * s.n*s.ng*sigma_i*vi(s.Tg);
    double Pg3 = kappa*(s.Tg-Tg0)/lambda_0 * A/V;
    double r4 = Pg1 + Pg2 - Pg3;

    double P2 = Eiz  * s.n*s.ng*Kiz(s.Te);
    double P3 = Eexc * s.n*s.ng*Kex(s.Te);
    double P4 = 3.0*me/M * kB*(s.Te*conv-s.Tg) * s.n*s.ng*Kel(s.Te);
    double P5 = 7.0*kB*s.Te*conv * s.n*uB(s.Te)*Aeff(lambda_i(s.ng))/V;
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

struct StationaryPowerSolveResult {
    bool converged = false;
    bool hit_limit = false;
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
// 4D-Solver mit aeusserer Power-Bisection
//
// Architektur (zurueck zum physikalisch korrekten Ansatz):
//   Innerer Solver: stationary_solve_lm loest ALLE 4 Gleichungen
//     (r1, r2, r3, r4) gleichzeitig fuer gegebenes P_RFG.
//     r3 koppelt P_RFG an den Plasmazustand — das ist entscheidend.
//   Aeusserer Solver: Bisection auf P_RFG bis I_mA = I_soll.
//
// Die vorherige ng-Dekomposition war physikalisch falsch, weil
// r3 die Plasmadichte ueber die verfuegbare Leistung begrenzt.
// Ohne r3 im inneren Solver gibt es keine obere Schranke fuer n.
// ============================================================

static StationaryPowerSolveResult stationary_solve_for_target_current(const StationaryPlasmaState& guess) {
    StationaryPowerSolveResult out;

    debug_emit(2, "TARGET_CURRENT_BEGIN",
        "Q0sccm=" + std::to_string(Q0 / 4.477962312e17) +
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
    double p_lo = 5.0;
    double p_hi = 100.0;
    StationaryPlasmaState warm = guess;
    if (!stationary_state_in_bounds(warm))
        warm = stationary_safe_defaults_for_q(Q0);

    // Evaluiere untere Grenze
    StationarySolveResult s_lo = solve_at_power(p_lo, warm);
    double I_lo = (s_lo.converged && stationary_state_finite_positive(s_lo.state))
                  ? stationary_beam_current_mA(s_lo.state) : 0.0;
    double f_lo = I_lo - I_soll;

    // Evaluiere obere Grenze
    StationarySolveResult s_hi = solve_at_power(p_hi, warm);
    double I_hi = (s_hi.converged && stationary_state_finite_positive(s_hi.state))
                  ? stationary_beam_current_mA(s_hi.state) : 0.0;
    double f_hi = I_hi - I_soll;

    debug_emit(2, "POWER_BRACKET_INIT",
        "p_lo=" + std::to_string(p_lo) + " I_lo=" + std::to_string(I_lo) +
        " p_hi=" + std::to_string(p_hi) + " I_hi=" + std::to_string(I_hi) +
        " lo_conv=" + std::to_string(s_lo.converged) +
        " hi_conv=" + std::to_string(s_hi.converged));

    // Bracket erweitern falls noetig (P_hi verdoppeln, max 5000 W)
    const double P_expand_max = 5000.0;
    int expand_count = 0;
    while (f_lo * f_hi > 0.0 && p_hi < P_expand_max && expand_count < 15) {
        p_hi = std::min(p_hi * 2.0, P_expand_max);
        s_hi = solve_at_power(p_hi, s_hi.converged ? s_hi.state : warm);
        I_hi = (s_hi.converged && stationary_state_finite_positive(s_hi.state))
               ? stationary_beam_current_mA(s_hi.state) : I_hi;
        f_hi = I_hi - I_soll;

        std::cout << "POWER_BRACKET " << std::fixed << std::setprecision(4)
                  << p_lo << " " << p_hi << " " << f_lo << " " << f_hi << std::endl;
        expand_count++;

        // Warm-start fuer naechste Iteration
        if (s_hi.converged) warm = s_hi.state;
    }

    // Pruefen: Bracket gefunden?
    if (f_lo * f_hi > 0.0) {
        out.hit_limit = true;
        out.reason = "no power bracket (I(P) stays on one side, I_lo="
                     + std::to_string(I_lo) + " I_hi=" + std::to_string(I_hi)
                     + " at P_hi=" + std::to_string(p_hi) + ")";
        out.P_trial_last = p_hi;
        StationarySolveResult& best_s = (std::fabs(f_lo) < std::fabs(f_hi)) ? s_lo : s_hi;
        if (stationary_state_finite_positive(best_s.state)) {
            out.state = best_s.state;
            out.rf = best_s.rf;
            out.I_mA = stationary_beam_current_mA(best_s.state);
            out.err_mA = I_soll - out.I_mA;
        }
        debug_emit(1, "TARGET_CURRENT_FAIL", out.reason, true);
        return out;
    }

    // --- Bisection + Regula Falsi auf P_RFG ---
    StationaryPlasmaState last_good = s_lo.converged ? s_lo.state : warm;
    for (int iter = 0; iter < power_max_iter; ++iter) {
        // Secant-Schritt mit Bisection-Fallback
        double p_mid;
        double denom = f_hi - f_lo;
        if (std::fabs(denom) > 1e-12) {
            p_mid = p_lo - f_lo * (p_hi - p_lo) / denom;
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
                out.reason = "inner solve failed at P=" + std::to_string(p_mid)
                             + " (" + s_mid.reason + ")";
                out.P_trial_last = p_mid;
                return out;
            }
        }

        double I_mid = stationary_beam_current_mA(s_mid.state);
        double error = I_soll - I_mid;
        double f_mid = I_mid - I_soll;

        std::cout << "PID_DONE " << std::fixed << std::setprecision(4)
                  << I_mid << " " << error << " " << p_mid << " "
                  << s_mid.state.Te << " " << s_mid.state.Tg << std::endl;

        last_good = s_mid.state;
        out.iterations = iter;
        out.state = s_mid.state;
        out.rf = s_mid.rf;
        out.inner_resid_norm = s_mid.resid_norm;
        out.P_RFG_sol = p_mid;
        out.P_trial_last = p_mid;
        out.I_mA = I_mid;
        out.err_mA = error;

        if (std::fabs(error) < power_tol_mA) {
            out.converged = true;
            out.reason = "ok";
            debug_emit(2, "TARGET_CURRENT_OK",
                "iter=" + std::to_string(iter) +
                " P_sol=" + std::to_string(p_mid) +
                " I_mA=" + std::to_string(I_mid));
            std::cout << "CONVERGED " << iter << std::endl;
            return out;
        }

        if (f_lo * f_mid <= 0.0) { p_hi = p_mid; f_hi = f_mid; }
        else                      { p_lo = p_mid; f_lo = f_mid; }
    }

    out.reason = "power bisection max iter";
    out.P_trial_last = out.P_RFG_sol;
    return out;
}

// ============================================================
// Ableitungen (auf Energiedichte-Basis)
// ============================================================
double d_n(double n, double ng, double Te) {
    return n*ng*Kiz(Te) - n*uB(Te)*Aeff(lambda_i(ng))/V;
}
double d_ng(double n, double ng, double Te, double Tg) {
    return Q0/V + n*uB(Te)*Aeff1(lambda_i(ng))/V - n*ng*Kiz(Te) - Gamma_g(ng,vg(Tg))*Ag/V;
}
// Gibt dU_g/dt in J/(m^3 s) zurueck
double d_energy_gas(double Te, double Tg, double ng, double n) {
    double Pg1 = 3.0*me/M * kB*(Te*conv-Tg) * n*ng*Kel(Te);
    double Pg2 = 0.25*M*uB(Te)*uB(Te) * n*ng*sigma_i*vi(Tg);
    double Pg3 = kappa*(Tg-Tg0)/lambda_0 * A/V;
    return Pg1 + Pg2 - Pg3;
}
// Gibt dU_e/dt in J/(m^3 s) zurueck
double d_energy_elec(double Te, double Tg, double ng, double n, double P_vol) {
    double P2 = Eiz * n*ng*Kiz(Te);
    double P3 = Eexc * n*ng*Kex(Te);
    double P4 = 3.0*me/M * kB*(Te*conv-Tg) * n*ng*Kel(Te);
    double P5 = 7.0*kB*Te*conv * n*uB(Te)*Aeff(lambda_i(ng))/V;
    return P_vol - (P2+P3+P4+P5);
}

// ============================================================
// Zustandsvektor-Hilfsstruct fuer RK-Verfahren
// ============================================================
struct State {
    double n, ng, Ug, Ue; // Ug=energy_gas, Ue=energy_elec
};

State operator+(const State& a, const State& b) { return {a.n+b.n, a.ng+b.ng, a.Ug+b.Ug, a.Ue+b.Ue}; }
State operator-(const State& a, const State& b) { return {a.n-b.n, a.ng-b.ng, a.Ug-b.Ug, a.Ue-b.Ue}; }
State operator*(double s, const State& a)        { return {s*a.n,   s*a.ng,   s*a.Ug,   s*a.Ue}; }
State operator*(const State& a, double s)        { return s*a; }
State operator-(const State& a)                  { return {-a.n, -a.ng, -a.Ug, -a.Ue}; }

// Rechte Seite des DGL-Systems als Funktion des Zustandsvektors
State rhs(const State& s, double P_vol) {
    double Tg_ = (2.0/3.0)*s.Ug/(s.ng*kB);
    double Te_ = (2.0/3.0)*s.Ue/(s.n *kB*conv);
    return {
        d_n (s.n, s.ng, Te_),
        d_ng(s.n, s.ng, Te_, Tg_),
        d_energy_gas (Te_, Tg_, s.ng, s.n),
        d_energy_elec(Te_, Tg_, s.ng, s.n, P_vol)
    };
}

static bool invalid_updated_state(double new_n, double new_ng, double new_energy_gas, double new_energy_elec,
                                  double new_Te, double new_Tg) {
    return bad_value(new_n) || bad_value(new_ng) || bad_value(new_energy_gas) || bad_value(new_energy_elec) ||
           bad_value(new_Te) || bad_value(new_Tg) ||
           new_n <= 0.0 || new_ng <= 0.0 || new_energy_gas <= 0.0 || new_energy_elec <= 0.0 ||
           new_Te <= 0.0 || new_Tg <= 0.0;
}

static std::string diagnose_updated_state(double new_n, double new_ng, double new_energy_gas, double new_energy_elec,
                                          double new_Te, double new_Tg, double P_vol, int step, int meth,
                                          double old_n, double old_ng, double old_Te, double old_Tg,
                                          double old_Ug, double old_Ue) {
    std::ostringstream os;
    os << "method=" << meth << " step=" << step << " P_vol=" << std::scientific << std::setprecision(6) << P_vol << " ";
    if (bad_value(new_n)) os << "bad=new_n ";
    if (bad_value(new_ng)) os << "bad=new_ng ";
    if (bad_value(new_energy_gas)) os << "bad=new_Ug ";
    if (bad_value(new_energy_elec)) os << "bad=new_Ue ";
    if (bad_value(new_Te)) os << "bad=new_Te ";
    if (bad_value(new_Tg)) os << "bad=new_Tg ";
    if (new_n <= 0.0) os << "bad=new_n<=0 ";
    if (new_ng <= 0.0) os << "bad=new_ng<=0 ";
    if (new_energy_gas <= 0.0) os << "bad=new_Ug<=0 ";
    if (new_energy_elec <= 0.0) os << "bad=new_Ue<=0 ";
    if (new_Te <= 0.0) os << "bad=new_Te<=0 ";
    if (new_Tg <= 0.0) os << "bad=new_Tg<=0 ";
    os << "old{" << state_summary(old_n, old_ng, old_Te, old_Tg, old_Ug, old_Ue) << "} ";
    os << "new{" << state_summary(new_n, new_ng, new_Te, new_Tg, new_energy_gas, new_energy_elec) << "}";
    return os.str();
}

// ============================================================
// Integrator: Euler(1), RK4(2), RK45 Dormand-Prince(3)
// ============================================================
void runGM(double& n, double& ng, double& Te, double& Tg,
           double& energy_gas, double& energy_elec,
           double* n_ptr, double* ng_ptr, double* Tg_ptr, double* Te_ptr,
           int* it_ptr, double* tol_ptr,
           double* P_abs_ptr, double* R_induktiv_ptr, double* I_coil_ptr,
           int meth) {

    debug_emit(2, "RUNGM_BEGIN", "method=" + std::to_string(meth) + " P_RFG=" + std::to_string(P_RFG) + " " + state_summary(n, ng, Te, Tg, energy_gas, energy_elec));

    const int    Nmax = 600000000;
    // Euler stabil bis ~dt=5e-8, RK4 erlaubt ~10x groesseres dt
    const double dt   = (meth == 1) ? 1e-8 : 1e-7;
    const double tol  = 1e-10;
    const int    check_interval = (meth == 1) ? 10000 : 2000;
    const int    print_interval = (meth == 1) ? 5000000 : 1000000;

    if (meth == 1) {
        // ── Euler ──────────────────────────────────────────────────────────
        // Erste Konvergenzprüfung erst nach check_interval Schritten
        // (verhindert Sofort-Konvergenz bei Warmstart)
        for (int i = 0; i < Nmax; ++i) {
            double n_help=n, ng_help=ng, Tg_help=Tg, Te_help=Te;

            do_the_RF_magic(n, ng, Te, I_coil_ptr, P_abs_ptr, R_induktiv_ptr);
            double P_vol = *P_abs_ptr / V;

            double new_n           = n  + d_n(n,ng,Te)          * dt;
            double new_ng          = ng + d_ng(n,ng,Te,Tg)       * dt;
            double new_energy_gas  = energy_gas  + d_energy_gas(Te,Tg,ng,n)       * dt;
            double new_energy_elec = energy_elec + d_energy_elec(Te,Tg,ng,n,P_vol)* dt;
            double new_Tg = (2.0/3.0) * new_energy_gas  / (ng  * kB);
            double new_Te = (2.0/3.0) * new_energy_elec / (n   * kB * conv);

            n=new_n; ng=new_ng; Tg=new_Tg; Te=new_Te;
            energy_gas=new_energy_gas; energy_elec=new_energy_elec;

            if (invalid_updated_state(new_n, new_ng, new_energy_gas, new_energy_elec, new_Te, new_Tg)) {
                std::string diag = diagnose_updated_state(new_n, new_ng, new_energy_gas, new_energy_elec, new_Te, new_Tg, P_vol, i, meth,
                                                          n_help, ng_help, Te_help, Tg_help, energy_gas, energy_elec);
                debug_emit(1, "RUNGM_FAIL", diag, true);
                *it_ptr = -1;
                break;
            }

            if (i % check_interval == 0) {
                double R_ng = abs((new_ng-ng_help)/new_ng);
                double R_n  = abs((new_n -n_help) /new_n);
                double R_Tg = abs((new_Tg-Tg_help)/new_Tg);
                double R_Te = abs((new_Te-Te_help) /new_Te);
                double mean_tol = 0.25*(R_ng+R_n+R_Tg+R_Te);
                if (i % print_interval == 0) {
                    // Strukturierte Zeile fuer GUI-Parser
                    std::cout << "STEP "
                              << std::fixed << std::setprecision(3) << Te << " "
                              << std::scientific << std::setprecision(3) << n << " "
                              << ng << " " << dt << std::endl;
                    std::cout.flush();
                }
                if (R_ng<=tol && R_n<=tol && R_Tg<=tol && R_Te<=tol) {
                    *it_ptr=i; *tol_ptr=mean_tol; break;
                }
            }
        }
        std::cout << endl;
        *n_ptr=n; *ng_ptr=ng; *Tg_ptr=Tg; *Te_ptr=Te;
        debug_emit((*it_ptr == -1) ? 1 : 2, (*it_ptr == -1) ? "RUNGM_END_FAIL" : "RUNGM_END_OK", "method=1 it=" + std::to_string(*it_ptr) + " tol=" + std::to_string(*tol_ptr) + " " + state_summary(n, ng, Te, Tg, energy_gas, energy_elec), *it_ptr == -1);
    }

    if (meth == 2) {
        // ── Vollstaendig gekoppeltes RK4 ────────────────────────────────────
        // Ersten check_interval Schritte immer durchlaufen (kein Sofort-Konvergenz)
        double n_help0=n*0.5, ng_help0=ng*0.5, Tg_help0=Tg*0.5, Te_help0=Te*0.5;
        for (int i = 0; i < Nmax; ++i) {
            double n_help  = (i==0) ? n_help0  : n;
            double ng_help = (i==0) ? ng_help0 : ng;
            double Tg_help = (i==0) ? Tg_help0 : Tg;
            double Te_help = (i==0) ? Te_help0 : Te;

            do_the_RF_magic(n, ng, Te, I_coil_ptr, P_abs_ptr, R_induktiv_ptr);
            double P_vol = *P_abs_ptr / V;

            // Hilfslambdas: Temperatur aus Energiedichte
            auto get_Tg = [&](double ng_, double Ug_) { return (2.0/3.0)*Ug_/(ng_*kB); };
            auto get_Te = [&](double n_,  double Ue_) { return (2.0/3.0)*Ue_/(n_*kB*conv); };

            // k1: Ableitungen am aktuellen Punkt
            double k1_n  = d_n (n,  ng, Te)*dt;
            double k1_ng = d_ng(n,  ng, Te, Tg)*dt;
            double k1_Ug = d_energy_gas (Te, Tg, ng, n)*dt;
            double k1_Ue = d_energy_elec(Te, Tg, ng, n, P_vol)*dt;

            // k2: Ableitungen am Halbschritt (alle Variablen gleichzeitig)
            double n2  = n  + k1_n /2.0,  ng2 = ng + k1_ng/2.0;
            double Ug2 = energy_gas  + k1_Ug/2.0;
            double Ue2 = energy_elec + k1_Ue/2.0;
            double Tg2 = get_Tg(ng2, Ug2), Te2 = get_Te(n2, Ue2);
            double k2_n  = d_n (n2,  ng2, Te2)*dt;
            double k2_ng = d_ng(n2,  ng2, Te2, Tg2)*dt;
            double k2_Ug = d_energy_gas (Te2, Tg2, ng2, n2)*dt;
            double k2_Ue = d_energy_elec(Te2, Tg2, ng2, n2, P_vol)*dt;

            // k3: Ableitungen am zweiten Halbschritt
            double n3  = n  + k2_n /2.0,  ng3 = ng + k2_ng/2.0;
            double Ug3 = energy_gas  + k2_Ug/2.0;
            double Ue3 = energy_elec + k2_Ue/2.0;
            double Tg3 = get_Tg(ng3, Ug3), Te3 = get_Te(n3, Ue3);
            double k3_n  = d_n (n3,  ng3, Te3)*dt;
            double k3_ng = d_ng(n3,  ng3, Te3, Tg3)*dt;
            double k3_Ug = d_energy_gas (Te3, Tg3, ng3, n3)*dt;
            double k3_Ue = d_energy_elec(Te3, Tg3, ng3, n3, P_vol)*dt;

            // k4: Ableitungen am vollen Schritt
            double n4  = n  + k3_n,   ng4 = ng + k3_ng;
            double Ug4 = energy_gas  + k3_Ug;
            double Ue4 = energy_elec + k3_Ue;
            double Tg4 = get_Tg(ng4, Ug4), Te4 = get_Te(n4, Ue4);
            double k4_n  = d_n (n4,  ng4, Te4)*dt;
            double k4_ng = d_ng(n4,  ng4, Te4, Tg4)*dt;
            double k4_Ug = d_energy_gas (Te4, Tg4, ng4, n4)*dt;
            double k4_Ue = d_energy_elec(Te4, Tg4, ng4, n4, P_vol)*dt;

            // Gewichtete Summe
            double new_n           = n           + (k1_n  + 2*k2_n  + 2*k3_n  + k4_n )/6.0;
            double new_ng          = ng          + (k1_ng + 2*k2_ng + 2*k3_ng + k4_ng)/6.0;
            double new_energy_gas  = energy_gas  + (k1_Ug + 2*k2_Ug + 2*k3_Ug + k4_Ug)/6.0;
            double new_energy_elec = energy_elec + (k1_Ue + 2*k2_Ue + 2*k3_Ue + k4_Ue)/6.0;
            double new_Tg = get_Tg(new_ng, new_energy_gas);
            double new_Te = get_Te(new_n,  new_energy_elec);

            n=new_n; ng=new_ng; Tg=new_Tg; Te=new_Te;
            energy_gas=new_energy_gas; energy_elec=new_energy_elec;

            if (invalid_updated_state(new_n, new_ng, new_energy_gas, new_energy_elec, new_Te, new_Tg)) {
                std::string diag = diagnose_updated_state(new_n, new_ng, new_energy_gas, new_energy_elec, new_Te, new_Tg, P_vol, i, meth,
                                                          n_help, ng_help, Te_help, Tg_help, energy_gas, energy_elec);
                debug_emit(1, "RUNGM_FAIL", diag, true);
                *it_ptr = -1;
                break;
            }

            if (i % check_interval == 0) {
                double R_ng = abs((new_ng-ng_help)/new_ng);
                double R_n  = abs((new_n -n_help) /new_n);
                double R_Tg = abs((new_Tg-Tg_help)/new_Tg);
                double R_Te = abs((new_Te-Te_help) /new_Te);
                double mean_tol = 0.25*(R_ng+R_n+R_Tg+R_Te);
                if (i % print_interval == 0) {
                    // Strukturierte Zeile fuer GUI-Parser
                    std::cout << "STEP "
                              << std::fixed << std::setprecision(3) << Te << " "
                              << std::scientific << std::setprecision(3) << n << " "
                              << ng << " " << dt << std::endl;
                    std::cout.flush();
                }
                if (R_ng<=tol && R_n<=tol && R_Tg<=tol && R_Te<=tol) {
                    *it_ptr=i; *tol_ptr=mean_tol; break;
                }
            }
        }
        std::cout << endl;
        *n_ptr=n; *ng_ptr=ng; *Tg_ptr=Tg; *Te_ptr=Te;
        debug_emit((*it_ptr == -1) ? 1 : 2, (*it_ptr == -1) ? "RUNGM_END_FAIL" : "RUNGM_END_OK", "method=2 it=" + std::to_string(*it_ptr) + " tol=" + std::to_string(*tol_ptr) + " " + state_summary(n, ng, Te, Tg, energy_gas, energy_elec), *it_ptr == -1);
    }

    if (meth == 3) {
        // ── RK45 Dormand-Prince mit adaptiver Schrittweite ──────────────────
        const double atol      = 1e-8;
        const double rtol      = 1e-6;
        const double dt_min    = 1e-12;
        const double dt_max    = 1e-5;
        // Konvergenzkriterium: relative Aenderung ueber ein Zeitfenster.
        // Fenstergroesse = check_int * dt_avg ~ 1000 * 2e-6 = 2ms
        const double tol       = 1e-6;   // lockerer als Euler/RK4: passt zu adaptivem dt
        const int    check_int = 1000;   // alle 1000 akzeptierten Schritte pruefen
        const int    print_int = 5000;   // alle 5000 Schritte ausgeben

        do_the_RF_magic(n, ng, Te, I_coil_ptr, P_abs_ptr, R_induktiv_ptr);
        double P_vol = *P_abs_ptr / V;

        State s = {n, ng, energy_gas, energy_elec};
        double dt_rk45 = 1e-9;
        int accepted = 0;
        const int Nmax_rk45 = 5000000; // 5M akzeptierte Schritte max

        // Absichtlich schlechte Vorwerte setzen damit mindestens ein volles
        // check_int-Fenster integriert wird (verhindert Sofort-Konvergenz nach Warmstart)
        double n_prev  = n  * 0.5;
        double ng_prev = ng * 0.5;
        double Tg_prev = Tg * 0.5;
        double Te_prev = Te * 0.5;

        for (int i = 0; i < Nmax_rk45; ++i) {
            // k1-k6 (Dormand-Prince)
            State k1 = rhs(s,                                                          P_vol);
            State k2 = rhs(s + (dt_rk45*1.0/5.0)*k1,                                 P_vol);
            State k3 = rhs(s + dt_rk45*(3.0/40.0*k1  + 9.0/40.0*k2),                 P_vol);
            State k4 = rhs(s + dt_rk45*(44.0/45.0*k1 - 56.0/15.0*k2 + 32.0/9.0*k3), P_vol);
            State k5 = rhs(s + dt_rk45*(19372.0/6561.0*k1 - 25360.0/2187.0*k2
                                       + 64448.0/6561.0*k3 - 212.0/729.0*k4),         P_vol);
            State k6 = rhs(s + dt_rk45*(9017.0/3168.0*k1  - 355.0/33.0*k2
                                       + 46732.0/5247.0*k3 + 49.0/176.0*k4
                                       - 5103.0/18656.0*k5),                           P_vol);

            // RK4-Loesung (Ordnung 4)
            State s4 = s + dt_rk45*(35.0/384.0*k1 + 500.0/1113.0*k3
                                   + 125.0/192.0*k4 - 2187.0/6784.0*k5
                                   + 11.0/84.0*k6);
            // k7 fuer Fehlerschaetzer
            State k7 = rhs(s4, P_vol);
            // Fehlerschaetzer (Differenz RK4-RK5)
            State err = dt_rk45*(71.0/57600.0*k1 - 71.0/16695.0*k3
                                + 71.0/1920.0*k4  - 17253.0/339200.0*k5
                                + 22.0/525.0*k6   - 1.0/40.0*k7);

            // Normierter Fehler
            double sc_n  = atol + rtol*std::max(std::abs(s.n),  std::abs(s4.n));
            double sc_ng = atol + rtol*std::max(std::abs(s.ng), std::abs(s4.ng));
            double sc_Ug = atol + rtol*std::max(std::abs(s.Ug), std::abs(s4.Ug));
            double sc_Ue = atol + rtol*std::max(std::abs(s.Ue), std::abs(s4.Ue));
            double err_norm = sqrt(0.25*(
                pow(err.n /sc_n, 2) + pow(err.ng/sc_ng, 2) +
                pow(err.Ug/sc_Ug,2) + pow(err.Ue/sc_Ue, 2)));

            if (err_norm <= 1.0) {
                // Schritt akzeptieren
                s = s4;
                accepted++;

                double new_Tg = (2.0/3.0)*s.Ug/(s.ng*kB);
                double new_Te = (2.0/3.0)*s.Ue/(s.n *kB*conv);

                if (std::isnan(new_Te)||std::isnan(s.n)||s.n<=0||s.ng<=0) {
                    std::cerr << "\n[RK45] NaN bei Schritt " << accepted << std::endl;
                    *it_ptr = -1; break;
                }

                if (accepted % check_int == 0) {
                    double R_ng = std::abs((s.ng-ng_prev)/s.ng);
                    double R_n  = std::abs((s.n -n_prev) /s.n);
                    double R_Tg = std::abs((new_Tg-Tg_prev)/new_Tg);
                    double R_Te = std::abs((new_Te-Te_prev) /new_Te);
                    if (accepted % print_int == 0) {
                        std::cout << "  [" << accepted << " steps, dt="
                                  << std::scientific << std::setprecision(1) << dt_rk45
                                  << ", Te=" << std::fixed << std::setprecision(2) << new_Te << "eV"
                                  << ", n=" << std::scientific << std::setprecision(2) << s.n << "]"
                                  << std::endl;
                        std::cout.flush();
                    }
                    if (R_ng<=tol && R_n<=tol && R_Tg<=tol && R_Te<=tol) {
                        n=s.n; ng=s.ng; Tg=new_Tg; Te=new_Te;
                        energy_gas=s.Ug; energy_elec=s.Ue;
                        *it_ptr=accepted; *tol_ptr=0.25*(R_ng+R_n+R_Tg+R_Te);
                        break;
                    }
                    n_prev=s.n; ng_prev=s.ng; Tg_prev=new_Tg; Te_prev=new_Te;
                }
            }

            // Schrittweite anpassen (PI-Regler)
            double factor = 0.9 * pow(err_norm, -0.2);
            factor = std::min(std::max(factor, 0.1), 5.0);
            dt_rk45 = std::min(std::max(dt_rk45 * factor, dt_min), dt_max);
        }

        if (accepted >= Nmax_rk45)
            std::cerr << "[RK45] Nmax erreicht ohne Konvergenz!" << std::endl;
        n=s.n; ng=s.ng;
        Tg=(2.0/3.0)*s.Ug/(s.ng*kB);
        Te=(2.0/3.0)*s.Ue/(s.n *kB*conv);
        energy_gas=s.Ug; energy_elec=s.Ue;
        std::cout << endl;
        *n_ptr=n; *ng_ptr=ng; *Tg_ptr=Tg; *Te_ptr=Te;
    }
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

    std::string red   = "\033[31m";
    std::string green = "\033[32m";
    std::string cyan  = "\033[36m";
    std::string reset = "\033[0m";
    // HINWEIS: system("reset") entfernt – stoert subprocess-Ausgaben

    try {
        cout << green
             << "####################################################\n"
             << "#              Global Xenon Model                  #\n"
             << "#   Solver: "
             << (method==1 ? "Euler (RK1)          "
                 : method==2 ? "RK4                  "
                 : method==3 ? "RK45 (Dormand-Prince)"
                             : "Newton (stationaer) ")
             << "         #\n"
             << "####################################################\n"
             << reset << endl;

        // Live-Ergebnis-Datei leeren (neu starten)
        { std::ofstream live("live_results.txt", std::ios::trunc); }

        std::ofstream datei("output_kh.txt");
        if (!datei) { cerr << "Fehler: Ausgabedatei nicht oeffenbar!" << endl; return 1; }

        datei << "Method, Q0sccm, Te, Tg, n, ng, iondeg, P_RFG, P_abs, I_extr_mA, "
              << "collision_freq, R_induktiv, I_coil, epsilon_p_real, epsilon_p_imag, "
              << "u_Bohm, J_i, zeta, gamma, xi, eta, plasmafrequenz, frequency_MHz\n";

        // Startwerte (Warmstart: n, Te, Tg werden zwischen jj-Iterationen beibehalten)
        double Te = 3.75, Tg = 300.0, n = 1.0e17, ng = 1.0e19;
        double P_abs = 0.95*P_RFG, R_induktiv = 0.0, I_coil = 0.0;
        double tolerance = 0.0;
        int    it = 0;
        double P_RFG_start = P_RFG;
        bool   first_run = true;  // erstes jj: feste Startwerte, danach Warmstart

        double* n_ptr          = &n;
        double* ng_ptr         = &ng;
        double* Tg_ptr         = &Tg;
        double* Te_ptr         = &Te;
        int*    it_ptr         = &it;
        double* tol_ptr        = &tolerance;
        double* P_abs_ptr      = &P_abs;
        double* R_induktiv_ptr = &R_induktiv;
        double* I_coil_ptr     = &I_coil;

        StationaryPlasmaState stationary_prev = stationary_safe_defaults_for_q(Q0sccm_start * 4.477962312e17);
        bool stationary_have_prev = false;
        StationaryPlasmaState stationary_last_good{};
        bool stationary_have_good = false;
        double stationary_last_good_q0 = 0.0;

        // BUG-FIX: crash_of_radius-Logik war im Original invertiert.
        // Rc = Spulenradius, R = Kammerradius.
        // Spule muss AUSSERHALB liegen (Rc >= R).
        if (Rc < R) {
            cout << red << "FEHLER: Spulenradius Rc=" << Rc << " < Kammerradius R=" << R
                 << " – Spule liegt innerhalb der Kammer! Programm stoppt." << reset << endl;
            return 1;
        }

        // Aeussere Schleife: Q0sccm variieren
        for (int jj = 0; jj < jjmax; ++jj) {
            Q0sccm = Q0sccm_start + jj * Q0sccm_step;
            Q0     = Q0sccm * 4.477962312e17;
            cout << "Q0_STEP " << fixed << setprecision(4) << Q0sccm
                 << " " << (jj+1) << " " << jjmax << endl;
            debug_emit(2, "Q0_STEP", "jj=" + std::to_string(jj) + " Q0sccm=" + std::to_string(Q0sccm) + " Q0=" + std::to_string(Q0));

            if (method == 4) {
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

                StationaryPowerSolveResult ps = stationary_solve_for_target_current(guess);
                if (!stationary_power_result_valid(ps)) {
                    if (ps.hit_limit) {
                        cout << "LIMIT_HIT " << jj << " " << fixed << setprecision(4)
                             << Q0sccm << " "
                             << (std::isfinite(ps.P_trial_last) ? ps.P_trial_last : -1.0) << " "
                             << P_RFG_max << " "
                             << (std::isfinite(ps.I_mA) ? ps.I_mA : -1.0) << " "
                             << (std::isfinite(ps.err_mA) ? ps.err_mA : -1.0) << endl;
                    }

                    cout << "SOLVER_FAIL " << jj << " " << fixed << setprecision(4)
                         << Q0sccm << " " << ps.reason
                         << " P_try=" << (std::isfinite(ps.P_trial_last) ? ps.P_trial_last : -1.0)
                         << " I_mA=" << (std::isfinite(ps.I_mA) ? ps.I_mA : -1.0)
                         << " err_mA=" << (std::isfinite(ps.err_mA) ? ps.err_mA : -1.0)
                         << endl;
                    debug_emit(1, "SOLVER_FAIL", "jj=" + std::to_string(jj) + " Q0sccm=" + std::to_string(Q0sccm) + " reason=" + ps.reason, true);
                    continue;
                }

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

                double lam = lambda_i(ng);
                double I_extr = stationary_beam_current_mA(ps.state);
                double iondeg = n/ng*100.0;
                complex<double> eps_p = my_calc_eps_p(n, ng, Te);
                double cf = coll_freq(ng, Te);
                double u_Bohm = uB(Te);
                double J_i = e * Gamma_i_func(lam, Te, n);
                double zeta = R_induktiv/(R_induktiv+R_ohm);
                double pf = plasma_freq(n);

                if (std::isfinite(ps.inner_resid_norm) && ps.inner_resid_norm >= newton_tol) {
                    cout << "SOFT_ACCEPT " << fixed << setprecision(4)
                         << Q0sccm << " " << scientific << ps.inner_resid_norm << endl;
                }

                cout << "RESULT " << scientific << setprecision(4)
                     << n << " " << ng << " " << fixed << setprecision(3)
                     << Te << " " << Tg << " " << I_extr << " " << P_RFG << endl;

                datei << "Stationary" << ", " << Q0sccm << ", " << Te << ", " << Tg << ", "
                      << scientific << n << ", " << ng << ", "
                      << fixed << iondeg << ", " << P_RFG << ", "
                      << P_abs << ", " << I_extr << ", "
                      << cf << ", " << R_induktiv << ", " << I_coil << ", "
                      << real(eps_p) << ", " << imag(eps_p) << ", "
                      << u_Bohm << ", " << J_i << ", "
                      << zeta << ", " << 0.0 << ", "
                      << 0.0 << ", " << 0.0 << ", "
                      << pf << ", " << frequency/1e6 << "\n";

                {
                    std::ofstream live("live_results.txt", std::ios::app);
                    live << std::fixed << std::setprecision(6)
                         << Q0sccm << " " << P_RFG << " " << Te << " " << Tg << " "
                         << n << " " << ng << " " << I_extr << "\n";
                    live.flush();
                }

                P_RFG = P_RFG_start;
                continue;
            }

            // ng aus stationaerem Gleichgewicht ohne Plasma
            double p_null = 4 * kB * Tg0 * Q0 / (vg(Tg0) * Ag);
            ng = p_null / kB / Tg0;

            // Warmstart: ab 2. Iteration werden n, Te, Tg der vorigen Loesung wiederverwendet
            if (first_run) {
                n  = 1.0e17;
                Te = 3.75;
                Tg = 300.0;
                first_run = false;
            }
            // Energiedichten aus aktuellen Werten (konsistent initialisieren)
            double energy_gas  = 1.5 * ng * kB * Tg;
            double energy_elec = 1.5 * n  * kB * Te * conv;
            Tg = (2.0/3.0) * energy_gas  / (ng * kB);
            Te = (2.0/3.0) * energy_elec / (n  * kB * conv);

            // PID-Regler: P_RFG -> I_beam = I_soll
            P_RFG = P_RFG_start;

            // Feste Startwerte fuer jeden PID-Durchlauf speichern.
            // WICHTIG: Kein Warmstart zwischen PID-Iterationen!
            // Grund: Das System hat bei manchen P_RFG-Werten zwei stabile Zustaende
            // (Bifurkation: hoher Te-Modus ~8eV und niederer ~5eV). Ein Warmstart
            // springt zwischen diesen hin und her. Feste Startbedingungen erzwingen
            // immer denselben Ast der Loesung.
            const double n_pid0          = n;
            const double ng_pid0         = ng;
            const double Te_pid0         = Te;
            const double Tg_pid0         = Tg;
            const double energy_gas_pid0  = energy_gas;
            const double energy_elec_pid0 = energy_elec;

            // Einfacher P-Regler (kein I/D): stabiler bei nichtlinearen Systemen
            double K_p = 0.3;

            for (int iter = 0; iter < 2000; ++iter) {
                // Zustand vor jeder PID-Iteration auf feste Startwerte zuruecksetzen
                n           = n_pid0;
                ng          = ng_pid0;
                Te          = Te_pid0;
                Tg          = Tg_pid0;
                energy_gas  = energy_gas_pid0;
                energy_elec = energy_elec_pid0;

                cout << "PID_START " << iter << " " << fixed << setprecision(4) << P_RFG << endl;
                cout.flush();
                debug_emit(3, "PID_START", "iter=" + std::to_string(iter) + " P_RFG=" + std::to_string(P_RFG));

                runGM(n,ng,Te,Tg,energy_gas,energy_elec,
                      n_ptr,ng_ptr,Tg_ptr,Te_ptr,
                      it_ptr,tol_ptr,P_abs_ptr,R_induktiv_ptr,I_coil_ptr,
                      method);

                double lam    = lambda_i(ng);
                double I_extr = Ai * e * Gamma_i_func(lam,Te,n) * 1000.0;
                double error  = I_soll - I_extr;

                // Einfacher P-Schritt: P_RFG proportional zum Fehler anpassen
                P_RFG += K_p * error;
                P_RFG  = max(1.0, min(P_RFG, 200.0));

                cout << "PID_DONE " << fixed << setprecision(4) << I_extr
                     << " " << error << " " << P_RFG << " " << Te << " " << Tg << endl;
                debug_emit((it < 0) ? 1 : 3, (it < 0) ? "PID_FAIL" : "PID_DONE", "iter=" + std::to_string(iter) + " I_extr=" + std::to_string(I_extr) + " error=" + std::to_string(error) + " P_RFG=" + std::to_string(P_RFG) + " runGM_it=" + std::to_string(it), it < 0);

                if (it < 0) {
                    debug_emit(1, "PID_ABORT", "iter=" + std::to_string(iter) + " reason=runGM_failed", true);
                    break;
                }
                if (abs(error) < 0.005) {
                    cout << "CONVERGED " << iter << endl;
                    debug_emit(2, "PID_CONVERGED", "iter=" + std::to_string(iter) + " I_extr=" + std::to_string(I_extr) + " P_RFG=" + std::to_string(P_RFG));
                    break;
                }
                if (iter == 1999)
                    cout << "PID_MAXITER" << endl;
            }

            // Abgeleitete Groessen
            double lam_f     = lambda_i(ng);
            double Gi_f      = Gamma_i_func(lam_f,Te,n);
            double vg_f      = vg(Tg);
            double Gn_f      = 0.25*ng*vg_f;
            double T_i       = Thrust(Gi_f);
            double T_n       = Gn_f * M * vg_f * Ag;
            double u_Bohm    = uB(Te);
            double J_i       = e * Gi_f;
            double I_extr    = Ai * J_i;
            double iondeg    = n/ng*100.0;
            double zeta      = R_induktiv/(R_induktiv+R_ohm);
            double v_extr    = sqrt(2*e*Vgrid/M);
            double pow_i     = 0.5*M*v_extr*v_extr*Gi_f*Ai;
            double pow_n     = 0.5*M*vg_f*vg_f*0.25*ng*vg_f*Ag;
            double gamma_eff = (pow_i+pow_n)/(pow_i+pow_n+P_RFG);
            double xi        = 1000.0*(T_i+T_n)/P_RFG;
            double eta       = Gi_f*Ai/Q0;
            complex<double> eps_p = my_calc_eps_p(n,ng,Te);
            double pf        = plasma_freq(n);
            double cf        = coll_freq(ng,Te);

            cout << "RESULT " << scientific << setprecision(4)
                 << n << " " << ng << " " << fixed << setprecision(3)
                 << Te << " " << Tg << " " << I_extr*1000 << " " << P_RFG << endl;

            string meth_str = (method==1) ? "Euler" : (method==2 ? "RK4" : (method==3 ? "RK45" : "Stationary"));
            datei << meth_str  << ", " << Q0sccm    << ", " << Te         << ", " << Tg    << ", "
                  << scientific << n   << ", "       << ng                << ", "
                  << fixed      << iondeg << ", "    << P_RFG             << ", "
                  << P_abs      << ", " << I_extr*1000 << ", "
                  << cf         << ", " << R_induktiv  << ", " << I_coil  << ", "
                  << real(eps_p)<< ", " << imag(eps_p) << ", "
                  << u_Bohm     << ", " << J_i          << ", "
                  << zeta       << ", " << gamma_eff    << ", "
                  << xi         << ", " << eta          << ", "
                  << pf         << ", " << frequency/1e6 << "\n";

            // Live-Ergebnis fuer Python-Plot schreiben (gefiltert)
            {
                static bool first_point = true;
                static double P_prev = -1.0;

                bool valid = true;

                // physikalische Grenzen
                if (P_RFG <= 0.0 || P_RFG > 1e4) valid = false;

                // erster Punkt verwerfen
                if (first_point) {
                    valid = false;
                    first_point = false;
                }

                // Sprungfilter
                if (P_prev > 0.0) {
                    if (fabs(P_RFG - P_prev) > 0.5 * P_prev) valid = false;
                }

                if (valid) {
                    std::ofstream live("live_results.txt", std::ios::app);
                    live << std::fixed << std::setprecision(6)
                         << Q0sccm << " " << P_RFG << " " << Te << " " << Tg << " "
                         << n << " " << ng << " " << I_extr*1000 << "\n";
                    live.flush();
                    P_prev = P_RFG;
                }
            }

            P_RFG = P_RFG_start;
        }

    } catch (const exception& ex) {
        cerr << "Fehler: " << ex.what() << endl;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    cout << "\nDer Code hat " << fixed << setprecision(3)
         << std::chrono::duration<double>(end_time - start_time).count()
         << " Sekunden gebraucht." << endl;
    return 0;
}
