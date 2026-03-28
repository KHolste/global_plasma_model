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

    double A    = 2 * pi * R * R + pi * R * L;
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
    double newton_tol        = 1e-7;
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
        A        = 2 * pi * R * R + pi * R * L;
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
    double hR = 0.80 * pow(4 + L/lambda,     -0.5);
    return 2*hR*pi*R*L + 2*hL*pi*R*R;
}
double Aeff1(double lambda) {
    double hL = 0.86 * pow(3 + L/(2*lambda), -0.5);
    double hR = 0.80 * pow(4 + L/lambda,     -0.5);
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

    if (!std::isfinite(denom_c.real()) || !std::isfinite(denom_c.imag()) || std::abs(denom_c) == 0.0) return out;

    complex<double> result = (1i * kR * j1) / denom_c;
    double R_ind = 2*pi*Nw*Nw / (epsilon0*L*omega) * result.real();
    if (!std::isfinite(R_ind) || R_ind <= 0.0) return out;

    double Ic = sqrt(2.0 * P_RFG_local / (R_ind + R_ohm));
    double P_abs = 0.5 * R_ind * Ic * Ic;
    if (!std::isfinite(Ic) || !std::isfinite(P_abs) || Ic < 0.0 || P_abs < 0.0) return out;

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
    double scale1 = max(1e10, s.n);
    double scale2 = max(1e14, s.ng);
    double scale3 = max(1.0, fabs(rf.P_abs / V));
    double scale4 = max(1.0, fabs(kappa*(s.Tg-Tg0)/lambda_0 * A/V));
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


static StationarySolveResult stationary_solve_newton(double P_RFG_local, const StationaryPlasmaState& initial) {
    StationarySolveResult out;
    out.state = initial;
    if (!stationary_state_in_bounds(initial)) {
        out.reason = "invalid initial state";
        return out;
    }

    std::array<double,4> x = {log(initial.n), log(initial.ng), log(initial.Te), log(initial.Tg)};

    for (int iter = 0; iter < newton_max_iter; ++iter) {
        StationaryPlasmaState s{exp(x[0]), exp(x[1]), exp(x[2]), exp(x[3])};
        StationaryRFState rf;
        auto r = stationary_residual_scaled(s, P_RFG_local, &rf);
        double merit0 = stationary_merit(r, s);

        std::cout << "NEWTON_IT " << iter << " "
                  << std::fixed << std::setprecision(4) << P_RFG_local << " "
                  << std::scientific << std::setprecision(6)
                  << merit0 << " " << s.n << " " << s.ng << " "
                  << std::fixed << std::setprecision(6)
                  << s.Te << " " << s.Tg << std::endl;

        if (!std::isfinite(merit0) || !rf.valid) {
            out.reason = "invalid residual/RF";
            out.state = s;
            out.rf = rf;
            out.iterations = iter;
            return out;
        }

        if (merit0 < newton_tol) {
            out.converged = true;
            out.state = s;
            out.rf = rf;
            out.iterations = iter;
            out.resid_norm = merit0;
            out.reason = "ok";
            debug_emit(2, "NEWTON_OK", "iter=" + std::to_string(iter) + " merit=" + std::to_string(out.resid_norm));
            return out;
        }

        double J[4][4];
        for (int j = 0; j < 4; ++j) {
            std::array<double,4> xp = x, xm = x;
            double h = newton_fd_eps * std::max(1.0, std::fabs(x[j]));
            xp[j] += h;
            xm[j] -= h;

            StationaryPlasmaState sp{exp(xp[0]), exp(xp[1]), exp(xp[2]), exp(xp[3])};
            StationaryPlasmaState sm{exp(xm[0]), exp(xm[1]), exp(xm[2]), exp(xm[3])};
            auto rp = stationary_residual_scaled(sp, P_RFG_local, nullptr);
            auto rm = stationary_residual_scaled(sm, P_RFG_local, nullptr);

            for (int i = 0; i < 4; ++i) {
                if (!std::isfinite(rp[i]) || !std::isfinite(rm[i])) {
                    out.reason = "jacobian nan";
                    out.state = s;
                    out.rf = rf;
                    out.iterations = iter;
                    out.resid_norm = merit0;
                    return out;
                }
                J[i][j] = (rp[i] - rm[i]) / (2.0*h);
            }
        }

        double minus_r[4] = {-r[0], -r[1], -r[2], -r[3]};
        double dx[4];
        if (!stationary_solve_linear_4x4(J, minus_r, dx)) {
            out.reason = "linear solve failed";
            out.state = s;
            out.rf = rf;
            out.iterations = iter;
            out.resid_norm = merit0;
            return out;
        }

        for (int i = 0; i < 4; ++i) {
            if (!std::isfinite(dx[i])) {
                out.reason = "dx nan";
                out.state = s;
                out.rf = rf;
                out.iterations = iter;
                out.resid_norm = merit0;
                return out;
            }
            dx[i] = std::max(-newton_max_log_step, std::min(newton_max_log_step, dx[i]));
        }

        bool accepted = false;
        std::array<double,4> x_trial{};
        double best_trial_merit = std::numeric_limits<double>::infinity();
        double best_trial_alpha = 0.0;
        int best_trial_in_bounds = 0;
        int best_trial_finite = 0;
        for (double alpha : {1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625}) {
            for (int i = 0; i < 4; ++i) x_trial[i] = x[i] + alpha * dx[i];
            StationaryPlasmaState st{exp(x_trial[0]), exp(x_trial[1]), exp(x_trial[2]), exp(x_trial[3])};
            if (!stationary_state_in_bounds(st)) {
                DBG(3, "LINESEARCH_TRY iter=" << iter << " alpha=" << alpha << " status=out_of_bounds");
                continue;
            }
            best_trial_in_bounds++;
            auto rt = stationary_residual_scaled(st, P_RFG_local, nullptr);
            double mt = stationary_merit(rt, st);
            if (std::isfinite(mt)) {
                best_trial_finite++;
                if (mt < best_trial_merit) {
                    best_trial_merit = mt;
                    best_trial_alpha = alpha;
                }
                DBG(3, "LINESEARCH_TRY iter=" << iter << " alpha=" << alpha << " merit=" << mt << " merit0=" << merit0);
            } else {
                DBG(3, "LINESEARCH_TRY iter=" << iter << " alpha=" << alpha << " status=nonfinite");
            }
            if (std::isfinite(mt) && mt < merit0) {
                DBG(2, "LINESEARCH_ACCEPT iter=" << iter << " alpha=" << alpha << " merit0=" << merit0 << " merit=" << mt);
                x = x_trial;
                accepted = true;
                break;
            }
        }

        if (!accepted) {
            bool fallback_ok = false;
            int fallback_tries = 0;
            double best_fallback_merit = std::numeric_limits<double>::infinity();
            int best_fallback_j = -1;
            double best_fallback_step = 0.0;
            for (int j = 0; j < 4 && !fallback_ok; ++j) {
                for (double step : {0.1, -0.1, 0.03, -0.03}) {
                    fallback_tries++;
                    x_trial = x;
                    x_trial[j] += step;
                    StationaryPlasmaState st{exp(x_trial[0]), exp(x_trial[1]), exp(x_trial[2]), exp(x_trial[3])};
                    if (!stationary_state_in_bounds(st)) {
                        DBG(3, "FALLBACK_TRY iter=" << iter << " j=" << j << " step=" << step << " status=out_of_bounds");
                        continue;
                    }
                    auto rt = stationary_residual_scaled(st, P_RFG_local, nullptr);
                    double mt = stationary_merit(rt, st);
                    if (std::isfinite(mt) && mt < best_fallback_merit) {
                        best_fallback_merit = mt;
                        best_fallback_j = j;
                        best_fallback_step = step;
                    }
                    DBG(3, "FALLBACK_TRY iter=" << iter << " j=" << j << " step=" << step << " merit=" << mt << " merit0=" << merit0);
                    if (std::isfinite(mt) && mt < merit0) {
                        DBG(2, "FALLBACK_ACCEPT iter=" << iter << " j=" << j << " step=" << step << " merit0=" << merit0 << " merit=" << mt);
                        x = x_trial;
                        fallback_ok = true;
                        break;
                    }
                }
            }
            if (!fallback_ok) {
                DBG(1, "LINESEARCH_FAIL_DETAIL iter=" << iter
                    << " merit0=" << merit0
                    << " best_alpha=" << best_trial_alpha
                    << " best_trial_merit=" << best_trial_merit
                    << " in_bounds=" << best_trial_in_bounds
                    << " finite=" << best_trial_finite
                    << " best_fallback_j=" << best_fallback_j
                    << " best_fallback_step=" << best_fallback_step
                    << " best_fallback_merit=" << best_fallback_merit
                    << " dx=[" << dx[0] << "," << dx[1] << "," << dx[2] << "," << dx[3] << "]");
                out.reason = "line search failed";
                out.state = s;
                out.rf = rf;
                out.iterations = iter;
                out.resid_norm = merit0;
                return out;
            }
        }
    }

    StationaryPlasmaState s{exp(x[0]), exp(x[1]), exp(x[2]), exp(x[3])};
    StationaryRFState rf;
    auto r = stationary_residual_scaled(s, P_RFG_local, &rf);
    out.converged = stationary_merit(r, s) < newton_tol && rf.valid && stationary_state_in_bounds(s);
    out.state = s;
    out.rf = rf;
    out.iterations = newton_max_iter;
    out.resid_norm = stationary_merit(r, s);
    out.reason = out.converged ? "ok" : "newton max iter";
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

    // 1) Erst direkter Newton-Versuch
    best = stationary_solve_newton(P_RFG_local, initial);
    if (best.converged) return best;
    if (stationary_soft_accept(best, merit_initial)) {
        best.converged = true;
        best.reason = "soft-ok-direct";
        return best;
    }

    // 2) Pseudo-transient continuation in Log-Variablen
    std::array<double,4> x = x0;
    StationaryPlasmaState s = initial;
    StationaryRFState rf;
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
            StationarySolveResult polished = stationary_solve_newton(P_RFG_local, s);
            if (polished.converged) return polished;
            if (stationary_soft_accept(polished, merit_initial)) {
                polished.converged = true;
                polished.reason = "soft-ok-polished";
                return polished;
            }
            if (polished.resid_norm < best_local.resid_norm) best_local = polished;
        }
    }

    // 3) Letzter Newton-Finish vom besten gefundenen Zustand
    StationarySolveResult final_try = stationary_solve_newton(P_RFG_local, best_local.state);
    if (final_try.converged) return final_try;
    if (stationary_soft_accept(final_try, merit_initial)) {
        final_try.converged = true;
        final_try.reason = "soft-ok-final";
        return final_try;
    }
    if (stationary_soft_accept(best_local, merit_initial)) {
        best_local.converged = true;
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
           ps.P_RFG_sol <= P_RFG_max + 1e-9 &&
           std::isfinite(ps.I_mA);
}

static StationaryPowerSolveResult stationary_solve_for_target_current(const StationaryPlasmaState& guess) {
    StationaryPowerSolveResult out;

    debug_emit(2, "TARGET_CURRENT_BEGIN", "Q0=" + std::to_string(Q0) + " I_soll=" + std::to_string(I_soll) + " P_RFG_seed=" + std::to_string(P_RFG));

    double p_center = std::max(power_min, std::min(P_RFG, P_RFG_max));
    double p_lo = std::max(power_min, std::min(p_center * 0.6, P_RFG_max));
    double p_hi = std::max(p_lo + 0.5, std::min(P_RFG_max, std::max(p_center, p_lo + 2.0)));

    std::vector<StationaryPlasmaState> starts;
    StationaryPlasmaState safe = stationary_safe_defaults_for_q(Q0);
    starts.push_back(guess);
    starts.push_back(safe);
    starts.push_back(StationaryPlasmaState{
        std::sqrt(std::max(1.0, guess.n * safe.n)),
        std::sqrt(std::max(1.0, guess.ng * safe.ng)),
        0.5*(guess.Te + safe.Te),
        0.5*(guess.Tg + safe.Tg)
    });

    StationarySolveResult s_center = stationary_solve_power_robust(p_center, starts);
    bool center_ok = s_center.converged;
    if (!center_ok) {
        out.reason = "center power solve failed: " + s_center.reason;
        debug_emit(1, "TARGET_CURRENT_FAIL", "stage=center P_try=" + std::to_string(p_center) + " reason=" + out.reason + " resid=" + std::to_string(s_center.resid_norm), true);
        out.P_trial_last = p_center;
        out.state = s_center.state;
        out.rf = s_center.rf;
        out.I_mA = stationary_state_finite_positive(s_center.state) ? stationary_beam_current_mA(s_center.state)
                                                                    : std::numeric_limits<double>::quiet_NaN();
        out.err_mA = std::isfinite(out.I_mA) ? (I_soll - out.I_mA) : std::numeric_limits<double>::quiet_NaN();
        return out;
    }

    StationarySolveResult s_lo = stationary_solve_power_robust(p_lo, std::vector<StationaryPlasmaState>{s_center.state, stationary_safe_defaults_for_q(Q0)});
    if (!s_lo.converged) s_lo = s_center;
    double f_lo = stationary_beam_current_mA(s_lo.state) - I_soll;

    StationarySolveResult s_hi = stationary_solve_power_robust(p_hi, std::vector<StationaryPlasmaState>{s_center.state, stationary_safe_defaults_for_q(Q0)});
    if (!s_hi.converged) {
        bool hi_ok = false;
        for (double ph : {std::min(P_RFG_max, p_center + 5.0), std::min(P_RFG_max, p_center + 10.0), P_RFG_max}) {
            s_hi = stationary_solve_power_robust(ph, std::vector<StationaryPlasmaState>{s_center.state, stationary_safe_defaults_for_q(Q0)});
            if (s_hi.converged) { p_hi = ph; hi_ok = true; break; }
        }
        if (!hi_ok) {
            out.reason = "high power initial solve failed";
            out.P_trial_last = p_hi;
            return out;
        }
    }
    double f_hi = stationary_beam_current_mA(s_hi.state) - I_soll;

    {
        double I_lo = stationary_beam_current_mA(s_lo.state);
        double I_hi = stationary_beam_current_mA(s_hi.state);
        std::ostringstream os;
        os << std::fixed << std::setprecision(6)
           << "p_lo=" << p_lo << " p_hi=" << p_hi
           << " I_lo=" << I_lo << " I_hi=" << I_hi
           << " dI=" << (I_hi - I_lo)
           << " f_lo=" << f_lo << " f_hi=" << f_hi;
        debug_emit(2, "POWER_RESPONSE", os.str());
        if (std::isfinite(I_lo) && std::isfinite(I_hi) && std::fabs(I_hi - I_lo) < 1e-3 && std::fabs(p_hi - p_lo) > 1.0) {
            debug_emit(1, "POWER_RESPONSE_FLAT", os.str() + " reason=stationary solver returns nearly identical current for different power", true);
        }
    }

    int bracket_expand = 0;
    while (f_lo * f_hi > 0.0 && p_hi < P_RFG_max && bracket_expand < 12) {
        double new_hi = std::min(P_RFG_max, std::max(p_hi * 1.35, p_hi + 5.0));
        std::cout << "POWER_BRACKET " << std::fixed << std::setprecision(4)
                  << p_lo << " " << new_hi << " " << f_lo << " " << f_hi << std::endl;
        StationarySolveResult s_try = stationary_solve_power_robust(new_hi, std::vector<StationaryPlasmaState>{s_hi.state, s_center.state, stationary_safe_defaults_for_q(Q0)});
        if (s_try.converged) {
            p_hi = new_hi;
            s_hi = s_try;
            f_hi = stationary_beam_current_mA(s_hi.state) - I_soll;
        } else {
            p_hi = new_hi;
        }
        bracket_expand++;
    }

    if (f_lo * f_hi > 0.0) {
        out.hit_limit = true;
        out.reason = "no bracket within power limit (I(P) stayed on one side or was nearly flat)";
        out.P_trial_last = p_hi;
        out.state = s_hi.state;
        out.rf = s_hi.rf;
        out.I_mA = stationary_state_finite_positive(s_hi.state) ? stationary_beam_current_mA(s_hi.state)
                                                                : std::numeric_limits<double>::quiet_NaN();
        out.err_mA = std::isfinite(out.I_mA) ? (I_soll - out.I_mA) : std::numeric_limits<double>::quiet_NaN();
        return out;
    }

    StationaryPlasmaState last_guess = s_center.state;
    for (int iter = 0; iter < power_max_iter; ++iter) {
        double denom = (f_hi - f_lo);
        double p_mid = (std::fabs(denom) > 1e-10) ? p_hi - f_hi * (p_hi - p_lo) / denom : 0.5 * (p_lo + p_hi);
        if (!(p_mid > p_lo && p_mid < p_hi) || !std::isfinite(p_mid)) p_mid = 0.5 * (p_lo + p_hi);

        std::cout << "PID_START " << iter << " " << std::fixed << std::setprecision(4) << p_mid << std::endl;

        StationarySolveResult s_mid = stationary_solve_power_robust(p_mid, std::vector<StationaryPlasmaState>{last_guess, s_center.state, stationary_safe_defaults_for_q(Q0)});
        if (!s_mid.converged) {
            StationaryPlasmaState safe2 = stationary_safe_defaults_for_q(Q0);
            p_mid = 0.5 * (p_lo + p_hi);
            s_mid = stationary_solve_power_robust(p_mid, std::vector<StationaryPlasmaState>{safe2, s_center.state});
            if (!s_mid.converged) {
                out.reason = "mid solve failed";
                out.P_trial_last = p_mid;
                return out;
            }
        }

        double I_mA = stationary_beam_current_mA(s_mid.state);
        double error = I_soll - I_mA;

        std::cout << "PID_DONE " << std::fixed << std::setprecision(4)
                  << I_mA << " " << error << " " << p_mid << " "
                  << s_mid.state.Te << " " << s_mid.state.Tg << std::endl;

        last_guess = s_mid.state;
        out.iterations = iter;
        out.state = s_mid.state;
        out.rf = s_mid.rf;
        out.P_RFG_sol = p_mid;
        out.P_trial_last = p_mid;
        out.I_mA = I_mA;
        out.err_mA = error;

        if (std::fabs(error) < power_tol_mA) {
            out.converged = true;
            out.reason = "ok";
            debug_emit(2, "TARGET_CURRENT_OK", "iter=" + std::to_string(iter) + " P_sol=" + std::to_string(out.P_RFG_sol) + " I_mA=" + std::to_string(out.I_mA) + " err_mA=" + std::to_string(out.err_mA));
            std::cout << "CONVERGED " << iter << std::endl;
            return out;
        }

        double f_mid = I_mA - I_soll;
        if (f_lo * f_mid <= 0.0) {
            p_hi = p_mid;
            f_hi = f_mid;
        } else {
            p_lo = p_mid;
            f_lo = f_mid;
        }
    }

    out.reason = "power max iter";
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
                StationaryPlasmaState guess = stationary_have_prev ? stationary_prev : stationary_safe_defaults_for_q(Q0);
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
