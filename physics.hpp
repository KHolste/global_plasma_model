// physics.hpp -- Plasmaphysik, RF-Kopplung, Residuen, abgeleitete Groessen.
//
// Hot-Path-Funktionen sind inline. Groessere Funktionen in physics.cpp.
#ifndef PHYSICS_HPP
#define PHYSICS_HPP

#include "sim_context.hpp"
#include "rates.hpp"
#include "beam_extraction_cpp.hpp"
#include <cmath>
#include <complex>
#include <array>

// ═══ Kleine Hilfsfunktionen (inline, hot path) ═══════════
using namespace PhysConst;

inline double vg(double M, double Tg)    { return std::sqrt(8*kB*Tg/(pi*M)); }
inline double vi(double M, double Ti)    { return std::sqrt(8*kB*Ti/(pi*M)); }
inline double uB(double M, double Te)    { return std::sqrt(kB*Te*conv/M); }
inline double lambda_i(double ng, double sigma) { return 1.0/(ng*sigma); }
inline double plasma_freq(double n)      { return std::sqrt(n*e*e/(me*epsilon0)); }
inline double coll_freq(const SimContext& ctx, double ng, double Te) { return Kel(ctx, Te)*ng; }

inline double Aeff(const Thruster& t, double lambda) {
    double hL = 0.86 * std::pow(3 + t.L/(2*lambda), -0.5);
    double hR = 0.80 * std::pow(4 + t.R/lambda,     -0.5);
    return 2*hR*pi*t.R*t.L + 2*hL*pi*t.R*t.R;
}

inline double Aeff1(const Thruster& t, double lambda) {
    double hL = 0.86 * std::pow(3 + t.L/(2*lambda), -0.5);
    double hR = 0.80 * std::pow(4 + t.R/lambda,     -0.5);
    return 2*hR*pi*t.R*t.L + (2 - t.betai)*hL*pi*t.R*t.R;
}

inline double Gamma_i_func(const Thruster& t, double lambda, double Te, double n, double M) {
    double hL = 0.86 * std::pow(3 + t.L/(2*lambda), -0.5);
    return hL * n * uB(M, Te);
}

inline double beam_current_mA(const SimContext& ctx, const PlasmaState& s) {
    // Volles Extraktionsmodell: Bohm + Child-Langmuir + eta_opt
    return beam_current_mA_full(ctx, s);
}

template <typename T> inline int signum(T val) { return (T(0)<val)-(val<T(0)); }

// ═══ Zustandspruefungen ════════════════════════════════════

constexpr double iondeg_max = 0.95;

inline bool state_finite_positive(const PlasmaState& s) {
    return std::isfinite(s.n) && std::isfinite(s.ng) && std::isfinite(s.Te) && std::isfinite(s.Tg)
        && s.n > 0 && s.ng > 0 && s.Te > 0 && s.Tg > 0;
}

inline bool state_in_bounds(const SolverParams& sp, const PlasmaState& s) {
    if (!state_finite_positive(s)) return false;
    if (s.n < sp.n_min || s.n > sp.n_max) return false;
    if (s.ng < sp.ng_min || s.ng > sp.ng_max) return false;
    if (s.Te < sp.Te_min || s.Te > sp.Te_max) return false;
    if (s.Tg < sp.Tg_min || s.Tg > sp.Tg_max) return false;
    double id = s.n / std::max(1.0, s.ng);
    return std::isfinite(id) && id >= 0 && id <= iondeg_max;
}

inline PlasmaState safe_defaults(const SimContext& ctx) {
    PlasmaState s;
    double p0 = 4*kB*Tg0*ctx.thruster.Q0 / (vg(ctx.gas.M, Tg0)*ctx.thruster.Ag);
    s.ng = std::max(1e16, p0/(kB*Tg0));
    s.n = 1e17; s.Te = 3.75; s.Tg = 300.0;
    return s;
}

// ═══ Groessere Funktionen (deklariert, impl in physics.cpp) ═

RFState compute_rf(const SimContext& ctx, double n, double ng, double Te, double P_RFG);
// Gleiche Rechnung, aber mit von aussen bestimmter Stossfrequenz. Der
// generische Weg bildet sie aus den Reaktionen des Chemiepakets; der feste
// Weg laesst sie aus dem einen elastischen Ratenkoeffizienten bilden.
RFState compute_rf_nu(const SimContext& ctx, double n, double nu_m, double P_RFG);
std::array<double,4> residual_raw(const SimContext& ctx, const PlasmaState& s, double P_RFG, RFState* rf_out = nullptr);
std::array<double,4> residual_scaled(const SimContext& ctx, const PlasmaState& s, double P_RFG, RFState* rf_out = nullptr);
DerivedQuantities compute_derived(const SimContext& ctx, double n, double ng, double Te, double Tg, double R_ind, double I_coil, double P_abs);

// Gleiche Groessen, aber mit bereits gerechneter Extraktion. Der generische
// Weg reicht hier die Rechnung ueber alle Ionensorten herein; der feste Weg
// laesst sie aus dem Einsortenzustand bilden. Strahlstrom, Schub und
// Wirkungsgrade stammen so in beiden Faellen aus derselben Quelle.
DerivedQuantities compute_derived_from(const SimContext& ctx, const ExtractionResult& ex,
                                        double n, double ng, double Te, double Tg,
                                        double R_ind, double I_coil, double P_abs);

void emit_csv_row(std::ofstream& datei, const std::string& method, double Q0sccm,
                   double n, double ng, double Te, double Tg,
                   double P_RFG, double P_abs, double R_ind, double I_coil,
                   const SimContext& ctx, const DerivedQuantities& d);

#endif
