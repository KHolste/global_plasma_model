// rates.hpp -- Ratenkoeffizienten und Tabellenzugriff.
//
// Hot-Path-Funktionen (Kiz, Kel, Kex) sind inline fuer Performance.
// Tabellen-I/O lebt in rates.cpp.
#ifndef RATES_HPP
#define RATES_HPP

#include "sim_context.hpp"
#include <cmath>
#include <string>
#include <iostream>

// Cold path: Tabellen laden (implementiert in rates.cpp)
bool load_kiz_table(RateConfig& r, const std::string& path);
bool load_kel_table(RateConfig& r, const std::string& path);
bool load_kex_table(RateConfig& r, const std::string& path);

// ═══ Inline-Interpolation (Hot Path) ═══════════════════════

namespace detail {
    template<typename E, typename F>
    inline double interp(const std::vector<E>& tbl, double Te, F getter) {
        if (tbl.empty() || Te <= 0.0) return 0.0;
        if (Te <= tbl.front().Te_eV) return getter(tbl.front());
        if (Te >= tbl.back().Te_eV)  return getter(tbl.back());
        size_t lo = 0, hi = tbl.size() - 1;
        while (hi - lo > 1) { size_t m = (lo+hi)/2; if (tbl[m].Te_eV <= Te) lo = m; else hi = m; }
        double t = (Te - tbl[lo].Te_eV) / (tbl[hi].Te_eV - tbl[lo].Te_eV);
        return getter(tbl[lo]) * (1.0-t) + getter(tbl[hi]) * t;
    }
}

// ═══ Ratenkoeffizienten ════════════════════════════════════

inline double Kiz(const SimContext& ctx, double Te) {
    if (ctx.rates.ionization_model == 1 && !ctx.rates.kiz.empty())
        return detail::interp(ctx.rates.kiz, Te, [](const KizEntry& e){ return e.Kiz; });
    using namespace PhysConst;
    double TeV = kB * Te * conv / e;
    double K1 = 6.73e-15 * std::sqrt(TeV) * (3.97 + 0.643*TeV - 0.0368*TeV*TeV) * std::exp(-12.127/TeV);
    double K2 = 6.73e-15 * std::sqrt(TeV) * (-0.0001031*TeV*TeV + 6.386*std::exp(-12.127/TeV));
    return 0.5 * (K1 + K2);
}

inline double Kel(const SimContext& ctx, double Te) {
    double val;
    if (ctx.rates.elastic_model == 1 && !ctx.rates.kel.empty())
        val = detail::interp(ctx.rates.kel, Te, [](const KelEntry& e){ return e.Kel; });
    else
        val = ctx.rates.kel_constant;
    if (val <= 0.0) { val = 1e-15; }
    return val;
}

inline double Kex(const SimContext& ctx, double Te) {
    using namespace PhysConst;
    return ctx.rates.Kex_scale * 1.2921e-13 * std::exp(-e * 11.6 / (kB * Te * conv));
}

inline double Pexc_coeff(const SimContext& ctx, double Te) {
    if (ctx.rates.excitation_model == 1 && !ctx.rates.kex.empty())
        return detail::interp(ctx.rates.kex, Te, [](const KexEntry& e){ return e.Pexc_coeff; });
    return 0.0;
}

#endif
