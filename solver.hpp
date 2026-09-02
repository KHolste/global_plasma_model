// solver.hpp -- Numerische Solver: LM, PTC, Multi-Start, I-fix, SC-Modus.
//
// Alle Funktionen nehmen SimContext per Referenz. Kein globaler Zustand.
#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "physics.hpp"

// Merit-Funktion (inline, hot path)
inline double merit(const std::array<double,4>& r, const SolverParams& sp, const PlasmaState& s) {
    if (!state_in_bounds(sp, s)) return 1e30;
    double m = 0;
    for (double v : r) { if (!std::isfinite(v)) return 1e30; m = std::max(m, std::fabs(v)); }
    return m;
}

// Solver-Algorithmen (impl in solver.cpp)
SolveResult solve_lm(const SimContext& ctx, double P_RFG, const PlasmaState& initial);
SolveResult solve_ptc_then_lm(const SimContext& ctx, double P_RFG, const PlasmaState& initial);
PowerSolveResult solve_at_fixed_power(const SimContext& ctx, double P_RFG, const PlasmaState& guess);
PowerSolveResult solve_for_target_current(const SimContext& ctx, const PlasmaState& guess);

bool power_result_valid(const SimContext& ctx, const PowerSolveResult& ps);

#endif
