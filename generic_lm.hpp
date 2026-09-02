// generic_lm.hpp -- N-dimensionaler Levenberg-Marquardt Solver.
//
// Arbeitet mit generischem ChemSystem und assemble_residual().
// Zustandsvektor hat variable Laenge: N_species + 2 (Te, Tg).
// Ersetzt den bisherigen festen 4D-Solver fuer den produktiven Pfad.
//
// Architektur:
//   GenSolveResult   -- Solver-Ergebnis mit generischem Zustand
//   gen_solve_lm()   -- LM-Solver fuer beliebige Zustandsdimension
//   gen_solve_sc()   -- SC-Modus (feste Leistung)
//   gen_solve_ifix() -- I-fix-Modus (Zielsrom-Bisection)
//
// Der Xenon-Fall (N=4: Xe, Xe+, Te, Tg) ist der erste Nutzer.
#ifndef GENERIC_LM_HPP
#define GENERIC_LM_HPP

#include "chem_system.hpp"
#include "physics.hpp"
#include "beam_extraction_cpp.hpp"
#include "sim_logging.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string>
#include <limits>

// ═════════════════════════════════════════════════════════════
// Zustandslayout
// ═════════════════════════════════════════════════════════════
struct StateLayout {
    int N;          // Gesamtdimension des Zustandsvektors
    int n_species;  // Anzahl Spezies (ohne Te, Tg)
    int Te_idx;     // Index von Te im Vektor
    int Tg_idx;     // Index von Tg im Vektor

    StateLayout() : N(0), n_species(0), Te_idx(0), Tg_idx(0) {}
    explicit StateLayout(const ChemSystem& sys)
        : N(sys.state_size()), n_species((int)sys.species.size()),
          Te_idx(sys.Te_idx()), Tg_idx(sys.Tg_idx()) {}
};

// ═════════════════════════════════════════════════════════════
// Solver-Ergebnis (generisch)
// ═════════════════════════════════════════════════════════════
struct GenSolveResult {
    bool converged = false;
    bool soft_ok = false;
    std::vector<double> state;     // Generischer Zustandsvektor
    RFState rf;
    int iterations = 0;
    double resid_norm = 1e30;
    std::string reason;

    // Convenience: PlasmaState fuer Xenon-Kompatibilitaet extrahieren
    PlasmaState to_plasma_state(const ChemSystem& sys) const {
        PlasmaState ps{};
        if ((int)state.size() >= sys.state_size()) {
            // Erste positive Ionenspezies = n (Elektronendichte)
            for (int i = 0; i < (int)sys.species.size(); ++i) {
                if (sys.species[i].is_positive_ion()) { ps.n = state[i]; break; }
            }
            // Feedstock = ng
            for (int i = 0; i < (int)sys.species.size(); ++i) {
                if (sys.species[i].is_feedstock) { ps.ng = state[i]; break; }
            }
            ps.Te = state[sys.Te_idx()];
            ps.Tg = state[sys.Tg_idx()];
        }
        return ps;
    }
};

// ═════════════════════════════════════════════════════════════
// NxN Gauss mit Pivotisierung (generisch)
// ═════════════════════════════════════════════════════════════
static bool solve_NxN(int N, std::vector<std::vector<double>>& A,
                       std::vector<double>& b, std::vector<double>& x) {
    // Augmentierte Matrix [A|b]
    std::vector<std::vector<double>> M(N, std::vector<double>(N+1));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) M[i][j] = A[i][j];
        M[i][N] = b[i];
    }
    for (int c = 0; c < N; ++c) {
        int piv = c;
        for (int r = c+1; r < N; ++r)
            if (std::fabs(M[r][c]) > std::fabs(M[piv][c])) piv = r;
        if (std::fabs(M[piv][c]) < 1e-20 || !std::isfinite(M[piv][c])) return false;
        if (piv != c) std::swap(M[piv], M[c]);
        double d = M[c][c];
        for (int j = c; j <= N; ++j) M[c][j] /= d;
        for (int r = 0; r < N; ++r) {
            if (r == c) continue;
            double f = M[r][c];
            for (int j = c; j <= N; ++j) M[r][j] -= f * M[c][j];
        }
    }
    x.resize(N);
    for (int i = 0; i < N; ++i) x[i] = M[i][N];
    return true;
}

// ═════════════════════════════════════════════════════════════
// Merit-Funktion (generisch)
// ═════════════════════════════════════════════════════════════
static double gen_merit(const std::vector<double>& F) {
    double m = 0;
    for (double v : F) {
        if (!std::isfinite(v)) return 1e30;
        m = std::max(m, std::fabs(v));
    }
    return m;
}

// ═════════════════════════════════════════════════════════════
// Zustandsgrenzen-Check (generisch)
// ═════════════════════════════════════════════════════════════
static bool gen_state_valid(const std::vector<double>& state, const StateLayout& lay,
                             const SolverParams& sp) {
    for (int i = 0; i < lay.n_species; ++i) {
        if (!std::isfinite(state[i]) || state[i] <= 0) return false;
    }
    double Te = state[lay.Te_idx], Tg = state[lay.Tg_idx];
    if (Te < sp.Te_min || Te > sp.Te_max) return false;
    if (Tg < sp.Tg_min || Tg > sp.Tg_max) return false;
    return true;
}

// ═════════════════════════════════════════════════════════════
// Residual-Wrapper: ChemSystem → skalierter Residualvektor
// ═════════════════════════════════════════════════════════════
static std::vector<double> gen_residual_scaled(
    const ChemSystem& sys, const std::vector<double>& state,
    double P_RFG, const SimContext& ctx, RFState* rf_out)
{
    const auto& t = ctx.thruster;
    int N = sys.state_size();

    // Physikalischen Zustand extrahieren fuer RF
    double n_ion = 0;
    for (int i = 0; i < (int)sys.species.size(); ++i)
        if (sys.species[i].is_positive_ion()) n_ion += state[i];
    double ng = 0;
    for (int i = 0; i < (int)sys.species.size(); ++i)
        if (sys.species[i].is_neutral()) ng += state[i];
    double Te = state[sys.Te_idx()];

    RFState rf = compute_rf(ctx, n_ion, ng, Te, P_RFG);
    if (rf_out) *rf_out = rf;
    if (!rf.valid) return std::vector<double>(N, std::numeric_limits<double>::quiet_NaN());

    double P_abs_V = rf.P_abs * ctx.solver.P_abs_scale / t.V;

    auto raw = assemble_residual(sys, state, P_abs_V, t.Q0, t, ctx,
                                  ctx.solver.alpha_e_wall, ctx.solver.density_profile_factor);

    // Skalierung: physikalisch motiviert (analog zu Legacy residual_scaled)
    // Spezies-Bilanzen: skaliert mit dominantem Produktionsterm
    // Energiebilanzen: skaliert mit Leistungsdichte oder Waermeleitung
    std::vector<double> scaled(N);
    double ne = electron_density(sys, state);
    for (int i = 0; i < N; ++i) {
        if (!std::isfinite(raw[i])) return std::vector<double>(N, std::numeric_limits<double>::quiet_NaN());
        double sc;
        if (i == sys.Te_idx()) {
            // Elektronenenergie: P_abs/V als Skala
            sc = std::max(1e-6, rf.P_abs / t.V);
        } else if (i == sys.Tg_idx()) {
            // Gasenergie: Waermeleitung + elastische Heizung
            double Te = state[sys.Te_idx()], Tg = state[sys.Tg_idx()];
            double ng_total = 0;
            for (int j = 0; j < (int)sys.species.size(); ++j)
                if (sys.species[j].is_neutral()) ng_total += state[j];
            sc = std::max(1e-6, std::fabs(sys.species[0].thermal_cond * (Tg - PhysConst::Tg0) / t.lambda_0 * t.A / t.V)
                        + std::fabs(3.0*PhysConst::me/sys.species[0].mass_kg*PhysConst::kB*(Te*PhysConst::conv-Tg)*ne*ng_total*1e-13));
        } else {
            // Spezies-Bilanz: ne*ng*Kiz als Skala (wie Legacy scale1/scale2)
            double Te = state[sys.Te_idx()];
            double ng_main = 0;
            for (int j = 0; j < (int)sys.species.size(); ++j)
                if (sys.species[j].is_feedstock) ng_main = state[j];
            // Ionisationsrate (erster CTX_KIZ-Reaktionskoeffizient)
            double Kiz_val = 1e-15;
            for (auto& rxn : sys.reactions) {
                if (rxn.rate.type == RateType::CTX_KIZ) {
                    Kiz_val = rxn.rate.evaluate(ctx, Te); break;
                }
            }
            double ion_rate = std::fabs(ne * ng_main * Kiz_val);
            // Fuer Ionen: Ionisationsrate; fuer Neutrale: Q0/V
            if (sys.species[i].is_positive_ion())
                sc = std::max(1e-20, ion_rate);
            else
                sc = std::max(1e-20, std::max(ion_rate, t.Q0 / t.V));
        }
        scaled[i] = raw[i] / sc;
    }
    return scaled;
}

// ═════════════════════════════════════════════════════════════
// Generischer LM-Solver
// ═════════════════════════════════════════════════════════════
inline GenSolveResult gen_solve_lm(
    const ChemSystem& sys, const SimContext& ctx,
    double P_RFG, const std::vector<double>& initial)
{
    const auto& sp = ctx.solver;
    StateLayout lay(sys);
    int N = lay.N;

    GenSolveResult out;
    out.state = initial;
    if (!gen_state_valid(initial, lay, sp)) { out.reason = "invalid initial"; return out; }

    // Log-Koordinaten fuer Dichten, linear fuer Te/Tg
    std::vector<double> x(N);
    for (int i = 0; i < lay.n_species; ++i) x[i] = std::log(initial[i]);
    x[lay.Te_idx] = std::log(initial[lay.Te_idx]);
    x[lay.Tg_idx] = std::log(initial[lay.Tg_idx]);

    auto from_log = [&](const std::vector<double>& u) -> std::vector<double> {
        std::vector<double> s(N);
        for (int i = 0; i < N; ++i) s[i] = std::exp(std::clamp(u[i], -50.0, 50.0));
        s[lay.Te_idx] = std::clamp(s[lay.Te_idx], sp.Te_min, sp.Te_max);
        s[lay.Tg_idx] = std::clamp(s[lay.Tg_idx], sp.Tg_min, sp.Tg_max);
        return s;
    };

    auto state = from_log(x);
    RFState rf;
    auto F = gen_residual_scaled(sys, state, P_RFG, ctx, &rf);
    if (!rf.valid || gen_merit(F) > 1e29) {
        out.reason = "invalid initial residual"; return out;
    }

    double cost = 0; for (double v : F) cost += v*v; cost *= 0.5;
    double lambda = 1e-2;
    const double lam_min = 1e-12, lam_max = 1e8;
    const int max_iter = 120, max_tries = 15, max_stag = 15;
    int stag = 0; double prev_cost = cost;

    for (int iter = 0; iter < max_iter; ++iter) {
        double m = gen_merit(F);
        if (m < sp.newton_tol) {
            out.converged = true; out.state = state; out.rf = rf;
            out.iterations = iter; out.resid_norm = m; out.reason = "ok";
            return out;
        }
        if (!std::isfinite(m)) { out.reason = "nan"; out.state = state; out.iterations = iter; out.resid_norm = m; return out; }

        // Jacobian (zentrale FD)
        std::vector<std::vector<double>> J(N, std::vector<double>(N, 0));
        bool jok = true;
        for (int j = 0; j < N && jok; ++j) {
            auto xp = x, xm = x;
            double h = sp.newton_fd_eps * std::max(1.0, std::fabs(x[j]));
            xp[j] += h; xm[j] -= h;
            auto rp = gen_residual_scaled(sys, from_log(xp), P_RFG, ctx, nullptr);
            auto rm = gen_residual_scaled(sys, from_log(xm), P_RFG, ctx, nullptr);
            for (int i = 0; i < N; ++i) {
                if (!std::isfinite(rp[i]) || !std::isfinite(rm[i])) { jok = false; break; }
                J[i][j] = (rp[i] - rm[i]) / (2*h);
            }
        }
        if (!jok) { out.reason = "jacobian nan"; out.state = state; out.iterations = iter; out.resid_norm = m; return out; }

        // J^T J und J^T F
        std::vector<std::vector<double>> JtJ(N, std::vector<double>(N, 0));
        std::vector<double> JtF(N, 0);
        double dmax = 0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                double s = 0; for (int k = 0; k < N; ++k) s += J[k][i]*J[k][j]; JtJ[i][j] = s;
            }
            double sf = 0; for (int k = 0; k < N; ++k) sf += J[k][i]*F[k]; JtF[i] = sf;
            dmax = std::max(dmax, JtJ[i][i]);
        }

        // LM-Schritt
        bool accepted = false;
        for (int tr = 0; tr < max_tries; ++tr) {
            auto A = JtJ;
            for (int i = 0; i < N; ++i)
                A[i][i] += lambda * std::max(JtJ[i][i], 1e-10*std::max(dmax, 1e-20));
            std::vector<double> neg_JtF(N), dx(N);
            for (int i = 0; i < N; ++i) neg_JtF[i] = -JtF[i];
            if (!solve_NxN(N, A, neg_JtF, dx)) {
                lambda = std::min(lambda*4, lam_max); if (lambda >= lam_max) break; continue;
            }
            for (int i = 0; i < N; ++i)
                dx[i] = std::clamp(dx[i], -sp.newton_max_log_step, sp.newton_max_log_step);

            std::vector<double> xt(N);
            for (int i = 0; i < N; ++i) xt[i] = x[i] + dx[i];
            auto st = from_log(xt);
            if (!gen_state_valid(st, lay, sp)) {
                lambda = std::min(lambda*4, lam_max); if (lambda >= lam_max) break; continue;
            }
            RFState rft;
            auto Ft = gen_residual_scaled(sys, st, P_RFG, ctx, &rft);
            double ct = 0; bool ok = rft.valid;
            for (double v : Ft) { if (!std::isfinite(v)) { ok=false; break; } ct += v*v; } ct *= 0.5;
            if (!ok || !std::isfinite(ct)) { lambda = std::min(lambda*4, lam_max); if (lambda >= lam_max) break; continue; }

            if (cost - ct > 0) {
                x = xt; state = st; rf = rft; F = Ft; cost = ct; accepted = true;
                double rho = (cost - ct) / std::max(1e-30, cost);
                if (rho > 0.75) lambda = std::max(lambda/3, lam_min);
                else if (rho < 0.25) lambda = std::min(lambda*2, lam_max);
                break;
            }
            lambda = std::min(lambda*4, lam_max); if (lambda >= lam_max) break;
        }
        if (!accepted) { out.reason = "step rejected"; out.state = state; out.rf = rf; out.iterations = iter; out.resid_norm = m; return out; }

        double rr = (prev_cost - cost) / std::max(prev_cost, 1e-30);
        if (rr < 1e-10) { stag++; if (stag >= max_stag) { out.reason = "stagnation"; out.state = state; out.rf = rf; out.iterations = iter; out.resid_norm = gen_merit(F); return out; } } else stag = 0;
        prev_cost = cost;
    }
    double fm = gen_merit(F);
    out.converged = fm < sp.newton_tol; out.state = state; out.rf = rf; out.iterations = max_iter; out.resid_norm = fm;
    out.reason = out.converged ? "ok" : "max iter";
    return out;
}

// ═════════════════════════════════════════════════════════════
// Multi-Start Wrapper
// ═════════════════════════════════════════════════════════════
inline GenSolveResult gen_solve_multistart(
    const ChemSystem& sys, const SimContext& ctx, double P_RFG,
    const std::vector<std::vector<double>>& starts)
{
    GenSolveResult best; best.reason = "all starts failed";
    for (auto& s0 : starts) {
        if (!gen_state_valid(s0, StateLayout(sys), ctx.solver)) continue;
        auto r = gen_solve_lm(sys, ctx, P_RFG, s0);
        if (r.converged) return r;
        if (r.resid_norm < best.resid_norm) best = r;
    }
    // Soft-accept: resid < 5*tol und signifikante Verbesserung
    if (best.resid_norm < 5.0 * ctx.solver.newton_tol) {
        best.soft_ok = true;
        best.converged = true;
        best.reason = "soft-converged (" + std::to_string(best.resid_norm) + ")";
    }
    return best;
}

// ═════════════════════════════════════════════════════════════
// Default-Startzustand aus ChemSystem
// ═════════════════════════════════════════════════════════════
inline std::vector<double> gen_default_state(const ChemSystem& sys, double Q0_pps, double V) {
    std::vector<double> state(sys.state_size(), 1e16);
    for (int i = 0; i < (int)sys.species.size(); ++i) {
        if (sys.species[i].is_feedstock)
            state[i] = std::max(1e16, Q0_pps * 4.0 / (V * std::sqrt(8*PhysConst::kB*300/(PhysConst::pi*sys.species[i].mass_kg))));
        else if (sys.species[i].is_positive_ion())
            state[i] = 1e17;
        else
            state[i] = 1e16;
    }
    state[sys.Te_idx()] = 3.75;
    state[sys.Tg_idx()] = 300.0;
    return state;
}

// ═════════════════════════════════════════════════════════════
// Beam-Strom aus generischem Zustand
// ═════════════════════════════════════════════════════════════
inline double gen_beam_current_mA(const ChemSystem& sys, const SimContext& ctx,
                                    const std::vector<double>& state) {
    PlasmaState ps;
    // Summe aller positiven Ionen als n
    ps.n = 0;
    for (int i = 0; i < (int)sys.species.size(); ++i)
        if (sys.species[i].is_positive_ion()) ps.n += state[i];
    // Feedstock als ng
    for (int i = 0; i < (int)sys.species.size(); ++i)
        if (sys.species[i].is_feedstock) { ps.ng = state[i]; break; }
    ps.Te = state[sys.Te_idx()];
    ps.Tg = state[sys.Tg_idx()];
    return beam_current_mA(ctx, ps);
}

#endif // GENERIC_LM_HPP
