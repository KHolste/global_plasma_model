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

    // Stossfrequenz aus dem Reaktionsnetz; ist dort nichts gekennzeichnet,
    // bleibt es beim bisherigen Weg ueber den einen elastischen Koeffizienten.
    double nu_m = chem_nu_m(sys, state, ctx, Te);
    RFState rf = (nu_m >= 0) ? compute_rf_nu(ctx, n_ion, nu_m, P_RFG)
                             : compute_rf(ctx, n_ion, ng, Te, P_RFG);
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

            // Guete des Schritts: tatsaechliche gegen vorhergesagte Senkung.
            // Beides vor dem Uebernehmen bilden -- sonst ist das Verhaeltnis
            // immer null und die Daempfung waechst bei jedem Schritt.
            double act = cost - ct;
            if (act > 0) {
                double quad = 0, damp = 0;
                for (int i = 0; i < N; ++i) {
                    double jd = 0;
                    for (int j = 0; j < N; ++j) jd += JtJ[i][j] * dx[j];
                    quad += dx[i] * jd;
                    damp += lambda * std::max(JtJ[i][i], 1e-10*std::max(dmax, 1e-20)) * dx[i] * dx[i];
                }
                double pred = 0.5*quad + damp;
                if (pred <= 0) pred = 1e-30;
                double rho = act / pred;
                x = xt; state = st; rf = rft; F = Ft; cost = ct; accepted = true;
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
// Extraktion aus dem generischen Zustand
//
// Jede positive Ionensorte geht mit ihrer eigenen Masse und Ladungszahl ein,
// jede neutrale Sorte mit ihrer eigenen Masse in den Neutralschub. Der
// zugefuehrte Massenstrom ist der Teilchenstrom mal der Masse der zugefuehrten
// Sorte -- bei molekularem Treibstoff also die Molekuelmasse.
// ═════════════════════════════════════════════════════════════
inline ExtractionInput gen_extraction_input(const ChemSystem& sys, const SimContext& ctx,
                                             const std::vector<double>& state) {
    ExtractionInput in;
    for (int i = 0; i < (int)sys.species.size(); ++i) {
        const ChemSpecies& sp = sys.species[i];
        if (sp.is_positive_ion())
            in.ions.push_back({sp.id, state[i], sp.mass_kg, sp.charge});
        else if (sp.is_neutral())
            in.neutrals.push_back({sp.id, state[i], sp.mass_kg});
    }
    in.Te_eV = state[sys.Te_idx()];
    in.Tg_K = state[sys.Tg_idx()];
    in.sigma_i = sys.sigma_i;
    const ChemSpecies* feed = sys.feedstock();
    in.mdot_in_kg_s = feed ? ctx.thruster.Q0 * feed->mass_kg : 0.0;
    return in;
}

inline ExtractionResult gen_extraction(const ChemSystem& sys, const SimContext& ctx,
                                        const std::vector<double>& state) {
    return compute_extraction(gen_extraction_input(sys, ctx, state), ctx.thruster);
}

inline double gen_beam_current_mA(const ChemSystem& sys, const SimContext& ctx,
                                    const std::vector<double>& state) {
    return gen_extraction(sys, ctx, state).I_beam_mA;
}

// ═════════════════════════════════════════════════════════════
// Betriebsarten: feste Leistung und fester Strahlstrom
//
// Gegenstueck zu solve_at_fixed_power/solve_for_target_current des fest
// verdrahteten Pfades, aber ueber dem generischen Zustandsvektor. Die
// Ausgabezeilen sind dieselben, damit Oberflaeche und Protokoll den
// Unterschied nicht bemerken muessen.
// ═════════════════════════════════════════════════════════════

struct GenPowerResult {
    bool converged = false;
    SolveFailType fail_type = SolveFailType::NONE;
    std::vector<double> state;
    RFState rf;
    double P_RFG_sol = 0, P_trial_last = 0;
    double I_mA = 0, err_mA = 0;
    double inner_resid_norm = 1e30;
    int iterations = 0;
    std::string reason;
};

// Startzustaende: uebergebener Warmstart, Vorgabezustand, und deren
// geometrisches Mittel. Dieselbe Staffelung wie im festen Pfad.
inline std::vector<std::vector<double>> gen_starts(
    const ChemSystem& sys, const SimContext& ctx, const std::vector<double>& guess)
{
    const int N = sys.state_size();
    auto def = gen_default_state(sys, ctx.thruster.Q0, ctx.thruster.V);
    std::vector<std::vector<double>> starts;
    if ((int)guess.size() == N) starts.push_back(guess);
    starts.push_back(def);
    if ((int)guess.size() == N) {
        std::vector<double> mid(N);
        for (int i = 0; i < (int)sys.species.size(); ++i)
            mid[i] = std::sqrt(std::max(1.0, guess[i] * def[i]));
        mid[sys.Te_idx()] = 0.5 * (guess[sys.Te_idx()] + def[sys.Te_idx()]);
        mid[sys.Tg_idx()] = 0.5 * (guess[sys.Tg_idx()] + def[sys.Tg_idx()]);
        starts.push_back(mid);
    }
    return starts;
}

inline GenPowerResult gen_solve_at_fixed_power(
    const ChemSystem& sys, const SimContext& ctx,
    double P_RFG, const std::vector<double>& guess)
{
    GenPowerResult out;
    GenSolveResult r = gen_solve_multistart(sys, ctx, P_RFG, gen_starts(sys, ctx, guess));
    out.state = r.state;
    out.rf = r.rf;
    out.iterations = r.iterations;
    out.inner_resid_norm = r.resid_norm;
    out.P_RFG_sol = P_RFG;
    out.P_trial_last = P_RFG;
    if (r.converged && (int)r.state.size() == sys.state_size()) {
        out.converged = true;
        out.I_mA = gen_beam_current_mA(sys, ctx, r.state);
        out.reason = r.reason;
    } else {
        out.fail_type = SolveFailType::NUMERICAL_FAIL;
        out.reason = "generisch gescheitert: " + r.reason;
    }
    return out;
}

inline GenPowerResult gen_solve_for_target_current(
    const ChemSystem& sys, const SimContext& ctx, const std::vector<double>& guess)
{
    const auto& sp = ctx.solver;
    const double P_max = std::max(sp.P_RFG_max, 200.0);
    GenPowerResult out;

    // Untere Schranke
    double p_lo = std::max(1.0, sp.power_min);
    GenPowerResult lo = gen_solve_at_fixed_power(sys, ctx, p_lo, guess);
    double f_lo = (lo.converged ? lo.I_mA : 0.0) - sp.I_soll;
    std::vector<double> warm = lo.converged ? lo.state : guess;

    // Obere Schranke suchen: verdoppeln, bis der Zielstrom ueberschritten ist
    double p_hi = 2.0 * p_lo, f_hi = -sp.I_soll, p_last_good = p_lo;
    GenPowerResult hi;
    bool bracket = false;
    for (int k = 0; k < 20 && p_hi <= P_max; ++k) {
        hi = gen_solve_at_fixed_power(sys, ctx, p_hi, warm);
        std::cout << "POWER_BRACKET " << std::fixed << std::setprecision(4)
                  << p_lo << " " << p_hi << " " << f_lo << " "
                  << (hi.converged ? hi.I_mA - sp.I_soll : 0.0)
                  << " hi_conv=" << hi.converged << std::endl;
        if (hi.converged) {
            f_hi = hi.I_mA - sp.I_soll;
            warm = hi.state;
            p_last_good = p_hi;
            if (f_lo * f_hi < 0) { bracket = true; break; }
            p_hi = std::min(p_hi * 2.0, P_max + 1.0);
        } else {
            // Oberhalb liegt keine Loesung mehr: zwischen dem letzten guten
            // Punkt und hier weiter einschachteln, sonst aufgeben.
            double mitte = 0.5 * (p_last_good + p_hi);
            if (mitte - p_last_good < 0.25) break;
            p_hi = mitte;
        }
    }

    if (!bracket) {
        const GenPowerResult& best = (std::fabs(f_lo) <= std::fabs(f_hi)) ? lo : hi;
        out = best;
        out.converged = false;
        out.P_trial_last = p_hi;
        out.I_mA = best.converged ? best.I_mA : 0.0;
        out.err_mA = sp.I_soll - out.I_mA;
        out.fail_type = (lo.converged && hi.converged)
                        ? SolveFailType::NO_PHYSICAL_SOLUTION
                        : SolveFailType::NUMERICAL_FAIL;
        out.reason = "keine Einschachtelung (I_lo=" + std::to_string(lo.I_mA) +
                     " I_hi=" + std::to_string(hi.I_mA) + ")";
        return out;
    }

    // Regula falsi mit Rueckfall auf Bisektion
    GenPowerResult best; double best_err = 1e30;
    int inner_fail = 0;
    for (int it = 0; it < 60; ++it) {
        double pm, ds = f_hi - f_lo;
        if (std::fabs(ds) > 1e-12) {
            pm = p_lo - f_lo * (p_hi - p_lo) / ds;
            if (!std::isfinite(pm) || pm <= p_lo || pm >= p_hi) pm = 0.5 * (p_lo + p_hi);
        } else {
            pm = 0.5 * (p_lo + p_hi);
        }
        std::cout << "PID_START " << it << " " << std::fixed << std::setprecision(4) << pm << std::endl;

        GenPowerResult m = gen_solve_at_fixed_power(sys, ctx, pm, warm);
        if (!m.converged) {
            if (++inner_fail >= 5) {
                out = (best_err < 1e30) ? best : m;
                out.converged = false;
                out.fail_type = SolveFailType::NUMERICAL_FAIL;
                out.reason = "innerer Loeser " + std::to_string(inner_fail) + "x gescheitert";
                return out;
            }
            p_hi = pm; f_hi = f_lo;  // Fenster verkleinern und erneut versuchen
            continue;
        }
        warm = m.state;
        double fm = m.I_mA - sp.I_soll, ae = std::fabs(fm);
        std::cout << "PID_DONE " << std::fixed << std::setprecision(4) << m.I_mA << " "
                  << -fm << " " << pm << " " << m.state[sys.Te_idx()] << " "
                  << m.state[sys.Tg_idx()] << std::endl;

        m.P_RFG_sol = pm; m.P_trial_last = pm; m.err_mA = -fm; m.iterations = it;
        if (ae < best_err) { best_err = ae; best = m; }

        if (ae < sp.power_tol_mA) {
            best.converged = true; best.reason = "ok";
            std::cout << "CONVERGED " << it << std::endl;
            return best;
        }
        double fenster = p_hi - p_lo;
        if (fenster < 0.01) {
            if (best_err < 1.0) {
                best.converged = true; best.reason = "Plateau angenommen";
                std::cout << "CONVERGED " << it << std::endl;
            } else {
                best.converged = false;
                best.fail_type = SolveFailType::NO_PHYSICAL_SOLUTION;
                best.reason = "Plateau";
            }
            return best;
        }
        if (f_lo * fm <= 0) { p_hi = pm; f_hi = fm; } else { p_lo = pm; f_lo = fm; }
    }

    if (best_err < 10 * sp.power_tol_mA) {
        best.converged = true; best.reason = "nach Hoechstzahl angenommen";
        std::cout << "CONVERGED 60" << std::endl;
        return best;
    }
    best.converged = false;
    best.fail_type = SolveFailType::NUMERICAL_FAIL;
    best.reason = "Bisektion ohne Treffer";
    return best;
}

#endif // GENERIC_LM_HPP
