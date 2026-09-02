// solver.cpp -- LM-Solver, PTC, Multi-Start, I-fix, SC-Modus.
#include "solver.hpp"
#include "sim_logging.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <limits>

using namespace PhysConst;

// ═══ 4x4 Gauss mit Pivotisierung ══════════════════════════

static bool solve_4x4(double A[4][4], double b[4], double x[4]) {
    double M[4][5];
    for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) M[i][j] = A[i][j]; M[i][4] = b[i]; }
    for (int c = 0; c < 4; ++c) {
        int piv = c;
        for (int r = c+1; r < 4; ++r) if (std::fabs(M[r][c]) > std::fabs(M[piv][c])) piv = r;
        if (std::fabs(M[piv][c]) < 1e-20 || !std::isfinite(M[piv][c])) return false;
        if (piv != c) for (int j = c; j < 5; ++j) std::swap(M[piv][j], M[c][j]);
        double d = M[c][c];
        for (int j = c; j < 5; ++j) M[c][j] /= d;
        for (int r = 0; r < 4; ++r) { if (r==c) continue; double f = M[r][c]; for (int j = c; j < 5; ++j) M[r][j] -= f*M[c][j]; }
    }
    for (int i = 0; i < 4; ++i) x[i] = M[i][4];
    return true;
}

// ═══ Soft-Accept ═══════════════════════════════════════════

static bool soft_accept(const SimContext& ctx, const SolveResult& r, double initial_merit) {
    const auto& sp = ctx.solver;
    if (!state_finite_positive(r.state) || !state_in_bounds(sp, r.state) || !r.rf.valid) return false;
    if (!std::isfinite(r.resid_norm) || !std::isfinite(initial_merit) || initial_merit <= 0) return false;
    double rel = sp.soft_rel_improve * initial_merit;
    return r.resid_norm <= rel && r.resid_norm <= sp.soft_abs_resid && r.iterations >= 2;
}

// ═══ Levenberg-Marquardt ═══════════════════════════════════

SolveResult solve_lm(const SimContext& ctx, double P_RFG_local, const PlasmaState& initial) {
    const auto& sp = ctx.solver;
    SolveResult out; out.state = initial;
    if (!state_in_bounds(sp, initial)) { out.reason = "invalid initial state"; return out; }

    std::array<double,4> x = {std::log(initial.n), std::log(initial.ng), std::log(initial.Te), std::log(initial.Tg)};
    PlasmaState s{std::exp(x[0]), std::exp(x[1]), std::exp(x[2]), std::exp(x[3])};
    RFState rf;
    auto F = residual_scaled(ctx, s, P_RFG_local, &rf);
    if (!rf.valid || !std::isfinite(F[0]) || !std::isfinite(F[1]) || !std::isfinite(F[2]) || !std::isfinite(F[3])) {
        out.reason = "invalid initial residual"; out.resid_norm = 1e30; return out;
    }

    double cost = 0; for (double v : F) cost += v*v; cost *= 0.5;
    double lambda = 1e-2;
    const double lam_min = 1e-12, lam_max = 1e8;
    const int max_iter = 80, max_tries = 15, max_stag = 8;
    int stag = 0; double prev_cost = cost;

    for (int iter = 0; iter < max_iter; ++iter) {
        double m = merit(F, sp, s);
        if (m < sp.newton_tol) {
            out.converged = true; out.state = s; out.rf = rf;
            out.iterations = iter; out.resid_norm = m; out.reason = "ok";
            debug_emit(ctx, 2, "LM_OK", "iter=" + std::to_string(iter) + " merit=" + std::to_string(m));
            return out;
        }
        std::cout << "NEWTON_IT " << iter << " " << std::fixed << std::setprecision(4) << P_RFG_local << " "
                  << std::scientific << std::setprecision(6) << m << " " << s.n << " " << s.ng << " "
                  << std::fixed << std::setprecision(6) << s.Te << " " << s.Tg << std::endl;
        if (!std::isfinite(m)) { out.reason = "invalid residual"; out.state = s; out.rf = rf; out.iterations = iter; out.resid_norm = m; return out; }

        // Jacobian (zentrale FD)
        double J[4][4]; bool jok = true;
        for (int j = 0; j < 4 && jok; ++j) {
            std::array<double,4> xp = x, xm = x;
            double h = sp.newton_fd_eps * std::max(1.0, std::fabs(x[j]));
            xp[j] += h; xm[j] -= h;
            PlasmaState sp_{std::exp(xp[0]),std::exp(xp[1]),std::exp(xp[2]),std::exp(xp[3])};
            PlasmaState sm_{std::exp(xm[0]),std::exp(xm[1]),std::exp(xm[2]),std::exp(xm[3])};
            auto rp = residual_scaled(ctx, sp_, P_RFG_local, nullptr);
            auto rm = residual_scaled(ctx, sm_, P_RFG_local, nullptr);
            for (int i = 0; i < 4; ++i) { if (!std::isfinite(rp[i])||!std::isfinite(rm[i])) { jok=false; break; } J[i][j] = (rp[i]-rm[i])/(2*h); }
        }
        if (!jok) { out.reason = "jacobian nan"; out.state = s; out.rf = rf; out.iterations = iter; out.resid_norm = m; return out; }

        double JtJ[4][4]={}, JtF[4]={}, dmax = 0;
        for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) { double sm = 0; for (int k = 0; k < 4; ++k) sm += J[k][i]*J[k][j]; JtJ[i][j] = sm; }
            double sf = 0; for (int k = 0; k < 4; ++k) sf += J[k][i]*F[k]; JtF[i] = sf; dmax = std::max(dmax, JtJ[i][i]); }

        double gn = 0; for (int i = 0; i < 4; ++i) gn += JtF[i]*JtF[i];
        if (gn < 1e-28*std::max(1.0, cost*cost)) { out.reason = "local minimum"; out.state = s; out.rf = rf; out.iterations = iter; out.resid_norm = m; return out; }

        bool accepted = false;
        for (int tr = 0; tr < max_tries; ++tr) {
            double A[4][4];
            for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) A[i][j] = JtJ[i][j]; A[i][i] += lambda*std::max(JtJ[i][i], 1e-10*std::max(dmax, 1e-20)); }
            double nb[4] = {-JtF[0],-JtF[1],-JtF[2],-JtF[3]}, dx[4];
            if (!solve_4x4(A, nb, dx)) { lambda = std::min(lambda*4, lam_max); if (lambda >= lam_max) break; continue; }
            for (int i = 0; i < 4; ++i) dx[i] = std::max(-sp.newton_max_log_step, std::min(sp.newton_max_log_step, dx[i]));
            std::array<double,4> xt; for (int i = 0; i < 4; ++i) xt[i] = x[i]+dx[i];
            PlasmaState st{std::exp(xt[0]),std::exp(xt[1]),std::exp(xt[2]),std::exp(xt[3])};
            if (!state_in_bounds(sp, st)) { lambda = std::min(lambda*4, lam_max); if (lambda >= lam_max) break; continue; }
            RFState rft; auto Ft = residual_scaled(ctx, st, P_RFG_local, &rft);
            double ct = 0; bool ok = rft.valid;
            for (double v : Ft) { if (!std::isfinite(v)) { ok=false; break; } ct += v*v; } ct *= 0.5;
            if (!ok || !std::isfinite(ct)) { lambda = std::min(lambda*4, lam_max); if (lambda >= lam_max) break; continue; }
            double act = cost - ct;
            double pred = 0; { double jj=0,dp=0; for (int i=0;i<4;++i){double jd=0;for(int j=0;j<4;++j)jd+=JtJ[i][j]*dx[j];jj+=dx[i]*jd;dp+=lambda*std::max(JtJ[i][i],1e-10*std::max(dmax,1e-20))*dx[i]*dx[i];}pred=0.5*jj+dp;if(pred<=0)pred=1e-30;}
            double rho = act/pred;
            if (act > 0) { x = xt; s = st; rf = rft; F = Ft; cost = ct; accepted = true;
                if (rho > 0.75) lambda = std::max(lambda/3, lam_min); else if (rho < 0.25) lambda = std::min(lambda*2, lam_max); break;
            } else { lambda = std::min(lambda*4, lam_max); if (lambda >= lam_max) break; }
        }
        if (!accepted) { out.reason = "lm step rejected"; out.state = s; out.rf = rf; out.iterations = iter; out.resid_norm = m; return out; }
        double rr = (prev_cost-cost)/std::max(prev_cost, 1e-30);
        if (rr < 1e-10) { stag++; if (stag >= max_stag) { out.reason = "stagnation"; out.state = s; out.rf = rf; out.iterations = iter; out.resid_norm = merit(F,sp,s); return out; } } else stag = 0;
        prev_cost = cost;
    }
    double fm = merit(F, sp, s);
    out.converged = fm < sp.newton_tol; out.state = s; out.rf = rf; out.iterations = max_iter; out.resid_norm = fm;
    out.reason = out.converged ? "ok" : "lm max iter";
    return out;
}

// ═══ PTC + LM ═════════════════════════════════════════════

SolveResult solve_ptc_then_lm(const SimContext& ctx, double P_RFG_local, const PlasmaState& initial) {
    const auto& sp = ctx.solver;
    SolveResult best; best.state = initial;
    if (!state_in_bounds(sp, initial)) { best.reason = "invalid initial state"; return best; }

    RFState rf0; auto r0 = residual_scaled(ctx, initial, P_RFG_local, &rf0);
    double mi = merit(r0, sp, initial);

    best = solve_lm(ctx, P_RFG_local, initial);
    if (best.converged) return best;
    if (soft_accept(ctx, best, mi)) { best.soft_ok = true; best.reason = "soft-ok-direct"; return best; }

    std::array<double,4> x = {std::log(initial.n), std::log(initial.ng), std::log(initial.Te), std::log(initial.Tg)};
    PlasmaState s = initial; RFState rf = rf0; auto r = r0; double m = mi;
    if (!std::isfinite(m) || !rf.valid) return best;

    double gain = sp.ptc_start_gain;
    SolveResult bl = best;
    if (std::isfinite(m) && m < bl.resid_norm) { bl.state = s; bl.rf = rf; bl.resid_norm = m; bl.iterations = 0; }

    for (int it = 0; it < sp.ptc_max_iter; ++it) {
        bool acc = false; std::array<double,4> xt = x;
        for (double al : {1.0, 0.5, 0.25, 0.125, 0.0625}) {
            for (int j = 0; j < 4; ++j) { double st = -al*gain*r[j]; xt[j] = x[j] + std::max(-0.12, std::min(0.12, st)); }
            PlasmaState st{std::exp(xt[0]),std::exp(xt[1]),std::exp(xt[2]),std::exp(xt[3])};
            if (!state_in_bounds(sp, st)) continue;
            RFState rft; auto rt = residual_scaled(ctx, st, P_RFG_local, &rft);
            double mt = merit(rt, sp, st);
            if (!std::isfinite(mt) || !rft.valid) continue;
            if (mt < sp.ptc_accept_ratio * m) {
                x = xt; s = st; r = rt; rf = rft; m = mt; acc = true;
                if (m < bl.resid_norm) { bl.state = s; bl.rf = rf; bl.resid_norm = m; bl.iterations = it+1; }
                break;
            }
        }
        if (!acc) { gain *= 0.5; if (gain < sp.ptc_min_gain) break; continue; }
        if (m < sp.ptc_switch_merit) {
            SolveResult pol = solve_lm(ctx, P_RFG_local, s);
            if (pol.converged) return pol;
            if (soft_accept(ctx, pol, mi)) { pol.soft_ok = true; pol.reason = "soft-ok-polished"; return pol; }
            if (pol.resid_norm < bl.resid_norm) bl = pol;
        }
    }
    SolveResult fin = solve_lm(ctx, P_RFG_local, bl.state);
    if (fin.converged) return fin;
    if (soft_accept(ctx, fin, mi)) { fin.soft_ok = true; fin.reason = "soft-ok-final"; return fin; }
    if (soft_accept(ctx, bl, mi)) { bl.soft_ok = true; return bl; }
    return (fin.resid_norm < bl.resid_norm) ? fin : bl;
}

// ═══ Hilfsfunktionen fuer Multi-Start ══════════════════════

static std::vector<PlasmaState> make_starts(const SimContext& ctx, const PlasmaState& guess) {
    PlasmaState safe = safe_defaults(ctx);
    std::vector<PlasmaState> starts = {guess, safe};
    if (state_in_bounds(ctx.solver, guess))
        starts.push_back(PlasmaState{
            std::sqrt(std::max(1.0, guess.n*safe.n)),
            std::sqrt(std::max(1.0, guess.ng*safe.ng)),
            0.5*(guess.Te+safe.Te), 0.5*(guess.Tg+safe.Tg)});
    return starts;
}

static SolveResult multi_start_solve(const SimContext& ctx, double P_RFG, const PlasmaState& guess, bool use_ptc) {
    auto starts = make_starts(ctx, guess);
    SolveResult best; best.reason = "all starts failed";
    for (const auto& st : starts) {
        if (!state_in_bounds(ctx.solver, st)) continue;
        SolveResult cur = use_ptc ? solve_ptc_then_lm(ctx, P_RFG, st) : solve_lm(ctx, P_RFG, st);
        if (cur.converged) return cur;
        if (cur.resid_norm < best.resid_norm) best = cur;
    }
    return best;
}

// ═══ SC-Modus ═════════════════════════════════════════════

PowerSolveResult solve_at_fixed_power(const SimContext& ctx, double P_RFG_fixed, const PlasmaState& guess) {
    PowerSolveResult out;
    SolveResult best = multi_start_solve(ctx, P_RFG_fixed, guess, true);
    out.iterations = best.iterations; out.inner_resid_norm = best.resid_norm;
    out.P_RFG_sol = P_RFG_fixed; out.P_trial_last = P_RFG_fixed;
    if (best.converged && state_finite_positive(best.state)) {
        out.converged = true; out.state = best.state; out.rf = best.rf;
        out.I_mA = beam_current_mA(ctx, best.state); out.err_mA = 0; out.reason = "ok";
        std::cout << "PID_DONE " << std::fixed << std::setprecision(4)
                  << out.I_mA << " 0.0000 " << P_RFG_fixed << " " << best.state.Te << " " << best.state.Tg << std::endl;
        std::cout << "CONVERGED " << best.iterations << std::endl;
    } else {
        out.fail_type = SolveFailType::NUMERICAL_FAIL;
        out.state = best.state; out.rf = best.rf;
        out.I_mA = state_finite_positive(best.state) ? beam_current_mA(ctx, best.state) : 0;
        out.reason = "self-consistent failed: " + best.reason;
    }
    return out;
}

// ═══ I-fix-Modus ══════════════════════════════════════════

PowerSolveResult solve_for_target_current(const SimContext& ctx, const PlasmaState& guess) {
    const auto& sp = ctx.solver;
    PowerSolveResult out;

    auto solve_at = [&](double p, const PlasmaState& start) { return multi_start_solve(ctx, p, start, false); };

    PlasmaState warm = state_in_bounds(sp, guess) ? guess : safe_defaults(ctx);
    double p_lo = 5.0;
    SolveResult s_lo = solve_at(p_lo, warm);
    double I_lo = (s_lo.converged && state_finite_positive(s_lo.state)) ? beam_current_mA(ctx, s_lo.state) : 0;
    double f_lo = I_lo - sp.I_soll;
    if (s_lo.converged) warm = s_lo.state;

    double p_hi = 10, I_hi = 0, f_hi = -sp.I_soll;
    double p_last_good = p_lo; PlasmaState st_good = warm; double I_good = I_lo;
    SolveResult s_hi;

    const double P_max = std::max(sp.P_RFG_max, 200.0);
    bool bracket = false;
    for (int ex = 0; p_hi <= P_max && ex < 20; ++ex) {
        s_hi = solve_at(p_hi, st_good);
        bool ok = s_hi.converged && state_finite_positive(s_hi.state);
        I_hi = ok ? beam_current_mA(ctx, s_hi.state) : 0; f_hi = I_hi - sp.I_soll;
        std::cout << "POWER_BRACKET " << std::fixed << std::setprecision(4) << p_lo << " " << p_hi << " " << f_lo << " " << f_hi << " hi_conv=" << ok << std::endl;
        if (ok) {
            p_last_good = p_hi; st_good = s_hi.state; I_good = I_hi;
            if (f_lo*f_hi < 0) { bracket = true; break; }
            p_hi = std::min(p_hi*2, P_max+1);
        } else {
            double pp = 0.5*(p_last_good+p_hi);
            for (int r = 0; r < 8; ++r) {
                auto sp_ = solve_at(pp, st_good);
                bool pok = sp_.converged && state_finite_positive(sp_.state);
                if (pok) { p_last_good = pp; st_good = sp_.state; I_good = beam_current_mA(ctx, sp_.state); pp = 0.5*(pp+p_hi); }
                else { p_hi = pp; pp = 0.5*(p_last_good+pp); }
            }
            p_hi = p_last_good; s_hi = solve_at(p_hi, st_good);
            I_hi = (s_hi.converged && state_finite_positive(s_hi.state)) ? beam_current_mA(ctx, s_hi.state) : I_good;
            f_hi = I_hi - sp.I_soll;
            if (f_lo*f_hi < 0) bracket = true;
            break;
        }
    }

    if (!bracket) {
        out.hit_limit = true; out.P_trial_last = p_hi;
        auto& bs = (std::fabs(f_lo) < std::fabs(f_hi)) ? s_lo : s_hi;
        if (state_finite_positive(bs.state)) { out.state = bs.state; out.rf = bs.rf; out.I_mA = beam_current_mA(ctx, bs.state); out.err_mA = sp.I_soll - out.I_mA; }
        out.fail_type = (s_lo.converged && s_hi.converged) ? SolveFailType::NO_PHYSICAL_SOLUTION : SolveFailType::NUMERICAL_FAIL;
        out.reason = "no bracket (I_lo=" + std::to_string(I_lo) + " I_hi=" + std::to_string(I_hi) + ")";
        return out;
    }

    // Bisection + Regula Falsi
    PlasmaState lg = s_lo.converged ? s_lo.state : warm;
    double be = 1e30; SolveResult bs; double bp = 0, bi = 0;
    int ifc = 0;
    for (int it = 0; it < 60; ++it) {
        double pm; double ds = f_hi - f_lo;
        if (std::fabs(ds) > 1e-12) { pm = p_lo - f_lo*(p_hi-p_lo)/ds; if (!(pm>p_lo&&pm<p_hi)||!std::isfinite(pm)) pm = 0.5*(p_lo+p_hi); } else pm = 0.5*(p_lo+p_hi);
        std::cout << "PID_START " << it << " " << std::fixed << std::setprecision(4) << pm << std::endl;
        auto sm = solve_at(pm, lg);
        if (!sm.converged) { pm = 0.5*(p_lo+p_hi); sm = solve_at(pm, safe_defaults(ctx));
            if (!sm.converged) { ifc++; if (ifc >= 5) { out.fail_type = SolveFailType::NUMERICAL_FAIL; out.reason = "inner fail x" + std::to_string(ifc); if (be<1e30){out.state=bs.state;out.rf=bs.rf;out.P_RFG_sol=bp;out.I_mA=bi;out.err_mA=sp.I_soll-bi;} return out; }
                p_hi = 0.5*(p_lo+p_hi); f_hi = f_lo; continue; } }
        double Im = beam_current_mA(ctx, sm.state), err = sp.I_soll - Im, fm = Im - sp.I_soll, ae = std::fabs(err);
        std::cout << "PID_DONE " << std::fixed << std::setprecision(4) << Im << " " << err << " " << pm << " " << sm.state.Te << " " << sm.state.Tg << std::endl;
        lg = sm.state;
        if (ae < be) { be = ae; bs = sm; bp = pm; bi = Im; }
        out.iterations = it; out.state = sm.state; out.rf = sm.rf; out.inner_resid_norm = sm.resid_norm;
        out.P_RFG_sol = pm; out.P_trial_last = pm; out.I_mA = Im; out.err_mA = err;
        if (ae < sp.power_tol_mA) { out.converged = true; out.reason = "ok"; std::cout << "CONVERGED " << it << std::endl; return out; }
        double pw = p_hi - p_lo;
        if (ae < 10*sp.power_tol_mA && pw < 0.1) { out.converged = true; out.state = bs.state; out.rf = bs.rf; out.P_RFG_sol = bp; out.I_mA = bi; out.err_mA = sp.I_soll-bi; out.reason = "near-hit"; std::cout << "CONVERGED " << it << std::endl; return out; }
        if (pw < 0.01) { if (be < 1) { out.converged = true; out.state = bs.state; out.rf = bs.rf; out.P_RFG_sol = bp; out.I_mA = bi; out.err_mA = sp.I_soll-bi; out.reason = "plateau-accept"; std::cout << "CONVERGED " << it << std::endl; } else { out.fail_type = SolveFailType::NO_PHYSICAL_SOLUTION; out.reason = "plateau"; } return out; }
        if (f_lo*fm <= 0) { p_hi = pm; f_hi = fm; } else { p_lo = pm; f_lo = fm; }
    }
    if (be < 10*sp.power_tol_mA) { out.converged = true; out.state = bs.state; out.rf = bs.rf; out.P_RFG_sol = bp; out.I_mA = bi; out.err_mA = sp.I_soll-bi; out.reason = "max-iter-accept"; std::cout << "CONVERGED 60" << std::endl; return out; }
    out.fail_type = SolveFailType::NUMERICAL_FAIL; out.reason = "bisection max iter"; return out;
}

bool power_result_valid(const SimContext& ctx, const PowerSolveResult& ps) {
    return ps.converged && state_finite_positive(ps.state) && state_in_bounds(ctx.solver, ps.state)
        && ps.rf.valid && std::isfinite(ps.P_RFG_sol) && ps.P_RFG_sol > 0 && std::isfinite(ps.I_mA);
}
