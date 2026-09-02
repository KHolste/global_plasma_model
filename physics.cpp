// physics.cpp -- RF-Kopplung, Residuen, abgeleitete Groessen.
#include "physics.hpp"
#include "bessel_wrapper.hpp"
#include <limits>
#include <fstream>
#include <iomanip>

using namespace std::literals::complex_literals;
using namespace PhysConst;

// ═══ RF-Kopplung (Bessel-BVP) ══════════════════════════════

RFState compute_rf(const SimContext& ctx, double n, double ng, double Te, double P_RFG_local) {
    const auto& t = ctx.thruster;
    RFState out{};
    if (!(std::isfinite(n) && std::isfinite(ng) && std::isfinite(Te) && std::isfinite(P_RFG_local))
        || n<=0 || ng<=0 || Te<=0 || P_RFG_local<=0)
        return out;

    double a = plasma_freq(n), b = t.omega, cf = coll_freq(ctx, ng, Te);
    double den = b*b + cf*cf;
    if (!std::isfinite(a) || !std::isfinite(cf) || den <= 0) return out;

    double it = a*a*cf/(b*den);
    double rt = 1.0 - a*a/den;
    double rabs = std::sqrt(it*it + rt*rt);
    double o_rt = std::sqrt(std::max(0.0, (rabs+rt)/2.0));
    double o_it = signum(it) * std::sqrt(std::max(0.0, (rabs-rt)/2.0));

    std::complex<double> k1(o_rt*t.k_0, -o_it*t.k_0);
    std::complex<double> kR = t.R * k1;
    std::complex<double> eps_p(rt, -it);
    std::complex<double> j0 = bessel_cyl_j(0, kR);
    std::complex<double> j1 = bessel_cyl_j(1, kR);
    std::complex<double> dc = eps_p * j0;

    if (!std::isfinite(dc.real()) || !std::isfinite(dc.imag())) return out;
    double dca = std::abs(dc);
    if (dca < 1e-30) dc = dc / std::max(dca, 1e-300) * 1e-30;

    std::complex<double> res = (1i * kR * j1) / dc;
    double R_ind = std::max(1e-4, 2*pi*t.Nw*t.Nw / (epsilon0*t.L*t.omega) * res.real());
    if (!std::isfinite(R_ind)) return out;

    double Ic = std::sqrt(2.0 * P_RFG_local / (R_ind + t.R_ohm));
    double Pabs = 0.5 * R_ind * Ic * Ic;
    if (!std::isfinite(Ic) || !std::isfinite(Pabs)) return out;

    out.P_abs = Pabs; out.R_ind = R_ind; out.I_coil = Ic; out.valid = true;
    return out;
}

// ═══ Residuen ══════════════════════════════════════════════

std::array<double,4> residual_raw(const SimContext& ctx, const PlasmaState& s,
                                   double P_RFG_local, RFState* rf_out) {
    const auto& t = ctx.thruster;
    const auto& g = ctx.gas;
    const auto& sp = ctx.solver;

    RFState rf = compute_rf(ctx, s.n, s.ng, s.Te, P_RFG_local);
    if (rf_out) *rf_out = rf;
    if (!rf.valid) { double nan = std::numeric_limits<double>::quiet_NaN(); return {nan,nan,nan,nan}; }

    double P_vol = rf.P_abs * sp.P_abs_scale / t.V;
    double n_eff = sp.density_profile_factor * s.n;
    double lam = lambda_i(s.ng, g.sigma_i);
    double ub = uB(g.M, s.Te);
    double aeff = Aeff(t, lam);

    // r1: Ionenbilanz
    double r1 = n_eff*s.ng*Kiz(ctx, s.Te) - s.n*ub*aeff/t.V;
    // r2: Neutralbilanz
    double r2 = t.Q0/t.V + s.n*ub*Aeff1(t, lam)/t.V - s.n*s.ng*Kiz(ctx, s.Te) - 0.25*s.ng*vg(g.M, s.Tg)*t.Ag/t.V;
    // r4: Gasenergie
    double Pg1 = 3.0*me/g.M * kB*(s.Te*conv-s.Tg) * s.n*s.ng*Kel(ctx, s.Te);
    double Pg2 = 0.25*g.M*ub*ub * s.n*s.ng*g.sigma_i*vi(g.M, s.Tg);
    double Pg3 = g.kappa*(s.Tg-Tg0)/t.lambda_0 * t.A/t.V;
    double r4 = Pg1 + Pg2 - Pg3;
    // r3: Elektronenenergie
    double P2 = g.Eiz * n_eff*s.ng*Kiz(ctx, s.Te);
    double P3;
    if (ctx.rates.excitation_model == 1 && !ctx.rates.kex.empty())
        P3 = Pexc_coeff(ctx, s.Te) * n_eff * s.ng;
    else
        P3 = g.Eexc * n_eff*s.ng*Kex(ctx, s.Te);
    double P4 = 3.0*me/g.M * kB*(s.Te*conv-s.Tg) * n_eff*s.ng*Kel(ctx, s.Te);
    double P5 = sp.alpha_e_wall*kB*s.Te*conv * s.n*ub*aeff/t.V;
    double r3 = P_vol - (P2 + P3 + P4 + P5);

    return {r1, r2, r3, r4};
}

std::array<double,4> residual_scaled(const SimContext& ctx, const PlasmaState& s,
                                      double P_RFG_local, RFState* rf_out) {
    if (!state_in_bounds(ctx.solver, s)) {
        double nan = std::numeric_limits<double>::quiet_NaN();
        return {nan,nan,nan,nan};
    }
    auto raw = residual_raw(ctx, s, P_RFG_local, rf_out);
    for (double v : raw) if (!std::isfinite(v)) { double nan = std::numeric_limits<double>::quiet_NaN(); return {nan,nan,nan,nan}; }

    RFState rf = rf_out ? *rf_out : compute_rf(ctx, s.n, s.ng, s.Te, P_RFG_local);
    if (!rf.valid || !std::isfinite(rf.P_abs) || rf.P_abs < 0) {
        double nan = std::numeric_limits<double>::quiet_NaN();
        return {nan,nan,nan,nan};
    }
    const auto& t = ctx.thruster;
    const auto& g = ctx.gas;
    double sc1 = std::max(1e-20, std::fabs(s.n * s.ng * Kiz(ctx, s.Te)));
    double sc2 = std::max(1e-20, t.Q0 / t.V);
    double sc3 = std::max(1e-6,  rf.P_abs / t.V);
    double sc4 = std::max(1e-6,  std::fabs(g.kappa*(s.Tg-Tg0)/t.lambda_0*t.A/t.V)
                        + std::fabs(3.0*me/g.M*kB*(s.Te*conv-s.Tg)*s.n*s.ng*Kel(ctx, s.Te)));
    return {raw[0]/sc1, raw[1]/sc2, raw[2]/sc3, raw[3]/sc4};
}

// ═══ Abgeleitete Groessen ══════════════════════════════════

DerivedQuantities compute_derived(const SimContext& ctx, double n, double ng, double Te, double Tg,
                                   double R_ind, double I_coil_val, double P_abs_val) {
    const auto& t = ctx.thruster;
    const auto& g = ctx.gas;
    DerivedQuantities d{};

    double lam = lambda_i(ng, g.sigma_i);
    double Gi = Gamma_i_func(t, lam, Te, n, g.M);

    // Volles Extraktionsmodell (Bohm + CL + eta_opt)
    PlasmaState ps{n, ng, Te, Tg};
    ExtractionResult ex = compute_beam_extraction(ctx, ps);
    d.I_extr_mA = ex.I_beam_mA;
    d.I_CL_limit_mA = ex.I_CL_limit_mA;
    d.I_Bohm_limit_mA = ex.I_Bohm_limit_mA;
    d.eta_opt_used = ex.eta_opt;
    d.beam_limiting = ex.limiting;
    d.iondeg = n / std::max(ng, 1.0) * 100;
    d.cf = coll_freq(ctx, ng, Te);
    d.u_Bohm = uB(g.M, Te);
    d.J_i = e * Gi;
    d.pf = plasma_freq(n);

    double a = d.pf, b = t.omega, cf = d.cf;
    double den = b*b + cf*cf;
    d.eps_p_real = 1 - a*a/den;
    d.eps_p_imag = -a*a*cf/(b*den);

    d.zeta = (R_ind > 1e-10) ? R_ind/(R_ind + t.R_ohm) : 0;
    d.icp_eff = d.zeta;

    double vgf = vg(g.M, Tg);
    double Gn = 0.25 * ng * vgf;
    d.T_i_N = Gi * g.M * t.Ai * std::sqrt(2*e*t.Vgrid/g.M);
    d.T_n_N = Gn * g.M * vgf * t.Ag;
    d.T_total_N = d.T_i_N + d.T_n_N;
    d.P_RF = 0.5 * (R_ind + t.R_ohm) * I_coil_val * I_coil_val;

    double ve = std::sqrt(2*e*t.Vgrid/g.M);
    double pi_ = 0.5*g.M*ve*ve*Gi*t.Ai;
    double pn_ = 0.5*g.M*vgf*vgf*Gn*t.Ag;
    d.gamma_eff = (d.P_RF > 1e-10) ? (pi_+pn_)/(pi_+pn_+d.P_RF) : 0;
    d.xi_mN_kW = (d.P_RF > 1e-10) ? 1000*d.T_total_N/d.P_RF : 0;
    d.eta_mass = (t.Q0 > 0) ? Gi*t.Ai/t.Q0 : 0;

    return d;
}

void emit_csv_row(std::ofstream& f, const std::string& method, double Q0sccm,
                   double n, double ng, double Te, double Tg,
                   double P_RFG, double P_abs, double R_ind, double I_coil,
                   const SimContext& ctx, const DerivedQuantities& d) {
    f << method << ", " << Q0sccm << ", " << Te << ", " << Tg << ", "
      << std::scientific << n << ", " << ng << ", "
      << std::fixed << d.iondeg << ", " << P_RFG << ", "
      << P_abs << ", " << d.I_extr_mA << ", "
      << d.cf << ", " << R_ind << ", " << I_coil << ", "
      << d.eps_p_real << ", " << d.eps_p_imag << ", "
      << d.u_Bohm << ", " << d.J_i << ", "
      << d.zeta << ", " << d.gamma_eff << ", "
      << d.xi_mN_kW << ", " << d.eta_mass << ", "
      << d.pf << ", " << ctx.thruster.frequency/1e6 << ", "
      << d.T_i_N*1e3 << ", " << d.T_n_N*1e3 << ", " << d.T_total_N*1e3 << ", "
      << d.icp_eff << ", " << d.P_RF << ", "
      << ctx.solver.density_profile_factor * n << ", " << ctx.solver.density_profile_factor << "\n";
}
