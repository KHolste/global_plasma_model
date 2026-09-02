// sim_logging.cpp -- Masterlog-Writer.
#include "sim_logging.hpp"

using namespace PhysConst;

void write_masterlog(const SimContext& ctx, const std::string& config_file,
                      double elapsed, int count_ok, int count_nophys, int count_numfail,
                      const std::vector<SimLogRow>& rows,
                      const std::vector<SimLogEvent>& events) {
    const auto& t = ctx.thruster;
    const auto& sp = ctx.solver;
    std::string ts = make_timestamp_readable();
    std::string logname = "simulation_log_" + make_timestamp() + ".txt";
    std::ofstream lf(logname);
    if (!lf) return;
    std::cout << "LOG_FILE " << logname << std::endl;

    const char* rmn = (ctx.rate_model==2) ? "Full tabulated" : (ctx.rate_model==1) ? "Conservative tabulated" : "Legacy";
    lf << "==================================================\nGLOBAL PLASMA MODEL - SIMULATION LOG\n"
       << "gas_species:     " << ctx.gas.species << "\ntimestamp_end:   " << ts
       << "\nruntime_seconds: " << std::fixed << std::setprecision(3) << elapsed
       << "\nsolver_mode:     stationary (LM), solve_mode=" << sp.solve_mode
       << "\nrate_model:      " << ctx.rate_model << " (" << rmn
       << ")\nconfig_file:     " << config_file
       << "\n==================================================\n\n";

    lf << "--------------------------------------------------\nSIMULATION PARAMETERS\n--------------------------------------------------\n";
    auto pp = [&](const char* label, double val, const char* unit, const char* key) {
        lf << std::left << std::setw(35) << label << "| " << std::right << std::setw(14)
           << std::scientific << std::setprecision(6) << val << " | " << std::setw(6) << unit << " | " << key << "\n";
    };
    pp("Radius",t.R,"m","R"); pp("Laenge",t.L,"m","L");
    pp("betai",t.betai,"-","betai"); pp("betag",t.betag,"-","betag");
    pp("Frequenz",t.frequency,"Hz","frequency"); pp("Nw",t.Nw,"-","Nw");
    pp("R_ohm",t.R_ohm,"Ohm","R_ohm"); pp("Rc",t.Rc,"m","Rc"); pp("lc",t.lc,"m","lc");
    pp("Vgrid",t.Vgrid,"V","Vgrid"); pp("sgrid",t.sgrid,"m","sgrid");
    pp("P_RFG",t.P_RFG,"W","P_RFG"); pp("P_RFG_max",sp.P_RFG_max,"W","P_RFG_max");
    lf << "--------------------------------------------------\n";
    pp("Q0sccm_start",sp.Q0sccm_start,"sccm","Q0sccm_start");
    pp("Q0sccm_step",sp.Q0sccm_step,"sccm","Q0sccm_step");
    pp("jjmax",(double)sp.jjmax,"-","jjmax");
    pp("I_soll",sp.I_soll,"mA","I_soll");
    lf << "--------------------------------------------------\n";
    pp("V",t.V,"m^3","V"); pp("A",t.A,"m^2","A");
    pp("omega",t.omega,"rad/s","omega");
    pp("Ai",t.Ai,"m^2","Ai"); pp("Ag",t.Ag,"m^2","Ag");
    pp("J_CL",t.J_CL,"A/m^2","J_CL");
    pp("newton_tol",sp.newton_tol,"rel","newton_tol");
    lf << std::left << std::setw(35) << "gas_species" << "| " << ctx.gas.species << "\n";
    lf << std::left << std::setw(35) << "cs_database" << "| " << ctx.cs_database << "\n\n";

    // Ergebnistabelle
    lf << "--------------------------------------------------\nRESULT TABLE\n--------------------------------------------------\n";
    lf << std::left << std::setw(5) << "idx" << "| " << std::setw(9) << "Q0_sccm" << "| "
       << std::setw(22) << "status" << "| " << std::setw(14) << "P_sol_W" << "| "
       << std::setw(10) << "I_mA" << "| " << std::setw(8) << "Te_eV" << "| "
       << std::setw(8) << "Tg_K" << "| " << std::setw(12) << "n_m-3" << "| "
       << std::setw(12) << "ng_m-3" << "| " << std::setw(10) << "resid" << "| note\n";
    lf << std::string(130, '-') << "\n";
    for (const auto& r : rows) {
        lf << std::left << std::setw(5) << r.idx << "| " << std::fixed << std::setprecision(4) << std::setw(9) << r.Q0sccm << "| ";
        lf << std::left << std::setw(22) << r.status << "| ";
        if (r.has_data) {
            lf << std::right << std::fixed << std::setprecision(2) << std::setw(14) << r.P_sol << "| "
               << std::setprecision(3) << std::setw(10) << r.I_mA << "| " << std::setw(8) << r.Te << "| "
               << std::setprecision(1) << std::setw(8) << r.Tg << "| " << std::scientific << std::setprecision(2)
               << std::setw(12) << r.n << "| " << std::setw(12) << r.ng << "| " << std::setprecision(1) << std::setw(10) << r.resid << "| ";
        } else {
            lf << std::right << std::setw(14) << "-" << "| " << std::setw(10) << "-" << "| " << std::setw(8) << "-" << "| "
               << std::setw(8) << "-" << "| " << std::setw(12) << "-" << "| " << std::setw(12) << "-" << "| " << std::setw(10) << "-" << "| ";
        }
        lf << std::left << r.note << "\n";
    }
    lf << "\n";

    // Events
    lf << "--------------------------------------------------\nEVENT DETAILS\n--------------------------------------------------\n";
    int li = -1;
    for (const auto& ev : events) {
        if (ev.idx != li) { lf << "\n[Q0 idx=" << ev.idx << " | Q0_sccm=" << std::fixed << std::setprecision(4) << ev.Q0sccm << "]\n"; li = ev.idx; }
        lf << "  - " << ev.message << "\n";
    }
    if (events.empty()) lf << "(keine Ereignisse)\n";
    lf << "\n";

    // Summary
    lf << "--------------------------------------------------\nSUMMARY\n--------------------------------------------------\n";
    lf << std::left << std::setw(30) << "total_points" << "| " << sp.jjmax << "\n"
       << std::setw(30) << "converged" << "| " << count_ok << "\n"
       << std::setw(30) << "no_physical_solution" << "| " << count_nophys << "\n"
       << std::setw(30) << "numerical_fail" << "| " << count_numfail << "\n"
       << std::setw(30) << "runtime_seconds" << "| " << std::fixed << std::setprecision(3) << elapsed << "\n";
    lf << "--------------------------------------------------\n\n";

    // Machine-readable
    lf << "--------------------------------------------------\nMACHINE-READABLE DATA TABLE\n--------------------------------------------------\n";
    lf << "# CL_LIMIT_J_A_per_m2=" << std::scientific << std::setprecision(6) << t.J_CL << "\n";
    lf << "# CL_LIMIT_I_mA=" << std::fixed << std::setprecision(4) << t.J_CL*t.Ai*1000 << "\n";
    lf << "DATA_HEADER|idx|Q0sccm|status|P_RFG_W|P_abs_W|P_RF_W|I_mA|Te_eV|Tg_K"
       << "|n_m3|ng_m3|iondeg_pct|collision_freq|R_induktiv_Ohm|I_coil_A"
       << "|eps_p_real|eps_p_imag|u_Bohm_ms|J_i_A_m2|plasmafrequenz_rad_s"
       << "|thrust_ions_mN|thrust_atoms_mN|thrust_total_mN"
       << "|icp_power_efficiency|gamma_thrust_eff|xi_mN_per_kW|eta_mass_util"
       << "|n_eff|density_profile_factor|resid|note\n";
    for (const auto& r : rows) {
        lf << "DATA|" << r.idx << "|" << std::fixed << std::setprecision(4) << r.Q0sccm << "|" << r.status << "|";
        if (r.has_data) {
            lf << r.P_sol << "|" << r.P_abs << "|" << r.P_RF << "|" << r.I_mA << "|" << r.Te << "|"
               << std::setprecision(2) << r.Tg << "|" << std::scientific << std::setprecision(6)
               << r.n << "|" << r.ng << "|" << std::fixed << std::setprecision(4) << r.iondeg << "|"
               << std::scientific << r.collision_freq << "|" << r.R_induktiv << "|" << r.I_coil << "|"
               << r.eps_p_real << "|" << r.eps_p_imag << "|" << r.u_Bohm << "|" << r.J_i << "|" << r.plasmafrequenz << "|"
               << std::fixed << std::setprecision(6)
               << r.thrust_ions_mN << "|" << r.thrust_atoms_mN << "|" << r.thrust_total_mN << "|"
               << r.icp_eff << "|" << r.gamma_eff << "|" << r.xi_mN_kW << "|" << r.eta_mass << "|"
               << sp.density_profile_factor*r.n << "|" << sp.density_profile_factor << "|" << std::scientific << r.resid;
        } else { for (int f = 0; f < 28; ++f) lf << "|"; }
        lf << "|" << r.note << "\n";
    }
    lf << "--------------------------------------------------\n";
}
