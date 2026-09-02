// test_extraction.cpp -- Test fuer die gemeinsame Extraktionsrechnung.
//
// Geprueft wird vor allem das, was vorher nicht stimmte:
//   1. Strahlstrom und Schub stammen aus derselben Rechnung. Greift die
//      Raumladungsgrenze, sinken beide im gleichen Verhaeltnis.
//   2. Die Ladungszahl wirkt an allen drei Stellen: Strom je Teilchen,
//      Bohm-Geschwindigkeit und Austrittsgeschwindigkeit.
//   3. Der Massenwirkungsgrad ist ein Massenverhaeltnis und bleibt damit
//      auch bei molekularem Treibstoff sinnvoll.
//
// Build:
//   g++ -O3 -std=c++17 test_extraction.cpp sim_config.o rates.o physics.o \
//       solver.o sim_logging.o bessel_wrapper.o chem_loader.o -o test_extraction
#include "beam_extraction_cpp.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

static int passed = 0, failed = 0;

static void check(const string& name, bool cond, const string& detail = "") {
    if (cond) { ++passed; cout << "  PASS: " << name << endl; }
    else { ++failed; cout << "  FAIL: " << name << " -- " << detail << endl; }
}

static bool nahe(double a, double b, double rel = 1e-9) {
    double skala = max(fabs(a), fabs(b));
    return fabs(a - b) <= rel * max(skala, 1e-300);
}

static const double M_XE = 2.1801711e-25;
static const double M_I  = 2.1071711e-25;
static const double M_I2 = 4.2143422e-25;

static Thruster standard_thruster() {
    Thruster t;
    t.R = 0.02; t.L = 0.04; t.betai = 0.5; t.betag = 0.05145;
    t.Vgrid = 1500.0; t.sgrid = 0.001; t.eta_opt = 1.0;
    t.Q0sccm = 0.40;
    t.recompute(M_XE);
    return t;
}

static ExtractionInput xenon_input(double n_ion, double n_gas) {
    ExtractionInput in;
    in.ions.push_back({"Xe+", n_ion, M_XE, 1});
    in.neutrals.push_back({"Xe", n_gas, M_XE});
    in.Te_eV = 4.0;
    in.Tg_K = 300.0;
    in.sigma_i = 1e-18;
    in.mdot_in_kg_s = 0;
    return in;
}

int main() {
    const Thruster t = standard_thruster();
    const double A_open = t.betai * PhysConst::pi * t.R * t.R;

    // ── Test 1: eine Ionensorte, von Hand nachgerechnet ──────
    cout << "\n--- Test 1: eine Ionensorte ---" << endl;
    {
        ExtractionInput in = xenon_input(2.0e17, 2.2e19);
        ExtractionResult r = compute_extraction(in, t);

        double lam = 1.0 / (2.2e19 * 1e-18);
        double hL = 0.86 * pow(3.0 + t.L/(2*lam), -0.5);
        double uB = sqrt(PhysConst::kB * 4.0 * PhysConst::conv / M_XE);
        double J = PhysConst::e * hL * 2.0e17 * uB;
        double I = J * A_open * 1000.0;
        double v = sqrt(2 * PhysConst::e * t.Vgrid / M_XE);
        double mdot = (J * A_open / PhysConst::e) * M_XE;

        check("Randschichtfaktor", nahe(r.hL, hL), to_string(r.hL));
        check("Bohm-Geschwindigkeit", nahe(r.ions[0].u_B, uB), to_string(r.ions[0].u_B));
        check("Stromdichte", nahe(r.J_Bohm_total, J), to_string(r.J_Bohm_total));
        check("Strahlstrom", nahe(r.I_beam_mA, I), to_string(r.I_beam_mA));
        check("Austrittsgeschwindigkeit", nahe(r.ions[0].v_exhaust, v));
        check("Ionenschub", nahe(r.thrust_ions_N, mdot * v), to_string(r.thrust_ions_N));
        check("plasmabegrenzt", r.limiting == "plasma", r.limiting);
        check("keine Drosselung", nahe(r.throttle, 1.0));
    }

    // ── Test 2: Raumladungsgrenze wirkt auf alles gleich ─────
    // Das ist der Kern der Reparatur: frueher kam der Strahlstrom aus dem
    // begrenzten Modell, der Schub aus dem unbegrenzten Bohm-Fluss.
    cout << "\n--- Test 2: Raumladungsgrenze ---" << endl;
    {
        Thruster eng = t;
        eng.Vgrid = 60.0;      // niedrige Spannung, weiter Spalt
        eng.sgrid = 0.004;
        eng.recompute(M_XE);

        ExtractionInput in = xenon_input(5.0e17, 2.2e19);
        ExtractionResult frei = compute_extraction(in, t);
        ExtractionResult begrenzt = compute_extraction(in, eng);

        check("jetzt raumladungsbegrenzt", begrenzt.limiting == "space_charge",
              begrenzt.limiting);
        check("Verhaeltnis groesser eins", begrenzt.space_charge_ratio > 1.0,
              to_string(begrenzt.space_charge_ratio));
        check("Drosselung ist der Kehrwert",
              nahe(begrenzt.throttle, 1.0/begrenzt.space_charge_ratio));

        // Strom und ausgetragener Massenstrom sinken um denselben Faktor.
        double f_I = begrenzt.I_beam_mA / (frei.I_beam_mA);
        double f_m = begrenzt.mdot_out_kg_s / frei.mdot_out_kg_s;
        check("Strom und Massenstrom sinken gleich", nahe(f_I, f_m, 1e-12),
              to_string(f_I) + " vs " + to_string(f_m));
        check("Strom sinkt um die Drosselung", nahe(f_I, begrenzt.throttle),
              to_string(f_I));

        // Der Schub folgt dem Strom, weil beide aus demselben Teilchenstrom
        // kommen; die Austrittsgeschwindigkeit haengt nur an der Spannung.
        double v_eng = sqrt(2 * PhysConst::e * eng.Vgrid / M_XE);
        check("Schub aus demselben Teilchenstrom",
              nahe(begrenzt.thrust_ions_N, begrenzt.mdot_out_kg_s * v_eng),
              to_string(begrenzt.thrust_ions_N));
        check("Grenzstrom wird ausgewiesen",
              nahe(begrenzt.I_CL_limit_mA, begrenzt.I_beam_mA, 1e-12),
              to_string(begrenzt.I_CL_limit_mA) + " vs " + to_string(begrenzt.I_beam_mA));
    }

    // ── Test 3: Ladungszahl wirkt an allen drei Stellen ──────
    cout << "\n--- Test 3: zweifach geladene Ionen ---" << endl;
    {
        ExtractionInput einfach = xenon_input(1.0e17, 2.2e19);
        ExtractionInput zweifach = xenon_input(1.0e17, 2.2e19);
        zweifach.ions[0] = {"Xe2+", 1.0e17, M_XE, 2};

        ExtractionResult r1 = compute_extraction(einfach, t);
        ExtractionResult r2 = compute_extraction(zweifach, t);

        // u_B ~ sqrt(Z), Strom ~ Z * u_B ~ Z^1.5, v_aus ~ sqrt(Z),
        // Teilchenstrom ~ Z^1.5 / Z = sqrt(Z), Schub ~ sqrt(Z) * sqrt(Z) = Z
        check("Bohm-Geschwindigkeit mal Wurzel zwei",
              nahe(r2.ions[0].u_B, r1.ions[0].u_B * sqrt(2.0)),
              to_string(r2.ions[0].u_B / r1.ions[0].u_B));
        check("Strom mal zwei hoch drei halbe",
              nahe(r2.I_beam_mA, r1.I_beam_mA * pow(2.0, 1.5)),
              to_string(r2.I_beam_mA / r1.I_beam_mA));
        check("Austrittsgeschwindigkeit mal Wurzel zwei",
              nahe(r2.ions[0].v_exhaust, r1.ions[0].v_exhaust * sqrt(2.0)));
        check("Massenstrom mal Wurzel zwei",
              nahe(r2.mdot_out_kg_s, r1.mdot_out_kg_s * sqrt(2.0)),
              to_string(r2.mdot_out_kg_s / r1.mdot_out_kg_s));
        check("Schub mal zwei",
              nahe(r2.thrust_ions_N, r1.thrust_ions_N * 2.0),
              to_string(r2.thrust_ions_N / r1.thrust_ions_N));
    }

    // ── Test 4: Gemisch aus beiden Ladungszustaenden ─────────
    cout << "\n--- Test 4: Gemisch ---" << endl;
    {
        ExtractionInput mix;
        mix.ions.push_back({"Xe+",  8.0e16, M_XE, 1});
        mix.ions.push_back({"Xe2+", 2.0e16, M_XE, 2});
        mix.neutrals.push_back({"Xe", 2.2e19, M_XE});
        mix.Te_eV = 4.0; mix.Tg_K = 300.0; mix.sigma_i = 1e-18;
        mix.mdot_in_kg_s = 0;

        ExtractionInput nur1 = mix; nur1.ions.resize(1);
        ExtractionInput nur2 = mix; nur2.ions.erase(nur2.ions.begin());

        ExtractionResult rm = compute_extraction(mix, t);
        ExtractionResult r1 = compute_extraction(nur1, t);
        ExtractionResult r2 = compute_extraction(nur2, t);

        check("zwei Anteile ausgewiesen", rm.ions.size() == 2, to_string(rm.ions.size()));
        check("Strom ist die Summe", nahe(rm.I_beam_mA, r1.I_beam_mA + r2.I_beam_mA, 1e-12),
              to_string(rm.I_beam_mA));
        check("Schub ist die Summe",
              nahe(rm.thrust_ions_N, r1.thrust_ions_N + r2.thrust_ions_N, 1e-12));
        check("zweifach geladen traegt mehr Strom als Dichte",
              rm.ions[1].I_mA / rm.I_beam_mA > 2.0e16 / 1.0e17,
              to_string(rm.ions[1].I_mA / rm.I_beam_mA));
    }

    // ── Test 5: Massenwirkungsgrad bei molekularem Treibstoff ─
    cout << "\n--- Test 5: Massenwirkungsgrad ---" << endl;
    {
        // Zugefuehrt wird I2, extrahiert wird I+. Ein Teilchenverhaeltnis
        // waere hier sinnlos, das Massenverhaeltnis nicht.
        ExtractionInput in;
        in.ions.push_back({"I+", 2.0e17, M_I, 1});
        in.neutrals.push_back({"I2", 1.0e19, M_I2});
        in.neutrals.push_back({"I", 5.0e18, M_I});
        in.Te_eV = 4.0; in.Tg_K = 350.0; in.sigma_i = 1e-18;
        const double Q0 = 1.0e18;              // zugefuehrte Molekuele je Sekunde
        in.mdot_in_kg_s = Q0 * M_I2;

        ExtractionResult r = compute_extraction(in, t);
        check("Massenwirkungsgrad ist ein Massenverhaeltnis",
              nahe(r.eta_mass, r.mdot_out_kg_s / in.mdot_in_kg_s),
              to_string(r.eta_mass));
        check("Wirkungsgrad zwischen null und eins",
              r.eta_mass > 0 && r.eta_mass < 1.0, to_string(r.eta_mass));
        check("beide Neutralsorten schieben mit", r.thrust_neutrals_N > 0,
              to_string(r.thrust_neutrals_N));
        check("Gesamtschub ist die Summe",
              nahe(r.thrust_total_N, r.thrust_ions_N + r.thrust_neutrals_N));
    }

    // ── Test 6: Gitteroptik wirkt auf Strom und Schub gleich ─
    cout << "\n--- Test 6: Gitteroptik ---" << endl;
    {
        Thruster halb = t;
        halb.eta_opt = 0.5;
        halb.recompute(M_XE);
        ExtractionInput in = xenon_input(2.0e17, 2.2e19);
        ExtractionResult voll = compute_extraction(in, t);
        ExtractionResult halbiert = compute_extraction(in, halb);

        check("Strom halbiert", nahe(halbiert.I_beam_mA, 0.5*voll.I_beam_mA));
        check("Schub halbiert", nahe(halbiert.thrust_ions_N, 0.5*voll.thrust_ions_N));
        check("Massenstrom halbiert", nahe(halbiert.mdot_out_kg_s, 0.5*voll.mdot_out_kg_s));
    }

    // ── Test 7: Randfaelle ───────────────────────────────────
    cout << "\n--- Test 7: Randfaelle ---" << endl;
    {
        ExtractionInput leer;
        leer.Te_eV = 4.0; leer.Tg_K = 300.0;
        ExtractionResult r = compute_extraction(leer, t);
        check("ohne Ionen kein Strom", nahe(r.I_beam_mA, 0.0));
        check("ohne Ionen kein Schub", nahe(r.thrust_total_N, 0.0));

        ExtractionInput kalt = xenon_input(2.0e17, 2.2e19);
        kalt.Te_eV = 0.0;
        ExtractionResult rk = compute_extraction(kalt, t);
        check("ohne Elektronentemperatur kein Strom", nahe(rk.I_beam_mA, 0.0));
    }

    cout << "\n========================================" << endl;
    cout << "  Bestanden: " << passed << "   Fehlgeschlagen: " << failed << endl;
    cout << "========================================" << endl;
    return failed == 0 ? 0 : 1;
}
