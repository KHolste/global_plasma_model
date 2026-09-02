// rates.cpp -- Ratentabellen-I/O (Cold Path).
#include "rates.hpp"
#include <fstream>
#include <sstream>

static bool load_table_2col(std::vector<KizEntry>& tbl, const std::string& path) {
    std::ifstream f(path); if (!f) return false;
    tbl.clear(); std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == 'T') continue;
        std::istringstream ss(line); KizEntry e{}; char c;
        if (ss >> e.Te_eV >> c >> e.Kiz) tbl.push_back(e);
    }
    return !tbl.empty();
}

bool load_kiz_table(RateConfig& r, const std::string& path) { return load_table_2col(r.kiz, path); }

bool load_kel_table(RateConfig& r, const std::string& path) {
    std::ifstream f(path); if (!f) return false;
    r.kel.clear(); std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == 'T') continue;
        std::istringstream ss(line); KelEntry e{}; char c;
        if (ss >> e.Te_eV >> c >> e.Kel) r.kel.push_back(e);
    }
    return !r.kel.empty();
}

bool load_kex_table(RateConfig& r, const std::string& path) {
    std::ifstream f(path); if (!f) return false;
    r.kex.clear(); std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == 'T') continue;
        std::istringstream ss(line); KexEntry e{}; char c;
        if (ss >> e.Te_eV >> c >> e.Kex_total >> c >> e.Pexc_coeff)
            r.kex.push_back(e);
    }
    return !r.kex.empty();
}
