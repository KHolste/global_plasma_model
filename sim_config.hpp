// sim_config.hpp -- Konfigurationsdatei-Parser und Anwendung auf SimContext.
#ifndef SIM_CONFIG_HPP
#define SIM_CONFIG_HPP

#include "sim_context.hpp"
#include <string>
#include <map>

struct ConfigData {
    std::map<std::string, double> numeric;
    std::map<std::string, std::string> strings;
};

ConfigData loadConfig(const std::string& filename);
void applyConfig(SimContext& ctx, const ConfigData& cd);

#endif
