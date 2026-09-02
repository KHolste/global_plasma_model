// sim_logging.hpp -- Debug-Ausgabe und Masterlog.
#ifndef SIM_LOGGING_HPP
#define SIM_LOGGING_HPP

#include "sim_context.hpp"
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <algorithm>

// ═══ Debug ═════════════════════════════════════════════════

inline void debug_emit(const SimContext& ctx, int lvl, const std::string& tag,
                        const std::string& msg, bool also_stderr = false) {
    if (ctx.debug_level < lvl) return;
    std::string line = "DBG " + tag + " " + msg;
    std::cout << line << std::endl; std::cout.flush();
    if (also_stderr || lvl <= 1) { std::cerr << line << std::endl; std::cerr.flush(); }
    std::ofstream log("chabert_debug.log", std::ios::app);
    if (log) log << line << '\n';
}

inline std::string make_timestamp() {
    auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    char buf[64]; std::strftime(buf, sizeof(buf), "%Y%m%d_%H%M%S", std::localtime(&t));
    return buf;
}

inline std::string make_timestamp_readable() {
    auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    char buf[64]; std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
    return buf;
}

// ═══ Masterlog (impl in sim_logging.cpp) ═══════════════════

void write_masterlog(const SimContext& ctx, const std::string& config_file,
                      double elapsed, int count_ok, int count_nophys, int count_numfail,
                      const std::vector<SimLogRow>& rows,
                      const std::vector<SimLogEvent>& events);

#endif
