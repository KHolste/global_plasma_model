// phys_const.hpp -- Unveraenderliche physikalische Naturkonstanten.
//
// NUR constexpr-Werte. Kein mutabler Zustand.
// Gasspezifische Groessen (M, Eiz, ...) leben im SimContext.
#ifndef PHYS_CONST_HPP
#define PHYS_CONST_HPP

namespace PhysConst {
    constexpr double me       = 9.10938215e-31;
    constexpr double e        = 1.602176487e-19;
    constexpr double kB       = 1.3806504e-23;
    constexpr double epsilon0 = 8.854187e-12;
    constexpr double pi       = 3.141592653589793;
    constexpr double c        = 299792458.0;
    constexpr double mu_0     = 1.256637061e-6;
    constexpr double SCCM_TO_PPS = 4.477962312e17;
    constexpr double Tg0      = 293.00;
    constexpr double conv     = 11604.5250061657;
}

#endif
