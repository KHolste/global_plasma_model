"""
rf_diagnostics.py -- RF-Kopplungsdiagnostik fuer den generischen Plasma-Solver.

Berechnet aus dem Plasma-Zustand und der Triebwerksgeometrie physikalisch
interpretierbare RF-Kopplungsgroessen. Angelehnt an das ICP-Modell im
C++-Solver (Chabert 2012, Eq.14-20), vereinfacht fuer den Python-Pfad.

Annahmen (explizit):
- Zylindrische ICP-Geometrie mit Solenoid-Spule
- Uniform-Plasma-Naeherung (0D)
- Kollisionsdominiertes Regime (nu_m >> omega oder vergleichbar)
- Bessel-Funktionen J0, J1 aus scipy
- Keine Skin-Effekt-Korrektur (duennes Plasma)

Referenz: Chabert et al., Phys. Plasmas 19, 073512 (2012), Eq.14-20
          Grondein et al., Phys. Plasmas 23, 033514 (2016), Eq.14-20
"""
from __future__ import annotations

import math
import cmath
from dataclasses import dataclass

import numpy as np
from scipy.special import jv as bessel_j  # J_n(x) Bessel function

from plasma_chemistry import ThrusterGeometry, KB, CONV, E_CH, ME, PI, EPS0


MU0 = 4 * PI * 1e-7   # Vakuumpermeabilitaet [H/m]
C0 = 299792458.0       # Lichtgeschwindigkeit [m/s]


@dataclass
class RFDiagnostics:
    """RF-Kopplungsdiagnosegroessen."""
    # Primaere Transportgroessen
    nu_m: float = 0.0          # Elektron-Neutral-Kollisionsfrequenz [s^-1]
    omega_pe: float = 0.0      # Plasmafrequenz [rad/s]
    omega: float = 0.0         # RF-Kreisfrequenz [rad/s]
    nu_m_over_omega: float = 0.0  # Normierter Kollisionsparameter [-]

    # Komplexe Plasma-Dielektrizitaet
    eps_p_real: float = 0.0    # Re(epsilon_p) [-]
    eps_p_imag: float = 0.0    # Im(epsilon_p) [-]

    # ICP-Kopplungsgroessen
    R_ind: float = 0.0         # Induktiver Widerstand [Ohm]
    L_ind: float = 0.0         # Induktivitaet [H]
    R_coil: float = 0.0        # Ohmscher Spulenwiderstand [Ohm]
    zeta: float = 0.0          # ICP-Kopplungseffizienz [-]
    P_abs_frac: float = 0.0    # Absorbierte Leistungsfraktion [-]

    # Abgeleitete Groessen
    skin_depth: float = 0.0    # Klassische Eindringtiefe [m]
    ne: float = 0.0            # Elektronendichte [m^-3] (Eingabe)

    def summary(self) -> str:
        return (f"nu_m={self.nu_m:.2e} s-1, nu_m/omega={self.nu_m_over_omega:.3f}, "
                f"eps_p=({self.eps_p_real:.1f}, {self.eps_p_imag:.1f}), "
                f"R_ind={self.R_ind:.4f} Ohm, zeta={self.zeta:.4f}, "
                f"P_abs_frac={self.P_abs_frac:.4f}")


def compute_rf_diagnostics(
    ne: float,
    nu_m: float,
    geom: ThrusterGeometry,
) -> RFDiagnostics:
    """Berechne RF-Kopplungsdiagnostik aus Plasma-Zustand und Geometrie.

    Args:
        ne: Elektronendichte [m^-3]
        nu_m: Elektron-Neutral-Kollisionsfrequenz [s^-1]
        geom: Triebwerksgeometrie

    Returns:
        RFDiagnostics mit allen berechneten Groessen
    """
    d = RFDiagnostics()
    d.ne = ne
    d.nu_m = nu_m
    d.omega = geom.omega
    d.R_coil = geom.R_ohm

    # Normierter Kollisionsparameter
    d.nu_m_over_omega = nu_m / geom.omega if geom.omega > 0 else 0

    # Plasmafrequenz
    if ne > 0:
        d.omega_pe = math.sqrt(ne * E_CH**2 / (ME * EPS0))
    else:
        d.omega_pe = 0.0

    # Komplexe Plasma-Dielektrizitaet (Drude-Modell)
    # eps_p = 1 - omega_pe^2 / (omega^2 + nu_m^2) - j * omega_pe^2 * nu_m / (omega * (omega^2 + nu_m^2))
    omega = geom.omega
    if omega > 0 and d.omega_pe > 0:
        denom = omega**2 + nu_m**2
        d.eps_p_real = 1.0 - d.omega_pe**2 / denom
        d.eps_p_imag = -d.omega_pe**2 * nu_m / (omega * denom)
    else:
        d.eps_p_real = 1.0
        d.eps_p_imag = 0.0

    eps_p = complex(d.eps_p_real, d.eps_p_imag)

    # Wellenzahl und Bessel-Argument
    k0 = omega / C0  # Vakuum-Wellenzahl
    kR = k0 * geom.R  # Argument fuer Bessel-Funktionen

    # Komplexes Bessel-Argument mit Plasma-Permittivitaet
    # k_p = k0 * sqrt(eps_p)
    eps_p_c = complex(d.eps_p_real, d.eps_p_imag)
    k_p_R = k0 * cmath.sqrt(eps_p_c) * geom.R

    # ICP-Kopplungsmodell (Chabert Eq.17-18)
    # R_ind = (2*pi*N^2) / (L*omega*eps0) * Re[ j*k_p*R * J1(k_p*R) / (eps_p * J0(k_p*R)) ]
    # L_ind = L_coil * (1 - R^2/Rc^2 + (2*pi*N^2)/(L*omega^2*eps0) * Im[...])

    N = geom.Nw
    L = geom.L
    Rc = geom.Rc

    # Bessel-Funktionen mit komplexem Argument
    try:
        J0_val = complex(bessel_j(0, k_p_R.real + 1j * k_p_R.imag))
        J1_val = complex(bessel_j(1, k_p_R.real + 1j * k_p_R.imag))
    except (ValueError, OverflowError):
        # Fallback: kleine Argumente
        J0_val = complex(1.0, 0.0)
        J1_val = complex(k_p_R / 2.0)

    # Vermeide Division durch 0
    if abs(J0_val) < 1e-30:
        J0_val = complex(1e-30, 0)
    if abs(eps_p_c) < 1e-30:
        eps_p_c = complex(1e-30, 0)

    # Bessel-Quotient
    bessel_ratio = 1j * k_p_R * J1_val / (eps_p_c * J0_val)

    prefactor = 2 * PI * N**2 / (L * omega * EPS0)

    d.R_ind = max(prefactor * bessel_ratio.real, 1e-6)

    # Spuleninduktivitaet
    L_coil = MU0 * PI * Rc**2 * N**2 / geom.lc
    d.L_ind = L_coil * (1.0 - geom.R**2 / Rc**2) + prefactor / omega * bessel_ratio.imag

    # ICP-Kopplungseffizienz
    d.zeta = d.R_ind / (d.R_ind + d.R_coil) if (d.R_ind + d.R_coil) > 0 else 0

    # Absorbierte Leistungsfraktion (= zeta fuer einfaches Modell)
    d.P_abs_frac = d.zeta

    # Klassische Eindringtiefe (kollisionsdominiert)
    if d.omega_pe > 0 and nu_m > 0:
        d.skin_depth = C0 / d.omega_pe * math.sqrt(2 * omega / nu_m) if omega > 0 else 0
    else:
        d.skin_depth = 0.0

    return d
