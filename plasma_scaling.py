"""
plasma_scaling.py -- Leistungsabhaengige Plasma-Skalierung fuer kleine ICP-Thruster.

Korrigiert die 0D-Homogenitaetsannahme fuer kleine Geometrien, in denen
die Elektronendichte und damit das Plasmaangebot am Grid staerker mit
der absorbierten Leistung skaliert als im reinen 0D-Modell.

Physikalische Grundlage:
- Im 0D-Modell: ne(P) ~ P / (Kiz·ng·Eiz·V) -> schwach, da ng sinkt
- In realen kleinen ICPs: Skin-Effekt lokalisiert die Leistungsdeposition
- Lokale Produktionsrate steigt staerker als mittlere Dichte
- Effektiver Korrekturfaktor: f_power = (P_abs/P_ref)^alpha

Parameter:
  P_ref: Referenzleistung bei der f_power = 1 [W]
         Skaliert mit Kammergroesse: P_ref ~ n0·V·Eiz·Kiz(Te0)·uB(Te0)·Aeff
         Typisch: P_ref ~ 1-5 W fuer NPT30-I2, ~ 50-200 W fuer grosse ICP

  alpha: Skalierungsexponent [-]
         alpha = 0: kein Effekt (reines 0D)
         alpha = 0.5: Wurzel-Skalierung (konservativ)
         alpha = 0.7: moderat (empfohlen fuer kleine ICP)
         alpha = 1.0: linear (obere Grenze)

Referenzen:
  Lieberman & Lichtenberg 2005, Ch.10 (ICP power scaling)
  Chabert & Braithwaite 2011, Ch.11 (power balance in ICPs)
  Lafleur et al. 2022 (NPT30-I2 experimental scaling)

Annahme: Maxwell-EEDF, elektropositive Sheath, keine Profilaufloesung.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

KB = 1.3806504e-23
CONV = 11604.5250061657
E_CH = 1.602176487e-19
PI = math.pi


@dataclass
class PlasmaScaling:
    """Leistungsabhaengige Skalierungskorrektur fuer kleine ICP-Thruster."""

    # Konfiguration
    alpha: float = 0.7          # Skalierungsexponent
    P_ref_W: float = 3.0        # Referenzleistung [W]
    enabled: bool = True         # Korrektur aktiv?

    # Diagnostik (nach Berechnung)
    f_power: float = 1.0        # Berechneter Korrekturfaktor
    P_abs_used: float = 0.0     # Verwendete P_abs [W]

    def compute_factor(self, P_abs_W: float) -> float:
        """Berechne den Leistungsskalierungsfaktor.

        f_power = (P_abs / P_ref)^alpha    fuer P_abs >= P_ref
        f_power = P_abs / P_ref            fuer P_abs < P_ref (linear im Startup)

        Returns: f_power >= 0
        """
        if not self.enabled or P_abs_W <= 0:
            self.f_power = 1.0
            self.P_abs_used = P_abs_W
            return 1.0

        P_ratio = P_abs_W / max(self.P_ref_W, 0.1)
        self.P_abs_used = P_abs_W

        if P_ratio <= 1.0:
            # Unterhalb Referenz: linearer Anlauf (vermeidet f < 1)
            self.f_power = P_ratio
        else:
            # Oberhalb Referenz: Power-Law-Skalierung
            self.f_power = P_ratio ** self.alpha

        return self.f_power

    @staticmethod
    def estimate_P_ref(R: float, L: float, Te_eV: float = 4.0,
                       n_ref: float = 1e17, sigma_i: float = 1e-18) -> float:
        """Schaetze P_ref aus Geometrie und typischen Plasmaparametern.

        P_ref ~ n_ref * uB * Aeff * Eiz / V * V
              ~ n_ref * uB * Aeff * Eiz

        Fuer kleine Thruster: P_ref ~ 1-5 W
        Fuer grosse Thruster: P_ref ~ 50-200 W
        """
        M_I = 2.107e-25
        uB = math.sqrt(KB * Te_eV * CONV / M_I)
        V = PI * R**2 * L
        Aeff = 2 * PI * R**2 * 0.5 + 2 * PI * R * L * 0.3  # hL~0.5, hR~0.3
        Eiz = 10.45 * E_CH  # Iod Ionisierungsenergie

        # P_ref als Leistung, die noetig ist um n_ref zu halten
        P_ref = n_ref * uB * Aeff * Eiz
        return P_ref


def apply_scaling_to_extraction(
    I_beam_0D_mA: float,
    P_abs_W: float,
    scaling: PlasmaScaling,
) -> tuple[float, float]:
    """Wende Leistungsskalierung auf den 0D-Beamstrom an.

    Der 0D-Beamstrom wird mit dem Leistungsfaktor korrigiert:
    I_beam_corrected = I_beam_0D * f_power

    Returns: (I_beam_corrected_mA, f_power)
    """
    f = scaling.compute_factor(P_abs_W)
    return I_beam_0D_mA * f, f
