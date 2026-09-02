"""
plasma_chemistry.py -- Generisches Chemie-Framework fuer 0D-Globalmodelle.

Definiert Spezies, Reaktionen und Chemie-Pakete als datengetriebene Strukturen.
Der Bilanz-Assembler baut aus einem Chemiepaket automatisch den Zustandsvektor
und die Residualfunktion auf.

Architektur-Annahmen (Version 1):
- Eine Elektronentemperatur Te (Maxwell-EEDF)
- Eine gemeinsame Gastemperatur Tg fuer alle schweren Spezies
- Elektropositive Sheath-Annahme (erweiterbar)
- Keine Vibrationskinetik, keine Metastabilen-Bilanz (architektonisch vorbereitet)

Spezies-Klassen:
  electron, neutral_atom, neutral_molecule, positive_ion, negative_ion

Reaktionstypen:
  ionization, molecular_ionization, dissociation, dissociative_ionization,
  dissociative_attachment, excitation, elastic, recombination,
  charge_exchange, surface_recombination

Referenz: Grondein et al., Phys. Plasmas 23, 033514 (2016)
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional

import numpy as np


# ═══════════════════════════════════════════════════════════════
# Spezies
# ═══════════════════════════════════════════════════════════════

class SpeciesType(str, Enum):
    ELECTRON = "electron"
    NEUTRAL_ATOM = "neutral_atom"
    NEUTRAL_MOLECULE = "neutral_molecule"
    POSITIVE_ION = "positive_ion"
    NEGATIVE_ION = "negative_ion"


@dataclass
class Species:
    """Eine Plasma-Spezies."""
    id: str                          # z.B. "Xe", "I2", "I+", "e"
    name: str                        # Anzeigename
    species_type: SpeciesType
    mass_kg: float                   # Teilchenmasse [kg]
    charge: int = 0                  # Ladung in Elementarladungen
    is_feedstock: bool = False       # Wird von aussen zugefuehrt?
    is_beam_extracted: bool = False  # Wird ueber Grid extrahiert?
    thermal_conductivity: float = 0.0  # kappa [W/(m*K)], fuer Neutrale

    @property
    def is_electron(self) -> bool:
        return self.species_type == SpeciesType.ELECTRON

    @property
    def is_neutral(self) -> bool:
        return self.species_type in (SpeciesType.NEUTRAL_ATOM, SpeciesType.NEUTRAL_MOLECULE)

    @property
    def is_ion(self) -> bool:
        return self.species_type in (SpeciesType.POSITIVE_ION, SpeciesType.NEGATIVE_ION)

    @property
    def is_positive_ion(self) -> bool:
        return self.species_type == SpeciesType.POSITIVE_ION

    @property
    def is_negative_ion(self) -> bool:
        return self.species_type == SpeciesType.NEGATIVE_ION

    @staticmethod
    def from_dict(d: dict) -> Species:
        return Species(
            id=d["id"], name=d.get("name", d["id"]),
            species_type=SpeciesType(d["type"]),
            mass_kg=float(d["mass_kg"]), charge=int(d.get("charge", 0)),
            is_feedstock=d.get("is_feedstock", False),
            is_beam_extracted=d.get("is_beam_extracted", False),
            thermal_conductivity=float(d.get("thermal_conductivity", 0.0)),
        )


# ═══════════════════════════════════════════════════════════════
# Reaktionen
# ═══════════════════════════════════════════════════════════════

class ReactionType(str, Enum):
    IONIZATION = "ionization"
    MOLECULAR_IONIZATION = "molecular_ionization"
    DISSOCIATION = "dissociation"
    DISSOCIATIVE_IONIZATION = "dissociative_ionization"
    DISSOCIATIVE_ATTACHMENT = "dissociative_attachment"
    EXCITATION = "excitation"
    ELASTIC = "elastic"
    RECOMBINATION = "recombination"
    CHARGE_EXCHANGE = "charge_exchange"
    SURFACE_RECOMBINATION = "surface_recombination"


class RateModel(str, Enum):
    """Wie K(Te) bestimmt wird."""
    ARRHENIUS = "arrhenius"         # K = A * exp(-E_a / Te)
    CONSTANT = "constant"           # K = const
    TABULATED = "tabulated"         # Lookup aus CSV
    POLYNOMIAL = "polynomial"       # Polynomfit


@dataclass
class RateCoefficient:
    """Ratenkoeffizient K(Te) [m^3/s] oder K [m^3/s] (konstant)."""
    model: RateModel
    value: float = 0.0              # Fuer CONSTANT
    A: float = 0.0                  # Praeexponentiell (ARRHENIUS)
    E_a_eV: float = 0.0            # Aktivierungsenergie [eV] (ARRHENIUS)
    table_file: str = ""            # CSV-Pfad (TABULATED)
    poly_coeffs: list[float] = field(default_factory=list)  # POLYNOMIAL

    # Geladene Tabelle (lazy)
    _table_Te: Optional[np.ndarray] = field(default=None, repr=False)
    _table_K: Optional[np.ndarray] = field(default=None, repr=False)

    def evaluate(self, Te_eV: float) -> float:
        """Berechne K(Te) [m^3/s]."""
        if self.model == RateModel.CONSTANT:
            return self.value
        elif self.model == RateModel.ARRHENIUS:
            if Te_eV <= 0:
                return 0.0
            return self.A * math.exp(-self.E_a_eV / Te_eV)
        elif self.model == RateModel.TABULATED:
            if self._table_Te is None:
                return 0.0
            return float(np.interp(Te_eV, self._table_Te, self._table_K))
        elif self.model == RateModel.POLYNOMIAL:
            return sum(c * Te_eV**i for i, c in enumerate(self.poly_coeffs))
        return 0.0

    def load_table(self, base_dir: Path):
        """Lade Tabelle aus CSV."""
        if self.model != RateModel.TABULATED or not self.table_file:
            return
        path = base_dir / self.table_file
        if not path.exists():
            return
        Te_vals, K_vals = [], []
        with open(path) as f:
            for line in f:
                if line.startswith("#") or line.startswith("Te"):
                    continue
                parts = line.strip().split(",")
                if len(parts) >= 2:
                    Te_vals.append(float(parts[0]))
                    K_vals.append(float(parts[1]))
        if Te_vals:
            self._table_Te = np.array(Te_vals)
            self._table_K = np.array(K_vals)

    @staticmethod
    def from_dict(d: dict) -> RateCoefficient:
        model = RateModel(d["model"])
        rc = RateCoefficient(model=model)
        if model == RateModel.CONSTANT:
            rc.value = float(d["value"])
        elif model == RateModel.ARRHENIUS:
            rc.A = float(d["A"])
            rc.E_a_eV = float(d["E_a_eV"])
        elif model == RateModel.TABULATED:
            rc.table_file = d.get("file", "")
        elif model == RateModel.POLYNOMIAL:
            rc.poly_coeffs = [float(c) for c in d["coeffs"]]
        return rc


@dataclass
class Reaction:
    """Eine Plasmareaktion."""
    id: str                          # z.B. "iz_Xe", "diss_I2"
    name: str                        # Anzeigename
    reaction_type: ReactionType
    reactants: dict[str, int]        # {species_id: stoechiometrie} (verbraucht)
    products: dict[str, int]         # {species_id: stoechiometrie} (erzeugt)
    rate: RateCoefficient
    energy_eV: float = 0.0          # Energieverlust pro Ereignis [eV]
    is_electron_impact: bool = True  # Benoetiigt Elektron als Reaktionspartner?

    # Flags fuer spezielle Behandlung
    contributes_to_elastic_heating: bool = False  # Beitrag zur Gasheizung
    contributes_to_nu_m: bool = False  # Beitrag zur Kollisionsfrequenz
    surface_gamma: float = 0.0        # Oberflaechenkoeffizient

    def net_stoichiometry(self) -> dict[str, int]:
        """Netto-Stoechiometrie (Produkte - Edukte)."""
        net = {}
        for sp, n in self.reactants.items():
            net[sp] = net.get(sp, 0) - n
        for sp, n in self.products.items():
            net[sp] = net.get(sp, 0) + n
        return {k: v for k, v in net.items() if v != 0}

    @staticmethod
    def from_dict(d: dict) -> Reaction:
        return Reaction(
            id=d["id"], name=d.get("name", d["id"]),
            reaction_type=ReactionType(d["type"]),
            reactants=d.get("reactants", {}),
            products=d.get("products", {}),
            rate=RateCoefficient.from_dict(d["rate"]),
            energy_eV=float(d.get("energy_eV", 0.0)),
            is_electron_impact=d.get("is_electron_impact", True),
            contributes_to_elastic_heating=d.get("elastic_heating", False),
            contributes_to_nu_m=d.get("nu_m", False),
            surface_gamma=float(d.get("surface_gamma", 0.0)),
        )


# ═══════════════════════════════════════════════════════════════
# Chemie-Paket
# ═══════════════════════════════════════════════════════════════

@dataclass
class ChemistryPackage:
    """Vollstaendige Chemie-Definition fuer ein Gas oder Gasgemisch."""
    name: str
    description: str = ""
    species: dict[str, Species] = field(default_factory=dict)
    reactions: list[Reaction] = field(default_factory=list)

    # Physikalische Konstanten
    wall_temperature_K: float = 293.0
    sigma_i: float = 1e-18            # Ion-neutral cross-section [m^2]

    def validate(self) -> list[str]:
        """Pruefe Konsistenz. Gibt Liste von Fehlern zurueck."""
        errors = []
        sp_ids = set(self.species.keys())

        # Muss Elektronen enthalten
        electrons = [s for s in self.species.values() if s.is_electron]
        if len(electrons) != 1:
            errors.append(f"Genau 1 Elektron-Spezies erwartet, {len(electrons)} gefunden")

        # Muss mindestens ein Neutral und ein Ion enthalten
        neutrals = [s for s in self.species.values() if s.is_neutral]
        ions = [s for s in self.species.values() if s.is_positive_ion]
        if not neutrals:
            errors.append("Mindestens eine neutrale Spezies erforderlich")
        if not ions:
            errors.append("Mindestens ein positives Ion erforderlich")

        # Feedstock-Spezies?
        feedstocks = [s for s in self.species.values() if s.is_feedstock]
        if not feedstocks:
            errors.append("Mindestens eine Feedstock-Spezies erforderlich")

        # Reaktionen: alle Spezies muessen definiert sein
        for rxn in self.reactions:
            for sp in list(rxn.reactants.keys()) + list(rxn.products.keys()):
                if sp not in sp_ids and sp != "e":
                    errors.append(f"Reaktion '{rxn.id}': Spezies '{sp}' nicht definiert")

        return errors

    @property
    def heavy_species(self) -> list[Species]:
        """Alle Spezies ausser Elektronen (Zustandsvariablen)."""
        return [s for s in self.species.values() if not s.is_electron]

    @property
    def neutral_species(self) -> list[Species]:
        return [s for s in self.species.values() if s.is_neutral]

    @property
    def ion_species(self) -> list[Species]:
        return [s for s in self.species.values() if s.is_ion]

    @property
    def positive_ions(self) -> list[Species]:
        return [s for s in self.species.values() if s.is_positive_ion]

    @property
    def negative_ions(self) -> list[Species]:
        return [s for s in self.species.values() if s.is_negative_ion]

    @property
    def feedstock_species(self) -> list[Species]:
        return [s for s in self.species.values() if s.is_feedstock]

    def load_rate_tables(self, base_dir: Path):
        """Lade alle tabulierten Raten."""
        for rxn in self.reactions:
            rxn.rate.load_table(base_dir)

    @staticmethod
    def from_json(path: Path) -> ChemistryPackage:
        """Lade Chemiepaket aus JSON-Datei."""
        with open(path, encoding="utf-8") as f:
            d = json.load(f)

        pkg = ChemistryPackage(
            name=d["name"],
            description=d.get("description", ""),
            wall_temperature_K=float(d.get("wall_temperature_K", 293.0)),
            sigma_i=float(d.get("sigma_i", 1e-18)),
        )

        for sd in d["species"]:
            sp = Species.from_dict(sd)
            pkg.species[sp.id] = sp

        for rd in d["reactions"]:
            pkg.reactions.append(Reaction.from_dict(rd))

        return pkg


# ═══════════════════════════════════════════════════════════════
# Bilanz-Assembler
# ═══════════════════════════════════════════════════════════════

# Physikalische Konstanten
ME = 9.10938215e-31
E_CH = 1.602176487e-19
KB = 1.3806504e-23
EPS0 = 8.854187817e-12
PI = math.pi
CONV = 11604.5250061657  # eV -> K


@dataclass
class ThrusterGeometry:
    """Triebwerksgeometrie."""
    R: float = 0.02       # Kammerradius [m]
    L: float = 0.04       # Kammerlaenge [m]
    betai: float = 0.5    # Ionentransparenz
    betag: float = 0.05   # Gastransparenz
    Vgrid: float = 1500.0 # Gitterspannung [V]
    sgrid: float = 0.001  # Gitterabstand [m]
    frequency: float = 2.5e6  # RF-Frequenz [Hz]
    Nw: float = 6.0       # Windungen
    R_ohm: float = 0.36   # Spulenwiderstand
    Rc: float = 0.02      # Spulenradius
    lc: float = 0.04      # Spulenlaenge
    eta_opt: float = 1.0  # Grid-Optik-Effizienz [-] (1.0=Legacy, <1 bei explizitem Preset)

    @property
    def V(self): return PI * self.R**2 * self.L

    @property
    def A(self): return 2*PI*self.R**2 + 2*PI*self.R*self.L

    @property
    def Ai(self): return self.betai * PI * self.R**2

    @property
    def Ag(self): return self.betag * PI * self.R**2

    @property
    def omega(self): return 2 * PI * self.frequency

    @property
    def lambda_0(self): return self.R / 2.405 + self.L / PI


class BalanceAssembler:
    """Baut aus einem ChemistryPackage die Residualfunktion zusammen.

    Zustandsvektor: [n_species_0, n_species_1, ..., Te, Tg]
    Alle schweren Spezies (nicht Elektron) bilden den Dichtenteil.
    Te und Tg sind die letzten zwei Eintraege.

    Elektrondichte wird aus Quasineutralitaet berechnet:
      n_e = sum(Z_i * n_i) fuer alle Ionen
    """

    def __init__(self, chem: ChemistryPackage, geom: ThrusterGeometry):
        self.chem = chem
        self.geom = geom

        # Zustandsvektor-Index
        self._heavy = chem.heavy_species
        self._n_heavy = len(self._heavy)
        self._idx = {sp.id: i for i, sp in enumerate(self._heavy)}
        self._Te_idx = self._n_heavy
        self._Tg_idx = self._n_heavy + 1
        self.state_size = self._n_heavy + 2

        # State labels
        self.state_labels = [sp.id for sp in self._heavy] + ["Te", "Tg"]

    def default_state(self, Q0_pps: float = 2e17) -> np.ndarray:
        """Erzeuge einen plausiblen Startzustand."""
        x = np.zeros(self.state_size)
        for sp in self._heavy:
            i = self._idx[sp.id]
            if sp.is_feedstock:
                x[i] = Q0_pps * 4 / (self.geom.V * self._mean_speed(sp, 300))
            elif sp.is_positive_ion:
                x[i] = 1e16
            elif sp.is_negative_ion:
                x[i] = 1e14
            else:
                x[i] = 1e16
        x[self._Te_idx] = 3.0
        x[self._Tg_idx] = 300.0
        return x

    def electron_density(self, state: np.ndarray) -> float:
        """Quasineutralitaet: n_e = sum(|Z| * n) fuer positive Ionen."""
        ne = 0.0
        for sp in self._heavy:
            if sp.is_positive_ion:
                ne += abs(sp.charge) * state[self._idx[sp.id]]
            elif sp.is_negative_ion:
                ne -= abs(sp.charge) * state[self._idx[sp.id]]
        return max(ne, 0.0)

    def residual(self, state: np.ndarray, P_abs_V: float, Q0_pps: float,
                 alpha_e_wall: float = 7.0, density_profile_factor: float = 1.0
                 ) -> np.ndarray:
        """Berechne Residualvektor.

        state: [n_sp0, ..., n_spN, Te, Tg]
        P_abs_V: absorbierte Leistungsdichte [W/m^3]
        Q0_pps: Teilchenfluss [particles/s]

        Returns: Residualvektor gleicher Groesse
        """
        Te = state[self._Te_idx]
        Tg = state[self._Tg_idx]
        ne = self.electron_density(state)
        V = self.geom.V
        geom = self.geom

        resid = np.zeros(self.state_size)

        # ── Spezies-Bilanzen ─────────────────────────────────
        for rxn in self.chem.reactions:
            K = rxn.rate.evaluate(Te)
            if K <= 0:
                continue

            # Reaktionsrate R = K * Produkt(n_reactant^stoech)
            rate = K
            for sp_id, stoech in rxn.reactants.items():
                if sp_id == "e":
                    rate *= ne ** stoech
                elif sp_id in self._idx:
                    rate *= state[self._idx[sp_id]] ** stoech

            # Netto-Stoechiometrie -> Produktions-/Verbrauchsterme
            net = rxn.net_stoichiometry()
            for sp_id, delta in net.items():
                if sp_id == "e":
                    continue  # Elektronen aus Quasineutralitaet
                if sp_id in self._idx:
                    resid[self._idx[sp_id]] += delta * rate

            # Energieverlust -> Elektronenenergiebilanz
            if rxn.energy_eV > 0 and rxn.is_electron_impact:
                E_loss_J = rxn.energy_eV * E_CH
                n_eff = ne * density_profile_factor
                # Rate fuer Energiebilanz: ne_eff * n_neutral * K
                rate_e = K
                for sp_id, stoech in rxn.reactants.items():
                    if sp_id == "e":
                        rate_e *= n_eff ** stoech
                    elif sp_id in self._idx:
                        rate_e *= state[self._idx[sp_id]] ** stoech
                resid[self._Te_idx] -= E_loss_J * rate_e

            # Elastische Heizung -> Gastemperatur
            if rxn.contributes_to_elastic_heating:
                for sp_id in rxn.reactants:
                    if sp_id != "e" and sp_id in self._idx:
                        sp = self.chem.species[sp_id]
                        n_sp = state[self._idx[sp_id]]
                        Pg = 3 * ME / sp.mass_kg * KB * (Te * CONV - Tg) * ne * n_sp * K
                        resid[self._Tg_idx] += Pg

        # ── Massenfluss-Zufuhr (Feedstock) ────────────────────
        for sp in self.chem.feedstock_species:
            resid[self._idx[sp.id]] += Q0_pps / V

        # ── Oberflaechen-Prozesse (z.B. Wandrekombination) ────
        # Surface reactions: rate = gamma * 0.25 * v_th * n * A_wall / V
        for rxn in self.chem.reactions:
            if rxn.surface_gamma > 0:
                for sp_id in rxn.reactants:
                    if sp_id in self._idx:
                        sp = self.chem.species[sp_id]
                        n_sp = state[self._idx[sp_id]]
                        v_th = self._mean_speed(sp, Tg)
                        # Oberflaechenrate: gamma * (1/4) * n * v * A_n / V
                        # A_n = gesamte Neutralverlustflaeche (Waende + Gitter)
                        A_n = geom.A  # Vollstaendige Wandflaeche
                        surf_rate = rxn.surface_gamma * 0.25 * v_th * n_sp * A_n / V
                        resid[self._idx[sp_id]] -= surf_rate
                        # Produkte
                        for prod_id, stoech in rxn.products.items():
                            if prod_id in self._idx:
                                # Stoechiometrie: 2 I -> 1 I2, also pro I-Verlust 0.5 I2 produziert
                                resid[self._idx[prod_id]] += stoech * surf_rate

        # ── Wandverluste (Ionen + Neutrale) ───────────────────
        sigma_i = self.chem.sigma_i
        total_neutral_density = sum(
            state[self._idx[sp.id]] for sp in self.chem.neutral_species
        )

        for sp in self._heavy:
            idx = self._idx[sp.id]
            n_sp = state[idx]

            if sp.is_positive_ion:
                uB = math.sqrt(KB * Te * CONV / sp.mass_kg)
                lam = 1.0 / max(total_neutral_density * sigma_i, 1e-10)
                hL = 0.86 * (3 + geom.L / (2 * lam))**(-0.5)
                hR = 0.80 * (4 + geom.R / lam)**(-0.5)
                Aeff = 2*hR*PI*geom.R*geom.L + 2*hL*PI*geom.R**2
                wall_loss = n_sp * uB * Aeff / V
                resid[idx] -= wall_loss

                # Ion neutralisation at wall -> returns neutral atom
                # I+ -> I, I2+ -> I2 (simplified: product is atom with same mass)
                for neut in self.chem.neutral_species:
                    if abs(neut.mass_kg - sp.mass_kg) / sp.mass_kg < 0.01:
                        if neut.id in self._idx:
                            resid[self._idx[neut.id]] += wall_loss
                        break

                # Wandverlust-Energie (Elektronen verlieren Energie an Wand)
                resid[self._Te_idx] -= alpha_e_wall * KB * Te * CONV * wall_loss

            elif sp.is_neutral:
                # Alle Neutrale: Austritt durch Gitter (nicht nur Feedstock)
                v_mean = self._mean_speed(sp, Tg)
                outflow = 0.25 * n_sp * v_mean * geom.Ag / V
                resid[idx] -= outflow

        # ── RF-Leistungseintrag ───────────────────────────────
        resid[self._Te_idx] += P_abs_V

        # ── Gasheizung: Waermeleitung ─────────────────────────
        # Effektive kappa gewichtet nach Molenbruch
        n_total = max(total_neutral_density, 1e10)
        kappa_eff = 0.0
        for sp in self.chem.neutral_species:
            if sp.thermal_conductivity > 0 and sp.id in self._idx:
                x_frac = state[self._idx[sp.id]] / n_total
                kappa_eff += sp.thermal_conductivity * x_frac
        if kappa_eff > 0:
            Pg_cond = kappa_eff * (Tg - self.chem.wall_temperature_K) / geom.lambda_0 * geom.A / V
            resid[self._Tg_idx] -= Pg_cond

        return resid

    def _mean_speed(self, sp: Species, T_K: float) -> float:
        """Mittlere thermische Geschwindigkeit [m/s]."""
        return math.sqrt(8 * KB * T_K / (PI * sp.mass_kg))

    def collision_frequency(self, state: np.ndarray) -> float:
        """Effektive Elektron-Neutral-Stossfrequenz fuer RF-Kopplung.
        nu_m = sum_j (K_el_j * n_j) ueber alle Neutrale mit nu_m-Flag."""
        Te = state[self._Te_idx]
        nu = 0.0
        for rxn in self.chem.reactions:
            if rxn.contributes_to_nu_m:
                K = rxn.rate.evaluate(Te)
                for sp_id in rxn.reactants:
                    if sp_id != "e" and sp_id in self._idx:
                        nu += K * state[self._idx[sp_id]]
        return nu


# ═══════════════════════════════════════════════════════════════
# Convenience: Chemiepaket laden
# ═══════════════════════════════════════════════════════════════

def load_chemistry(path: str | Path) -> ChemistryPackage:
    """Lade und validiere ein Chemiepaket."""
    path = Path(path)
    pkg = ChemistryPackage.from_json(path)
    errors = pkg.validate()
    if errors:
        raise ValueError(f"Chemiepaket '{pkg.name}' ungueltig:\n" +
                          "\n".join(f"  - {e}" for e in errors))
    # Lade Tabellen relativ zum JSON-Verzeichnis
    pkg.load_rate_tables(path.parent)
    return pkg
