"""
run_config.py -- Einheitliche Runtime-Konfiguration fuer C++ und Python Backend.

Primaere Quelle der Wahrheit fuer alle Laufparameter. Wird von der GUI erzeugt,
als run_config.json serialisiert und von beiden Backends direkt gelesen.

Schema:
  geometry.*    -- Kammergeometrie (R, L, betai, betag)
  grid.*        -- Extraktionsgitter (Vgrid, sgrid, eta_opt)
  coil.*        -- RF-Spule (frequency, Nw, R_ohm, Rc, lc)
  operation.*   -- Betriebspunkt (solve_mode, P_max, I_soll, ...)
  sweep.*       -- Sweep-Parameter (Q0_start, Q0_step, N)
  rates.*       -- Ratenmodell (rate_model)
  meta.*        -- Lauf-Metadaten (backend, package, preset, gas, ...)

Primaerer Pfad:
  GUI -> RunConfig -> run_config.json -> Backend (C++ oder Python)

Kompatibilitaet:
  to_params_txt()    -- Legacy-Export fuer alte C++-Workflows
  from_params_txt()  -- Legacy-Import
  save_json()        -- Primaere Serialisierung (run_config.json)
  load_json()        -- Primaeres Laden
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional


@dataclass
class GeometryConfig:
    R: float = 0.02           # Kammerradius [m]
    L: float = 0.04           # Kammerlaenge [m]
    betai: float = 0.5        # Ionentransparenz [-]
    betag: float = 0.05145    # Gastransparenz [-]


@dataclass
class GridConfig:
    Vgrid: float = 1500.0     # Gitterspannung [V]
    sgrid: float = 0.001      # Gitterabstand [m]
    eta_opt: float = 1.0      # Grid-Optik-Effizienz [-] (1.0=Legacy, <1 bei Preset)


@dataclass
class CoilConfig:
    frequency: float = 2.5e6  # RF-Frequenz [Hz]
    Nw: float = 6.0           # Windungszahl [-]
    R_ohm: float = 0.36       # Spulenwiderstand [Ohm]
    Rc: float = 0.02          # Spulenradius [m]
    lc: float = 0.04          # Spulenlaenge [m]


@dataclass
class OperationConfig:
    solve_mode: int = 1       # 1=I-fix, 2=SC
    P_max: float = 80.0       # Max. RF-Leistung [W]
    I_soll: float = 15.0      # Sollstrom [mA] (nur mode=1)
    density_profile_factor: float = 1.0  # Dichtekorrekturfaktor [-]
    alpha_e_wall: float = 7.0 # Wand-Energieverlustfaktor [-]
    rf_coupling: bool = False  # True: P_input = P_RFG, P_abs via RF model; False: P_input = P_abs


@dataclass
class SweepConfig:
    Q0_start: float = 0.27    # Start-Massenfluss [sccm]
    Q0_step: float = 0.01     # Schrittweite [sccm]
    N: int = 73               # Anzahl Schritte


@dataclass
class RateConfig:
    rate_model: int = 0       # 0=legacy, 1=conservative, 2=full tabulated


@dataclass
class MetaConfig:
    backend: str = ""         # "cpp" oder "python"
    package_id: str = ""      # z.B. "xenon_biagi", "py_iodine_lafleur_v1"
    gas: str = "xenon"        # Gas-Spezies
    cs_database: str = "biagi"  # Cross-Section-Datenbank
    preset_id: str = "custom" # Thruster-Preset-ID


@dataclass
class RunConfig:
    """Einheitliche Runtime-Konfiguration fuer beide Backends."""
    geometry: GeometryConfig = field(default_factory=GeometryConfig)
    grid: GridConfig = field(default_factory=GridConfig)
    coil: CoilConfig = field(default_factory=CoilConfig)
    operation: OperationConfig = field(default_factory=OperationConfig)
    sweep: SweepConfig = field(default_factory=SweepConfig)
    rates: RateConfig = field(default_factory=RateConfig)
    meta: MetaConfig = field(default_factory=MetaConfig)

    # ── Serialisierung ───────────────────────────────────────

    def to_json(self) -> dict:
        """Gibt die gesamte Konfiguration als strukturiertes Dict zurueck."""
        return asdict(self)

    def save_json(self, path: str | Path = "run_config.json"):
        """Schreibe primaere Konfigurationsdatei (JSON)."""
        with open(path, "w", encoding="utf-8") as f:
            json.dump(asdict(self), f, indent=2)

    @classmethod
    def load_json(cls, path: str | Path = "run_config.json") -> "RunConfig":
        """Lade RunConfig aus JSON-Datei (primaerer Pfad)."""
        with open(path, encoding="utf-8") as f:
            d = json.load(f)
        cfg = cls()
        if "geometry" in d:
            g = d["geometry"]
            cfg.geometry = GeometryConfig(**{k: g[k] for k in g if hasattr(cfg.geometry, k)})
        if "grid" in d:
            g = d["grid"]
            cfg.grid = GridConfig(**{k: g[k] for k in g if hasattr(cfg.grid, k)})
        if "coil" in d:
            c = d["coil"]
            cfg.coil = CoilConfig(**{k: c[k] for k in c if hasattr(cfg.coil, k)})
        if "operation" in d:
            o = d["operation"]
            cfg.operation = OperationConfig(**{k: o[k] for k in o if hasattr(cfg.operation, k)})
        if "sweep" in d:
            s = d["sweep"]
            cfg.sweep = SweepConfig(**{k: s[k] for k in s if hasattr(cfg.sweep, k)})
        if "rates" in d:
            r = d["rates"]
            cfg.rates = RateConfig(**{k: r[k] for k in r if hasattr(cfg.rates, k)})
        if "meta" in d:
            m = d["meta"]
            cfg.meta = MetaConfig(**{k: m[k] for k in m if hasattr(cfg.meta, k)})
        return cfg

    def to_flat_dict(self) -> dict:
        """Flache Schluessel-Wert-Darstellung (fuer Logging und Vergleich)."""
        flat = {}
        for section_name, section in asdict(self).items():
            if isinstance(section, dict):
                for k, v in section.items():
                    flat[k] = v
            else:
                flat[section_name] = section
        return flat

    def to_params_txt(self, path: str | Path = "params.txt"):
        """Schreibe params.txt (C++-kompatibel, flaches Key-Value-Format)."""
        # Mapping: RunConfig-Feld -> params.txt-Schluessel
        FIELD_MAP = {
            # Geometrie
            "R": self.geometry.R, "L": self.geometry.L,
            "betai": self.geometry.betai, "betag": self.geometry.betag,
            # Grid
            "Vgrid": self.grid.Vgrid, "sgrid": self.grid.sgrid,
            "eta_opt": self.grid.eta_opt,
            # Spule
            "frequency": self.coil.frequency, "Nw": self.coil.Nw,
            "R_ohm": self.coil.R_ohm, "Rc": self.coil.Rc, "lc": self.coil.lc,
            # Betrieb
            "solve_mode": self.operation.solve_mode,
            "P_RFG_max": self.operation.P_max,
            "I_soll": self.operation.I_soll,
            "density_profile_factor": self.operation.density_profile_factor,
            "alpha_e_wall": self.operation.alpha_e_wall,
            # Sweep
            "Q0sccm_start": self.sweep.Q0_start,
            "Q0sccm_step": self.sweep.Q0_step,
            "jjmax": self.sweep.N,
            # Raten
            "rate_model": self.rates.rate_model,
            # Meta (String-Werte)
            "gas_species": self.meta.gas,
            "cs_database": self.meta.cs_database,
        }

        with open(path, "w", encoding="utf-8") as f:
            f.write(f"# RunConfig generated — preset={self.meta.preset_id}\n")
            for key, val in FIELD_MAP.items():
                f.write(f"{key} {val}\n")
            f.write("use_paper_kel 1\n")

    @classmethod
    def from_params_txt(cls, path: str | Path = "params.txt") -> "RunConfig":
        """Lese params.txt und erzeuge RunConfig."""
        cfg = cls()
        p = Path(path)
        if not p.exists():
            return cfg

        params: dict[str, str] = {}
        for line in p.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            if len(parts) == 2:
                params[parts[0]] = parts[1]

        def _f(key, default=0.0):
            try:
                return float(params[key])
            except (KeyError, ValueError):
                return default

        def _i(key, default=0):
            try:
                return int(float(params[key]))
            except (KeyError, ValueError):
                return default

        def _s(key, default=""):
            return params.get(key, default)

        cfg.geometry = GeometryConfig(
            R=_f("R", 0.02), L=_f("L", 0.04),
            betai=_f("betai", 0.5), betag=_f("betag", 0.05145))
        cfg.grid = GridConfig(
            Vgrid=_f("Vgrid", 1500), sgrid=_f("sgrid", 0.001),
            eta_opt=_f("eta_opt", 1.0))
        cfg.coil = CoilConfig(
            frequency=_f("frequency", 2.5e6), Nw=_f("Nw", 6),
            R_ohm=_f("R_ohm", 0.36), Rc=_f("Rc", 0.02), lc=_f("lc", 0.04))
        cfg.operation = OperationConfig(
            solve_mode=_i("solve_mode", 1),
            P_max=_f("P_RFG_max", 80),
            I_soll=_f("I_soll", 15),
            density_profile_factor=_f("density_profile_factor", 1.0),
            alpha_e_wall=_f("alpha_e_wall", 7.0))
        cfg.sweep = SweepConfig(
            Q0_start=_f("Q0sccm_start", 0.27),
            Q0_step=_f("Q0sccm_step", 0.01),
            N=_i("jjmax", 73))
        cfg.rates = RateConfig(rate_model=_i("rate_model", 0))
        cfg.meta = MetaConfig(
            gas=_s("gas_species", "xenon"),
            cs_database=_s("cs_database", "biagi"))

        return cfg

    # ── Validierung ──────────────────────────────────────────

    def validate(self) -> list[str]:
        """Pruefe Konsistenz. Gibt Liste von Fehlern zurueck (leer = OK)."""
        errs = []
        g = self.geometry
        if g.R <= 0:   errs.append("R muss > 0 sein")
        if g.L <= 0:   errs.append("L muss > 0 sein")
        if g.betai <= 0 or g.betai > 1: errs.append("betai muss in (0, 1] liegen")
        if g.betag < 0 or g.betag > 1:  errs.append("betag muss in [0, 1] liegen")

        gr = self.grid
        if gr.Vgrid <= 0:   errs.append("Vgrid muss > 0 sein")
        if gr.sgrid <= 0:   errs.append("sgrid muss > 0 sein")
        if gr.eta_opt <= 0 or gr.eta_opt > 1: errs.append("eta_opt muss in (0, 1] liegen")

        c = self.coil
        if c.frequency <= 0: errs.append("frequency muss > 0 sein")
        if c.Rc < g.R:       errs.append(f"Rc ({c.Rc}) < R ({g.R}) — Spule innerhalb Kammer")

        o = self.operation
        if o.solve_mode not in (1, 2): errs.append(f"solve_mode={o.solve_mode} ungueltig (1 oder 2)")
        if o.P_max <= 0:     errs.append("P_max muss > 0 sein")
        if o.solve_mode == 1 and o.I_soll <= 0:
            errs.append("I-fix-Modus: I_soll muss > 0 sein")

        s = self.sweep
        if s.Q0_step <= 0:   errs.append("Q0_step muss > 0 sein")
        if s.N <= 0:         errs.append("N muss > 0 sein")

        return errs

    # ── Preset-Integration ───────────────────────────────────

    @classmethod
    def from_preset(cls, preset: dict, package_gas: str = "xenon",
                    package_id: str = "", cs_database: str = "biagi") -> "RunConfig":
        """Erzeuge RunConfig aus einem Thruster-Preset-Dict."""
        cfg = cls()
        p = preset.get("params", {})
        ext = preset.get("extraction", {})

        cfg.geometry = GeometryConfig(
            R=p.get("R", 0.02), L=p.get("L", 0.04),
            betai=p.get("betai", 0.5), betag=p.get("betag", 0.05145))
        cfg.grid = GridConfig(
            Vgrid=p.get("Vgrid", 1500), sgrid=p.get("sgrid", 0.001),
            eta_opt=ext.get("eta_opt", 1.0))
        cfg.coil = CoilConfig(
            frequency=p.get("frequency", 2.5e6), Nw=p.get("Nw", 6),
            R_ohm=p.get("R_ohm", 0.36), Rc=p.get("Rc", 0.02),
            lc=p.get("lc", 0.04))
        cfg.operation = OperationConfig(
            P_max=p.get("P_RFG_max", 80),
            I_soll=p.get("I_soll", 15),
            density_profile_factor=p.get("density_profile_factor", 1.0))
        cfg.sweep = SweepConfig(
            Q0_start=p.get("Q0sccm", 0.475))  # Preset-Default als Startwert
        cfg.meta = MetaConfig(
            preset_id=preset.get("id", "custom"),
            gas=package_gas, package_id=package_id,
            cs_database=cs_database)

        return cfg
