"""
gui.py – PyQt6-GUI fuer das Global Plasma Model mit Echtzeit-Streaming.

Konzept:
- Die GUI liest die strukturierte C++-Ausgabe direkt aus stdout/stderr
- RESULT-Zeilen werden sofort in eingebettete Live-Plots uebernommen

Voraussetzungen:
    pip install PyQt6 pyqtgraph
"""

from __future__ import annotations

import os
import sys
import shutil
import subprocess
from pathlib import Path

from PyQt6.QtCore import Qt, QProcess
from PyQt6.QtGui import QAction, QDesktopServices, QTextCursor
from PyQt6.QtWidgets import (
    QApplication,
    QButtonGroup,
    QComboBox,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QProgressBar,
    QRadioButton,
    QScrollArea,
    QSizePolicy,
    QSplitter,
    QStatusBar,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

import pyqtgraph as pg

SCRIPT_DIR = Path(__file__).resolve().parent
os.chdir(SCRIPT_DIR)

CPP_SOURCE = "chabert_modified.cpp"
CONFIG_FILE = "params.txt"
OUTPUT_FILE = "output_kh.txt"

# ─── Optionen ────────────────────────────────────────────────

SOLVE_MODE_OPTIONS = [
    (1, "Fester Strahlstrom"),
    (2, "Selbstkonsistent"),
]
DEFAULT_SOLVE_MODE = 1

RATE_MODEL_OPTIONS = [
    (0, "Legacy"),
    (1, "Conservative"),
    (2, "Full tabulated"),
]
DEFAULT_RATE_MODEL = 0

GAS_SPECIES_OPTIONS = [
    ("xenon",   "Xenon (Xe)"),
    ("krypton", "Krypton (Kr)"),
    ("argon",   "Argon (Ar)"),
]

# ─── Parameter-Definitionen ─────────────────────────────────
# Jeder Eintrag: (key, kurz-label, unit, default)

# Gruppe 1: System / Geometrie  (ohne Gas/Raten – die sitzen als Combos drin)
GEOM_PARAMS = [
    # Zeile 1: Kammer
    ("R",     "Radius R",   "m",   0.02),
    ("L",     "Laenge L",   "m",   0.04),
    # Zeile 2: Spule
    ("Nw",    "Windungen",  "--",  6.0),
    ("R_ohm", "R_ohm",      "Ohm", 0.36),
    # Zeile 3: Spule cont.
    ("frequency", "RF-Frequenz", "Hz", 2.5e6),
    ("lc",        "L_Spule",     "m",  0.04),
    # Zeile 4: Spule radius + Gitter
    ("Rc",    "R_Spule",  "m", 0.02),
    ("Vgrid", "V_grid",   "V", 1500.0),
    # Zeile 5: Gitter + Transparenzen
    ("sgrid", "s_grid",   "m",  0.001),
    ("betai", "betai",    "--", 0.5),
    # Zeile 6: Transparenzen
    ("betag", "betag",    "--", 0.05145),
]

# Gruppe 2: Betriebsbedingungen
OPERATE_PARAMS = [
    ("P_RFG",     "P_RFG",      "W",    18.0),
    ("P_RFG_max", "P_max",      "W",    80.0),
    ("Q0sccm",    "Massenfluss","sccm", 0.475),
    ("I_soll",    "I_soll",     "mA",   15.0),
    ("density_profile_factor", "Profil-Fakt.", "--", 1.0),
]

# Gruppe 4: Sweep
SWEEP_PARAMS = [
    ("Q0sccm_start", "Q0 Start",  "sccm", 0.27),
    ("Q0sccm_step",  "Schritt",   "sccm", 0.01),
    ("jjmax",        "Schritte",  "--",    73),
]

# Flache Lookup-Dicts
ALL_PRIMARY: dict[str, tuple[str, str, float]] = {}
for key, label, unit, default in GEOM_PARAMS + OPERATE_PARAMS:
    ALL_PRIMARY[key] = (label, unit, default)

ALL_SWEEP: dict[str, tuple[str, str, float]] = {}
for key, label, unit, default in SWEEP_PARAMS:
    ALL_SWEEP[key] = (label, unit, default)

# ─── Validierung ─────────────────────────────────────────────

VALIDATION_RULES: list[tuple[str, object, str]] = [
    ("R",         lambda v: v > 0,       "R muss > 0"),
    ("L",         lambda v: v > 0,       "L muss > 0"),
    ("betai",     lambda v: 0 < v <= 1,  "betai in (0,1]"),
    ("betag",     lambda v: 0 <= v <= 1, "betag in [0,1]"),
    ("frequency", lambda v: v > 0,       "Frequenz > 0"),
    ("Nw",        lambda v: v >= 1,      "Nw >= 1"),
    ("R_ohm",     lambda v: v > 0,       "R_ohm > 0"),
    ("Rc",        lambda v: v > 0,       "Rc > 0"),
    ("lc",        lambda v: v > 0,       "lc > 0"),
    ("Vgrid",     lambda v: v > 0,       "Vgrid > 0"),
    ("sgrid",     lambda v: v > 0,       "sgrid > 0"),
    ("P_RFG",     lambda v: v > 0,       "P_RFG > 0"),
    ("P_RFG_max", lambda v: v > 0,       "P_max > 0"),
    ("Q0sccm",    lambda v: v >= 0,      "Q0sccm >= 0"),
    ("I_soll",    lambda v: v > 0,       "I_soll > 0"),
    ("Q0sccm_start", lambda v: v >= 0,   "Q0 Start >= 0"),
    ("Q0sccm_step",  lambda v: v > 0,    "Schritt > 0"),
    ("jjmax",     lambda v: v >= 1,      "Schritte >= 1"),
    ("density_profile_factor", lambda v: 0 < v <= 1.0, "Profil-Fakt. in (0,1]"),
]


def cross_validate(vals: dict[str, float]) -> list[str]:
    errors = []
    if "Rc" in vals and "R" in vals and vals["Rc"] < vals["R"]:
        errors.append(f"Rc ({vals['Rc']:.4f}) < R ({vals['R']:.4f})")
    if "P_RFG" in vals and "P_RFG_max" in vals and vals["P_RFG"] > vals["P_RFG_max"]:
        errors.append(f"P_RFG ({vals['P_RFG']:.1f}) > P_max ({vals['P_RFG_max']:.1f})")
    return errors


# ─── Hilfsfunktion: kompaktes Eingabefeld ────────────────────

def _add_field(grid: QGridLayout, row: int, col: int,
               label: str, unit: str, default, entries: dict,
               key: str, width: int = 80):
    """Fuegt Label | Edit | Unit in 3 aufeinanderfolgende Spalten ein."""
    base = col * 3
    lbl = QLabel(label)
    lbl.setObjectName("ParamLabel")
    grid.addWidget(lbl, row, base)
    edit = QLineEdit(str(default))
    edit.setAlignment(Qt.AlignmentFlag.AlignRight)
    edit.setFixedWidth(width)
    grid.addWidget(edit, row, base + 1)
    u = QLabel(unit)
    u.setObjectName("MutedLabel")
    u.setFixedWidth(30)
    grid.addWidget(u, row, base + 2)
    entries[key] = edit


# ─── MetricCard ──────────────────────────────────────────────

class MetricCard(QFrame):
    def __init__(self, title: str, unit: str):
        super().__init__()
        self.setObjectName("MetricCard")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(4, 3, 4, 3)
        lay.setSpacing(1)
        t = QLabel(f"{title} [{unit}]")
        t.setObjectName("MetricTitle")
        self.value_lbl = QLabel("--")
        self.value_lbl.setObjectName("MetricValue")
        self.value_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        lay.addWidget(t)
        lay.addWidget(self.value_lbl)

    def set_value(self, v: str):
        self.value_lbl.setText(v)


# ─── LivePlotGrid ────────────────────────────────────────────

class LivePlotGrid(QGroupBox):
    def __init__(self):
        super().__init__("Live-Plots")
        layout = QGridLayout(self)
        layout.setHorizontalSpacing(4)
        layout.setVerticalSpacing(4)
        layout.setContentsMargins(2, 12, 2, 2)

        pg.setConfigOptions(antialias=True, background="#0b1220", foreground="#d7e3ff")

        self.data: dict = {
            "first_point": None,
            "Q0sccm": [], "P_RFG": [], "Te": [], "Tg": [],
            "n": [], "ng": [], "I_mA": [],
        }

        defs = [
            ("P_RFG", "P (W)"),  ("Te", "Te (eV)"),   ("Tg", "Tg (K)"),
            ("n", "n (m-3)"),    ("ng", "ng (m-3)"),   ("I_mA", "I (mA)"),
        ]

        self.widgets = {}
        self.curves = {}
        self.markers = {}
        self.first_markers = {}

        for idx, (key, title) in enumerate(defs):
            w = pg.PlotWidget()
            w.setMinimumHeight(100)
            w.showGrid(x=True, y=True, alpha=0.18)
            w.setTitle(title, size="7pt")
            w.setLabel("bottom", "Q0", units="sccm")
            w.getPlotItem().getAxis("left").enableAutoSIPrefix(True)
            w.getPlotItem().getAxis("bottom").enableAutoSIPrefix(False)

            pen = pg.mkPen(width=2)
            curve = w.plot([], [], pen=pen, symbol="o", symbolSize=4)
            marker = pg.ScatterPlotItem(size=7, brush=pg.mkBrush("w"), pen=pg.mkPen("w"))
            first_m = pg.ScatterPlotItem(size=10, brush=pg.mkBrush("#ff6b6b"), pen=pg.mkPen("#ff6b6b"))
            w.addItem(marker)
            w.addItem(first_m)

            self.widgets[key] = w
            self.curves[key] = curve
            self.markers[key] = marker
            self.first_markers[key] = first_m
            layout.addWidget(w, idx // 3, idx % 3)

    def clear(self):
        self.data["first_point"] = None
        for k in ("Q0sccm", "P_RFG", "Te", "Tg", "n", "ng", "I_mA"):
            self.data[k].clear()
        self.redraw()

    def add_point(self, q0, p_rfg, te, tg, n, ng, i_ma):
        if self.data["first_point"] is None:
            self.data["first_point"] = (q0, p_rfg, te, tg, n, ng, i_ma)
        else:
            self.data["Q0sccm"].append(q0)
            self.data["P_RFG"].append(p_rfg)
            self.data["Te"].append(te)
            self.data["Tg"].append(tg)
            self.data["n"].append(n)
            self.data["ng"].append(ng)
            self.data["I_mA"].append(i_ma)
        self.redraw()

    def redraw(self):
        x = self.data["Q0sccm"]
        fp = self.data["first_point"]
        fp_map = None
        if fp:
            q, p, te, tg, n, ng, i = fp
            fp_map = {"P_RFG": (q, p), "Te": (q, te), "Tg": (q, tg),
                      "n": (q, n), "ng": (q, ng), "I_mA": (q, i)}

        for key in self.curves:
            y = self.data[key]
            self.curves[key].setData(x, y)
            self.markers[key].setData([{"pos": (x[-1], y[-1])}] if x else [])
            if fp_map:
                qv, yv = fp_map[key]
                self.first_markers[key].setData([{"pos": (qv, yv)}])
            else:
                self.first_markers[key].setData([])


# ═════════════════════════════════════════════════════════════
# Hauptfenster
# ═════════════════════════════════════════════════════════════

class SimulatorWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Global Plasma Model")
        screen = QApplication.primaryScreen()
        if screen:
            avail = screen.availableGeometry()
            w = min(1580, int(avail.width() * 0.93))
            h = min(980, int(avail.height() * 0.91))
            self.resize(w, h)
            self.move(avail.x() + (avail.width() - w) // 2,
                      avail.y() + (avail.height() - h) // 2)
        else:
            self.resize(1580, 980)

        self.entries: dict[str, QLineEdit] = {}
        self.sweep_entries: dict[str, QLineEdit] = {}
        self.metric_cards: dict[str, MetricCard] = {}
        self.process: QProcess | None = None
        self.cancel_requested = False
        self._stdout_buffer = ""
        self.current_q0 = None

        self._build_ui()
        self._apply_styles()

    # ── Haupt-UI ─────────────────────────────────────────────

    def _build_ui(self):
        root = QWidget()
        self.setCentralWidget(root)
        outer = QVBoxLayout(root)
        outer.setContentsMargins(6, 4, 6, 4)
        outer.setSpacing(4)

        # ── Header ───────────────────────────────────────────
        header = QFrame()
        header.setObjectName("Header")
        hl = QHBoxLayout(header)
        hl.setContentsMargins(12, 6, 12, 6)
        title = QLabel("Global Plasma Model")
        title.setObjectName("HeaderTitle")
        subtitle = QLabel("0D Solver  --  Echtzeit-Streaming")
        subtitle.setObjectName("HeaderSubtitle")
        hl.addWidget(title)
        hl.addStretch(1)

        # Header-rechts: kompakte Status-Info
        self._hdr_gas = QLabel("Xe")
        self._hdr_gas.setObjectName("HeaderTag")
        self._hdr_rate = QLabel("Legacy")
        self._hdr_rate.setObjectName("HeaderTag")
        self._hdr_mode = QLabel("I-fix")
        self._hdr_mode.setObjectName("HeaderTag")
        for w in (self._hdr_gas, self._hdr_rate, self._hdr_mode):
            hl.addWidget(w)

        hl.addWidget(subtitle)
        outer.addWidget(header)

        # ── 2-Spalten-Splitter ───────────────────────────────
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setChildrenCollapsible(False)
        outer.addWidget(splitter, 1)

        # LINKS: scrollbar
        left_scroll = QScrollArea()
        left_scroll.setWidgetResizable(True)
        left_scroll.setFrameShape(QFrame.Shape.NoFrame)
        left_panel = QWidget()
        left_scroll.setWidget(left_panel)
        left_lay = QVBoxLayout(left_panel)
        left_lay.setContentsMargins(0, 0, 4, 0)
        left_lay.setSpacing(4)

        left_lay.addWidget(self._build_group_geometry())
        left_lay.addWidget(self._build_group_operate())
        left_lay.addWidget(self._build_group_model())
        left_lay.addWidget(self._build_group_sweep())
        left_lay.addWidget(self._build_group_actions())
        left_lay.addStretch(1)

        # RECHTS
        right_w = QWidget()
        right_lay = QVBoxLayout(right_w)
        right_lay.setContentsMargins(4, 0, 0, 0)
        right_lay.setSpacing(3)

        right_scroll = QScrollArea()
        right_scroll.setWidgetResizable(True)
        right_scroll.setFrameShape(QFrame.Shape.NoFrame)
        right_scroll.setWidget(right_w)

        splitter.addWidget(left_scroll)
        splitter.addWidget(right_scroll)
        splitter.setStretchFactor(0, 33)
        splitter.setStretchFactor(1, 67)

        # Fortschritt
        right_lay.addWidget(self._build_progress_row())

        # Metriken 3x3
        right_lay.addWidget(self._build_metrics_grid())

        # Plots + Log (vertikal gesplittet)
        self.plot_grid = LivePlotGrid()
        log_box = self._build_log_group()
        log_box.setMinimumHeight(80)

        vsplit = QSplitter(Qt.Orientation.Vertical)
        vsplit.setChildrenCollapsible(False)
        vsplit.addWidget(self.plot_grid)
        vsplit.addWidget(log_box)
        vsplit.setStretchFactor(0, 4)
        vsplit.setStretchFactor(1, 1)
        right_lay.addWidget(vsplit, 1)

        self.setStatusBar(QStatusBar())
        self.set_status("Bereit.")

        # Keyboard shortcuts
        for text, handler in [
            ("Defaults", self.reset_defaults),
            ("Kompilieren", self.compile_cpp),
            ("Simulation", self.run_sim),
            ("Abbrechen", self.cancel_sim),
        ]:
            a = QAction(text, self)
            a.triggered.connect(handler)
            self.addAction(a)

    # ── Gruppe 1: System / Geometrie ─────────────────────────

    def _build_group_geometry(self) -> QGroupBox:
        box = QGroupBox("System / Geometrie")
        g = QGridLayout(box)
        g.setContentsMargins(6, 14, 6, 4)
        g.setHorizontalSpacing(8)
        g.setVerticalSpacing(3)

        # Zeile 0: Gas + Ratenmodell (Combos)
        g.addWidget(QLabel("Gas"), 0, 0)
        self.cmb_gas_species = QComboBox()
        for val, txt in GAS_SPECIES_OPTIONS:
            self.cmb_gas_species.addItem(txt, val)
        self.cmb_gas_species.setCurrentIndex(0)
        self.cmb_gas_species.setToolTip("Gasspezies (M, Eiz, Cross-Sections)")
        self.cmb_gas_species.currentIndexChanged.connect(self._sync_header_tags)
        g.addWidget(self.cmb_gas_species, 0, 1, 1, 2)

        g.addWidget(QLabel("Ratenmodell"), 0, 3)
        self.cmb_rate_model = QComboBox()
        for val, txt in RATE_MODEL_OPTIONS:
            self.cmb_rate_model.addItem(txt, val)
        self.cmb_rate_model.setCurrentIndex(DEFAULT_RATE_MODEL)
        self.cmb_rate_model.setToolTip(
            "Legacy: Chabert 2012 Fits, Kel=const\n"
            "Conservative: Kiz+Kex tabelliert, Kel=const\n"
            "Full: Alle Raten aus Biagi/LXCat")
        self.cmb_rate_model.currentIndexChanged.connect(self._sync_header_tags)
        g.addWidget(self.cmb_rate_model, 0, 4, 1, 2)

        # Zeile 1-6: Geometrie-Parameter paarweise
        pairs = [
            (GEOM_PARAMS[0], GEOM_PARAMS[1]),   # R, L
            (GEOM_PARAMS[2], GEOM_PARAMS[3]),   # Nw, R_ohm
            (GEOM_PARAMS[4], GEOM_PARAMS[5]),   # freq, lc
            (GEOM_PARAMS[6], GEOM_PARAMS[7]),   # Rc, Vgrid
            (GEOM_PARAMS[8], GEOM_PARAMS[9]),   # sgrid, betai
            (GEOM_PARAMS[10],),                  # betag
        ]
        for row_idx, pair in enumerate(pairs, start=1):
            for col_idx, (key, label, unit, default) in enumerate(pair):
                _add_field(g, row_idx, col_idx, label, unit, default, self.entries, key)

        return box

    # ── Gruppe 2: Betriebsbedingungen ────────────────────────

    def _build_group_operate(self) -> QGroupBox:
        box = QGroupBox("Betriebsbedingungen")
        g = QGridLayout(box)
        g.setContentsMargins(6, 14, 6, 4)
        g.setHorizontalSpacing(8)
        g.setVerticalSpacing(3)

        # Zeile 0: Betriebsmodus (Radio)
        self.mode_group = QButtonGroup(self)
        self.mode_buttons: dict[int, QRadioButton] = {}
        g.addWidget(QLabel("Modus"), 0, 0)
        for i, (val, txt) in enumerate(SOLVE_MODE_OPTIONS):
            btn = QRadioButton(txt)
            if val == DEFAULT_SOLVE_MODE:
                btn.setChecked(True)
            self.mode_group.addButton(btn, val)
            self.mode_buttons[val] = btn
            g.addWidget(btn, 0, 1 + i, 1, 2)
        self.mode_group.idToggled.connect(self._on_mode_changed)

        # Zeile 1: P_RFG, P_max
        _add_field(g, 1, 0, "P_RFG", "W", 18.0, self.entries, "P_RFG")
        _add_field(g, 1, 1, "P_max", "W", 80.0, self.entries, "P_RFG_max")

        # Zeile 2: Q0, I_soll
        _add_field(g, 2, 0, "Massenfluss", "sccm", 0.475, self.entries, "Q0sccm")
        _add_field(g, 2, 1, "I_soll", "mA", 15.0, self.entries, "I_soll")

        # Zeile 3: density_profile_factor
        _add_field(g, 3, 0, "Profil-Fakt.", "--", 1.0, self.entries, "density_profile_factor")

        # I_soll disable bei Selbstkonsistent
        self._on_mode_changed(self.mode_group.checkedId(), True)

        return box

    def _on_mode_changed(self, mode_id: int, checked: bool = True):
        if not checked:
            return
        is_sc = (mode_id == 2)
        if "I_soll" in self.entries:
            self.entries["I_soll"].setEnabled(not is_sc)
        self._sync_header_tags()

    # ── Gruppe 3: Modelloptionen ─────────────────────────────

    def _build_group_model(self) -> QGroupBox:
        box = QGroupBox("Modelloptionen")
        g = QGridLayout(box)
        g.setContentsMargins(6, 14, 6, 4)
        g.setHorizontalSpacing(8)
        g.setVerticalSpacing(3)

        # Zeile 0: Info-Labels (read-only Zusammenfassung)
        g.addWidget(QLabel("Ratenmodell:"), 0, 0)
        self._lbl_rate_info = QLabel("--")
        self._lbl_rate_info.setObjectName("MutedLabel")
        g.addWidget(self._lbl_rate_info, 0, 1, 1, 3)

        g.addWidget(QLabel("Gas:"), 1, 0)
        self._lbl_gas_info = QLabel("--")
        self._lbl_gas_info.setObjectName("MutedLabel")
        g.addWidget(self._lbl_gas_info, 1, 1, 1, 3)

        # Wird aktualisiert ueber _sync_header_tags
        self._sync_header_tags()

        return box

    # ── Gruppe 4: Sweep / Studie ─────────────────────────────

    def _build_group_sweep(self) -> QGroupBox:
        box = QGroupBox("Sweep / Studie")
        g = QGridLayout(box)
        g.setContentsMargins(6, 14, 6, 4)
        g.setHorizontalSpacing(8)
        g.setVerticalSpacing(3)

        # Zeile 0: Q0 Start, Schritt, Schritte
        for i, (key, label, unit, default) in enumerate(SWEEP_PARAMS):
            _add_field(g, 0, i, label, unit, default, self.sweep_entries, key, width=65)

        # Zeile 1: Q0 Ende (berechnet, read-only)
        g.addWidget(QLabel("Q0 Ende"), 1, 0)
        self._lbl_q0_end = QLineEdit("--")
        self._lbl_q0_end.setReadOnly(True)
        self._lbl_q0_end.setFixedWidth(65)
        self._lbl_q0_end.setAlignment(Qt.AlignmentFlag.AlignRight)
        self._lbl_q0_end.setObjectName("ReadOnlyField")
        g.addWidget(self._lbl_q0_end, 1, 1)
        g.addWidget(QLabel("sccm"), 1, 2)

        # Verbinde Felder fuer automatische Berechnung
        for key in ("Q0sccm_start", "Q0sccm_step", "jjmax"):
            if key in self.sweep_entries:
                self.sweep_entries[key].textChanged.connect(self._update_q0_end)
        self._update_q0_end()

        return box

    def _update_q0_end(self):
        try:
            start = float(self.sweep_entries["Q0sccm_start"].text())
            step = float(self.sweep_entries["Q0sccm_step"].text())
            n = int(float(self.sweep_entries["jjmax"].text()))
            end = start + step * (n - 1)
            self._lbl_q0_end.setText(f"{end:.4f}")
        except (ValueError, KeyError):
            self._lbl_q0_end.setText("--")

    # ── Gruppe 5: Aktionen ───────────────────────────────────

    def _build_group_actions(self) -> QGroupBox:
        box = QGroupBox("Aktionen")
        lay = QVBoxLayout(box)
        lay.setContentsMargins(6, 14, 6, 4)
        lay.setSpacing(3)

        # Hauptreihe
        r1 = QHBoxLayout()
        r1.setSpacing(4)
        r1.addWidget(self._btn("Kompilieren", self.compile_cpp))
        r1.addWidget(self._btn("Simulation starten", self.run_sim, primary=True))
        r1.addWidget(self._btn("Abbrechen", self.cancel_sim, danger=True))
        lay.addLayout(r1)

        # Sekundaer
        r2 = QHBoxLayout()
        r2.setSpacing(4)
        r2.addWidget(self._btn("Defaults", self.reset_defaults, small=True))
        r2.addWidget(self._btn("Plots leeren", self.clear_plots, small=True))
        r2.addWidget(self._btn("Output", self.open_output, small=True))
        r2.addWidget(self._btn("Log-Viewer", self.open_log_viewer, small=True))
        lay.addLayout(r2)

        return box

    def _btn(self, text, slot, primary=False, danger=False, small=False):
        b = QPushButton(text)
        b.clicked.connect(slot)
        if primary:
            b.setProperty("primary", True)
        if danger:
            b.setProperty("danger", True)
        b.setMinimumHeight(24 if small else 32)
        return b

    # ── Rechte Seite: Fortschritt ────────────────────────────

    def _build_progress_row(self) -> QWidget:
        w = QWidget()
        lay = QHBoxLayout(w)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(6)
        lay.addWidget(QLabel("Scan:"))
        self.scan_progress = QProgressBar()
        self.scan_progress.setRange(0, 100)
        self.scan_progress.setFixedHeight(16)
        lay.addWidget(self.scan_progress, 3)
        self.scan_label = QLabel("--")
        self.scan_label.setFixedWidth(55)
        lay.addWidget(self.scan_label)
        lay.addWidget(QLabel("Iter:"))
        self.pid_progress = QProgressBar()
        self.pid_progress.setRange(0, 100)
        self.pid_progress.setFixedHeight(16)
        lay.addWidget(self.pid_progress, 2)
        self.pid_label = QLabel("--")
        self.pid_label.setFixedWidth(55)
        lay.addWidget(self.pid_label)
        return w

    # ── Rechte Seite: Metriken 3x3 Grid ─────────────────────

    def _build_metrics_grid(self) -> QWidget:
        w = QWidget()
        g = QGridLayout(w)
        g.setContentsMargins(0, 0, 0, 0)
        g.setSpacing(3)

        defs = [
            # Reihe 1
            ("Q0",  "Q0sccm", "sccm"), ("P_sol", "P_RFG", "W"),    ("P_set", "P_set", "W"),
            # Reihe 2
            ("I",   "I_beam", "mA"),    ("Te",    "Te",    "eV"),   ("Tg",    "Tg",   "K"),
            # Reihe 3
            ("n",   "n",      "m-3"),   ("ng",    "ng",    "m-3"),  ("iondeg","iondeg","%"),
        ]

        for i, (title, key, unit) in enumerate(defs):
            card = MetricCard(title, unit)
            g.addWidget(card, i // 3, i % 3)
            self.metric_cards[key] = card

        return w

    # ── Log ──────────────────────────────────────────────────

    def _build_log_group(self) -> QGroupBox:
        box = QGroupBox("Log")
        lay = QVBoxLayout(box)
        lay.setContentsMargins(2, 12, 2, 2)
        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self.log.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        lay.addWidget(self.log)
        return box

    # ── Styles ───────────────────────────────────────────────

    def _apply_styles(self):
        p = SCRIPT_DIR / "style.qss"
        try:
            self.setStyleSheet(p.read_text(encoding="utf-8"))
        except FileNotFoundError:
            pass

    # ── Header-Tags synchronisieren ──────────────────────────

    def _sync_header_tags(self, *_args):
        gas_text = self.cmb_gas_species.currentText() if hasattr(self, 'cmb_gas_species') else "?"
        rate_text = self.cmb_rate_model.currentText() if hasattr(self, 'cmb_rate_model') else "?"
        mode_text = "SC" if self.solve_mode() == 2 else "I-fix"
        if hasattr(self, '_hdr_gas'):
            self._hdr_gas.setText(gas_text)
            self._hdr_rate.setText(rate_text)
            self._hdr_mode.setText(mode_text)
        if hasattr(self, '_lbl_rate_info'):
            self._lbl_rate_info.setText(rate_text)
        if hasattr(self, '_lbl_gas_info'):
            self._lbl_gas_info.setText(gas_text)

    # ═════════════════════════════════════════════════════════
    # Backend-Logik (unveraendert)
    # ═════════════════════════════════════════════════════════

    def solve_mode(self) -> int:
        return self.mode_group.checkedId() if hasattr(self, 'mode_group') else 1

    def log_msg(self, msg: str):
        self.log.append(msg)
        c = self.log.textCursor()
        c.movePosition(QTextCursor.MoveOperation.End)
        self.log.setTextCursor(c)

    def set_status(self, msg: str):
        self.statusBar().showMessage(msg)

    def update_metric(self, key: str, value: str):
        c = self.metric_cards.get(key)
        if c:
            c.set_value(value)

    def update_solution_power(self, v):
        self.update_metric("P_RFG", "--" if v is None else f"{v:.2f}")

    def update_set_power(self, v):
        self.update_metric("P_set", "--" if v is None else f"{v:.2f}")

    def solver_method(self) -> int:
        return 4

    def solver_label(self) -> str:
        return "Newton (stationaer)"

    def wsl_available(self) -> bool:
        return sys.platform == "win32" and shutil.which("wsl") is not None

    def win_to_wsl(self, path: str) -> str:
        path = path.replace("\\", "/")
        if len(path) >= 2 and path[1] == ":":
            return f"/mnt/{path[0].lower()}{path[2:]}"
        return path

    def read_float(self, le: QLineEdit, name: str) -> float:
        try:
            return float(le.text().strip())
        except ValueError as e:
            raise ValueError(f"{name}: ungueltiger Wert") from e

    # ── Validierung + Config ─────────────────────────────────

    def validate_inputs(self) -> bool:
        errors: list[str] = []
        vals: dict[str, float] = {}
        all_e = {**self.entries, **self.sweep_entries}
        for key, fn, msg in VALIDATION_RULES:
            if key not in all_e:
                continue
            try:
                v = float(all_e[key].text().strip())
                vals[key] = v
                if not fn(v):
                    errors.append(f"{msg} (={v})")
            except ValueError:
                errors.append(f"{key}: ungueltig")
        errors.extend(cross_validate(vals))
        if errors:
            QMessageBox.warning(self, "Eingabefehler",
                                "\n".join(f"- {e}" for e in errors))
            return False
        return True

    def write_config(self) -> bool:
        if not self.validate_inputs():
            return False
        try:
            with open(CONFIG_FILE, "w", encoding="utf-8") as f:
                f.write("# Automatisch generiert\n")
                for key, edit in self.entries.items():
                    f.write(f"{key} {self.read_float(edit, key)}\n")
                for key, edit in self.sweep_entries.items():
                    f.write(f"{key} {self.read_float(edit, key)}\n")
                f.write(f"method {self.solver_method()}\n")
                f.write(f"solve_mode {self.solve_mode()}\n")
                f.write(f"gas_species {self.cmb_gas_species.currentData()}\n")
                f.write(f"rate_model {self.cmb_rate_model.currentData()}\n")
                f.write("use_paper_kel 1\n")
            g = self.cmb_gas_species.currentText()
            r = self.cmb_rate_model.currentText()
            m = "SC" if self.solve_mode() == 2 else "I-fix"
            self.log_msg(f"Config: {g}, {m}, {r}")
            return True
        except ValueError as e:
            QMessageBox.critical(self, "Fehler", str(e))
            return False
        except OSError as e:
            QMessageBox.critical(self, "Fehler", str(e))
            return False

    # ── Kompilierung ─────────────────────────────────────────

    def compile_cpp(self):
        if hasattr(self, '_compile_proc') and self._compile_proc and \
           self._compile_proc.state() != QProcess.ProcessState.NotRunning:
            QMessageBox.information(self, "Info", "Kompilierung laeuft.")
            return
        self.set_status("Kompiliere ...")
        self.log_msg("Kompilierung gestartet ...")
        if self.wsl_available():
            cwd = self.win_to_wsl(str(SCRIPT_DIR))
            prog = "wsl"
            args = ["bash", "-c",
                    f'cd "{cwd}" && '
                    f'g++ -O3 -march=native -std=c++17 -c bessel_wrapper.cpp -o bessel_wrapper.o 2>&1 && '
                    f'g++ -O3 -march=native -std=c++17 {CPP_SOURCE} bessel_wrapper.o -o chabert 2>&1']
        else:
            cc = shutil.which("g++")
            if not cc:
                self.log_msg("FEHLER: g++ nicht gefunden")
                self.set_status("Kein Compiler.")
                QMessageBox.critical(self, "Fehler", "g++ nicht gefunden.")
                return
            out = "chabert.exe" if sys.platform == "win32" else "chabert"
            obj = "bessel_wrapper.obj" if sys.platform == "win32" else "bessel_wrapper.o"
            if sys.platform == "win32":
                prog = "cmd"
                args = ["/c",
                        f'"{cc}" -O3 -march=native -std=c++17 -c bessel_wrapper.cpp -o {obj} && '
                        f'"{cc}" -O3 -march=native -std=c++17 {CPP_SOURCE} {obj} -o {out}']
            else:
                prog = "bash"
                args = ["-c",
                        f'"{cc}" -O3 -march=native -std=c++17 -c bessel_wrapper.cpp -o {obj} && '
                        f'"{cc}" -O3 -march=native -std=c++17 {CPP_SOURCE} {obj} -o {out}']
        self._compile_proc = QProcess(self)
        self._compile_proc.setWorkingDirectory(str(SCRIPT_DIR))
        self._compile_proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self._compile_proc.readyReadStandardOutput.connect(self._on_compile_out)
        self._compile_proc.finished.connect(self._on_compile_done)
        self._compile_proc.start(prog, args)

    def _on_compile_out(self):
        if self._compile_proc:
            t = bytes(self._compile_proc.readAllStandardOutput()).decode("utf-8", errors="replace")
            for l in t.splitlines():
                if l.strip():
                    self.log_msg(f"  [gcc] {l.strip()}")

    def _on_compile_done(self, code, _):
        if code == 0:
            self.log_msg("Kompilierung OK.")
            self.set_status("Kompiliert.")
        else:
            self.log_msg(f"Kompilierung fehlgeschlagen ({code}).")
            self.set_status("Kompilierfehler.")
            QMessageBox.critical(self, "Fehler", f"Kompilierung fehlgeschlagen ({code}).")

    # ── Simulation ───────────────────────────────────────────

    def run_sim(self):
        if self.process and self.process.state() != QProcess.ProcessState.NotRunning:
            QMessageBox.information(self, "Info", "Simulation laeuft.")
            return
        bn = "chabert.exe" if sys.platform == "win32" else "chabert"
        if self.wsl_available():
            cwd = self.win_to_wsl(str(SCRIPT_DIR))
            chk = subprocess.run(["wsl", "bash", "-c", f'test -f "{cwd}/chabert"'])
            if chk.returncode != 0:
                QMessageBox.warning(self, "Fehler", "Binary nicht gefunden.")
                return
        else:
            if not (SCRIPT_DIR / bn).exists():
                QMessageBox.warning(self, "Fehler", f"'{bn}' nicht gefunden.")
                return
        if not self.write_config():
            return
        self.cancel_requested = False
        self._stdout_buffer = ""
        self.current_q0 = None
        self.scan_progress.setValue(0)
        self.pid_progress.setValue(0)
        self.scan_label.setText("--")
        self.pid_label.setText("--")
        self.clear_plots()
        self.process = QProcess(self)
        self.process.setWorkingDirectory(str(SCRIPT_DIR))
        self.process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.process.readyReadStandardOutput.connect(self._on_sim_out)
        self.process.finished.connect(self._on_sim_done)
        g = self.cmb_gas_species.currentText()
        r = self.cmb_rate_model.currentText()
        self.log_msg(f"\nSimulation: {g}, {r}")
        self.set_status("Simulation laeuft ...")
        if self.wsl_available():
            cwd = self.win_to_wsl(str(SCRIPT_DIR))
            self.process.start("wsl", ["bash", "-c", f'cd "{cwd}" && ./chabert {CONFIG_FILE} 2>&1'])
        else:
            self.process.start(str(SCRIPT_DIR / bn), [CONFIG_FILE])
        if not self.process.waitForStarted(3000):
            QMessageBox.critical(self, "Fehler", "Start fehlgeschlagen.")
            self.set_status("Startfehler.")

    def _on_sim_out(self):
        if not self.process:
            return
        self._stdout_buffer += bytes(self.process.readAllStandardOutput()).decode("utf-8", errors="replace")
        while "\n" in self._stdout_buffer:
            line, self._stdout_buffer = self._stdout_buffer.split("\n", 1)
            self._parse_line(line.strip())

    def _parse_line(self, line: str):
        if not line:
            return
        p = line.split()
        tag = p[0]
        try:
            if tag == "Q0_STEP" and len(p) >= 4:
                self.current_q0 = float(p[1])
                pct = int(p[2]) / int(p[3]) * 100
                self.scan_progress.setValue(int(pct))
                self.scan_label.setText(f"{p[2]}/{p[3]}")
                self.update_metric("Q0sccm", f"{float(p[1]):.4f}")
                self.pid_progress.setValue(0)
                self.pid_label.setText("--")
                self.log_msg(f"Q0={p[1]} ({p[2]}/{p[3]})")
                return
            if tag == "PID_START" and len(p) >= 3:
                self.pid_progress.setValue(int(min(int(p[1])/50*100, 100)))
                self.pid_label.setText(f"It {p[1]}")
                self.update_set_power(float(p[2]))
                return
            if tag == "PID_DONE" and len(p) >= 6:
                self.update_metric("I_beam", f"{float(p[1]):.3f}")
                self.update_metric("Te", f"{float(p[4]):.2f}")
                self.update_metric("Tg", f"{float(p[5]):.1f}")
                self.update_set_power(float(p[3]))
                return
            if tag == "CONVERGED":
                self.pid_progress.setValue(100)
                self.pid_label.setText("OK")
                return
            if tag == "RESULT" and len(p) >= 7:
                n, ng, te, tg, i, pf = [float(x) for x in p[1:7]]
                self.update_metric("n", f"{n:.2e}")
                self.update_metric("ng", f"{ng:.2e}")
                self.update_metric("Te", f"{te:.2f}")
                self.update_metric("Tg", f"{tg:.1f}")
                self.update_metric("I_beam", f"{i:.3f}")
                if ng > 0:
                    self.update_metric("iondeg", f"{n/ng*100:.2f}")
                self.update_solution_power(pf)
                self.update_set_power(pf)
                if self.current_q0 is not None:
                    self.plot_grid.add_point(self.current_q0, pf, te, tg, n, ng, i)
                return
            if tag == "PID_MAXITER":
                self.log_msg("  PID: max iter")
                return
            if tag == "P_LIMIT_REACHED" and len(p) >= 6:
                self.log_msg(f"  RF-Limit Q0={float(p[1]):.4f} P={float(p[2]):.1f}>{float(p[3]):.1f}W")
                self.update_set_power(float(p[3]))
                return
            if tag in ("SWEEP_RECOVERY", "SWEEP_RECOVERY_FAIL"):
                self.log_msg(f"  {tag}: {' '.join(p[1:])}")
                return
            if tag == "SOLVER_FAIL":
                self.update_solution_power(None)
                self.log_msg(f"  FAIL: {line}")
                return
            if not line.startswith("STEP "):
                self.log_msg(line)
        except Exception as e:
            self.log_msg(f"Parse: {line} | {e}")

    def _on_sim_done(self, code, _):
        if self.cancel_requested:
            self.log_msg("\nAbgebrochen.")
            self.set_status("Abgebrochen.")
            self.cancel_requested = False
            return
        if code == 0:
            self.log_msg("\nFertig.")
            self.set_status("Fertig.")
            self.scan_progress.setValue(100)
            if Path(OUTPUT_FILE).exists():
                try:
                    n = max(0, len(Path(OUTPUT_FILE).read_text(encoding="utf-8").splitlines()) - 1)
                    self.log_msg(f"{OUTPUT_FILE}: {n} Punkte")
                except Exception:
                    pass
        else:
            self.log_msg(f"\nFehler ({code}).")
            self.set_status(f"Fehler ({code}).")

    # ── Aktionen ─────────────────────────────────────────────

    def cancel_sim(self):
        if self.process and self.process.state() != QProcess.ProcessState.NotRunning:
            self.cancel_requested = True
            self.process.kill()
            self.set_status("Abbruch ...")
        else:
            self.log_msg("Kein Prozess.")

    def reset_defaults(self):
        for key, (_, _, d) in ALL_PRIMARY.items():
            if key in self.entries:
                self.entries[key].setText(str(d))
        for key, (_, _, d) in ALL_SWEEP.items():
            if key in self.sweep_entries:
                self.sweep_entries[key].setText(str(d))
        self.mode_buttons[DEFAULT_SOLVE_MODE].setChecked(True)
        self.cmb_gas_species.setCurrentIndex(0)
        self.cmb_rate_model.setCurrentIndex(DEFAULT_RATE_MODEL)
        self.log_msg("Defaults geladen.")

    def clear_plots(self):
        self.plot_grid.clear()

    def open_output(self):
        p = SCRIPT_DIR / OUTPUT_FILE
        if not p.exists():
            QMessageBox.information(self, "Info", f"'{OUTPUT_FILE}' nicht vorhanden.")
            return
        QDesktopServices.openUrl(p.as_uri())

    def open_log_viewer(self):
        vs = SCRIPT_DIR / "log_viewer.py"
        if not vs.exists():
            QMessageBox.warning(self, "Fehler", "log_viewer.py nicht gefunden.")
            return
        logs = sorted(SCRIPT_DIR.glob("simulation_log_*.txt"),
                      key=lambda x: x.stat().st_mtime, reverse=True)
        args = [sys.executable, str(vs)]
        if logs:
            args.append(str(logs[0]))
        subprocess.Popen(args, cwd=str(SCRIPT_DIR))
        self.log_msg("Log-Viewer gestartet" + (f" ({logs[0].name})" if logs else ""))

    def closeEvent(self, event):
        try:
            if self.process and self.process.state() != QProcess.ProcessState.NotRunning:
                self.cancel_requested = True
                self.process.kill()
                self.process.waitForFinished(1500)
        except Exception:
            pass
        super().closeEvent(event)


# ═════════════════════════════════════════════════════════════

def main() -> int:
    app = QApplication(sys.argv)
    w = SimulatorWindow()
    w.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
