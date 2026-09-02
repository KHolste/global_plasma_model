"""
gui.py -- PyQt6-GUI fuer das Global Plasma Model mit Echtzeit-Streaming.

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
    QCheckBox,
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

CPP_SOURCE = "main.cpp"
CONFIG_FILE = "params.txt"
OUTPUT_FILE = "output_kh.txt"

# ─── UI-Skalierung ──────────────────────────────────────────
# Erkennt DPI / Skalierungsfaktor und stellt Hilfsfunktionen bereit.
# Bei 96 DPI (100%) ist der Faktor 1.0; bei 120 DPI (125%) ist er 1.25 usw.
# Manueller Override: UI_SCALE_OVERRIDE auf gewuenschten Wert setzen (z.B. 1.0).

UI_SCALE_OVERRIDE: float | None = None  # None = automatisch, float = fester Faktor


class UI:
    """Zentrale UI-Skalierung. Einmal init() aufrufen nach QApplication-Erstellung."""
    _factor: float = 1.0

    @classmethod
    def init(cls):
        if UI_SCALE_OVERRIDE is not None:
            cls._factor = max(0.5, min(3.0, UI_SCALE_OVERRIDE))
            return
        screen = QApplication.primaryScreen()
        if screen:
            dpi = screen.logicalDotsPerInch()
            cls._factor = max(0.75, min(2.5, dpi / 96.0))
        else:
            cls._factor = 1.0

    @classmethod
    def factor(cls) -> float:
        return cls._factor

    @classmethod
    def px(cls, base: int) -> int:
        """Skaliere einen Pixel-Wert (fuer setFixedWidth, setMinimumHeight etc.)."""
        return max(1, round(base * cls._factor))

    @classmethod
    def pt(cls, base: float) -> float:
        """Skaliere einen pt-Wert fuer Stylesheets. Qt skaliert pt bereits mit DPI,
        daher nur bei manuellem Override oder extremen DPI-Werten korrigieren."""
        return round(base * cls._factor, 1)

    @classmethod
    def qss(cls) -> str:
        """Generiert das Stylesheet mit skalierten Schriftgroessen."""
        f = cls.pt  # Kurzform
        return f"""
/* Auto-scaled stylesheet (factor={cls._factor:.2f}) */

QMainWindow, QWidget {{
    background: #0b1020;
    color: #e6edf7;
    font-size: {f(9)}pt;
}}

QFrame#Header {{
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1, stop:0 #15203b, stop:1 #1f3b73);
    border: 1px solid #2a4b8d;
    border-radius: {cls.px(6)}px;
}}
QLabel#HeaderTitle {{
    font-size: {f(14)}pt;
    font-weight: 700;
    color: white;
}}
QLabel#HeaderTag {{
    background: #1a2744;
    border: 1px solid #3a5080;
    border-radius: {cls.px(3)}px;
    padding: {cls.px(1)}px {cls.px(5)}px;
    color: #8cf0d8;
    font-size: {f(8)}pt;
    font-weight: 600;
}}

QGroupBox {{
    border: 1px solid #24314f;
    border-radius: {cls.px(6)}px;
    margin-top: {cls.px(6)}px;
    padding-top: {cls.px(8)}px;
    background: #0f172a;
    font-weight: 600;
    font-size: {f(9)}pt;
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    left: {cls.px(6)}px;
    padding: 0 {cls.px(3)}px;
    color: #b0c0e0;
    font-size: {f(8)}pt;
}}

QLabel#ParamLabel {{
    font-size: {f(8)}pt;
    color: #b0bdd0;
}}
QLabel#MutedLabel {{
    color: #6070a0;
    font-size: {f(8)}pt;
}}

QLineEdit {{
    background: #111a31;
    border: 1px solid #2a3d60;
    border-radius: {cls.px(4)}px;
    padding: {cls.px(2)}px {cls.px(4)}px;
    color: white;
    font-size: {f(9)}pt;
}}
QLineEdit#ReadOnlyField {{
    background: #0a1020;
    color: #506080;
    border-color: #1a2840;
}}

QPushButton {{
    background: #17233f;
    border: 1px solid #31466d;
    border-radius: {cls.px(5)}px;
    padding: {cls.px(3)}px {cls.px(7)}px;
    font-weight: 600;
    font-size: {f(8)}pt;
}}
QPushButton:hover {{ background: #1e2f56; }}
QPushButton[primary="true"] {{ background: #2553d8; border-color: #4f7eff; }}
QPushButton[primary="true"]:hover {{ background: #2d62fb; }}
QPushButton[danger="true"] {{ background: #7e2431; border-color: #b44b5b; }}
QPushButton[danger="true"]:hover {{ background: #9a2c3c; }}

QProgressBar {{
    border: 1px solid #2a3d60;
    border-radius: {cls.px(4)}px;
    background: #111a31;
    text-align: center;
    font-size: {f(7)}pt;
}}
QProgressBar::chunk {{
    background: #2dd4bf;
    border-radius: {cls.px(3)}px;
}}

QTextEdit {{
    background: #050812;
    border: 1px solid #1e2e4a;
    border-radius: {cls.px(4)}px;
    color: #8ff7c8;
    font-family: Consolas, Menlo, monospace;
    font-size: {f(8)}pt;
    padding: {cls.px(2)}px;
}}

QFrame#MetricCard {{
    background: #111a31;
    border: 1px solid #1e2e4a;
    border-radius: {cls.px(4)}px;
    max-height: {cls.px(42)}px;
}}
QLabel#MetricTitle {{
    color: #607898;
    font-size: {f(7)}pt;
}}
QLabel#MetricValue {{
    background: #030712;
    border-radius: {cls.px(3)}px;
    padding: {cls.px(1)}px {cls.px(2)}px;
    color: #7fffd4;
    font-size: {f(10)}pt;
    font-weight: 700;
    font-family: Consolas, Menlo, monospace;
}}

QRadioButton {{
    spacing: {cls.px(3)}px;
    color: #e6edf7;
    font-size: {f(8)}pt;
    padding: {cls.px(1)}px {cls.px(2)}px;
}}
QRadioButton::indicator {{
    width: {cls.px(10)}px; height: {cls.px(10)}px;
    border-radius: {cls.px(5)}px;
    border: 2px solid #7aa2ff;
    background: #09101f;
}}
QRadioButton::indicator:checked {{ border-color: #2dd4bf; background: #2dd4bf; }}
QRadioButton:checked {{ color: #2dd4bf; font-weight: 700; }}

QComboBox {{
    background: #111a31;
    border: 1px solid #2a3d60;
    border-radius: {cls.px(4)}px;
    padding: {cls.px(2)}px {cls.px(4)}px;
    color: white;
    font-size: {f(8)}pt;
}}
QComboBox::drop-down {{ border: none; width: {cls.px(14)}px; }}
QComboBox QAbstractItemView {{
    background: #111a31;
    border: 1px solid #2a3d60;
    color: white;
    selection-background-color: #2553d8;
}}

QSplitter::handle {{ background: #1a2844; width: {cls.px(3)}px; height: {cls.px(3)}px; }}
QSplitter::handle:hover {{ background: #3a5888; }}

QStatusBar {{ font-size: {f(8)}pt; color: #6080a0; }}
"""

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

# Paket-Registry (importiert beim ersten Zugriff)
from package_registry import (
    PackageInfo, BackendType, PackageStatus,
    discover_packages, resolve_backend, get_default_package,
    generate_cpp_config,
)

# ─── Parameter ───────────────────────────────────────────────
# (key, label, unit, default)

GEOM_PARAMS = [
    ("R",     "R",       "m",   0.02),
    ("L",     "L",       "m",   0.04),
    ("Nw",    "Nw",      "--",  6.0),
    ("R_ohm", "R_ohm",   "Ohm", 0.36),
    ("frequency", "f_RF", "Hz", 2.5e6),
    ("lc",    "l_coil",  "m",   0.04),
    ("Rc",    "R_coil",  "m",   0.02),
    ("Vgrid", "V_grid",  "V",   1500.0),
    ("sgrid", "s_grid",  "m",   0.001),
    ("betai", "betai",   "--",  0.5),
    ("betag", "betag",   "--",  0.05145),
]

OPERATE_PARAMS = [
    ("P_RFG",     "P_RFG",  "W",    18.0),
    ("P_RFG_max", "P_max",  "W",    80.0),
    ("Q0sccm",    "Q0",     "sccm", 0.475),
    ("I_soll",    "I_soll", "mA",   15.0),
    ("density_profile_factor", "alpha_p", "--", 1.0),
]

SWEEP_PARAMS = [
    ("Q0sccm_start", "Start", "sccm", 0.27),
    ("Q0sccm_step",  "Step",  "sccm", 0.01),
    ("jjmax",        "N",     "--",   73),
]

ALL_PRIMARY: dict[str, tuple] = {}
for key, label, unit, default in GEOM_PARAMS + OPERATE_PARAMS:
    ALL_PRIMARY[key] = (label, unit, default)

ALL_SWEEP: dict[str, tuple] = {}
for key, label, unit, default in SWEEP_PARAMS:
    ALL_SWEEP[key] = (label, unit, default)

# ─── Validierung ─────────────────────────────────────────────

VALIDATION_RULES = [
    ("R",         lambda v: v > 0,       "R > 0"),
    ("L",         lambda v: v > 0,       "L > 0"),
    ("betai",     lambda v: 0 < v <= 1,  "betai in (0,1]"),
    ("betag",     lambda v: 0 <= v <= 1, "betag in [0,1]"),
    ("frequency", lambda v: v > 0,       "f_RF > 0"),
    ("Nw",        lambda v: v >= 1,      "Nw >= 1"),
    ("R_ohm",     lambda v: v > 0,       "R_ohm > 0"),
    ("Rc",        lambda v: v > 0,       "Rc > 0"),
    ("lc",        lambda v: v > 0,       "lc > 0"),
    ("Vgrid",     lambda v: v > 0,       "Vgrid > 0"),
    ("sgrid",     lambda v: v > 0,       "sgrid > 0"),
    ("P_RFG",     lambda v: v > 0,       "P_RFG > 0"),
    ("P_RFG_max", lambda v: v > 0,       "P_max > 0"),
    ("Q0sccm",    lambda v: v >= 0,      "Q0 >= 0"),
    ("I_soll",    lambda v: v > 0,       "I_soll > 0"),
    ("Q0sccm_start", lambda v: v >= 0,   "Start >= 0"),
    ("Q0sccm_step",  lambda v: v > 0,    "Step > 0"),
    ("jjmax",     lambda v: v >= 1,      "N >= 1"),
    ("density_profile_factor", lambda v: 0 < v <= 1.0, "alpha_p in (0,1]"),
]


def cross_validate(vals):
    errs = []
    if "Rc" in vals and "R" in vals and vals["Rc"] < vals["R"]:
        errs.append(f"Rc < R ({vals['Rc']:.4f} < {vals['R']:.4f})")
    if "P_RFG" in vals and "P_RFG_max" in vals and vals["P_RFG"] > vals["P_RFG_max"]:
        errs.append(f"P_RFG > P_max")
    return errs


# ─── Helpers ─────────────────────────────────────────────────


def _field(grid, row, col, label, unit, default, entries, key, width=None):
    if width is None:
        width = UI.px(68)
    base = col * 3
    lbl = QLabel(label)
    lbl.setObjectName("ParamLabel")
    grid.addWidget(lbl, row, base)
    ed = QLineEdit(str(default))
    ed.setAlignment(Qt.AlignmentFlag.AlignRight)
    ed.setFixedWidth(width)
    grid.addWidget(ed, row, base + 1)
    u = QLabel(unit)
    u.setObjectName("MutedLabel")
    u.setFixedWidth(UI.px(28))
    grid.addWidget(u, row, base + 2)
    entries[key] = ed


class MetricCard(QFrame):
    def __init__(self, title, unit):
        super().__init__()
        self.setObjectName("MetricCard")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(3, 2, 3, 2)
        lay.setSpacing(0)
        t = QLabel(f"{title} [{unit}]")
        t.setObjectName("MetricTitle")
        self.val = QLabel("--")
        self.val.setObjectName("MetricValue")
        self.val.setAlignment(Qt.AlignmentFlag.AlignCenter)
        lay.addWidget(t)
        lay.addWidget(self.val)

    def set_value(self, v):
        self.val.setText(v)


# ─── LivePlotGrid ────────────────────────────────────────────

class LivePlotGrid(QWidget):
    def __init__(self):
        super().__init__()
        layout = QGridLayout(self)
        layout.setHorizontalSpacing(3)
        layout.setVerticalSpacing(3)
        layout.setContentsMargins(0, 0, 0, 0)

        pg.setConfigOptions(antialias=True, background="#0b1220", foreground="#d7e3ff")

        self.data = {
            "first_point": None,
            "Q0sccm": [], "P_RFG": [], "P_abs_rf": [], "Te": [], "Tg": [],
            "n": [], "ng": [], "I_mA": [],
            # Iod-Spezies (gefuellt wenn IODINE_EXT vorhanden)
            "diss": [], "fIp": [], "alpha": [],
            "nI": [], "nI2": [], "nIp": [], "nI2p": [], "nIm": [],
        }
        self._fail_q0 = []
        self._fail_P = []
        self._fail_I = []
        self._I_target = None

        # 3x3 Grid: Kern-Plots + Iod-Diagnostik
        defs = [
            # Zeile 1: Kern
            ("P_RFG", "P [W]"),       ("Te", "Te [eV]"),     ("I_mA", "I [mA]"),
            # Zeile 2: Dichten + Temperatur
            ("n", "n_e [m-3]"),       ("ng", "n_g [m-3]"),   ("Tg", "Tg [K]"),
            # Zeile 3: Iod-Diagnostik (leer fuer Xenon)
            ("diss", "Dissoziation"), ("fIp", "f(I+) / f(I2+)"), ("alpha", "I-/n_e"),
        ]

        self.widgets = {}
        self.curves = {}
        self.markers = {}
        self.first_markers = {}

        self._fail_scatters = {}
        # Zusaetzliche Kurven fuer Multi-Spezies-Plots
        self._extra_curves = {}

        for idx, (key, title) in enumerate(defs):
            w = pg.PlotWidget()
            w.setMinimumHeight(UI.px(55))  # Kompakter fuer 3x3
            w.showGrid(x=True, y=True, alpha=0.15)
            w.setTitle(title, size=f"{UI.pt(8)}pt")
            w.setLabel("bottom", "Q0 [sccm]")
            ax_l = w.getPlotItem().getAxis("left")
            ax_b = w.getPlotItem().getAxis("bottom")
            ax_l.enableAutoSIPrefix(True)
            ax_b.enableAutoSIPrefix(False)
            ax_l.setStyle(tickTextWidth=UI.px(45))
            tick_font = pg.QtGui.QFont("Consolas", round(UI.pt(7)))
            for ax in (ax_l, ax_b):
                ax.setStyle(tickFont=tick_font)

            pen = pg.mkPen(color="#5cb8ff", width=2)
            # P-Plot Hauptkurve bekommt Namen fuer Legende
            curve_name = "P_RFG" if key == "P_RFG" else None
            curve = w.plot([], [], pen=pen, symbol="o", symbolSize=3,
                           symbolBrush="#5cb8ff", symbolPen=None, name=curve_name)
            marker = pg.ScatterPlotItem(size=5, brush=pg.mkBrush("w"), pen=pg.mkPen("w"))
            fm = pg.ScatterPlotItem(size=7, brush=pg.mkBrush("#ff6b6b"), pen=pg.mkPen("#ff6b6b"))
            w.addItem(marker)
            w.addItem(fm)
            fail_scatter = pg.ScatterPlotItem(
                size=7, symbol="x",
                brush=pg.mkBrush("#ff4444"), pen=pg.mkPen("#ff4444", width=2))
            w.addItem(fail_scatter)

            self.widgets[key] = w
            self.curves[key] = curve
            self.markers[key] = marker
            self.first_markers[key] = fm
            self._fail_scatters[key] = fail_scatter
            layout.addWidget(w, idx // 3, idx % 3)

            # P-Plot: zweite Kurve fuer P_abs (RF-gekoppelt)
            if key == "P_RFG":
                pen_pabs = pg.mkPen(color="#ff9f43", width=2.5, style=pg.QtCore.Qt.PenStyle.DashLine)
                c_pabs = w.plot([], [], pen=pen_pabs, symbol="s", symbolSize=5,
                                symbolBrush="#ff9f43", symbolPen=None, name="P_abs")
                self._extra_curves["P_abs_on_P"] = c_pabs
                # Legende im P-Plot
                w.addLegend(offset=(5, 5), labelTextSize=f"{UI.pt(7)}pt")

            # fIp-Plot: zweite Kurve fuer f(I2+) in anderer Farbe
            if key == "fIp":
                pen2 = pg.mkPen(color="#ff9f43", width=2, style=pg.QtCore.Qt.PenStyle.DashLine)
                c2 = w.plot([], [], pen=pen2, symbol="s", symbolSize=3,
                            symbolBrush="#ff9f43", symbolPen=None, name="f(I2+)")
                self._extra_curves["fI2p"] = c2
                # Legende
                w.addLegend(offset=(5, 5), labelTextSize=f"{UI.pt(7)}pt")
                self.curves[key].setData([], [], name="f(I+)")

            # ng-Plot: zweite Kurve fuer Iod-Neutraldichten
            if key == "ng":
                pen_i = pg.mkPen(color="#34d399", width=1.5, style=pg.QtCore.Qt.PenStyle.DashLine)
                c_i = w.plot([], [], pen=pen_i, name="n(I)")
                self._extra_curves["nI_on_ng"] = c_i
                pen_i2 = pg.mkPen(color="#ff9f43", width=1.5, style=pg.QtCore.Qt.PenStyle.DashLine)
                c_i2 = w.plot([], [], pen=pen_i2, name="n(I2)")
                self._extra_curves["nI2_on_ng"] = c_i2

        # I_target Referenzlinie im I_mA-Plot
        self._I_target_line = None
        if "I_mA" in self.widgets:
            self._I_target_line = pg.InfiniteLine(
                angle=0, pen=pg.mkPen("#ffaa44", width=1, style=pg.QtCore.Qt.PenStyle.DashLine))
            self._I_target_line.setVisible(False)
            self.widgets["I_mA"].addItem(self._I_target_line)

    def clear(self):
        self.data["first_point"] = None
        for k in self.data:
            if k != "first_point":
                self.data[k].clear()
        self._fail_q0.clear()
        self._fail_P.clear()
        self._fail_I.clear()
        self._I_target = None
        if self._I_target_line:
            self._I_target_line.setVisible(False)
        self._redraw()

    def set_I_target(self, I_target):
        """Setze I_soll-Referenzlinie im I_mA-Plot."""
        self._I_target = I_target
        if self._I_target_line:
            self._I_target_line.setValue(I_target)
            self._I_target_line.setVisible(True)

    def add_point(self, q0, p, te, tg, n, ng, i, iodine=None, P_abs_rf=0.0):
        """Fuge Sweep-Punkt hinzu. P_abs_rf fuer RF-gekoppelten Modus."""
        vals = {"Q0sccm": q0, "P_RFG": p, "P_abs_rf": P_abs_rf, "Te": te, "Tg": tg, "n": n, "ng": ng, "I_mA": i}
        if iodine:
            for k in ("diss", "fIp", "alpha", "nI", "nI2", "nIp", "nI2p", "nIm"):
                vals[k] = iodine.get(k, 0)
            vals["fIp"] = iodine.get("fIp", 0)

        if self.data["first_point"] is None:
            self.data["first_point"] = vals
        else:
            for k, v in vals.items():
                if k in self.data:
                    self.data[k].append(v)
        self._redraw()

    def update_last_iodine(self, iod):
        """Aktualisiere den letzten Punkt mit Iod-Daten (kommt nach RESULT)."""
        for k in ("diss", "fIp", "alpha", "nI", "nI2", "nIp", "nI2p", "nIm"):
            if k in self.data and k in iod:
                if self.data["Q0sccm"]:  # Nicht-erster Punkt
                    self.data[k].append(iod[k])
                elif self.data["first_point"] is not None:
                    self.data["first_point"][k] = iod[k]
        self._redraw()

    def set_power_label(self, label: str):
        """Update the P-plot title (e.g. 'P_abs [W]' or 'P_RFG [W]')."""
        if "P_RFG" in self.widgets:
            self.widgets["P_RFG"].setTitle(label, size=f"{UI.pt(8)}pt")

    def set_flow_unit(self, unit: str):
        """Update x-axis label on all plots (e.g. 'Q0 [sccm]' or 'Q0 [mg/s]')."""
        label = f"Q0 [{unit}]"
        for w in self.widgets.values():
            w.setLabel("bottom", label)

    def add_fail_point(self, q0, p_best, i_best):
        """Fehlgeschlagener I-fix-Punkt (wird als X-Marker geplottet)."""
        self._fail_q0.append(q0)
        self._fail_P.append(p_best)
        self._fail_I.append(i_best)
        self._redraw()

    def _redraw(self):
        x = self.data["Q0sccm"]
        fp = self.data["first_point"]  # Dict oder None
        for key in self.curves:
            yd = self.data.get(key, [])
            # Laengenpruefung: nur plotten wenn Arrays gleich lang
            if len(x) == len(yd):
                self.curves[key].setData(x, yd)
                self.markers[key].setData([{"pos": (x[-1], yd[-1])}] if x and yd else [])
            else:
                self.curves[key].setData([], [])
                self.markers[key].setData([])
            if fp and key in fp:
                q0f = fp.get("Q0sccm", 0)
                self.first_markers[key].setData([{"pos": (q0f, fp[key])}])
            else:
                self.first_markers[key].setData([])

            # Fail-Marker
            fs = self._fail_scatters.get(key)
            if fs:
                if key == "P_RFG" and self._fail_q0:
                    fs.setData([{"pos": (fq, fp_)} for fq, fp_ in zip(self._fail_q0, self._fail_P)])
                elif key == "I_mA" and self._fail_q0:
                    fs.setData([{"pos": (fq, fi)} for fq, fi in zip(self._fail_q0, self._fail_I)])
                else:
                    fs.setData([])

        # Extra-Kurve: P_abs im P-Plot
        if "P_abs_on_P" in self._extra_curves:
            pabs_data = self.data.get("P_abs_rf", [])
            if len(x) == len(pabs_data) and any(v > 0 for v in pabs_data):
                self._extra_curves["P_abs_on_P"].setData(x, pabs_data)
            else:
                self._extra_curves["P_abs_on_P"].setData([], [])

        # Extra-Kurven: Iod f(I2+) im fIp-Plot
        if "fI2p" in self._extra_curves:
            fIp_data = self.data.get("fIp", [])
            if len(x) == len(fIp_data):
                fI2p_data = [1.0 - v for v in fIp_data]
                self._extra_curves["fI2p"].setData(x, fI2p_data)

        # Extra-Kurven: n(I) und n(I2) im ng-Plot
        nI_data = self.data.get("nI", [])
        nI2_data = self.data.get("nI2", [])
        if "nI_on_ng" in self._extra_curves and len(x) == len(nI_data):
            self._extra_curves["nI_on_ng"].setData(x, nI_data)
        if "nI2_on_ng" in self._extra_curves and len(x) == len(nI2_data):
            self._extra_curves["nI2_on_ng"].setData(x, self.data.get("nI2", []))


# ═════════════════════════════════════════════════════════════
# Hauptfenster
# ═════════════════════════════════════════════════════════════

class SimulatorWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Global Plasma Model")
        screen = QApplication.primaryScreen()
        if screen:
            av = screen.availableGeometry()
            w = min(1600, int(av.width() * 0.94))
            h = min(1000, int(av.height() * 0.92))
            self.resize(w, h)
            self.move(av.x() + (av.width() - w) // 2, av.y() + (av.height() - h) // 2)
        else:
            self.resize(1600, 1000)

        self.entries: dict[str, QLineEdit] = {}
        self.sweep_entries: dict[str, QLineEdit] = {}
        self.metric_cards: dict[str, MetricCard] = {}
        self.process: QProcess | None = None
        self.cancel_requested = False
        self._stdout_buffer = ""
        self.current_q0 = None
        self._q0_unit_mode = "sccm"
        self._q0_internal_sccm = 0.475

        self._build_ui()
        self._apply_styles()

    def _build_ui(self):
        root = QWidget()
        self.setCentralWidget(root)
        outer = QVBoxLayout(root)
        outer.setContentsMargins(4, 3, 4, 3)
        outer.setSpacing(3)

        # Header
        hdr = QFrame()
        hdr.setObjectName("Header")
        hl = QHBoxLayout(hdr)
        hl.setContentsMargins(10, 5, 10, 5)
        t = QLabel("Global Plasma Model")
        t.setObjectName("HeaderTitle")
        hl.addWidget(t)
        hl.addStretch(1)
        self._hdr_pkg = QLabel("--")
        self._hdr_pkg.setObjectName("HeaderTag")
        self._hdr_rate = QLabel("Legacy")
        self._hdr_rate.setObjectName("HeaderTag")
        self._hdr_mode = QLabel("I-fix")
        self._hdr_mode.setObjectName("HeaderTag")
        for w in (self._hdr_pkg, self._hdr_rate, self._hdr_mode):
            hl.addWidget(w)
        outer.addWidget(hdr)

        # Main splitter
        self.main_split = QSplitter(Qt.Orientation.Horizontal)
        self.main_split.setChildrenCollapsible(False)
        outer.addWidget(self.main_split, 1)

        # ── LEFT PANEL ───────────────────────────────────────
        left_scroll = QScrollArea()
        left_scroll.setWidgetResizable(True)
        left_scroll.setFrameShape(QFrame.Shape.NoFrame)
        left_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        left_panel = QWidget()
        left_panel.setMaximumWidth(UI.px(380))
        left_scroll.setWidget(left_panel)
        ll = QVBoxLayout(left_panel)
        ll.setContentsMargins(0, 0, 2, 0)
        ll.setSpacing(3)

        ll.addWidget(self._build_geometry())
        ll.addWidget(self._build_operate())
        ll.addWidget(self._build_diagnosis())
        ll.addWidget(self._build_sweep())
        ll.addWidget(self._build_actions())
        ll.addWidget(self._build_progress())
        ll.addWidget(self._build_metrics())
        ll.addStretch(1)

        # ── RIGHT PANEL (plots + log) ────────────────────────
        right_w = QWidget()
        rl = QVBoxLayout(right_w)
        rl.setContentsMargins(2, 0, 0, 0)
        rl.setSpacing(2)

        self.plot_grid = LivePlotGrid()
        self.plot_grid.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        log_box = QGroupBox("Log")
        log_lay = QVBoxLayout(log_box)
        log_lay.setContentsMargins(2, 10, 2, 2)
        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self.log.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.log.setMaximumHeight(UI.px(120))
        log_lay.addWidget(self.log)

        vsplit = QSplitter(Qt.Orientation.Vertical)
        vsplit.setChildrenCollapsible(False)
        vsplit.addWidget(self.plot_grid)
        vsplit.addWidget(log_box)
        vsplit.setStretchFactor(0, 8)   # Mehr Platz fuer Plots
        vsplit.setStretchFactor(1, 1)   # Log kompakter
        rl.addWidget(vsplit, 1)

        self.main_split.addWidget(left_scroll)
        self.main_split.addWidget(right_w)
        self.main_split.setStretchFactor(0, 0)  # left: no stretch
        self.main_split.setStretchFactor(1, 1)  # right: all stretch
        self.main_split.setSizes([UI.px(320), UI.px(900)])

        self.setStatusBar(QStatusBar())
        self.set_status("Bereit.")

        for txt, fn in [("Defaults", self.reset_defaults), ("Compile", self.compile_cpp),
                         ("Run", self.run_sim), ("Cancel", self.cancel_sim)]:
            a = QAction(txt, self)
            a.triggered.connect(fn)
            self.addAction(a)

    # ── System / Geometrie ───────────────────────────────────

    def _build_geometry(self):
        box = QGroupBox("System / Geometrie")
        g = QGridLayout(box)
        g.setContentsMargins(4, 12, 4, 3)
        g.setHorizontalSpacing(4)
        g.setVerticalSpacing(2)

        # ── Zeile 0: Chemiepaket ──────────────────────────────
        g.addWidget(QLabel("Paket"), 0, 0)
        self.cmb_package = QComboBox()
        self.cmb_package.setToolTip("Chemiepaket (Gas + Datenquelle + Backend)")
        self.cmb_package.currentIndexChanged.connect(self._on_package_changed)
        g.addWidget(self.cmb_package, 0, 1, 1, 5)

        # ── Zeile 1: Thruster-Preset ─────────────────────────
        g.addWidget(QLabel("Thruster"), 1, 0)
        self.cmb_thruster = QComboBox()
        self.cmb_thruster.setToolTip("Thruster-Konfiguration (Literatur-Preset oder Custom)")
        self.cmb_thruster.currentIndexChanged.connect(self._on_thruster_changed)
        g.addWidget(self.cmb_thruster, 1, 1, 1, 5)
        self._load_thruster_presets()

        # ── Zeile 2b: Preset-Info ────────────────────────────
        self._lbl_preset_info = QLabel("")
        self._lbl_preset_info.setObjectName("MutedLabel")
        self._lbl_preset_info.setWordWrap(True)
        g.addWidget(self._lbl_preset_info, 2, 0, 1, 6)

        # ── Zeile 3: Ratenmodell + Info ───────────────────────
        g.addWidget(QLabel("Raten"), 3, 0)
        self.cmb_rate_model = QComboBox()
        for v, t in RATE_MODEL_OPTIONS:
            self.cmb_rate_model.addItem(t, v)
        self.cmb_rate_model.setCurrentIndex(DEFAULT_RATE_MODEL)
        self.cmb_rate_model.setToolTip("Legacy / Conservative / Full tabulated")
        self.cmb_rate_model.currentIndexChanged.connect(self._sync_hdr)
        g.addWidget(self.cmb_rate_model, 3, 1, 1, 2)

        self._lbl_backend = QLabel("")
        self._lbl_backend.setObjectName("MutedLabel")
        g.addWidget(self._lbl_backend, 3, 3, 1, 3)

        # ── Zeile 4+: Geometrie ──────────────────────────────
        pairs = [
            (GEOM_PARAMS[0], GEOM_PARAMS[1]),    # R, L
            (GEOM_PARAMS[2], GEOM_PARAMS[3]),    # Nw, R_ohm
            (GEOM_PARAMS[4], GEOM_PARAMS[5]),    # f_RF, l_coil
            (GEOM_PARAMS[6], GEOM_PARAMS[7]),    # R_coil, Vgrid
            (GEOM_PARAMS[8], GEOM_PARAMS[9]),    # sgrid, betai
            (GEOM_PARAMS[10],),                   # betag
        ]
        for ri, pair in enumerate(pairs, 4):
            for ci, (k, lb, un, df) in enumerate(pair):
                _field(g, ri, ci, lb, un, df, self.entries, k)

        # Pakete laden
        self._packages = discover_packages(SCRIPT_DIR)
        self._populate_packages()

        return box

    def _populate_packages(self):
        """Befuelle den Paket-Combo aus der Registry."""
        self.cmb_package.blockSignals(True)
        self.cmb_package.clear()
        for pkg in self._packages:
            label = pkg.display_name + pkg.status_label
            self.cmb_package.addItem(label, pkg.id)
        # Default: xenon_biagi
        default = get_default_package(self._packages)
        if default:
            for i in range(self.cmb_package.count()):
                if self.cmb_package.itemData(i) == default.id:
                    self.cmb_package.setCurrentIndex(i)
                    break
        self.cmb_package.blockSignals(False)
        self._on_package_changed()

    def _current_package(self) -> PackageInfo | None:
        """Aktuell gewaehltes Paket."""
        pkg_id = self.cmb_package.currentData() if hasattr(self, "cmb_package") else None
        if pkg_id:
            for p in self._packages:
                if p.id == pkg_id:
                    return p
        return None

    def _on_package_changed(self, _idx=None):
        """Reagiert auf Paketwechsel: Backend-Info, Rate-Model-Verfuegbarkeit."""
        pkg = self._current_package()
        if not pkg:
            return

        # Backend-Info anzeigen
        backend_str = "C++" if pkg.is_cpp else "Python"
        status_str = {
            PackageStatus.PRODUCTION: "",
            PackageStatus.EXPERIMENTAL: " (experimentell)",
            PackageStatus.DEMO: " (Demo)",
        }.get(pkg.status, "")
        self._lbl_backend.setText(f"{backend_str}{status_str}")

        # Rate-Model nur fuer C++-Pakete sinnvoll
        self.cmb_rate_model.setEnabled(pkg.is_cpp)
        if not pkg.is_cpp:
            self.cmb_rate_model.setToolTip("Ratenmodell nur fuer C++-Pakete relevant")

        self._sync_hdr()
        if hasattr(self, "_lbl_diag_compat"):
            self._update_diagnosis()
        if hasattr(self, "_q0_alt_lbl"):
            self._update_q0_display()
        self._update_power_label()

    # ── Thruster-Presets ────────────────────────────────────

    def _load_thruster_presets(self):
        """Lade Thruster-Presets aus JSON."""
        import json
        self._thruster_presets = []
        preset_path = SCRIPT_DIR / "thruster_presets.json"
        if preset_path.exists():
            try:
                with open(preset_path, encoding="utf-8") as f:
                    data = json.load(f)
                self._thruster_presets = data.get("presets", [])
            except (json.JSONDecodeError, OSError):
                pass

        self.cmb_thruster.blockSignals(True)
        self.cmb_thruster.clear()
        for p in self._thruster_presets:
            label = p["name"]
            if p.get("description"):
                label += f"  ({p['description'][:40]})"
            self.cmb_thruster.addItem(label, p["id"])
        # Default: custom (letzter Eintrag)
        for i in range(self.cmb_thruster.count()):
            if self.cmb_thruster.itemData(i) == "custom":
                self.cmb_thruster.setCurrentIndex(i)
                break
        self.cmb_thruster.blockSignals(False)

    def _on_thruster_changed(self, _idx=None):
        """Lade Parameter aus dem gewaehlten Thruster-Preset."""
        preset_id = self.cmb_thruster.currentData()
        if not preset_id or preset_id == "custom":
            return  # Custom = keine Aenderung

        preset = next((p for p in self._thruster_presets if p["id"] == preset_id), None)
        if not preset:
            return

        params = preset.get("params", {})
        for key, value in params.items():
            if key in self.entries:
                self.entries[key].setText(str(value))

        # eta_opt in Metadaten speichern (fuer Extraktionsmodell)
        ext = preset.get("extraction", {})
        self._eta_opt = ext.get("eta_opt", 0.25)

        # Preset-Info anzeigen
        plaus = preset.get("plausibility", {})
        rec = preset.get("recommended_gas", [])
        info_parts = [f"eta_opt={self._eta_opt}"]
        if rec:
            info_parts.append("Gas: " + ", ".join(r.capitalize() for r in rec))
        I_min = plaus.get("I_beam_typ_min_mA")
        I_max = plaus.get("I_beam_typ_max_mA")
        if I_min is not None and I_max is not None:
            info_parts.append(f"I_beam ~{I_min}-{I_max} mA")
        if hasattr(self, "_lbl_preset_info"):
            self._lbl_preset_info.setText("  ".join(info_parts))

        self.log_msg(f"Preset: {preset['name']} geladen (eta_opt={self._eta_opt})")
        self._sync_hdr()
        if hasattr(self, "_lbl_diag_compat"):
            self._update_diagnosis()

    # ── Betrieb ──────────────────────────────────────────────

    def _build_operate(self):
        box = QGroupBox("Betrieb")
        g = QGridLayout(box)
        g.setContentsMargins(4, 12, 4, 3)
        g.setHorizontalSpacing(4)
        g.setVerticalSpacing(2)

        self.mode_group = QButtonGroup(self)
        self.mode_buttons = {}
        for i, (v, t) in enumerate(SOLVE_MODE_OPTIONS):
            btn = QRadioButton(t)
            if v == DEFAULT_SOLVE_MODE:
                btn.setChecked(True)
            self.mode_group.addButton(btn, v)
            self.mode_buttons[v] = btn
            g.addWidget(btn, 0, i * 3, 1, 3)
        self.mode_group.idToggled.connect(self._on_mode)

        # P-Feld mit dynamischem Label (P_abs fuer Iod direct, P_RFG fuer Xe/RF)
        self._p_label = QLabel("P_RFG")
        self._p_label.setObjectName("ParamLabel")
        g.addWidget(self._p_label, 1, 0)
        ed_p = QLineEdit("18.0")
        ed_p.setAlignment(Qt.AlignmentFlag.AlignRight)
        ed_p.setFixedWidth(UI.px(68))
        g.addWidget(ed_p, 1, 1)
        ul_p = QLabel("W")
        ul_p.setObjectName("MutedLabel")
        ul_p.setFixedWidth(UI.px(28))
        g.addWidget(ul_p, 1, 2)
        self.entries["P_RFG"] = ed_p

        _field(g, 1, 1, "P_max", "W", 80.0, self.entries, "P_RFG_max")

        # RF coupling checkbox (nur relevant fuer Python/Iod)
        self._chk_rf_coupling = QCheckBox("RF-Kopplung")
        self._chk_rf_coupling.setToolTip(
            "Aktiviert: P = P_RFG (Generator), P_abs wird via RF-Modell abgeleitet.\n"
            "Deaktiviert: P = P_abs (direkt absorbierte Leistung).")
        self._chk_rf_coupling.toggled.connect(self._on_rf_coupling_changed)
        g.addWidget(self._chk_rf_coupling, 1, 5)

        # ── Q0 mit Einheitenumschaltung ──────────────────────
        g.addWidget(QLabel("Q0"), 2, 0)
        self._q0_edit = QLineEdit("0.475")
        self._q0_edit.setAlignment(Qt.AlignmentFlag.AlignRight)
        self._q0_edit.setFixedWidth(UI.px(68))
        g.addWidget(self._q0_edit, 2, 1)
        self.entries["Q0sccm"] = self._q0_edit  # Intern immer sccm

        self._q0_unit_combo = QComboBox()
        self._q0_unit_combo.addItems(["sccm", "mg/s"])
        self._q0_unit_combo.setFixedWidth(UI.px(48))
        self._q0_unit_combo.currentIndexChanged.connect(self._on_q0_unit_changed)
        g.addWidget(self._q0_unit_combo, 2, 2)

        self._q0_alt_lbl = QLabel("")
        self._q0_alt_lbl.setObjectName("MutedLabel")
        g.addWidget(self._q0_alt_lbl, 2, 3, 1, 3)
        self._q0_edit.textChanged.connect(self._update_q0_display)
        self._q0_unit_mode = "sccm"  # Aktive Einheit
        self._q0_internal_sccm = 0.475  # Source of truth

        _field(g, 3, 0, "I_soll", "mA", 15.0, self.entries, "I_soll")
        _field(g, 3, 1, "alpha_p", "--", 1.0, self.entries, "density_profile_factor")

        # I_soll-Aenderung -> Diagnose aktualisieren
        if "I_soll" in self.entries:
            self.entries["I_soll"].textChanged.connect(
                lambda: self._update_ifix_range_hint() if hasattr(self, "_lbl_diag_range") else None)

        self._on_mode(self.mode_group.checkedId(), True)
        self._update_q0_display()  # Initialen Korrespondenzwert berechnen
        return box

    def _on_rf_coupling_changed(self, checked):
        """Update P-Label when RF coupling checkbox changes."""
        self._update_power_label()

    def _update_power_label(self):
        """Set P-field label based on package type and RF coupling."""
        if not hasattr(self, "_p_label"):
            return
        pkg = self._current_package()
        is_python = pkg and pkg.is_python if pkg else False
        rf_on = self._chk_rf_coupling.isChecked() if hasattr(self, "_chk_rf_coupling") else False

        if is_python and not rf_on:
            self._p_label.setText("P_abs")
            self._p_label.setToolTip("Direkt absorbierte Leistung im Plasma [W]")
            if hasattr(self, "plot_grid"):
                self.plot_grid.set_power_label("P_abs [W]")
        else:
            self._p_label.setText("P_RFG")
            self._p_label.setToolTip("RF-Generator-Leistung [W]")
            if hasattr(self, "plot_grid"):
                self.plot_grid.set_power_label("P_RFG [W]")

    def _get_gas(self):
        pkg = self._current_package()
        return pkg.gas if pkg else "xenon"

    def _on_q0_unit_changed(self, idx):
        """Einheitenumschaltung sccm <-> mg/s."""
        try:
            from flow_units import sccm_to_mg_per_s, mg_per_s_to_sccm
        except ImportError:
            return
        new_unit = "sccm" if idx == 0 else "mg/s"
        if new_unit == self._q0_unit_mode:
            return
        # Aktuellen Wert konvertieren
        try:
            val = float(self._q0_edit.text())
        except ValueError:
            self._q0_unit_mode = new_unit
            return
        gas = self._get_gas()
        self._q0_edit.blockSignals(True)
        if new_unit == "mg/s":
            # sccm -> mg/s: Anzeige umstellen
            self._q0_internal_sccm = val
            mg = sccm_to_mg_per_s(val, gas)
            self._q0_edit.setText(f"{mg:.4f}")
        else:
            # mg/s -> sccm: zurueckrechnen
            sccm = mg_per_s_to_sccm(val, gas)
            self._q0_internal_sccm = sccm
            self._q0_edit.setText(f"{sccm:.4f}")
        self._q0_edit.blockSignals(False)
        self._q0_unit_mode = new_unit
        self._update_q0_display()
        # Sweep-Felder ebenfalls umstellen
        self._update_sweep_units(new_unit, gas)

    def _update_q0_display(self):
        """Aktualisiere Korrespondenzanzeige und interne sccm-Source-of-Truth."""
        try:
            from flow_units import sccm_to_mg_per_s, mg_per_s_to_sccm
        except ImportError:
            return
        gas = self._get_gas()
        try:
            val = float(self._q0_edit.text())
        except ValueError:
            self._q0_alt_lbl.setText("")
            return
        if self._q0_unit_mode == "sccm":
            self._q0_internal_sccm = val
            mg = sccm_to_mg_per_s(val, gas)
            self._q0_alt_lbl.setText(f"= {mg:.3f} mg/s")
        else:
            sccm = mg_per_s_to_sccm(val, gas)
            self._q0_internal_sccm = sccm
            self._q0_alt_lbl.setText(f"= {sccm:.4f} sccm")
        # Interne sccm-Source-of-Truth fuer den Solver bereitstellen
        self.entries["Q0sccm"].blockSignals(True)
        if self._q0_unit_mode != "sccm":
            # Wenn mg/s aktiv, das interne sccm-Feld heimlich aktualisieren
            pass  # entries["Q0sccm"] IST self._q0_edit, Solver liest sccm
        self.entries["Q0sccm"].blockSignals(False)

    def _update_sweep_units(self, unit, gas):
        """Sweep-Felder Start/Step/Ende bei Einheitenwechsel konvertieren."""
        try:
            from flow_units import sccm_to_mg_per_s, mg_per_s_to_sccm
        except ImportError:
            return
        for key in ("Q0sccm_start", "Q0sccm_step"):
            if key not in self.sweep_entries:
                continue
            ed = self.sweep_entries[key]
            try:
                val = float(ed.text())
            except ValueError:
                continue
            ed.blockSignals(True)
            if unit == "mg/s":
                ed.setText(f"{sccm_to_mg_per_s(val, gas):.4f}")
            else:
                ed.setText(f"{mg_per_s_to_sccm(val, gas):.4f}")
            ed.blockSignals(False)
        self._upd_q0end()
        # Update Sweep-Labels
        if hasattr(self, "_sweep_unit_labels"):
            for lbl in self._sweep_unit_labels:
                lbl.setText(unit)
        # Update Plot x-Achsen
        if hasattr(self, "plot_grid"):
            self.plot_grid.set_flow_unit(unit)

    def _on_mode(self, mid, checked=True):
        if not checked:
            return
        if "I_soll" in self.entries:
            self.entries["I_soll"].setEnabled(mid != 2)
        self._sync_hdr()
        if hasattr(self, "_lbl_diag_range"):
            self._update_diagnosis()

    # ── Diagnose / Plausibilisierung ────────────────────────

    def _build_diagnosis(self):
        """Kompaktes Diagnose-Panel: Preset-Kompatibilitaet + erreichbarer I_beam-Bereich."""
        box = QGroupBox("Diagnose")
        lay = QVBoxLayout(box)
        lay.setContentsMargins(4, 10, 4, 3)
        lay.setSpacing(2)

        self._lbl_diag_compat = QLabel("")
        self._lbl_diag_compat.setWordWrap(True)
        lay.addWidget(self._lbl_diag_compat)

        self._lbl_diag_range = QLabel("")
        self._lbl_diag_range.setWordWrap(True)
        lay.addWidget(self._lbl_diag_range)

        return box

    def _update_diagnosis(self):
        """Aktualisiere Diagnose-Panel: Kompatibilitaet + I_beam-Bereich."""
        self._update_gas_compatibility()
        self._update_ifix_range_hint()

    def _get_current_preset(self):
        preset_id = self.cmb_thruster.currentData() if hasattr(self, "cmb_thruster") else None
        if not preset_id or preset_id == "custom":
            return None
        return next((p for p in self._thruster_presets if p["id"] == preset_id), None)

    def _update_gas_compatibility(self):
        """Pruefe ob Paket-Gas und Preset-Empfehlung zusammenpassen."""
        if not hasattr(self, "_lbl_diag_compat"):
            return
        pkg = self._current_package()
        preset = self._get_current_preset()

        if not pkg or not preset:
            self._lbl_diag_compat.setText("")
            return

        rec = preset.get("recommended_gas", [])
        plaus = preset.get("plausibility", {})
        gas_warnings = plaus.get("gas_warnings", {})
        pkg_gas = pkg.gas  # "xenon", "iodine", etc.

        # Gasspezifische Warnung?
        if pkg_gas in gas_warnings and gas_warnings[pkg_gas]:
            self._lbl_diag_compat.setText(f"⚠ {gas_warnings[pkg_gas]}")
            self._lbl_diag_compat.setStyleSheet("color: #e8a838; font-size: 11px;")
            return

        if rec and pkg_gas not in rec:
            rec_str = ", ".join(r.capitalize() for r in rec)
            self._lbl_diag_compat.setText(
                f"⚠ Preset '{preset['name']}' empfohlen fuer {rec_str}, "
                f"aber Paket ist {pkg_gas.capitalize()}.")
            self._lbl_diag_compat.setStyleSheet("color: #e8a838; font-size: 11px;")
            return

        # Alles OK
        self._lbl_diag_compat.setText("")
        self._lbl_diag_compat.setStyleSheet("")

    def _update_ifix_range_hint(self):
        """Zeige den ungefaehr erreichbaren Strahlstrombereich im I-fix-Modus."""
        if not hasattr(self, "_lbl_diag_range"):
            return
        if self.solve_mode() != 1:
            self._lbl_diag_range.setText("")
            return

        preset = self._get_current_preset()
        if not preset:
            self._lbl_diag_range.setText("")
            return

        plaus = preset.get("plausibility", {})
        I_min = plaus.get("I_beam_typ_min_mA")
        I_max = plaus.get("I_beam_typ_max_mA")

        if I_min is None or I_max is None:
            self._lbl_diag_range.setText("")
            return

        try:
            I_soll = float(self.entries["I_soll"].text())
        except (ValueError, KeyError):
            self._lbl_diag_range.setText("")
            return

        range_str = f"Typischer Bereich: {I_min} - {I_max} mA"

        if I_soll < I_min:
            self._lbl_diag_range.setText(
                f"{range_str}\n"
                f"⚠ I_soll = {I_soll:.0f} mA liegt unterhalb des typischen Bereichs.")
            self._lbl_diag_range.setStyleSheet("color: #e8a838; font-size: 11px;")
        elif I_soll > I_max:
            self._lbl_diag_range.setText(
                f"{range_str}\n"
                f"⚠ I_soll = {I_soll:.0f} mA liegt oberhalb des typischen Bereichs. "
                f"Erreichbarkeit abhaengig von P_max und Q0.")
            self._lbl_diag_range.setStyleSheet("color: #e8a838; font-size: 11px;")
        else:
            self._lbl_diag_range.setText(
                f"{range_str}  —  I_soll = {I_soll:.0f} mA liegt im Bereich.")
            self._lbl_diag_range.setStyleSheet("color: #5a8a5a; font-size: 11px;")

    # ── Sweep ────────────────────────────────────────────────

    def _build_sweep(self):
        box = QGroupBox("Sweep")
        g = QGridLayout(box)
        g.setContentsMargins(4, 12, 4, 3)
        g.setHorizontalSpacing(4)
        g.setVerticalSpacing(2)

        self._sweep_unit_labels = []
        for i, (k, lb, un, df) in enumerate(SWEEP_PARAMS):
            _field(g, 0, i, lb, un, df, self.sweep_entries, k, width=UI.px(55))
            # Merke die Unit-Labels fuer spaetere Aktualisierung (nur sccm-Felder)
            if un == "sccm":
                # Unit-Label ist das 3. Widget in der Spalte (base+2)
                base = i * 3
                w = g.itemAtPosition(0, base + 2)
                if w and w.widget():
                    self._sweep_unit_labels.append(w.widget())

        g.addWidget(QLabel("Ende"), 1, 0)
        self._q0end = QLineEdit("--")
        self._q0end.setReadOnly(True)
        self._q0end.setFixedWidth(UI.px(55))
        self._q0end.setAlignment(Qt.AlignmentFlag.AlignRight)
        self._q0end.setObjectName("ReadOnlyField")
        g.addWidget(self._q0end, 1, 1)
        self._q0end_unit = QLabel("sccm")
        self._q0end_unit.setObjectName("MutedLabel")
        g.addWidget(self._q0end_unit, 1, 2)
        self._sweep_unit_labels.append(self._q0end_unit)

        for k in ("Q0sccm_start", "Q0sccm_step", "jjmax"):
            if k in self.sweep_entries:
                self.sweep_entries[k].textChanged.connect(self._upd_q0end)
        self._upd_q0end()
        return box

    def _upd_q0end(self):
        try:
            s = float(self.sweep_entries["Q0sccm_start"].text())
            d = float(self.sweep_entries["Q0sccm_step"].text())
            n = int(float(self.sweep_entries["jjmax"].text()))
            self._q0end.setText(f"{s + d * (n - 1):.3f}")
        except (ValueError, KeyError):
            self._q0end.setText("--")

    # ── Aktionen ─────────────────────────────────────────────

    def _build_actions(self):
        box = QGroupBox("Aktionen")
        lay = QVBoxLayout(box)
        lay.setContentsMargins(4, 12, 4, 3)
        lay.setSpacing(2)

        r1 = QHBoxLayout()
        r1.setSpacing(3)
        r1.addWidget(self._btn("Kompilieren", self.compile_cpp))
        r1.addWidget(self._btn("Start", self.run_sim, primary=True))
        r1.addWidget(self._btn("Stop", self.cancel_sim, danger=True))
        lay.addLayout(r1)

        r2 = QHBoxLayout()
        r2.setSpacing(3)
        r2.addWidget(self._btn("Defaults", self.reset_defaults, small=True))
        r2.addWidget(self._btn("Plots", self.clear_plots, small=True))
        r2.addWidget(self._btn("Output", self.open_output, small=True))
        r2.addWidget(self._btn("Logs", self.open_log_viewer, small=True))
        lay.addLayout(r2)

        r3 = QHBoxLayout()
        r3.setSpacing(3)
        r3.addWidget(self._btn("Physik-Daten", self._open_physics_data, small=True))
        r3.addStretch()
        lay.addLayout(r3)
        return box

    def _btn(self, text, slot, primary=False, danger=False, small=False):
        b = QPushButton(text)
        b.clicked.connect(slot)
        if primary:
            b.setProperty("primary", True)
        if danger:
            b.setProperty("danger", True)
        b.setMinimumHeight(UI.px(22) if small else UI.px(28))
        return b

    # ── Fortschritt (links) ──────────────────────────────────

    def _build_progress(self):
        box = QGroupBox("Status")
        g = QGridLayout(box)
        g.setContentsMargins(4, 12, 4, 3)
        g.setHorizontalSpacing(4)
        g.setVerticalSpacing(2)

        g.addWidget(QLabel("Scan"), 0, 0)
        self.scan_progress = QProgressBar()
        self.scan_progress.setRange(0, 100)
        self.scan_progress.setFixedHeight(UI.px(14))
        g.addWidget(self.scan_progress, 0, 1)
        self.scan_label = QLabel("--")
        self.scan_label.setFixedWidth(UI.px(45))
        self.scan_label.setObjectName("MutedLabel")
        g.addWidget(self.scan_label, 0, 2)

        g.addWidget(QLabel("Iter"), 1, 0)
        self.pid_progress = QProgressBar()
        self.pid_progress.setRange(0, 100)
        self.pid_progress.setFixedHeight(UI.px(14))
        g.addWidget(self.pid_progress, 1, 1)
        self.pid_label = QLabel("--")
        self.pid_label.setFixedWidth(UI.px(45))
        self.pid_label.setObjectName("MutedLabel")
        g.addWidget(self.pid_label, 1, 2)

        return box

    # ── Live-Werte (links, kompakt) ──────────────────────────

    def _build_metrics(self):
        box = QGroupBox("Live")
        g = QGridLayout(box)
        g.setContentsMargins(3, 10, 3, 3)
        g.setSpacing(2)

        defs = [
            ("Q0", "Q0sccm", "sccm"), ("P",  "P_RFG", "W"),    ("P_set", "P_set", "W"),
            ("I",  "I_beam", "mA"),    ("Te", "Te",    "eV"),   ("Tg",    "Tg",   "K"),
            ("n",  "n",      "m-3"),   ("ng", "ng",    "m-3"),  ("ion%",  "iondeg", "%"),
        ]
        # I-fix + RF-Diagnostik Metriken
        ifix_defs = [
            ("I_soll", "I_target", "mA"), ("dI", "delta_I", "mA"), ("Stat", "ifix_status", ""),
            ("P_abs", "P_abs", "W"),      ("zeta", "zeta", ""),     ("", "_spacer", ""),
        ]
        all_defs = defs + ifix_defs
        for i, (title, key, unit) in enumerate(all_defs):
            card = MetricCard(title, unit)
            g.addWidget(card, i // 3, i % 3)
            self.metric_cards[key] = card
        return box

    # ── Styles ───────────────────────────────────────────────

    def _apply_styles(self):
        self.setStyleSheet(UI.qss())

    # ── Header sync ──────────────────────────────────────────

    def _sync_hdr(self, *_):
        pkg = self._current_package() if hasattr(self, "_packages") else None
        pkg_label = pkg.display_name if pkg else "?"
        r = self.cmb_rate_model.currentText() if hasattr(self, "cmb_rate_model") else "?"
        m = "SC" if self.solve_mode() == 2 else "I-fix"
        if hasattr(self, "_hdr_pkg"):
            self._hdr_pkg.setText(pkg_label)
            self._hdr_rate.setText(r)
            self._hdr_mode.setText(m)

    # ═════════════════════════════════════════════════════════
    # Backend (unveraendert)
    # ═════════════════════════════════════════════════════════

    def solve_mode(self):
        return self.mode_group.checkedId() if hasattr(self, "mode_group") else 1

    def log_msg(self, msg):
        self.log.append(msg)
        c = self.log.textCursor()
        c.movePosition(QTextCursor.MoveOperation.End)
        self.log.setTextCursor(c)

    def set_status(self, msg):
        self.statusBar().showMessage(msg)

    def update_metric(self, key, value):
        c = self.metric_cards.get(key)
        if c:
            c.set_value(value)

    def update_solution_power(self, v):
        self.update_metric("P_RFG", "--" if v is None else f"{v:.2f}")

    def update_set_power(self, v):
        self.update_metric("P_set", "--" if v is None else f"{v:.2f}")

    def solver_method(self):
        return 4

    def solver_label(self):
        return "Newton (stationaer)"

    def wsl_available(self):
        return sys.platform == "win32" and shutil.which("wsl") is not None

    def win_to_wsl(self, path):
        path = path.replace("\\", "/")
        if len(path) >= 2 and path[1] == ":":
            return f"/mnt/{path[0].lower()}{path[2:]}"
        return path

    def read_float(self, le, name):
        try:
            return float(le.text().strip())
        except ValueError as e:
            raise ValueError(f"{name}: ungueltig") from e

    def validate_inputs(self):
        errs = []
        vals = {}
        ae = {**self.entries, **self.sweep_entries}
        for key, fn, msg in VALIDATION_RULES:
            if key not in ae:
                continue
            try:
                v = float(ae[key].text().strip())
                vals[key] = v
                if not fn(v):
                    errs.append(f"{msg} (={v})")
            except ValueError:
                errs.append(f"{key}: ungueltig")
        errs.extend(cross_validate(vals))
        if errs:
            QMessageBox.warning(self, "Fehler", "\n".join(f"- {e}" for e in errs))
            return False
        return True

    def _build_run_config(self):
        """Erzeuge RunConfig aus aktuellem GUI-Zustand."""
        from run_config import RunConfig, GeometryConfig, GridConfig, CoilConfig, OperationConfig, SweepConfig, RateConfig, MetaConfig
        pkg = self._current_package()
        preset_id = self.cmb_thruster.currentData() if hasattr(self, "cmb_thruster") else "custom"
        rf = self.read_float

        cfg = RunConfig()
        cfg.geometry = GeometryConfig(
            R=rf(self.entries["R"], "R"),
            L=rf(self.entries["L"], "L"),
            betai=rf(self.entries["betai"], "betai"),
            betag=rf(self.entries["betag"], "betag"))
        cfg.grid = GridConfig(
            Vgrid=rf(self.entries["Vgrid"], "Vgrid"),
            sgrid=rf(self.entries["sgrid"], "sgrid"),
            eta_opt=getattr(self, "_eta_opt", 1.0))
        cfg.coil = CoilConfig(
            frequency=rf(self.entries["frequency"], "frequency"),
            Nw=rf(self.entries["Nw"], "Nw"),
            R_ohm=rf(self.entries["R_ohm"], "R_ohm"),
            Rc=rf(self.entries["Rc"], "Rc"),
            lc=rf(self.entries["lc"], "lc"))
        cfg.operation = OperationConfig(
            solve_mode=self.solve_mode(),
            P_max=rf(self.entries["P_RFG_max"], "P_RFG_max"),
            I_soll=rf(self.entries["I_soll"], "I_soll"),
            density_profile_factor=rf(self.entries["density_profile_factor"], "density_profile_factor"),
            rf_coupling=self._chk_rf_coupling.isChecked() if hasattr(self, "_chk_rf_coupling") else False)
        # Sweep: Intern immer in sccm. Bei mg/s-Modus zurueckrechnen.
        q0_start_raw = rf(self.sweep_entries["Q0sccm_start"], "Q0sccm_start")
        q0_step_raw = rf(self.sweep_entries["Q0sccm_step"], "Q0sccm_step")
        if getattr(self, "_q0_unit_mode", "sccm") == "mg/s":
            try:
                from flow_units import mg_per_s_to_sccm
                gas = pkg.gas if pkg else "xenon"
                q0_start_raw = mg_per_s_to_sccm(q0_start_raw, gas)
                q0_step_raw = mg_per_s_to_sccm(q0_step_raw, gas)
            except ImportError:
                pass
        cfg.sweep = SweepConfig(
            Q0_start=q0_start_raw,
            Q0_step=q0_step_raw,
            N=int(rf(self.sweep_entries["jjmax"], "jjmax")))
        cfg.rates = RateConfig(
            rate_model=self.cmb_rate_model.currentData() if hasattr(self, "cmb_rate_model") else 0)
        cfg.meta = MetaConfig(
            backend="cpp" if (pkg and pkg.is_cpp) else "python",
            package_id=pkg.id if pkg else "",
            gas=pkg.gas if pkg else "xenon",
            cs_database=pkg.database if (pkg and pkg.database) else "biagi",
            preset_id=preset_id or "custom")
        return cfg

    def write_config(self):
        if not self.validate_inputs():
            return False
        pkg = self._current_package()
        if not pkg:
            QMessageBox.critical(self, "Fehler", "Kein Chemiepaket gewaehlt.")
            return False
        try:
            cfg = self._build_run_config()

            # Validierung via RunConfig
            config_errs = cfg.validate()
            if config_errs:
                QMessageBox.warning(self, "Konfigurationsfehler",
                    "\n".join(f"- {e}" for e in config_errs))
                return False

            # Primaer: run_config.json (beide Backends)
            cfg.save_json(SCRIPT_DIR / "run_config.json")
            # Kompatibilitaet: params.txt (C++-Legacy-Pfad)
            cfg.to_params_txt(CONFIG_FILE)

            backend = "C++" if pkg.is_cpp else "Python"
            self.log_msg(f"Config: {pkg.display_name} [{backend}], "
                         f"{'SC' if self.solve_mode() == 2 else 'I-fix'}, "
                         f"{self.cmb_rate_model.currentText()}, "
                         f"Preset={cfg.meta.preset_id}")
            return True
        except (ValueError, OSError) as e:
            QMessageBox.critical(self, "Fehler", str(e))
            return False

    def compile_cpp(self):
        if hasattr(self, "_cp") and self._cp and self._cp.state() != QProcess.ProcessState.NotRunning:
            return
        self.set_status("Kompiliere ...")
        self.log_msg("Kompilierung ...")
        # Modulare Kompilierung: 6 Objekte + main -> chabert
        CPP_MODULES = ["bessel_wrapper", "sim_config", "rates", "physics", "solver", "sim_logging"]
        if self.wsl_available():
            cwd = self.win_to_wsl(str(SCRIPT_DIR))
            compile_steps = " && ".join(
                f'g++ -O3 -march=native -std=c++17 -c {m}.cpp -o {m}.o 2>&1'
                for m in CPP_MODULES)
            link_objs = " ".join(f'{m}.o' for m in CPP_MODULES)
            prog, args = "wsl", ["bash", "-c",
                f'cd "{cwd}" && {compile_steps} && '
                f'g++ -O3 -march=native -std=c++17 {CPP_SOURCE} {link_objs} -o chabert 2>&1']
        else:
            cc = shutil.which("g++")
            if not cc:
                QMessageBox.critical(self, "Fehler", "g++ nicht gefunden.")
                return
            out = "chabert.exe" if sys.platform == "win32" else "chabert"
            ext = ".obj" if sys.platform == "win32" else ".o"
            shell = "cmd" if sys.platform == "win32" else "bash"
            flag = "/c" if sys.platform == "win32" else "-c"
            compile_steps = " && ".join(
                f'"{cc}" -O3 -march=native -std=c++17 -c {m}.cpp -o {m}{ext}'
                for m in CPP_MODULES)
            link_objs = " ".join(f'{m}{ext}' for m in CPP_MODULES)
            prog, args = shell, [flag,
                f'{compile_steps} && '
                f'"{cc}" -O3 -march=native -std=c++17 {CPP_SOURCE} {link_objs} -o {out}']
        self._cp = QProcess(self)
        self._cp.setWorkingDirectory(str(SCRIPT_DIR))
        self._cp.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self._cp.readyReadStandardOutput.connect(
            lambda: self.log_msg("  " + bytes(self._cp.readAllStandardOutput()).decode("utf-8", errors="replace").strip()) if self._cp else None)
        self._cp.finished.connect(lambda c, _: (
            self.log_msg("OK." if c == 0 else f"Fehler ({c})."),
            self.set_status("Kompiliert." if c == 0 else "Fehler.")))
        self._cp.start(prog, args)

    def _open_physics_data(self):
        """Open physics data viewer for the currently selected package."""
        pkg = self._current_package()
        if not pkg:
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Fehler", "Kein Paket ausgewaehlt.")
            return
        try:
            from physics_data_viewer import PhysicsDataWindow
            pkg_info = {
                "display_name": pkg.display_name,
                "gas": pkg.gas,
                "backend": "cpp" if pkg.is_cpp else "python",
                "database": pkg.database or "",
                "path": str(pkg.path),
                "id": pkg.id,
            }
            self._phys_win = PhysicsDataWindow(pkg_info, parent=self)
            self._phys_win.show()
        except Exception as e:
            self.log_msg(f"Physik-Daten: {e}")

    def _ifix_preflight_check(self):
        """Pre-Flight-Pruefung fuer I-fix: Warnung bei unplausiblem Sollstrom.

        Prueft nur gegen Preset-Metadaten (kein Solver-Aufruf).
        Gibt True zurueck wenn Lauf starten soll, False bei Abbruch.
        """
        try:
            I_soll = float(self.entries["I_soll"].text())
        except (ValueError, KeyError):
            return True  # Ungueltige Eingabe wird von validate_inputs gefangen

        preset = self._get_current_preset()
        if not preset:
            return True  # Custom oder kein Preset -> keine Pruefung

        plaus = preset.get("plausibility", {})
        I_max = plaus.get("I_beam_typ_max_mA")
        I_min = plaus.get("I_beam_typ_min_mA")

        if I_max is not None and I_soll > I_max * 1.5:
            # Deutlich ausserhalb -> Bestaetigungsdialog
            ans = QMessageBox.warning(
                self, "I_soll ausserhalb Preset-Bereich",
                f"Der Sollstrom I_soll = {I_soll:.0f} mA liegt deutlich oberhalb "
                f"des typischen Bereichs fuer '{preset['name']}' ({I_min}-{I_max} mA).\n\n"
                f"Der Sweep wird voraussichtlich keine Konvergenz erreichen.\n"
                f"Trotzdem starten?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No)
            if ans != QMessageBox.StandardButton.Yes:
                self.log_msg(f"Abgebrochen: I_soll={I_soll:.0f} mA ausserhalb Preset-Bereich ({I_min}-{I_max} mA)")
                return False
            self.log_msg(f"Warnung: I_soll={I_soll:.0f} mA oberhalb Preset-Bereich ({I_min}-{I_max} mA)")

        # Gas-Kompatibilitaet pruefen
        pkg = self._current_package()
        if pkg:
            rec = preset.get("recommended_gas", [])
            gas_warnings = plaus.get("gas_warnings", {})
            if pkg.gas in gas_warnings and gas_warnings[pkg.gas]:
                self.log_msg(f"Hinweis: {gas_warnings[pkg.gas]}")
            elif rec and pkg.gas not in rec:
                self.log_msg(f"Hinweis: Preset '{preset['name']}' ist fuer "
                             f"{', '.join(r.capitalize() for r in rec)} empfohlen, "
                             f"Paket ist {pkg.gas.capitalize()}.")

        return True

    def run_sim(self):
        if self.process and self.process.state() != QProcess.ProcessState.NotRunning:
            return
        pkg = self._current_package()
        if not pkg:
            QMessageBox.warning(self, "Fehler", "Kein Chemiepaket gewaehlt.")
            return

        # Backend-spezifische Pruefung
        if pkg.is_cpp:
            bn = "chabert.exe" if sys.platform == "win32" else "chabert"
            if self.wsl_available():
                cwd = self.win_to_wsl(str(SCRIPT_DIR))
                if subprocess.run(["wsl", "bash", "-c", f'test -f "{cwd}/chabert"']).returncode != 0:
                    QMessageBox.warning(self, "Fehler", "C++ Binary nicht gefunden. Bitte kompilieren.")
                    return
            elif not (SCRIPT_DIR / bn).exists():
                QMessageBox.warning(self, "Fehler", f"'{bn}' nicht gefunden. Bitte kompilieren.")
                return
        elif pkg.is_python:
            if not (SCRIPT_DIR / "generic_solver.py").exists():
                QMessageBox.warning(self, "Fehler", "generic_solver.py nicht gefunden.")
                return

        if not self.write_config():
            return

        # I-fix Pre-Flight: Warnung bei unplausiblem I_soll
        if self.solve_mode() == 1 and not self._ifix_preflight_check():
            return

        self.cancel_requested = False
        self._stdout_buffer = ""
        self.current_q0 = None
        self._ifix_target_set = False
        self._last_delta_I = 0.0
        self._last_P_abs_rf = 0.0
        self.scan_progress.setValue(0)
        self.pid_progress.setValue(0)
        self.scan_label.setText("--")
        self.pid_label.setText("--")
        self.clear_plots()

        self.process = QProcess(self)
        self.process.setWorkingDirectory(str(SCRIPT_DIR))
        self.process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.process.readyReadStandardOutput.connect(self._on_out)
        self.process.finished.connect(self._on_done)

        backend_label = "C++" if pkg.is_cpp else "Python"
        self.log_msg(f"\nSim [{backend_label}]: {pkg.display_name}")
        self.set_status(f"Simulation ({backend_label}) ...")

        if pkg.is_cpp:
            # C++-Backend: chabert mit run_config.json (primaer) oder params.txt (Fallback)
            cfg_file = "run_config.json" if (SCRIPT_DIR / "run_config.json").exists() else CONFIG_FILE
            if self.wsl_available():
                cwd = self.win_to_wsl(str(SCRIPT_DIR))
                self.process.start("wsl", ["bash", "-c", f'cd "{cwd}" && ./chabert {cfg_file} 2>&1'])
            else:
                self.process.start(str(SCRIPT_DIR / bn), [cfg_file])
        else:
            # Python-Backend: generic_solver.py mit Sweep-Parametern
            # Bei mg/s-Modus intern nach sccm konvertieren fuer den Solver
            sm = self.solve_mode()  # 1=I-fix, 2=SC
            Q0_start = float(self.sweep_entries["Q0sccm_start"].text())
            Q0_step = float(self.sweep_entries["Q0sccm_step"].text())
            if getattr(self, "_q0_unit_mode", "sccm") == "mg/s":
                try:
                    from flow_units import mg_per_s_to_sccm
                    gas = pkg.gas if pkg else "xenon"
                    Q0_start = mg_per_s_to_sccm(Q0_start, gas)
                    Q0_step = mg_per_s_to_sccm(Q0_step, gas)
                except ImportError:
                    pass
            N_steps = int(float(self.sweep_entries["jjmax"].text()))
            if sm == 1:
                # I-fix: param = I_soll [mA]
                param = float(self.entries["I_soll"].text())
            else:
                # SC: param = P_RFG [W]
                param = float(self.entries["P_RFG"].text())
            self.process.start(sys.executable, [
                str(SCRIPT_DIR / "generic_solver.py"),
                pkg.chemistry_json,
                str(sm),
                str(param),
                str(Q0_start),
                str(Q0_step),
                str(N_steps),
            ])

        if not self.process.waitForStarted(3000):
            self.set_status("Startfehler.")

    def _on_out(self):
        if not self.process:
            return
        self._stdout_buffer += bytes(self.process.readAllStandardOutput()).decode("utf-8", errors="replace")
        while "\n" in self._stdout_buffer:
            line, self._stdout_buffer = self._stdout_buffer.split("\n", 1)
            self._parse(line.strip())

    def _parse(self, line):
        if not line:
            return
        p = line.split()
        tag = p[0]
        try:
            if tag == "Q0_STEP" and len(p) >= 4:
                q0_sccm = float(p[1])
                self.current_q0 = q0_sccm
                self.scan_progress.setValue(int(int(p[2]) / int(p[3]) * 100))
                self.scan_label.setText(f"{p[2]}/{p[3]}")
                # Q0-Anzeige in aktiver Einheit
                if getattr(self, "_q0_unit_mode", "sccm") == "mg/s":
                    try:
                        from flow_units import sccm_to_mg_per_s
                        q0_display = sccm_to_mg_per_s(q0_sccm, self._get_gas())
                        self.update_metric("Q0sccm", f"{q0_display:.4f}")
                    except ImportError:
                        self.update_metric("Q0sccm", f"{q0_sccm:.4f}")
                else:
                    self.update_metric("Q0sccm", f"{q0_sccm:.4f}")
                self.pid_progress.setValue(0)
                self.pid_label.setText("--")
                self.log_msg(f"Q0={p[1]} sccm ({p[2]}/{p[3]})")
                return
            if tag == "PID_START" and len(p) >= 3:
                self.pid_progress.setValue(int(min(int(p[1]) / 50 * 100, 100)))
                self.pid_label.setText(f"It {p[1]}")
                self.update_set_power(float(p[2]))
                return
            if tag == "PID_DONE" and len(p) >= 6:
                self.update_metric("I_beam", f"{float(p[1]):.3f}")
                self.update_metric("Te", f"{float(p[4]):.2f}")
                self.update_metric("Tg", f"{float(p[5]):.1f}")
                self.update_set_power(float(p[3]))
                return
            if tag == "IFIX_RESULT" and len(p) >= 7:
                # IFIX_RESULT Q0 P I_target I_beam delta_I status
                q0_if = float(p[1])
                P_if = float(p[2])
                I_target = float(p[3])
                I_beam_if = float(p[4])
                delta_I = float(p[5])
                ifix_status = p[6]

                self.update_metric("I_target", f"{I_target:.1f}")
                self.update_metric("delta_I", f"{delta_I:+.2f}")
                self.update_metric("ifix_status", ifix_status)

                # I_target-Referenzlinie beim ersten Punkt setzen
                if not hasattr(self, "_ifix_target_set") or not self._ifix_target_set:
                    self.plot_grid.set_I_target(I_target)
                    self._ifix_target_set = True

                # Fehlpunkte als rote Marker (kein normaler Plot-Punkt)
                if ifix_status != "converged":
                    q0_fail = q0_if
                    if getattr(self, "_q0_unit_mode", "sccm") == "mg/s":
                        try:
                            from flow_units import sccm_to_mg_per_s
                            q0_fail = sccm_to_mg_per_s(q0_if, self._get_gas())
                        except ImportError:
                            pass
                    self.plot_grid.add_fail_point(q0_fail, P_if, I_beam_if)
                    self.log_msg(f"  I-fix: Q0={q0_if:.3f} P={P_if:.1f}W "
                                 f"I={I_beam_if:.1f}/{I_target:.0f}mA "
                                 f"dI={delta_I:+.1f} [{ifix_status}]")
                else:
                    self.log_msg(f"  I-fix: Q0={q0_if:.3f} P={P_if:.1f}W "
                                 f"I={I_beam_if:.2f}/{I_target:.0f}mA dI={delta_I:+.2f}")

                # delta_I fuer nachfolgenden RESULT-Punkt merken
                self._last_delta_I = delta_I
                return
            if tag == "CONVERGED":
                self.pid_progress.setValue(100)
                self.pid_label.setText("OK")
                return
            if tag == "RESULT" and len(p) >= 7:
                n, ng, te, tg, i, pf = (float(x) for x in p[1:7])
                self.update_metric("n", f"{n:.2e}")
                self.update_metric("ng", f"{ng:.2e}")
                self.update_metric("Te", f"{te:.2f}")
                self.update_metric("Tg", f"{tg:.1f}")
                self.update_metric("I_beam", f"{i:.3f}")
                if ng > 0:
                    self.update_metric("iondeg", f"{n / ng * 100:.2f}")
                self.update_solution_power(pf)
                self.update_set_power(pf)
                if self.current_q0 is not None:
                    # Iod-Daten werden im naechsten IODINE_EXT nachgeliefert
                    # Q0 fuer Plot in aktiver Einheit
                    q0_plot = self.current_q0
                    if getattr(self, "_q0_unit_mode", "sccm") == "mg/s":
                        try:
                            from flow_units import sccm_to_mg_per_s
                            q0_plot = sccm_to_mg_per_s(self.current_q0, self._get_gas())
                        except ImportError:
                            pass
                    pabs_rf = getattr(self, "_last_P_abs_rf", 0.0)
                    self.plot_grid.add_point(q0_plot, pf, te, tg, n, ng, i, P_abs_rf=pabs_rf)
                    self._last_P_abs_rf = 0.0
                    self._last_delta_I = 0.0
                return
            if tag == "IODINE_EXT" and len(p) >= 10:
                # nI nI2 nIp nI2p nIm diss fIp fI2p alpha
                nI, nI2, nIp, nI2p, nIm = (float(x) for x in p[1:6])
                diss, fIp, fI2p, alpha = (float(x) for x in p[6:10])
                self.update_metric("iondeg", f"{fIp*100:.1f}")
                iod = {"nI": nI, "nI2": nI2, "nIp": nIp, "nI2p": nI2p,
                       "nIm": nIm, "diss": diss, "fIp": fIp, "fI2p": fI2p, "alpha": alpha}
                self._last_iodine = iod
                # Iod-Daten direkt in Plot-Daten einfuegen (letzter Punkt)
                self.plot_grid.update_last_iodine(iod)
                self.log_msg(f"  f_I+={fIp:.3f} diss={diss:.3f} alpha={alpha:.4f}")
                return
            if tag == "PID_MAXITER":
                self.log_msg("  PID max iter")
                return
            if tag == "P_LIMIT_REACHED" and len(p) >= 6:
                self.log_msg(f"  RF-Limit Q0={float(p[1]):.3f}")
                self.update_set_power(float(p[3]))
                return
            if tag in ("SWEEP_RECOVERY", "SWEEP_RECOVERY_FAIL"):
                self.log_msg(f"  {tag}: {' '.join(p[1:])}")
                return
            if tag == "RF_COUPLING":
                self.log_msg(f"  RF-Kopplung: {' '.join(p[1:])}")
                return
            if tag == "RF_DIAG" and len(p) >= 4:
                # RF_DIAG P_RFG=400.0 P_abs=278.8 zeta=0.6971
                for part in p[1:]:
                    if part.startswith("P_abs="):
                        pabs_val = part.split("=")[1]
                        self.update_metric("P_abs", pabs_val)
                        self._last_P_abs_rf = float(pabs_val)
                    elif part.startswith("zeta="):
                        self.update_metric("zeta", part.split("=")[1])
                self.log_msg(f"  RF: {' '.join(p[1:])}")
                return
            if tag == "SOLVER_FAIL":
                self.update_solution_power(None)
                # delta_I zuruecksetzen (verhindert Bleed zum naechsten Punkt)
                self._last_delta_I = 0.0
                # Menschenlesbare Diagnose extrahieren
                diag_idx = line.find("diag=")
                if diag_idx >= 0:
                    diag_text = line[diag_idx + 5:]
                    self.log_msg(f"  ⚠ {diag_text}")
                else:
                    self.log_msg(f"  FAIL: {line}")
                return
            if not line.startswith("STEP "):
                self.log_msg(line)
        except Exception as e:
            self.log_msg(f"Parse: {line} | {e}")

    def _on_done(self, code, _):
        if self.cancel_requested:
            self.log_msg("\nAbgebrochen.")
            self.set_status("Abgebrochen.")
            self.cancel_requested = False
            return
        if code == 0:
            self.log_msg("\nFertig.")
            self.set_status("Fertig.")
            self.scan_progress.setValue(100)
        else:
            self.log_msg(f"\nFehler ({code}).")
            self.set_status(f"Fehler ({code}).")

    def cancel_sim(self):
        if self.process and self.process.state() != QProcess.ProcessState.NotRunning:
            self.cancel_requested = True
            self.process.kill()
            self.set_status("Abbruch ...")

    def reset_defaults(self):
        for k, (_, _, d) in ALL_PRIMARY.items():
            if k in self.entries:
                self.entries[k].setText(str(d))
        for k, (_, _, d) in ALL_SWEEP.items():
            if k in self.sweep_entries:
                self.sweep_entries[k].setText(str(d))
        self.mode_buttons[DEFAULT_SOLVE_MODE].setChecked(True)
        self._populate_packages()  # Reset to default package
        self.cmb_rate_model.setCurrentIndex(DEFAULT_RATE_MODEL)
        self.log_msg("Defaults.")

    def clear_plots(self):
        self.plot_grid.clear()

    def open_output(self):
        p = SCRIPT_DIR / OUTPUT_FILE
        if p.exists():
            QDesktopServices.openUrl(p.as_uri())

    def open_log_viewer(self):
        vs = SCRIPT_DIR / "log_viewer.py"
        if not vs.exists():
            return
        logs = sorted(SCRIPT_DIR.glob("simulation_log_*.txt"), key=lambda x: x.stat().st_mtime, reverse=True)
        args = [sys.executable, str(vs)] + ([str(logs[0])] if logs else [])
        subprocess.Popen(args, cwd=str(SCRIPT_DIR))
        self.log_msg("Log-Viewer" + (f" ({logs[0].name})" if logs else ""))

    def closeEvent(self, event):
        try:
            if self.process and self.process.state() != QProcess.ProcessState.NotRunning:
                self.cancel_requested = True
                self.process.kill()
                self.process.waitForFinished(1500)
        except Exception:
            pass
        super().closeEvent(event)


def main():
    app = QApplication(sys.argv)
    UI.init()
    w = SimulatorWindow()
    w.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
