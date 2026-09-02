"""
physics_data_viewer.py -- Physics Data Inspection Window.

Shows cross sections sigma(E) and rate coefficients K(Te) for any
established model/package. Supports multi-curve selection, axis scaling
control, tooltips, and a process details panel.

Features:
  - Separate plots for cross sections and rate coefficients
  - Multi-curve selection with distinct colors and legend
  - Manual axis control: log/linear for both x and y on each plot
  - Tooltips on process list entries with compact metadata
  - Details panel showing rich metadata for selected process
  - Handles missing data gracefully per package type
"""
from __future__ import annotations

import json
import math
import numpy as np
from pathlib import Path
from typing import Optional

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QSplitter,
    QGroupBox, QLabel, QListWidget, QListWidgetItem,
    QAbstractItemView, QCheckBox, QTextEdit, QStatusBar,
)
import pyqtgraph as pg

SCRIPT_DIR = Path(__file__).resolve().parent

COLORS = [
    "#5cb8ff", "#ff6b6b", "#2dd4bf", "#fbbf24", "#a78bfa",
    "#34d399", "#f472b6", "#60a5fa", "#fb923c", "#c084fc",
    "#ff9f43", "#26de81", "#778ca3", "#fc5c65", "#45aaf2",
]

# ═══════════════════════════════════════════════════════════════
# Data Loading (unchanged logic, richer metadata extraction)
# ═══════════════════════════════════════════════════════════════

def _read_csv_2col(path: Path) -> Optional[tuple[list[float], list[float], dict]]:
    """Read a 2-column CSV. Returns (x, y, metadata) or None."""
    if not path.exists():
        return None
    meta = {"file": path.name}
    x, y = [], []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                for key in ("Process:", "Threshold:", "Source:", "Type:",
                            "Database:", "Species:", "Comment:", "Unit:"):
                    if key in s:
                        meta[key.rstrip(":")] = s.split(key, 1)[1].strip()
                continue
            parts = s.split(",")
            if len(parts) >= 2:
                try:
                    x.append(float(parts[0].strip()))
                    y.append(float(parts[1].strip()))
                except ValueError:
                    continue
    if not x:
        return None
    meta["n_points"] = len(x)
    meta["x_min"] = min(x)
    meta["x_max"] = max(x)
    return x, y, meta


def discover_cross_sections(pkg_path: Path) -> list[dict]:
    """Find all cross section CSV files for a C++ package."""
    items = []
    for name, label, proc_type in [
        ("elastic.csv", "Elastic (momentum transfer)", "elastic"),
        ("ionization.csv", "Ionization", "ionization"),
    ]:
        p = pkg_path / name
        if p.exists():
            result = _read_csv_2col(p)
            meta = result[2] if result else {}
            items.append({
                "path": str(p), "label": label, "process": proc_type,
                "data_type": "cross_section", "unit_x": "eV", "unit_y": "m^2",
                "threshold_eV": meta.get("Threshold", ""),
                "source": meta.get("Source", meta.get("Database", "")),
                "n_points": meta.get("n_points", 0),
                "E_range": f"{meta.get('x_min', '?')}-{meta.get('x_max', '?')} eV",
            })

    exc_dir = pkg_path / "excitation"
    if exc_dir.is_dir():
        for f in sorted(exc_dir.glob("excitation_*.csv")):
            state = f.stem.replace("excitation_", "").replace("_", " ").upper()
            result = _read_csv_2col(f)
            meta = result[2] if result else {}
            items.append({
                "path": str(f), "label": f"Excitation {state}",
                "process": "excitation", "state": state,
                "data_type": "cross_section", "unit_x": "eV", "unit_y": "m^2",
                "threshold_eV": meta.get("Threshold", ""),
                "source": meta.get("Source", ""),
                "n_points": meta.get("n_points", 0),
                "E_range": f"{meta.get('x_min', '?')}-{meta.get('x_max', '?')} eV",
            })
    return items


def discover_rate_coefficients(pkg_path: Path, pkg_info: dict) -> list[dict]:
    """Find all rate coefficient data for a package."""
    items = []
    backend = pkg_info.get("backend", "")

    if backend == "cpp":
        for name, label, proc in [
            ("kel_table.csv", "K_el (elastic)", "elastic"),
            ("kiz_table.csv", "K_iz (ionization)", "ionization"),
            ("kex_table.csv", "K_ex (excitation, total)", "excitation"),
        ]:
            p = pkg_path / name
            if p.exists():
                result = _read_csv_2col(p)
                meta = result[2] if result else {}
                items.append({
                    "path": str(p), "label": label, "process": proc,
                    "data_type": "rate_coefficient", "source": "Maxwell-Boltzmann integration",
                    "representation": "tabulated", "unit_x": "eV", "unit_y": "m^3/s",
                    "n_points": meta.get("n_points", 0),
                    "Te_range": f"{meta.get('x_min', '?')}-{meta.get('x_max', '?')} eV",
                })

    elif backend == "python":
        rates_dir = pkg_path / "rates"
        # Try to load reaction metadata from chemistry.json
        rxn_meta = {}
        chem_json = pkg_path / "chemistry.json"
        if chem_json.exists():
            try:
                with open(chem_json, encoding="utf-8") as fj:
                    chem = json.load(fj)
                for rxn in chem.get("reactions", []):
                    rid = rxn.get("id", "")
                    rxn_meta[rid] = rxn
            except (json.JSONDecodeError, OSError):
                pass

        if rates_dir.is_dir():
            for f in sorted(rates_dir.glob("*.csv")):
                stem = f.stem
                label = stem.replace("_", " ")
                result = _read_csv_2col(f)
                meta = result[2] if result else {}
                # Match to reaction metadata
                rxn = {}
                for rid, rm in rxn_meta.items():
                    rate_file = rm.get("rate", {}).get("file", "")
                    if stem in rate_file or rid in stem.lower():
                        rxn = rm
                        break
                reaction_str = ""
                if rxn:
                    reactants = _format_reaction_side(rxn.get("reactants", {}))
                    products = _format_reaction_side(rxn.get("products", {}))
                    if reactants and products:
                        reaction_str = f"{reactants} -> {products}"

                items.append({
                    "path": str(f), "label": label,
                    "process": rxn.get("type", stem),
                    "data_type": "rate_coefficient",
                    "representation": "tabulated",
                    "reaction": reaction_str,
                    "energy_eV": rxn.get("energy_eV", ""),
                    "source": rxn.get("source", meta.get("Source", "")),
                    "unit_x": "eV", "unit_y": "m^3/s",
                    "n_points": meta.get("n_points", 0),
                    "Te_range": f"{meta.get('x_min', '?')}-{meta.get('x_max', '?')} eV",
                })

        # Arrhenius-only fallback
        if not items and chem_json.exists():
            for rid, rxn in rxn_meta.items():
                rate = rxn.get("rate", {})
                model = rate.get("model", "")
                if model in ("arrhenius", "constant"):
                    A = rate.get("A", rate.get("value", 0))
                    E_a = rate.get("E_a_eV", 0)
                    reactants = _format_reaction_side(rxn.get("reactants", {}))
                    products = _format_reaction_side(rxn.get("products", {}))
                    items.append({
                        "label": f"{rxn.get('name', rid)} ({model})",
                        "process": rxn.get("type", "?"),
                        "data_type": "rate_coefficient",
                        "representation": model,
                        "reaction": f"{reactants} -> {products}" if reactants else "",
                        "energy_eV": rxn.get("energy_eV", ""),
                        "A": A, "E_a_eV": E_a, "model": model,
                        "source": "Arrhenius fit",
                        "unit_x": "eV", "unit_y": "m^3/s",
                    })
    return items


def _format_reaction_side(species_dict: dict) -> str:
    """Format a reactants or products dict with stoichiometry.

    Example: {"e": 2, "I+": 1, "I": 1} -> "2e + I+ + I"
    """
    parts = []
    for sp, n in species_dict.items():
        if n > 1:
            parts.append(f"{n}{sp}")
        elif n == 1:
            parts.append(sp)
    return " + ".join(parts)


def compute_arrhenius_curve(A: float, E_a_eV: float, model: str = "arrhenius"):
    Te = np.linspace(0.5, 20, 200)
    if model == "constant":
        K = np.full_like(Te, A)
    else:
        K = A * np.exp(-E_a_eV / np.maximum(Te, 0.01))
    return Te.tolist(), K.tolist()


def _build_tooltip(item: dict) -> str:
    """Build a compact tooltip string from item metadata."""
    parts = []
    if item.get("reaction"):
        parts.append(item["reaction"])
    if item.get("process"):
        parts.append(f"Type: {item['process']}")
    if item.get("threshold_eV"):
        parts.append(f"Threshold: {item['threshold_eV']} eV")
    if item.get("energy_eV"):
        parts.append(f"Energy loss: {item['energy_eV']} eV")
    if item.get("representation"):
        parts.append(f"Repr: {item['representation']}")
    if item.get("n_points"):
        parts.append(f"Points: {item['n_points']}")
    r = item.get("E_range") or item.get("Te_range")
    if r and "?" not in str(r):
        parts.append(f"Range: {r}")
    if item.get("source"):
        parts.append(f"Source: {item['source']}")
    return "\n".join(parts) if parts else "(no metadata)"


def _build_details(item: dict) -> str:
    """Build a detailed multi-line info text."""
    lines = [f"<b>{item.get('label', '?')}</b>"]
    if item.get("reaction"):
        lines.append(f"Reaction: <code>{item['reaction']}</code>")
    lines.append(f"Data type: {item.get('data_type', '?')}")
    lines.append(f"Process: {item.get('process', '?')}")
    if item.get("representation"):
        lines.append(f"Representation: {item['representation']}")
    if item.get("threshold_eV"):
        lines.append(f"Threshold energy: {item['threshold_eV']} eV")
    if item.get("energy_eV"):
        lines.append(f"Energy loss per event: {item['energy_eV']} eV")
    if item.get("state"):
        lines.append(f"Excited state: {item['state']}")
    if item.get("source"):
        lines.append(f"Source: {item['source']}")
    if item.get("n_points"):
        lines.append(f"Data points: {item['n_points']}")
    r = item.get("E_range") or item.get("Te_range")
    if r:
        lines.append(f"Range: {r}")
    if item.get("unit_x") and item.get("unit_y"):
        lines.append(f"Units: {item['unit_x']} vs {item['unit_y']}")
    if item.get("A"):
        lines.append(f"A = {item['A']:.4e}")
    if item.get("E_a_eV"):
        lines.append(f"E_a = {item['E_a_eV']} eV")
    if item.get("path"):
        lines.append(f"<small>File: {Path(item['path']).name}</small>")
    return "<br>".join(lines)


# ═══════════════════════════════════════════════════════════════
# Physics Data Viewer Window
# ═══════════════════════════════════════════════════════════════

class PhysicsDataWindow(QMainWindow):

    def __init__(self, pkg_info: dict, parent=None):
        super().__init__(parent)
        self.pkg = pkg_info
        self.pkg_path = Path(pkg_info.get("path", ""))
        self.setWindowTitle(f"Physics Data -- {pkg_info.get('display_name', '?')}")
        self.resize(1300, 800)

        self._cs_items = discover_cross_sections(self.pkg_path) if pkg_info.get("backend") == "cpp" else []
        self._rate_items = discover_rate_coefficients(self.pkg_path, pkg_info)

        self._build_ui()
        self._apply_style()

    def _build_ui(self):
        root = QWidget()
        self.setCentralWidget(root)
        outer = QVBoxLayout(root)
        outer.setContentsMargins(6, 6, 6, 6)

        hdr = QLabel(f"<b>{self.pkg.get('display_name', '?')}</b> -- "
                     f"Gas: {self.pkg.get('gas', '?')}, "
                     f"Backend: {self.pkg.get('backend', '?')}, "
                     f"Database: {self.pkg.get('database', '-')}")
        hdr.setStyleSheet("color: #2dd4bf; font-size: 11pt; padding: 4px;")
        outer.addWidget(hdr)

        main_split = QSplitter(Qt.Orientation.Horizontal)
        outer.addWidget(main_split, 1)

        # ── Left: Selection + Details ────────────────────────
        left = QWidget()
        ll = QVBoxLayout(left)
        ll.setContentsMargins(0, 0, 4, 0)

        # Cross sections list
        grp_cs = QGroupBox(f"Cross Sections ({len(self._cs_items)})")
        cl = QVBoxLayout(grp_cs)
        self._cs_list = QListWidget()
        self._cs_list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        if self._cs_items:
            for item in self._cs_items:
                wi = QListWidgetItem(item["label"])
                wi.setToolTip(_build_tooltip(item))
                self._cs_list.addItem(wi)
        else:
            wi = QListWidgetItem("(no cross section data)")
            self._cs_list.addItem(wi)
            self._cs_list.setEnabled(False)
        self._cs_list.itemSelectionChanged.connect(self._on_cs_selection)
        self._cs_list.currentRowChanged.connect(self._on_cs_focus)
        cl.addWidget(self._cs_list)
        ll.addWidget(grp_cs)

        # Rate coefficients list
        grp_rate = QGroupBox(f"Rate Coefficients ({len(self._rate_items)})")
        rl = QVBoxLayout(grp_rate)
        self._rate_list = QListWidget()
        self._rate_list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        if self._rate_items:
            for item in self._rate_items:
                src = item.get("representation", item.get("source", "?"))
                wi = QListWidgetItem(f"{item['label']}  [{src}]")
                wi.setToolTip(_build_tooltip(item))
                self._rate_list.addItem(wi)
        else:
            wi = QListWidgetItem("(no rate data)")
            self._rate_list.addItem(wi)
            self._rate_list.setEnabled(False)
        self._rate_list.itemSelectionChanged.connect(self._on_rate_selection)
        self._rate_list.currentRowChanged.connect(self._on_rate_focus)
        rl.addWidget(self._rate_list)
        ll.addWidget(grp_rate)

        # Details panel
        grp_det = QGroupBox("Process Details")
        dl = QVBoxLayout(grp_det)
        self._details = QTextEdit()
        self._details.setReadOnly(True)
        self._details.setMaximumHeight(140)
        self._details.setStyleSheet("background: #050812; border: 1px solid #24314f; "
                                     "border-radius: 6px; color: #8ff7c8; font-size: 9pt;")
        self._details.setHtml("<i>Select a process to see details</i>")
        dl.addWidget(self._details)
        ll.addWidget(grp_det)

        main_split.addWidget(left)

        # ── Right: Plots with axis controls ──────────────────
        right = QWidget()
        rv = QVBoxLayout(right)
        rv.setContentsMargins(0, 0, 0, 0)
        rv.setSpacing(4)

        pg.setConfigOptions(antialias=True, background="#0b1220", foreground="#d7e3ff")

        # Cross section plot + controls
        cs_ctrl = QHBoxLayout()
        cs_ctrl.addWidget(QLabel("Cross Sections:"))
        self._cs_log_x = QCheckBox("Log X")
        self._cs_log_y = QCheckBox("Log Y")
        self._cs_log_y.setChecked(True)  # Default: log y for cross sections
        self._cs_log_x.toggled.connect(self._update_cs_axes)
        self._cs_log_y.toggled.connect(self._update_cs_axes)
        cs_ctrl.addWidget(self._cs_log_x)
        cs_ctrl.addWidget(self._cs_log_y)
        cs_ctrl.addStretch()
        rv.addLayout(cs_ctrl)

        self._cs_plot = pg.PlotWidget()
        self._cs_plot.setLabel("bottom", "Energy [eV]")
        self._cs_plot.setLabel("left", "Cross Section [m^2]")
        self._cs_plot.showGrid(x=True, y=True, alpha=0.2)
        self._cs_plot.setLogMode(x=False, y=True)
        rv.addWidget(self._cs_plot, 1)

        # Rate coefficient plot + controls
        rate_ctrl = QHBoxLayout()
        rate_ctrl.addWidget(QLabel("Rate Coefficients:"))
        self._rate_log_x = QCheckBox("Log X")
        self._rate_log_y = QCheckBox("Log Y")
        self._rate_log_y.setChecked(True)  # Default: log y for rates
        self._rate_log_x.toggled.connect(self._update_rate_axes)
        self._rate_log_y.toggled.connect(self._update_rate_axes)
        rate_ctrl.addWidget(self._rate_log_x)
        rate_ctrl.addWidget(self._rate_log_y)
        rate_ctrl.addStretch()
        rv.addLayout(rate_ctrl)

        self._rate_plot = pg.PlotWidget()
        self._rate_plot.setLabel("bottom", "Electron Temperature [eV]")
        self._rate_plot.setLabel("left", "Rate Coefficient [m^3/s]")
        self._rate_plot.showGrid(x=True, y=True, alpha=0.2)
        self._rate_plot.setLogMode(x=False, y=True)
        rv.addWidget(self._rate_plot, 1)

        main_split.addWidget(right)
        main_split.setSizes([320, 980])
        self.setStatusBar(QStatusBar())

    def _apply_style(self):
        self.setStyleSheet("""
            QMainWindow, QWidget { background: #0b1020; color: #e6edf7; font-size: 10pt; }
            QGroupBox { border: 1px solid #24314f; border-radius: 8px; margin-top: 8px;
                        padding-top: 10px; background: #0f172a; font-weight: 600; }
            QGroupBox::title { subcontrol-origin: margin; left: 8px; padding: 0 4px; color: #d9e4ff; }
            QListWidget { background: #111a31; border: 1px solid #24314f; border-radius: 6px;
                          color: white; font-size: 9pt; }
            QListWidget::item:selected { background: #2553d8; }
            QLabel { color: #9fb4dd; }
            QCheckBox { color: #e6edf7; spacing: 3px; }
        """)

    # ── Axis Control ─────────────────────────────────────────

    def _update_cs_axes(self):
        self._cs_plot.setLogMode(x=self._cs_log_x.isChecked(), y=self._cs_log_y.isChecked())
        lx = " [log]" if self._cs_log_x.isChecked() else ""
        ly = " [log]" if self._cs_log_y.isChecked() else ""
        self._cs_plot.setLabel("bottom", f"Energy [eV]{lx}")
        self._cs_plot.setLabel("left", f"Cross Section [m^2]{ly}")
        self._plot_cross_sections()

    def _update_rate_axes(self):
        self._rate_plot.setLogMode(x=self._rate_log_x.isChecked(), y=self._rate_log_y.isChecked())
        lx = " [log]" if self._rate_log_x.isChecked() else ""
        ly = " [log]" if self._rate_log_y.isChecked() else ""
        self._rate_plot.setLabel("bottom", f"Electron Temperature [eV]{lx}")
        self._rate_plot.setLabel("left", f"Rate Coefficient [m^3/s]{ly}")
        self._plot_rate_coefficients()

    # ── Details Panel ────────────────────────────────────────

    def _on_cs_focus(self, row):
        if 0 <= row < len(self._cs_items):
            self._details.setHtml(_build_details(self._cs_items[row]))

    def _on_rate_focus(self, row):
        if 0 <= row < len(self._rate_items):
            self._details.setHtml(_build_details(self._rate_items[row]))

    # ── Selection -> Plot ────────────────────────────────────

    def _on_cs_selection(self):
        self._plot_cross_sections()
        sel = self._cs_list.selectedIndexes()
        if len(sel) == 1 and sel[0].row() < len(self._cs_items):
            self._details.setHtml(_build_details(self._cs_items[sel[0].row()]))

    def _on_rate_selection(self):
        self._plot_rate_coefficients()
        sel = self._rate_list.selectedIndexes()
        if len(sel) == 1 and sel[0].row() < len(self._rate_items):
            self._details.setHtml(_build_details(self._rate_items[sel[0].row()]))

    # ── Plotting ─────────────────────────────────────────────

    @staticmethod
    def _filter_positive(x, y):
        """Filter out non-positive values for log-safe plotting."""
        xf, yf = [], []
        for xi, yi in zip(x, y):
            if xi > 0 and yi > 0:
                xf.append(xi)
                yf.append(yi)
        return xf, yf

    def _plot_cross_sections(self):
        self._cs_plot.clear()
        selected = self._cs_list.selectedIndexes()
        if not selected or not self._cs_items:
            return

        legend = self._cs_plot.addLegend(offset=(10, 10))
        plotted = 0
        for sel in selected:
            idx = sel.row()
            if idx >= len(self._cs_items):
                continue
            item = self._cs_items[idx]
            result = _read_csv_2col(Path(item["path"]))
            if result is None:
                continue
            x, y, meta = result
            xf, yf = self._filter_positive(x, y)
            if not xf:
                continue
            color = COLORS[plotted % len(COLORS)]
            self._cs_plot.plot(xf, yf, pen=pg.mkPen(color, width=2), name=item["label"])
            plotted += 1

        self.statusBar().showMessage(f"{plotted} cross section(s) plotted")

    def _plot_rate_coefficients(self):
        self._rate_plot.clear()
        selected = self._rate_list.selectedIndexes()
        if not selected or not self._rate_items:
            return

        legend = self._rate_plot.addLegend(offset=(10, 10))
        plotted = 0
        for sel in selected:
            idx = sel.row()
            if idx >= len(self._rate_items):
                continue
            item = self._rate_items[idx]
            if item.get("representation") in ("arrhenius", "constant"):
                x, y = compute_arrhenius_curve(item.get("A", 0), item.get("E_a_eV", 0),
                                                item.get("model", "arrhenius"))
            else:
                result = _read_csv_2col(Path(item.get("path", "")))
                if result is None:
                    continue
                x, y, meta = result
            xf, yf = self._filter_positive(x, y)
            if not xf:
                continue
            color = COLORS[plotted % len(COLORS)]
            self._rate_plot.plot(xf, yf, pen=pg.mkPen(color, width=2), name=item["label"])
            plotted += 1

        self.statusBar().showMessage(f"{plotted} rate coefficient(s) plotted")
