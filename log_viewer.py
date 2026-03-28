"""
log_viewer.py – Standalone Log-Viewer fuer Global Xenon Model Masterlogs.

Features:
  - Laedt die strukturierte Masterlog-Datei (simulation_log_*.txt)
  - Vordefinierte Einzel- und Multi-Kurven-Plots
  - Frei waehlbare x/y-Achsen mit Multi-Select fuer y
  - CL-Limit-Marker bei RF-Power-Plots
  - Legende bei Multi-Kurven
  - Log-Achsen (x/y unabhaengig umschaltbar)
  - Benutzerdefinierte Kurven-Stile (Linienstaerke, Marker, Farbe)
  - Standalone oder aus der Haupt-GUI startbar

Start:
    python log_viewer.py [logfile.txt]

Voraussetzungen:
    pip install PyQt6 pyqtgraph
"""
from __future__ import annotations

import math
import sys
import re
from pathlib import Path
from typing import Optional

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QSplitter, QPushButton, QComboBox, QLabel, QFileDialog,
    QGroupBox, QTextEdit, QStatusBar, QCheckBox, QListWidget,
    QAbstractItemView, QScrollArea, QFrame, QSpinBox,
    QTabWidget, QTableWidget, QTableWidgetItem, QHeaderView,
)
import pyqtgraph as pg

# ─── Farben fuer Multi-Kurven ────────────────────────────────

CURVE_COLORS = [
    "#2dd4bf", "#7aa2ff", "#ff6b6b", "#fbbf24", "#a78bfa",
    "#34d399", "#f472b6", "#60a5fa", "#fb923c", "#c084fc",
]

LINE_STYLES = {
    "Solid":    Qt.PenStyle.SolidLine,
    "Dash":     Qt.PenStyle.DashLine,
    "Dot":      Qt.PenStyle.DotLine,
    "DashDot":  Qt.PenStyle.DashDotLine,
}

# ─── Log Parser ───────────────────────────────────────────────

def parse_master_log(path: str) -> dict:
    """Parse eine Masterlog-Datei."""
    result: dict = {
        "metadata": {}, "params": {}, "columns": [], "data": [],
        "events": [], "summary": {},
        "cl_limit_I_mA": None, "cl_limit_J": None,
    }

    text = Path(path).read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    columns: list[str] = []

    for line in lines:
        s = line.strip()

        for key in ("timestamp_start:", "timestamp_end:", "runtime_seconds:",
                     "solver_mode:", "config_file:"):
            if s.startswith(key):
                result["metadata"][key.rstrip(":")] = s[len(key):].strip()

        if s.startswith("# CL_LIMIT_I_mA="):
            try: result["cl_limit_I_mA"] = float(s.split("=", 1)[1])
            except ValueError: pass
        if s.startswith("# CL_LIMIT_J_A_per_m2="):
            try: result["cl_limit_J"] = float(s.split("=", 1)[1])
            except ValueError: pass

        m = re.match(r"^(.+?)\|\s*([\d.eE+\-]+)\s*\|\s*(\S+)\s*\|\s*(\S+)", s)
        if m and not s.startswith("idx") and not s.startswith("-"):
            result["params"][m.group(4).strip()] = m.group(2).strip()

        if s.startswith("DATA_HEADER|"):
            columns = [c.strip() for c in s.split("|")[1:]]
            result["columns"] = columns

        if s.startswith("DATA|") and columns:
            parts = s.split("|")[1:]
            row: dict = {}
            for i, col in enumerate(columns):
                if i < len(parts):
                    val = parts[i].strip()
                    try: row[col] = float(val)
                    except (ValueError, TypeError): row[col] = val
                else:
                    row[col] = None
            result["data"].append(row)

        if s.startswith("- ") and not s.startswith("--"):
            result["events"].append(s[2:])

        for skey in ("total_points", "converged", "no_physical_solution",
                     "numerical_fail", "runtime_seconds"):
            if s.startswith(skey):
                sm = re.match(rf"^{re.escape(skey)}\s*\|\s*(.+)", s)
                if sm: result["summary"][skey] = sm.group(1).strip()

    return result


def get_numeric_columns(parsed: dict) -> list[str]:
    if not parsed["data"]:
        return []
    numeric = []
    for col in parsed["columns"]:
        for row in parsed["data"]:
            v = row.get(col)
            if isinstance(v, (int, float)) and v != 0:
                numeric.append(col)
                break
    return numeric


def get_column_data(parsed: dict, col: str, only_converged: bool = True) -> list[float]:
    vals = []
    for row in parsed["data"]:
        if only_converged and row.get("status") != "CONVERGED":
            continue
        v = row.get(col)
        vals.append(float(v) if isinstance(v, (int, float)) else float("nan"))
    return vals


# ─── Datenfilter fuer Log-Achsen ──────────────────────────────

def filter_for_log(xvals: list[float], yvals: list[float],
                   log_x: bool, log_y: bool) -> tuple[list[float], list[float], int]:
    """Entferne Punkte, die bei Log-Darstellung ungueltig waeren (<=0).

    Returns (x_filtered, y_filtered, n_removed).
    """
    xf, yf = [], []
    removed = 0
    for x, y in zip(xvals, yvals):
        x_ok = (not log_x) or (math.isfinite(x) and x > 0)
        y_ok = (not log_y) or (math.isfinite(y) and y > 0)
        if x_ok and y_ok:
            xf.append(x)
            yf.append(y)
        else:
            removed += 1
    return xf, yf, removed


# ─── Vordefinierte Plot-Gruppen ───────────────────────────────

PLOT_GROUPS = [
    ("Q0sccm", ["thrust_total_mN", "thrust_ions_mN", "thrust_atoms_mN"],
     "Thrust vs Flow (all)"),
    ("Q0sccm", ["thrust_total_mN"],  "Total Thrust vs Flow"),
    ("Q0sccm", ["Te_eV"],            "Te vs Flow"),
    ("Q0sccm", ["Tg_K"],             "Tg vs Flow"),
    ("Q0sccm", ["n_m3", "ng_m3"],    "Densities vs Flow"),
    ("Q0sccm", ["I_mA"],             "Beam Current vs Flow"),
    ("Q0sccm", ["P_RFG_W"],          "RF Power vs Flow"),
    ("Q0sccm", ["iondeg_pct"],       "Ionization Degree vs Flow"),
    ("Q0sccm", ["icp_power_efficiency", "gamma_thrust_eff", "eta_mass_util"],
     "Efficiencies vs Flow"),
    ("P_RFG_W", ["n_m3", "ng_m3"],   "Densities vs RF Power"),
    ("P_RFG_W", ["Te_eV", "Tg_K"],   "Temperatures vs RF Power"),
    ("P_RFG_W", ["Te_eV"],           "Te vs RF Power"),
    ("P_RFG_W", ["Tg_K"],            "Tg vs RF Power"),
    ("P_RFG_W", ["I_mA"],            "Beam Current vs RF Power"),
    ("P_RFG_W", ["icp_power_efficiency", "gamma_thrust_eff", "eta_mass_util"],
     "Efficiencies vs RF Power"),
    ("P_RFG_W", ["xi_mN_per_kW"],    "Thrust Efficiency vs RF Power"),
    ("P_RFG_W", ["thrust_total_mN", "thrust_ions_mN", "thrust_atoms_mN"],
     "Thrust vs RF Power"),
]


# ─── GUI ──────────────────────────────────────────────────────

class LogViewerWindow(QMainWindow):
    def __init__(self, initial_file: Optional[str] = None):
        super().__init__()
        self.setWindowTitle("Global Xenon Model – Log Viewer")
        self.resize(1300, 850)
        self.parsed: Optional[dict] = None
        self._build_ui()
        self._apply_style()
        if initial_file and Path(initial_file).exists():
            self._load_file(initial_file)

    def _build_ui(self):
        root = QWidget()
        self.setCentralWidget(root)
        outer = QVBoxLayout(root)
        outer.setContentsMargins(8, 8, 8, 8)

        # Toolbar
        toolbar = QHBoxLayout()
        btn_open = QPushButton("Open Log File...")
        btn_open.clicked.connect(self._on_open)
        toolbar.addWidget(btn_open)
        self.file_label = QLabel("No file loaded")
        toolbar.addWidget(self.file_label, 1)
        outer.addLayout(toolbar)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        outer.addWidget(splitter, 1)

        # --- Linke Seite (scrollbar) ---
        left_scroll = QScrollArea()
        left_scroll.setWidgetResizable(True)
        left_scroll.setFrameShape(QFrame.Shape.NoFrame)
        left_inner = QWidget()
        left = QVBoxLayout(left_inner)
        left.setContentsMargins(0, 0, 4, 0)
        left_scroll.setWidget(left_inner)

        # Vordefinierte Plots
        grp_pre = QGroupBox("Predefined Plots")
        gl = QVBoxLayout(grp_pre)
        self.combo_predef = QComboBox()
        self.combo_predef.addItem("-- select --")
        for xcol, ycols, label in PLOT_GROUPS:
            tag = f" [{len(ycols)}]" if len(ycols) > 1 else ""
            self.combo_predef.addItem(f"{label}{tag}", (xcol, ycols, label))
        self.combo_predef.currentIndexChanged.connect(self._on_predef_changed)
        gl.addWidget(self.combo_predef)
        left.addWidget(grp_pre)

        # Custom Multi-Plot
        grp_custom = QGroupBox("Custom Multi-Plot")
        cl = QVBoxLayout(grp_custom)
        cl.addWidget(QLabel("X axis:"))
        self.combo_x = QComboBox()
        cl.addWidget(self.combo_x)
        cl.addWidget(QLabel("Y axes (Ctrl+Click):"))
        self.list_y = QListWidget()
        self.list_y.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.list_y.setMaximumHeight(160)
        cl.addWidget(self.list_y)
        btn_plot = QPushButton("Plot Selected")
        btn_plot.clicked.connect(self._on_custom_plot)
        cl.addWidget(btn_plot)
        self.chk_converged = QCheckBox("Only converged points")
        self.chk_converged.setChecked(True)
        cl.addWidget(self.chk_converged)
        left.addWidget(grp_custom)

        # ── Plot Settings ─────────────────────────────────────
        grp_style = QGroupBox("Plot Settings")
        sl = QVBoxLayout(grp_style)

        # Log-Achsen
        row_log = QHBoxLayout()
        self.chk_log_x = QCheckBox("Log X")
        self.chk_log_y = QCheckBox("Log Y")
        row_log.addWidget(self.chk_log_x)
        row_log.addWidget(self.chk_log_y)
        row_log.addStretch()
        sl.addLayout(row_log)

        # Linienstil
        row_line = QHBoxLayout()
        row_line.addWidget(QLabel("Line width:"))
        self.spin_linewidth = QSpinBox()
        self.spin_linewidth.setRange(1, 8)
        self.spin_linewidth.setValue(2)
        row_line.addWidget(self.spin_linewidth)
        row_line.addWidget(QLabel("Style:"))
        self.combo_linestyle = QComboBox()
        for name in LINE_STYLES:
            self.combo_linestyle.addItem(name)
        row_line.addWidget(self.combo_linestyle)
        sl.addLayout(row_line)

        # Marker
        row_marker = QHBoxLayout()
        self.chk_markers = QCheckBox("Show markers")
        self.chk_markers.setChecked(True)
        row_marker.addWidget(self.chk_markers)
        row_marker.addWidget(QLabel("Size:"))
        self.spin_markersize = QSpinBox()
        self.spin_markersize.setRange(1, 20)
        self.spin_markersize.setValue(4)
        row_marker.addWidget(self.spin_markersize)
        row_marker.addStretch()
        sl.addLayout(row_marker)

        # Farbe (Override fuer alle Kurven oder "Auto")
        row_color = QHBoxLayout()
        row_color.addWidget(QLabel("Color:"))
        self.combo_color = QComboBox()
        self.combo_color.addItem("Auto (per curve)")
        for c in CURVE_COLORS:
            self.combo_color.addItem(c)
        row_color.addWidget(self.combo_color, 1)
        sl.addLayout(row_color)

        left.addWidget(grp_style)

        # Run Info
        grp_info = QGroupBox("Run Info")
        il = QVBoxLayout(grp_info)
        self.info_text = QTextEdit()
        self.info_text.setReadOnly(True)
        self.info_text.setMaximumHeight(180)
        il.addWidget(self.info_text)
        left.addWidget(grp_info)
        left.addStretch()

        splitter.addWidget(left_scroll)

        # --- Rechte Seite: Tabs (Plot + Datentabelle) ---
        self.right_tabs = QTabWidget()

        # Tab 1: Plot
        pg.setConfigOptions(antialias=True, background="#0b1220", foreground="#d7e3ff")
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.showGrid(x=True, y=True, alpha=0.2)
        self.right_tabs.addTab(self.plot_widget, "Plot")

        # Tab 2: Datentabelle
        self.data_table = QTableWidget()
        self.data_table.setAlternatingRowColors(True)
        self.data_table.setSortingEnabled(True)
        self.data_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.data_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.right_tabs.addTab(self.data_table, "Datentabelle")

        splitter.addWidget(self.right_tabs)
        splitter.setSizes([340, 960])
        self.setStatusBar(QStatusBar())

    def _apply_style(self):
        self.setStyleSheet("""
            QMainWindow, QWidget { background: #0b1020; color: #e6edf7; font-size: 10pt; }
            QGroupBox { border: 1px solid #24314f; border-radius: 10px; margin-top: 8px;
                        padding-top: 10px; background: #0f172a; font-weight: 600; }
            QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 4px;
                               color: #d9e4ff; }
            QPushButton { background: #17233f; border: 1px solid #31466d; border-radius: 8px;
                          padding: 6px 12px; font-weight: 600; }
            QPushButton:hover { background: #1e2f56; }
            QComboBox { background: #111a31; border: 1px solid #31466d; border-radius: 8px;
                        padding: 4px 8px; color: white; }
            QTextEdit { background: #050812; border: 1px solid #24314f; border-radius: 8px;
                        color: #8ff7c8; font-family: Consolas, monospace; font-size: 9pt; }
            QLabel { color: #9fb4dd; }
            QCheckBox { color: #e6edf7; }
            QSpinBox { background: #111a31; border: 1px solid #31466d; border-radius: 6px;
                       padding: 2px 6px; color: white; }
            QListWidget { background: #111a31; border: 1px solid #31466d; border-radius: 8px;
                          color: white; font-size: 9pt; }
            QListWidget::item:selected { background: #2553d8; }
            QTableWidget { background: #0b1220; color: #d7e3ff; gridline-color: #24314f;
                           font-family: Consolas, monospace; font-size: 9pt; }
            QTableWidget::item { padding: 2px 4px; }
            QTableWidget::item:selected { background: #2553d8; }
            QHeaderView::section { background: #111a31; color: #9fb4dd; border: 1px solid #24314f;
                                   padding: 3px 6px; font-weight: 600; }
            QTabWidget::pane { border: 1px solid #24314f; }
            QTabBar::tab { background: #111a31; color: #9fb4dd; padding: 6px 16px;
                           border: 1px solid #24314f; border-bottom: none; border-radius: 6px 6px 0 0; }
            QTabBar::tab:selected { background: #0b1220; color: #2dd4bf; font-weight: 700; }
        """)

    # ── Datei laden ───────────────────────────────────────────

    def _on_open(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open Master Log", "",
            "Log files (simulation_log_*.txt);;All (*)")
        if path:
            self._load_file(path)

    def _load_file(self, path: str):
        try:
            self.parsed = parse_master_log(path)
        except Exception as exc:
            self.statusBar().showMessage(f"Fehler: {exc}")
            return

        self.file_label.setText(Path(path).name)
        n_rows = len(self.parsed["data"])
        n_conv = sum(1 for r in self.parsed["data"] if r.get("status") == "CONVERGED")
        self.statusBar().showMessage(f"{n_rows} Punkte geladen ({n_conv} konvergiert)")

        cols = get_numeric_columns(self.parsed)
        self.combo_x.clear()
        self.list_y.clear()
        for c in cols:
            self.combo_x.addItem(c)
            self.list_y.addItem(c)

        info = []
        for k, v in self.parsed["metadata"].items():
            info.append(f"{k}: {v}")
        for k, v in self.parsed["summary"].items():
            info.append(f"{k}: {v}")
        cl = self.parsed.get("cl_limit_I_mA")
        if cl:
            info.append(f"CL limit: {cl:.2f} mA")
        self.info_text.setPlainText("\n".join(info))

        # Datentabelle befuellen
        self._populate_table()

    # ── Datentabelle ──────────────────────────────────────────

    def _populate_table(self):
        """Befuelle die QTableWidget mit den geparsten DATA-Zeilen."""
        if not self.parsed or not self.parsed["data"]:
            self.data_table.clear()
            self.data_table.setRowCount(0)
            self.data_table.setColumnCount(0)
            return

        columns = self.parsed["columns"]
        data = self.parsed["data"]

        self.data_table.setSortingEnabled(False)  # Waehrend Befuellung deaktivieren
        self.data_table.setRowCount(len(data))
        self.data_table.setColumnCount(len(columns))
        self.data_table.setHorizontalHeaderLabels(columns)

        for row_idx, row in enumerate(data):
            for col_idx, col_name in enumerate(columns):
                val = row.get(col_name)
                if val is None:
                    text = ""
                elif isinstance(val, float):
                    # Kompakte Darstellung: wissenschaftlich fuer grosse/kleine Werte
                    if val == 0.0:
                        text = "0"
                    elif abs(val) < 0.01 or abs(val) >= 1e6:
                        text = f"{val:.4e}"
                    else:
                        text = f"{val:.4f}"
                else:
                    text = str(val)

                item = QTableWidgetItem(text)
                # Fuer Sortierung: numerische Werte als UserRole speichern
                if isinstance(val, (int, float)):
                    item.setData(Qt.ItemDataRole.UserRole, val)
                self.data_table.setItem(row_idx, col_idx, item)

        self.data_table.setSortingEnabled(True)
        self.data_table.resizeColumnsToContents()

        # Spaltenbreite begrenzen (sonst werden breite Spalten unleserlich)
        header = self.data_table.horizontalHeader()
        for i in range(len(columns)):
            if self.data_table.columnWidth(i) > 150:
                self.data_table.setColumnWidth(i, 150)
        header.setStretchLastSection(True)

    # ── Plot-Helfer ───────────────────────────────────────────

    def _get_pen(self, color: str, index: int) -> pg.mkPen:
        """Erzeuge Pen mit aktuellen UI-Einstellungen."""
        # Farbe: Auto oder Override
        color_override = self.combo_color.currentText()
        if color_override.startswith("#"):
            c = color_override
        else:
            c = color  # Auto: Farbe aus Zyklusliste

        width = self.spin_linewidth.value()
        style_name = self.combo_linestyle.currentText()
        style = LINE_STYLES.get(style_name, Qt.PenStyle.SolidLine)
        return pg.mkPen(c, width=width, style=style)

    def _get_symbol(self, color: str) -> dict:
        """Symbol-Optionen fuer Plot-Aufruf."""
        if not self.chk_markers.isChecked():
            return {}
        color_override = self.combo_color.currentText()
        c = color_override if color_override.startswith("#") else color
        return {
            "symbol": "o",
            "symbolSize": self.spin_markersize.value(),
            "symbolBrush": c,
        }

    def _safe_clear_plot(self):
        """Plot und Legende sicher zuruecksetzen."""
        plot_item = self.plot_widget.getPlotItem()
        self.plot_widget.clear()
        legend = plot_item.legend
        if legend is not None:
            try:
                scene = legend.scene()
                if scene is not None:
                    scene.removeItem(legend)
            except Exception:
                pass
            plot_item.legend = None

    # ── Haupt-Plot-Funktion ───────────────────────────────────

    def _plot_multi(self, xcol: str, ycols: list[str], title: str = ""):
        """Plotte mehrere y-Spalten gegen eine x-Spalte."""
        if not self.parsed or not ycols:
            return

        only_conv = self.chk_converged.isChecked()
        xdata_raw = get_column_data(self.parsed, xcol, only_conv)
        if not xdata_raw:
            self.statusBar().showMessage("Keine Daten fuer X-Achse")
            return

        log_x = self.chk_log_x.isChecked()
        log_y = self.chk_log_y.isChecked()

        self._safe_clear_plot()
        plot_item = self.plot_widget.getPlotItem()
        plot_item.addLegend(offset=(10, 10))

        # Log-Modus setzen (pyqtgraph transformiert intern)
        plot_item.setLogMode(x=log_x, y=log_y)

        self.plot_widget.setTitle(title or " / ".join(ycols) + f" vs {xcol}")
        self.plot_widget.setLabel("bottom", xcol + (" [log]" if log_x else ""))
        if len(ycols) == 1:
            self.plot_widget.setLabel("left", ycols[0] + (" [log]" if log_y else ""))
        else:
            self.plot_widget.setLabel("left", "[log]" if log_y else "")

        total_plotted = 0
        total_filtered = 0

        for i, ycol in enumerate(ycols):
            ydata_raw = get_column_data(self.parsed, ycol, only_conv)
            if len(ydata_raw) != len(xdata_raw):
                continue

            # Nichtpositive Werte fuer Log-Achsen herausfiltern
            xdata, ydata, n_removed = filter_for_log(
                xdata_raw, ydata_raw, log_x, log_y)
            total_filtered += n_removed

            if not xdata:
                self.statusBar().showMessage(
                    f"{ycol}: alle Punkte ungueltig fuer Log-Darstellung")
                continue

            color = CURVE_COLORS[i % len(CURVE_COLORS)]
            pen = self._get_pen(color, i)
            sym = self._get_symbol(color)

            self.plot_widget.plot(
                xdata, ydata, pen=pen, name=ycol, **sym)
            total_plotted += 1

        # CL-Limit-Markierung bei RF-Power-x-Achse
        if xcol.startswith("P_R") and self.parsed.get("cl_limit_I_mA"):
            cl_i = self.parsed["cl_limit_I_mA"]
            p_cl = None
            for row in self.parsed["data"]:
                if row.get("status") != "CONVERGED":
                    continue
                i_val = row.get("I_mA")
                p_val = row.get(xcol)
                if isinstance(i_val, (int, float)) and isinstance(p_val, (int, float)):
                    if i_val >= cl_i:
                        p_cl = p_val
                        break
            if p_cl is not None:
                # Bei Log-X muss die Position log10-transformiert werden
                pos = math.log10(p_cl) if (log_x and p_cl > 0) else p_cl
                cl_line = pg.InfiniteLine(
                    pos=pos, angle=90,
                    pen=pg.mkPen("#ff6b6b", width=2, style=Qt.PenStyle.DashLine),
                    label=f"CL limit ({p_cl:.0f} W)",
                    labelOpts={"color": "#ff6b6b", "position": 0.95},
                )
                self.plot_widget.addItem(cl_line)

        # Statusmeldung
        msg = f"{total_plotted} Kurve(n), {len(xdata_raw)} Punkte"
        if total_filtered > 0:
            msg += f" ({total_filtered} Punkte gefiltert fuer Log-Skala)"
        self.statusBar().showMessage(msg)

    # ── Event-Handler ─────────────────────────────────────────

    def _on_predef_changed(self, index):
        data = self.combo_predef.itemData(index)
        if data:
            xcol, ycols, label = data
            self._plot_multi(xcol, ycols, label)

    def _on_custom_plot(self):
        xcol = self.combo_x.currentText()
        ycols = [item.text() for item in self.list_y.selectedItems()]
        if xcol and ycols:
            self._plot_multi(xcol, ycols)
        elif xcol:
            self.statusBar().showMessage("Bitte Y-Spalte(n) auswaehlen")


# ─── Main ─────────────────────────────────────────────────────

def main():
    app = QApplication(sys.argv)
    initial = sys.argv[1] if len(sys.argv) > 1 else None
    win = LogViewerWindow(initial)
    win.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
