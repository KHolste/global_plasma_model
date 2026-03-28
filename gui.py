"""
gui_realtime_stream_pyqt.py – moderne PyQt6-GUI mit Echtzeit-Streaming über stdout

Konzept:
- Kein live_results.txt für Live-Plot mehr
- Keine externe live_plot.py nötig
- Die GUI liest die strukturierte C++-Ausgabe direkt aus stdout/stderr
- RESULT-Zeilen werden sofort in eingebettete Live-Plots übernommen

Voraussetzungen:
    pip install PyQt6 pyqtgraph

Optional:
- g++ oder WSL zum Kompilieren des C++-Solvers
"""

from __future__ import annotations

import os
import sys
import shutil
import signal
import subprocess
from pathlib import Path

from PyQt6.QtCore import Qt, QProcess
from PyQt6.QtGui import QAction, QDesktopServices, QTextCursor
from PyQt6.QtWidgets import (
    QApplication,
    QButtonGroup,
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

# Solver-IDs fuer params.txt
# WICHTIG: Newton ist hier auf Methode 4 gemappt.
# Falls dein C++-Backend einen anderen Code erwartet, nur hier anpassen.
SOLVER_OPTIONS = [
    (4, "Newton (stationär)"),
    (1, "Euler (RK1)"),
    (2, "RK4"),
    (3, "RK45 adaptiv"),
]
DEFAULT_SOLVER_ID = 4

PRIMARY_PARAMS = {
    "R":         ("Radius Entladungsgefäß", "m", 0.02),
    "L":         ("Länge Entladungsgefäß", "m", 0.04),
    "betai":     ("Transparenz Ionen", "–", 0.5),
    "betag":     ("Transparenz Neutralgas", "–", 0.05145),
    "frequency": ("Anregefrequenz", "Hz", 2.5e6),
    "Nw":        ("Anzahl Spulenwindungen", "–", 6.0),
    "R_ohm":     ("Ohmscher Widerstand Spule", "Ω", 0.36),
    "Rc":        ("Spulenradius", "m", 0.02),
    "lc":        ("Spulenlänge", "m", 0.04),
    "Vgrid":     ("Screengitter-Spannung", "V", 1500.0),
    "sgrid":     ("Gitterabstand", "m", 0.001),
    "P_RFG":     ("DC-Leistung RFG", "W", 18.0),
    "P_RFG_max": ("Max. RF-Leistung", "W", 80.0),
    "Q0sccm":    ("Massenfluss", "sccm", 0.475),
}

SWEEP_PARAMS = {
    "Q0sccm_start": ("Q0sccm Start", "sccm", 0.27),
    "Q0sccm_step":  ("Q0sccm Schritt", "sccm", 0.01),
    "jjmax":        ("Anzahl Schritte", "–", 73),
    "I_soll":       ("Ziel-Strahlstrom", "mA", 15.0),
}


class MetricCard(QFrame):
    def __init__(self, title: str, unit: str):
        super().__init__()
        self.setObjectName("MetricCard")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 8, 10, 8)
        layout.setSpacing(4)

        title_lbl = QLabel(f"{title} [{unit}]")
        title_lbl.setObjectName("MetricTitle")
        self.value_lbl = QLabel("–")
        self.value_lbl.setObjectName("MetricValue")
        self.value_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)

        layout.addWidget(title_lbl)
        layout.addWidget(self.value_lbl)

    def set_value(self, value: str) -> None:
        self.value_lbl.setText(value)


class LivePlotGrid(QGroupBox):
    def __init__(self):
        super().__init__("Echtzeit-Streaming – Live-Plots")
        layout = QGridLayout(self)
        layout.setHorizontalSpacing(10)
        layout.setVerticalSpacing(10)

        pg.setConfigOptions(antialias=True, background="#0b1220", foreground="#d7e3ff")

        self.data = {
            "first_point": None,
            "Q0sccm": [],
            "P_RFG": [],
            "Te": [],
            "Tg": [],
            "n": [],
            "ng": [],
            "I_mA": [],
        }

        defs = [
            ("P_RFG", "P_solution (W)"),
            ("Te", "Tₑ (eV)"),
            ("Tg", "Tg (K)"),
            ("n", "n (m⁻³)"),
            ("ng", "ng (m⁻³)"),
            ("I_mA", "I_Strahl (mA)"),
        ]

        self.widgets = {}
        self.curves = {}
        self.markers = {}
        self.first_markers = {}

        for idx, (key, title) in enumerate(defs):
            w = pg.PlotWidget()
            w.setMinimumHeight(220)
            w.showGrid(x=True, y=True, alpha=0.20)
            w.setTitle(title)
            w.setLabel("bottom", "Q₀ (sccm)")
            w.getPlotItem().getAxis("left").enableAutoSIPrefix(True)
            w.getPlotItem().getAxis("bottom").enableAutoSIPrefix(False)

            pen = pg.mkPen(width=2)
            curve = w.plot([], [], pen=pen, symbol="o", symbolSize=6)
            marker = pg.ScatterPlotItem(size=10, brush=pg.mkBrush("w"), pen=pg.mkPen("w"))
            first_marker = pg.ScatterPlotItem(size=14, brush=pg.mkBrush("#ff6b6b"), pen=pg.mkPen("#ff6b6b"))
            w.addItem(marker)
            w.addItem(first_marker)

            if key in ("n", "ng"):
                w.getPlotItem().getAxis("left").setLabel(text=title)

            self.widgets[key] = w
            self.curves[key] = curve
            self.markers[key] = marker
            self.first_markers[key] = first_marker
            layout.addWidget(w, idx // 3, idx % 3)

    def clear(self):
        self.data["first_point"] = None
        for key in ("Q0sccm", "P_RFG", "Te", "Tg", "n", "ng", "I_mA"):
            self.data[key].clear()
        self.redraw()

    def add_point(self, q0: float, p_rfg: float, te: float, tg: float, n: float, ng: float, i_ma: float):
        # first point special handling
        if self.data.get("first_point") is None:
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
        fp = self.data.get("first_point")

        fp_map = None
        if fp is not None:
            q0_fp, p_rfg_fp, te_fp, tg_fp, n_fp, ng_fp, i_ma_fp = fp
            fp_map = {
                "P_RFG": (q0_fp, p_rfg_fp),
                "Te": (q0_fp, te_fp),
                "Tg": (q0_fp, tg_fp),
                "n": (q0_fp, n_fp),
                "ng": (q0_fp, ng_fp),
                "I_mA": (q0_fp, i_ma_fp),
            }

        for key in self.curves:
            y = self.data[key]

            # main sweep curve (without first point)
            self.curves[key].setData(x, y)

            # last sweep point
            if x:
                self.markers[key].setData([{"pos": (x[-1], y[-1])}])
            else:
                self.markers[key].setData([])

            # highlighted first point
            if fp_map is not None:
                q0v, yv = fp_map[key]
                self.first_markers[key].setData([{"pos": (q0v, yv)}])
            else:
                self.first_markers[key].setData([])


class SimulatorWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Xenon Plasma Simulator – Realtime Streaming UI")
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

    def _build_ui(self) -> None:
        root = QWidget()
        self.setCentralWidget(root)
        outer = QVBoxLayout(root)
        outer.setContentsMargins(14, 14, 14, 14)
        outer.setSpacing(12)

        header = QFrame()
        header.setObjectName("Header")
        header_layout = QVBoxLayout(header)
        header_layout.setContentsMargins(18, 16, 18, 16)
        title = QLabel("🛰 Global Xenon Model")
        title.setObjectName("HeaderTitle")
        subtitle = QLabel("Echtzeit-Streaming direkt aus stdout – keine Live-Dateien mehr")
        subtitle.setObjectName("HeaderSubtitle")
        header_layout.addWidget(title)
        header_layout.addWidget(subtitle)
        outer.addWidget(header)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setChildrenCollapsible(False)
        outer.addWidget(splitter, 1)

        left_scroll = QScrollArea()
        left_scroll.setWidgetResizable(True)
        left_panel = QWidget()
        left_scroll.setWidget(left_panel)
        left_layout = QVBoxLayout(left_panel)
        left_layout.setContentsMargins(0, 0, 8, 0)
        left_layout.setSpacing(12)

        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        right_layout.setContentsMargins(8, 0, 0, 0)
        right_layout.setSpacing(12)

        splitter.addWidget(left_scroll)
        splitter.addWidget(right_panel)
        splitter.setSizes([720, 860])

        left_layout.addWidget(self._build_params_group())
        left_layout.addWidget(self._build_solver_group())
        left_layout.addWidget(self._build_sweep_group())
        left_layout.addWidget(self._build_button_row())
        left_layout.addStretch(1)

        right_layout.addWidget(self._build_progress_group())
        right_layout.addWidget(self._build_metrics_group())
        self.plot_grid = LivePlotGrid()
        right_layout.addWidget(self.plot_grid, 1)
        right_layout.addWidget(self._build_log_group(), 1)

        status = QStatusBar()
        self.setStatusBar(status)
        self.set_status("Bereit.")

        for text, handler in [
            ("Defaults", self.reset_defaults),
            ("Kompilieren", self.compile_cpp),
            ("Simulation starten", self.run_sim),
            ("Abbrechen", self.cancel_sim),
            ("Output öffnen", self.open_output),
            ("Plots leeren", self.clear_plots),
        ]:
            action = QAction(text, self)
            action.triggered.connect(handler)
            self.addAction(action)

    def _build_params_group(self) -> QGroupBox:
        box = QGroupBox("Simulationsparameter")
        grid = QGridLayout(box)
        grid.setHorizontalSpacing(12)
        grid.setVerticalSpacing(8)

        for col, header in enumerate(["Parameter", "Wert", "Einheit", "Schlüssel"]):
            lbl = QLabel(header)
            lbl.setObjectName("SectionHeader")
            grid.addWidget(lbl, 0, col)

        for i, (key, (label, unit, default)) in enumerate(PRIMARY_PARAMS.items(), start=1):
            grid.addWidget(QLabel(label), i, 0)
            edit = QLineEdit(str(default))
            edit.setAlignment(Qt.AlignmentFlag.AlignRight)
            grid.addWidget(edit, i, 1)
            grid.addWidget(QLabel(unit), i, 2)
            key_lbl = QLabel(f"({key})")
            key_lbl.setObjectName("MutedLabel")
            grid.addWidget(key_lbl, i, 3)
            self.entries[key] = edit
        return box

    def _build_solver_group(self) -> QGroupBox:
        box = QGroupBox("Solver / Verfahren")
        layout = QHBoxLayout(box)

        self.solver_group = QButtonGroup(self)
        self.solver_buttons: dict[int, QRadioButton] = {}

        for value, text in SOLVER_OPTIONS:
            btn = QRadioButton(text)
            if value == DEFAULT_SOLVER_ID:
                btn.setChecked(True)
            self.solver_group.addButton(btn, value)
            self.solver_buttons[value] = btn
            layout.addWidget(btn)

        layout.addStretch(1)
        return box

    def _build_sweep_group(self) -> QGroupBox:
        box = QGroupBox("Scan / Sweep")
        grid = QGridLayout(box)
        grid.setHorizontalSpacing(12)
        grid.setVerticalSpacing(8)

        for i, (key, (label, unit, default)) in enumerate(SWEEP_PARAMS.items()):
            grid.addWidget(QLabel(label), 0, i * 3)
            edit = QLineEdit(str(default))
            edit.setAlignment(Qt.AlignmentFlag.AlignRight)
            grid.addWidget(edit, 0, i * 3 + 1)
            grid.addWidget(QLabel(unit), 0, i * 3 + 2)
            self.sweep_entries[key] = edit
        return box

    def _make_button(self, text: str, slot, primary: bool = False, danger: bool = False) -> QPushButton:
        btn = QPushButton(text)
        btn.clicked.connect(slot)
        if primary:
            btn.setProperty("primary", True)
        if danger:
            btn.setProperty("danger", True)
        btn.setMinimumHeight(40)
        return btn

    def _build_button_row(self) -> QWidget:
        row = QWidget()
        layout = QGridLayout(row)
        layout.setHorizontalSpacing(10)
        layout.setVerticalSpacing(10)

        buttons = [
            self._make_button("🔧 Kompilieren", self.compile_cpp),
            self._make_button("▶ Simulation starten", self.run_sim, primary=True),
            self._make_button("⏹ Abbrechen", self.cancel_sim, danger=True),
            self._make_button("↺ Defaults", self.reset_defaults),
            self._make_button("🧹 Plots leeren", self.clear_plots),
            self._make_button("📄 Output öffnen", self.open_output),
        ]
        for idx, button in enumerate(buttons):
            layout.addWidget(button, 0, idx)
        return row

    def _build_progress_group(self) -> QGroupBox:
        box = QGroupBox("Simulation – Status")
        grid = QGridLayout(box)

        self.scan_progress = QProgressBar()
        self.scan_progress.setRange(0, 100)
        self.scan_label = QLabel("–")

        self.pid_progress = QProgressBar()
        self.pid_progress.setRange(0, 100)
        self.pid_label = QLabel("–")

        grid.addWidget(QLabel("Scan-Fortschritt"), 0, 0)
        grid.addWidget(self.scan_progress, 0, 1)
        grid.addWidget(self.scan_label, 0, 2)

        grid.addWidget(QLabel("PID-Iteration"), 1, 0)
        grid.addWidget(self.pid_progress, 1, 1)
        grid.addWidget(self.pid_label, 1, 2)

        return box

    def _build_metrics_group(self) -> QGroupBox:
        box = QGroupBox("Live-Metriken")
        grid = QGridLayout(box)

        defs = [
            ("Q₀", "Q0sccm", "sccm"),
            ("P_solution", "P_RFG", "W"),
            ("P_set", "P_set", "W"),
            ("I_Strahl", "I_beam", "mA"),
            ("Tₑ", "Te", "eV"),
            ("T_g", "Tg", "K"),
            ("n", "n", "m⁻³"),
        ]

        for i, (title, key, unit) in enumerate(defs):
            card = MetricCard(title, unit)
            grid.addWidget(card, i // 3, i % 3)
            self.metric_cards[key] = card
        return box

    def _build_log_group(self) -> QGroupBox:
        box = QGroupBox("Log")
        layout = QVBoxLayout(box)

        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self.log.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(self.log)
        return box

    def _apply_styles(self) -> None:
        self.setStyleSheet("""
            QMainWindow, QWidget {
                background: #0b1020;
                color: #e6edf7;
                font-size: 13px;
            }
            QFrame#Header {
                background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                    stop:0 #15203b, stop:1 #1f3b73);
                border: 1px solid #2a4b8d;
                border-radius: 16px;
            }
            QLabel#HeaderTitle {
                font-size: 28px;
                font-weight: 700;
                color: white;
            }
            QLabel#HeaderSubtitle {
                color: #c9d6ee;
                font-size: 14px;
            }
            QGroupBox {
                border: 1px solid #24314f;
                border-radius: 14px;
                margin-top: 10px;
                padding-top: 12px;
                background: #0f172a;
                font-weight: 600;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 12px;
                padding: 0 6px;
                color: #d9e4ff;
            }
            QLabel#SectionHeader {
                color: #9fb4dd;
                font-weight: 700;
            }
            QLabel#MutedLabel {
                color: #7f91b2;
            }
            QLineEdit {
                background: #111a31;
                border: 1px solid #31466d;
                border-radius: 10px;
                padding: 8px 10px;
                color: white;
            }
            QPushButton {
                background: #17233f;
                border: 1px solid #31466d;
                border-radius: 12px;
                padding: 8px 14px;
                font-weight: 600;
            }
            QPushButton:hover {
                background: #1e2f56;
            }
            QPushButton[primary="true"] {
                background: #2553d8;
                border-color: #4f7eff;
            }
            QPushButton[primary="true"]:hover {
                background: #2d62fb;
            }
            QPushButton[danger="true"] {
                background: #7e2431;
                border-color: #b44b5b;
            }
            QPushButton[danger="true"]:hover {
                background: #9a2c3c;
            }
            QProgressBar {
                border: 1px solid #31466d;
                border-radius: 9px;
                background: #111a31;
                text-align: center;
                min-height: 22px;
            }
            QProgressBar::chunk {
                background: #2dd4bf;
                border-radius: 8px;
            }
            QTextEdit {
                background: #050812;
                border: 1px solid #24314f;
                border-radius: 12px;
                color: #8ff7c8;
                font-family: Consolas, Menlo, monospace;
                font-size: 12px;
                padding: 6px;
            }
            QFrame#MetricCard {
                background: #111a31;
                border: 1px solid #24314f;
                border-radius: 14px;
            }
            QLabel#MetricTitle {
                color: #8aa2ce;
                font-size: 11px;
            }
            QLabel#MetricValue {
                background: #030712;
                border-radius: 10px;
                padding: 8px;
                color: #7fffd4;
                font-size: 18px;
                font-weight: 700;
                font-family: Consolas, Menlo, monospace;
            }
            QRadioButton {
                spacing: 10px;
                color: #e6edf7;
                font-weight: 500;
                padding: 4px 8px 4px 4px;
            }
            QRadioButton::indicator {
                width: 18px;
                height: 18px;
                border-radius: 9px;
                border: 2px solid #7aa2ff;
                background: #09101f;
            }
            QRadioButton::indicator:hover {
                border: 2px solid #a9c1ff;
                background: #0d1830;
            }
            QRadioButton::indicator:checked {
                border: 2px solid #2dd4bf;
                background: #2dd4bf;
            }
            QRadioButton::indicator:checked:hover {
                border: 2px solid #59f0d7;
                background: #59f0d7;
            }
            QRadioButton:checked {
                color: #2dd4bf;
                font-weight: 700;
            }
        """)

    def log_msg(self, msg: str) -> None:
        self.log.append(msg)
        cursor = self.log.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        self.log.setTextCursor(cursor)

    def set_status(self, msg: str) -> None:
        self.statusBar().showMessage(msg)

    def update_metric(self, key: str, value: str) -> None:
        card = self.metric_cards.get(key)
        if card:
            card.set_value(value)

    def update_solution_power(self, value: float | None) -> None:
        if value is None:
            self.update_metric("P_RFG", "—")
        else:
            self.update_metric("P_RFG", f"{value:.2f}")

    def update_set_power(self, value: float | None) -> None:
        if value is None:
            self.update_metric("P_set", "—")
        else:
            self.update_metric("P_set", f"{value:.2f}")

    def solver_method(self) -> int:
        return self.solver_group.checkedId()

    def solver_label(self) -> str:
        for value, text in SOLVER_OPTIONS:
            if value == self.solver_method():
                return text
        return f"Unbekannt ({self.solver_method()})"

    def wsl_available(self) -> bool:
        return sys.platform == "win32" and shutil.which("wsl") is not None

    def win_to_wsl(self, path: str) -> str:
        path = path.replace("\\", "/")
        if len(path) >= 2 and path[1] == ":":
            drive = path[0].lower()
            rest = path[2:]
            return f"/mnt/{drive}{rest}"
        return path

    def read_float(self, line_edit: QLineEdit, name: str) -> float:
        try:
            return float(line_edit.text().strip())
        except ValueError as exc:
            raise ValueError(f"{name}: ungültiger numerischer Wert") from exc

    def write_config(self) -> bool:
        try:
            with open(CONFIG_FILE, "w", encoding="utf-8") as f:
                f.write("# Automatisch generiert von gui_realtime_stream_pyqt.py\n")
                for key, edit in self.entries.items():
                    value = self.read_float(edit, key)
                    f.write(f"{key} {value}\n")
                for key, edit in self.sweep_entries.items():
                    value = self.read_float(edit, key)
                    f.write(f"{key} {value}\n")
                f.write(f"method {self.solver_method()}\n")
            self.log_msg(f"📝 {CONFIG_FILE} geschrieben. Solver={self.solver_label()} (method={self.solver_method()})")
            return True
        except ValueError as exc:
            QMessageBox.critical(self, "Eingabefehler", str(exc))
            return False
        except OSError as exc:
            QMessageBox.critical(self, "Dateifehler", str(exc))
            return False

    def compile_cpp(self) -> None:
        self.set_status("Kompiliere …")
        self.log_msg("⏳ Kompilierung gestartet …")

        if self.wsl_available():
            cwd_wsl = self.win_to_wsl(str(SCRIPT_DIR))
            cmd = (
                f'cd "{cwd_wsl}" && '
                f'g++ -O3 -march=native -std=c++17 {CPP_SOURCE} -o chabert -lboost_system 2>&1'
            )
            result = subprocess.run(["wsl", "bash", "-c", cmd], capture_output=True, text=True)
            output = (result.stdout or "") + (result.stderr or "")
            if result.returncode == 0:
                self.log_msg("✅ Kompilierung erfolgreich (WSL).")
                self.set_status("Kompilierung erfolgreich.")
            else:
                self.log_msg("❌ Kompilierfehler:\n" + output)
                self.set_status("Fehler beim Kompilieren.")
                QMessageBox.critical(self, "Kompilierfehler", output[:2000] or "Unbekannter Fehler")
            return

        compiler = shutil.which("g++")
        if not compiler:
            msg = "g++ wurde nicht gefunden und WSL ist nicht verfügbar."
            self.log_msg("❌ " + msg)
            self.set_status("Kein Compiler gefunden.")
            QMessageBox.critical(self, "Compiler nicht gefunden", msg)
            return

        out_name = "chabert.exe" if sys.platform == "win32" else "chabert"
        cmd = [compiler, "-O3", "-march=native", "-std=c++17", CPP_SOURCE, "-o", out_name, "-lboost_system"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            self.log_msg("✅ Kompilierung erfolgreich.")
            self.set_status("Kompilierung erfolgreich.")
        else:
            self.log_msg("❌ Kompilierfehler:\n" + (result.stderr or result.stdout))
            self.set_status("Fehler beim Kompilieren.")
            QMessageBox.critical(self, "Kompilierfehler", (result.stderr or result.stdout)[:2000])

    def run_sim(self) -> None:
        if self.process and self.process.state() != QProcess.ProcessState.NotRunning:
            QMessageBox.information(self, "Läuft bereits", "Es läuft bereits eine Simulation.")
            return

        binary_name = "chabert.exe" if sys.platform == "win32" else "chabert"

        if self.wsl_available():
            cwd_wsl = self.win_to_wsl(str(SCRIPT_DIR))
            check = subprocess.run(["wsl", "bash", "-c", f'test -f "{cwd_wsl}/chabert"'])
            if check.returncode != 0:
                QMessageBox.warning(self, "Nicht gefunden", "WSL-Binary 'chabert' nicht gefunden. Bitte zuerst kompilieren.")
                return
        else:
            binary_path = SCRIPT_DIR / binary_name
            if not binary_path.exists():
                QMessageBox.warning(self, "Nicht gefunden", f"'{binary_path.name}' nicht gefunden. Bitte zuerst kompilieren.")
                return

        if not self.write_config():
            return

        self.cancel_requested = False
        self._stdout_buffer = ""
        self.current_q0 = None
        self.scan_progress.setValue(0)
        self.pid_progress.setValue(0)
        self.scan_label.setText("–")
        self.pid_label.setText("–")
        self.clear_plots()

        self.process = QProcess(self)
        self.process.setWorkingDirectory(str(SCRIPT_DIR))
        self.process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.process.readyReadStandardOutput.connect(self.on_process_output)
        self.process.finished.connect(self.on_process_finished)

        self.log_msg(f"\n🚀 Starte Simulation …  Solver={self.solver_label()} (method={self.solver_method()})")
        self.set_status("Simulation läuft …")

        if self.wsl_available():
            cwd_wsl = self.win_to_wsl(str(SCRIPT_DIR))
            bash_cmd = f'cd "{cwd_wsl}" && ./chabert {CONFIG_FILE} 2>&1'
            self.process.start("wsl", ["bash", "-c", bash_cmd])
        else:
            program = str(SCRIPT_DIR / binary_name)
            self.process.start(program, [CONFIG_FILE])

        if not self.process.waitForStarted(3000):
            QMessageBox.critical(self, "Startfehler", "Simulation konnte nicht gestartet werden.")
            self.set_status("Startfehler.")
            self.log_msg("❌ Prozessstart fehlgeschlagen.")

    def on_process_output(self) -> None:
        if not self.process:
            return
        text = bytes(self.process.readAllStandardOutput()).decode("utf-8", errors="replace")
        self._stdout_buffer += text

        while "\n" in self._stdout_buffer:
            line, self._stdout_buffer = self._stdout_buffer.split("\n", 1)
            self.parse_cpp_line(line.strip())

    def parse_cpp_line(self, line: str) -> None:
        if not line:
            return
        parts = line.split()
        tag = parts[0]

        try:
            if tag == "Q0_STEP" and len(parts) >= 4:
                q0, step, total = parts[1], parts[2], parts[3]
                self.current_q0 = float(q0)
                pct = int(step) / int(total) * 100
                self.scan_progress.setValue(int(pct))
                self.scan_label.setText(f"Schritt {step}/{total}")
                self.update_metric("Q0sccm", f"{float(q0):.4f}")
                self.pid_progress.setValue(0)
                self.pid_label.setText("–")
                self.log_msg(f"▶ Q0={q0} sccm ({step}/{total})")
                return

            if tag == "PID_START" and len(parts) >= 3:
                it, prfg = parts[1], parts[2]
                pct = min(int(it) / 50 * 100, 100)
                self.pid_progress.setValue(int(pct))
                self.pid_label.setText(f"Iter {it}")
                self.update_set_power(float(prfg))
                return

            if tag == "PID_DONE" and len(parts) >= 6:
                I, err, prfg, te, tg = parts[1], parts[2], parts[3], parts[4], parts[5]
                self.update_metric("I_beam", f"{float(I):.3f}")
                self.update_metric("Te", f"{float(te):.2f}")
                self.update_metric("Tg", f"{float(tg):.1f}")
                self.update_set_power(float(prfg))
                color = "🟢" if abs(float(err)) < 0.1 else ("🟡" if abs(float(err)) < 1 else "🔴")
                self.log_msg(f"  {color} I={float(I):.3f} mA  err={float(err):+.3f} mA  P={float(prfg):.2f} W")
                return

            if tag == "CONVERGED" and len(parts) >= 2:
                self.pid_progress.setValue(100)
                self.pid_label.setText("✅ konvergiert")
                self.log_msg(f"  ✅ Konvergiert nach {parts[1]} PID-Iterationen")
                return

            if tag == "RESULT" and len(parts) >= 7:
                n, ng, te, tg, i_beam_mA, p_rfg = parts[1], parts[2], parts[3], parts[4], parts[5], parts[6]
                n_f = float(n)
                ng_f = float(ng)
                te_f = float(te)
                tg_f = float(tg)
                i_f = float(i_beam_mA)
                p_f = float(p_rfg)

                self.update_metric("n", f"{n_f:.2e}")
                self.update_metric("Te", f"{te_f:.2f}")
                self.update_metric("Tg", f"{tg_f:.1f}")
                self.update_metric("I_beam", f"{i_f:.3f}")
                self.update_solution_power(p_f)
                self.update_set_power(p_f)

                if self.current_q0 is not None:
                    self.plot_grid.add_point(self.current_q0, p_f, te_f, tg_f, n_f, ng_f, i_f)
                return

            if tag == "PID_MAXITER":
                self.log_msg("  ⚠ PID: max. Iterationen ohne Konvergenz")
                return

            if tag == "P_LIMIT_REACHED" and len(parts) >= 6:
                q0, p_now, p_max, i_beam, err = parts[1], parts[2], parts[3], parts[4], parts[5]
                self.log_msg(
                    f"  ⛔ RF-Limit erreicht bei Q0={float(q0):.4f} sccm | "
                    f"P={float(p_now):.2f} W > Pmax={float(p_max):.2f} W | "
                    f"I={float(i_beam):.3f} mA | err={float(err):+.3f} mA"
                )
                self.update_set_power(float(p_max))
                return

            if tag == "SWEEP_RECOVERY" and len(parts) >= 2:
                self.log_msg("  ♻ " + " ".join(parts[1:]))
                return

            if tag == "SWEEP_RECOVERY_FAIL" and len(parts) >= 2:
                self.log_msg("  ❌ " + " ".join(parts[1:]))
                return

            if tag == "SOLVER_FAIL" and len(parts) >= 4:
                self.update_solution_power(None)
                self.log_msg("  ❌ " + line)
                return

            if not line.startswith("STEP "):
                self.log_msg(line)

        except Exception as exc:
            self.log_msg(f"⚠ Parse-Fehler bei Zeile: {line} | {exc}")

    def on_process_finished(self, exit_code: int, _status) -> None:
        if self.cancel_requested:
            self.log_msg("\n⏹ Simulation wurde abgebrochen.")
            self.set_status("Abgebrochen.")
            self.cancel_requested = False
            return

        if exit_code == 0:
            self.log_msg("\n✅ Simulation abgeschlossen.")
            self.set_status("✅ Fertig.")
            self.scan_progress.setValue(100)
            if Path(OUTPUT_FILE).exists():
                try:
                    lines = Path(OUTPUT_FILE).read_text(encoding="utf-8").splitlines()
                    self.log_msg(f"\n📄 {OUTPUT_FILE} ({max(0, len(lines)-1)} Datenpunkte gespeichert)")
                except Exception:
                    self.log_msg(f"\n📄 {OUTPUT_FILE} gespeichert")
        else:
            self.log_msg(f"\n❌ Simulation mit Fehler beendet (Code {exit_code}).")
            self.set_status(f"❌ Fehler (Code {exit_code}).")

    def cancel_sim(self) -> None:
        if self.process and self.process.state() != QProcess.ProcessState.NotRunning:
            self.cancel_requested = True
            self.process.kill()
            self.log_msg("⏹ Abbruch angefordert …")
            self.set_status("Abbruch läuft …")
        else:
            self.log_msg("ℹ Kein laufender Prozess.")

    def reset_defaults(self) -> None:
        for key, (_, _, default) in PRIMARY_PARAMS.items():
            self.entries[key].setText(str(default))
        for key, (_, _, default) in SWEEP_PARAMS.items():
            self.sweep_entries[key].setText(str(default))
        self.solver_buttons[DEFAULT_SOLVER_ID].setChecked(True)
        self.log_msg("↺ Standardwerte geladen.")

    def clear_plots(self) -> None:
        self.plot_grid.clear()
        self.log_msg("🧹 Live-Plots geleert.")

    def open_output(self) -> None:
        path = SCRIPT_DIR / OUTPUT_FILE
        if not path.exists():
            QMessageBox.information(self, "Nicht vorhanden", f"'{OUTPUT_FILE}' existiert noch nicht.")
            return
        QDesktopServices.openUrl(path.as_uri())

    def closeEvent(self, event) -> None:
        try:
            if self.process and self.process.state() != QProcess.ProcessState.NotRunning:
                self.cancel_requested = True
                self.process.kill()
                self.process.waitForFinished(1500)
        except Exception:
            pass
        super().closeEvent(event)


def main() -> int:
    app = QApplication(sys.argv)
    window = SimulatorWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
