"""
run_logger.py -- Strukturiertes Run-Log fuer den Python-Solver.

Erzeugt pro Lauf zwei Dateien:
  1. JSONL-Log (maschinenlesbar): python_run_<ts>.jsonl
  2. Masterlog (C++-kompatibel): simulation_log_<ts>.txt

Verwendung:
    logger = RunLogger()
    logger.set_metadata(backend="python", package="iodine_lafleur_v1", ...)
    logger.add_point(q0=0.5, P=14.3, I_beam=12.05, status="converged", ...)
    logger.add_fail(q0=0.3, status="above_P_max", diag="...")
    logger.finalize(elapsed=12.3, count_ok=5, count_fail=2)
"""
from __future__ import annotations
import json, time
from pathlib import Path
from dataclasses import dataclass, field, asdict


def _ts():
    return time.strftime("%Y%m%d_%H%M%S")


def _ts_readable():
    return time.strftime("%Y-%m-%d %H:%M:%S")


@dataclass
class SweepPoint:
    """Ein Sweep-Punkt (Erfolg oder Fehler)."""
    idx: int = 0
    Q0_sccm: float = 0.0
    status: str = ""              # converged, above_P_max, below_P_min, ...
    P_W: float = 0.0
    I_beam_mA: float = 0.0
    Te_eV: float = 0.0
    Tg_K: float = 0.0
    ne: float = 0.0
    ng: float = 0.0
    # I-fix spezifisch
    I_target_mA: float = 0.0
    delta_I_mA: float = 0.0
    # Iod-spezifisch (optional)
    diss: float = 0.0
    fIp: float = 0.0
    fI2p: float = 0.0
    alpha: float = 0.0
    nI: float = 0.0
    nI2: float = 0.0
    nIp: float = 0.0
    nI2p: float = 0.0
    nIm: float = 0.0
    # Diagnostik
    converged: bool = False
    merit: float = 0.0
    iterations: int = 0
    diag: str = ""


class RunLogger:
    """Sammelt Sweep-Punkte und schreibt strukturierte Logs."""

    def __init__(self, output_dir: str | Path = "."):
        self._dir = Path(output_dir)
        self._ts = _ts()
        self._meta: dict = {}
        self._params: dict = {}
        self._points: list[SweepPoint] = []
        self._events: list[dict] = []

    # ── Metadaten ────────────────────────────────────────────

    def set_metadata(self, **kwargs):
        """Setze Laufmetadaten (backend, package, mode, ...)."""
        self._meta.update(kwargs)
        self._meta.setdefault("timestamp_start", _ts_readable())
        self._meta.setdefault("backend", "python")

    def set_params(self, params: dict):
        """Setze die Eingabeparameter (Geometrie, Sweep, etc.)."""
        self._params = dict(params)

    # ── Sweep-Punkte ─────────────────────────────────────────

    def add_point(self, **kwargs) -> SweepPoint:
        """Fuge einen Sweep-Punkt hinzu (Erfolg oder Fehler)."""
        pt = SweepPoint(**kwargs)
        self._points.append(pt)
        return pt

    def add_event(self, idx: int, q0: float, message: str):
        """Freiform-Event (Warnung, Diagnose, etc.)."""
        self._events.append({"idx": idx, "Q0_sccm": q0, "message": message})

    # ── Finalisierung ────────────────────────────────────────

    def finalize(self, elapsed: float = 0.0,
                 count_ok: int = 0, count_fail: int = 0):
        """Schreibe JSONL und Masterlog."""
        self._meta["timestamp_end"] = _ts_readable()
        self._meta["runtime_seconds"] = round(elapsed, 3)

        summary = {
            "total_points": len(self._points),
            "converged": count_ok,
            "failed": count_fail,
            "runtime_seconds": round(elapsed, 3),
        }

        # 1. JSONL-Log
        jsonl_path = self._dir / f"python_run_{self._ts}.jsonl"
        with open(jsonl_path, "w", encoding="utf-8") as f:
            f.write(json.dumps({"type": "metadata", **self._meta}) + "\n")
            if self._params:
                f.write(json.dumps({"type": "params", **self._params}) + "\n")
            for pt in self._points:
                f.write(json.dumps({"type": "point", **asdict(pt)}) + "\n")
            for ev in self._events:
                f.write(json.dumps({"type": "event", **ev}) + "\n")
            f.write(json.dumps({"type": "summary", **summary}) + "\n")

        # 2. C++-kompatibles Masterlog
        txt_path = self._dir / f"simulation_log_{self._ts}.txt"
        self._write_masterlog(txt_path, summary)

        return str(jsonl_path), str(txt_path)

    # ── Masterlog (C++-kompatibles Textformat) ───────────────

    def _write_masterlog(self, path: Path, summary: dict):
        with open(path, "w", encoding="utf-8") as f:
            # Header
            f.write("==================================================\n")
            f.write("GLOBAL PLASMA MODEL - SIMULATION LOG\n")
            for k, v in self._meta.items():
                f.write(f"{k + ':':<20} {v}\n")
            f.write("==================================================\n\n")

            # Parameters
            f.write("--------------------------------------------------\n")
            f.write("SIMULATION PARAMETERS\n")
            f.write("--------------------------------------------------\n")
            for k, v in self._params.items():
                f.write(f"{k:<35}| {v}\n")
            f.write("\n")

            # Result Table
            f.write("--------------------------------------------------\n")
            f.write("RESULT TABLE\n")
            f.write("--------------------------------------------------\n")
            hdr = (f"{'idx':<5}| {'Q0_sccm':<9}| {'status':<18}| {'P_W':>10}| "
                   f"{'I_mA':>8}| {'Te_eV':>7}| {'Tg_K':>7}| "
                   f"{'ne':>12}| {'ng':>12}| {'note'}\n")
            f.write(hdr)
            f.write("-" * 110 + "\n")
            for pt in self._points:
                s = pt.status[:18]
                if pt.converged:
                    f.write(f"{pt.idx:<5}| {pt.Q0_sccm:<9.4f}| {s:<18}| "
                            f"{pt.P_W:>10.2f}| {pt.I_beam_mA:>8.3f}| "
                            f"{pt.Te_eV:>7.3f}| {pt.Tg_K:>7.1f}| "
                            f"{pt.ne:>12.3e}| {pt.ng:>12.3e}| {pt.diag}\n")
                else:
                    f.write(f"{pt.idx:<5}| {pt.Q0_sccm:<9.4f}| {s:<18}| "
                            f"{pt.P_W:>10.1f}| {pt.I_beam_mA:>8.2f}| "
                            f"{'--':>7}| {'--':>7}| {'--':>12}| {'--':>12}| "
                            f"{pt.diag[:60]}\n")
            f.write("\n")

            # Events
            if self._events:
                f.write("--------------------------------------------------\n")
                f.write("EVENT DETAILS\n")
                f.write("--------------------------------------------------\n")
                for ev in self._events:
                    f.write(f"  [{ev['idx']}] Q0={ev['Q0_sccm']:.4f}: {ev['message']}\n")
                f.write("\n")

            # Summary
            f.write("--------------------------------------------------\n")
            f.write("SUMMARY\n")
            f.write("--------------------------------------------------\n")
            for k, v in summary.items():
                f.write(f"{k:<30}| {v}\n")
            f.write("--------------------------------------------------\n\n")

            # Machine-readable data (C++-kompatibel)
            f.write("--------------------------------------------------\n")
            f.write("MACHINE-READABLE DATA TABLE\n")
            f.write("--------------------------------------------------\n")
            cols = ("idx|Q0sccm|status|P_W|I_mA|Te_eV|Tg_K|ne|ng"
                    "|I_target_mA|delta_I_mA"
                    "|diss|fIp|fI2p|alpha|merit|iterations|note")
            f.write(f"DATA_HEADER|{cols}\n")
            for pt in self._points:
                vals = (f"{pt.idx}|{pt.Q0_sccm:.4f}|{pt.status}|"
                        f"{pt.P_W:.2f}|{pt.I_beam_mA:.3f}|"
                        f"{pt.Te_eV:.3f}|{pt.Tg_K:.1f}|"
                        f"{pt.ne:.6e}|{pt.ng:.6e}|"
                        f"{pt.I_target_mA:.2f}|{pt.delta_I_mA:.2f}|"
                        f"{pt.diss:.4f}|{pt.fIp:.4f}|{pt.fI2p:.4f}|"
                        f"{pt.alpha:.4f}|{pt.merit:.2e}|{pt.iterations}|"
                        f"{pt.diag}")
                f.write(f"DATA|{vals}\n")
            f.write("--------------------------------------------------\n")

    # ── Hilfsmethoden ────────────────────────────────────────

    @property
    def point_count(self):
        return len(self._points)

    @property
    def jsonl_filename(self):
        return f"python_run_{self._ts}.jsonl"
