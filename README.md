# Global Plasma Model

A stationary 0D global discharge model for inductively coupled RF ion thrusters.

Based on **Chabert et al., Phys. Plasmas 19, 073512 (2012)**, extended with tabulated cross-section data from the **Biagi/LXCat database (MagBoltz V8.97)**.

## Features

- **4-equation 0D model**: ion density, neutral density, electron temperature, gas temperature
- **3 rate model presets**: Legacy (paper-compatible), Conservative tabulated, Full tabulated
- **2 operating modes**: fixed beam current (I-fix) and self-consistent (SC)
- **Multi-gas architecture**: Xenon fully integrated, Krypton/Argon structurally prepared
- **PyQt6 GUI** with real-time streaming, live plots, and parameter sweeps
- **Levenberg-Marquardt solver** with Pseudo-Transient Continuation (PTC) warm-start
- **Adaptive power bracketing** for robust fixed-current operation across all rate models

## Model Overview

The solver finds the stationary equilibrium of four coupled balance equations:

| Equation | Balance | Primary unknown |
|----------|---------|-----------------|
| r1 | Ion production = wall loss | Te |
| r2 | Gas inflow = ionization + outflow | ng |
| r3 | RF absorption = ionization + excitation + elastic + wall losses | n (via P_RFG) |
| r4 | Gas heating = thermal conduction | Tg |

RF coupling is modeled via the complex plasma permittivity with cylindrical Bessel functions for the ICP geometry. The collision frequency nu_m = Kel(Te) * ng is the central coupling parameter between plasma state and RF absorption.

## Rate Models

| Preset | rate_model | Kiz | Kex | Kel |
|--------|-----------|-----|-----|-----|
| Legacy | 0 | Chabert polynomial | Chabert Arrhenius | Constant 1e-13 m3/s |
| Conservative | 1 | Biagi tabulated | Biagi tabulated | Constant 1e-13 m3/s |
| Full tabulated | 2 | Biagi tabulated | Biagi tabulated | Biagi tabulated |

**Key finding**: The Biagi elastic momentum transfer cross-section yields Kel values 2-3x higher than the Legacy constant in the operating range (Te = 2-8 eV). This significantly affects RF coupling and shifts Te by ~1.2 eV upward in self-consistent mode.

## Supported Gases

| Gas | Status | Cross-sections |
|-----|--------|---------------|
| Xenon (Xe) | Fully integrated | Biagi/LXCat (elastic, ionization, 50 excitation processes) |
| Krypton (Kr) | Structurally prepared | Physical constants set, cross-section data needed |
| Argon (Ar) | Structurally prepared | Physical constants set, cross-section data needed |

## Project Structure

```
global-xenon-model/
  chabert_modified.cpp      C++ solver (main program, ~2400 lines)
  bessel_wrapper.cpp/hpp    Bessel functions (separate compilation unit)
  bessel-library.hpp        Template Bessel library
  gui.py                    PyQt6 GUI with real-time streaming
  style.qss                 Dark theme stylesheet
  log_viewer.py             Standalone log viewer with pyqtgraph
  rate_coefficients.py      Maxwell-Boltzmann integration, legacy fits
  convert_lxcat.py          LXCat raw data to CSV converter
  precompute_rates.py       Generate lookup tables from cross-sections
  test_solver.py            Test suite (5 tests)
  dietz_validation.py       RIT-4 benchmark validation (Dietz 2021)
  config_dietz_rit4.txt     RIT-4 configuration file
  param_study.py            P x Q0 parameter study
  rate_model_comparison.py  Compare all 3 rate model presets
  audit_kel.py              Kel implementation audit
  cross_sections/
    xenon/                  Complete Biagi cross-section data
      elastic.csv           199 data points (momentum transfer)
      ionization.csv        54 data points
      excitation/           50 individual process files
      kel_table.csv         Pre-integrated Kel(Te), 391 points
      kiz_table.csv         Pre-integrated Kiz(Te), 391 points
      kex_table.csv         Pre-integrated Kex(Te), 391 points
      metadata.json         Process overview
    krypton/                Placeholder (README only)
    argon/                  Placeholder (README only)
    tests/
      xenon_all.txt         LXCat source file (Biagi V8.97)
  docs/
    Global_Plasma_Model_Documentation_Current.docx
    validation_report_rit4_kk.docx
    validation_summary_global_xenon_model.docx
```

## Prerequisites

- **C++ compiler**: g++ with C++17 support
- **Python 3.10+** with packages: `PyQt6`, `pyqtgraph`, `numpy`, `scipy`
- **Windows**: WSL recommended for compilation (Linux g++ via WSL)
- **Linux/macOS**: native g++ works directly

```bash
pip install PyQt6 pyqtgraph numpy scipy python-docx
```

## Build

```bash
# Compile Bessel wrapper (separate unit for faster rebuilds)
g++ -O3 -march=native -std=c++17 -c bessel_wrapper.cpp -o bessel_wrapper.o

# Compile solver
g++ -O3 -march=native -std=c++17 chabert_modified.cpp bessel_wrapper.o -o chabert
```

On Windows with WSL:
```bash
wsl g++ -O3 -march=native -std=c++17 -c bessel_wrapper.cpp -o bessel_wrapper.o
wsl g++ -O3 -march=native -std=c++17 chabert_modified.cpp bessel_wrapper.o -o chabert
```

## Run

### GUI
```bash
python gui.py
```
The GUI provides compilation, parameter editing, sweep configuration, and real-time visualization.

### Command line
```bash
./chabert params.txt
```

### Example configuration (params.txt)
```
R 0.02
L 0.04
betai 0.5
betag 0.05145
frequency 2.5e6
Nw 6.0
R_ohm 0.36
Rc 0.02
lc 0.04
Vgrid 1500.0
sgrid 0.001
P_RFG 45.0
P_RFG_max 200.0
Q0sccm_start 0.30
Q0sccm_step 0.01
jjmax 50
I_soll 15.0
solve_mode 2
gas_species xenon
rate_model 2
use_paper_kel 1
```

## Example Workflow

1. Start the GUI: `python gui.py`
2. Select gas (Xenon), rate model (Full tabulated), and operating mode
3. Set thruster geometry and operating parameters
4. Click "Kompilieren" to build the solver
5. Configure sweep range (Q0 start, step, number of steps)
6. Click "Simulation starten" to run
7. Watch real-time plots update as each Q0 point converges
8. Open log viewer for detailed post-analysis

## Cross-Section Data

The `cross_sections/xenon/` directory contains Biagi/LXCat data for Xenon. To add a new gas:

1. Download cross-section data from [LXCat](https://www.lxcat.net/) (Biagi database)
2. Convert: `python convert_lxcat.py <lxcat_file> <gas_name>`
3. Precompute rate tables: adapt `precompute_rates.py` for the new gas
4. The solver automatically resolves paths via `cross_sections/<gas_species>/`

## Validation

- **Chabert 2012**: Legacy mode reproduces paper results
- **Dietz et al. 2021**: RIT-4 benchmark comparison available via `dietz_validation.py`
- **Test suite**: `python test_solver.py` runs 5 automated tests (reference, convergence, unphysical params, regression, self-consistent mode)

## Known Limitations

- 0D assumption: homogeneous density profile (correctable via `density_profile_factor`)
- No doubly-charged ions (Xe++)
- No separate metastable state balance
- Solver convergence limit at very high power (>60-70 W for Full tabulated Xe)
- Krypton/Argon: physical constants available but no tabulated cross-sections yet
- Legacy rate fits (Kiz, Kex) are calibrated for Xenon only

## References

- P. Chabert, J. Arancibia Monreal, J. Bredin, L. Popelier, A. Aanesland, *Phys. Plasmas* **19**, 073512 (2012)
- T. Dietz, S. Seeger, C. Volkmar, U. Scholz, *Proc. IEPC* (2021)
- S. F. Biagi, MagBoltz, via [LXCat](https://www.lxcat.net/)
- M. A. Lieberman, A. J. Lichtenberg, *Principles of Plasma Discharges and Materials Processing*, Wiley (2005)
