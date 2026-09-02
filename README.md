# Global Plasma Model

A stationary 0D global discharge model for inductively coupled RF ion thrusters.

Based on **Chabert et al., Phys. Plasmas 19, 073512 (2012)**, extended with tabulated cross-section data from the **Biagi/LXCat database (MagBoltz V8.97)**.

## Features

- **4-equation 0D model**: ion density, neutral density, electron temperature, gas temperature
- **3 rate model presets**: Legacy (paper-compatible), Conservative tabulated, Full tabulated
- **2 operating modes**: fixed beam current (I-fix) and self-consistent (SC)
- **Multi-gas architecture**: Xenon fully integrated, Krypton/Argon structurally prepared
- **Chemistry packages**: species/reactions/rate tables read from JSON by both the Python and the C++ side; select one with `chemistry_package` to run the generic solver over it
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
  C++ core (seven translation units + entry point)
    main.cpp                Entry point, sweep loop
    sim_context.hpp         Single mutable context, replaces all globals
    sim_config.cpp/hpp      Config file parsing
    phys_const.hpp          Physical constants
    rates.cpp/hpp           Rate coefficients, table lookup
    physics.cpp/hpp         RF coupling, residuals, derived quantities
    solver.cpp/hpp          LM, PTC, multi-start, I-fix, SC mode
    sim_logging.cpp/hpp     Run logging
    chem_loader.cpp/hpp     Loads chemistry packages (JSON + rate tables)
    beam_extraction_cpp.hpp Bohm / Child-Langmuir extraction model
    bessel_wrapper.cpp/hpp  Bessel functions (separate compilation unit)
    bessel-library.hpp      Template Bessel library
  build.py                  Build definition, used by CLI and GUI alike
  run_tests.py              Runs the whole test suite, reports pass/fail
  gui.py                    PyQt6 GUI with real-time streaming
  style.qss                 Dark theme stylesheet
  log_viewer.py             Standalone log viewer with pyqtgraph
  physics_data_viewer.py    Cross sections and rate coefficients inspector
  run_config.py             Single source of truth for run parameters
  package_registry.py       Discovers chemistry packages, routes backends
  Generic chemistry layer (not yet wired into the production path)
    plasma_chemistry.py     Species, reactions, balance assembler
    generic_solver.py       LM solver over the dynamic state vector
    chem_system.hpp         C++ counterpart
    generic_lm.hpp          C++ generic LM solver
  third_party/
    json.hpp                nlohmann/json 3.11.3 (MIT), used only by chem_loader
  rate_coefficients.py      Maxwell-Boltzmann integration, legacy fits
  convert_lxcat.py          LXCat raw data to CSV converter
  precompute_rates.py       Generate lookup tables from cross-sections
  dietz_validation.py       RIT-4 benchmark validation (Dietz 2021)
  config_dietz_rit4.txt     RIT-4 configuration file
  param_study.py            P x Q0 parameter study
  rate_model_comparison.py  Compare all 3 rate model presets
  audit_kel.py              Kel implementation audit
  cross_sections/
    xenon/
      biagi/                Biagi/LXCat data (elastic, ionization, 50 excitations)
        elastic.csv         Momentum transfer
        ionization.csv
        excitation/         Individual process files
        kel_table.csv       Pre-integrated Kel(Te)
        kiz_table.csv       Pre-integrated Kiz(Te)
        kex_table.csv       Pre-integrated Kex(Te)
        db_info.json        Database identification
        metadata.json       Process overview
      hayashi/              Second database, same layout
    krypton/                Placeholder (README only)
    argon/                  Placeholder (README only)
    tests/                  LXCat source files
  chemistry/                Chemistry packages (xenon_simple, two iodine sets)
  Iodmodell/                Legacy iodine model (2019), reference only
  archive/                  Pre-split monolith, for reference
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
python build.py            # compile the C++ core
python build.py --clean    # remove objects and binary first
python build.py --tests    # also build the C++ test programs
```

`build.py` holds the single definition of which translation units are built
with which flags. The GUI's "Kompilieren" button uses the same definition, so
there is no second place that can drift. On Windows the GUI prefers WSL when
available and falls back to a native g++.

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

Add `chemistry_package xenon_simple` (or any directory under `chemistry/`) to
run the generic solver over that package instead of the hard-wired physics.
If the package cannot be loaded, the run falls back to the hard-wired path and
says so on stderr.

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
- **Test suite**: `python run_tests.py` runs the whole test bestand — 19 Python
  tests plus the three compiled C++ test programs — and reports pass/fail per test.

```bash
python run_tests.py              # everything
python run_tests.py --build      # compile the core first, then test
python run_tests.py --only ifix  # only tests whose name contains "ifix"
python run_tests.py --list       # just list what would run
```

The runner uses the interpreter it was started with and names it in the header.
Several tests import PyQt6; if the interpreter lacks it, those tests are
reported as an incomplete environment rather than as a failure. Last full run:
22 of 22 pass.

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
