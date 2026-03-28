# Cross Section Data for Global Xenon Plasma Model

## Directory Structure

```
cross_sections/
  xenon/
    ionization.csv      # Xe + e -> Xe+ + 2e
    excitation.csv      # Xe + e -> Xe* + e (placeholder)
    elastic.csv         # Xe + e -> Xe + e (momentum transfer)
    gas_properties.csv  # Mass, threshold energies, etc.
  krypton/              # (prepared, data not yet populated)
  argon/                # (prepared, data not yet populated)
```

## CSV Format

Each cross section file uses a simple CSV format:

```csv
# Comment lines start with #
# Threshold: 12.127 eV
# Process: Xe + e -> Xe+ + 2e
energy_eV,cross_section_m2
12.1298,3.665e-21
13.0,1.547e-20
...
```

- Column 1: `energy_eV` – Electron energy in electronvolts
- Column 2: `cross_section_m2` – Cross section in m^2
- Lines starting with `#` are comments (metadata)
- Header line `energy_eV,cross_section_m2` is skipped

## Rate Coefficient Computation

Rate coefficients K(Te) are computed by Maxwell-Boltzmann integration:

```
K(Te) = integral_0^inf sigma(E) * v(E) * f_MB(E, Te) dE
```

See `rate_coefficients.py` for the implementation.

## Usage

```python
from rate_coefficients import GasRateModel

model = GasRateModel.load("cross_sections/xenon")
Kiz = model.Kiz(3.5)  # Te = 3.5 eV -> Kiz in m^3/s
Kel = model.Kel(3.5)  # -> Kel in m^3/s
```

## Configuration

In params.txt or config files:
```
rate_model legacy      # Use Chabert polynomial fits (default)
rate_model tabulated   # Use cross section tables from files
gas_species xenon      # Which gas directory to load
```

## Adding New Gases

1. Create a subdirectory: `cross_sections/krypton/`
2. Add CSV files: `ionization.csv`, `elastic.csv`, `excitation.csv`
3. Add `gas_properties.csv` with mass, threshold energies
4. Set `gas_species krypton` in config
