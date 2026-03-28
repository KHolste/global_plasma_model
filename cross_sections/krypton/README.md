# Krypton Cross Sections

Placeholder for Krypton (Kr) cross-section data.

To populate, download Krypton data from LXCat (Biagi database) and run:

    python convert_lxcat.py <lxcat_file> krypton

Then precompute rate tables:

    python precompute_rates.py  # (requires adaptation for krypton)

Gas properties:
- Atomic mass: 83.798 u (1.3915e-25 kg)
- Ionization energy: 13.9996 eV
- First excitation energy: ~10.0 eV
