# Argon Cross Sections

Placeholder for Argon (Ar) cross-section data.

To populate, download Argon data from LXCat (Biagi database) and run:

    python convert_lxcat.py <lxcat_file> argon

Then precompute rate tables:

    python precompute_rates.py  # (requires adaptation for argon)

Gas properties:
- Atomic mass: 39.948 u (6.6335e-26 kg)
- Ionization energy: 15.7596 eV
- First excitation energy: ~11.55 eV
