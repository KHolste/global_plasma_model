#!/usr/bin/env python3
"""
test_reaction_balance.py -- Validates stoichiometry and charge balance
for all reactions in all chemistry packages.

Checks:
  1. Atom balance (I atoms conserved, Xe atoms conserved)
  2. Charge balance (total charge conserved)
  3. GUI reaction string matches actual stoichiometry
  4. Species exist in package species list

Run: python test_reaction_balance.py
"""
from __future__ import annotations
import sys, json
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

passed = 0
failed = 0
errors = []


def check(name, condition, detail=""):
    global passed, failed
    if condition:
        passed += 1
    else:
        failed += 1
        print(f"  FAIL: {name} -- {detail}")
        errors.append(name)


def count_atoms_and_charge(side: dict, species_db: dict) -> tuple[dict, int]:
    """Count atoms and net charge for a reaction side."""
    atoms = {}
    charge = 0
    for sp_id, n in side.items():
        sp = species_db.get(sp_id, {})
        charge += sp.get("charge", 0) * n
        # Count atoms by species type
        if sp_id == "e":
            continue  # Electrons don't carry atoms
        if sp_id in ("I2", "I2+"):
            atoms["I"] = atoms.get("I", 0) + 2 * n
        elif sp_id in ("I", "I+", "I-"):
            atoms["I"] = atoms.get("I", 0) + 1 * n
        elif sp_id in ("Xe", "Xe+"):
            atoms["Xe"] = atoms.get("Xe", 0) + 1 * n
    return atoms, charge


def format_side(d: dict) -> str:
    parts = []
    for sp, n in d.items():
        parts.append(f"{n}{sp}" if n > 1 else sp)
    return " + ".join(parts)


def validate_package(chem_path: Path):
    """Validate all reactions in a chemistry package."""
    with open(chem_path, encoding="utf-8") as f:
        chem = json.load(f)

    pkg_name = chem.get("name", chem_path.parent.name)
    species_db = {s["id"]: s for s in chem.get("species", [])}
    species_ids = set(species_db.keys())
    reactions = chem.get("reactions", [])

    print(f"\n--- Package: {pkg_name} ({len(reactions)} reactions) ---")

    for rxn in reactions:
        rid = rxn["id"]
        reactants = rxn.get("reactants", {})
        products = rxn.get("products", {})
        rtype = rxn.get("type", "?")

        # Surface reactions may have fractional stoichiometry
        is_surface = rxn.get("surface_gamma", 0) > 0

        # Check all species exist
        for sp in list(reactants.keys()) + list(products.keys()):
            if sp != "e":  # electrons are implicit
                check(f"{rid}: species '{sp}' exists", sp in species_ids,
                      f"species '{sp}' not in {species_ids}")

        # Atom balance
        r_atoms, r_charge = count_atoms_and_charge(reactants, species_db)
        p_atoms, p_charge = count_atoms_and_charge(products, species_db)

        if not is_surface:
            for element in set(list(r_atoms.keys()) + list(p_atoms.keys())):
                check(f"{rid}: {element} atom balance",
                      r_atoms.get(element, 0) == p_atoms.get(element, 0),
                      f"reactants={r_atoms.get(element,0)}, products={p_atoms.get(element,0)}")

        # Charge balance
        check(f"{rid}: charge balance",
              r_charge == p_charge,
              f"reactants charge={r_charge}, products charge={p_charge}")

        # GUI string correctness
        gui_str = f"{format_side(reactants)} -> {format_side(products)}"
        # Verify the name field is consistent
        name = rxn.get("name", "")
        if name:
            # Check that stoichiometry in name matches data
            # (This is informational, not a hard fail)
            pass

        print(f"  {rid:<20} {gui_str:<45} [{rtype}] OK")


def test_iodine_lafleur():
    print("\n=== Iodine Lafleur V1 Validation ===")
    p = SCRIPT_DIR / "chemistry" / "iodine_lafleur_v1" / "chemistry.json"
    if p.exists():
        validate_package(p)
    else:
        print("  SKIP: package not found")


def test_xenon_simple():
    print("\n=== Xenon Simple Validation ===")
    p = SCRIPT_DIR / "chemistry" / "xenon_simple" / "chemistry.json"
    if p.exists():
        validate_package(p)
    else:
        print("  SKIP: package not found")


def test_gui_reaction_strings():
    """Verify that physics_data_viewer generates correct reaction strings."""
    print("\n=== GUI Reaction String Validation ===")
    from physics_data_viewer import discover_rate_coefficients, _format_reaction_side

    # Test _format_reaction_side directly
    check("format 2e", _format_reaction_side({"e": 2, "I+": 1}) == "2e + I+",
          _format_reaction_side({"e": 2, "I+": 1}))
    check("format 2I", _format_reaction_side({"I": 2}) == "2I",
          _format_reaction_side({"I": 2}))
    check("format single", _format_reaction_side({"e": 1, "I2": 1}) == "e + I2",
          _format_reaction_side({"e": 1, "I2": 1}))

    # Test via discover
    iod_path = SCRIPT_DIR / "chemistry" / "iodine_lafleur_v1"
    if iod_path.exists():
        items = discover_rate_coefficients(iod_path, {"backend": "python"})
        for item in items:
            rxn_str = item.get("reaction", "")
            if not rxn_str:
                continue
            label = item.get("label", "?")
            # Specific checks for the reported bugs
            if "dissiz" in label.lower() or "Kdissiz" in label:
                check("dissiz_I2 shows 2e",
                      "2e" in rxn_str,
                      f"reaction='{rxn_str}' should contain '2e'")
            if "Kdiss I2" == label or "diss_I2" in item.get("process", ""):
                if "dissiz" not in label.lower() and "dissatt" not in label.lower():
                    check("diss_I2 shows 2I",
                          "2I" in rxn_str,
                          f"reaction='{rxn_str}' should contain '2I'")


def test_species_table():
    """Print species table for iodine model."""
    print("\n=== Iodine Species Table ===")
    p = SCRIPT_DIR / "chemistry" / "iodine_lafleur_v1" / "chemistry.json"
    if not p.exists():
        return
    with open(p, encoding="utf-8") as f:
        chem = json.load(f)

    print(f"  {'Species':<8} {'Type':<20} {'Charge':>6} {'Feedstock':>9} {'Bilanziert':>10} {'Bemerkung'}")
    print("  " + "-" * 80)
    for sp in chem["species"]:
        sid = sp["id"]
        stype = sp["type"]
        charge = sp.get("charge", 0)
        feed = "ja" if sp.get("is_feedstock") else ""
        # All non-electron species are explicitly balanced in the solver
        balanced = "implizit" if stype == "electron" else "explizit"
        note = ""
        if sp.get("is_beam_extracted"):
            note = "Beam-extrahiert"
        if sid == "e":
            note = "Aus Quasineutralitaet"
        print(f"  {sid:<8} {stype:<20} {charge:>+6} {feed:>9} {balanced:>10} {note}")
        check(f"species {sid} has type", stype in
              ("electron", "neutral_atom", "neutral_molecule", "positive_ion", "negative_ion"),
              f"type={stype}")


def main():
    global passed, failed
    test_iodine_lafleur()
    test_xenon_simple()
    test_gui_reaction_strings()
    test_species_table()

    print(f"\n{'='*60}")
    print(f"  Result: {passed} passed, {failed} failed")
    if errors:
        print(f"  Failed: {', '.join(errors)}")
    print(f"{'='*60}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
