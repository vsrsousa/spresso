#!/usr/bin/env python3
"""
ASE Database Usage Examples

This script demonstrates how to use ASE databases to store and retrieve atomic structures.
ASE (Atomic Simulation Environment) provides a convenient way to store atomic structures
with metadata in SQLite databases.
"""

import os
from ase.db import connect as ase_db_connect
from ase import Atoms

def main():
    """Demonstrate ASE database operations."""

    # Default database path (same as used in the GUI)
    db_path = os.path.expanduser("~/.xespresso/structures.db")

    print("ASE Database Usage Examples")
    print("=" * 40)
    print(f"Database path: {db_path}")
    print()

    # 1. Connect to the database
    print("1. Connecting to database...")
    db = ase_db_connect(db_path)
    print("   ✓ Connected successfully")
    print()

    # 2. List all structures
    print("2. Listing all structures:")
    rows = list(db.select())
    print(f"   Found {len(rows)} structure(s)")
    print()

    if rows:
        print("   Available structures:")
        print("   ID | Formula | Atoms | Energy | Tags")
        print("   " + "-" * 50)

        for row in rows:
            # Get energy if available
            energy = f"{row.energy:.3f}" if hasattr(row, 'energy') and row.energy is not None else "N/A"

            # Get tags (key-value pairs)
            tags = ", ".join(row.key_value_pairs.keys()) if row.key_value_pairs else ""

            print(f"   {row.id:2d} | {row.formula:8s} | {row.natoms:3d} | {energy:>8s} | {tags}")

        print()

        # 3. Query specific structures
        print("3. Query examples:")

        # Find structures with specific elements
        print("   Structures containing Gd:")
        gd_rows = list(db.select("Gd"))
        for row in gd_rows:
            print(f"      ID {row.id}: {row.formula}")

        # Find structures by tag
        print("   Structures tagged as 'bulk':")
        bulk_rows = list(db.select(bulk=True))
        for row in bulk_rows:
            print(f"      ID {row.id}: {row.formula}")

        # Find structures by formula pattern (using manual filtering)
        print("   Structures with 'Gd' in formula:")
        gd_in_formula = [row for row in rows if 'Gd' in row.formula]
        for row in gd_in_formula:
            print(f"      ID {row.id}: {row.formula}")

        print()

        # 4. Load a specific structure
        print("4. Loading structure ID 1:")
        try:
            row = db.get(id=1)
            atoms = row.toatoms()
            print(f"   Formula: {atoms.get_chemical_formula()}")
            print(f"   Number of atoms: {len(atoms)}")
            print(f"   Chemical symbols: {atoms.get_chemical_symbols()}")
            print(f"   Cell: {atoms.cell}")
            print(f"   Positions shape: {atoms.positions.shape}")

            # Show metadata
            if row.key_value_pairs:
                print(f"   Metadata: {row.key_value_pairs}")

        except Exception as e:
            print(f"   Error loading structure: {e}")

        print()

        # 5. Advanced queries
        print("5. Advanced query examples:")

        # Count structures by element
        print("   Element counts across all structures:")
        element_counts = {}
        for row in rows:
            symbols = row.toatoms().get_chemical_symbols()
            for symbol in set(symbols):  # unique elements per structure
                element_counts[symbol] = element_counts.get(symbol, 0) + 1

        for element, count in sorted(element_counts.items()):
            print(f"      {element}: {count} structure(s)")

        print()

        # Find largest structures
        print("   Largest structures (by atom count):")
        sorted_rows = sorted(rows, key=lambda r: r.natoms, reverse=True)
        for row in sorted_rows[:3]:  # top 3
            print(f"      {row.formula} ({row.natoms} atoms)")

    else:
        print("   Database is empty. Here's how to add structures:")
        print()
        print("   # Example: Add a simple structure")
        print("   from ase import Atoms")
        print("   from ase.db import connect")
        print("   ")
        print("   # Create atoms object")
        print("   atoms = Atoms('H2O', positions=[[0,0,0], [0.758,0,0], [0.379,0.652,0]])")
        print("   ")
        print("   # Connect to database")
        print("   db = connect('~/.xespresso/structures.db')")
        print("   ")
        print("   # Write with metadata")
        print("   db.write(atoms, name='water', molecule=True, optimized=False)")
        print()

    print("Database operations completed!")

if __name__ == "__main__":
    main()