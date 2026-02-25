#!/usr/bin/env python3
"""
Simple ASE Database Structure Lister

This script shows how to connect to an ASE database and list all available structures.
"""

import os
from ase.db import connect as ase_db_connect

def list_structures(db_path=None):
    """List all structures in the ASE database."""

    if db_path is None:
        db_path = os.path.expanduser("~/.xespresso/structures.db")

    print(f"Connecting to database: {db_path}")

    if not os.path.exists(db_path):
        print("Database file does not exist.")
        return

    try:
        db = ase_db_connect(db_path)
        rows = list(db.select())

        print(f"\nFound {len(rows)} structure(s) in database:\n")

        if rows:
            print("ID  | Formula  | Atoms | Tags")
            print("-" * 50)

            for row in rows:
                tags = ", ".join(row.key_value_pairs.keys()) if row.key_value_pairs else ""
                print(f"{row.id:3d} | {row.formula:9s} | {row.natoms:5d} | {tags}")

        else:
            print("Database is empty.")

    except Exception as e:
        print(f"Error accessing database: {e}")

if __name__ == "__main__":
    list_structures()