#!/usr/bin/env python3
"""
Initialize DatabaseWorkflow and verify databases are created
"""

from pathlib import Path
import os

print("="*70)
print("DATABASE INITIALIZATION")
print("="*70)

# Check current directory
cwd = os.getcwd()
print(f"\nCurrent directory: {cwd}")

# Import and initialize
print("\nImporting DatabaseWorkflow...")
from xespresso.db import DatabaseWorkflow

print("Creating DatabaseWorkflow instance...")
wf = DatabaseWorkflow()

print("✅ DatabaseWorkflow initialized")
print(f"   Structures DB: ~/.xespresso/structures.db")
print(f"   Provenance DB: ~/.xespresso/provenance.db")

# Check if files were created
print("\nChecking database files...")
for db_file in ['~/.xespresso/structures.db', '~/.xespresso/provenance.db']:
    path = Path(db_file).expanduser()
    if path.exists():
        size = path.stat().st_size
        print(f"✅ {db_file} exists ({size} bytes)")
    else:
        print(f"❌ {db_file} NOT FOUND")

# Try a test write
print("\nTesting database write...")
try:
    from ase.build import bulk
    si = bulk('Si', 'diamond', a=5.43)
    id = wf.db.write(si)
    print(f"✅ Successfully wrote test structure with ID: {id}")
    
    # Check again
    print("\nChecking database files again...")
    for db_file in ['~/.xespresso/structures.db', '~/.xespresso/provenance.db']:
        path = Path(db_file).expanduser()
        if path.exists():
            size = path.stat().st_size
            print(f"✅ {db_file} exists ({size} bytes)")
except Exception as e:
    print(f"❌ Error: {e}")

print("\n" + "="*70)
print("Database initialization complete!")
print("="*70)
