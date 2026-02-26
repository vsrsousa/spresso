#!/usr/bin/env python3
"""
Silicon SCF calculation with caching and provenance tracking
"""

from ase.build import bulk
from xespresso.db import DatabaseWorkflow
import logging
import os
from pathlib import Path

# Clean old databases (optional - remove if you want to keep history)
for db_file in ['~/.xespresso/structures.db', '~/.xespresso/provenance.db']:
    db_path = Path(db_file).expanduser()
    if db_path.exists():
        db_path.unlink()
        print(f"Cleaned {db_file}")

# Reduce logging noise (optional)
logging.getLogger('xespresso').setLevel(logging.WARNING)

# Create silicon structure
si = bulk('Si', "diamond", a=5.43)

# Initialize database workflow (uses defaults: ~/.xespresso/structures.db and ~/.xespresso/provenance.db)
wf = DatabaseWorkflow()

print("="*70)
print("SILICON SCF CALCULATION WITH CACHING")
print("="*70)

# STEP 1: Store input structure and get its ID
print("\nSTEP 1: Store input structure")
print("-"*70)
input_structure_id = wf.db.write(si)
print(f"✅ Structure stored with ID: {input_structure_id}")
# Metadados são opcionais! Você também pode fazer:
# input_structure_id = wf.db.write(si, structure_name='Si_diamond', lattice_parameter=5.43)

# STEP 2: Run SCF calculation (will be cached)
print("\nSTEP 2: First SCF calculation")
print("-"*70)
print("Running SCF calculation...")
result1, from_cache1 = wf.get_or_calculate(
    atoms=si,
    calculation_params={
        'protocol': 'fast',
        'machine': 'medusa',
        'code_version' : '7.4.1',
        'pseudopotentials_config' : 'default'
    },
    calculation_method='scf',
    input_structure_id=input_structure_id  # ← Track input structure
)

print(f"✅ First run - From cache: {from_cache1}")
if not from_cache1:
    print("   (Calculation executed and stored)")

# STEP 3: Run the SAME calculation again (should hit cache)
print("\nSTEP 3: Second SCF calculation (identical parameters)")
print("-"*70)
print("Running same calculation again...")
result2, from_cache2 = wf.get_or_calculate(
    atoms=si,
    calculation_params={
        'protocol': 'fast',
        'machine': 'medusa',
        'code_version' : '7.4.1',
        'pseudopotentials_config' : 'default'
    },
    calculation_method='scf',
    input_structure_id=input_structure_id
)

print(f"✅ Second run - From cache: {from_cache2}")
if from_cache2:
    print("   🎯 CACHE HIT! Retrieved previous result instantly")
else:
    print("   (Calculation executed again)")

# STEP 4: Show provenance info
print("\nSTEP 4: Provenance summary")
print("-"*70)

# Count total structures in database
total_structures = len(wf.db)
print(f"✅ Total structures in database: {total_structures}")

# Get calculation hash to query provenance
calc_params = {
    'protocol': 'fast',
    'machine': 'medusa',
    'code_version' : '7.4.1',
    'pseudopotentials_config' : 'default'
}
calc_hash = wf._compute_calculation_hash(si, calc_params)

# Query provenance
calc_record = wf.provenance.query_by_hash(calc_hash)
if calc_record:
    print(f"✅ Found in provenance database:")
    print(f"   Calculation ID: {calc_record['id']}")
    print(f"   Hash: {calc_record['calculation_hash'][:16]}...")
    print(f"   Input structure ID: {calc_record['input_structure_id']}")
    print(f"   Success: {calc_record['success']}")
    
    # Show execution history
    executions = wf.provenance.query_executions(calc_record['id'])
    print(f"   Executions: {len(executions)}")
    for i, exec_rec in enumerate(executions, 1):
        print(f"     {i}. Machine: {exec_rec['machine']}, "
              f"Time: {exec_rec['execution_timestamp']}")

print("\n" + "="*70)
print("✅ Provenance test completed!")
print("="*70)
