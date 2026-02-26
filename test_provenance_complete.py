#!/usr/bin/env python3
"""
Complete provenance example with structure tracking
"""

from ase.build import bulk
from xespresso.db import DatabaseWorkflow
import json

# Create Silicon structure
si = bulk('Si', 'diamond', a=5.43)
print("Silicon structure created")

# Initialize workflow
wf = DatabaseWorkflow(db_path='provenance_complete.db', ase_db_path='ase_complete.db')

# STEP 1: Store input structure in ASE database to get its ID
print("\n" + "="*70)
print("STEP 1: Store input structure")
print("="*70)

input_structure_id = wf.db.write(
    si,
    structure_name='Si_diamond',
    lattice_parameter=5.43
)
print(f"✅ Input structure stored with ID: {input_structure_id}")

# STEP 2: Run SCF calculation with provenance tracking
print("\n" + "="*70)
print("STEP 2: First SCF calculation (not cached)")
print("="*70)

result1, from_cache1 = wf.get_or_calculate(
    atoms=si,
    calculation_params={
        'protocol': 'fast',
        'machine': 'medusa',
        'code_version': '7.4.1',
        'pseudopotentials_config': 'default'
    },
    calculation_method='scf',
    input_structure_id=input_structure_id  # Link to input structure
)

print(f"✅ First run - From cache: {from_cache1}")

# STEP 3: Query provenance to see the calculation
print("\n" + "="*70)
print("STEP 3: Query provenance database")
print("="*70)

# Get the calculation from provenance
calc_hash = wf._compute_calculation_hash(si, {
    'protocol': 'fast',
    'machine': 'medusa',
    'code_version': '7.4.1',
    'pseudopotentials_config': 'default'
})

calc = wf.provenance.query_by_hash(calc_hash)
if calc:
    print(f"✅ Found calculation in provenance:")
    print(f"   ID: {calc['id']}")
    print(f"   Hash: {calc['calculation_hash'][:8]}...")
    print(f"   Input structure ID: {calc['input_structure_id']}")
    print(f"   Output structure ID: {calc['output_structure_id']}")
    print(f"   Energy: {calc['energy']} eV" if calc['energy'] else "   Energy: Not stored")
    print(f"   Success: {calc['success']}")

# STEP 4: Run same calculation again (should hit cache)
print("\n" + "="*70)
print("STEP 4: Same calculation again (should hit cache)")
print("="*70)

result2, from_cache2 = wf.get_or_calculate(
    atoms=si,
    calculation_params={
        'protocol': 'fast',
        'machine': 'medusa',
        'code_version': '7.4.1',
        'pseudopotentials_config': 'default'
    },
    calculation_method='scf',
    input_structure_id=input_structure_id
)

print(f"✅ Second run - From cache: {from_cache2}")

# STEP 5: Show provenance summary
print("\n" + "="*70)
print("STEP 5: Provenance Summary")
print("="*70)

# Query execution history
exec_records = wf.provenance.query_executions(calc['id'])
print(f"✅ Total executions: {len(exec_records)}")
for i, exec_rec in enumerate(exec_records, 1):
    print(f"   Execution {i}:")
    print(f"      Machine: {exec_rec['machine']}")
    print(f"      Success: {exec_rec['success']}")
    print(f"      Timestamp: {exec_rec['execution_timestamp']}")

print("\n" + "="*70)
print("🎉 Complete provenance test finished!")
print("="*70)
print(f"Total structures in database: {len(wf.db)}")
print(f"Cache test: Second run used cache = {from_cache2}")
