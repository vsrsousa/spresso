#!/usr/bin/env python3
"""
How to retrieve structures from the ASE database for future calculations.
Demonstra as diferentes formas de buscar e reutilizar estruturas.
"""

from ase.build import bulk
from xespresso.db import DatabaseWorkflow
import logging

# Reduce logging noise
logging.getLogger('xespresso').setLevel(logging.WARNING)

print("="*70)
print("RETRIEVING STRUCTURES FROM DATABASE")
print("="*70)

# Initialize workflow
wf = DatabaseWorkflow(db_path='retrieve_example.db', ase_db_path='retrieve_example_ase.db')
print("\n✅ DatabaseWorkflow initialized")

# STEP 1: Store some structures
print("\n" + "-"*70)
print("STEP 1: Store initial structures")
print("-"*70)

si = bulk('Si', 'diamond', a=5.43)
id_si = wf.db.write(si, structure_name='Si_diamond', lattice_param=5.43)
print(f"✅ Stored Si: ID = {id_si}")

fe = bulk('Fe', 'bcc', a=2.87)
id_fe = wf.db.write(fe, structure_name='Fe_bcc', lattice_param=2.87)
print(f"✅ Stored Fe: ID = {id_fe}")

# STEP 2: Retrieve by ID (most straightforward)
print("\n" + "-"*70)
print("STEP 2: Retrieve by ID")
print("-"*70)

# Recover Si structure
si_retrieved = wf.db.get_atoms(id_si)
print(f"✅ Retrieved structure ID {id_si}:")
print(f"   Formula: {si_retrieved.get_chemical_formula()}")
print(f"   Lattice: {si_retrieved.cell[0,0]:.3f} Å")

# STEP 3: Query with filters
print("\n" + "-"*70)
print("STEP 3: Query structures with filters")
print("-"*70)

print("\n3a) Get all structures in database:")
print(f"   Total structures: {len(wf.db)}")

print("\n3b) Query by key-value pairs:")
# Get structures with specific lattice parameter
try:
    for row in wf.db.select('structure_name=Si_diamond'):
        atoms = row.toatoms()
        print(f"   ✅ Found: {row['structure_name']}")
        print(f"      ID: {row.id}")
        print(f"      Formula: {atoms.get_chemical_formula()}")
except Exception as e:
    print(f"   Note: Need to use different query syntax")

print("\n3c) Iterate through all structures:")
for i, row in enumerate(wf.db.select()):
    atoms = row.toatoms()
    print(f"   {i+1}. ID {row.id}: {atoms.get_chemical_formula()} "
          f"({row.get('structure_name', 'unnamed')})")

# STEP 4: Store metadata for easy retrieval
print("\n" + "-"*70)
print("STEP 4: Use metadata to find structures")
print("-"*70)

# Store with meaningful metadata
mo = bulk('Mo', 'bcc', a=3.15)
id_mo = wf.db.write(
    mo,
    structure_name='Mo_bcc',
    lattice_param=3.15,
    element='Mo',
    crystal_system='bcc',
    calculation_type='scf'
)
print(f"✅ Stored Mo with rich metadata: ID = {id_mo}")

print("\n   Retrieving by different queries:")
for row in wf.db.select():
    if row.get('element') == 'Mo':
        print(f"   ✅ Found Mo structure: ID = {row.id}")
        atoms = row.toatoms()
        print(f"      Formula: {atoms.get_chemical_formula()}")
        print(f"      Cell: {atoms.cell[0,0]:.3f} Å")

# STEP 5: Use retrieved structure in calculation
print("\n" + "-"*70)
print("STEP 5: Use retrieved structure for new calculation")
print("-"*70)

print(f"""
Example workflow:

# Step A: Retrieve previous structure
retrieved_atoms = wf.db.get_atoms(structure_id)

# Step B: Use in new calculation
result, cached = wf.get_or_calculate(
    atoms=retrieved_atoms,
    calculation_params={{'protocol': 'fast', ...}},
    calculation_method='scf',
    input_structure_id=structure_id  # Link to original
)

# Step C: Store new result
result_id = wf.db.write(result, structure_name='Si_scf_result')

# Step D: Query the lineage
lineage = wf.provenance.query_by_structure(structure_id)
→ Shows all calculations using this structure
""")

# STEP 6: Practical retrieval example
print("\n" + "-"*70)
print("STEP 6: Complete retrieval + calculation example")
print("-"*70)

# Get the Si structure we stored earlier
si_for_calc = wf.db.get_atoms(id_si)
print(f"✅ Retrieved Si structure (ID={id_si}):")
print(f"   Atoms: {si_for_calc}")
print(f"   Lattice parameter: {si_for_calc.cell[0,0]:.3f} Å")

print("\nNow this structure can be used for:")
print("✅ SCF calculation: wf.get_or_calculate(atoms=si_for_calc, ...)")
print("✅ Relaxation:     wf.get_or_calculate(atoms=si_for_calc, ..., method='relax')")
print("✅ Band structure: wf.get_or_calculate(atoms=si_for_calc, ..., method='band')")
print("✅ Phonons:        wf.get_or_calculate(atoms=si_for_calc, ..., method='phonon')")

# STEP 7: Show database info
print("\n" + "-"*70)
print("STEP 7: Database summary")
print("-"*70)

# Count groups
print(f"Total structures stored: {len(wf.db)}")

# Show structure IDs available
print("\nStructure IDs available for retrieval:")
for row in wf.db.select():
    name = row.get('structure_name', f'structure_{row.id}')
    print(f"  → retrieve_structure = wf.db.get_atoms({row.id})  # {name}")

print("\n" + "="*70)
print("🎉 Examples complete!")
print("="*70)
print("""
SUMMARY:

To retrieve a structure:
  1. By ID (direct):
     atoms = wf.db.get_atoms(structure_id)

  2. By query:
     for row in wf.db.select():
         atoms = row.toatoms()

  3. By metadata filter:
     for row in wf.db.select('element=Si'):
         atoms = row.toatoms()

With retrieved structure, use for new calculation:
  result, cached = wf.get_or_calculate(
      atoms=atoms,
      calculation_params={...},
      calculation_method='scf',
      input_structure_id=structure_id  # Optional but recommended
  )
""")
