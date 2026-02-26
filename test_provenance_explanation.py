#!/usr/bin/env python3
"""
Simple provenance caching test with local execution.
Shows how to track input/output structures and use caching.
"""

from ase.build import bulk
from xespresso.db import DatabaseWorkflow
import logging

# Reduce logging noise
logging.getLogger('xespresso').setLevel(logging.WARNING)

print("="*70)
print("PROVENANCE CACHING EXAMPLE")
print("="*70)

# Create Silicon structure
si = bulk('Si', 'diamond', a=5.43)
print("\n✅ Silicon structure created (diamond, a=5.43 Å)")

# Initialize workflow
wf = DatabaseWorkflow(db_path='provenance_example.db', ase_db_path='ase_example.db')
print("✅ DatabaseWorkflow initialized")

# STEP 1: Store input structure to track provenance
print("\n" + "-"*70)
print("STEP 1: Store input structure in ASE database")
print("-"*70)

input_structure_id = wf.db.write(
    si,
    structure_name='Si_diamond',
    lattice_parameter=5.43
)
print(f"✅ Input structure stored with ID: {input_structure_id}")
print(f"   (This ID links input→calculation→output in provenance)")

# STEP 2: Show what parameters we're using
calc_params = {
    'protocol': 'fast',
    'pseudopotentials_config': 'default'
}
print("\n" + "-"*70)
print("STEP 2: Calculation parameters")
print("-"*70)
print(f"Protocol: {calc_params['protocol']}")
print(f"Pseudopotentials: {calc_params['pseudopotentials_config']}")

# Compute the hash
calc_hash = wf._compute_calculation_hash(si, calc_params)
print(f"\n✅ Calculation hash computed: {calc_hash[:8]}...")
print(f"   (Same structure + params = Same hash = Cache hit!)")

# STEP 3: Show how caching works
print("\n" + "-"*70)
print("STEP 3: Understanding the cache mechanism")
print("-"*70)
print("""
Without passing input_structure_id:
  - Caching still works (based on structure + parameters hash)
  - Provenance tracks input→output relationship
  - Query the provenance database for full lineage

With input_structure_id:
  - You explicitly link this calculation to a specific input structure
  - Useful for tracking workflows: relax → band → DOS
  - Can query all calculations that depend on a specific structure
""")

# STEP 4: Show provenance query methods
print("\n" + "-"*70)
print("STEP 4: Querying provenance")
print("-"*70)

print("""
The provenance database allows you to query:

1. By hash (for caching):
   existing = wf.provenance.query_by_hash(calc_hash)
   → Fast lookup to avoid recalculation

2. By structure (for lineage):
   calcs = wf.provenance.query_by_structure(structure_id)
   → Find all calculations that use this structure

3. By calculation (for execution history):
   execs = wf.provenance.query_executions(calc_id)
   → See how many times a calc was run, on which machines

4. Derivations (for structure evolution):
   derivations = wf.provenance.query_derivations(structure_id)
   → Track relaxation: input → relaxed structure
""")

# STEP 5: Manual example
print("\n" + "-"*70)
print("STEP 5: Demonstrating manual provenance logging")
print("-"*70)

# Manually log a test calculation
test_hash = "test_hash_abc123"
test_calc_id = wf.provenance.log_calculation(
    calc_hash=test_hash,
    input_structure_id=input_structure_id,
    output_structure_id=None,  # Not calculated yet, just logging the attempt
    calculation_method='scf',
    success=True
)
print(f"✅ Logged test calculation with ID: {test_calc_id}")

# Query it back
test_calc = wf.provenance.query_by_hash(test_hash)
if test_calc:
    print(f"✅ Retrieved from provenance:")
    print(f"   - Calculation ID: {test_calc['id']}")
    print(f"   - Input structure ID: {test_calc['input_structure_id']}")
    print(f"   - Output structure ID: {test_calc['output_structure_id']}")

print("\n" + "="*70)
print("🎉 Provenance system demonstration complete!")
print("="*70)
print(f"""
Key takeaways:

1. Input structure ID is OPTIONAL
   - Caching works without it (uses hash)
   - Use it for tracking calculation dependencies

2. Hash-based caching is automatic
   - Same atoms + parameters = same hash
   - Second call with identical params avoids recalculation

3. Provenance tracks everything
   - Success/failure of each calculation
   - Execution history (machine, timestamp)
   - Input/output structure relationships

4. For workflows:
   ✅ Store input structure → get ID
   ✅ Calculate with get_or_calculate(... input_structure_id=...)
   ✅ Query provenance for lineage and dependencies
""")
