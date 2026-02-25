"""
Example demonstrating DatabaseWorkflow with ASE Database and Provenance tracking.

This example shows:
1. Storing structures in ASE database
2. Caching calculations to avoid redundancy
3. Tracking provenance of calculations
4. Inheriting convergence parameters
5. Analyzing large datasets
"""

import os
from ase.build import bulk, add_adsorbate
from ase.io import read

from xespresso.db import DatabaseWorkflow, get_structure_history, validate_consistency
from xespresso.db.queries import (
    get_all_derivatives,
    export_to_dataframe,
    get_calculation_stats,
)

# Setup environment
os.environ['ASE_ESPRESSO_PSEUDO'] = '/home/vinicius/pseudos/SSSP_1.3.0_PBE_efficiency'

print("\n" + "="*80)
print("DatabaseWorkflow Example: Structure Evolution and Calculation Caching")
print("="*80)

# ============================================================================
# Step 1: Initialize DatabaseWorkflow
# ============================================================================
print("\n[Step 1] Initialize DatabaseWorkflow")
print("-" * 80)

db_workflow = DatabaseWorkflow(
    db_path='~/.xespresso/database.db',
    provenance_path='~/.xespresso/provenance.db'
)

# ============================================================================
# Step 2: Store original structure and run SCF
# ============================================================================
print("\n[Step 2] Store original Si structure and run SCF")
print("-" * 80)

# Create structure
atoms_original = bulk('Si', 'diamond', a=5.431)

# Define calculation parameters
scf_params = {
    'protocol': 'moderate',
    'pseudopotentials': {'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'},
    'machine': 'medusa'
}

# Run SCF (first time - will calculate)
print("Running SCF on original structure...")
scf_result, from_cache = db_workflow.get_or_calculate(
    atoms_original,
    scf_params,
    calculation_method='scf'
)
print(f"SCF energy: {scf_result.get_potential_energy():.6f} eV")
print(f"Retrieved from cache: {from_cache}")

# Get the structure ID for this result
original_struct_id = len(db_workflow.db) - 1  # Last added structure

# ============================================================================
# Step 3: Run relaxation on original structure
# ============================================================================
print("\n[Step 3] Run structural relaxation")
print("-" * 80)

print("Running vc-relax on original structure...")
relax_params = scf_params.copy()

relax_result, from_cache = db_workflow.get_or_calculate(
    atoms_original,
    relax_params,
    calculation_method='relax'
)
print(f"Relaxed energy: {relax_result.get_potential_energy():.6f} eV")

relaxed_struct_id = len(db_workflow.db) - 1

# Log the derivation
energy_change = relax_result.get_potential_energy() - scf_result.get_potential_energy()
db_workflow.log_structure_derivation(
    source_structure_id=original_struct_id,
    derived_structure_id=relaxed_struct_id,
    derivation_method='vc-relax',
    energy_change=energy_change
)
print(f"Energy change from relaxation: {energy_change:.6f} eV")

# ============================================================================
# Step 4: Calculate properties on relaxed structure
# ============================================================================
print("\n[Step 4] Calculate phonon properties on relaxed structure")
print("-" * 80)

print("Running phonon calculation on relaxed structure...")
phonon_result, from_cache = db_workflow.run_properties_on_structure(
    structure_id=relaxed_struct_id,
    calc_type='phonon'
)
print(f"Phonon calculation complete (from cache: {from_cache})")

# ============================================================================
# Step 5: Re-run same SCF calculation (should be cached!)
# ============================================================================
print("\n[Step 5] Re-run SCF on original structure (should be cached)")
print("-" * 80)

atoms_original_copy = bulk('Si', 'diamond', a=5.431)
print("Running SCF again on same structure...")
scf_result_2, from_cache = db_workflow.get_or_calculate(
    atoms_original_copy,
    scf_params,
    calculation_method='scf'
)
print(f"✅ Retrieved from cache: {from_cache}")
print(f"Energy matches: {abs(scf_result_2.get_potential_energy() - scf_result.get_potential_energy()) < 1e-6}")

# ============================================================================
# Step 6: Analyze structure evolution
# ============================================================================
print("\n[Step 6] Analyze structure evolution")
print("-" * 80)

# Get lineage of relaxed structure
history = get_structure_history(db_workflow.provenance, relaxed_struct_id)
print(f"Structure evolution history ({len(history)} steps):")
for i, step in enumerate(history):
    print(f"  Step {i+1}: {step['method']} → Energy: {step['energy']:.6f} eV (ID: {step['id']})")

# ============================================================================
# Step 7: Check consistency
# ============================================================================
print("\n[Step 7] Verify calculation consistency")
print("-" * 80)

is_consistent = validate_consistency(db_workflow.provenance, relaxed_struct_id)
if is_consistent:
    print("✅ All calculations on this structure used consistent parameters")
else:
    print("⚠️  Inconsistent parameters detected")

# ============================================================================
# Step 8: Export data for large-scale analysis
# ============================================================================
print("\n[Step 8] Export data for analysis")
print("-" * 80)

# Get statistics
stats = get_calculation_stats(db_workflow.provenance)
print(f"Total calculations: {stats['total_calculations']}")
print(f"By method: {stats['by_method']}")
print(f"Average convergence steps: {stats['avg_convergence_steps']:.1f}")

# Export to DataFrame
df = export_to_dataframe(db_workflow.db, db_workflow.provenance)
print(f"\nExported to DataFrame:")
print(df[['formula', 'energy', 'method', 'protocol', 'ecutwfc']])

# ============================================================================
# Step 9: Get derivatives
# ============================================================================
print("\n[Step 9] Find all structures derived from original")
print("-" * 80)

derivatives = get_all_derivatives(db_workflow.provenance, original_struct_id)
print(f"Found {len(derivatives)} derived structures:")
for deriv in derivatives:
    print(f"  ID: {deriv['derived_structure_id']}, Method: {deriv['method']}, "
          f"Energy change: {deriv['energy_change']:.6f} eV")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "="*80)
print("Summary")
print("="*80)
print("""
Key features demonstrated:

✅ Automatic Caching
   - First SCF calculation executed and stored
   - Second SCF calculation retrieved from cache (zero overhead)

✅ Structure Derivation Tracking
   - Relaxation logged as derivation of original structure
   - Energy changes recorded

✅ Parameter Inheritance
   - Phonon calculation inherited convergence params from SCF
   - Automatic consistency checking

✅ Provenance & Auditability
   - Complete calculation history stored
   - Machine information recorded
   - Dependencies tracked

✅ Data Analysis
   - Easy export to pandas for statistical analysis
   - Query by protocol, element, or method
   - Compare structures and calculations

For large-scale studies with thousands of calculations:
- Avoid redundant computations using cache
- Maintain parameter consistency across related calculations
- Track complete provenance for scientific reproducibility
- Analyze results systematically with database queries
""")

# Cleanup
db_workflow.close()
print("\nDatabaseWorkflow closed successfully")
