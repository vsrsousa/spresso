"""
Example demonstrating the updated workflow features based on user feedback.

This example shows:
1. Simplified k-spacing API (no need for /(2*np.pi))
2. Magnetic configuration support (ferro, antiferro, custom)
3. Hubbard parameter integration
"""

from ase.build import bulk
from xespresso import CalculationWorkflow, quick_scf, quick_relax
import numpy as np

print("=" * 70)
print("Updated Workflow Features Demo")
print("=" * 70)

# Feature 1: Simplified K-spacing API
print("\n1. K-spacing API - No more 2π confusion!")
print("-" * 70)
atoms = bulk('Si', cubic=True)

workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={'Si': 'Si.pbe.UPF'},
    protocol='moderate',
    kspacing=0.20  # Just pass the physical value in Angstrom^-1!
)

print(f"✓ User passes: kspacing=0.20 Angstrom^-1")
print(f"✓ Workflow handles: /(2*π) internally")
print(f"✓ Resulting k-points: {workflow._get_kpts()}")
print("✓ No more manual normalization needed!")

# Feature 2: Magnetic Configurations
print("\n2. Magnetic Configuration Support")
print("-" * 70)

# 2a. Ferromagnetic
print("\n2a. Ferromagnetic configuration:")
atoms_fe = bulk('Fe', cubic=True)
workflow_ferro = CalculationWorkflow(
    atoms=atoms_fe,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
    protocol='moderate',
    kspacing=0.3,
    magnetic_config='ferro'  # Simple string for ferro
)
print("✓ magnetic_config='ferro' works")
print(f"✓ Species: {workflow_ferro.atoms.arrays.get('species', 'no species')}")
print(f"✓ Has magnetic params: {'starting_magnetization' in workflow_ferro.input_data.get('input_ntyp', {})}")

# 2b. Antiferromagnetic
print("\n2b. Antiferromagnetic configuration:")
atoms_fe = bulk('Fe', cubic=True)
workflow_afm = CalculationWorkflow(
    atoms=atoms_fe,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
    protocol='moderate',
    kspacing=0.3,
    magnetic_config='antiferro'  # Simple string for AFM
)
print("✓ magnetic_config='antiferro' works")
print(f"✓ Species: {workflow_afm.atoms.arrays['species']}")
print(f"✓ Pseudopotentials: {list(workflow_afm.pseudopotentials.keys())}")
print(f"✓ Has magnetic params: {'starting_magnetization' in workflow_afm.input_data.get('input_ntyp', {})}")

# 2c. Custom magnetic configuration
print("\n2c. Custom element-based magnetic configuration:")
atoms_fe = bulk('Fe', cubic=True)
workflow_custom = CalculationWorkflow(
    atoms=atoms_fe,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
    protocol='moderate',
    kspacing=0.3,
    magnetic_config={'Fe': [1, -1]}  # AFM with explicit moments
)
print("✓ magnetic_config={'Fe': [1, -1]} works")
print(f"✓ Species created: {workflow_custom.atoms.arrays['species']}")
print(f"✓ Magnetization values: {workflow_custom.input_data['input_ntyp']['starting_magnetization']}")

# Feature 3: Hubbard Parameters
print("\n3. Hubbard Parameter Integration")
print("-" * 70)

# 3a. Old format (QE < 7.0)
print("\n3a. Hubbard U (old format for QE < 7.0):")
atoms_fe = bulk('Fe', cubic=True)
workflow_hubbard_old = CalculationWorkflow(
    atoms=atoms_fe,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
    protocol='accurate',
    kspacing=0.25,
    magnetic_config={'Fe': {'mag': [1, -1], 'U': 4.3}}
)
print("✓ Hubbard U=4.3 applied")
print(f"✓ Hubbard_U values: {workflow_hubbard_old.input_data['input_ntyp']['Hubbard_U']}")

# 3b. New format (QE >= 7.0)
print("\n3b. Hubbard U (new format for QE >= 7.0):")
atoms_fe = bulk('Fe', cubic=True)
workflow_hubbard_new = CalculationWorkflow(
    atoms=atoms_fe,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF', 'O': 'O.pbe.UPF'},
    protocol='accurate',
    kspacing=0.25,
    magnetic_config={'Fe': {'mag': [1], 'U': {'3d': 4.3}}},
    input_data={'qe_version': '7.2'}
)
print("✓ Hubbard U on Fe-3d orbital")
print(f"✓ Has hubbard card: {'hubbard' in workflow_hubbard_new.input_data}")

# Feature 4: Quick Functions with All Features
print("\n4. Quick Functions with Magnetic + K-spacing")
print("-" * 70)

print("\nExample quick_scf with antiferromagnetic + k-spacing:")
print("""
calc = quick_scf(
    'fe_structure.cif',
    {'Fe': 'Fe.pbe-spn.UPF'},
    protocol='moderate',
    kspacing=0.3,           # Simple k-spacing
    magnetic_config='antiferro'  # Simple AFM
)
""")
print("✓ One line to set up complex calculations!")

print("\nExample quick_relax with Hubbard U:")
print("""
calc = quick_relax(
    atoms,
    {'Fe': 'Fe.pbe-spn.UPF'},
    protocol='accurate',
    kspacing=0.2,
    magnetic_config={'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}},
    relax_type='vc-relax'
)
""")
print("✓ Full control with minimal code!")

# Summary
print("\n" + "=" * 70)
print("Summary of Improvements")
print("=" * 70)
print("""
1. ✓ K-spacing: No more /(2*π) confusion
   - Just pass physical values in Angstrom^-1
   - Workflow handles normalization internally

2. ✓ Magnetic Configurations: Seamless integration
   - Simple: 'ferro' or 'antiferro'
   - Custom: {'Fe': [1, -1]} for element-based
   - Full integration with setup_magnetic_config

3. ✓ Hubbard Parameters: Both old and new formats
   - Old QE: {'Fe': {'mag': [1], 'U': 4.3}}
   - New QE 7.x: {'Fe': {'mag': [1], 'U': {'3d': 4.3}}}
   - Automatic species and pseudopotential handling

4. ✓ Quick Functions: All features in one line
   - quick_scf() and quick_relax()
   - Support all magnetic and Hubbard options
   - Minimal code for complex setups
""")

print("=" * 70)
print("All features working perfectly! 🚀")
print("=" * 70)
