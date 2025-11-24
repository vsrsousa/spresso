"""
Example demonstrating the kpts_from_spacing utility function.

This example shows how to use kpts_from_spacing to convert k-spacing
to k-point grids without worrying about 2π normalization.
"""

from ase.build import bulk
from xespresso import kpts_from_spacing
import numpy as np

print("=" * 70)
print("kpts_from_spacing - Utility Function Demo")
print("=" * 70)

# Example 1: Basic usage
print("\n1. Basic Usage")
print("-" * 70)
atoms = bulk('Si', cubic=True)

# Simple usage - no need for /(2*π)
kpts = kpts_from_spacing(atoms, 0.20)
print(f"atoms = bulk('Si', cubic=True)")
print(f"kpts = kpts_from_spacing(atoms, 0.20)")
print(f"Result: kpts = {kpts}")
print("\n✓ Clean and simple - no manual 2π normalization needed!")

# Example 2: Comparison with manual approach
print("\n2. Comparison with Manual Approach")
print("-" * 70)

# Old way (manual normalization)
from ase.io.espresso import kspacing_to_grid
kpts_manual = kspacing_to_grid(atoms, 0.20 / (2 * np.pi))

print("Old way (confusing):")
print(f"  kpts = kspacing_to_grid(atoms, 0.20/(2*np.pi))")
print(f"  Result: {tuple(kpts_manual)}")

print("\nNew way (clean):")
print(f"  kpts = kpts_from_spacing(atoms, 0.20)")
print(f"  Result: {kpts}")

print(f"\n✓ Both give the same result: {kpts == tuple(kpts_manual)}")

# Example 3: Different k-spacing values
print("\n3. Testing Different K-spacing Values")
print("-" * 70)

kspacing_values = [0.5, 0.3, 0.2, 0.15]
print(f"Structure: Si (cubic, a=5.43 Å)")
print()
for ksp in kspacing_values:
    kpts = kpts_from_spacing(atoms, ksp)
    print(f"  kspacing = {ksp:.2f} Å⁻¹  →  kpts = {kpts}")

# Example 4: Different structures
print("\n4. Different Crystal Structures")
print("-" * 70)

structures = {
    'Si (diamond)': bulk('Si', cubic=True),
    'Fe (bcc)': bulk('Fe', cubic=True),
    'Al (fcc)': bulk('Al', cubic=True),
}

ksp = 0.3
print(f"Using kspacing = {ksp} Å⁻¹ for all structures:\n")
for name, struct in structures.items():
    kpts = kpts_from_spacing(struct, ksp)
    a = struct.cell.cellpar()[0]
    print(f"  {name:15s} (a={a:.2f} Å)  →  kpts = {kpts}")

# Example 5: Using in a calculation setup
print("\n5. Using in a Calculation Setup")
print("-" * 70)

print("""
# Example calculation setup
from ase.build import bulk
from xespresso import Espresso, kpts_from_spacing

atoms = bulk('Si', cubic=True)
kpts = kpts_from_spacing(atoms, 0.25)  # Clean k-spacing specification

calc = Espresso(
    pseudopotentials={'Si': 'Si.pbe.UPF'},
    ecutwfc=50,
    kpts=kpts,  # Use the calculated k-points
    label='scf/silicon'
)

atoms.calc = calc
energy = atoms.get_potential_energy()
""")

print("✓ Simple and intuitive API for k-point specification")

# Example 6: Integration with workflows
print("\n6. Integration with Workflows")
print("-" * 70)

print("""
The kpts_from_spacing function is also used internally by the
CalculationWorkflow class:

from xespresso import CalculationWorkflow

# Workflow automatically uses kpts_from_spacing internally
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={'Si': 'Si.pbe.UPF'},
    quality='moderate',
    kspacing=0.3  # Workflow uses kpts_from_spacing internally
)

# Or use it manually for custom setups
from xespresso import Espresso, kpts_from_spacing

kpts = kpts_from_spacing(atoms, 0.3)
calc = Espresso(..., kpts=kpts, ...)
""")

print("✓ Works seamlessly with both manual and workflow-based setups")

# Summary
print("\n" + "=" * 70)
print("Summary")
print("=" * 70)
print("""
✓ kpts_from_spacing(atoms, kspacing) simplifies k-point specification
✓ No need to manually apply /(2*π) normalization
✓ Works with any ASE Atoms structure
✓ Can be used in custom calculation setups or with workflows
✓ Returns standard tuple of k-points (kx, ky, kz)

Usage:
  from xespresso import kpts_from_spacing
  kpts = kpts_from_spacing(atoms, 0.20)  # Clean and simple!
""")

print("=" * 70)
print("Example complete!")
print("=" * 70)
