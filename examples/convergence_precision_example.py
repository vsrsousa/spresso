#!/usr/bin/env python3
"""
Example demonstrating the simplified ConvergenceWorkflow interface.

This example shows how to use the new precision-based parameter optimization
with minimal user input.
"""

from ase.build import bulk
from xespresso.workflow import ConvergenceWorkflow

def main():
    # Create a simple silicon structure
    atoms = bulk('Si', 'diamond', a=5.43)

    # Define pseudopotentials (you would use actual UPF files)
    pseudopotentials = {'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'}

    print("=== Simplified Convergence Workflow Example ===\n")

    # Method 1: Use the optimize_parameters convenience method
    print("1. Using optimize_parameters (simplest interface):")
    print("   optimal = ConvergenceWorkflow.optimize_parameters(")
    print("       atoms=atoms,")
    print("       pseudopotentials=pseudopotentials,")
    print("       precision='medium'")
    print("   )")
    print("   print(f'Optimal ecutwfc: {optimal[\"ecutwfc\"]} Ry')")
    print("   print(f'Optimal kspacing: {optimal[\"kspacing\"]} Å⁻¹')")
    print()

    # Method 2: Create workflow with precision, then run manually
    print("2. Creating workflow with precision, then running manually:")
    print("   workflow = ConvergenceWorkflow(atoms, pseudopotentials, precision='high')")
    print("   results = workflow.run_convergence_study()")
    print("   recommendations = workflow.get_recommendations()")
    print()

    # Method 3: Create from CIF file with precision
    print("3. Creating from CIF file with precision:")
    print("   workflow = ConvergenceWorkflow.from_cif('structure.cif', pseudopotentials, precision='low')")
    print("   results = workflow.run_convergence_study()")
    print()

    # Method 4: Advanced usage with custom ranges still available
    print("\n4. Advanced usage (custom ranges still supported):")
    print("   workflow = ConvergenceWorkflow(atoms, pseudopotentials,")
    print("                                 ecutwfc_range=[40, 60, 80],")
    print("                                 kspacing_range=[0.3, 0.2, 0.1])")
    print("   results = workflow.run_convergence_study()")
    print()

    print("=== Precision Levels ===")
    print("• 'low': Quick calculations (ecutwfc: 30-50 Ry, kspacing: 0.5-0.3 Å⁻¹)")
    print("• 'medium': Balanced speed/accuracy (ecutwfc: 40-70 Ry, kspacing: 0.4-0.2 Å⁻¹)")
    print("• 'high': High accuracy (ecutwfc: 50-90 Ry, kspacing: 0.3-0.12 Å⁻¹)")
    print("• 'ultra': Maximum accuracy (ecutwfc: 60-140 Ry, kspacing: 0.25-0.1 Å⁻¹)")
    print("• Ranges are automatically adjusted based on pseudopotential requirements!")

    print("\n=== Benefits of optimize_parameters ===")
    print("✓ Minimal user input required")
    print("✓ Independent parameter convergence (efficient)")
    print("✓ Automatic parameter optimization")
    print("✓ Precision-based defaults")
    print("✓ Returns optimal parameters directly")
    print()

    print("=== How Independent Convergence Works ===")
    print("1. Converge ecutwfc using coarse kspacing (fast)")
    print("2. Converge kspacing using lower ecutwfc (efficient)")
    print("3. Return optimal parameter combination")
    print()

if __name__ == '__main__':
    main()
    print("✓ Automatic pseudopotential-aware range adjustment")
    print("✓ Backward compatibility maintained")
    print("✓ Advanced options still available")

if __name__ == '__main__':
    main()