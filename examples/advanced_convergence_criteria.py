#!/usr/bin/env python3
"""
Example: Advanced Convergence Criteria

This example demonstrates the new convergence criteria system that allows
specifying which physical quantities to check for convergence.
"""

from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow

def demonstrate_convergence_criteria():
    print("=== Advanced Convergence Criteria Demo ===\n")

    # Create a silicon structure
    atoms = bulk('Si', 'diamond', a=5.43)

    # Define pseudopotentials
    pseudopotentials = {'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'}

    print("Available convergence criteria:")
    print("• 'energy': Total energy convergence")
    print("• 'forces': Maximum force convergence")
    print("• 'geometry': Atomic position/geometry convergence")
    print("• 'magnetic_moments': Magnetic moment convergence")
    print()

    # Example 1: Energy only (fastest)
    print("1. Energy-only convergence (fastest):")
    workflow1 = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        precision='medium',
        convergence_criteria_list=['energy']
    )
    print(f"   Criteria: {workflow1.convergence_criteria_list}")
    print(f"   Tolerances: Energy = {workflow1.convergence_criteria['energy_tolerance']*1000:.1f} meV/atom")
    print()

    # Example 2: Energy + Forces
    print("2. Energy + Forces convergence:")
    workflow2 = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        precision='medium',
        convergence_criteria_list=['energy', 'forces']
    )
    print(f"   Criteria: {workflow2.convergence_criteria_list}")
    print(f"   Tolerances: Energy = {workflow2.convergence_criteria['energy_tolerance']*1000:.1f} meV/atom")
    print(f"               Forces = {workflow2.convergence_criteria['force_tolerance']:.1f} eV/Å")
    print()

    # Example 3: Full convergence (slowest but most accurate)
    print("3. Full convergence (energy + forces + geometry + magnetic):")
    workflow3 = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        precision='high',
        convergence_criteria_list=['energy', 'forces', 'geometry', 'magnetic_moments']
    )
    print(f"   Criteria: {workflow3.convergence_criteria_list}")
    print(f"   Tolerances: Energy = {workflow3.convergence_criteria['energy_tolerance']*1000:.1f} meV/atom")
    print(f"               Forces = {workflow3.convergence_criteria['force_tolerance']:.1f} eV/Å")
    print(f"               Geometry = {workflow3.convergence_criteria['geometry_tolerance']*1000:.1f} Å")
    print(f"               Magnetic = {workflow3.convergence_criteria['magnetic_tolerance']*1000:.1f} μB")
    print()

    # Example 4: Custom tolerances
    print("4. Custom convergence criteria and tolerances:")
    custom_criteria = {
        'energy_tolerance': 1e-3,      # 1 meV/atom
        'force_tolerance': 0.05,       # 0.05 eV/Å
        'geometry_tolerance': 0.005,   # 0.005 Å
        'magnetic_tolerance': 0.0001,  # 0.0001 μB
    }
    workflow4 = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        convergence_criteria_list=['energy', 'forces'],
        convergence_criteria=custom_criteria
    )
    print(f"   Criteria: {workflow4.convergence_criteria_list}")
    print(f"   Custom tolerances: Energy = {workflow4.convergence_criteria['energy_tolerance']*1000:.1f} meV/atom")
    print(f"                       Forces = {workflow4.convergence_criteria['force_tolerance']:.1f} eV/Å")
    print()

    print("=== Running Convergence Study ===")
    print("To run an actual convergence study:")
    print("""
# Run with energy + forces convergence
results = workflow2.run_convergence_study(verbose=True)
print(f"Converged at ecutwfc = {results['ecutwfc'].iloc[-1]} Ry")
    """)

if __name__ == '__main__':
    demonstrate_convergence_criteria()