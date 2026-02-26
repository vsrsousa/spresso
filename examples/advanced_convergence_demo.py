#!/usr/bin/env python3
"""
Example: Advanced Convergence Criteria

This example demonstrates the enhanced ConvergenceWorkflow with multiple
convergence criteria: energy, forces, geometry, and magnetic moments.
"""

from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow

def demonstrate_convergence_criteria():
    print("=== Advanced Convergence Criteria Demo ===\n")

    # Create a simple silicon structure
    atoms = bulk('Si', 'diamond', a=5.43)

    # Define pseudopotentials
    pseudopotentials = {'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'}

    print("Available convergence criteria:")
    print("• 'energy': Energy convergence (meV/atom)")
    print("• 'forces': Force convergence (eV/Å)")
    print("• 'geometry': Geometry optimization convergence")
    print("• 'magnetic_moments': Magnetic moment convergence")
    print()

    # Example 1: Energy-only convergence (fastest)
    print("1. Energy-only convergence (precision='low'):")
    workflow1 = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        precision='low'
    )
    print(f"   Criteria: {workflow1.convergence_criteria_list}")
    print(f"   Tolerances: {workflow1.convergence_criteria}")
    print()

    # Example 2: Energy + Forces convergence
    print("2. Energy + Forces convergence (precision='medium'):")
    workflow2 = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        precision='medium'
    )
    print(f"   Criteria: {workflow2.convergence_criteria_list}")
    print(f"   Tolerances: {workflow2.convergence_criteria}")
    print()

    # Example 3: Custom criteria
    print("3. Custom convergence criteria:")
    custom_criteria = ['energy', 'forces', 'geometry']
    custom_tolerances = {
        'energy_tolerance': 5e-4,    # 0.5 meV/atom
        'force_tolerance': 0.05,     # eV/Å
        'geometry_tolerance': 0.005, # Å
    }
    workflow3 = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        convergence_criteria_list=custom_criteria,
        convergence_criteria=custom_tolerances
    )
    print(f"   Criteria: {workflow3.convergence_criteria_list}")
    print(f"   Tolerances: {workflow3.convergence_criteria}")
    print()

    print("=== Key Points ===")
    print("• ecutwfc and kspacing convergence are INDEPENDENT")
    print("• For each ecutwfc, we test if it converges with kspacing refinement")
    print("• User can choose ANY combination of criteria, regardless of precision level")
    print("• All precision levels define tolerances for ALL criteria types")
    print()

    # Example 4: Energy-only with ultra precision
    print("4. Energy-only with ultra precision:")
    workflow4 = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        precision='ultra',
        convergence_criteria_list=['energy']  # Only energy, but ultra tolerance
    )
    print(f"   Criteria: {workflow4.convergence_criteria_list}")
    print(f"   Energy tolerance: {workflow4.convergence_criteria['energy_tolerance']*1000:.1f} meV/atom")
    print(f"   (Note: ultra precision gives 0.5 meV/atom even for energy-only)")
    print()

    # Example 5: Forces-only with low precision
    print("5. Forces-only with low precision:")
    workflow5 = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        precision='low',
        convergence_criteria_list=['forces']  # Only forces, but low tolerance
    )
    print(f"   Criteria: {workflow5.convergence_criteria_list}")
    print(f"   Force tolerance: {workflow5.convergence_criteria['force_tolerance']:.1f} eV/Å")
    print(f"   (Note: low precision gives 0.5 eV/Å even for forces-only)")
    print()

    print("=== Running Convergence Study ===")
    print("The workflow will:")
    print("1. Start with ecutwfc = 30 Ry")
    print("2. For each ecutwfc, test kspacing values from coarse to fine")
    print("3. Check if energy/forces/geometry converge with kspacing")
    print("4. If converged, stop; if not, increase ecutwfc by 10 Ry and repeat")
    print("5. Stop when ALL specified criteria converge with kspacing")
    print()

    # Uncomment to run actual convergence study
    # results = workflow2.run_convergence_study(verbose=True)
    # print(f"Converged at ecutwfc={results.converged_ecutwfc} Ry")

if __name__ == '__main__':
    demonstrate_convergence_criteria()