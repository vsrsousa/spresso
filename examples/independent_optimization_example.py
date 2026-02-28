#!/usr/bin/env python3
"""
Example: Independent Parameter Optimization

This example demonstrates the new optimize_parameters method that performs
two independent convergence runs for maximum efficiency:
1. Converge ecutwfc using coarse kspacing
2. Converge kspacing using lower ecutwfc
"""

from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow

def demonstrate_independent_optimization():
    print("=== Independent Parameter Optimization Demo ===\n")

    # Create a silicon structure
    atoms = bulk('Si', 'diamond', a=5.43)

    # Define pseudopotentials
    pseudopotentials = {'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'}

    print("This method performs TWO independent convergence studies:")
    print("1. Phase 1: Converge ecutwfc using coarse kspacing (efficient)")
    print("2. Phase 2: Converge kspacing using lower ecutwfc (fast)")
    print()
    print("KEY: Precision controls PARAMETER RANGES, criteria control CONVERGENCE CHECKS")
    print("Ranges are automatically adjusted based on pseudopotential requirements!")
    print("You can mix any precision level with any convergence criteria!")
    print()

    # Example 1: Simple optimization with medium precision
    print("Example 1: Medium precision optimization")
    print("-" * 40)

    try:
        optimal_params = ConvergenceWorkflow.optimize_parameters(
            atoms=atoms,
            pseudopotentials=pseudopotentials,
            precision='medium',
            verbose=True
        )

        print("Results:")
        print(f"  Optimal ecutwfc: {optimal_params['ecutwfc']} Ry")
        print(f"  Optimal kspacing: {optimal_params['kspacing']} Å⁻¹")
        print(f"  Precision level: {optimal_params['precision']}")
        print(f"  Convergence criteria: {optimal_params['convergence_criteria']}")
        print()

        # Show convergence data
        print("Phase 1 data (ecutwfc convergence):")
        for result in optimal_params['ecutwfc_convergence_data']:
            print(".1f"
        print()
        print("Phase 2 data (kspacing convergence):")
        for result in optimal_params['kspacing_convergence_data']:
            print(".3f"
        print()

    except Exception as e:
        print(f"Error during optimization: {e}")
        print("Note: This is expected if QE/xespresso is not properly configured")
        print()

    # Example 2: Low precision with strict criteria (unusual combination!)
    print("Example 2: Low precision parameters + strict convergence criteria")
    print("-" * 60)

    try:
        optimal_params = ConvergenceWorkflow.optimize_parameters(
            atoms=atoms,
            pseudopotentials=pseudopotentials,
            precision='low',  # Coarse parameter ranges
            convergence_criteria_list=['energy', 'forces', 'geometry'],  # But strict criteria
            verbose=True
        )

        print("Results:")
        print(f"  Optimal ecutwfc: {optimal_params['ecutwfc']} Ry")
        print(f"  Optimal kspacing: {optimal_params['kspacing']} Å⁻¹")
        print(f"  Precision: {optimal_params['precision']} (coarse ranges)")
        print(f"  Criteria: {optimal_params['convergence_criteria']} (strict checks)")
        print()

    except Exception as e:
        print(f"Error during optimization: {e}")
        print()

    # Example 3: Ultra precision with minimal criteria
    print("Example 3: Ultra precision parameters + minimal convergence criteria")
    print("-" * 65)

    try:
        optimal_params = ConvergenceWorkflow.optimize_parameters(
            atoms=atoms,
            pseudopotentials=pseudopotentials,
            precision='ultra',  # Very fine parameter ranges
            convergence_criteria_list=['energy'],  # But only energy check
            verbose=True
        )

        print("Results:")
        print(f"  Optimal ecutwfc: {optimal_params['ecutwfc']} Ry")
        print(f"  Optimal kspacing: {optimal_params['kspacing']} Å⁻¹")
        print(f"  Precision: {optimal_params['precision']} (very fine ranges)")
        print(f"  Criteria: {optimal_params['convergence_criteria']} (minimal checks)")
        print()

    except Exception as e:
        print(f"Error during optimization: {e}")
        print()

    print("=== Key Advantages ===")
    print("✓ Independent convergence runs (more efficient)")
    print("✓ Phase 1: ecutwfc with coarse kspacing")
    print("✓ Phase 2: kspacing with lower ecutwfc")
    print("✓ Automatic parameter selection")
    print("✓ Returns optimal parameters directly")
    print("✓ Minimal computational cost")
    print()

    print("=== Usage in Production Code ===")
    print("""
# Get optimal parameters
optimal = ConvergenceWorkflow.optimize_parameters(
    atoms=your_atoms,
    pseudopotentials=your_pseudos,
    precision='medium'
)

# Use optimal parameters in production calculations
workflow = CalculationWorkflow(
    atoms=your_atoms,
    pseudopotentials=your_pseudos,
    ecutwfc=optimal['ecutwfc'],
    kspacing=optimal['kspacing']
)
    """)

if __name__ == '__main__':
    demonstrate_independent_optimization()