#!/usr/bin/env python3
"""
Example: Adaptive Precision System with Structural Complexity

This example demonstrates how the enhanced ConvergenceWorkflow now considers
both pseudopotential requirements AND structural complexity for truly adaptive
parameter selection.
"""

from ase.build import bulk, fcc111
from ase.cluster import Octahedron
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow

def demonstrate_adaptive_precision():
    print("=== Adaptive Precision System Demo ===\n")

    # Define pseudopotentials with mock suggested ecutwfc values
    pseudopotentials_si = {'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'}  # Mock: 30 Ry
    pseudopotentials_cu = {'Cu': 'Cu.pbe.UPF'}  # Mock: 35 Ry  
    pseudopotentials_al = {'Al': 'Al.pbe.UPF'}  # Mock: 25 Ry

    structures = [
        ("Bulk Silicon", bulk('Si', 'diamond', a=5.43)),
        ("Copper Cluster", Octahedron('Cu', 2, cutoff=0.5)),
        ("Aluminum Surface", fcc111('Al', size=(2, 2, 4), vacuum=10.0)),
    ]

    pseudos = [pseudopotentials_si, pseudopotentials_cu, pseudopotentials_al]
    mock_suggested = [30, 35, 25]  # Mock suggested ecutwfc values

    for (name, atoms), pseudo, suggested in zip(structures, pseudos, mock_suggested):
        print(f"Structure: {name} ({len(atoms)} atoms)")
        print(f"Formula: {atoms.get_chemical_formula()}")

        # Create workflow with precision
        workflow = ConvergenceWorkflow(
            atoms=atoms,
            pseudopotentials=pseudo,
            precision='medium'
        )

        # Show structural complexity analysis
        complexity = workflow._analyze_structural_complexity(atoms)
        structural_factor = (
            complexity['surface_factor'] *
            complexity['vacuum_factor'] *
            complexity['heterogeneity_factor']
        )

        print(f"  Structural factors: {complexity}")
        print(f"  Combined factor: {structural_factor:.2f}")
        
        # Demonstrate adaptive calculation
        adjusted_suggested = suggested * structural_factor
        print(f"  Mock suggested ecutwfc: {suggested} Ry")
        print(f"  Adjusted for structure: {suggested} × {structural_factor:.2f} = {adjusted_suggested:.1f} Ry")
        
        # Show what the range would be
        base_range = [40, 50, 60, 70]  # medium precision base
        if adjusted_suggested > max(base_range):
            extended_range = base_range[:-1] + [round(adjusted_suggested * 1.1, 1)]
            print(f"  Adaptive range: {sorted(extended_range)} (extended for structure)")
        else:
            print(f"  Base range: {base_range} (no extension needed)")
        
        print(f"  Actual ranges: ecutwfc={workflow.ecutwfc_range}, kspacing={workflow.kspacing_range}")
        print()

    print("=== Key Benefits ===")
    print("✓ Pseudopotential-aware: Ranges adjusted based on suggested_ecutwfc")
    print("✓ Structure-aware: Additional scaling for surface/cluster/vacuum systems")
    print("✓ Adaptive precision: Same 'medium' level gives different ranges for different systems")
    print("✓ Automatic optimization: No manual parameter tuning required")

if __name__ == '__main__':
    demonstrate_adaptive_precision()