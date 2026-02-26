#!/usr/bin/env python3
"""
Test script for structural complexity analysis in ConvergenceWorkflow
"""

from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow

def test_structural_analysis():
    """Test structural complexity analysis with different structure types."""
    print("Testing structural complexity analysis...")

    # Test 1: Bulk silicon (simple structure, 2 atoms)
    atoms_si = bulk('Si', 'diamond', a=5.43)
    workflow_si = ConvergenceWorkflow(
        atoms=atoms_si,
        pseudopotentials={'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'},
        precision='medium'
    )

    complexity_si = workflow_si._analyze_structural_complexity(atoms_si)
    print(f"Silicon bulk (2 atoms): {complexity_si}")

    # Verify bulk silicon factors
    assert complexity_si['surface_factor'] == 2.0, f"Expected surface_factor=2.0 for small bulk, got {complexity_si['surface_factor']}"
    assert complexity_si['vacuum_factor'] == 1.0, f"Expected vacuum_factor=1.0 for bulk, got {complexity_si['vacuum_factor']}"
    assert complexity_si['heterogeneity_factor'] == 1.0, f"Expected heterogeneity_factor=1.0 for pure element, got {complexity_si['heterogeneity_factor']}"

    # Test 2: Small cluster (should have high surface factor)
    from ase import Atoms
    atoms_cluster = Atoms('Cu4', positions=[[0, 0, 0], [2.5, 0, 0], [0, 2.5, 0], [0, 0, 2.5]])
    workflow_cluster = ConvergenceWorkflow(
        atoms=atoms_cluster,
        pseudopotentials={'Cu': 'Cu.pbe.UPF'},
        precision='medium'
    )

    complexity_cluster = workflow_cluster._analyze_structural_complexity(atoms_cluster)
    print(f"Copper cluster ({len(atoms_cluster)} atoms): {complexity_cluster}")

    # Verify cluster factors
    assert complexity_cluster['surface_factor'] == 2.0, f"Expected surface_factor=2.0 for small cluster, got {complexity_cluster['surface_factor']}"
    assert complexity_cluster['vacuum_factor'] == 1.0, f"Expected vacuum_factor=1.0 for cluster, got {complexity_cluster['vacuum_factor']}"
    assert complexity_cluster['heterogeneity_factor'] == 1.0, f"Expected heterogeneity_factor=1.0 for pure element, got {complexity_cluster['heterogeneity_factor']}"

    # Test 3: Surface with vacuum (should have vacuum factor > 1)
    from ase.build import fcc111
    atoms_surface = fcc111('Al', size=(2, 2, 4), vacuum=10.0)
    workflow_surface = ConvergenceWorkflow(
        atoms=atoms_surface,
        pseudopotentials={'Al': 'Al.pbe.UPF'},
        precision='medium'
    )

    complexity_surface = workflow_surface._analyze_structural_complexity(atoms_surface)
    print(f"Aluminum surface ({len(atoms_surface)} atoms): {complexity_surface}")

    # Verify surface factors
    assert complexity_surface['surface_factor'] == 1.5, f"Expected surface_factor=1.5 for medium system, got {complexity_surface['surface_factor']}"
    assert complexity_surface['vacuum_factor'] > 1.0, f"Expected vacuum_factor>1.0 for surface with vacuum, got {complexity_surface['vacuum_factor']}"
    assert complexity_surface['heterogeneity_factor'] == 1.0, f"Expected heterogeneity_factor=1.0 for pure element, got {complexity_surface['heterogeneity_factor']}"

    # Test 4: Multi-element system (should have higher heterogeneity factor)
    atoms_compound = bulk('GaAs', 'zincblende', a=5.65)
    workflow_compound = ConvergenceWorkflow(
        atoms=atoms_compound,
        pseudopotentials={'Ga': 'Ga.pbe.UPF', 'As': 'As.pbe.UPF'},
        precision='medium'
    )

    complexity_compound = workflow_compound._analyze_structural_complexity(atoms_compound)
    print(f"GaAs compound ({len(atoms_compound)} atoms): {complexity_compound}")

    # Verify compound factors
    assert complexity_compound['heterogeneity_factor'] == 1.2, f"Expected heterogeneity_factor=1.2 for binary compound, got {complexity_compound['heterogeneity_factor']}"

    print("✓ All structural analysis tests passed!")
    print("✓ Complexity factors are calculated correctly for different structure types")

if __name__ == '__main__':
    test_structural_analysis()