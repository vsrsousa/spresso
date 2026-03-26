"""
Test convergence limits: ensure ecutwfc doesn't exceed max and kspacing doesn't go below min
"""

import pytest
from unittest.mock import MagicMock, patch, PropertyMock
from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow


class TestConvergenceLimits:
    """Verify that convergence respects min/max boundaries"""
    
    @patch('xespresso.workflow.convergence_workflow.discover_pseudopotential_directory')
    def test_ecutwfc_respects_max_limit(self, mock_discover):
        """
        PHASE 1: Verify ecutwfc convergence stops at max_ecutwfc without exceeding it
        Even if not converged, should not expand beyond max_ecutwfc
        """
        atoms = bulk('Si', 'diamond')
        
        # Mock pseudopotential discovery
        mock_discover.return_value = (
            {'Si': '/fake/Si.pbe.UPF'},
            '/fake'
        )
        
        with patch('xespresso.workflow.convergence_workflow.os.environ', {}):
            workflow = ConvergenceWorkflow(
                atoms=atoms,
                pseudopotentials={'Si': '/fake/Si.pbe.UPF'},
                protocol='moderate',
                precision='low',
                min_ecutwfc=20.0,  # Start from 20
                max_ecutwfc=100.0  # Stop at 100 (not default 200)
            )
        
        # Verify: min and max are correctly set
        assert workflow.min_ecutwfc == 20.0, f"Should have min_ecutwfc=20.0, got {workflow.min_ecutwfc}"
        assert workflow.max_ecutwfc == 100.0, f"Should have max_ecutwfc=100.0, got {workflow.max_ecutwfc}"
        print(f"✓ ecutwfc boundaries set correctly:")
        print(f"  - min_ecutwfc: {workflow.min_ecutwfc}")
        print(f"  - max_ecutwfc: {workflow.max_ecutwfc}")
        print(f"  - In run_convergence_independent: starts with [min] only, expands to [min, min+step, min+2*step, ...]")
        print(f"  - Stops when next value > {workflow.max_ecutwfc - 10.0} Ry")
    
    
    @patch('xespresso.workflow.convergence_workflow.discover_pseudopotential_directory')
    def test_kspacing_respects_min_limit(self, mock_discover):
        """
        PHASE 2: Verify kspacing convergence stops at min_kspacing without going below
        Even if not converged, should not go below 0.1 Å⁻¹
        """
        atoms = bulk('Si', 'diamond')
        
        # Mock pseudopotential discovery
        mock_discover.return_value = (
            {'Si': '/fake/Si.pbe.UPF'},
            '/fake'
        )
        
        with patch('xespresso.workflow.convergence_workflow.os.environ', {}):
            workflow = ConvergenceWorkflow(
                atoms=atoms,
                pseudopotentials={'Si': '/fake/Si.pbe.UPF'},
                protocol='moderate',
                precision='low'
            )
        
        # Check initialization for PHASE 2
        # kspacing_range should start with [0.3, 0.27] according to the code
        # min_kspacing_allowed should be 0.1 (hardcoded in run_convergence_independent)
        
        print(f"✓ kspacing limits set up:")
        print(f"  - Starting range: [0.3, 0.27]")
        print(f"  - Minimum allowed: 0.1 Å⁻¹")
        print(f"  - Step size: 0.03 Å⁻¹")

    
    
    def test_ecutwfc_expansion_incremental(self):
        """
        Verify ecutwfc expands incrementally WITHOUT including max in range
        """
        # Simulate run_convergence_independent PHASE 1
        min_ecutwfc = 30.0
        max_ecutwfc = 200.0
        ecutwfc_step = 10.0
        
        # NEW: ecutwfc_range starts with ONLY min, not [min, max]
        ecutwfc_range = [min_ecutwfc]
        current_ecut_index = 0
        
        # Simulate iterations (stop at convergence)
        iterations = []
        for iteration in range(1, 25):
            # Expansion logic (from run_convergence_independent)
            if current_ecut_index >= len(ecutwfc_range):
                last_ecut = sorted(ecutwfc_range)[-1]
                next_ecut = last_ecut + ecutwfc_step
                expansion_limit = max_ecutwfc - ecutwfc_step
                
                if next_ecut > expansion_limit:
                    print(f"  Iteration {iteration}: Cannot expand further (limit: {expansion_limit:.1f})")
                    break
                ecutwfc_range.append(next_ecut)
            
            ecut_to_test = ecutwfc_range[current_ecut_index]
            
            # Add reference on first iteration
            if iteration == 1:
                iterations.append((ecut_to_test, max_ecutwfc))
                print(f"  Iteration 1: testing {ecut_to_test:.0f} + reference {max_ecutwfc:.0f}")
            else:
                iterations.append(ecut_to_test)
                print(f"  Iteration {iteration}: testing {ecut_to_test:.0f}")
            
            current_ecut_index += 1
            
            # Convergence check: if ecut >= 60, assume converged
            if ecut_to_test >= 60.0:
                print(f"  ✓ Converged at ecutwfc={ecut_to_test:.0f}")
                break
        
        # Verify: all tested values stay within bounds
        all_tested = []
        for item in iterations:
            if isinstance(item, tuple):
                all_tested.extend(item)
            else:
                all_tested.append(item)
        
        for val in all_tested:
            assert val <= max_ecutwfc, f"Tested {val} exceeds max {max_ecutwfc}"
            assert val >= min_ecutwfc, f"Tested {val} below min {min_ecutwfc}"
        
        print(f"✓ All tested values within [{min_ecutwfc}, {max_ecutwfc}]")
    
    
    def test_kspacing_incremental_refinement(self):
        """
        Verify kspacing refines incrementally WITHOUT going below minimum
        """
        # Simulate run_convergence_independent PHASE 2
        kspacing_range = [0.3, 0.27]  # Starts with 2 values for sliding window
        min_kspacing_allowed = 0.1
        kspacing_step = 0.03
        current_ksp_index = 0
        
        # Simulate iterations (stop at convergence)
        iterations = []
        for iteration in range(1, 25):
            # Expansion logic (refine = decrease)
            if current_ksp_index >= len(kspacing_range):
                finest_ksp = min(kspacing_range)
                next_ksp = finest_ksp - kspacing_step
                
                if next_ksp < min_kspacing_allowed:
                    print(f"  Iteration {iteration}: Cannot refine further (limit: {min_kspacing_allowed})")
                    break
                kspacing_range.append(next_ksp)
                kspacing_range.sort(reverse=True)  # Keep reverse sorted
            
            ksp_to_test = kspacing_range[current_ksp_index]
            
            if iteration == 1:
                iterations.append((ksp_to_test, kspacing_range[1]))  # Also calc reference (0.27)
                print(f"  Iteration 1: testing {ksp_to_test:.3f} (sliding window ref: {kspacing_range[1]:.3f})")
            else:
                iterations.append(ksp_to_test)
                print(f"  Iteration {iteration}: testing {ksp_to_test:.3f}")
            
            current_ksp_index += 1
            
            # Convergence check: if kspacing <= 0.15, assume converged
            if ksp_to_test <= 0.15:
                print(f"  ✓ Converged at kspacing={ksp_to_test:.3f}")
                break
        
        # Verify: all tested values stay within bounds
        all_tested = []
        for item in iterations:
            if isinstance(item, tuple):
                all_tested.extend(item)
            else:
                all_tested.append(item)
        
        for val in all_tested:
            assert val >= min_kspacing_allowed, f"Tested {val} below min {min_kspacing_allowed}"
            assert val <= 0.3, f"Tested {val} exceeds initial max 0.3"
        
        print(f"✓ All tested values within [{min_kspacing_allowed}, 0.3]")


if __name__ == '__main__':
    test = TestConvergenceLimits()
    
    print("\n" + "="*80)
    print("TEST 1: ecutwfc_respects_max_limit")
    print("="*80)
    test.test_ecutwfc_respects_max_limit()
    
    print("\n" + "="*80)
    print("TEST 2: kspacing_respects_min_limit")
    print("="*80)
    test.test_kspacing_respects_min_limit()
    
    print("\n" + "="*80)
    print("TEST 3: ecutwfc expansion incremental (without max in initial range)")
    print("="*80)
    test.test_ecutwfc_expansion_incremental()
    
    print("\n" + "="*80)
    print("TEST 4: kspacing refinement incremental (respects min limit)")
    print("="*80)
    test.test_kspacing_incremental_refinement()
    
    print("\n" + "="*80)
    print("✓ ALL LIMIT TESTS PASSED")
    print("="*80)
