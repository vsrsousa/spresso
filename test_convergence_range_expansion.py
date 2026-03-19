"""
Test convergence workflow range expansion logic with mocking.

Tests that the ecutwfc and kspacing ranges expand correctly,
especially that the reference values don't block expansion.
"""

import pytest
import numpy as np
from unittest.mock import Mock, MagicMock, patch, call
from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow


class TestConvergenceRangeExpansion:
    """Test suite for convergence range expansion logic."""
    
    @pytest.fixture
    def si_atoms(self):
        """Create a simple Si bulk structure."""
        return bulk('Si', 'diamond', a=5.43)
    
    @pytest.fixture
    def pseudo_dict(self):
        """Dummy pseudopotentials dict."""
        return {'Si': '/path/to/Si.pbe.UPF'}
    
    def create_mock_completion(self, ecutwfc, energy, forces=None):
        """Create a mock completion dict with given parameters."""
        forces_arr = forces if forces is not None else np.zeros((len(bulk('Si', 'diamond', a=5.43)), 3))
        return {
            'success': True,
            'energy': energy,
            'forces': forces_arr,
        }
    
    @patch('xespresso.workflow.convergence_workflow.CalculationWorkflow')
    def test_ecutwfc_range_expands_toward_200_without_blocking_at_200(
        self, mock_calc_wf_class, si_atoms, pseudo_dict
    ):
        """
        Test that ecutwfc range can expand up to 190 (one step below reference 200).
        
        This was the bug: range contained 200, blocking expansion because
        max(current_range) >= expansion_limit.
        
        Now: reference (200) is added to batch only, not to current_range.
        So expansion should work: [30, 40, 50] → [30, 40, 50, 60, 70, ...] → [60, 70, 80, ...] etc.
        """
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low',
            convergence_criteria_list=['energy']
        )
        
        # Mock the CalculationWorkflow
        mock_wf_instance = MagicMock()
        mock_calc_wf_class.return_value = mock_wf_instance
        
        # Track which ecutwfc values are requested
        requested_ecutwfc = []
        last_batch_params = {}
        
        def mock_submit(*args, **kwargs):
            """Mock submit_scf_batch_multiple - just return batch IDs."""
            return {'batch_id': 'test_batch_1'}
        
        def mock_wait(batch_results, *args, **kwargs):
            """Mock wait_for_batch_jobs to return completions."""
            # Get batch_params from the _last_submit which we saved
            batch_params = last_batch_params.get('batch_params', [])
            
            completions = []
            for param in batch_params:
                ecut = param['ecutwfc']
                requested_ecutwfc.append(ecut)
                
                # Return converged energies for ecutwfc >= 60, not converged for < 60
                if ecut == 200:
                    # Reference value
                    energy = -100.0  # Reference energy (per atom)
                elif ecut >= 60:
                    # Converged: very close to reference
                    energy = -100.0 + 0.0001
                else:
                    # Not converged: far from reference
                    energy = -100.0 + 0.05
                
                completions.append(self.create_mock_completion(ecut, energy))
            
            return completions
        
        def mock_submit_with_save(batch_params, *args, **kwargs):
            """Wrapper to save batch_params for wait_for_batch_jobs."""
            last_batch_params['batch_params'] = batch_params
            return mock_submit(batch_params, *args, **kwargs)
        
        mock_wf_instance.submit_scf_batch_multiple = mock_submit_with_save
        mock_wf_instance.wait_for_batch_jobs = mock_wait
        
        # Run convergence
        results = workflow.run_convergence_independent(
            max_ecutwfc=200.0,
            ecutwfc_step=10.0,
            verbose=False
        )
        
        # Verify that we tested values beyond initial range
        print(f"\n✓ Requested ecutwfc values: {sorted(set(requested_ecutwfc))}")
        
        # Should have tested initial [30, 40, 50] plus reference [200]
        assert 30 in requested_ecutwfc, "Should test ecutwfc=30"
        assert 40 in requested_ecutwfc, "Should test ecutwfc=40"
        assert 50 in requested_ecutwfc, "Should test ecutwfc=50"
        assert 200 in requested_ecutwfc, "Should test reference ecutwfc=200"
        
        # Should have expanded incrementally (one value per iteration)
        # With incremental expansion, should have gone through multiple iterations
        # Each iteration adds just ONE new value
        assert 60 in requested_ecutwfc, "Should expand to ecutwfc=60"
        
        # May or may not reach higher values depending on when convergence happens
        # The important thing is that it expands ONE value at a time, not all at once
        
        # Should NOT try to go beyond 190 (one step below reference)
        assert all(e <= 190 for e in requested_ecutwfc if e < 200), \
            f"Should not expand beyond 190, but got: {sorted(set(requested_ecutwfc))}"
        
        print(f"\n✓✓ Range expansion working correctly: one value per iteration, respects limit at 200")
    
    
    def test_expand_range_respects_max_val(self, si_atoms, pseudo_dict):
        """Test that _expand_range respects the max_val limit and adds ONE value at a time."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        # Test ecutwfc expansion
        current = [30, 40, 50]
        step = 10.0
        max_val = 190.0  # One step below reference (200)
        
        expansion = workflow._expand_range(current, step, max_val)
        
        # Should expand with exactly ONE new value
        assert len(expansion) == 1, f"Should generate exactly 1 new value, got {len(expansion)}: {expansion}"
        assert expansion[0] == 60.0, f"Should add next value 60 (from max 50 + step 10), got {expansion}"
        
        # Test incremental expansion
        current2 = [30, 40, 50, 60]
        expansion2 = workflow._expand_range(current2, step, max_val)
        assert expansion2[0] == 70.0, f"Should add 70 next, got {expansion2}"
        
        # Test at limit
        current_at_limit = [30, 40, 50, 180]  # Next would be 190
        expansion3 = workflow._expand_range(current_at_limit, step, max_val)
        assert expansion3[0] == 190.0, f"Should allow up to max_val 190, got {expansion3}"
        
        # Test would exceed limit
        current_exceeds = [30, 40, 50, 190]  # Next would be 200, which exceeds
        expansion4 = workflow._expand_range(current_exceeds, step, max_val)
        assert len(expansion4) == 0, f"Should return empty when next value would exceed max_val, got {expansion4}"


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
