"""
Test that convergence stops immediately when first value converges.
"""

import pytest
import numpy as np
from unittest.mock import Mock, MagicMock, patch
from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow


class TestConvergenceStopsAtFirst:
    """Test that convergence halts when first value converges."""
    
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
    def test_stops_at_first_convergence(
        self, mock_calc_wf_class, si_atoms, pseudo_dict
    ):
        """
        Test that convergence stops when first value converges.
        
        Scenario:
        - ecutwfc < 60: not converged (ΔE > tolerance)
        - ecutwfc = 60: CONVERGED (ΔE < tolerance)
        - ecutwfc > 60: not tested (should stop at 60)
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
        
        def mock_submit_with_save(batch_params, *args, **kwargs):
            """Save batch_params for wait_for_batch_jobs."""
            last_batch_params['batch_params'] = batch_params
            return {'batch_id': 'test_batch_1'}
        
        def mock_wait(batch_results, *args, **kwargs):
            """Mock wait_for_batch_jobs to return completions."""
            batch_params = last_batch_params.get('batch_params', [])
            
            completions = []
            for param in batch_params:
                ecut = param['ecutwfc']
                requested_ecutwfc.append(ecut)
                
                # Reference value
                if ecut == 200:
                    energy = -100.0
                # ecutwfc=60 converges!
                elif ecut == 60:
                    energy = -100.0 + 0.0001  # ΔE = 0.0001 < tolerance(0.003) ✓
                # Others don't converge
                else:
                    energy = -100.0 + 0.05  # ΔE = 0.05 > tolerance(0.003) ✗
                
                completions.append(self.create_mock_completion(ecut, energy))
            
            return completions
        
        mock_wf_instance.submit_scf_batch_multiple = mock_submit_with_save
        mock_wf_instance.wait_for_batch_jobs = mock_wait
        
        # Run convergence
        results = workflow.run_convergence_independent(
            max_ecutwfc=200.0,
            ecutwfc_step=10.0,
            verbose=True
        )
        
        # Verify
        print(f"\n✓ Requested ecutwfc values: {sorted(set(requested_ecutwfc))}")
        
        # Should have tested initial [30, 40, 50] + reference [200]
        # Then expanded to [60] which converges → STOPS
        # Should NOT have gone to [70], [80], etc.
        
        assert 30 in requested_ecutwfc, "Should test initial ecutwfc=30"
        assert 40 in requested_ecutwfc, "Should test initial ecutwfc=40"
        assert 50 in requested_ecutwfc, "Should test initial ecutwfc=50"
        assert 200 in requested_ecutwfc, "Should test reference ecutwfc=200"
        assert 60 in requested_ecutwfc, "Should expand to ecutwfc=60"
        
        # Critical assertion: should NOT go beyond 60
        max_tested = max(e for e in requested_ecutwfc if e < 200)
        assert max_tested == 60, \
            f"Should stop at ecutwfc=60 (first convergence), but tested up to {max_tested}. Full list: {sorted(set(requested_ecutwfc))}"
        
        # Should have selected ecutwfc = 60 (the converged one)
        assert workflow.optimal_ecutwfc == 60, \
            f"Should select ecutwfc=60 (first convergence), got {workflow.optimal_ecutwfc}"
        
        print(f"\n✓✓ SUCCESS: Stops at first convergence (ecutwfc=60), doesn't expand further!")


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
