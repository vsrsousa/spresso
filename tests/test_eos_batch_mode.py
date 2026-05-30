#!/usr/bin/env python
"""
Test batch mode detection and job submission with mocked SLURM machine.
This test verifies that when a queue with scheduler='slurm' is provided,
the EOS workflow uses SLURM batch mode (_run_eos_slurm_batch) instead of
parallel worker mode (_run_eos_parallel).
"""
import pytest
import numpy as np
from unittest.mock import patch, MagicMock
from ase.build import bulk
from xespresso.workflow import EOSWorkflow


class TestEOSWorkflowBatchMode:
    """Test SLURM batch mode detection and execution."""
    
    @pytest.fixture
    def atoms_si(self):
        """Create Si atoms object."""
        return bulk('Si', 'diamond', a=5.43)
    
    @pytest.fixture
    def slurm_queue(self):
        """Create mocked SLURM queue dict."""
        return {
            'execution': 'remote',
            'scheduler': 'slurm',
            'remote_host': 'cluster.example.com',
            'remote_user': 'user',
            'remote_dir': '/scratch/xespresso'
        }
    
    def test_batch_mode_detection_with_slurm_queue(self, atoms_si, slurm_queue):
        """Verify that SLURM queue is properly detected and batch mode is enabled."""
        
        eos = EOSWorkflow(
            atoms=atoms_si,
            pseudopotentials_config='default',
            protocol='moderate',
            queue=slurm_queue
        )
        
        # Verify queue is set and contains SLURM scheduler
        assert eos.queue is not None
        assert eos.queue.get('scheduler') == 'slurm'
    
    def test_batch_mode_executes_all_jobs_at_once(self, atoms_si, slurm_queue):
        """Verify that batch mode submits all 11 jobs at once instead of using workers."""
        
        eos = EOSWorkflow(
            atoms=atoms_si,
            pseudopotentials_config='default',
            protocol='moderate',
            queue=slurm_queue
        )
        
        # Mock scaled structures
        scale_factors = np.linspace(0.9, 1.1, 11)
        scaled_structures = {f: atoms_si.copy() for f in scale_factors}
        
        with patch.object(eos, 'create_scaled_structures') as mock_scaled:
            mock_scaled.return_value = scaled_structures
            
            # Mock _run_eos_slurm_batch to return proper results
            batch_call_count = [0]
            
            def mock_batch_method(label):
                batch_call_count[0] += 1
                # Return list of tuples: (factor, energy, volume, error)
                factors = sorted(list(scaled_structures.keys()))
                results = []
                for f in factors:
                    energy = -150.0 - 0.5*(f-1.0)**2  # Parabolic curve
                    volume = scaled_structures[f].get_volume()
                    results.append((f, energy, volume, None))  # None = success
                return results
            
            # Mock _run_eos_parallel to verify it's NOT called
            parallel_call_count = [0]
            
            def mock_parallel_method(label, max_workers=4):
                parallel_call_count[0] += 1
                return []
            
            with patch.object(eos, '_run_eos_slurm_batch', side_effect=mock_batch_method):
                with patch.object(eos, '_run_eos_parallel', side_effect=mock_parallel_method):
                    # Run EOS study
                    df = eos.run_eos_study(
                        volume_range=(0.90, 1.10),
                        n_points=11,
                        label='test_batch_mode'
                    )
                    
                    # Verify batch mode was used
                    assert batch_call_count[0] == 1, "Batch method should be called once"
                    assert parallel_call_count[0] == 0, "Parallel method should NOT be called"
                    
                    # Verify all 11 points were collected
                    assert len(df) == 11, f"Should have 11 points, got {len(df)}"
                    assert df['factor'].min() >= 0.9
                    assert df['factor'].max() <= 1.1
    
    def test_batch_mode_not_used_without_slurm(self, atoms_si):
        """Verify that batch mode is NOT used when no queue is provided."""
        
        eos = EOSWorkflow(
            atoms=atoms_si,
            pseudopotentials_config='default',
            protocol='moderate',
            queue=None  # No queue → no SLURM
        )
        
        # Verify queue is None
        assert eos.queue is None
        
        # When running with no queue, it should use parallel or sequential mode
        # (not batch mode) by default
        scale_factors = np.linspace(0.9, 1.1, 5)
        
        batch_called = [False]
        parallel_called = [False]
        sequential_called = [False]
        
        def mock_batch_method(label):
            batch_called[0] = True
            return []
        
        def mock_parallel_method(label, max_workers=4):
            parallel_called[0] = True
            factors = sorted(list(eos.eos_structures.keys()))
            results = [
                (f, -150.0 - 0.5*(f-1.0)**2, eos.eos_structures[f].get_volume(), None)
                for f in factors
            ]
            return results
        
        with patch.object(eos, '_run_eos_slurm_batch', side_effect=mock_batch_method):
            with patch.object(eos, '_run_eos_parallel', side_effect=mock_parallel_method):
                df = eos.run_eos_study(
                    scale_factors=scale_factors,
                    label='test_no_batch'
                )
                
                # Batch should NOT be called
                assert not batch_called[0], "Batch method should NOT be called when queue is None"
                # Parallel should be called
                assert parallel_called[0], "Parallel method should be called when queue is None"
                # Results should be returned
                assert len(df) == 5



if __name__ == '__main__':
    pytest.main([__file__, '-v'])
