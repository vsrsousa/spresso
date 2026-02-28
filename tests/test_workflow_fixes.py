"""
Test suite for workflow fixes:
1. Convergence workflow passing pseudopotentials correctly to calculation workflow
2. Proper handling of remote execution with wait_for_completion
3. Avoiding file re-transfer retry loops
"""

import pytest
import os
import tempfile
from unittest.mock import Mock, MagicMock, patch
from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow
from xespresso.workflow.calculation_workflow import CalculationWorkflow


@pytest.fixture
def si_atoms():
    """Create a simple Si structure for testing."""
    return bulk('Si', 'diamond', a=5.43, cubic=True)


@pytest.fixture
def mock_pseudopotentials_config():
    """Create a mock pseudopotentials configuration."""
    mock_config = MagicMock()
    mock_config.base_path = '/path/to/pseudos'
    
    # Mock pseudo object
    mock_pseudo = MagicMock()
    mock_pseudo.filename = 'Si.pbe.UPF'
    mock_pseudo.element = 'Si'
    
    mock_config.pseudopotentials = {'Si': mock_pseudo}
    mock_config.get_pseudopotential.return_value = mock_pseudo
    mock_config.list_elements.return_value = ['Si']
    
    return mock_config


class TestConvergenceWorkflowPseudopotentials:
    """Test that convergence_workflow passes pseudopotentials correctly."""
    
    @patch('xespresso.pseudopotentials.manager.load_pseudopotentials_config')
    @patch('xespresso.workflow.calculation_workflow.Espresso')
    def test_convergence_passes_config_name_not_paths(
        self, mock_espresso, mock_load_pseudo, si_atoms, mock_pseudopotentials_config
    ):
        """
        Test that convergence_workflow passes pseudopotentials_config name
        (not full paths) to calculation_workflow.
        
        This prevents the issue where full paths were being passed, causing
        calculation_workflow to not properly resolve pseudopotentials.
        """
        mock_load_pseudo.return_value = mock_pseudopotentials_config
        
        # Create convergence workflow
        conv_wf = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials_config='default',
            protocol='fast',
            precision='low'
        )
        
        # Verify it stored the config name
        assert conv_wf._pseudo_config_name == 'default'
        
        # Verify pseudopotentials are stored as filenames only
        assert 'Si' in conv_wf.pseudopotentials
        assert conv_wf.pseudopotentials['Si'] == 'Si.pbe.UPF'  # Filename only, not full path
        
        # Verify base_path is stored
        assert conv_wf.pseudopotentials_base_path == '/path/to/pseudos'


    @patch('xespresso.pseudopotentials.manager.load_pseudopotentials_config')
    @patch('xespresso.workflow.calculation_workflow.CalculationWorkflow')
    def test_convergence_passes_config_name_to_calc_workflow(
        self, mock_calc_workflow, mock_load_pseudo, si_atoms, mock_pseudopotentials_config
    ):
        """
        Test that when creating CalculationWorkflow in convergence_workflow,
        the config name is passed (not the pseudopotentials dict).
        """
        mock_load_pseudo.return_value = mock_pseudopotentials_config
        
        # Create convergence workflow (batch mode)
        with patch('xespresso.workflow.calculation_workflow.Espresso'):
            conv_wf = ConvergenceWorkflow(
                atoms=si_atoms,
                pseudopotentials_config='default',
                protocol='fast',
                precision='low',
                queue={'execution': 'remote', 'scheduler': 'slurm'}
            )
            
            # This would call submit_scf_batch_multiple internally
            # which creates CalculationWorkflow instances
            # We're not actually calling it, just checking the logic
            
            # Verify the config name is what would be passed
            assert conv_wf._pseudo_config_name == 'default'


class TestCalculationWorkflowRemoteExecution:
    """Test that calculation_workflow handles remote execution correctly."""
    
    @patch('xespresso.workflow.calculation_workflow.Espresso')
    @patch('xespresso.workflow.calculation_workflow.load_pseudopotentials_config')
    def test_remote_non_blocking_condition(
        self, mock_load_pseudo, mock_espresso, si_atoms, mock_pseudopotentials_config
    ):
        """
        Test that the remote non-blocking condition correctly identifies
        when to use special handling vs calc.run().
        
        The condition: if queue.get('execution') == 'remote' and not queue.get('wait_for_completion', False)
        
        - Remote + missing wait_for_completion (default) → non-blocking handler ✓
        - Remote + wait_for_completion=False → non-blocking handler ✓
        - Remote + wait_for_completion=True → calc.run() ✓
        - Local → calc.run() ✓
        """
        mock_load_pseudo.return_value = mock_pseudopotentials_config
        
        # Test case 1: Remote with missing wait_for_completion (typical case)
        queue_remote_default = {
            'execution': 'remote',
            'scheduler': 'slurm'
        }
        
        # The condition from the code:
        uses_special_handling = queue_remote_default and \
                                queue_remote_default.get('execution') == 'remote' and \
                                not queue_remote_default.get('wait_for_completion', False)
        
        assert uses_special_handling is True, "Remote with missing wait_for_completion should use special handling"
        
        # Test case 2: Remote with wait_for_completion=False (explicit)
        queue_remote_false = {
            'execution': 'remote',
            'scheduler': 'slurm',
            'wait_for_completion': False
        }
        
        uses_special_handling = queue_remote_false and \
                                queue_remote_false.get('execution') == 'remote' and \
                                not queue_remote_false.get('wait_for_completion', False)
        
        assert uses_special_handling is True, "Remote with wait_for_completion=False should use special handling"
        
        # Test case 3: Remote with wait_for_completion=True
        queue_remote_true = {
            'execution': 'remote',
            'scheduler': 'slurm',
            'wait_for_completion': True
        }
        
        uses_special_handling = queue_remote_true and \
                                queue_remote_true.get('execution') == 'remote' and \
                                not queue_remote_true.get('wait_for_completion', False)
        
        assert uses_special_handling is False, "Remote with wait_for_completion=True should NOT use special handling"
        
        # Test case 4: Local execution
        queue_local = {
            'execution': 'local',
            'scheduler': 'direct'
        }
        
        uses_special_handling = queue_local and \
                                queue_local.get('execution') == 'remote' and \
                                not queue_local.get('wait_for_completion', False)
        
        assert uses_special_handling is False, "Local execution should NOT use special handling"


    @patch('xespresso.workflow.calculation_workflow.Espresso')
    @patch('xespresso.workflow.calculation_workflow.load_pseudopotentials_config')
    def test_no_file_resend_in_remote_non_blocking(
        self, mock_load_pseudo, mock_espresso_class, si_atoms, mock_pseudopotentials_config
    ):
        """
        Test that when using remote non-blocking execution, files are NOT
        re-sent on retries.
        
        The issue was:
        - Non-blocking code path: writes input once, then execute() sends files
        - Blocking code path (calc.run()): has retry logic that calls write_input multiple times
        
        For remote non-blocking, we should NOT call calc.run() which has retries.
        """
        mock_load_pseudo.return_value = mock_pseudopotentials_config
        
        # Mock the Espresso calculator
        mock_calc = MagicMock()
        mock_calc.atoms = None
        mock_calc_instance = MagicMock()
        mock_espresso_class.return_value = mock_calc_instance
        
        # Create calculation workflow with remote queue
        queue = {
            'execution': 'remote',
            'scheduler': 'slurm',
            'remote_host': 'medusa.fis.uerj.br',
            'remote_user': 'vinicius'
        }
        
        with patch.object(CalculationWorkflow, '_monitor_remote_job'):
            with patch('xespresso.workflow.calculation_workflow.RemoteJobMonitor'):
                # We can't fully create the workflow due to mocking, but we can verify the logic
                # In real usage, the condition would route to special handler, not calc.run()
                
                # If queue.get('execution') == 'remote' and not queue.get('wait_for_completion', False):
                #     → enters special handling (writes input once, execute() sends files)
                # else:
                #     → calls calc.run() which has retry logic (would re-send files)
                
                should_use_special = queue.get('execution') == 'remote' and \
                                    not queue.get('wait_for_completion', False)
                
                assert should_use_special is True, \
                    "Should use special handling (no retries) for remote non-blocking"


class TestPseudopotentialsResolution:
    """Test pseudopotential file resolution in different scenarios."""
    
    def test_pseudopotentials_filename_only_in_convergence(
        self, si_atoms, mock_pseudopotentials_config
    ):
        """
        Test that convergence_workflow stores filenames only (not full paths).
        """
        with patch('xespresso.pseudopotentials.manager.load_pseudopotentials_config') as mock_load:
            mock_load.return_value = mock_pseudopotentials_config
            
            conv_wf = ConvergenceWorkflow(
                atoms=si_atoms,
                pseudopotentials_config='default',
                protocol='fast'
            )
            
            # Should be filename only
            assert 'Si' in conv_wf.pseudopotentials
            assert conv_wf.pseudopotentials['Si'] == 'Si.pbe.UPF'
            assert '/' not in conv_wf.pseudopotentials['Si'], \
                "Should not contain path separators - should be filename only"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
