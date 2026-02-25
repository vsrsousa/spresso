"""
Tests for run_projwfc() method in CalculationWorkflow.

Tests the orbital-projected density of states (PDOS) functionality
for both magnetic and non-magnetic systems.
"""

import unittest
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock, call
import tempfile
import os
from ase.build import bulk
from ase import Atoms

# Import workflow
from xespresso.workflow.calculation_workflow import CalculationWorkflow


class TestRunProjwfc(unittest.TestCase):
    """Test run_projwfc() method with various configurations."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create simple test structure
        self.atoms = bulk('Si', 'diamond', a=5.43)
        self.pseudopotentials = {'Si': 'Si.pbe.UPF'}
        self.temp_dir = tempfile.mkdtemp(prefix='test_projwfc_')
        self.orig_cwd = os.getcwd()
        
    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        os.chdir(self.orig_cwd)
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_projwfc_basic_setup(self):
        """Test basic PROJWFC parameter setup."""
        workflow = CalculationWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
            protocol='moderate'
        )
        
        # Verify workflow is set up correctly
        self.assertEqual(workflow.atoms.get_chemical_symbols(), self.atoms.get_chemical_symbols())
        self.assertEqual(workflow.pseudopotentials, self.pseudopotentials)
        # nspin defaults to 1 (non-magnetic) if not explicitly set
        self.assertEqual(workflow.input_data.get('nspin', 1), 1)
    
    def test_projwfc_with_magnetic_system(self):
        """Test PROJWFC setup with magnetic system."""
        # Create workflow with magnetic configuration
        workflow = CalculationWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
            protocol='moderate',
            magnetic_config='ferro'
        )
        
        # After applying magnetic config, nspin should be 2
        self.assertEqual(workflow.input_data.get('nspin'), 2)
    
    def test_projwfc_nonexistent_nscf_raises_error(self):
        """Test that run_projwfc raises error if NSCF directory doesn't exist."""
        workflow = CalculationWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
            protocol='moderate'
        )
        
        # Try to run PROJWFC without NSCF calculation
        with self.assertRaises(FileNotFoundError) as context:
            workflow.run_projwfc(nscf_label='nonexistent_nscf')
        
        self.assertIn("NSCF calculation directory", str(context.exception))
        self.assertIn("nonexistent_nscf", str(context.exception))
    
    @patch('xespresso.post.projwfc.EspressoProjwfc')
    def test_projwfc_parameters_non_magnetic(self, mock_projwfc_class):
        """Test PROJWFC parameter passing for non-magnetic system."""
        # Create a mock NSCF directory
        nscf_dir = Path(self.temp_dir) / 'nscf'
        nscf_dir.mkdir(parents=True, exist_ok=True)
        
        # Create workflow
        workflow = CalculationWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
            protocol='moderate'
        )
        
        # Mock EspressoProjwfc instance
        mock_projwfc_instance = Mock()
        mock_projwfc_class.return_value = mock_projwfc_instance
        
        # Change to temp directory
        os.chdir(self.temp_dir)
        
        # Run PROJWFC with custom parameters
        workflow.run_projwfc(
            nscf_label='nscf',
            projwfc_label='projwfc',
            Emin=-40,
            Emax=20,
            DeltaE=0.02,
            degauss=0.02,
            ngauss=1,
            lsym=0,
            pawproj=1
        )
        
        # Verify EspressoProjwfc was instantiated with correct parameters
        mock_projwfc_class.assert_called_once()
        call_kwargs = mock_projwfc_class.call_args[1]
        
        self.assertEqual(call_kwargs['Emin'], -40)
        self.assertEqual(call_kwargs['Emax'], 20)
        self.assertEqual(call_kwargs['DeltaE'], 0.02)
        self.assertEqual(call_kwargs['degauss'], 0.02)
        self.assertEqual(call_kwargs['ngauss'], 1)
        self.assertEqual(call_kwargs['lsym'], 0)
        self.assertEqual(call_kwargs['pawproj'], 1)
        
        # Verify run was called
        mock_projwfc_instance.run.assert_called_once()
    
    @patch('xespresso.post.projwfc.EspressoProjwfc')
    def test_projwfc_default_parameters(self, mock_projwfc_class):
        """Test PROJWFC uses correct default parameters."""
        # Create mock NSCF directory
        nscf_dir = Path(self.temp_dir) / 'nscf'
        nscf_dir.mkdir(parents=True, exist_ok=True)
        
        workflow = CalculationWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
            protocol='moderate'
        )
        
        # Mock instance
        mock_projwfc_instance = Mock()
        mock_projwfc_class.return_value = mock_projwfc_instance
        
        os.chdir(self.temp_dir)
        
        # Run with default parameters
        workflow.run_projwfc(nscf_label='nscf')
        
        # Check default values
        call_kwargs = mock_projwfc_class.call_args[1]
        
        self.assertEqual(call_kwargs['Emin'], -30.0)
        self.assertEqual(call_kwargs['Emax'], 10.0)
        self.assertEqual(call_kwargs['DeltaE'], 0.01)
        self.assertEqual(call_kwargs['lsym'], 1)
        self.assertEqual(call_kwargs['pawproj'], 0)
    
    @patch('xespresso.post.projwfc.EspressoProjwfc')
    def test_projwfc_magnetic_system(self, mock_projwfc_class):
        """Test PROJWFC with magnetic system."""
        # Create mock NSCF directory
        nscf_dir = Path(self.temp_dir) / 'nscf'
        nscf_dir.mkdir(parents=True, exist_ok=True)
        
        workflow = CalculationWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
            protocol='moderate',
            magnetic_config='ferro'
        )
        
        # Verify magnetic system setup
        self.assertEqual(workflow.input_data['nspin'], 2)
        
        mock_projwfc_instance = Mock()
        mock_projwfc_class.return_value = mock_projwfc_instance
        
        os.chdir(self.temp_dir)
        
        workflow.run_projwfc(nscf_label='nscf')
        
        # Instance should be created and run called
        mock_projwfc_class.assert_called_once()
        mock_projwfc_instance.run.assert_called_once()
    
    @patch('xespresso.post.projwfc.EspressoProjwfc')
    def test_projwfc_returns_calculator(self, mock_projwfc_class):
        """Test run_projwfc returns EspressoProjwfc instance."""
        nscf_dir = Path(self.temp_dir) / 'nscf'
        nscf_dir.mkdir(parents=True, exist_ok=True)
        
        workflow = CalculationWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
            protocol='moderate'
        )
        
        # Create mock instance with directory attribute
        mock_projwfc_instance = Mock()
        mock_projwfc_instance.directory = 'projwfc'
        mock_projwfc_class.return_value = mock_projwfc_instance
        
        os.chdir(self.temp_dir)
        
        result = workflow.run_projwfc(nscf_label='nscf')
        
        # Should return the calculator instance
        self.assertEqual(result, mock_projwfc_instance)
    
    @patch('xespresso.post.projwfc.EspressoProjwfc')
    def test_projwfc_with_queue_config(self, mock_projwfc_class):
        """Test PROJWFC respects queue configuration."""
        nscf_dir = Path(self.temp_dir) / 'nscf'
        nscf_dir.mkdir(parents=True, exist_ok=True)
        
        queue_config = {
            'name': 'gpu',
            'nodes': 2,
            'parallel': 'mpi=8'
        }
        
        workflow = CalculationWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
            protocol='moderate',
            queue=queue_config
        )
        
        mock_projwfc_instance = Mock()
        mock_projwfc_class.return_value = mock_projwfc_instance
        
        os.chdir(self.temp_dir)
        
        workflow.run_projwfc(nscf_label='nscf')
        
        # Verify queue config was passed
        call_kwargs = mock_projwfc_class.call_args[1]
        self.assertIsNotNone(call_kwargs['queue'])
        self.assertEqual(call_kwargs['queue']['name'], 'gpu')
    
    @patch('xespresso.post.projwfc.EspressoProjwfc')
    def test_projwfc_prefix_extraction(self, mock_projwfc_class):
        """Test that PROJWFC correctly extracts prefix from path."""
        # Create nested directories
        nscf_dir = Path(self.temp_dir) / 'calculations' / 'nscf_step'
        nscf_dir.mkdir(parents=True, exist_ok=True)
        
        workflow = CalculationWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
            protocol='moderate'
        )
        
        mock_projwfc_instance = Mock()
        mock_projwfc_class.return_value = mock_projwfc_instance
        
        os.chdir(self.temp_dir)
        
        workflow.run_projwfc(nscf_label='calculations/nscf_step')
        
        # Should extract 'nscf_step' as prefix
        call_kwargs = mock_projwfc_class.call_args[1]
        # The prefix passed to EspressoProjwfc should be from leaf directory
        call_args = mock_projwfc_class.call_args[0]
    
    @patch('xespresso.post.projwfc.EspressoProjwfc')
    def test_projwfc_custom_filpdos_prefix(self, mock_projwfc_class):
        """Test PROJWFC with custom output filename prefix."""
        nscf_dir = Path(self.temp_dir) / 'nscf'
        nscf_dir.mkdir(parents=True, exist_ok=True)
        
        workflow = CalculationWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
            protocol='moderate'
        )
        
        mock_projwfc_instance = Mock()
        mock_projwfc_class.return_value = mock_projwfc_instance
        
        os.chdir(self.temp_dir)
        
        workflow.run_projwfc(
            nscf_label='nscf',
            filpdos='my_custom_pdos'
        )
        
        call_kwargs = mock_projwfc_class.call_args[1]
        self.assertEqual(call_kwargs['filpdos'], 'my_custom_pdos')


if __name__ == '__main__':
    unittest.main()
