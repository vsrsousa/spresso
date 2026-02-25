"""Test run_dos() method with spin polarization support."""
import unittest
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
from ase.build import bulk
from xespresso.workflow.calculation_workflow import CalculationWorkflow


class TestRunDos(unittest.TestCase):
    """Test DOS post-processing calculation with spin support."""

    def setUp(self):
        """Set up test fixtures."""
        self.atoms = bulk('Al', 'fcc', a=4.05)
        self.pseudopotentials = {'Al': 'Al.pbe.UPF'}
        
        # Non-magnetic workflow
        self.workflow = CalculationWorkflow(
            self.atoms,
            protocol='fast',
            pseudopotentials=self.pseudopotentials
        )
        
        # Magnetic workflow
        self.atoms_fe = bulk('Fe', 'bcc', a=2.87)
        self.workflow_mag = CalculationWorkflow(
            self.atoms_fe,
            protocol='fast',
            pseudopotentials={'Fe': 'Fe.pbe.UPF'},
            magnetic_config='ferro'
        )

    @patch('xespresso.post.dos.EspressoDos')
    def test_dos_basic_parameters(self, mock_dos_class):
        """Test DOS initialization with default parameters."""
        # Mock EspressoDos instance
        mock_dos_instance = MagicMock()
        mock_dos_instance.directory = 'nscf/dos'
        mock_dos_class.return_value = mock_dos_instance
        
        # Create mock NSCF directory
        Path('nscf').mkdir(exist_ok=True)
        Path('nscf/.keep').touch()
        
        try:
            result = self.workflow.run_dos(nscf_label='nscf')
            
            # Verify EspressoDos was called with correct parameters
            mock_dos_class.assert_called_once()
            call_args = mock_dos_class.call_args
            
            # Check positional/keyword arguments
            self.assertEqual(call_args.kwargs['parent_directory'], 'nscf')
            self.assertEqual(call_args.kwargs['prefix'], 'nscf')
            self.assertEqual(call_args.kwargs['Emin'], -30.0)  # default
            self.assertEqual(call_args.kwargs['Emax'], 10.0)   # default
            self.assertEqual(call_args.kwargs['DeltaE'], 0.01)
            
            # Verify run() was called
            mock_dos_instance.run.assert_called_once()
            
        finally:
            # Cleanup
            import shutil
            if Path('nscf').exists():
                shutil.rmtree('nscf')

    @patch('xespresso.post.dos.EspressoDos')
    def test_dos_custom_energy_windows(self, mock_dos_class):
        """Test DOS with custom energy windows."""
        mock_dos_instance = MagicMock()
        mock_dos_instance.directory = 'nscf/dos'
        mock_dos_class.return_value = mock_dos_instance
        
        Path('nscf').mkdir(exist_ok=True)
        Path('nscf/.keep').touch()
        
        try:
            result = self.workflow.run_dos(
                nscf_label='nscf',
                Emin=-50.0,
                Emax=20.0,
                DeltaE=0.05
            )
            
            call_args = mock_dos_class.call_args
            self.assertEqual(call_args.kwargs['Emin'], -50.0)
            self.assertEqual(call_args.kwargs['Emax'], 20.0)
            self.assertEqual(call_args.kwargs['DeltaE'], 0.05)
            
        finally:
            import shutil
            if Path('nscf').exists():
                shutil.rmtree('nscf')

    @patch('xespresso.post.dos.EspressoDos')
    def test_dos_with_degauss(self, mock_dos_class):
        """Test DOS with custom degauss parameter."""
        mock_dos_instance = MagicMock()
        mock_dos_instance.directory = 'nscf/dos'
        mock_dos_class.return_value = mock_dos_instance
        
        Path('nscf').mkdir(exist_ok=True)
        Path('nscf/.keep').touch()
        
        try:
            result = self.workflow.run_dos(
                nscf_label='nscf',
                degauss=0.05  # Custom from workflow's preset: 0.02 (fast)
            )
            
            call_args = mock_dos_class.call_args
            # Should use custom value, not preset's 0.02
            self.assertEqual(call_args.kwargs['degauss'], 0.05)
            
        finally:
            import shutil
            if Path('nscf').exists():
                shutil.rmtree('nscf')

    def test_dos_missing_nscf_raises_error(self):
        """Test that missing NSCF directory raises FileNotFoundError."""
        with self.assertRaises(FileNotFoundError) as context:
            self.workflow.run_dos(nscf_label='nonexistent_nscf')
        
        self.assertIn('NSCF calculation directory', str(context.exception))

    @patch('xespresso.post.dos.EspressoDos')
    def test_dos_magnetic_system(self, mock_dos_class):
        """Test DOS with magnetic system (nspin=2)."""
        # Verify magnetic workflow has nspin=2
        self.assertTrue(self.workflow_mag.input_data.get('nspin', 1) > 1)
        
        mock_dos_instance = MagicMock()
        mock_dos_instance.directory = 'nscf/dos'
        mock_dos_class.return_value = mock_dos_instance
        
        Path('nscf').mkdir(exist_ok=True)
        Path('nscf/.keep').touch()
        
        try:
            # Call run_dos on magnetic workflow
            result = self.workflow_mag.run_dos(nscf_label='nscf')
            
            # Verify it was called with correct params
            mock_dos_class.assert_called_once()
            mock_dos_instance.run.assert_called_once()
            
        finally:
            import shutil
            if Path('nscf').exists():
                shutil.rmtree('nscf')

    @patch('xespresso.post.dos.EspressoDos')
    def test_dos_with_pdos_flag(self, mock_dos_class):
        """Test DOS with PDOS (projected DOS) flag."""
        mock_dos_instance = MagicMock()
        mock_dos_instance.directory = 'nscf/dos'
        mock_dos_class.return_value = mock_dos_instance
        
        Path('nscf').mkdir(exist_ok=True)
        Path('nscf/.keep').touch()
        
        try:
            result = self.workflow_mag.run_dos(
                nscf_label='nscf',
                pdos=True  # Enable projected DOS
            )
            
            # Verify PDOS flag was handled (logged)
            mock_dos_instance.run.assert_called_once()
            
        finally:
            import shutil
            if Path('nscf').exists():
                shutil.rmtree('nscf')

    @patch('xespresso.post.dos.EspressoDos')
    def test_dos_with_queue_config(self, mock_dos_class):
        """Test DOS with queue configuration."""
        # Create workflow with queue
        queue_config = {
            'execution': 'local',
            'parallel': '-np 4'
        }
        workflow = CalculationWorkflow(
            self.atoms,
            protocol='fast',
            pseudopotentials=self.pseudopotentials,
            queue=queue_config
        )
        
        mock_dos_instance = MagicMock()
        mock_dos_instance.directory = 'nscf/dos'
        mock_dos_class.return_value = mock_dos_instance
        
        Path('nscf').mkdir(exist_ok=True)
        Path('nscf/.keep').touch()
        
        try:
            result = workflow.run_dos(nscf_label='nscf')
            
            call_args = mock_dos_class.call_args
            self.assertEqual(call_args.kwargs['queue'], queue_config)
            self.assertEqual(call_args.kwargs['parallel'], '-np 4')
            
        finally:
            import shutil
            if Path('nscf').exists():
                shutil.rmtree('nscf')


if __name__ == '__main__':
    unittest.main()
