"""Test run_bands() method with automatic bandpath generation."""
import unittest
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
from ase.build import bulk
from xespresso.workflow.calculation_workflow import CalculationWorkflow


class TestRunBands(unittest.TestCase):
    """Test band structure calculation with high-symmetry k-paths."""

    def setUp(self):
        """Set up test fixtures."""
        self.atoms_al = bulk('Al', 'fcc', a=4.05)
        self.pseudopotentials_al = {'Al': 'Al.pbe.UPF'}
        
        # Non-magnetic workflow
        self.workflow_al = CalculationWorkflow(
            self.atoms_al,
            protocol='fast',
            pseudopotentials=self.pseudopotentials_al
        )
        
        # Magnetic workflow
        self.atoms_fe = bulk('Fe', 'bcc', a=2.87)
        self.workflow_fe = CalculationWorkflow(
            self.atoms_fe,
            protocol='fast',
            pseudopotentials={'Fe': 'Fe.pbe.UPF'},
            magnetic_config='ferro'
        )

    @patch('xespresso.workflow.calculation_workflow.Espresso')
    def test_bands_auto_bandpath_generation(self, mock_espresso_class):
        """Test automatic high-symmetry bandpath generation."""
        mock_calc = MagicMock()
        mock_calc.last_job_id = '12345'
        mock_espresso_class.return_value = mock_calc
        
        result = self.workflow_al.run_bands(label='bands')
        
        # Verify Espresso was called
        mock_espresso_class.assert_called_once()
        call_args = mock_espresso_class.call_args
        
        # Check calculation type
        self.assertEqual(call_args.kwargs['calculation'], 'bands')
        
        # Check that kpts was a bandpath (has .kpts attribute)
        kpts = call_args.kwargs['kpts']
        self.assertTrue(hasattr(kpts, 'kpts'))  # BandPath object
        self.assertTrue(hasattr(kpts, 'path'))  # Path string
        self.assertTrue(hasattr(kpts, 'special_points'))  # Special k-points dict
        
        # Verify at least some k-points
        self.assertGreater(len(kpts.kpts), 0)
        
        # Verify run was called
        mock_calc.run.assert_called_once()

    @patch('xespresso.workflow.calculation_workflow.Espresso')
    def test_bands_magnetic_system(self, mock_espresso_class):
        """Test band structure for magnetic system (nspin=2)."""
        self.assertTrue(self.workflow_fe.input_data.get('nspin', 1) > 1)
        
        mock_calc = MagicMock()
        mock_calc.last_job_id = '12345'
        mock_espresso_class.return_value = mock_calc
        
        result = self.workflow_fe.run_bands(label='bands_fe')
        
        # Verify it was called with magnetic config
        mock_espresso_class.assert_called_once()
        call_args = mock_espresso_class.call_args
        
        # nspin should be preserved in input_data
        input_data = call_args.kwargs['input_data']
        self.assertEqual(input_data.get('nspin', 1), 2)
        
        # Verify magnetization is set
        self.assertTrue(input_data.get('nspin', 1) > 1)

    @patch('xespresso.workflow.calculation_workflow.Espresso')
    def test_bands_parameters_consistency(self, mock_espresso_class):
        """Test that bandpath parameters are consistent with protocol."""
        mock_calc = MagicMock()
        mock_calc.last_job_id = '12345'
        mock_espresso_class.return_value = mock_calc
        
        result = self.workflow_al.run_bands(label='bands')
        
        call_args = mock_espresso_class.call_args
        
        # Check ecutwfc and ecutrho are set
        self.assertEqual(call_args.kwargs['ecutwfc'], 30.0)  # fast preset
        self.assertEqual(call_args.kwargs['ecutrho'], 240.0)  # fast preset
        
        # Check label
        self.assertEqual(call_args.kwargs['label'], 'bands')

    @patch('xespresso.workflow.calculation_workflow.Espresso')
    def test_bands_custom_label(self, mock_espresso_class):
        """Test bands with custom label."""
        mock_calc = MagicMock()
        mock_calc.last_job_id = '12345'
        mock_espresso_class.return_value = mock_calc
        
        result = self.workflow_al.run_bands(label='custom_bands_path')
        
        call_args = mock_espresso_class.call_args
        self.assertEqual(call_args.kwargs['label'], 'custom_bands_path')

    @patch('xespresso.workflow.calculation_workflow.Espresso')
    def test_bands_with_queue_config(self, mock_espresso_class):
        """Test bands with queue configuration."""
        queue_config = {
            'execution': 'remote',
            'wait_for_completion': True,
            'job_timeout': 3600
        }
        
        workflow = CalculationWorkflow(
            self.atoms_al,
            protocol='fast',
            pseudopotentials=self.pseudopotentials_al,
            queue=queue_config
        )
        
        mock_calc = MagicMock()
        mock_calc.last_job_id = '12345'
        mock_espresso_class.return_value = mock_calc
        
        result = workflow.run_bands(label='bands')
        
        call_args = mock_espresso_class.call_args
        self.assertEqual(call_args.kwargs['queue'], queue_config)

    def test_bands_invalid_bandpath_type_raises_error(self):
        """Test that invalid bandpath_type raises NotImplementedError."""
        with self.assertRaises(NotImplementedError) as context:
            self.workflow_al.run_bands(
                label='bands',
                bandpath_type='custom'  # Not yet implemented
            )
        
        self.assertIn('not yet implemented', str(context.exception))
        self.assertIn('auto', str(context.exception))

    @patch('xespresso.workflow.calculation_workflow.Espresso')
    def test_bands_bandpath_attributes(self, mock_espresso_class):
        """Test that bandpath has correct attributes."""
        mock_calc = MagicMock()
        mock_calc.last_job_id = '12345'
        mock_espresso_class.return_value = mock_calc
        
        result = self.workflow_al.run_bands(label='bands')
        
        call_args = mock_espresso_class.call_args
        kpts = call_args.kwargs['kpts']
        
        # BandPath should have path string
        self.assertIsInstance(kpts.path, str)
        self.assertTrue(len(kpts.path) > 0)
        
        # Should have special points (like G, X, W, L)
        self.assertIsInstance(kpts.special_points, dict)
        self.assertGreater(len(kpts.special_points), 0)


if __name__ == '__main__':
    unittest.main()
