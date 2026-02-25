"""
Tests for ConvergenceWorkflow - parameter optimization workflow.

Tests the convergence study functionality including:
- Parameter ranges and setup
- Results collection and analysis
- Recommendation generation
"""

import unittest
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import tempfile
import os
import pandas as pd
from ase.build import bulk
from ase import Atoms

from xespresso.workflow.convergence_workflow import ConvergenceWorkflow


class TestConvergenceWorkflow(unittest.TestCase):
    """Test ConvergenceWorkflow class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.atoms = bulk('Si', 'diamond', a=5.43)
        self.pseudopotentials = {'Si': 'Si.pbe.UPF'}
        self.temp_dir = tempfile.mkdtemp(prefix='test_convergence_')
    
    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_initialization_default_ranges(self):
        """Test initialization with default parameter ranges."""
        conv = ConvergenceWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
        )
        
        # Check defaults
        self.assertEqual(len(conv.ecutwfc_range), 5)  # [30, 40, 50, 60, 70]
        self.assertEqual(len(conv.kspacing_range), 4)  # [0.5, 0.3, 0.2, 0.15]
        self.assertEqual(conv.protocol, 'moderate')
        self.assertIsNone(conv.results)
    
    def test_initialization_custom_ranges(self):
        """Test initialization with custom parameter ranges."""
        ecutwfc = [50, 60, 70]
        kspacing = [0.3, 0.2]
        
        conv = ConvergenceWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
            ecutwfc_range=ecutwfc,
            kspacing_range=kspacing,
        )
        
        self.assertEqual(conv.ecutwfc_range, ecutwfc)
        self.assertEqual(conv.kspacing_range, kspacing)
    
    def test_ranges_are_sorted(self):
        """Test that parameter ranges are sorted correctly."""
        conv = ConvergenceWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
            ecutwfc_range=[70, 40, 60],  # Unsorted
            kspacing_range=[0.2, 0.5, 0.15],  # Unsorted
        )
        
        # ecutwfc should be sorted ascending
        self.assertEqual(conv.ecutwfc_range, [40, 60, 70])
        # kspacing should be sorted descending (for convenience)
        self.assertEqual(conv.kspacing_range, [0.5, 0.2, 0.15])
    
    def test_from_cif_classmethod(self):
        """Test creating workflow from CIF file."""
        # Create a temporary CIF file
        cif_path = Path(self.temp_dir) / 'test.cif'
        self.atoms.write(str(cif_path), format='cif')
        
        conv = ConvergenceWorkflow.from_cif(
            str(cif_path),
            pseudopotentials=self.pseudopotentials,
        )
        
        self.assertEqual(conv.atoms.get_chemical_formula(), self.atoms.get_chemical_formula())
        self.assertEqual(conv.pseudopotentials, self.pseudopotentials)
    
    def test_results_dataframe_structure(self):
        """Test mock results DataFrame structure."""
        conv = ConvergenceWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
        )
        
        # Create mock results
        mock_results = pd.DataFrame({
            'ecutwfc': [40, 40, 50, 50],
            'kspacing': [0.3, 0.2, 0.3, 0.2],
            'energy': [-100.5, -100.6, -100.6, -100.61],
            'energy_per_atom': [-100.5, -100.6, -100.6, -100.61],
            'max_force': [0.05, 0.02, 0.02, 0.01],
            'n_kpoints': [8, 16, 8, 16],
            'label': ['conv/ecut40_ksp0.30', 'conv/ecut40_ksp0.20', 
                     'conv/ecut50_ksp0.30', 'conv/ecut50_ksp0.20'],
        })
        
        conv.results = mock_results
        
        # Verify structure
        self.assertEqual(len(conv.results), 4)
        expected_columns = {'ecutwfc', 'kspacing', 'energy', 'energy_per_atom',
                           'max_force', 'n_kpoints', 'label'}
        self.assertTrue(expected_columns.issubset(set(conv.results.columns)))
    
    def test_recommendations_without_results_raises_error(self):
        """Test that getting recommendations without results raises error."""
        conv = ConvergenceWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
        )
        
        with self.assertRaises(ValueError) as context:
            conv.get_recommendations()
        
        self.assertIn("No convergence results", str(context.exception))
    
    def test_get_recommendations_with_mock_data(self):
        """Test recommendation generation with mock data."""
        conv = ConvergenceWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
            ecutwfc_range=[40, 50, 60],
            kspacing_range=[0.3, 0.2, 0.15],
        )
        
        # Create mock results with clear convergence pattern
        mock_results = pd.DataFrame({
            'ecutwfc': [40, 40, 40, 50, 50, 50, 60, 60, 60],
            'kspacing': [0.3, 0.2, 0.15, 0.3, 0.2, 0.15, 0.3, 0.2, 0.15],
            'energy': [-100.4, -100.55, -100.59, -100.55, -100.595, -100.605,
                      -100.55, -100.595, -100.610],
            'energy_per_atom': [-100.4, -100.55, -100.59, -100.55, -100.595, -100.605,
                               -100.55, -100.595, -100.610],
            'max_force': [0.1, 0.05, 0.02, 0.05, 0.02, 0.01, 0.05, 0.02, 0.01],
            'n_kpoints': [4, 8, 16, 4, 8, 16, 4, 8, 16],
            'label': [f'conv/ecut{e}_ksp{k:.2f}' 
                     for e, k in zip([40]*3 + [50]*3 + [60]*3,
                                    [0.3, 0.2, 0.15]*3)],
        })
        
        conv.results = mock_results
        
        # Get recommendations with loose tolerance to get all options
        recs = conv.get_recommendations(
            energy_tolerance=1.0,  # Very loose tolerance (1 eV/atom)
            force_tolerance=0.5,   # Very loose tolerance (0.5 eV/Å)
            verbose=False
        )
        
        # Verify recommendation structure
        self.assertIn('fast', recs)
        self.assertIn('accurate', recs)
        self.assertIn('convergence_summary', recs)
        
        # Verify each recommendation has required keys
        for key in ['fast', 'accurate']:
            if key in recs and recs[key] is not None:
                rec = recs[key]
                self.assertIn('ecutwfc', rec)
                self.assertIn('kspacing', rec)
                self.assertIn('energy_per_atom', rec)
                self.assertIn('n_kpoints', rec)
    
    def test_convergence_summary_accuracy(self):
        """Test convergence summary statistics."""
        conv = ConvergenceWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
            ecutwfc_range=[40, 50, 60],
            kspacing_range=[0.3, 0.15],
        )
        
        mock_results = pd.DataFrame({
            'ecutwfc': [40, 40, 50, 50, 60, 60],
            'kspacing': [0.3, 0.15, 0.3, 0.15, 0.3, 0.15],
            'energy': [-100.4, -100.59, -100.55, -100.605, -100.55, -100.605],
            'energy_per_atom': [-100.4, -100.59, -100.55, -100.605, -100.55, -100.605],
            'max_force': [0.1, 0.02, 0.05, 0.01, 0.05, 0.01],
            'n_kpoints': [4, 16, 4, 16, 4, 16],
            'label': ['c/e40_k0.30', 'c/e40_k0.15', 'c/e50_k0.30', 
                     'c/e50_k0.15', 'c/e60_k0.30', 'c/e60_k0.15'],
        })
        
        conv.results = mock_results
        recs = conv.get_recommendations(
            energy_tolerance=1e-3,
            force_tolerance=0.05,
            verbose=False
        )
        
        summary = recs['convergence_summary']
        self.assertEqual(summary['n_total_tests'], 6)
        self.assertEqual(summary['energy_tolerance_eV_atom'], 1e-3)
        self.assertEqual(summary['force_tolerance_eV_A'], 0.05)
    
    def test_to_csv_export(self):
        """Test exporting results to CSV."""
        conv = ConvergenceWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
        )
        
        mock_results = pd.DataFrame({
            'ecutwfc': [40, 50],
            'kspacing': [0.3, 0.2],
            'energy': [-100.4, -100.59],
            'energy_per_atom': [-100.4, -100.59],
            'max_force': [0.1, 0.02],
            'n_kpoints': [4, 8],
            'label': ['c/e40_k0.30', 'c/e50_k0.20'],
        })
        
        conv.results = mock_results
        
        csv_path = Path(self.temp_dir) / 'convergence.csv'
        conv.to_csv(str(csv_path))
        
        self.assertTrue(csv_path.exists())
        loaded = pd.read_csv(csv_path)
        self.assertEqual(len(loaded), 2)
        self.assertIn('ecutwfc', loaded.columns)
    
    def test_from_csv_import(self):
        """Test importing results from CSV."""
        conv = ConvergenceWorkflow(
            atoms=self.atoms,
            pseudopotentials=self.pseudopotentials,
        )
        
        # Create and save CSV
        mock_results = pd.DataFrame({
            'ecutwfc': [40, 50],
            'kspacing': [0.3, 0.2],
            'energy': [-100.4, -100.59],
            'energy_per_atom': [-100.4, -100.59],
            'max_force': [0.1, 0.02],
            'n_kpoints': [4, 8],
            'label': ['c/e40_k0.30', 'c/e50_k0.20'],
        })
        
        csv_path = Path(self.temp_dir) / 'convergence.csv'
        mock_results.to_csv(csv_path, index=False)
        
        # Load from CSV
        conv.from_csv(str(csv_path))
        
        self.assertIsNotNone(conv.results)
        self.assertEqual(len(conv.results), 2)
        self.assertListEqual(list(conv.results['ecutwfc']), [40, 50])


if __name__ == '__main__':
    unittest.main()
