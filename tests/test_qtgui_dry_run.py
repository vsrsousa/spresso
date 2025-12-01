"""
Tests for qtgui dry run functionality.

These tests verify that the qtgui properly generates input files
with atomic positions and cell parameters, and creates job_file.
"""

import pytest
import os
import tempfile
from ase.build import bulk


@pytest.fixture(scope="module")
def mock_qapplication():
    """Create a QApplication for testing widgets."""
    import sys
    from PySide6.QtWidgets import QApplication
    
    # Check if QApplication already exists
    app = QApplication.instance()
    if app is None:
        # Create headless app for testing
        app = QApplication(sys.argv)
    
    yield app


class TestJobSubmissionHelpers:
    """Tests for job submission helper methods."""
    
    def test_build_input_data_basic(self, mock_qapplication):
        """Test _build_input_data creates proper input_data dict."""
        from qtgui.pages.job_submission import JobSubmissionPage
        from qtgui.main_app import SessionState
        
        # Reset singleton
        SessionState._instance = None
        state = SessionState()
        
        page = JobSubmissionPage(state)
        
        config = {
            'calc_type': 'scf',
            'ecutwfc': 50.0,
            'ecutrho': 400.0,
            'occupations': 'smearing',
            'smearing': 'gaussian',
            'degauss': 0.02,
            'conv_thr': 1.0e-8,
        }
        
        input_data = page._build_input_data(config, 'test')
        
        assert 'CONTROL' in input_data
        assert 'SYSTEM' in input_data
        assert 'ELECTRONS' in input_data
        assert input_data['CONTROL']['calculation'] == 'scf'
        assert input_data['SYSTEM']['ecutwfc'] == 50.0
        assert input_data['SYSTEM']['smearing'] == 'gaussian'
        assert input_data['ELECTRONS']['conv_thr'] == 1.0e-8
    
    def test_build_input_data_with_magnetism(self, mock_qapplication):
        """Test _build_input_data with magnetic configuration."""
        from qtgui.pages.job_submission import JobSubmissionPage
        from qtgui.main_app import SessionState
        
        # Reset singleton
        SessionState._instance = None
        state = SessionState()
        
        page = JobSubmissionPage(state)
        
        config = {
            'calc_type': 'scf',
            'ecutwfc': 50.0,
            'ecutrho': 400.0,
            'occupations': 'smearing',
            'conv_thr': 1.0e-8,
            'enable_magnetism': True,
            'magnetic_config': {'Fe': [2.0]},
        }
        
        input_data = page._build_input_data(config, 'test')
        
        assert input_data['SYSTEM']['nspin'] == 2
        assert 'INPUT_NTYP' in input_data
        assert 'starting_magnetization' in input_data['INPUT_NTYP']
        assert input_data['INPUT_NTYP']['starting_magnetization']['Fe'] == 2.0
    
    def test_build_input_data_with_hubbard(self, mock_qapplication):
        """Test _build_input_data with Hubbard U configuration."""
        from qtgui.pages.job_submission import JobSubmissionPage
        from qtgui.main_app import SessionState
        
        # Reset singleton
        SessionState._instance = None
        state = SessionState()
        
        page = JobSubmissionPage(state)
        
        config = {
            'calc_type': 'scf',
            'ecutwfc': 50.0,
            'ecutrho': 400.0,
            'occupations': 'smearing',
            'conv_thr': 1.0e-8,
            'enable_hubbard': True,
            'hubbard_u': {'Fe': 4.0},
        }
        
        input_data = page._build_input_data(config, 'test')
        
        assert input_data['SYSTEM']['lda_plus_u'] == True
        assert 'INPUT_NTYP' in input_data
        assert 'Hubbard_U' in input_data['INPUT_NTYP']
        assert input_data['INPUT_NTYP']['Hubbard_U']['Fe'] == 4.0
    
    def test_create_job_file(self, mock_qapplication):
        """Test _create_job_file creates proper job script."""
        from qtgui.pages.job_submission import JobSubmissionPage
        from qtgui.main_app import SessionState
        
        # Reset singleton
        SessionState._instance = None
        state = SessionState()
        
        page = JobSubmissionPage(state)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            job_file_path = os.path.join(tmpdir, 'job_file')
            config = {
                'calc_type': 'scf',
                'nprocs': 4,
            }
            
            page._create_job_file(job_file_path, 'test', config)
            
            assert os.path.exists(job_file_path)
            
            with open(job_file_path, 'r') as f:
                content = f.read()
            
            assert '#!/bin/bash' in content
            assert 'mpirun -np 4' in content
            assert 'pw.x' in content
            assert 'test.pwi' in content
            assert 'test.pwo' in content
    
    def test_create_job_file_single_proc(self, mock_qapplication):
        """Test _create_job_file with single processor (no mpirun)."""
        from qtgui.pages.job_submission import JobSubmissionPage
        from qtgui.main_app import SessionState
        
        # Reset singleton
        SessionState._instance = None
        state = SessionState()
        
        page = JobSubmissionPage(state)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            job_file_path = os.path.join(tmpdir, 'job_file')
            config = {
                'calc_type': 'scf',
                'nprocs': 1,
            }
            
            page._create_job_file(job_file_path, 'test', config)
            
            with open(job_file_path, 'r') as f:
                content = f.read()
            
            assert 'mpirun' not in content
            assert 'pw.x -in test.pwi > test.pwo' in content


class TestWriteEspressoInIntegration:
    """Tests for write_espresso_in integration with GUI."""
    
    def test_write_espresso_in_includes_structure(self):
        """Test that write_espresso_in properly includes atomic positions and cell."""
        from xespresso.xio import write_espresso_in
        
        atoms = bulk('Fe', 'bcc', a=2.87)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, 'test.pwi')
            
            input_data = {
                'CONTROL': {
                    'calculation': 'scf',
                    'prefix': 'test',
                },
                'SYSTEM': {
                    'ecutwfc': 50.0,
                },
                'ELECTRONS': {},
            }
            
            pseudopotentials = {'Fe': 'Fe.UPF'}
            
            write_espresso_in(
                input_path,
                atoms,
                input_data=input_data,
                pseudopotentials=pseudopotentials,
                kpts=(4, 4, 4)
            )
            
            with open(input_path, 'r') as f:
                content = f.read()
            
            # Verify structure is included
            assert 'ATOMIC_POSITIONS' in content
            assert 'CELL_PARAMETERS' in content
            assert 'Fe' in content
            
            # Verify the file is not just a preview with comments
            assert 'Note: This is a preview' not in content


class TestLabelInversion:
    """Tests for label inversion in Calculation Setup and Workflow Builder."""
    
    def test_calculation_setup_label_format(self, mock_qapplication):
        """Test that Calculation Setup creates label in structure/calc_type format."""
        from qtgui.pages.calculation_setup import CalculationSetupPage
        from qtgui.main_app import SessionState
        from ase.build import bulk
        
        # Reset singleton
        SessionState._instance = None
        state = SessionState()
        
        # Set a structure
        atoms = bulk('Al', 'fcc', a=4.05)
        state['current_structure'] = atoms
        
        page = CalculationSetupPage(state)
        page._update_structure_status()
        
        # Check the label format is structure/calc_type (e.g., Al/scf not scf/Al)
        label = page.label_edit.text()
        assert label == "Al/scf", f"Expected 'Al/scf' but got '{label}'"
    
    def test_job_submission_label_format(self, mock_qapplication):
        """Test that Job Submission creates label in structure/calc_type format."""
        from qtgui.pages.job_submission import JobSubmissionPage
        from qtgui.main_app import SessionState
        from ase.build import bulk
        
        # Reset singleton
        SessionState._instance = None
        state = SessionState()
        
        # Set a structure and config
        atoms = bulk('Fe', 'bcc', a=2.87)
        state['current_structure'] = atoms
        state['workflow_config'] = {'calc_type': 'relax'}
        
        page = JobSubmissionPage(state)
        page._update_dry_run_config()
        
        # Check the label format is structure/calc_type (e.g., Fe/relax not relax/Fe)
        label = page.label_edit.text()
        assert label == "Fe/relax", f"Expected 'Fe/relax' but got '{label}'"
    
    def test_placeholder_text_format(self, mock_qapplication):
        """Test that placeholder text shows structure/calc_type format."""
        from qtgui.pages.calculation_setup import CalculationSetupPage
        from qtgui.pages.job_submission import JobSubmissionPage
        from qtgui.main_app import SessionState
        
        # Reset singleton
        SessionState._instance = None
        state = SessionState()
        
        calc_page = CalculationSetupPage(state)
        job_page = JobSubmissionPage(state)
        
        # Check placeholders show inverted format (structure/calc_type)
        assert "Al/scf" in calc_page.label_edit.placeholderText()
        assert "Al/scf" in job_page.label_edit.placeholderText()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
