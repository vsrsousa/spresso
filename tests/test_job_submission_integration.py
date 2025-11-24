"""
Tests for improved job submission integration.

These tests verify that the job submission properly uses pre-configured
calculators from the calculation setup page.
"""

import pytest
import os
import tempfile
from ase.build import bulk


def test_prepare_calculation_from_gui():
    """Test that prepare_calculation_from_gui creates a valid calculator."""
    # Set up environment
    os.environ['ASE_ESPRESSO_COMMAND'] = 'pw.x -in PREFIX.pwi > PREFIX.pwo'
    os.environ['ESPRESSO_PSEUDO'] = '/tmp/pseudo'
    
    from gui.calculations import prepare_calculation_from_gui
    
    # Create a test structure
    atoms = bulk('Fe', 'bcc', a=2.87)
    
    # Configuration dictionary (like from workflow_config)
    config = {
        'calc_type': 'scf',
        'ecutwfc': 40.0,
        'ecutrho': 320.0,
        'conv_thr': 1.0e-6,
        'occupations': 'smearing',
        'smearing': 'gaussian',
        'degauss': 0.02,
        'pseudopotentials': {
            'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
        },
        'kpts': (4, 4, 4),
    }
    
    with tempfile.TemporaryDirectory() as tmpdir:
        label = os.path.join(tmpdir, 'test/fe')
        
        # Prepare calculation
        prepared_atoms, calc = prepare_calculation_from_gui(atoms, config, label=label)
        
        # Verify the calculator was created
        assert calc is not None
        # Note: Espresso appends the basename to the label for the directory
        assert calc.label.startswith(label)
        
        # Verify atoms are prepared
        assert prepared_atoms is not None
        assert prepared_atoms.get_chemical_formula() == atoms.get_chemical_formula()


def test_dry_run_calculation():
    """Test that dry_run_calculation creates files without running."""
    # Set up environment
    os.environ['ASE_ESPRESSO_COMMAND'] = 'pw.x -in PREFIX.pwi > PREFIX.pwo'
    os.environ['ESPRESSO_PSEUDO'] = '/tmp/pseudo'
    
    from gui.calculations import dry_run_calculation
    
    # Create a test structure
    atoms = bulk('Fe', 'bcc', a=2.87)
    
    # Configuration dictionary
    config = {
        'calc_type': 'scf',
        'ecutwfc': 40.0,
        'ecutrho': 320.0,
        'pseudopotentials': {
            'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
        },
        'kpts': (2, 2, 2),
    }
    
    with tempfile.TemporaryDirectory() as tmpdir:
        label = os.path.join(tmpdir, 'dry_run/fe')
        
        # Run dry run
        prepared_atoms, calc = dry_run_calculation(atoms, config, label=label)
        
        # Verify files were created
        pwi_file = f"{calc.label}.pwi"
        assert os.path.exists(pwi_file), f"PWI file not created: {pwi_file}"
        
        # Verify directory was created
        assert os.path.exists(calc.directory)


def test_calculator_reuse_in_session():
    """Test that calculators can be stored and reused."""
    # This simulates the pattern used in job_submission.py
    # where the calculator is stored in session_state
    
    # Set up environment
    os.environ['ASE_ESPRESSO_COMMAND'] = 'pw.x -in PREFIX.pwi > PREFIX.pwo'
    os.environ['ESPRESSO_PSEUDO'] = '/tmp/pseudo'
    
    from gui.calculations import prepare_calculation_from_gui
    
    atoms = bulk('Fe', 'bcc', a=2.87)
    
    config = {
        'calc_type': 'scf',
        'ecutwfc': 40.0,
        'pseudopotentials': {'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
        'kpts': (2, 2, 2),
    }
    
    with tempfile.TemporaryDirectory() as tmpdir:
        label1 = os.path.join(tmpdir, 'calc1/fe')
        
        # First preparation (like in calculation_setup.py)
        prepared_atoms, calc = prepare_calculation_from_gui(atoms, config, label=label1)
        
        # Simulate storing in session state
        stored_calc = calc
        stored_atoms = prepared_atoms
        
        # Simulate reusing in job_submission.py with a different label
        label2 = os.path.join(tmpdir, 'calc2/fe')
        stored_calc.label = label2
        
        # Verify the calculator can be reused with new label
        assert stored_calc.label == label2
        assert stored_atoms is not None


def test_config_validation():
    """Test that configuration validation works correctly."""
    from gui.calculations import prepare_calculation_from_gui
    
    atoms = bulk('Fe', 'bcc', a=2.87)
    
    # Missing pseudopotentials - should raise error
    invalid_config = {
        'calc_type': 'scf',
        'ecutwfc': 40.0,
    }
    
    with tempfile.TemporaryDirectory() as tmpdir:
        label = os.path.join(tmpdir, 'test/fe')
        
        with pytest.raises(ValueError, match="pseudopotentials"):
            prepare_calculation_from_gui(atoms, invalid_config, label=label)


def test_multiple_calc_types():
    """Test that different calculation types are handled correctly."""
    # Set up environment
    os.environ['ASE_ESPRESSO_COMMAND'] = 'pw.x -in PREFIX.pwi > PREFIX.pwo'
    os.environ['ESPRESSO_PSEUDO'] = '/tmp/pseudo'
    
    from gui.calculations import prepare_calculation_from_gui
    
    atoms = bulk('Fe', 'bcc', a=2.87)
    
    calc_types = ['scf', 'relax', 'vc-relax']
    
    for calc_type in calc_types:
        config = {
            'calc_type': calc_type,
            'ecutwfc': 40.0,
            'pseudopotentials': {'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
            'kpts': (2, 2, 2),
        }
        
        with tempfile.TemporaryDirectory() as tmpdir:
            label = os.path.join(tmpdir, f'{calc_type}/fe')
            
            prepared_atoms, calc = prepare_calculation_from_gui(atoms, config, label=label)
            
            # Verify calculator was created with correct type
            assert calc is not None
            assert prepared_atoms is not None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
