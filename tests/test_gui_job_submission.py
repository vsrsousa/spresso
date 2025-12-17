"""
Tests for GUI job submission functionality.

These tests verify that the GUI can properly generate input files
for dry run job submission using xespresso's write_input function.
"""

import pytest
import os
import tempfile
from pathlib import Path
from ase.build import bulk


def test_espresso_write_input_creates_files():
    """Test that Espresso.write_input creates the expected files."""
    # Set up environment
    os.environ['ASE_ESPRESSO_COMMAND'] = 'pw.x -in PREFIX.pwi > PREFIX.pwo'
    os.environ['ESPRESSO_PSEUDO'] = '/tmp/pseudo'
    
    from xespresso import Espresso
    
    # Create a test structure
    atoms = bulk('Fe', 'bcc', a=2.87)
    
    # Create temp directory
    with tempfile.TemporaryDirectory() as tmpdir:
        label = os.path.join(tmpdir, 'test/fe')
        
        # Setup calculator
        input_data = {
            'calculation': 'scf',
            'ecutwfc': 40.0,
            'ecutrho': 320.0,
            'conv_thr': 1.0e-6,
            'occupations': 'smearing',
            'smearing': 'gaussian',
            'degauss': 0.02,
        }
        
        pseudopotentials = {
            'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
        }
        
        calc = Espresso(
            label=label,
            pseudopotentials=pseudopotentials,
            input_data=input_data,
            kpts=(4, 4, 4),
            queue={},
        )
        
        atoms.calc = calc
        
        # Test write_input (dry run)
        calc.write_input(atoms)
        
        # Check that files were created
        pwi_file = f"{calc.label}.pwi"
        asei_file = f"{calc.label}.asei"
        
        assert os.path.exists(pwi_file), f"PWI file not created: {pwi_file}"
        assert os.path.exists(asei_file), f"ASEI file not created: {asei_file}"
        
        # Check PWI file content
        with open(pwi_file, 'r') as f:
            content = f.read()
            assert '&CONTROL' in content
            assert '&SYSTEM' in content
            assert "calculation      = 'scf'" in content
            assert 'ecutwfc          = 40.0' in content


def test_espresso_write_input_with_queue():
    """Test that Espresso.write_input creates job script with queue config."""
    # Set up environment
    os.environ['ASE_ESPRESSO_COMMAND'] = 'pw.x -in PREFIX.pwi > PREFIX.pwo'
    os.environ['ESPRESSO_PSEUDO'] = '/tmp/pseudo'
    
    from xespresso import Espresso
    
    # Create a test structure
    atoms = bulk('Fe', 'bcc', a=2.87)
    
    # Create temp directory
    with tempfile.TemporaryDirectory() as tmpdir:
        label = os.path.join(tmpdir, 'test/fe')
        
        # Setup calculator with queue
        input_data = {
            'calculation': 'scf',
            'ecutwfc': 40.0,
            'ecutrho': 320.0,
            'conv_thr': 1.0e-6,
            'occupations': 'smearing',
            'smearing': 'gaussian',
            'degauss': 0.02,
        }
        
        pseudopotentials = {
            'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
        }
        
        queue = {
            'nodes': 1,
            'ntasks-per-node': 4,
            'partition': 'normal',
            'time': '1:00:00',
        }
        
        calc = Espresso(
            label=label,
            pseudopotentials=pseudopotentials,
            input_data=input_data,
            kpts=(4, 4, 4),
            queue=queue,
        )
        
        atoms.calc = calc
        
        # Test write_input (dry run)
        calc.write_input(atoms)
        
        # Check for job script
        job_file = os.path.join(calc.directory, 'job_file')
        assert os.path.exists(job_file), f"Job script not created: {job_file}"
        
        # Check job script content
        with open(job_file, 'r') as f:
            content = f.read()
            assert '#!/bin/bash' in content
            assert 'pw.x' in content


def test_espresso_label_creates_directory():
    """Test that the label parameter correctly creates directories."""
    # Set up environment
    os.environ['ASE_ESPRESSO_COMMAND'] = 'pw.x -in PREFIX.pwi > PREFIX.pwo'
    os.environ['ESPRESSO_PSEUDO'] = '/tmp/pseudo'
    
    from xespresso import Espresso
    
    # Create a test structure
    atoms = bulk('Fe', 'bcc', a=2.87)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Test with nested label
        label = os.path.join(tmpdir, 'calc/scf/fe')
        
        calc = Espresso(
            label=label,
            pseudopotentials={'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
            input_data={'calculation': 'scf', 'ecutwfc': 40.0},
            kpts=(2, 2, 2),
            queue={},
        )
        
        atoms.calc = calc
        calc.write_input(atoms)
        
        # Check directory structure
        assert os.path.exists(calc.directory), f"Directory not created: {calc.directory}"
        assert os.path.basename(calc.directory) == 'fe'
        
        # Check files are in the correct location
        pwi_file = f"{calc.label}.pwi"
        assert os.path.exists(pwi_file)
        assert os.path.dirname(pwi_file) == calc.directory


def test_espresso_with_kspacing():
    """Test that Espresso works with kspacing parameter."""
    # Set up environment
    os.environ['ASE_ESPRESSO_COMMAND'] = 'pw.x -in PREFIX.pwi > PREFIX.pwo'
    os.environ['ESPRESSO_PSEUDO'] = '/tmp/pseudo'
    
    from xespresso import Espresso
    
    # Create a test structure
    atoms = bulk('Fe', 'bcc', a=2.87)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        label = os.path.join(tmpdir, 'test/fe')
        
        # Test with kspacing instead of kpts
        calc = Espresso(
            label=label,
            pseudopotentials={'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
            input_data={'calculation': 'scf', 'ecutwfc': 40.0},
            kspacing=0.3,  # Using kspacing
            queue={},
        )
        
        atoms.calc = calc
        calc.write_input(atoms)
        
        # Check that files were created
        pwi_file = f"{calc.label}.pwi"
        assert os.path.exists(pwi_file)
        
        # Check that k-points were generated from kspacing
        with open(pwi_file, 'r') as f:
            content = f.read()
            assert 'K_POINTS' in content


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
