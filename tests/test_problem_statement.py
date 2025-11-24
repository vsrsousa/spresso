"""
Final test reproducing the exact problem statement scenario
This test verifies that the bug where Hubbard parameters were 
silently removed from QE input files is fixed.
"""

import pytest
import tempfile
import os
import numpy as np
from ase.build import bulk
from xespresso.xio import write_espresso_in


def test_problem_statement_scenario():
    """
    Reproduces the exact scenario from the problem statement.
    
    The original issue was:
    - User creates an Espresso calculator with Hubbard parameters in QE 7.2 format
    - When write_input is called, the HUBBARD card is NOT generated
    - The input file is missing the Hubbard parameters completely
    
    This test verifies that the issue is fixed.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        atoms = bulk("Fe", cubic=True)
        atoms.new_array("species", np.array(atoms.get_chemical_symbols(), dtype="U20"))

        input_data_complete = {
            "ecutwfc": 30.0,
            "ecutrho": 240.0,
            "occupations": "smearing",
            "smearing": "gaussian",
            "degauss": 0.02,
            "nspin": 2,
            "lda_plus_u": True,
            "qe_version": "7.2",
            "hubbard": {
                "projector": "atomic",
                "u": {
                    "Fe-3d": 4.3,
                }
            }
        }

        pseudopotentials = {
            "Fe": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF",
        }

        input_file = os.path.join(tmpdir, "fe.pwi")
        write_espresso_in(
            input_file,
            atoms,
            input_data=input_data_complete,
            pseudopotentials=pseudopotentials,
            kpts=(4, 4, 4)
        )

        # Read the generated input file
        with open(input_file, 'r') as f:
            content = f.read()
        
        # Verify the fix: The HUBBARD card should be present
        assert "HUBBARD" in content, (
            "HUBBARD card is missing! This is the bug from the problem statement."
        )
        
        # Verify the specific parameters are correct
        assert "{atomic}" in content, "Projector type should be 'atomic'"
        assert "Fe-3d" in content, "Fe-3d orbital should be specified"
        assert "4.3" in content, "Hubbard U value should be 4.3"
        
        # Verify the HUBBARD card structure is correct
        lines = content.split('\n')
        hubbard_index = None
        for i, line in enumerate(lines):
            if 'HUBBARD' in line:
                hubbard_index = i
                break
        
        assert hubbard_index is not None, "HUBBARD card not found"
        
        # Check the lines around HUBBARD card
        hubbard_section = '\n'.join(lines[hubbard_index:hubbard_index+3])
        assert 'HUBBARD {atomic}' in hubbard_section
        assert 'U Fe-3d 4.3' in hubbard_section


def test_problem_statement_expected_output():
    """
    Verify the output matches the expected QE 7.x format
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        atoms = bulk("Fe", cubic=True)
        atoms.new_array("species", np.array(atoms.get_chemical_symbols(), dtype="U20"))

        input_data_complete = {
            "ecutwfc": 30.0,
            "ecutrho": 240.0,
            "occupations": "smearing",
            "smearing": "gaussian",
            "degauss": 0.02,
            "nspin": 2,
            "lda_plus_u": True,
            "qe_version": "7.2",
            "hubbard": {
                "projector": "atomic",
                "u": {
                    "Fe-3d": 4.3,
                }
            }
        }

        pseudopotentials = {
            "Fe": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF",
        }

        input_file = os.path.join(tmpdir, "fe.pwi")
        write_espresso_in(
            input_file,
            atoms,
            input_data=input_data_complete,
            pseudopotentials=pseudopotentials,
            kpts=(4, 4, 4)
        )

        with open(input_file, 'r') as f:
            content = f.read()
        
        # The input file should have standard sections
        assert "&CONTROL" in content
        assert "&SYSTEM" in content
        assert "&ELECTRONS" in content
        assert "ATOMIC_SPECIES" in content
        assert "K_POINTS automatic" in content
        assert "CELL_PARAMETERS angstrom" in content
        assert "ATOMIC_POSITIONS angstrom" in content
        
        # AND the HUBBARD card
        assert "HUBBARD {atomic}" in content
        assert "U Fe-3d 4.3" in content
        
        # Note: In QE 7.x with new HUBBARD card format, lda_plus_u is implicit
        # and doesn't need to be in the SYSTEM namelist


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
