"""
Integration test for Hubbard parameters with Espresso class
This test reproduces the issue from the problem statement.
"""

import pytest
import tempfile
import os
import numpy as np
from ase.build import bulk
from xespresso.xio import write_espresso_in


class TestHubbardIntegration:
    """Integration tests for Hubbard parameters"""
    
    def test_hubbard_parameters_in_write_input(self):
        """Test that Hubbard parameters are correctly written to input file"""
        # Create temporary directory for test
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

            # Read the input file
            with open(input_file, 'r') as f:
                content = f.read()
            
            # Check if HUBBARD card is present
            assert "HUBBARD" in content, "HUBBARD card should be present in the input file"
            assert "Fe-3d" in content, "Fe-3d should be specified in HUBBARD card"
            assert "4.3" in content, "Hubbard U value should be present"
            assert "{atomic}" in content, "Projector type should be specified"
    
    def test_hubbard_parameters_with_multiple_species(self):
        """Test Hubbard parameters with multiple species"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a simple FeO structure (just for testing, not realistic)
            atoms = bulk("Fe", cubic=True)
            atoms.append('O')
            atoms.new_array("species", np.array(['Fe', 'Fe', 'O'], dtype="U20"))

            input_data = {
                "ecutwfc": 30.0,
                "ecutrho": 240.0,
                "qe_version": "7.2",
                "hubbard": {
                    "projector": "ortho-atomic",
                    "u": {
                        "Fe-3d": 4.3,
                        "O-2p": 5.0,
                    }
                }
            }

            pseudopotentials = {
                "Fe": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF",
                "O": "O.pbe-n-kjpaw_psl.0.1.UPF",
            }

            input_file = os.path.join(tmpdir, "feo.pwi")
            write_espresso_in(
                input_file,
                atoms,
                input_data=input_data,
                pseudopotentials=pseudopotentials,
                kpts=(4, 4, 4)
            )

            # Read the input file
            with open(input_file, 'r') as f:
                content = f.read()
            
            # Check if both Hubbard U parameters are present
            assert "HUBBARD" in content
            assert "Fe-3d" in content
            assert "O-2p" in content
            assert "4.3" in content
            assert "5.0" in content
            assert "{ortho-atomic}" in content
    
    def test_old_format_hubbard_still_works(self):
        """Test that old format Hubbard parameters still work"""
        with tempfile.TemporaryDirectory() as tmpdir:
            atoms = bulk("Fe", cubic=True)
            atoms.new_array("species", np.array(atoms.get_chemical_symbols(), dtype="U20"))

            # Old format using input_ntyp
            input_data = {
                "ecutwfc": 30.0,
                "lda_plus_u": True,
                "qe_version": "6.8",  # Old QE version
                "input_ntyp": {
                    "Hubbard_U": {
                        "Fe": 4.3,
                    }
                }
            }

            pseudopotentials = {
                "Fe": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF",
            }

            input_file = os.path.join(tmpdir, "fe_old.pwi")
            write_espresso_in(
                input_file,
                atoms,
                input_data=input_data,
                pseudopotentials=pseudopotentials,
                kpts=(4, 4, 4)
            )

            # Read the input file
            with open(input_file, 'r') as f:
                content = f.read()
            
            # Old format should NOT have HUBBARD card
            # Instead, it should be in SYSTEM namelist
            assert "HUBBARD" not in content, "Old format should not have HUBBARD card"
            # QE input is case-insensitive, check for lowercase version
            assert "hubbard_u(1)" in content.lower(), "Old format should have Hubbard_U in SYSTEM"
            assert "4.3" in content


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
