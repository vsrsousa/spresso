"""
Test for read_convergence method to handle various pwo file formats
including files with warnings and extra headers.
"""

import os
import tempfile
import pytest
from _common_helpers import set_envs
from xespresso import Espresso


class TestReadConvergence:
    """Test read_convergence with various pwo file formats"""

    def test_read_convergence_with_mpi_warning(self):
        """
        Test that read_convergence handles pwo files with MPI warnings
        This reproduces the exact error from the problem statement.
        """
        set_envs()
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a mock pwo file with MPI warning (like in problem statement)
            pwo_content = """MPI startup(): PMI server not found. Please set I_MPI_PMI_LIBRARY variable if it is not a singleton case.

     Program PWSCF v.7.4.1 starts on 17Oct2025 at  1:39:50 

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
         "P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009);
         "P. Giannozzi et al., J. Phys.:Condens. Matter 29 465901 (2017);
         "P. Giannozzi et al., J. Chem. Phys. 152 154105 (2020);
          URL http://www.quantum-espresso.org", 
     in publications or presentations arising from this work. More details at
     http://www.quantum-espresso.org/quote

     Parallel version (MPI), running on     1 processors

     MPI processes distributed on     1 nodes
     JOB DONE.
"""
            label_dir = os.path.join(tmpdir, "test")
            os.makedirs(label_dir, exist_ok=True)
            pwo_file = os.path.join(label_dir, "test.pwo")
            with open(pwo_file, "w") as f:
                f.write(pwo_content)

            # Create a minimal Espresso calculator instance
            calc = Espresso(label=label_dir, debug=True, queue={})
            
            # This should not raise an IndexError
            convergence, msg = calc.read_convergence()
            
            # Should successfully detect the job completion
            assert convergence == 0, f"Expected convergence=0, got {convergence}"
            assert "JOB DONE" in msg

    def test_read_convergence_normal_output(self):
        """Test with normal pwo file (no warnings)"""
        set_envs()
        with tempfile.TemporaryDirectory() as tmpdir:
            pwo_content = """
     Program PWSCF v.7.4.1 starts on 17Oct2025 at  1:39:50 

     This program is part of the open-source Quantum ESPRESSO suite
     JOB DONE.
"""
            label_dir = os.path.join(tmpdir, "test")
            os.makedirs(label_dir, exist_ok=True)
            pwo_file = os.path.join(label_dir, "test.pwo")
            with open(pwo_file, "w") as f:
                f.write(pwo_content)

            calc = Espresso(label=label_dir, debug=True, queue={})
            convergence, msg = calc.read_convergence()
            
            assert convergence == 0
            assert "JOB DONE" in msg

    def test_read_convergence_empty_file(self):
        """Test with empty pwo file"""
        set_envs()
        with tempfile.TemporaryDirectory() as tmpdir:
            label_dir = os.path.join(tmpdir, "test")
            os.makedirs(label_dir, exist_ok=True)
            pwo_file = os.path.join(label_dir, "test.pwo")
            with open(pwo_file, "w") as f:
                f.write("")

            calc = Espresso(label=label_dir, debug=True, queue={})
            convergence, msg = calc.read_convergence()
            
            # Should return error code for empty file
            assert convergence == 1
            assert "pwo file has nothing" in msg

    def test_read_convergence_no_file(self):
        """Test when pwo file doesn't exist"""
        set_envs()
        with tempfile.TemporaryDirectory() as tmpdir:
            calc = Espresso(label=os.path.join(tmpdir, "nonexistent"), debug=True, queue={})
            convergence, msg = calc.read_convergence()
            
            # Should return error code for missing file
            assert convergence == 3
            assert "No pwo output file" in msg

    def test_read_convergence_not_converged(self):
        """Test with non-converged calculation"""
        set_envs()
        with tempfile.TemporaryDirectory() as tmpdir:
            pwo_content = """
     Program PWSCF v.7.4.1 starts on 17Oct2025 at  1:39:50 

     convergence NOT achieved after 100 iterations: stopping
"""
            label_dir = os.path.join(tmpdir, "test")
            os.makedirs(label_dir, exist_ok=True)
            pwo_file = os.path.join(label_dir, "test.pwo")
            with open(pwo_file, "w") as f:
                f.write(pwo_content)

            calc = Espresso(label=label_dir, debug=True, queue={})
            convergence, msg = calc.read_convergence()
            
            # Should detect non-convergence
            assert convergence == 1
            assert "convergence NOT achieved" in msg

    def test_read_convergence_max_cpu_time(self):
        """Test with maximum CPU time exceeded"""
        set_envs()
        with tempfile.TemporaryDirectory() as tmpdir:
            pwo_content = """
     Program PWSCF v.7.4.1 starts on 17Oct2025 at  1:39:50 

     Maximum CPU time exceeded
"""
            label_dir = os.path.join(tmpdir, "test")
            os.makedirs(label_dir, exist_ok=True)
            pwo_file = os.path.join(label_dir, "test.pwo")
            with open(pwo_file, "w") as f:
                f.write(pwo_content)

            calc = Espresso(label=label_dir, debug=True, queue={})
            convergence, msg = calc.read_convergence()
            
            # Should detect CPU time exceeded
            assert convergence == 2
            assert "Maximum CPU time exceeded" in msg

    def test_read_convergence_unknown_status(self):
        """Test with file that doesn't match any known pattern"""
        set_envs()
        with tempfile.TemporaryDirectory() as tmpdir:
            pwo_content = """
     Program PWSCF v.7.4.1 starts on 17Oct2025 at  1:39:50 

     Some calculation output
     More output
     Still running maybe?
"""
            label_dir = os.path.join(tmpdir, "test")
            os.makedirs(label_dir, exist_ok=True)
            pwo_file = os.path.join(label_dir, "test.pwo")
            with open(pwo_file, "w") as f:
                f.write(pwo_content)

            calc = Espresso(label=label_dir, debug=True, queue={})
            convergence, msg = calc.read_convergence()
            
            # Should return unknown error code
            assert convergence == 4


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
