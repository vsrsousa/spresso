"""
Tests for consecutive calls to get_potential_energy() with magnetic setup.

This test verifies that the calculator correctly detects when magnetic
parameters (INPUT_NTYP) have changed between consecutive calculations,
ensuring that the job is properly flagged to be rerun when needed.
"""

import pytest
import tempfile
import os
import shutil
import numpy as np
from ase.build import bulk
from xespresso import Espresso
from xespresso.tools import set_antiferromagnetic, set_magnetic_moments
from xespresso.xio import write_espresso_asei


class TestConsecutiveMagneticCalls:
    """Test consecutive get_potential_energy() calls with magnetic setup."""
    
    @pytest.fixture(autouse=True)
    def setup_env(self):
        """Set up environment for tests."""
        # Set required environment variables
        os.environ["ESPRESSO_PSEUDO"] = os.path.join(os.getcwd(), "tests/datas/pseudo")
        os.environ["ASE_ESPRESSO_COMMAND"] = "pw.x -in PREFIX.pwi > PREFIX.pwo"
        
        # Create temporary directory
        self.tmpdir = tempfile.mkdtemp()
        yield
        # Cleanup
        if os.path.exists(self.tmpdir):
            shutil.rmtree(self.tmpdir)
    
    def test_same_magnetic_parameters_no_rerun(self):
        """
        Test that check_state returns False when magnetic parameters are unchanged.
        
        When the same magnetic configuration is used in consecutive calculations,
        the calculator should detect that no parameters have changed and avoid
        unnecessarily rerunning the calculation.
        """
        atoms = bulk('Fe', cubic=True)
        
        # Set up antiferromagnetic configuration
        mag_config = set_antiferromagnetic(atoms, [[0], [1]], magnetic_moment=1.0)
        pseudopotentials = {
            'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
            'Fe1': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
        }
        mag_config['pseudopotentials'] = pseudopotentials
        
        # First calculation
        label = os.path.join(self.tmpdir, 'fe')
        
        calc1 = Espresso(
            label=label,
            pseudopotentials=mag_config['pseudopotentials'],
            input_data={'input_ntyp': mag_config['input_ntyp']},
            nspin=2,
            ecutwfc=30,
            kpts=(2, 2, 2),
            queue={},
        )
        
        # Create directory and write asei file to simulate a completed calculation
        os.makedirs(calc1.directory, exist_ok=True)
        asei_file = calc1.asei
        write_espresso_asei(asei_file, atoms, calc1.parameters)
        
        # Create a fake output file to simulate a converged calculation
        pwo_file = os.path.join(calc1.directory, f"{calc1.prefix}.pwo")
        with open(pwo_file, 'w') as f:
            f.write("     JOB DONE.\n")
        
        # Second calculation with same parameters
        calc2 = Espresso(
            label=label,
            pseudopotentials=mag_config['pseudopotentials'],
            input_data={'input_ntyp': mag_config['input_ntyp']},
            nspin=2,
            ecutwfc=30,
            kpts=(2, 2, 2),
            queue={},
        )
        
        # check_state should return False (no need to rerun)
        needs_run = calc2.check_state(atoms)
        
        assert not needs_run, (
            "check_state should return False when magnetic parameters are unchanged. "
            f"Changed parameters: {calc2.changed_parameters}"
        )
    
    def test_different_magnetic_moment_requires_rerun(self):
        """
        Test that check_state returns True when magnetic moments change.
        
        When magnetic moments are changed between calculations, the calculator
        must detect this change and flag that the job needs to be rerun.
        """
        atoms = bulk('Fe', cubic=True)
        
        # First calculation with magnetic moment = 1.0
        mag_config1 = set_antiferromagnetic(atoms, [[0], [1]], magnetic_moment=1.0)
        pseudopotentials = {
            'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
            'Fe1': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
        }
        mag_config1['pseudopotentials'] = pseudopotentials
        
        label = os.path.join(self.tmpdir, 'fe2')
        
        calc1 = Espresso(
            label=label,
            pseudopotentials=mag_config1['pseudopotentials'],
            input_data={'input_ntyp': mag_config1['input_ntyp']},
            nspin=2,
            ecutwfc=30,
            kpts=(2, 2, 2),
            queue={},
        )
        
        # Create directory and write asei file to simulate a completed calculation
        os.makedirs(calc1.directory, exist_ok=True)
        asei_file = calc1.asei
        write_espresso_asei(asei_file, atoms, calc1.parameters)
        
        # Second calculation with different magnetic moment = 2.0
        mag_config2 = set_antiferromagnetic(atoms, [[0], [1]], magnetic_moment=2.0)
        mag_config2['pseudopotentials'] = pseudopotentials
        
        calc2 = Espresso(
            label=label,
            pseudopotentials=mag_config2['pseudopotentials'],
            input_data={'input_ntyp': mag_config2['input_ntyp']},
            nspin=2,
            ecutwfc=30,
            kpts=(2, 2, 2),
            queue={},
        )
        
        # check_state should return True (needs rerun due to changed magnetic moment)
        needs_run = calc2.check_state(atoms)
        
        assert needs_run, (
            "check_state should return True when magnetic moments change. "
            f"Changed parameters: {calc2.changed_parameters}"
        )
        assert 'starting_magnetization' in calc2.changed_parameters, (
            "starting_magnetization should be in changed_parameters"
        )
    
    def test_adding_magnetic_setup_requires_rerun(self):
        """
        Test that adding magnetic setup to a non-magnetic calculation requires rerun.
        
        When magnetic parameters (INPUT_NTYP) are added to a calculation that
        previously didn't have them, the calculator must detect this and flag
        that the job needs to be rerun.
        """
        atoms = bulk('Fe', cubic=True)
        
        # First calculation without magnetic setup
        label = os.path.join(self.tmpdir, 'fe3')
        
        pseudopotentials = {'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'}
        
        calc1 = Espresso(
            label=label,
            pseudopotentials=pseudopotentials,
            ecutwfc=30,
            kpts=(2, 2, 2),
            queue={},
        )
        
        # Write asei file to simulate a completed non-magnetic calculation
        # Reset atoms species array to non-magnetic
        if 'species' in atoms.arrays:
            del atoms.arrays['species']
        
        os.makedirs(calc1.directory, exist_ok=True)
        asei_file = calc1.asei
        write_espresso_asei(asei_file, atoms, calc1.parameters)
        
        # Second calculation with magnetic setup added
        atoms2 = bulk('Fe', cubic=True)
        mag_config = set_antiferromagnetic(atoms2, [[0], [1]], magnetic_moment=1.0)
        pseudopotentials2 = {
            'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
            'Fe1': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
        }
        mag_config['pseudopotentials'] = pseudopotentials2
        
        calc2 = Espresso(
            label=label,
            pseudopotentials=mag_config['pseudopotentials'],
            input_data={'input_ntyp': mag_config['input_ntyp']},
            nspin=2,
            ecutwfc=30,
            kpts=(2, 2, 2),
            queue={},
        )
        
        # check_state should return True (needs rerun due to added magnetic setup)
        needs_run = calc2.check_state(atoms2)
        
        assert needs_run, (
            "check_state should return True when magnetic setup is added. "
            f"Changed parameters: {calc2.changed_parameters}"
        )
        # When INPUT_NTYP is added, either 'INPUT_NTYP' or 'starting_magnetization' 
        # should be in changed_parameters (both indicate magnetic setup changed)
        assert ('INPUT_NTYP' in calc2.changed_parameters or 
                'starting_magnetization' in calc2.changed_parameters), (
            "Either INPUT_NTYP or starting_magnetization should be in changed_parameters "
            "when magnetic setup is added"
        )
    
    def test_removing_magnetic_setup_requires_rerun(self):
        """
        Test that removing magnetic setup from a calculation requires rerun.
        
        When magnetic parameters (INPUT_NTYP) are removed from a calculation that
        previously had them, the calculator must detect this and flag that the
        job needs to be rerun.
        """
        atoms = bulk('Fe', cubic=True)
        
        # First calculation with magnetic setup
        mag_config = set_antiferromagnetic(atoms, [[0], [1]], magnetic_moment=1.0)
        pseudopotentials = {
            'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
            'Fe1': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
        }
        mag_config['pseudopotentials'] = pseudopotentials
        
        label = os.path.join(self.tmpdir, 'fe4')
        
        calc1 = Espresso(
            label=label,
            pseudopotentials=mag_config['pseudopotentials'],
            input_data={'input_ntyp': mag_config['input_ntyp']},
            nspin=2,
            ecutwfc=30,
            kpts=(2, 2, 2),
            queue={},
        )
        
        # Write asei file to simulate a completed magnetic calculation
        os.makedirs(calc1.directory, exist_ok=True)
        asei_file = calc1.asei
        write_espresso_asei(asei_file, atoms, calc1.parameters)
        
        # Second calculation without magnetic setup
        atoms2 = bulk('Fe', cubic=True)
        if 'species' in atoms2.arrays:
            del atoms2.arrays['species']
        
        calc2 = Espresso(
            label=label,
            pseudopotentials={'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
            ecutwfc=30,
            kpts=(2, 2, 2),
            queue={},
        )
        
        # check_state should return True (needs rerun due to removed magnetic setup)
        needs_run = calc2.check_state(atoms2)
        
        assert needs_run, (
            "check_state should return True when magnetic setup is removed. "
            f"Changed parameters: {calc2.changed_parameters}"
        )
    
    def test_complex_magnetic_configuration_changes(self):
        """
        Test detection of changes in complex magnetic configurations.
        
        Test with multiple elements and different magnetic moments to ensure
        the calculator properly detects any changes in the magnetic setup.
        """
        # Create a compound with Fe and Ni
        from ase import Atoms
        atoms = Atoms('Fe2Ni2', 
                     positions=[[0, 0, 0], [1.5, 0, 0], [0, 1.5, 0], [1.5, 1.5, 0]],
                     cell=[5, 5, 5])
        
        # First calculation: Fe AFM, Ni FM
        mag_config1 = set_magnetic_moments(atoms, [1.0, -1.0, 0.5, 0.5])
        pseudopotentials = {
            'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
            'Fe1': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
            'Ni': 'Ni.pbe-n-rrkjus_psl.1.0.0.UPF'
        }
        mag_config1['pseudopotentials'] = pseudopotentials
        
        label = os.path.join(self.tmpdir, 'feni')
        
        calc1 = Espresso(
            label=label,
            pseudopotentials=mag_config1['pseudopotentials'],
            input_data={'input_ntyp': mag_config1['input_ntyp']},
            nspin=2,
            ecutwfc=30,
            kpts=(2, 2, 2),
            queue={},
        )
        
        # Write asei file
        os.makedirs(calc1.directory, exist_ok=True)
        asei_file = calc1.asei
        write_espresso_asei(asei_file, atoms, calc1.parameters)
        
        # Second calculation: Change Fe magnetic moments
        atoms2 = Atoms('Fe2Ni2', 
                      positions=[[0, 0, 0], [1.5, 0, 0], [0, 1.5, 0], [1.5, 1.5, 0]],
                      cell=[5, 5, 5])
        mag_config2 = set_magnetic_moments(atoms2, [2.0, -2.0, 0.5, 0.5])
        mag_config2['pseudopotentials'] = pseudopotentials
        
        calc2 = Espresso(
            label=label,
            pseudopotentials=mag_config2['pseudopotentials'],
            input_data={'input_ntyp': mag_config2['input_ntyp']},
            nspin=2,
            ecutwfc=30,
            kpts=(2, 2, 2),
            queue={},
        )
        
        # check_state should return True (Fe magnetic moments changed)
        needs_run = calc2.check_state(atoms2)
        
        assert needs_run, (
            "check_state should return True when any magnetic moment changes in complex system. "
            f"Changed parameters: {calc2.changed_parameters}"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
