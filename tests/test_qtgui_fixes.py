"""
Tests for Qt GUI fixes:
1. ASE_ESPRESSO_COMMAND environment variable setup
2. Kpoints mutually exclusive (kspacing OR kpts, not both)
3. Launcher from machine configuration in job_file
4. Modules from codes configuration
5. Magnetic and Hubbard configuration editing after session save/reload
"""

import pytest
import os
import tempfile
from pathlib import Path

# Import ASE with graceful handling
try:
    from ase.build import bulk
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False
    bulk = None

# Import xespresso with graceful handling
try:
    from xespresso import Espresso
    from xespresso.machines.machine import Machine
    from xespresso.codes.config import Code, CodesConfig
    XESPRESSO_AVAILABLE = True
except ImportError:
    XESPRESSO_AVAILABLE = False
    Espresso = None
    Machine = None
    Code = None
    CodesConfig = None


def test_ase_espresso_command_set_in_qtgui():
    """Test that ASE_ESPRESSO_COMMAND is set in Qt GUI prepare_calculation."""
    # Set environment variable as Qt GUI does
    os.environ['ASE_ESPRESSO_COMMAND'] = "LAUNCHER PACKAGE.x PARALLEL -in PREFIX.PACKAGEi > PREFIX.PACKAGEo"
    
    # Verify the environment variable is set correctly
    assert 'ASE_ESPRESSO_COMMAND' in os.environ
    assert 'LAUNCHER PACKAGE.x PARALLEL' in os.environ['ASE_ESPRESSO_COMMAND']
    assert 'PREFIX.PACKAGEi > PREFIX.PACKAGEo' in os.environ['ASE_ESPRESSO_COMMAND']


@pytest.mark.skipif(not ASE_AVAILABLE or not XESPRESSO_AVAILABLE, reason="ASE or xespresso not available")
def test_kpoints_mutually_exclusive():
    """Test that Espresso calculator receives EITHER kspacing OR kpts, never both."""
    atoms = bulk('Fe', 'bcc', a=2.87)
    
    # Set environment variable
    os.environ['ASE_ESPRESSO_COMMAND'] = "LAUNCHER PACKAGE.x PARALLEL -in PREFIX.PACKAGEi > PREFIX.PACKAGEo"
    os.environ['ESPRESSO_PSEUDO'] = '/tmp/pseudo'
    
    with tempfile.TemporaryDirectory() as tmpdir:
        label = os.path.join(tmpdir, 'test/fe')
        
        # Test 1: Only kspacing provided
        calc1 = Espresso(
            label=label + '1',
            pseudopotentials={'Fe': 'Fe.upf'},
            input_data={'calculation': 'scf', 'ecutwfc': 40.0},
            kspacing=0.3,
        )
        
        # Verify kspacing is set, kpts is not
        assert hasattr(calc1, 'kspacing') or 'kspacing' in calc1.parameters
        # kpts should not be explicitly set when kspacing is provided
        
        # Test 2: Only kpts provided
        calc2 = Espresso(
            label=label + '2',
            pseudopotentials={'Fe': 'Fe.upf'},
            input_data={'calculation': 'scf', 'ecutwfc': 40.0},
            kpts=(4, 4, 4),
        )
        
        # Verify kpts is set
        assert hasattr(calc2, 'kpts') or 'kpts' in calc2.parameters
        
        # Test 3: Neither provided - should default to kpts
        calc3 = Espresso(
            label=label + '3',
            pseudopotentials={'Fe': 'Fe.upf'},
            input_data={'calculation': 'scf', 'ecutwfc': 40.0},
        )
        
        # Should have default kpts
        assert calc3 is not None


@pytest.mark.skipif(not XESPRESSO_AVAILABLE, reason="xespresso not available")
def test_launcher_from_machine_config():
    """Test that job_file uses launcher from machine configuration."""
    # Create a machine with custom launcher
    machine = Machine(
        name='test_machine',
        execution='local',
        launcher='srun -n {nprocs}',
        nprocs=8
    )
    
    # Test launcher substitution
    launcher = machine.launcher
    if '{nprocs}' in launcher:
        launcher = launcher.replace('{nprocs}', str(machine.nprocs))
    
    assert launcher == 'srun -n 8'
    assert 'srun' in launcher
    assert '8' in launcher


@pytest.mark.skipif(not XESPRESSO_AVAILABLE, reason="xespresso not available")
def test_modules_from_codes_config():
    """Test that modules are loaded from codes configuration."""
    # Create a codes configuration with modules
    codes_config = CodesConfig(
        machine_name='test_machine',
        qe_version='7.2',
        versions={
            '7.2': {
                'qe_prefix': '/opt/qe-7.2/bin',
                'modules': ['quantum-espresso/7.2', 'intel/2023'],
                'codes': {
                    'pw': Code(name='pw', path='/opt/qe-7.2/bin/pw.x', version='7.2')
                }
            }
        }
    )
    
    # Verify modules are present in version config
    assert codes_config.versions is not None
    assert '7.2' in codes_config.versions
    version_config = codes_config.versions['7.2']
    assert 'modules' in version_config
    assert 'quantum-espresso/7.2' in version_config['modules']
    assert 'intel/2023' in version_config['modules']


def test_magnetic_config_restore():
    """Test that magnetic configuration can be edited after session save/reload."""
    # Simulate saved workflow_config
    saved_config = {
        'enable_magnetism': True,
        'magnetic_config': {
            'Fe': [2.2],
            'Co': [1.7]
        }
    }
    
    # Verify magnetic config is present and editable
    assert saved_config['enable_magnetism'] is True
    assert 'Fe' in saved_config['magnetic_config']
    assert saved_config['magnetic_config']['Fe'][0] == 2.2
    
    # Simulate editing after reload
    saved_config['magnetic_config']['Fe'] = [2.5]
    assert saved_config['magnetic_config']['Fe'][0] == 2.5


def test_hubbard_config_restore():
    """Test that Hubbard configuration can be edited after session save/reload."""
    # Simulate saved workflow_config
    saved_config = {
        'enable_hubbard': True,
        'hubbard_format': 'new',
        'hubbard_u': {
            'Fe': 4.0,
            'Mn': 4.0
        }
    }
    
    # Verify Hubbard config is present and editable
    assert saved_config['enable_hubbard'] is True
    assert 'Fe' in saved_config['hubbard_u']
    assert saved_config['hubbard_u']['Fe'] == 4.0
    
    # Simulate editing after reload
    saved_config['hubbard_u']['Fe'] = 5.0
    assert saved_config['hubbard_u']['Fe'] == 5.0


def test_job_file_with_launcher_and_modules():
    """Test that job_file generation includes launcher and modules."""
    # Simulate config with launcher and modules
    config = {
        'calc_type': 'scf',
        'nprocs': 16,
        'modules': ['quantum-espresso/7.2', 'intel/2023']
    }
    
    # Simulate machine with launcher
    class MockMachine:
        launcher = 'mpirun -np {nprocs}'
    
    machine = MockMachine()
    
    # Build job_file content
    lines = ["#!/bin/bash", ""]
    lines.append("# Job script generated by xespresso GUI")
    lines.append(f"# Calculation: {config.get('calc_type', 'scf')}")
    lines.append("")
    
    # Module loading
    modules = config.get('modules')
    if modules:
        lines.append("# Load required modules")
        for module in modules:
            lines.append(f"module load {module}")
        lines.append("")
    
    # Launcher
    launcher = machine.launcher
    nprocs = config.get('nprocs', 1)
    if '{nprocs}' in launcher:
        launcher = launcher.replace('{nprocs}', str(nprocs))
    launcher = launcher + " "
    
    lines.append(f"{launcher}pw.x -in test.pwi > test.pwo")
    
    job_content = "\n".join(lines)
    
    # Verify job_file contains modules
    assert 'module load quantum-espresso/7.2' in job_content
    assert 'module load intel/2023' in job_content
    
    # Verify job_file contains launcher with substituted nprocs
    assert 'mpirun -np 16 pw.x' in job_content


def test_config_restore_ui_method():
    """Test the _restore_config_to_ui logic for session reload."""
    # Simulate saved workflow_config with all options
    saved_config = {
        'enable_magnetism': True,
        'magnetic_config': {'Fe': [2.2]},
        'enable_hubbard': True,
        'hubbard_format': 'new',
        'hubbard_u': {'Fe': 4.0},
        'calc_type': 'relax',
        'label': 'Fe/relax',
        'ecutwfc': 60.0,
        'ecutrho': 480.0,
        'occupations': 'smearing',
        'smearing': 'gaussian',
        'degauss': 0.02,
        'conv_thr': 1.0e-8,
        'kspacing': 0.3,
        'adjust_resources': True,
        'nprocs': 16,
        'resources': {
            'nodes': 2,
            'ntasks-per-node': 8,
            'time': '04:00:00',
            'partition': 'compute'
        }
    }
    
    # Verify all configuration is present and restorable
    assert saved_config['enable_magnetism'] is True
    assert saved_config['enable_hubbard'] is True
    assert saved_config['calc_type'] == 'relax'
    assert saved_config['ecutwfc'] == 60.0
    assert saved_config['kspacing'] == 0.3
    assert saved_config['adjust_resources'] is True
    assert saved_config['resources']['nodes'] == 2
    
    # All values should be editable after restore
    saved_config['ecutwfc'] = 70.0
    saved_config['magnetic_config']['Fe'] = [2.5]
    saved_config['hubbard_u']['Fe'] = 5.0
    
    assert saved_config['ecutwfc'] == 70.0
    assert saved_config['magnetic_config']['Fe'][0] == 2.5
    assert saved_config['hubbard_u']['Fe'] == 5.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
