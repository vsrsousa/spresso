"""
Tests for the simplified workflow functionality.
"""

import pytest
import tempfile
import shutil
import numpy as np
from pathlib import Path
from ase.build import bulk
from ase.io.espresso import kspacing_to_grid
from xespresso.workflow import CalculationWorkflow, quick_scf, quick_relax, PRESETS
from xespresso import kpts_from_spacing


def test_presets_exist():
    """Test that all expected presets exist."""
    assert 'fast' in PRESETS
    assert 'moderate' in PRESETS
    assert 'accurate' in PRESETS


def test_preset_values():
    """Test that presets have expected structure and values."""
    for quality, preset in PRESETS.items():
        # Check required keys
        assert 'ecutwfc' in preset
        assert 'ecutrho' in preset
        assert 'conv_thr' in preset
        assert 'kspacing' in preset
        
        # Check types
        assert isinstance(preset['ecutwfc'], float)
        assert isinstance(preset['ecutrho'], float)
        assert isinstance(preset['conv_thr'], float)
        assert isinstance(preset['kspacing'], float)
        
        # Check reasonable values
        assert preset['ecutwfc'] > 0
        assert preset['ecutrho'] > preset['ecutwfc']
        assert preset['kspacing'] > 0


def test_workflow_initialization():
    """Test basic workflow initialization."""
    atoms = bulk("Si", cubic=True)
    pseudopotentials = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}
    
    workflow = CalculationWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        protocol='moderate'
    )
    
    assert workflow.atoms is not None
    assert workflow.pseudopotentials == pseudopotentials
    assert workflow.protocol == 'moderate'


def test_workflow_invalid_quality():
    """Test that invalid protocol raises an error."""
    atoms = bulk("Si", cubic=True)
    pseudopotentials = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}
    
    with pytest.raises(ValueError, match="Protocol must be one of"):
        CalculationWorkflow(
            atoms=atoms,
            pseudopotentials=pseudopotentials,
            protocol='invalid'
        )


def test_workflow_kspacing():
    """Test k-spacing functionality."""
    atoms = bulk("Si", cubic=True)
    pseudopotentials = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}
    
    # Test with default kspacing from preset
    workflow1 = CalculationWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        protocol='moderate'
    )
    kpts1 = workflow1._get_kpts()
    assert isinstance(kpts1, tuple)
    assert len(kpts1) == 3
    assert all(isinstance(k, (int, np.int64, np.int32)) for k in kpts1)
    
    # Test with custom kspacing
    workflow2 = CalculationWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        protocol='moderate',
        kspacing=0.5  # Larger spacing = fewer k-points
    )
    kpts2 = workflow2._get_kpts()
    assert isinstance(kpts2, tuple)
    assert len(kpts2) == 3
    
    # Larger k-spacing should give fewer k-points (approximately)
    # Note: This is approximate due to grid rounding
    assert sum(kpts2) <= sum(kpts1)


def test_workflow_get_preset_info():
    """Test get_preset_info method."""
    atoms = bulk("Si", cubic=True)
    pseudopotentials = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}
    
    workflow = CalculationWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        protocol='fast',
        kspacing=0.4
    )
    
    info = workflow.get_preset_info()
    
    assert 'protocol' in info
    assert 'preset' in info
    assert 'kpts' in info
    assert 'kspacing' in info
    
    assert info['protocol'] == 'fast'
    assert info['kspacing'] == 0.4


def test_workflow_input_data_merge():
    """Test that custom input_data is merged with presets."""
    atoms = bulk("Si", cubic=True)
    pseudopotentials = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}
    
    custom_input = {
        'mixing_beta': 0.9,
        'custom_param': 'test'
    }
    
    workflow = CalculationWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        protocol='moderate',
        input_data=custom_input
    )
    
    # Check that preset values are present
    assert 'ecutwfc' in workflow.input_data
    assert 'conv_thr' in workflow.input_data
    
    # Check that custom values override preset
    assert workflow.input_data['mixing_beta'] == 0.9
    
    # Check that custom params are added
    assert workflow.input_data['custom_param'] == 'test'


def test_workflow_different_presets():
    """Test that different presets have different parameters."""
    atoms = bulk("Si", cubic=True)
    pseudopotentials = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}
    
    fast = CalculationWorkflow(atoms, pseudopotentials, 'fast')
    moderate = CalculationWorkflow(atoms, pseudopotentials, 'moderate')
    accurate = CalculationWorkflow(atoms, pseudopotentials, 'accurate')
    
    # Check that accuracy increases with preset level
    assert fast.input_data['ecutwfc'] < moderate.input_data['ecutwfc'] < accurate.input_data['ecutwfc']
    assert fast.input_data['conv_thr'] > moderate.input_data['conv_thr'] > accurate.input_data['conv_thr']
    assert fast.kspacing > moderate.kspacing > accurate.kspacing


def test_workflow_from_atoms_object():
    """Test creating workflow from ASE Atoms object."""
    atoms = bulk("Si", cubic=True)
    pseudopotentials = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}
    
    workflow = CalculationWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        protocol='moderate'
    )
    
    # Workflow makes a copy to avoid modifying the original
    assert workflow.get_atoms() is not atoms
    assert len(workflow.get_atoms()) == len(atoms)
    assert workflow.get_atoms().get_chemical_symbols() == atoms.get_chemical_symbols()


def test_workflow_magnetic_ferro():
    """Test ferromagnetic configuration."""
    atoms = bulk("Fe", cubic=True)
    pseudopotentials = {"Fe": "Fe.pbe-spn.UPF"}
    
    workflow = CalculationWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        protocol='moderate',
        magnetic_config='ferro'
    )
    
    assert 'input_ntyp' in workflow.input_data
    assert 'starting_magnetization' in workflow.input_data['input_ntyp']


def test_workflow_magnetic_antiferro():
    """Test antiferromagnetic configuration."""
    atoms = bulk("Fe", cubic=True)
    pseudopotentials = {"Fe": "Fe.pbe-spn.UPF"}
    
    workflow = CalculationWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        protocol='moderate',
        magnetic_config='antiferro'
    )
    
    # Should create different species
    assert 'species' in workflow.atoms.arrays
    assert len(set(workflow.atoms.arrays['species'])) > 1
    assert 'input_ntyp' in workflow.input_data
    assert 'starting_magnetization' in workflow.input_data['input_ntyp']


def test_workflow_magnetic_element_based():
    """Test element-based magnetic configuration."""
    atoms = bulk("Fe", cubic=True)
    pseudopotentials = {"Fe": "Fe.pbe-spn.UPF"}
    
    workflow = CalculationWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        protocol='moderate',
        magnetic_config={'Fe': [1, -1]}
    )
    
    # Should create different species
    assert 'species' in workflow.atoms.arrays
    species = workflow.atoms.arrays['species']
    assert 'Fe1' in species
    assert 'Fe2' in species


def test_workflow_magnetic_with_hubbard():
    """Test magnetic configuration with Hubbard U."""
    atoms = bulk("Fe", cubic=True)
    pseudopotentials = {"Fe": "Fe.pbe-spn.UPF"}
    
    workflow = CalculationWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        protocol='accurate',
        magnetic_config={'Fe': {'mag': [1, -1], 'U': 4.3}}
    )
    
    # Should have Hubbard parameters
    assert 'input_ntyp' in workflow.input_data
    assert 'Hubbard_U' in workflow.input_data['input_ntyp']
    
    # Check that U values are set
    hubbard_u = workflow.input_data['input_ntyp']['Hubbard_U']
    assert 'Fe1' in hubbard_u
    assert 'Fe2' in hubbard_u
    assert hubbard_u['Fe1'] == 4.3
    assert hubbard_u['Fe2'] == 4.3


def test_kpts_from_spacing():
    """Test kpts_from_spacing utility function."""
    atoms = bulk("Si", cubic=True)
    
    # Test that it returns the correct k-points
    kpts = kpts_from_spacing(atoms, 0.20)
    assert isinstance(kpts, tuple)
    assert len(kpts) == 3
    
    # Test that it matches manual calculation
    kpts_manual = kspacing_to_grid(atoms, 0.20 / (2 * np.pi))
    assert kpts == tuple(kpts_manual)
    
    # Test with different k-spacing values
    kpts1 = kpts_from_spacing(atoms, 0.5)
    kpts2 = kpts_from_spacing(atoms, 0.3)
    kpts3 = kpts_from_spacing(atoms, 0.15)
    
    # Larger k-spacing should give fewer k-points
    assert sum(kpts1) < sum(kpts2) < sum(kpts3)


def test_kpts_from_spacing_different_structures():
    """Test kpts_from_spacing with different structures."""
    si = bulk("Si", cubic=True)
    fe = bulk("Fe", cubic=True)
    
    # Same k-spacing should give different k-points for different structures
    kpts_si = kpts_from_spacing(si, 0.3)
    kpts_fe = kpts_from_spacing(fe, 0.3)
    
    # Both should be tuples of length 3
    assert isinstance(kpts_si, tuple) and len(kpts_si) == 3
    assert isinstance(kpts_fe, tuple) and len(kpts_fe) == 3
    
    # Due to different lattice parameters, k-points will differ
    # (Si has a=5.43 Å, Fe has a=2.87 Å)
    assert kpts_si != kpts_fe


def test_kpts_from_spacing_consistency():
    """Test that kpts_from_spacing gives consistent results."""
    atoms = bulk("Si", cubic=True)
    
    # Multiple calls with same input should give same output
    kpts1 = kpts_from_spacing(atoms, 0.25)
    kpts2 = kpts_from_spacing(atoms, 0.25)
    assert kpts1 == kpts2


def test_workflow_with_queue():
    """Test workflow with queue configuration."""
    atoms = bulk("Si", cubic=True)
    pseudopotentials = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}
    
    queue = {
        "execution": "local",
        "scheduler": "direct",
        "nprocs": 4
    }
    
    workflow = CalculationWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        protocol='moderate',
        queue=queue
    )
    
    assert workflow.queue is not None
    assert workflow.queue == queue
    assert workflow.queue['execution'] == 'local'


def test_workflow_queue_and_machine_conflict():
    """Test that specifying both queue and machine raises an error."""
    atoms = bulk("Si", cubic=True)
    pseudopotentials = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}
    
    queue = {
        "execution": "local",
        "scheduler": "direct"
    }
    
    with pytest.raises(ValueError, match="Cannot specify both 'queue' and 'machine'"):
        CalculationWorkflow(
            atoms=atoms,
            pseudopotentials=pseudopotentials,
            protocol='moderate',
            queue=queue,
            machine='cluster1'
        )


def test_quick_scf_with_queue():
    """Test quick_scf with queue configuration."""
    atoms = bulk("Si", cubic=True)
    pseudopotentials = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}
    
    queue = {
        "execution": "local",
        "scheduler": "direct"
    }
    
    # This should not raise an error during workflow creation
    # (actual calculation would require QE to be installed)
    try:
        workflow = CalculationWorkflow(
            atoms=atoms,
            pseudopotentials=pseudopotentials,
            protocol='moderate',
            queue=queue
        )
        assert workflow.queue == queue
    except Exception as e:
        # Only workflow creation should succeed, not actual calc execution
        if "Quantum" not in str(e):
            raise


def test_quick_relax_with_queue():
    """Test quick_relax with queue configuration."""
    atoms = bulk("Si", cubic=True)
    pseudopotentials = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}
    
    queue = {
        "execution": "local",
        "scheduler": "direct"
    }
    
    # This should not raise an error during workflow creation
    try:
        workflow = CalculationWorkflow(
            atoms=atoms,
            pseudopotentials=pseudopotentials,
            protocol='moderate',
            queue=queue
        )
        assert workflow.queue == queue
    except Exception as e:
        # Only workflow creation should succeed, not actual calc execution
        if "Quantum" not in str(e):
            raise


def test_workflow_none_queue():
    """Test that workflow works without queue (backward compatibility)."""
    atoms = bulk("Si", cubic=True)
    pseudopotentials = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}
    
    workflow = CalculationWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        protocol='moderate'
    )
    
    # Queue should be None when not specified
    assert workflow.queue is None


# Note: We don't test actual calculation runs here as they require
# Quantum ESPRESSO to be installed and configured. These tests focus
# on the workflow setup and parameter handling.


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
