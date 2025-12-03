"""
Test that calculators are properly attached to atoms objects.

This test verifies the fix for the issue:
"RuntimeError: Atoms object has no calculator."

The fix ensures that when prepare_calculation_from_gui() is called,
the calculator is attached to the atoms object following xespresso's pattern.
"""

import pytest
import tempfile
import os
from ase.build import bulk
from gui.calculations import prepare_calculation_from_gui


class TestCalculatorAttachment:
    """Test that calculators are properly attached to atoms."""
    
    def test_prepare_calculation_attaches_calculator(self):
        """
        Test that prepare_calculation_from_gui attaches calculator to atoms.
        
        This is the main fix - following xespresso's pattern where
        atoms.calc = calc is required before calling atoms.get_potential_energy()
        """
        atoms = bulk("Al", cubic=True)
        config = {
            'pseudopotentials': {'Al': 'Al.pbe.UPF'},
            'ecutwfc': 30.0,
            'kpts': (2, 2, 2),
            'calc_type': 'scf'
        }
        
        with tempfile.TemporaryDirectory() as tmpdir:
            label = os.path.join(tmpdir, 'test/Al')
            prepared_atoms, calc = prepare_calculation_from_gui(atoms, config, label=label)
            
            # Verify calculator is attached to atoms
            assert hasattr(prepared_atoms, 'calc'), "Atoms object should have calc attribute"
            assert prepared_atoms.calc is not None, "Calculator should be attached to atoms"
            assert prepared_atoms.calc == calc, "Attached calculator should be the same as returned calculator"
    
    def test_calculator_attached_for_scf(self):
        """Test calculator attachment for SCF calculation."""
        atoms = bulk("Fe", cubic=True)
        config = {
            'pseudopotentials': {'Fe': 'Fe.pbe.UPF'},
            'ecutwfc': 40.0,
            'kpts': (3, 3, 3),
            'calc_type': 'scf',
            'enable_magnetism': False,
            'enable_hubbard': False
        }
        
        with tempfile.TemporaryDirectory() as tmpdir:
            label = os.path.join(tmpdir, 'scf/Fe')
            prepared_atoms, calc = prepare_calculation_from_gui(atoms, config, label=label)
            
            # Verify calculator is attached
            assert prepared_atoms.calc is not None
            assert prepared_atoms.calc == calc
    
    def test_calculator_attached_for_relax(self):
        """Test calculator attachment for relax calculation."""
        atoms = bulk("Cu", cubic=True)
        config = {
            'pseudopotentials': {'Cu': 'Cu.pbe.UPF'},
            'ecutwfc': 35.0,
            'kpts': (4, 4, 4),
            'calc_type': 'relax',
            'forc_conv_thr': 1e-3
        }
        
        with tempfile.TemporaryDirectory() as tmpdir:
            label = os.path.join(tmpdir, 'relax/Cu')
            prepared_atoms, calc = prepare_calculation_from_gui(atoms, config, label=label)
            
            # Verify calculator is attached
            assert prepared_atoms.calc is not None
            assert prepared_atoms.calc == calc
    
    def test_calculator_attached_with_magnetism(self):
        """Test calculator attachment when magnetism is enabled."""
        atoms = bulk("Fe", cubic=True)
        config = {
            'pseudopotentials': {'Fe': 'Fe.pbe.UPF'},
            'ecutwfc': 40.0,
            'kpts': (3, 3, 3),
            'calc_type': 'scf',
            'enable_magnetism': True,
            'magnetic_config': {'Fe': [1.0]},
            'enable_hubbard': False
        }
        
        with tempfile.TemporaryDirectory() as tmpdir:
            label = os.path.join(tmpdir, 'scf/Fe-mag')
            prepared_atoms, calc = prepare_calculation_from_gui(atoms, config, label=label)
            
            # Verify calculator is attached even with magnetic configuration
            assert prepared_atoms.calc is not None
            assert prepared_atoms.calc == calc
    
    def test_calculator_ready_for_get_potential_energy(self):
        """
        Test that atoms.get_potential_energy() won't raise "no calculator" error.
        
        We don't actually run the calculation (no QE installed), but we verify
        that the calculator is properly attached so that when QE is available,
        get_potential_energy() will work.
        """
        atoms = bulk("Si", cubic=True)
        config = {
            'pseudopotentials': {'Si': 'Si.pbe.UPF'},
            'ecutwfc': 30.0,
            'kpts': (2, 2, 2),
            'calc_type': 'scf'
        }
        
        with tempfile.TemporaryDirectory() as tmpdir:
            label = os.path.join(tmpdir, 'test/Si')
            prepared_atoms, calc = prepare_calculation_from_gui(atoms, config, label=label)
            
            # Verify that get_potential_energy() can be called
            # (it will fail because QE is not installed, but it shouldn't fail 
            # because of missing calculator)
            assert hasattr(prepared_atoms, 'get_potential_energy')
            assert prepared_atoms.calc is not None
            
            # Try to call get_potential_energy - it should fail with a different error
            # not "Atoms object has no calculator"
            try:
                energy = prepared_atoms.get_potential_energy()
            except RuntimeError as e:
                # Should not be "Atoms object has no calculator"
                assert "Atoms object has no calculator" not in str(e)
            except Exception as e:
                # Other exceptions are fine (e.g., calculator not finding QE executable)
                pass


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
