"""
Integration test demonstrating the solution to the problem statement.

This test shows that users can now pass the entire config dict from
setup_magnetic_config directly to Espresso, without having to manually
extract and pass pseudopotentials separately.
"""

import pytest
import tempfile
import os
from ase.build import bulk
from xespresso.tools import setup_magnetic_config
from xespresso.xio import write_espresso_in


class TestProblemStatementSolution:
    """Integration tests showing the problem statement is solved."""
    
    def test_old_way_still_works(self):
        """Test that the old way of passing parameters still works."""
        atoms = bulk('Fe', cubic=True)
        base_pseudos = {'Fe': 'Fe.pbe-spn.UPF'}
        
        config = setup_magnetic_config(atoms, {'Fe': [1, -1]}, 
                                      pseudopotentials=base_pseudos)
        
        # Old way: manually extract and pass separately
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, 'test.pwi')
            
            write_espresso_in(
                input_file,
                config['atoms'],
                pseudopotentials=config['pseudopotentials'],  # Manual extraction
                input_data={'input_ntyp': config['input_ntyp']},
                ecutwfc=30.0,
                nspin=2,
                kpts=(4, 4, 4)
            )
            
            # Verify file was created
            assert os.path.exists(input_file)
            
            # Read and verify content
            with open(input_file, 'r') as f:
                content = f.read()
            
            assert 'ATOMIC_SPECIES' in content
            assert 'Fe1' in content
            assert 'Fe2' in content
    
    def test_new_way_with_input_data(self):
        """Test new way: pass pseudopotentials inside input_data."""
        atoms = bulk('Fe', cubic=True)
        base_pseudos = {'Fe': 'Fe.pbe-spn.UPF'}
        
        config = setup_magnetic_config(atoms, {'Fe': [1, -1]}, 
                                      pseudopotentials=base_pseudos)
        
        # New way: pass pseudopotentials inside input_data
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, 'test.pwi')
            
            write_espresso_in(
                input_file,
                config['atoms'],
                input_data={
                    'input_ntyp': config['input_ntyp'],
                    'pseudopotentials': config['pseudopotentials']  # Inside input_data
                },
                ecutwfc=30.0,
                nspin=2,
                kpts=(4, 4, 4)
            )
            
            # Verify file was created
            assert os.path.exists(input_file)
            
            # Read and verify content
            with open(input_file, 'r') as f:
                content = f.read()
            
            assert 'ATOMIC_SPECIES' in content
            assert 'Fe1' in content
            assert 'Fe2' in content
    
    def test_pass_entire_config_dict(self):
        """Test passing entire config dict as input_data (most convenient)."""
        atoms = bulk('Fe', cubic=True)
        base_pseudos = {'Fe': 'Fe.pbe-spn.UPF'}
        
        config = setup_magnetic_config(atoms, {'Fe': [1, -1]}, 
                                      pseudopotentials=base_pseudos)
        
        # Most convenient: pass entire config
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, 'test.pwi')
            
            write_espresso_in(
                input_file,
                config['atoms'],
                input_data=config,  # Pass entire config dict
                ecutwfc=30.0,
                nspin=2,
                kpts=(4, 4, 4)
            )
            
            # Verify file was created
            assert os.path.exists(input_file)
            
            # Read and verify content
            with open(input_file, 'r') as f:
                content = f.read()
            
            assert 'ATOMIC_SPECIES' in content
            assert 'Fe1' in content
            assert 'Fe2' in content
    
    def test_with_hubbard_parameters(self):
        """Test that Hubbard parameters work with new approach."""
        atoms = bulk('Fe', cubic=True)
        base_pseudos = {'Fe': 'Fe.pbe-spn.UPF'}
        
        config = setup_magnetic_config(atoms, {
            'Fe': {'mag': [1, -1], 'U': 4.3}
        }, pseudopotentials=base_pseudos)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, 'test.pwi')
            
            write_espresso_in(
                input_file,
                config['atoms'],
                input_data=config,  # Pass entire config
                ecutwfc=30.0,
                nspin=2,
                lda_plus_u=True,
                kpts=(4, 4, 4)
            )
            
            # Verify file was created
            assert os.path.exists(input_file)
            
            # Read and verify content
            with open(input_file, 'r') as f:
                content = f.read()
            
            assert 'ATOMIC_SPECIES' in content
            assert 'Fe1' in content
            assert 'Fe2' in content
            # Hubbard parameters should be present
            assert 'Hubbard' in content or 'hubbard' in content or '4.3' in content
    
    def test_auto_derive_from_base_element(self):
        """Test that derived species auto-inherit pseudopotentials from base element."""
        atoms = bulk('Fe', cubic=True)
        # Only provide base element pseudopotential
        base_pseudos = {'Fe': 'Fe.pbe-spn.UPF'}
        
        config = setup_magnetic_config(atoms, {'Fe': [1, -1]}, 
                                      pseudopotentials=base_pseudos)
        
        # Verify that derived species got the pseudopotential
        assert 'Fe1' in config['pseudopotentials']
        assert 'Fe2' in config['pseudopotentials']
        assert config['pseudopotentials']['Fe1'] == 'Fe.pbe-spn.UPF'
        assert config['pseudopotentials']['Fe2'] == 'Fe.pbe-spn.UPF'
        
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, 'test.pwi')
            
            write_espresso_in(
                input_file,
                config['atoms'],
                input_data=config,
                ecutwfc=30.0,
                nspin=2,
                kpts=(4, 4, 4)
            )
            
            assert os.path.exists(input_file)
            
            with open(input_file, 'r') as f:
                content = f.read()
            
            # Both species should be in ATOMIC_SPECIES
            assert 'Fe1' in content
            assert 'Fe2' in content
            assert content.count('Fe.pbe-spn.UPF') >= 2  # Both use same UPF
    
    def test_user_override_with_top_level_pseudopotentials(self):
        """Test that top-level pseudopotentials override input_data ones."""
        atoms = bulk('Fe', cubic=True)
        base_pseudos = {'Fe': 'Fe.pbe-spn.UPF'}
        
        config = setup_magnetic_config(atoms, {'Fe': [1, -1]}, 
                                      pseudopotentials=base_pseudos)
        
        # User wants to override one pseudopotential
        custom_pseudos = {
            'Fe1': 'Fe_custom.UPF'  # Override just Fe1
        }
        
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, 'test.pwi')
            
            write_espresso_in(
                input_file,
                config['atoms'],
                pseudopotentials=custom_pseudos,  # Top-level override
                input_data=config,  # Has Fe1 and Fe2 from config
                ecutwfc=30.0,
                nspin=2,
                kpts=(4, 4, 4)
            )
            
            assert os.path.exists(input_file)
            
            with open(input_file, 'r') as f:
                content = f.read()
            
            # Fe1 should use custom, Fe2 should use base
            assert 'Fe_custom.UPF' in content
            assert 'Fe.pbe-spn.UPF' in content


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
