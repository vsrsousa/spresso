"""
Tests for magnetic configuration helper functions.

These tests verify the new simplified API for setting up magnetic moments
for spin-polarized and antiferromagnetic calculations.
"""

import pytest
import numpy as np
from ase.build import bulk
from ase import Atoms
from xespresso.tools import (
    set_magnetic_moments,
    set_antiferromagnetic,
    set_ferromagnetic,
    setup_magnetic_config
)


class TestSetupMagneticConfig:
    """Tests for the new setup_magnetic_config function."""
    
    def test_simple_equivalent_atoms(self):
        """Test Fe=[1] with 2 Fe atoms - both should be equivalent."""
        atoms = bulk('Fe', cubic=True)  # 2 Fe atoms
        
        config = setup_magnetic_config(atoms, {'Fe': [1]})
        
        # Both atoms should have same species and magnetization
        assert 'input_ntyp' in config
        mag_dict = config['input_ntyp']['starting_magnetization']
        assert len(mag_dict) == 1  # Only one species
        assert list(mag_dict.values())[0] == 1.0
        
        # Check species - should be same for both atoms
        species = set(config['atoms'].arrays['species'])
        assert len(species) == 1  # Both atoms same species
    
    def test_simple_afm_non_equivalent(self):
        """Test Fe=[1, -1] with 2 Fe atoms - should create two species."""
        atoms = bulk('Fe', cubic=True)  # 2 Fe atoms
        
        config = setup_magnetic_config(atoms, {'Fe': [1, -1]})
        
        # Should have two species with opposite magnetization
        mag_dict = config['input_ntyp']['starting_magnetization']
        assert len(mag_dict) == 2
        values = list(mag_dict.values())
        assert 1.0 in values
        assert -1.0 in values
        
        # Check species
        species = set(config['atoms'].arrays['species'])
        assert len(species) == 2  # Two different species
    
    def test_multiple_elements(self):
        """Test multiple elements with different magnetic configs."""
        # Create FeMn compound (simplified)
        atoms = Atoms('Fe2Mn2', positions=[
            [0, 0, 0], [1.5, 0, 0],  # Fe
            [0, 1.5, 0], [1.5, 1.5, 0]  # Mn
        ])
        atoms.cell = [5, 5, 5]
        
        config = setup_magnetic_config(atoms, {
            'Fe': [1],       # Both Fe equivalent
            'Mn': [1, -1]    # Mn antiferromagnetic
        })
        
        mag_dict = config['input_ntyp']['starting_magnetization']
        
        # Should have 3 species: Fe, Mn, Mn1
        assert len(mag_dict) == 3
        
        # Check that we have the right magnetizations
        values = list(mag_dict.values())
        assert values.count(1.0) == 2  # Fe and one Mn
        assert values.count(-1.0) == 1  # One Mn
    
    def test_pattern_replication(self):
        """Test that pattern replicates when fewer moments than atoms."""
        atoms = bulk('Fe', cubic=True) * (2, 1, 1)  # 4 Fe atoms
        
        config = setup_magnetic_config(atoms, {'Fe': [1, -1]})
        
        # Should replicate pattern: 1, -1, 1, -1
        mag_dict = config['input_ntyp']['starting_magnetization']
        assert len(mag_dict) == 2  # Two species
        
        # Count occurrences
        species_list = list(config['atoms'].arrays['species'])
        species_counts = {}
        for sp in species_list:
            species_counts[sp] = species_counts.get(sp, 0) + 1
        
        # Each species should appear twice
        assert all(count == 2 for count in species_counts.values())
    
    def test_expansion_needed_error(self):
        """Test error when expansion needed but not allowed."""
        atoms = bulk('Fe', cubic=True)  # 2 Fe atoms
        
        # Try to specify 4 moments without allowing expansion
        with pytest.raises(ValueError, match="Set expand_cell=True"):
            setup_magnetic_config(atoms, {'Fe': [1, 1, -1, -1]}, expand_cell=False)
    
    def test_expansion_allowed(self):
        """Test cell expansion when needed."""
        atoms = bulk('Fe', cubic=True)  # 2 Fe atoms
        
        config = setup_magnetic_config(
            atoms, 
            {'Fe': [1, 1, -1, -1]},  # Need 4 Fe
            expand_cell=True
        )
        
        assert config['expanded'] is True
        assert len(config['atoms']) == 4  # Should have 4 atoms now
        
        mag_dict = config['input_ntyp']['starting_magnetization']
        assert len(mag_dict) == 2  # Two species (mag=1 and mag=-1)
    
    def test_with_hubbard_u_single_value(self):
        """Test adding Hubbard U with single value."""
        atoms = bulk('Fe', cubic=True)
        
        config = setup_magnetic_config(atoms, {
            'Fe': {'mag': [1, -1], 'U': 4.3}
        })
        
        # Check magnetization
        mag_dict = config['input_ntyp']['starting_magnetization']
        assert len(mag_dict) == 2
        
        # Check Hubbard U
        assert 'Hubbard_U' in config['input_ntyp']
        u_dict = config['input_ntyp']['Hubbard_U']
        assert len(u_dict) == 2  # Both species should have U
        assert all(u == 4.3 for u in u_dict.values())
    
    def test_with_hubbard_u_different_values(self):
        """Test adding different Hubbard U for each species."""
        atoms = bulk('Fe', cubic=True)
        
        config = setup_magnetic_config(atoms, {
            'Fe': {'mag': [1, -1], 'U': [4.3, 4.5]}
        })
        
        # Check Hubbard U
        u_dict = config['input_ntyp']['Hubbard_U']
        assert len(u_dict) == 2
        values = list(u_dict.values())
        assert 4.3 in values
        assert 4.5 in values
    
    def test_complex_femnal_system(self):
        """Test complex FeMnAl2 system as mentioned in the issue."""
        # Simulate FeMnAl2 with 2 Fe, 2 Mn, 4 Al
        atoms = Atoms('Fe2Mn2Al4', positions=[
            [0, 0, 0], [1, 0, 0],           # Fe
            [0, 1, 0], [1, 1, 0],           # Mn
            [0, 0, 1], [1, 0, 1],           # Al
            [0, 1, 1], [1, 1, 1]            # Al
        ])
        atoms.cell = [5, 5, 5]
        
        config = setup_magnetic_config(atoms, {
            'Fe': [1],        # Both Fe equivalent
            'Mn': [1, -1],    # Mn AFM
            'Al': [0]         # Al non-magnetic
        })
        
        mag_dict = config['input_ntyp']['starting_magnetization']
        
        # Only Fe and Mn should be in starting_magnetization (Al has 0)
        assert len(mag_dict) == 3  # Fe, Mn, Mn1
        
        # Verify magnetizations
        assert any(v == 1.0 for v in mag_dict.values())  # Fe
        assert any(v == -1.0 for v in mag_dict.values())  # One Mn
    
    def test_element_not_in_structure(self):
        """Test error when specified element not in structure."""
        atoms = bulk('Fe', cubic=True)
        
        with pytest.raises(ValueError, match="Element Mn not found"):
            setup_magnetic_config(atoms, {'Mn': [1]})
    
    def test_with_pseudopotentials(self):
        """Test with existing pseudopotentials."""
        atoms = bulk('Fe', cubic=True)
        
        base_pseudo = {'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'}
        
        config = setup_magnetic_config(
            atoms, 
            {'Fe': [1, -1]},
            pseudopotentials=base_pseudo
        )
        
        # Check pseudopotentials are updated
        pseudo = config['pseudopotentials']
        assert 'Fe1' in pseudo
        assert 'Fe2' in pseudo
        assert pseudo['Fe1'] == 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
        assert pseudo['Fe2'] == 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
    
    def test_simple_value_not_list(self):
        """Test that single values work (not in list)."""
        atoms = bulk('Fe', cubic=True)
        
        config = setup_magnetic_config(atoms, {'Fe': 1.5})
        
        mag_dict = config['input_ntyp']['starting_magnetization']
        assert len(mag_dict) == 1
        assert list(mag_dict.values())[0] == 1.5
    
    def test_new_hubbard_format_with_orbital(self):
        """Test new Hubbard format (QE 7.x) with orbital specification."""
        atoms = bulk('Fe', cubic=True)
        
        config = setup_magnetic_config(atoms, {
            'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}
        }, qe_version='7.2')
        
        # Check that hubbard dict is created (new format)
        assert 'hubbard' in config
        assert config['hubbard_format'] == 'new'
        
        # Check U parameters with orbital
        hubbard = config['hubbard']
        assert 'u' in hubbard
        assert 'Fe1-3d' in hubbard['u']
        assert 'Fe2-3d' in hubbard['u']
        assert hubbard['u']['Fe1-3d'] == 4.3
        assert hubbard['u']['Fe2-3d'] == 4.3
    
    def test_new_hubbard_format_different_u_per_species(self):
        """Test new format with different U for each species."""
        atoms = bulk('Fe', cubic=True)
        
        config = setup_magnetic_config(atoms, {
            'Fe': {'mag': [1, -1], 'U': {'3d': [4.3, 4.5]}}
        }, qe_version='7.2')
        
        hubbard = config['hubbard']
        assert hubbard['u']['Fe1-3d'] == 4.3
        assert hubbard['u']['Fe2-3d'] == 4.5
    
    def test_new_hubbard_format_with_v_parameter(self):
        """Test new format with V parameter (inter-site interaction)."""
        atoms = Atoms('FeO', positions=[[0, 0, 0], [1.5, 0, 0]])
        atoms.cell = [5, 5, 5]
        
        config = setup_magnetic_config(atoms, {
            'Fe': {
                'mag': [1],
                'U': {'3d': 4.3},
                'V': [{'species2': 'O', 'orbital1': '3d', 'orbital2': '2p', 'value': 1.0}]
            },
            'O': [0]
        }, qe_version='7.2')
        
        hubbard = config['hubbard']
        assert 'v' in hubbard
        assert len(hubbard['v']) == 1
        v_param = hubbard['v'][0]
        assert v_param['species1'] == 'Fe'
        assert v_param['species2'] == 'O'
        assert v_param['orbital1'] == '3d'
        assert v_param['orbital2'] == '2p'
        assert v_param['value'] == 1.0
    
    def test_old_format_with_dict_u_warns(self):
        """Test that using dict U format with old QE version extracts value."""
        atoms = bulk('Fe', cubic=True)
        
        # Should work but print warning
        config = setup_magnetic_config(atoms, {
            'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}
        }, qe_version='6.8')
        
        # Should fall back to old format
        assert config['hubbard_format'] == 'old'
        assert 'Hubbard_U' in config['input_ntyp']
        assert config['input_ntyp']['Hubbard_U']['Fe1'] == 4.3
        assert config['input_ntyp']['Hubbard_U']['Fe2'] == 4.3
    
    def test_new_format_error_without_orbital(self):
        """Test that new format requires orbital specification."""
        atoms = bulk('Fe', cubic=True)
        
        with pytest.raises(ValueError, match="requires orbital specification"):
            setup_magnetic_config(atoms, {
                'Fe': {'mag': [1, -1], 'U': 4.3}  # No orbital
            }, hubbard_format='new')
    
    def test_custom_projector_parameter(self):
        """Test that custom projector can be specified."""
        atoms = bulk('Fe', cubic=True)
        
        # Test with atomic projector
        config = setup_magnetic_config(atoms, {
            'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}
        }, qe_version='7.2', projector='atomic')
        
        assert config['hubbard']['projector'] == 'atomic'
        
        # Test with norm-atomic projector
        config = setup_magnetic_config(atoms, {
            'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}
        }, qe_version='7.2', projector='norm-atomic')
        
        assert config['hubbard']['projector'] == 'norm-atomic'
    
    def test_default_projector_is_ortho_atomic(self):
        """Test that default projector is ortho-atomic (QE recommended)."""
        atoms = bulk('Fe', cubic=True)
        
        config = setup_magnetic_config(atoms, {
            'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}
        }, qe_version='7.2')
        
        # Default should be ortho-atomic as recommended by QE
        assert config['hubbard']['projector'] == 'ortho-atomic'
    
    def test_invalid_projector_raises_error(self):
        """Test that invalid projector value raises an error."""
        atoms = bulk('Fe', cubic=True)
        
        with pytest.raises(ValueError, match="Invalid projector"):
            setup_magnetic_config(atoms, {
                'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}
            }, qe_version='7.2', projector='invalid-projector')


class TestSetMagneticMoments:
    """Tests for set_magnetic_moments function."""
    
    def test_simple_afm_list(self):
        """Test simple AFM configuration with list of magnetic moments."""
        atoms = bulk('Fe', cubic=True)
        
        mag_config = set_magnetic_moments(atoms, [1.0, -1.0])
        
        # Check that species were created
        assert 'species' in atoms.arrays
        assert len(set(atoms.arrays['species'])) == 2  # Two different species
        
        # Check input_ntyp structure
        assert 'input_ntyp' in mag_config
        assert 'starting_magnetization' in mag_config['input_ntyp']
        
        # Check that we have two species with opposite magnetization
        mag_dict = mag_config['input_ntyp']['starting_magnetization']
        assert len(mag_dict) == 2
        values = list(mag_dict.values())
        assert 1.0 in values
        assert -1.0 in values
        
        # Check pseudopotentials
        assert 'pseudopotentials' in mag_config
        assert len(mag_config['pseudopotentials']) == 2
    
    def test_afm_with_dict(self):
        """Test AFM configuration with dictionary of magnetic moments."""
        atoms = bulk('Fe', cubic=True)
        
        mag_config = set_magnetic_moments(atoms, {0: 1.0, 1: -1.0})
        
        # Check input_ntyp
        mag_dict = mag_config['input_ntyp']['starting_magnetization']
        assert len(mag_dict) == 2
        assert 1.0 in mag_dict.values()
        assert -1.0 in mag_dict.values()
    
    def test_mixed_elements(self):
        """Test with multiple elements."""
        from ase import Atoms
        atoms = Atoms('Fe2O', positions=[[0, 0, 0], [1.5, 0, 0], [0, 1.5, 0]])
        atoms.cell = [5, 5, 5]
        
        # Set magnetic moments only on Fe atoms
        mag_config = set_magnetic_moments(atoms, {0: 1.0, 1: -1.0, 2: 0.0})
        
        # Check that we have Fe and Fe1 species with magnetization
        mag_dict = mag_config['input_ntyp']['starting_magnetization']
        assert len(mag_dict) == 2  # Only Fe atoms with non-zero magnetization
    
    def test_with_pseudopotentials(self):
        """Test with existing pseudopotentials."""
        atoms = bulk('Fe', cubic=True)
        existing_pseudo = {'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'}
        
        mag_config = set_magnetic_moments(
            atoms, 
            [1.0, -1.0],
            pseudopotentials=existing_pseudo
        )
        
        # Check that pseudopotentials are updated correctly
        pseudo = mag_config['pseudopotentials']
        assert 'Fe1' in pseudo
        assert 'Fe2' in pseudo
        # Both should use the same pseudopotential file
        assert pseudo['Fe1'] == 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
        assert pseudo['Fe2'] == 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
    
    def test_zero_magnetic_moments(self):
        """Test atoms with zero magnetic moments are not included."""
        atoms = bulk('Fe', cubic=True) * (2, 1, 1)
        
        # Only first two atoms have magnetic moments
        mag_config = set_magnetic_moments(atoms, {0: 1.0, 1: -1.0, 2: 0.0, 3: 0.0})
        
        # Only non-zero magnetic moments should be in starting_magnetization
        mag_dict = mag_config['input_ntyp']['starting_magnetization']
        assert len(mag_dict) == 2


class TestSetAntiferromagnetic:
    """Tests for set_antiferromagnetic function."""
    
    def test_simple_afm(self):
        """Test simple antiferromagnetic configuration."""
        atoms = bulk('Fe', cubic=True)
        
        afm_config = set_antiferromagnetic(atoms, [[0], [1]])
        
        # Check magnetization
        mag_dict = afm_config['input_ntyp']['starting_magnetization']
        assert len(mag_dict) == 2
        values = list(mag_dict.values())
        assert 1.0 in values
        assert -1.0 in values
    
    def test_custom_magnetic_moment(self):
        """Test AFM with custom magnetic moment magnitude."""
        atoms = bulk('Fe', cubic=True)
        
        afm_config = set_antiferromagnetic(atoms, [[0], [1]], magnetic_moment=2.5)
        
        # Check magnetization magnitude
        mag_dict = afm_config['input_ntyp']['starting_magnetization']
        values = list(mag_dict.values())
        assert 2.5 in values
        assert -2.5 in values
    
    def test_larger_system(self):
        """Test AFM in larger system with multiple sublattices."""
        atoms = bulk('Fe', cubic=True) * (2, 2, 1)
        
        # Checkerboard AFM pattern
        afm_config = set_antiferromagnetic(
            atoms, 
            [[0, 3], [1, 2]],
            magnetic_moment=1.0
        )
        
        # Should have two species
        mag_dict = afm_config['input_ntyp']['starting_magnetization']
        assert len(mag_dict) == 2
    
    def test_invalid_sublattices(self):
        """Test that error is raised for invalid sublattice count."""
        atoms = bulk('Fe', cubic=True)
        
        with pytest.raises(ValueError, match="exactly 2 sublattices"):
            set_antiferromagnetic(atoms, [[0]])
        
        with pytest.raises(ValueError, match="exactly 2 sublattices"):
            set_antiferromagnetic(atoms, [[0], [1], [2]])


class TestSetFerromagnetic:
    """Tests for set_ferromagnetic function."""
    
    def test_simple_fm(self):
        """Test simple ferromagnetic configuration."""
        atoms = bulk('Fe', cubic=True)
        
        fm_config = set_ferromagnetic(atoms, magnetic_moment=2.0)
        
        # All atoms should have same magnetization
        mag_dict = fm_config['input_ntyp']['starting_magnetization']
        assert len(mag_dict) == 1  # Only one species
        assert list(mag_dict.values())[0] == 2.0
    
    def test_element_specific(self):
        """Test ferromagnetic configuration for specific element."""
        from ase import Atoms
        atoms = Atoms('Fe2O', positions=[[0, 0, 0], [1.5, 0, 0], [0, 1.5, 0]])
        atoms.cell = [5, 5, 5]
        
        # Only set magnetization for Fe
        fm_config = set_ferromagnetic(atoms, magnetic_moment=2.0, element='Fe')
        
        # Only Fe atoms should have magnetization
        mag_dict = fm_config['input_ntyp']['starting_magnetization']
        assert len(mag_dict) == 1
        assert list(mag_dict.values())[0] == 2.0
        
        # O should not have magnetization in input_ntyp
        species_with_mag = list(mag_dict.keys())
        assert all('Fe' in sp for sp in species_with_mag)


class TestIntegrationWithEspresso:
    """Integration tests to ensure helpers work with Espresso calculator."""
    
    def test_afm_fe_integration(self):
        """Test that AFM configuration works with Espresso calculator setup."""
        atoms = bulk('Fe', cubic=True)
        
        # Use helper function
        mag_config = set_antiferromagnetic(atoms, [[0], [1]], magnetic_moment=1.0)
        
        # Add actual pseudopotential names
        mag_config['pseudopotentials']['Fe1'] = 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
        mag_config['pseudopotentials']['Fe2'] = 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
        
        # Check that the configuration is correct
        assert 'input_ntyp' in mag_config
        assert 'starting_magnetization' in mag_config['input_ntyp']
        assert 'Fe1' in mag_config['pseudopotentials']
        assert 'Fe2' in mag_config['pseudopotentials']
        
        # Check magnetization values
        mag_dict = mag_config['input_ntyp']['starting_magnetization']
        assert len(mag_dict) == 2
        values = list(mag_dict.values())
        assert 1.0 in values
        assert -1.0 in values
    
    def test_comparison_with_old_method(self):
        """Compare new helper with old manual method."""
        atoms_old = bulk('Fe', cubic=True)
        atoms_new = bulk('Fe', cubic=True)
        
        # Old method (manual) - updated to match new numbering convention
        atoms_old.new_array('species', np.array(atoms_old.get_chemical_symbols(), dtype='U20'))
        atoms_old.arrays['species'][0] = 'Fe1'
        atoms_old.arrays['species'][1] = 'Fe2'
        input_ntyp_old = {
            'starting_magnetization': {
                'Fe1': 1.0,
                'Fe2': -1.0,
            }
        }
        
        # New method (helper)
        mag_config = set_magnetic_moments(atoms_new, [1.0, -1.0])
        
        # Both should have same number of species
        assert len(set(atoms_old.arrays['species'])) == len(set(atoms_new.arrays['species']))
        
        # Both should have same magnetization values
        old_mags = set(input_ntyp_old['starting_magnetization'].values())
        new_mags = set(mag_config['input_ntyp']['starting_magnetization'].values())
        assert old_mags == new_mags


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
