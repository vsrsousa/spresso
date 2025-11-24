"""
Tests for automatic pseudopotential extraction from input_data.

This test suite verifies that the Espresso calculator can extract
pseudopotentials from input_data when using setup_magnetic_config,
avoiding the need for users to manually pass pseudopotentials as a
separate parameter.
"""

import pytest
import tempfile
import os
import numpy as np
from ase.build import bulk
from ase import Atoms
from xespresso.tools import setup_magnetic_config
from xespresso.xio import sort_qe_input


class TestPseudopotentialExtraction:
    """Tests for pseudopotential extraction from input_data."""
    
    def test_extract_pseudopotentials_from_input_data(self):
        """Test that pseudopotentials are extracted from input_data."""
        parameters = {
            'input_data': {
                'pseudopotentials': {
                    'Fe1': 'Fe.pbe-spn.UPF',
                    'Fe2': 'Fe.pbe-spn.UPF'
                },
                'input_ntyp': {
                    'starting_magnetization': {
                        'Fe1': 1.0,
                        'Fe2': -1.0
                    }
                }
            },
            'ecutwfc': 30.0,
            'nspin': 2
        }
        
        sorted_params, unused = sort_qe_input(parameters)
        
        # Pseudopotentials should be extracted to top level
        assert 'pseudopotentials' in sorted_params
        assert sorted_params['pseudopotentials']['Fe1'] == 'Fe.pbe-spn.UPF'
        assert sorted_params['pseudopotentials']['Fe2'] == 'Fe.pbe-spn.UPF'
        
        # Should be removed from input_data
        assert 'pseudopotentials' not in sorted_params['input_data']
        
        # Other parameters should remain (normalized to uppercase)
        assert 'INPUT_NTYP' in sorted_params['input_data']
    
    def test_merge_pseudopotentials_from_multiple_sources(self):
        """Test that pseudopotentials from input_data merge with top-level ones."""
        parameters = {
            'pseudopotentials': {
                'Fe1': 'Fe1_custom.UPF'  # Top-level takes precedence
            },
            'input_data': {
                'pseudopotentials': {
                    'Fe1': 'Fe.pbe-spn.UPF',  # Should be overridden
                    'Fe2': 'Fe.pbe-spn.UPF'   # Should be added
                }
            }
        }
        
        sorted_params, unused = sort_qe_input(parameters)
        
        # Top-level should take precedence
        assert sorted_params['pseudopotentials']['Fe1'] == 'Fe1_custom.UPF'
        # Input_data should add missing ones
        assert sorted_params['pseudopotentials']['Fe2'] == 'Fe.pbe-spn.UPF'
    
    def test_species_map_auto_derivation(self):
        """Test that species_map helps auto-populate pseudopotentials."""
        parameters = {
            'pseudopotentials': {
                'Fe': 'Fe.pbe-spn.UPF'  # Only base element
            },
            'input_data': {
                'species_map': {
                    'Fe1': 'Fe',
                    'Fe2': 'Fe'
                },
                'input_ntyp': {
                    'starting_magnetization': {
                        'Fe1': 1.0,
                        'Fe2': -1.0
                    }
                }
            }
        }
        
        sorted_params, unused = sort_qe_input(parameters)
        
        # Derived species should inherit from base element
        assert sorted_params['pseudopotentials']['Fe'] == 'Fe.pbe-spn.UPF'
        assert sorted_params['pseudopotentials']['Fe1'] == 'Fe.pbe-spn.UPF'
        assert sorted_params['pseudopotentials']['Fe2'] == 'Fe.pbe-spn.UPF'
        
        # species_map should be removed from input_data
        assert 'species_map' not in sorted_params['input_data']
    
    def test_full_config_dict_from_setup_magnetic_config(self):
        """Test passing full config dict from setup_magnetic_config."""
        atoms = bulk('Fe', cubic=True)
        base_pseudos = {'Fe': 'Fe.pbe-spn.UPF'}
        
        config = setup_magnetic_config(atoms, {'Fe': [1, -1]}, 
                                      pseudopotentials=base_pseudos)
        
        # User passes entire config as input_data
        parameters = {
            'input_data': config,
            'ecutwfc': 30.0,
            'nspin': 2
        }
        
        sorted_params, unused = sort_qe_input(parameters)
        
        # Pseudopotentials should be extracted
        assert 'pseudopotentials' in sorted_params
        assert 'Fe1' in sorted_params['pseudopotentials']
        assert 'Fe2' in sorted_params['pseudopotentials']
        
        # Metadata fields should be removed
        assert 'expanded' not in sorted_params['input_data']
        assert 'hubbard_format' not in sorted_params['input_data']
        assert 'atoms' not in sorted_params['input_data']
        assert 'species_map' not in sorted_params['input_data']
        assert 'pseudopotentials' not in sorted_params['input_data']
        
        # input_ntyp should remain (normalized to uppercase)
        assert 'INPUT_NTYP' in sorted_params['input_data']
    
    def test_setup_magnetic_config_with_species_map(self):
        """Test that setup_magnetic_config properly creates species_map."""
        atoms = bulk('Fe', cubic=True)
        base_pseudos = {'Fe': 'Fe.pbe-spn.UPF'}
        
        config = setup_magnetic_config(atoms, {'Fe': [1, -1]},
                                      pseudopotentials=base_pseudos)
        
        # Verify config structure
        assert 'species_map' in config
        assert 'pseudopotentials' in config
        assert 'input_ntyp' in config
        
        # Verify species_map content
        assert config['species_map']['Fe1'] == 'Fe'
        assert config['species_map']['Fe2'] == 'Fe'
        
        # Verify pseudopotentials are properly derived
        assert config['pseudopotentials']['Fe1'] == 'Fe.pbe-spn.UPF'
        assert config['pseudopotentials']['Fe2'] == 'Fe.pbe-spn.UPF'
    
    def test_mixed_elements_with_base_pseudopotentials(self):
        """Test mixed elements where only base pseudopotentials are provided."""
        atoms = Atoms('Fe2Mn2', positions=[
            [0, 0, 0], [1.5, 0, 0],
            [0, 1.5, 0], [1.5, 1.5, 0]
        ])
        atoms.cell = [5, 5, 5]
        
        base_pseudos = {
            'Fe': 'Fe.pbe-spn.UPF',
            'Mn': 'Mn.pbe-spn.UPF'
        }
        
        config = setup_magnetic_config(atoms, {
            'Fe': [1],      # Both Fe equivalent
            'Mn': [1, -1]   # Mn AFM
        }, pseudopotentials=base_pseudos)
        
        # Now use with Espresso via input_data
        parameters = {
            'input_data': {
                'pseudopotentials': config['pseudopotentials'],
                'species_map': config['species_map'],
                'input_ntyp': config['input_ntyp']
            }
        }
        
        sorted_params, unused = sort_qe_input(parameters)
        
        # Should have all species pseudopotentials
        assert 'Fe' in sorted_params['pseudopotentials']
        assert 'Mn1' in sorted_params['pseudopotentials']
        assert 'Mn2' in sorted_params['pseudopotentials']
        
        # Verify they're properly derived from base elements
        assert sorted_params['pseudopotentials']['Fe'] == 'Fe.pbe-spn.UPF'
        assert sorted_params['pseudopotentials']['Mn1'] == 'Mn.pbe-spn.UPF'
        assert sorted_params['pseudopotentials']['Mn2'] == 'Mn.pbe-spn.UPF'
    
    def test_user_can_pass_only_base_pseudopotentials(self):
        """Test that user can pass only base element pseudopotentials."""
        parameters = {
            'pseudopotentials': {
                'Fe': 'Fe.pbe-spn.UPF'  # Only base element
            },
            'input_data': {
                'species_map': {
                    'Fe': 'Fe',   # Base maps to itself
                    'Fe1': 'Fe',  # Derived species
                    'Fe2': 'Fe'
                }
            }
        }
        
        sorted_params, unused = sort_qe_input(parameters)
        
        # All species should have pseudopotentials
        assert 'Fe' in sorted_params['pseudopotentials']
        assert 'Fe1' in sorted_params['pseudopotentials']
        assert 'Fe2' in sorted_params['pseudopotentials']
        
        # All should use the same file
        assert sorted_params['pseudopotentials']['Fe'] == 'Fe.pbe-spn.UPF'
        assert sorted_params['pseudopotentials']['Fe1'] == 'Fe.pbe-spn.UPF'
        assert sorted_params['pseudopotentials']['Fe2'] == 'Fe.pbe-spn.UPF'
    
    def test_hubbard_parameters_preserved(self):
        """Test that Hubbard parameters are preserved during extraction."""
        atoms = bulk('Fe', cubic=True)
        base_pseudos = {'Fe': 'Fe.pbe-spn.UPF'}
        
        config = setup_magnetic_config(atoms, {
            'Fe': {'mag': [1, -1], 'U': 4.3}
        }, pseudopotentials=base_pseudos)
        
        # Pass full config
        parameters = {
            'input_data': config
        }
        
        sorted_params, unused = sort_qe_input(parameters)
        
        # Pseudopotentials should be extracted
        assert 'pseudopotentials' in sorted_params
        
        # Hubbard parameters should remain in INPUT_NTYP (normalized to uppercase)
        assert 'INPUT_NTYP' in sorted_params['input_data']
        assert 'Hubbard_U' in sorted_params['input_data']['INPUT_NTYP']
    
    def test_no_pseudopotentials_in_input_data(self):
        """Test that sort_qe_input works when no pseudopotentials in input_data."""
        parameters = {
            'pseudopotentials': {
                'Fe': 'Fe.pbe-spn.UPF'
            },
            'input_data': {
                'input_ntyp': {
                    'starting_magnetization': {
                        'Fe': 1.0
                    }
                }
            }
        }
        
        sorted_params, unused = sort_qe_input(parameters)
        
        # Should work normally
        assert 'pseudopotentials' in sorted_params
        assert sorted_params['pseudopotentials']['Fe'] == 'Fe.pbe-spn.UPF'
    
    def test_empty_input_data(self):
        """Test that empty input_data doesn't cause errors."""
        parameters = {
            'pseudopotentials': {
                'Fe': 'Fe.pbe-spn.UPF'
            }
        }
        
        sorted_params, unused = sort_qe_input(parameters)
        
        # Should work and create input_data
        assert 'input_data' in sorted_params
        assert 'pseudopotentials' in sorted_params


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
