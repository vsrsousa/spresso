"""
Test that magnetic and Hubbard configurations are properly handled in preparation.py.

This test validates that:
1. Magnetic config is passed correctly to setup_magnetic_config
2. Hubbard parameters with orbitals are handled correctly
3. The input_data is properly constructed for both old and new formats
"""

import pytest
from ase.build import bulk
from gui.calculations.preparation import prepare_calculation_from_gui


class TestMagneticHubbardPreparation:
    """Test magnetic and Hubbard configuration in preparation."""
    
    def test_magnetic_config_as_list(self):
        """Test that magnetic config values are handled as lists."""
        atoms = bulk('Fe', 'bcc', a=2.87)
        
        config = {
            'pseudopotentials': {'Fe': 'Fe.pbe-spn.UPF'},
            'ecutwfc': 50.0,
            'ecutrho': 400.0,
            'occupations': 'smearing',
            'smearing': 'gaussian',
            'degauss': 0.02,
            'conv_thr': 1.0e-8,
            'enable_magnetism': True,
            'magnetic_config': {
                'Fe': [2.2]  # Should work with list
            },
            'enable_hubbard': False
        }
        
        # Should not raise an error
        prepared_atoms, calc = prepare_calculation_from_gui(atoms, config, label='test/scf')
        
        # Verify calculator was created
        assert calc is not None
        # Verify prepared_atoms has the calculator attached
        assert prepared_atoms.calc == calc
    
    def test_hubbard_orbitals_dict(self):
        """Test that hubbard_orbitals dict is correctly parsed."""
        atoms = bulk('Gd', 'bcc', a=3.6)
        
        config = {
            'pseudopotentials': {'Gd': 'Gd.pbe.UPF'},
            'ecutwfc': 50.0,
            'ecutrho': 400.0,
            'occupations': 'smearing',
            'smearing': 'gaussian',
            'degauss': 0.02,
            'conv_thr': 1.0e-8,
            'enable_magnetism': True,
            'magnetic_config': {
                'Gd': [7.0]
            },
            'enable_hubbard': True,
            'hubbard_format': 'new',
            'hubbard_u': {
                'Gd': 6.0
            },
            'hubbard_orbitals': {
                'Gd': '4f'  # Should use 4f orbital, not 3d
            },
            'qe_version': '7.4.1'
        }
        
        # Should not raise an error
        prepared_atoms, calc = prepare_calculation_from_gui(atoms, config, label='test/scf')
        
        # Verify calculator was created
        assert calc is not None
        
        # Verify Hubbard parameters are in parameters
        assert hasattr(calc, 'parameters')
        params = calc.parameters
        
        # Check that input_data exists
        assert 'input_data' in params, f"Expected input_data in parameters, got: {params.keys()}"
        
        input_data = params['input_data']
        
        # Check for Hubbard - it should be in either 'hubbard' key or 'INPUT_NTYP'/'input_ntyp'
        has_hubbard = ('hubbard' in input_data or 
                      'INPUT_NTYP' in input_data or 
                      'input_ntyp' in input_data)
        assert has_hubbard, \
            f"Expected hubbard data in input_data, got keys: {input_data.keys()}"
        
        # For new format, check the hubbard card
        if 'hubbard' in input_data:
            hubbard_dict = input_data['hubbard']
            if isinstance(hubbard_dict, dict) and 'u' in hubbard_dict:
                hubbard_u = hubbard_dict['u']
                # Check that Gd with 4f orbital is present
                has_correct_orbital = any('Gd' in str(key) and '4f' in str(key) for key in hubbard_u.keys())
                assert has_correct_orbital, \
                    f"Expected Gd-4f in hubbard U params, got: {hubbard_u}"
    
    def test_hubbard_without_magnetism(self):
        """Test that Hubbard works without magnetism enabled."""
        atoms = bulk('Gd', 'bcc', a=3.6)
        
        config = {
            'pseudopotentials': {'Gd': 'Gd.pbe.UPF'},
            'ecutwfc': 50.0,
            'ecutrho': 400.0,
            'occupations': 'smearing',
            'smearing': 'gaussian',
            'degauss': 0.02,
            'conv_thr': 1.0e-8,
            'enable_magnetism': False,  # Magnetism disabled
            'enable_hubbard': True,     # But Hubbard enabled
            'hubbard_format': 'new',
            'hubbard_u': {
                'Gd': 6.0
            },
            'hubbard_orbitals': {
                'Gd': '4f'
            },
            'qe_version': '7.4.1'
        }
        
        # Should not raise an error
        prepared_atoms, calc = prepare_calculation_from_gui(atoms, config, label='test/scf')
        
        # Verify calculator was created
        assert calc is not None
        
        # Verify Hubbard parameters are in parameters
        assert hasattr(calc, 'parameters')
        params = calc.parameters
        
        # Even without magnetism, Hubbard should be present
        input_data = params.get('input_data', {})
        
        # Check for Hubbard - could be in 'hubbard' key (new format) or in SYSTEM namelist (old format)
        has_hubbard = ('hubbard' in input_data or 
                      (isinstance(input_data.get('SYSTEM'), dict) and 
                       input_data['SYSTEM'].get('lda_plus_u') == True))
        assert has_hubbard, \
            f"Expected Hubbard to be configured (either in 'hubbard' key or SYSTEM.lda_plus_u), got: {input_data.keys()}"
    
    def test_magnetic_config_scalar_to_list_conversion(self):
        """Test that scalar magnetic values are converted to lists."""
        atoms = bulk('Ni', 'fcc', a=3.52)
        
        config = {
            'pseudopotentials': {'Ni': 'Ni.pbe.UPF'},
            'ecutwfc': 50.0,
            'ecutrho': 400.0,
            'occupations': 'smearing',
            'smearing': 'gaussian',
            'degauss': 0.02,
            'conv_thr': 1.0e-8,
            'enable_magnetism': True,
            'magnetic_config': {
                'Ni': 0.6  # Scalar value, should be converted to [0.6]
            },
            'enable_hubbard': False
        }
        
        # Should not raise an error - preparation should handle scalar conversion
        prepared_atoms, calc = prepare_calculation_from_gui(atoms, config, label='test/scf')
        
        # Verify calculator was created
        assert calc is not None
        # Verify prepared_atoms has the calculator attached
        assert prepared_atoms.calc == calc


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
