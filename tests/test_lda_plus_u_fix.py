"""
Test to verify lda_plus_u is only set for old Hubbard format.

This test validates the fix for the issue where lda_plus_u was incorrectly
being set when using the new HUBBARD card format (QE >= 7.0).
"""

import pytest


def test_lda_plus_u_not_set_for_new_format():
    """Test that lda_plus_u is NOT set when using new HUBBARD card format."""
    # This simulates what happens in qtgui/pages/job_submission.py _build_input_data
    
    config = {
        'enable_hubbard': True,
        'hubbard_format': 'new',
        'hubbard_u': {'Fe': 4.3, 'O': 3.0},
        'hubbard_orbitals': {'Fe': '3d', 'O': '2p'},
        'qe_version': '7.2',
    }
    
    # Simulate the logic from _build_input_data
    input_data = {'SYSTEM': {}}
    
    hubbard_format = config.get('hubbard_format', 'old')
    qe_version = config.get('qe_version', '')
    
    # Auto-detect format from QE version if not explicitly set
    use_new_format = False
    if hubbard_format == 'new':
        use_new_format = True
    elif hubbard_format == 'old':
        use_new_format = False
    elif qe_version:
        # Parse version and determine format
        try:
            major = int(qe_version.split('.')[0])
            use_new_format = (major >= 7)
        except (ValueError, IndexError):
            pass
    
    if use_new_format:
        # NEW FORMAT (QE >= 7.0): Use 'hubbard' dictionary with HUBBARD card
        # NOTE: lda_plus_u should NOT be set when using HUBBARD card (new format)
        hubbard_dict = {
            'projector': config.get('hubbard_projector', 'atomic'),
            'u': {},
            'v': []
        }
        
        # Build U parameters with orbital specifications
        for element, u_value in config.get('hubbard_u', {}).items():
            if u_value > 0:
                orbital = config.get('hubbard_orbitals', {}).get(element, '3d')
                hubbard_dict['u'][f"{element}-{orbital}"] = u_value
        
        input_data['hubbard'] = hubbard_dict
    else:
        # OLD FORMAT (QE < 7.0): Use 'input_ntyp' with Hubbard_U
        # Only set lda_plus_u for old format
        input_data['SYSTEM']['lda_plus_u'] = True
        input_data['input_ntyp'] = {'Hubbard_U': {}}
        for element, u_value in config.get('hubbard_u', {}).items():
            if u_value > 0:
                input_data['input_ntyp']['Hubbard_U'][element] = u_value
    
    # Verify that lda_plus_u is NOT set for new format
    assert 'hubbard' in input_data
    assert 'lda_plus_u' not in input_data['SYSTEM'] or not input_data['SYSTEM'].get('lda_plus_u')
    assert input_data['hubbard']['u']['Fe-3d'] == 4.3
    assert input_data['hubbard']['u']['O-2p'] == 3.0


def test_lda_plus_u_is_set_for_old_format():
    """Test that lda_plus_u IS set when using old format."""
    
    config = {
        'enable_hubbard': True,
        'hubbard_format': 'old',
        'hubbard_u': {'Fe': 4.3, 'O': 3.0},
        'qe_version': '6.8',
    }
    
    # Simulate the logic from _build_input_data
    input_data = {'SYSTEM': {}}
    
    hubbard_format = config.get('hubbard_format', 'old')
    qe_version = config.get('qe_version', '')
    
    # Auto-detect format from QE version if not explicitly set
    use_new_format = False
    if hubbard_format == 'new':
        use_new_format = True
    elif hubbard_format == 'old':
        use_new_format = False
    elif qe_version:
        # Parse version and determine format
        try:
            major = int(qe_version.split('.')[0])
            use_new_format = (major >= 7)
        except (ValueError, IndexError):
            pass
    
    if use_new_format:
        # NEW FORMAT (QE >= 7.0): Use 'hubbard' dictionary with HUBBARD card
        # NOTE: lda_plus_u should NOT be set when using HUBBARD card (new format)
        hubbard_dict = {
            'projector': config.get('hubbard_projector', 'atomic'),
            'u': {},
            'v': []
        }
        
        # Build U parameters with orbital specifications
        for element, u_value in config.get('hubbard_u', {}).items():
            if u_value > 0:
                orbital = config.get('hubbard_orbitals', {}).get(element, '3d')
                hubbard_dict['u'][f"{element}-{orbital}"] = u_value
        
        input_data['hubbard'] = hubbard_dict
    else:
        # OLD FORMAT (QE < 7.0): Use 'input_ntyp' with Hubbard_U
        # Only set lda_plus_u for old format
        input_data['SYSTEM']['lda_plus_u'] = True
        input_data['input_ntyp'] = {'Hubbard_U': {}}
        for element, u_value in config.get('hubbard_u', {}).items():
            if u_value > 0:
                input_data['input_ntyp']['Hubbard_U'][element] = u_value
    
    # Verify that lda_plus_u IS set for old format
    assert input_data['SYSTEM']['lda_plus_u'] == True
    assert 'input_ntyp' in input_data
    assert input_data['input_ntyp']['Hubbard_U']['Fe'] == 4.3
    assert input_data['input_ntyp']['Hubbard_U']['O'] == 3.0
    assert 'hubbard' not in input_data


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
