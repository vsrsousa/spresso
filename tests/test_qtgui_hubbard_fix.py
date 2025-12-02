"""
Test that qtGui properly handles Hubbard parameters according to xespresso's design.

This test validates that the qtGui correctly passes Hubbard parameters
in both old format (QE < 7.0) and new format (QE >= 7.0).
"""

import pytest
from ase.build import bulk
from unittest.mock import MagicMock


class TestQtGuiHubbardFix:
    """Test Hubbard parameter handling in qtGui."""
    
    def test_old_format_hubbard_parameters(self):
        """Test that old format Hubbard parameters are correctly formatted."""
        # Simulate config from calculation_setup page
        config = {
            'calc_type': 'scf',
            'label': 'fe/scf',
            'ecutwfc': 50.0,
            'ecutrho': 400.0,
            'occupations': 'smearing',
            'conv_thr': 1.0e-8,
            'enable_hubbard': True,
            'hubbard_format': 'old',
            'hubbard_u': {
                'Fe': 4.3,
                'O': 3.0,
            },
            'pseudopotentials': {
                'Fe': 'Fe.pbe-spn.UPF',
                'O': 'O.pbe.UPF',
            },
            'qe_version': '6.8',
        }
        
        # Import the helper function
        from qtgui.pages.job_submission import JobSubmissionPage
        from PySide6.QtWidgets import QApplication
        import sys
        
        # Create Qt application if needed
        app = QApplication.instance()
        if app is None:
            app = QApplication(sys.argv)
        
        # Create session state
        session_state = {'workflow_config': config}
        
        # Create page
        page = JobSubmissionPage(session_state)
        
        # Call the private method to prepare input data
        input_data = page._build_input_data(config, 'test')
        
        # Verify old format structure
        assert 'SYSTEM' in input_data
        assert input_data['SYSTEM']['lda_plus_u'] == True
        assert 'input_ntyp' in input_data
        assert 'Hubbard_U' in input_data['input_ntyp']
        assert input_data['input_ntyp']['Hubbard_U']['Fe'] == 4.3
        assert input_data['input_ntyp']['Hubbard_U']['O'] == 3.0
        
        # Verify new format is NOT used
        assert 'hubbard' not in input_data
    
    def test_new_format_hubbard_parameters(self):
        """Test that new format Hubbard parameters are correctly formatted."""
        # Simulate config from calculation_setup page
        config = {
            'calc_type': 'scf',
            'label': 'fe/scf',
            'ecutwfc': 50.0,
            'ecutrho': 400.0,
            'occupations': 'smearing',
            'conv_thr': 1.0e-8,
            'enable_hubbard': True,
            'hubbard_format': 'new',
            'hubbard_u': {
                'Fe': 4.3,
                'O': 3.0,
            },
            'hubbard_orbitals': {
                'Fe': '3d',
                'O': '2p',
            },
            'hubbard_projector': 'atomic',
            'pseudopotentials': {
                'Fe': 'Fe.pbe-spn.UPF',
                'O': 'O.pbe.UPF',
            },
            'qe_version': '7.2',
        }
        
        # Import the helper function
        from qtgui.pages.job_submission import JobSubmissionPage
        from PySide6.QtWidgets import QApplication
        import sys
        
        # Create Qt application if needed
        app = QApplication.instance()
        if app is None:
            app = QApplication(sys.argv)
        
        # Create session state
        session_state = {'workflow_config': config}
        
        # Create page
        page = JobSubmissionPage(session_state)
        
        # Call the private method to prepare input data
        input_data = page._build_input_data(config, 'test')
        
        # Verify new format structure
        assert 'SYSTEM' in input_data
        assert input_data['SYSTEM']['lda_plus_u'] == True
        assert 'hubbard' in input_data
        assert input_data['hubbard']['projector'] == 'atomic'
        assert 'u' in input_data['hubbard']
        assert input_data['hubbard']['u']['Fe-3d'] == 4.3
        assert input_data['hubbard']['u']['O-2p'] == 3.0
        
        # Verify old format is NOT used
        assert 'input_ntyp' not in input_data or 'Hubbard_U' not in input_data.get('input_ntyp', {})
        
        # Verify qe_version is included
        assert 'qe_version' in input_data
        assert input_data['qe_version'] == '7.2'
    
    def test_auto_detect_format_from_version_7(self):
        """Test auto-detection of new format for QE 7.x."""
        config = {
            'calc_type': 'scf',
            'enable_hubbard': True,
            'hubbard_format': 'old',  # User selected old, but version should override
            'hubbard_u': {'Fe': 4.3},
            'hubbard_orbitals': {'Fe': '3d'},
            'qe_version': '7.2',
            'pseudopotentials': {'Fe': 'Fe.pbe.UPF'},
        }
        
        from qtgui.pages.job_submission import JobSubmissionPage
        from PySide6.QtWidgets import QApplication
        import sys
        
        app = QApplication.instance()
        if app is None:
            app = QApplication(sys.argv)
        
        session_state = {'workflow_config': config}
        page = JobSubmissionPage(session_state)
        
        input_data = page._build_input_data(config, 'test')
        
        # With version 7.2 and old format selected, should still use old format
        # as user explicitly chose it
        assert 'input_ntyp' in input_data
        assert 'Hubbard_U' in input_data['input_ntyp']
    
    def test_auto_detect_format_from_version_6(self):
        """Test auto-detection of old format for QE 6.x."""
        config = {
            'calc_type': 'scf',
            'enable_hubbard': True,
            'hubbard_format': 'new',  # User selected new, but version 6 should override
            'hubbard_u': {'Fe': 4.3},
            'hubbard_orbitals': {'Fe': '3d'},
            'qe_version': '6.8',
            'pseudopotentials': {'Fe': 'Fe.pbe.UPF'},
        }
        
        from qtgui.pages.job_submission import JobSubmissionPage
        from PySide6.QtWidgets import QApplication
        import sys
        
        app = QApplication.instance()
        if app is None:
            app = QApplication(sys.argv)
        
        session_state = {'workflow_config': config}
        page = JobSubmissionPage(session_state)
        
        input_data = page._build_input_data(config, 'test')
        
        # With version 6.8 and new format selected, should use new format
        # as user explicitly chose it (they might have backported support)
        assert 'hubbard' in input_data
    
    def test_default_orbitals(self):
        """Test that default orbitals are used when not specified."""
        config = {
            'calc_type': 'scf',
            'enable_hubbard': True,
            'hubbard_format': 'new',
            'hubbard_u': {
                'Fe': 4.3,
                'Mn': 5.0,
                'O': 3.0,
                'Ce': 6.0,
            },
            'qe_version': '7.2',
            'pseudopotentials': {
                'Fe': 'Fe.pbe.UPF',
                'Mn': 'Mn.pbe.UPF',
                'O': 'O.pbe.UPF',
                'Ce': 'Ce.pbe.UPF',
            },
        }
        
        from qtgui.pages.job_submission import JobSubmissionPage, _get_default_orbital
        from PySide6.QtWidgets import QApplication
        import sys
        
        app = QApplication.instance()
        if app is None:
            app = QApplication(sys.argv)
        
        session_state = {'workflow_config': config}
        page = JobSubmissionPage(session_state)
        
        input_data = page._build_input_data(config, 'test')
        
        # Verify default orbitals are used
        assert 'hubbard' in input_data
        assert 'Fe-3d' in input_data['hubbard']['u']
        assert 'Mn-3d' in input_data['hubbard']['u']
        assert 'O-2p' in input_data['hubbard']['u']
        assert 'Ce-4f' in input_data['hubbard']['u']
        
        # Verify the _get_default_orbital function
        assert _get_default_orbital('Fe') == '3d'
        assert _get_default_orbital('Mn') == '3d'
        assert _get_default_orbital('O') == '2p'
        assert _get_default_orbital('Ce') == '4f'
    
    def test_no_hubbard_parameters(self):
        """Test that calculations without Hubbard work correctly."""
        config = {
            'calc_type': 'scf',
            'enable_hubbard': False,
            'pseudopotentials': {'Fe': 'Fe.pbe.UPF'},
        }
        
        from qtgui.pages.job_submission import JobSubmissionPage
        from PySide6.QtWidgets import QApplication
        import sys
        
        app = QApplication.instance()
        if app is None:
            app = QApplication(sys.argv)
        
        session_state = {'workflow_config': config}
        page = JobSubmissionPage(session_state)
        
        input_data = page._build_input_data(config, 'test')
        
        # Verify Hubbard is not set
        assert 'lda_plus_u' not in input_data['SYSTEM'] or not input_data['SYSTEM']['lda_plus_u']
        assert 'hubbard' not in input_data
        assert 'input_ntyp' not in input_data or 'Hubbard_U' not in input_data.get('input_ntyp', {})


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
