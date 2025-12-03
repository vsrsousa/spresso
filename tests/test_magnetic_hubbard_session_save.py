"""
Test that magnetic and Hubbard configurations are properly saved and restored in sessions.

This test validates the fixes for:
1. Saving magnetic/Hubbard configurations when user edits them after loading a session
2. Properly restoring checkbox states when loading a session with these configurations
"""

import pytest
import sys
import json
import tempfile
import os
from pathlib import Path
from PySide6.QtWidgets import QApplication


class TestMagneticHubbardSessionSave:
    """Test magnetic and Hubbard configuration session save/restore."""
    
    @pytest.fixture(scope="class")
    def qt_app(self):
        """Create Qt application for tests."""
        app = QApplication.instance()
        if app is None:
            app = QApplication(sys.argv)
        yield app
    
    def test_save_magnetic_config_after_session_load(self, qt_app):
        """Test that magnetic config can be saved after loading a session."""
        from qtgui.pages.calculation_setup import CalculationSetupPage
        from qtgui.main_app import SessionState
        from ase.build import bulk
        
        # Create session state
        session_state = SessionState()
        
        # Load a structure to enable magnetic inputs
        atoms = bulk('Fe', 'bcc', a=2.87)
        session_state['current_structure'] = atoms
        
        # Simulate a loaded session with a prepared calculation (has pseudopotentials)
        initial_config = {
            'calc_type': 'scf',
            'label': 'Fe/scf',
            'ecutwfc': 50.0,
            'ecutrho': 400.0,
            'occupations': 'smearing',
            'conv_thr': 1.0e-8,
            'pseudopotentials': {'Fe': 'Fe.pbe-spn.UPF'},
            'enable_magnetism': True,
            'magnetic_config': {'Fe': [2.2]},
            'enable_hubbard': False,
            'hubbard_u': {}
        }
        session_state['workflow_config'] = initial_config
        
        # Create calculation setup page
        page = CalculationSetupPage(session_state)
        
        # Trigger refresh to restore UI state from session
        page.refresh()
        
        # Verify magnetic checkbox is checked
        assert page.magnetic_group.isChecked() == True, "Magnetic group should be checked after loading session"
        
        # Verify magnetic value is restored
        assert 'Fe' in page.magnetic_edits
        assert abs(page.magnetic_edits['Fe'].value() - 2.2) < 0.01
        
        # User modifies the magnetic value
        page.magnetic_edits['Fe'].setValue(3.5)
        
        # Save state (as happens when user clicks Save Session)
        page.save_state()
        
        # Verify the config was properly saved with the new magnetic value
        saved_config = session_state.get('workflow_config')
        assert saved_config is not None
        assert saved_config.get('enable_magnetism') == True
        assert 'Fe' in saved_config.get('magnetic_config', {})
        # The value is stored as a list
        assert saved_config['magnetic_config']['Fe'] == [3.5]
        # Verify pseudopotentials were preserved
        assert 'Fe' in saved_config.get('pseudopotentials', {})
        assert saved_config['pseudopotentials']['Fe'] == 'Fe.pbe-spn.UPF'
    
    def test_save_hubbard_config_after_session_load(self, qt_app):
        """Test that Hubbard config can be saved after loading a session."""
        from qtgui.pages.calculation_setup import CalculationSetupPage
        from qtgui.main_app import SessionState
        from ase.build import bulk
        
        # Create session state
        session_state = SessionState()
        
        # Load a structure to enable Hubbard inputs
        atoms = bulk('Fe', 'bcc', a=2.87)
        session_state['current_structure'] = atoms
        
        # Simulate a loaded session with a prepared calculation (has pseudopotentials)
        initial_config = {
            'calc_type': 'scf',
            'label': 'Fe/scf',
            'ecutwfc': 50.0,
            'ecutrho': 400.0,
            'occupations': 'smearing',
            'conv_thr': 1.0e-8,
            'pseudopotentials': {'Fe': 'Fe.pbe-spn.UPF'},
            'enable_magnetism': False,
            'magnetic_config': {},
            'enable_hubbard': True,
            'hubbard_format': 'new',
            'hubbard_u': {'Fe': 4.0},
            'hubbard_orbitals': {'Fe': '3d'}
        }
        session_state['workflow_config'] = initial_config
        
        # Create calculation setup page
        page = CalculationSetupPage(session_state)
        
        # Trigger refresh to restore UI state from session
        page.refresh()
        
        # Verify Hubbard checkbox is checked
        assert page.hubbard_group.isChecked() == True, "Hubbard group should be checked after loading session"
        
        # Verify Hubbard value is restored
        assert 'Fe' in page.hubbard_edits
        assert abs(page.hubbard_edits['Fe'].value() - 4.0) < 0.01
        
        # Verify orbital is restored
        assert 'Fe' in page.hubbard_orbital_edits
        assert page.hubbard_orbital_edits['Fe'].currentText() == '3d'
        
        # User modifies the Hubbard U value
        page.hubbard_edits['Fe'].setValue(5.5)
        
        # User modifies the orbital
        idx = page.hubbard_orbital_edits['Fe'].findText('4s')
        if idx >= 0:
            page.hubbard_orbital_edits['Fe'].setCurrentIndex(idx)
        
        # Save state (as happens when user clicks Save Session)
        page.save_state()
        
        # Verify the config was properly saved with the new Hubbard value
        saved_config = session_state.get('workflow_config')
        assert saved_config is not None
        assert saved_config.get('enable_hubbard') == True
        assert 'Fe' in saved_config.get('hubbard_u', {})
        assert saved_config['hubbard_u']['Fe'] == 5.5
        # Verify orbital was saved
        if idx >= 0:
            assert 'Fe' in saved_config.get('hubbard_orbitals', {})
            assert saved_config['hubbard_orbitals']['Fe'] == '4s'
        # Verify pseudopotentials were preserved
        assert 'Fe' in saved_config.get('pseudopotentials', {})
        assert saved_config['pseudopotentials']['Fe'] == 'Fe.pbe-spn.UPF'
    
    def test_checkbox_state_restoration(self, qt_app):
        """Test that checkbox states are properly restored when loading a session."""
        from qtgui.pages.calculation_setup import CalculationSetupPage
        from qtgui.main_app import SessionState
        from ase.build import bulk
        
        # Create session state
        session_state = SessionState()
        
        # Load a structure
        atoms = bulk('Fe', 'bcc', a=2.87)
        session_state['current_structure'] = atoms
        
        # Test 1: Both magnetic and Hubbard enabled
        config1 = {
            'calc_type': 'scf',
            'pseudopotentials': {'Fe': 'Fe.pbe-spn.UPF'},
            'enable_magnetism': True,
            'magnetic_config': {'Fe': [2.2]},
            'enable_hubbard': True,
            'hubbard_u': {'Fe': 4.0}
        }
        session_state['workflow_config'] = config1
        
        page = CalculationSetupPage(session_state)
        page.refresh()
        
        assert page.magnetic_group.isChecked() == True
        assert page.hubbard_group.isChecked() == True
        
        # Test 2: Only magnetic enabled
        config2 = {
            'calc_type': 'scf',
            'pseudopotentials': {'Fe': 'Fe.pbe-spn.UPF'},
            'enable_magnetism': True,
            'magnetic_config': {'Fe': [2.2]},
            'enable_hubbard': False,
            'hubbard_u': {}
        }
        session_state['workflow_config'] = config2
        
        page2 = CalculationSetupPage(session_state)
        page2.refresh()
        
        assert page2.magnetic_group.isChecked() == True
        assert page2.hubbard_group.isChecked() == False
        
        # Test 3: Only Hubbard enabled
        config3 = {
            'calc_type': 'scf',
            'pseudopotentials': {'Fe': 'Fe.pbe-spn.UPF'},
            'enable_magnetism': False,
            'magnetic_config': {},
            'enable_hubbard': True,
            'hubbard_u': {'Fe': 4.0}
        }
        session_state['workflow_config'] = config3
        
        page3 = CalculationSetupPage(session_state)
        page3.refresh()
        
        assert page3.magnetic_group.isChecked() == False
        assert page3.hubbard_group.isChecked() == True
        
        # Test 4: Both disabled
        config4 = {
            'calc_type': 'scf',
            'pseudopotentials': {'Fe': 'Fe.pbe-spn.UPF'},
            'enable_magnetism': False,
            'magnetic_config': {},
            'enable_hubbard': False,
            'hubbard_u': {}
        }
        session_state['workflow_config'] = config4
        
        page4 = CalculationSetupPage(session_state)
        page4.refresh()
        
        assert page4.magnetic_group.isChecked() == False
        assert page4.hubbard_group.isChecked() == False
    
    def test_merge_configs_preserves_pseudopotentials(self, qt_app):
        """Test that _merge_configs properly preserves pseudopotentials."""
        from qtgui.pages.calculation_setup import CalculationSetupPage
        from qtgui.main_app import SessionState
        from ase.build import bulk
        
        session_state = SessionState()
        atoms = bulk('Fe', 'bcc', a=2.87)
        session_state['current_structure'] = atoms
        
        page = CalculationSetupPage(session_state)
        
        # Existing config with pseudopotentials
        existing_config = {
            'calc_type': 'scf',
            'pseudopotentials': {'Fe': 'Fe.pbe-spn.UPF'},
            'enable_magnetism': True,
            'magnetic_config': {'Fe': [2.2]}
        }
        
        # New config without pseudopotentials (as would happen when reading from UI)
        new_config = {
            'calc_type': 'scf',
            'enable_magnetism': True,
            'magnetic_config': {'Fe': [3.5]}
        }
        
        # Merge configs
        merged = page._merge_configs(new_config, existing_config)
        
        # Verify merge happened correctly
        assert merged is not None
        assert 'pseudopotentials' in merged
        assert merged['pseudopotentials']['Fe'] == 'Fe.pbe-spn.UPF'
        assert merged['magnetic_config']['Fe'] == [3.5]
        
        # Verify original configs were not modified
        assert 'pseudopotentials' not in new_config
        assert existing_config['magnetic_config']['Fe'] == [2.2]


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
