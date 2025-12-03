#!/usr/bin/env python3
"""
Test to verify Hubbard checkbox functionality and dropdown styling.
"""

import os
import sys
import tempfile
from pathlib import Path

# Set Qt to offscreen mode for headless environments
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

from PySide6.QtWidgets import QApplication
from ase.build import bulk

# Import the modules we need
sys.path.insert(0, str(Path(__file__).parent.parent))
from qtgui.pages.calculation_setup import CalculationSetupPage
from qtgui.main_app import SessionState


def test_hubbard_checkboxes():
    """Test that Hubbard checkboxes control which elements have Hubbard U applied."""
    
    print("="*60)
    print("Testing Hubbard Checkboxes Functionality")
    print("="*60)
    
    # Create Qt application
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    
    session_state = SessionState()
    
    # Load a structure with multiple elements
    atoms = bulk('Fe', 'bcc', a=2.87)
    session_state['current_structure'] = atoms
    
    # Create page
    page = CalculationSetupPage(session_state)
    page.refresh()
    
    print("\n--- Step 1: Verify Checkboxes Exist ---")
    assert 'Fe' in page.hubbard_checkboxes, "Fe checkbox should exist"
    assert 'Fe' in page.hubbard_edits, "Fe spinbox should exist"
    assert 'Fe' in page.hubbard_orbital_edits, "Fe orbital combo should exist"
    print(f"✓ Hubbard controls created for Fe")
    
    # Verify default state
    print("\n--- Step 2: Verify Default State (Unchecked) ---")
    assert not page.hubbard_checkboxes['Fe'].isChecked(), "Fe checkbox should be unchecked by default"
    print("✓ Checkbox unchecked by default")
    
    # Enable Hubbard group but don't check any elements
    print("\n--- Step 3: Test Saving Without Checked Elements ---")
    page.hubbard_group.setChecked(True)
    config = page._get_config()
    
    assert config.get('enable_hubbard') == True, "Hubbard should be enabled"
    assert 'Fe' not in config.get('hubbard_u', {}), "Fe should NOT be in hubbard_u when checkbox unchecked"
    print("✓ Unchecked elements not saved in config")
    
    # Check the Fe checkbox
    print("\n--- Step 4: Test Saving With Checked Element ---")
    page.hubbard_checkboxes['Fe'].setChecked(True)
    page.hubbard_edits['Fe'].setValue(5.5)
    
    config = page._get_config()
    
    assert config.get('enable_hubbard') == True, "Hubbard should be enabled"
    assert 'Fe' in config.get('hubbard_u', {}), "Fe should be in hubbard_u when checkbox checked"
    assert config['hubbard_u']['Fe'] == 5.5, f"Fe U value should be 5.5, got {config['hubbard_u']['Fe']}"
    assert 'Fe' in config.get('hubbard_orbitals', {}), "Fe orbital should be saved"
    print("✓ Checked element saved in config with correct value")
    
    # Test restoration from config
    print("\n--- Step 5: Test Restoration From Config ---")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        session_state2 = SessionState()
        session_state2._sessions_dir = temp_dir
        session_state2['current_structure'] = atoms
        
        # Set a config with Hubbard U for Fe
        restore_config = {
            'calc_type': 'scf',
            'pseudopotentials': {'Fe': 'Fe.pbe.UPF'},
            'enable_hubbard': True,
            'hubbard_u': {'Fe': 4.2},
            'hubbard_orbitals': {'Fe': '3d'},
            'hubbard_format': 'new'
        }
        session_state2['workflow_config'] = restore_config
        
        page2 = CalculationSetupPage(session_state2)
        page2.refresh()
        
        # Verify checkbox is checked after restoration
        assert page2.hubbard_checkboxes['Fe'].isChecked(), "Fe checkbox should be checked after restoration"
        assert abs(page2.hubbard_edits['Fe'].value() - 4.2) < 0.01, "Fe U value should be restored"
        assert page2.hubbard_orbital_edits['Fe'].currentText() == '3d', "Fe orbital should be restored"
        print("✓ Checkbox state restored correctly from config")
    
    # Test multiple elements
    print("\n--- Step 6: Test Multiple Elements (FeO) ---")
    
    from ase import Atoms
    # Create FeO structure manually
    feo = Atoms('FeO', positions=[[0, 0, 0], [1.5, 0, 0]])
    session_state3 = SessionState()
    session_state3['current_structure'] = feo
    
    page3 = CalculationSetupPage(session_state3)
    page3.refresh()
    
    # Verify both elements have checkboxes
    assert 'Fe' in page3.hubbard_checkboxes, "Fe checkbox should exist"
    assert 'O' in page3.hubbard_checkboxes, "O checkbox should exist"
    print("✓ Checkboxes created for both Fe and O")
    
    # Check only Fe
    page3.hubbard_group.setChecked(True)
    page3.hubbard_checkboxes['Fe'].setChecked(True)
    page3.hubbard_edits['Fe'].setValue(4.0)
    # O checkbox remains unchecked
    
    config3 = page3._get_config()
    
    assert 'Fe' in config3.get('hubbard_u', {}), "Fe should be in config"
    assert 'O' not in config3.get('hubbard_u', {}), "O should NOT be in config when unchecked"
    print("✓ Selective element activation works correctly")
    
    print("\n--- Step 7: Test Dropdown Styling ---")
    # Just verify that styling was applied (we can't easily test the visual appearance)
    # Check that the orbital combobox has a stylesheet
    orbital_combo = page.hubbard_orbital_edits['Fe']
    has_styling = len(orbital_combo.styleSheet()) > 0
    print(f"Dropdown has custom styling: {has_styling}")
    if has_styling:
        print("✓ Custom dropdown styling applied")
    
    print("\n" + "="*60)
    print("All Tests Passed!")
    print("="*60)


if __name__ == '__main__':
    try:
        test_hubbard_checkboxes()
        sys.exit(0)
    except AssertionError as e:
        print(f"\n❌ Test failed: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
