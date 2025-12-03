#!/usr/bin/env python3
"""
Manual test script to verify magnetic and Hubbard configuration save/restore behavior.

This script:
1. Creates a session with magnetic/Hubbard configurations
2. Saves it to disk
3. Loads it back
4. Modifies the configurations
5. Saves again
6. Verifies all changes were properly saved

Usage:
    python3 manual_test_magnetic_hubbard.py
"""

import os
import sys
import json
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


def print_header(text):
    """Print a formatted header."""
    print(f"\n{'='*60}")
    print(f"  {text}")
    print(f"{'='*60}")


def print_config_status(config, label):
    """Print the status of a configuration."""
    print(f"\n{label}:")
    print(f"  Magnetism enabled: {config.get('enable_magnetism', False)}")
    if config.get('magnetic_config'):
        print(f"  Magnetic values: {config['magnetic_config']}")
    print(f"  Hubbard enabled: {config.get('enable_hubbard', False)}")
    if config.get('hubbard_u'):
        print(f"  Hubbard U values: {config['hubbard_u']}")
    if config.get('hubbard_orbitals'):
        print(f"  Hubbard orbitals: {config['hubbard_orbitals']}")
    print(f"  Has pseudopotentials: {'pseudopotentials' in config}")


def test_magnetic_hubbard_session_workflow():
    """Test the complete workflow of creating, saving, loading, and modifying a session."""
    
    print_header("Starting Magnetic/Hubbard Session Test")
    
    # Create Qt application
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    
    # Create a temporary directory for sessions
    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"\nUsing temporary session directory: {temp_dir}")
        
        # Step 1: Create initial session with configurations
        print_header("Step 1: Create Session with Initial Config")
        
        session_state = SessionState()
        session_state._sessions_dir = temp_dir
        
        # Load a structure
        atoms = bulk('Fe', 'bcc', a=2.87)
        session_state['current_structure'] = atoms
        print("✓ Loaded Fe structure")
        
        # Create initial configuration
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
            'enable_hubbard': True,
            'hubbard_format': 'new',
            'hubbard_u': {'Fe': 4.0},
            'hubbard_orbitals': {'Fe': '3d'}
        }
        session_state['workflow_config'] = initial_config
        print_config_status(initial_config, "Initial configuration")
        
        # Create page and refresh to load config into UI
        page = CalculationSetupPage(session_state)
        page.refresh()
        print("✓ Created calculation setup page and loaded config")
        
        # Verify UI state
        assert page.magnetic_group.isChecked() == True, "Magnetic checkbox should be checked"
        assert page.hubbard_group.isChecked() == True, "Hubbard checkbox should be checked"
        print("✓ Checkboxes are correctly set")
        
        # Step 2: Save session to disk
        print_header("Step 2: Save Session to Disk")
        
        session_state['session_name'] = 'test_session'
        page.save_state()
        session_state.save_session()
        
        session_file = os.path.join(temp_dir, 'test_session.json')
        assert os.path.exists(session_file), "Session file should exist"
        print(f"✓ Session saved to {session_file}")
        
        # Read and verify saved data
        with open(session_file, 'r') as f:
            saved_data = json.load(f)
        print_config_status(saved_data.get('workflow_config', {}), "Saved configuration")
        
        # Step 3: Load session from disk (simulating app restart)
        print_header("Step 3: Load Session from Disk")
        
        # Create new session state and page (simulating fresh start)
        session_state2 = SessionState()
        session_state2._sessions_dir = temp_dir
        
        # Load the session
        success = session_state2.load_session_from_file(session_file)
        assert success, "Session should load successfully"
        print("✓ Session loaded successfully")
        
        # Restore structure
        session_state2['current_structure'] = atoms
        
        # Create new page and refresh
        page2 = CalculationSetupPage(session_state2)
        page2.refresh()
        print("✓ Created new page and refreshed")
        
        # Verify UI state was restored
        assert page2.magnetic_group.isChecked() == True, "Magnetic checkbox should be checked after load"
        assert page2.hubbard_group.isChecked() == True, "Hubbard checkbox should be checked after load"
        print("✓ Checkboxes correctly restored")
        
        # Verify values were restored
        assert 'Fe' in page2.magnetic_edits
        assert abs(page2.magnetic_edits['Fe'].value() - 2.2) < 0.01, "Magnetic value should be restored"
        assert 'Fe' in page2.hubbard_edits
        assert abs(page2.hubbard_edits['Fe'].value() - 4.0) < 0.01, "Hubbard U value should be restored"
        print("✓ Values correctly restored")
        
        # Step 4: Modify configurations
        print_header("Step 4: Modify Configurations")
        
        # User modifies magnetic value
        page2.magnetic_edits['Fe'].setValue(3.5)
        print("  Changed magnetic value: 2.2 → 3.5")
        
        # User modifies Hubbard U value
        page2.hubbard_edits['Fe'].setValue(5.5)
        print("  Changed Hubbard U value: 4.0 → 5.5")
        
        # Step 5: Save modified session
        print_header("Step 5: Save Modified Session")
        
        page2.save_state()
        session_state2.save_session()
        print("✓ Modified session saved")
        
        # Read and verify modified data
        with open(session_file, 'r') as f:
            modified_data = json.load(f)
        
        modified_config = modified_data.get('workflow_config', {})
        print_config_status(modified_config, "Modified configuration")
        
        # Verify changes were saved
        assert modified_config.get('enable_magnetism') == True, "Magnetism should still be enabled"
        assert 'Fe' in modified_config.get('magnetic_config', {}), "Magnetic config should have Fe"
        assert modified_config['magnetic_config']['Fe'] == [3.5], f"Magnetic value should be 3.5, got {modified_config['magnetic_config']['Fe']}"
        print("✓ Magnetic changes saved correctly")
        
        assert modified_config.get('enable_hubbard') == True, "Hubbard should still be enabled"
        assert 'Fe' in modified_config.get('hubbard_u', {}), "Hubbard U should have Fe"
        assert modified_config['hubbard_u']['Fe'] == 5.5, f"Hubbard U should be 5.5, got {modified_config['hubbard_u']['Fe']}"
        print("✓ Hubbard changes saved correctly")
        
        # Most importantly: verify pseudopotentials were preserved!
        assert 'pseudopotentials' in modified_config, "Pseudopotentials should be preserved"
        assert 'Fe' in modified_config['pseudopotentials'], "Fe pseudopotential should be preserved"
        assert modified_config['pseudopotentials']['Fe'] == 'Fe.pbe-spn.UPF', "Pseudopotential value should be preserved"
        print("✓ Pseudopotentials preserved correctly")
        
        # Step 6: Test disabling configurations
        print_header("Step 6: Test Disabling Configurations")
        
        # Uncheck magnetic
        page2.magnetic_group.setChecked(False)
        print("  Unchecked magnetic configuration")
        
        # Save again
        page2.save_state()
        session_state2.save_session()
        
        # Read and verify
        with open(session_file, 'r') as f:
            final_data = json.load(f)
        
        final_config = final_data.get('workflow_config', {})
        print_config_status(final_config, "Final configuration")
        
        assert final_config.get('enable_magnetism') == False, "Magnetism should be disabled"
        assert final_config.get('enable_hubbard') == True, "Hubbard should still be enabled"
        print("✓ Checkbox state changes saved correctly")
        
        # Load one more time to verify checkbox restoration
        session_state3 = SessionState()
        session_state3._sessions_dir = temp_dir
        session_state3.load_session_from_file(session_file)
        session_state3['current_structure'] = atoms
        
        page3 = CalculationSetupPage(session_state3)
        page3.refresh()
        
        assert page3.magnetic_group.isChecked() == False, "Magnetic checkbox should be unchecked"
        assert page3.hubbard_group.isChecked() == True, "Hubbard checkbox should still be checked"
        print("✓ Checkbox states restored correctly")
        
        print_header("All Tests Passed Successfully!")
        print("\nSummary:")
        print("  ✓ Session creation and initial configuration")
        print("  ✓ Session save to disk")
        print("  ✓ Session load from disk")
        print("  ✓ UI state restoration (checkboxes and values)")
        print("  ✓ Configuration modification")
        print("  ✓ Modified session save (preserving pseudopotentials)")
        print("  ✓ Checkbox state changes and restoration")
        print("\n" + "="*60)


if __name__ == '__main__':
    try:
        test_magnetic_hubbard_session_workflow()
        sys.exit(0)
    except AssertionError as e:
        print(f"\n❌ Test failed: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
