#!/usr/bin/env python3
"""
Test to verify pseudopotentials and smearing type are saved/restored correctly.
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


def test_pseudopotentials_and_smearing_save():
    """Test that pseudopotentials and smearing type are saved and restored correctly."""
    
    print("="*60)
    print("Testing Pseudopotentials and Smearing Type Save/Restore")
    print("="*60)
    
    # Create Qt application
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    
    # Create a temporary directory for sessions
    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"\nUsing temporary session directory: {temp_dir}")
        
        # Step 1: Create session with configuration
        print("\n--- Step 1: Create Initial Configuration ---")
        
        session_state = SessionState()
        session_state._sessions_dir = temp_dir
        
        # Load a structure
        atoms = bulk('Fe', 'bcc', a=2.87)
        session_state['current_structure'] = atoms
        
        # Create initial configuration with pseudopotentials and smearing
        initial_config = {
            'calc_type': 'scf',
            'label': 'Fe/scf',
            'ecutwfc': 50.0,
            'ecutrho': 400.0,
            'occupations': 'smearing',
            'smearing': 'marzari-vanderbilt',  # Specific smearing type
            'degauss': 0.02,
            'conv_thr': 1.0e-8,
            'pseudopotentials': {'Fe': 'Fe.pbe-spn-kjpaw_psl.1.0.0.UPF'},
            'enable_magnetism': True,
            'magnetic_config': {'Fe': [2.2]},
        }
        session_state['workflow_config'] = initial_config
        
        print(f"Initial config created:")
        print(f"  Pseudopotentials: {initial_config['pseudopotentials']}")
        print(f"  Smearing type: {initial_config['smearing']}")
        
        # Create page and refresh
        page = CalculationSetupPage(session_state)
        page.refresh()
        
        # Step 2: Save session
        print("\n--- Step 2: Save Session ---")
        
        session_state['session_name'] = 'test_pseudo_smearing'
        page.save_state()
        session_state.save_session()
        
        session_file = os.path.join(temp_dir, 'test_pseudo_smearing.json')
        assert os.path.exists(session_file), "Session file should exist"
        
        # Read and verify saved data
        with open(session_file, 'r') as f:
            saved_data = json.load(f)
        
        saved_config = saved_data.get('workflow_config', {})
        print(f"Saved config:")
        print(f"  Pseudopotentials: {saved_config.get('pseudopotentials', 'MISSING!')}")
        print(f"  Smearing type: {saved_config.get('smearing', 'MISSING!')}")
        
        # Verify pseudopotentials were saved
        assert 'pseudopotentials' in saved_config, "Pseudopotentials should be in saved config"
        assert 'Fe' in saved_config['pseudopotentials'], "Fe pseudopotential should be saved"
        assert saved_config['pseudopotentials']['Fe'] == 'Fe.pbe-spn-kjpaw_psl.1.0.0.UPF', \
            f"Fe pseudopotential should be 'Fe.pbe-spn-kjpaw_psl.1.0.0.UPF', got {saved_config['pseudopotentials']['Fe']}"
        
        # Verify smearing type was saved
        assert 'smearing' in saved_config, "Smearing type should be in saved config"
        assert saved_config['smearing'] == 'marzari-vanderbilt', \
            f"Smearing type should be 'marzari-vanderbilt', got {saved_config['smearing']}"
        
        print("✓ Pseudopotentials and smearing type saved correctly")
        
        # Step 3: Load session (simulating app restart)
        print("\n--- Step 3: Load Session ---")
        
        session_state2 = SessionState()
        session_state2._sessions_dir = temp_dir
        
        success = session_state2.load_session_from_file(session_file)
        assert success, "Session should load successfully"
        
        # Restore structure
        session_state2['current_structure'] = atoms
        
        # Create new page and refresh
        page2 = CalculationSetupPage(session_state2)
        page2.refresh()
        
        print("✓ Session loaded")
        
        # Step 4: Verify UI state was restored
        print("\n--- Step 4: Verify UI Restoration ---")
        
        # Check smearing combo box
        current_smearing = page2.smearing_combo.currentText()
        print(f"  UI smearing type: {current_smearing}")
        assert current_smearing == 'marzari-vanderbilt', \
            f"Smearing combo should show 'marzari-vanderbilt', got '{current_smearing}'"
        
        # Check pseudopotentials in UI
        if hasattr(page2, 'pseudo_selector'):
            restored_pseudos = page2.pseudo_selector.get_pseudopotentials()
        else:
            restored_pseudos = {}
            for element, edit in getattr(page2, 'pseudo_edits', {}).items():
                pseudo = edit.text().strip()
                if pseudo:
                    restored_pseudos[element] = pseudo
        
        print(f"  UI pseudopotentials: {restored_pseudos}")
        assert 'Fe' in restored_pseudos, "Fe pseudopotential should be in UI"
        assert restored_pseudos['Fe'] == 'Fe.pbe-spn-kjpaw_psl.1.0.0.UPF', \
            f"Fe pseudopotential in UI should be 'Fe.pbe-spn-kjpaw_psl.1.0.0.UPF', got {restored_pseudos['Fe']}"
        
        print("✓ UI state restored correctly")
        
        # Step 5: Save again and verify persistence
        print("\n--- Step 5: Save Again and Verify Persistence ---")
        
        page2.save_state()
        session_state2.save_session()
        
        with open(session_file, 'r') as f:
            final_data = json.load(f)
        
        final_config = final_data.get('workflow_config', {})
        print(f"Final saved config:")
        print(f"  Pseudopotentials: {final_config.get('pseudopotentials', 'MISSING!')}")
        print(f"  Smearing type: {final_config.get('smearing', 'MISSING!')}")
        
        # Verify everything is still there
        assert 'pseudopotentials' in final_config, "Pseudopotentials should still be in config"
        assert 'Fe' in final_config['pseudopotentials'], "Fe pseudopotential should still be saved"
        assert final_config['pseudopotentials']['Fe'] == 'Fe.pbe-spn-kjpaw_psl.1.0.0.UPF'
        assert 'smearing' in final_config, "Smearing type should still be in config"
        assert final_config['smearing'] == 'marzari-vanderbilt'
        
        print("✓ Pseudopotentials and smearing type persisted correctly")
        
        print("\n" + "="*60)
        print("All tests passed!")
        print("="*60)


if __name__ == '__main__':
    try:
        test_pseudopotentials_and_smearing_save()
        sys.exit(0)
    except AssertionError as e:
        print(f"\n❌ Test failed: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
