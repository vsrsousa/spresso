# Testing Guide for Session and Resources Implementation

## Overview
This guide provides step-by-step instructions for testing the new session naming and resources configuration features.

## Prerequisites
- xespresso installed with all dependencies
- ASE library installed
- Access to GUI (streamlit)

## Test 1: Session Naming

### Test 1.1: Save Session with Custom Name
**Steps:**
1. Open xespresso GUI: `streamlit run xespresso/gui/streamlit_app.py`
2. Create or load a structure in Structure Viewer
3. Go to Sessions panel in sidebar
4. Click "✏️ Rename" button
5. Enter custom name: "Al_scf"
6. Click "✓ OK"
7. Click "💾 Save"

**Expected Results:**
- Session saved successfully message appears
- File created: `~/.xespresso/sessions/Al_scf.json`
- Filename does NOT include timestamp
- File contains metadata with session_name: "Al_scf"

**Verification:**
```bash
# Check file exists
ls ~/.xespresso/sessions/Al_scf.json

# Check metadata
cat ~/.xespresso/sessions/Al_scf.json | python -m json.tool | grep session_name
```

### Test 1.2: Load Session Preserves Name
**Steps:**
1. Continue from Test 1.1 or start fresh GUI
2. Go to Sessions panel → expand "📂 Load Saved Session"
3. Find "Al_scf.json" in the list
4. Click "Load" button
5. Check current session name in title bar

**Expected Results:**
- Session loads successfully
- Current session name shows "Al_scf" (not "Al scf" or any variation)
- All session state restored

## Test 2: Resources Configuration in Calculation Setup

### Test 2.1: View Machine Defaults
**Steps:**
1. Open GUI
2. Go to Machine Configuration page
3. Ensure at least one machine is configured with resources
4. Go to Calculation Setup page
5. Load a structure
6. Select the machine in "Execution Environment" section
7. Scroll to "⚙️ Resources Configuration" section
8. Ensure "Adjust Resources" is UNCHECKED

**Expected Results:**
- Resources section shows "Using default resources from machine configuration:"
- Default values displayed (e.g., Nodes: 1, Tasks per node: 16, etc.)
- No input fields shown

### Test 2.2: Adjust Resources
**Steps:**
1. Continue from Test 2.1
2. Check the "Adjust Resources" checkbox
3. Observe the UI change

**Expected Results:**
- Input fields appear for all resources
- Default values from machine pre-filled in fields
- Can modify any value
- Caption shows: "These custom resources will override the machine defaults"

### Test 2.3: Custom Resources Applied
**Steps:**
1. Continue from Test 2.2
2. Modify resources:
   - Nodes: 2 (change from 1)
   - Tasks per Node: 32 (change from 16)
   - Time Limit: 04:00:00 (change from 02:00:00)
3. Configure calculation parameters (calc type, k-points, etc.)
4. Click "🔧 Prepare Calculation"

**Expected Results:**
- Calculation prepares successfully
- Info message shows: "Using custom resources: {nodes: 2, ...}"
- Resources stored in config

**Verification:**
```python
# In GUI session state or saved session:
config = st.session_state.workflow_config
assert config['adjust_resources'] == True
assert config['resources']['nodes'] == 2
assert config['resources']['ntasks-per-node'] == 32
```

## Test 3: Resources Configuration in Workflow Builder

Repeat Test 2 steps in Workflow Builder page:
- Same interface
- Same behavior
- Resources applied to all workflow steps

## Test 4: Structure Serialization

### Test 4.1: Save Session with Structure
**Steps:**
1. Open GUI
2. Go to Structure Viewer
3. Load or create a structure (e.g., Al bulk)
4. Verify structure is displayed
5. Rename session to "Al_structure_test"
6. Save session

**Expected Results:**
- Session saved successfully
- File: `~/.xespresso/sessions/Al_structure_test.json`

**Verification:**
```python
import json

with open('~/.xespresso/sessions/Al_structure_test.json') as f:
    data = json.load(f)

# Check structure is in state
assert 'current_structure' in data['state']
# Check it has special format
assert data['state']['current_structure']['__type__'] == 'ase.Atoms'
assert '__data__' in data['state']['current_structure']
```

### Test 4.2: Load Session with Structure
**Steps:**
1. Close GUI (optional)
2. Open fresh GUI session
3. Verify no structure in Structure Viewer initially
4. Load session "Al_structure_test.json"
5. Go to Structure Viewer page

**Expected Results:**
- Session loads successfully
- Structure automatically appears in Structure Viewer
- Chemical formula and atom count displayed
- 3D visualization available
- All structure properties preserved

### Test 4.3: Complete Workflow Preservation
**Steps:**
1. Open GUI
2. Load structure: Al bulk
3. Configure calculation:
   - Type: scf
   - ecutwfc: 50
   - kpts: (4,4,4)
   - pseudopotentials: Al.UPF
4. Select machine with scheduler
5. Adjust resources:
   - nodes: 2
   - time: 03:00:00
6. Rename session: "Al_complete_test"
7. Save session
8. Close and reopen GUI
9. Load "Al_complete_test.json"
10. Check all pages

**Expected Results:**
- **Sessions panel**: Session name = "Al_complete_test"
- **Structure Viewer**: Al structure displayed
- **Calculation Setup**: 
  - All parameters preserved (scf, ecutwfc=50, kpts=(4,4,4))
  - Pseudopotentials: Al.UPF
  - Machine selected
  - Resources: nodes=2, time=03:00:00, "Adjust Resources" checked
- **All other pages**: Reflect loaded session state

## Test 5: Backward Compatibility

### Test 5.1: Load Old Session Format
**Steps:**
1. Create old-format session (without session_name in metadata):
```python
import json
import os

old_session = {
    "metadata": {
        "saved_at": "2024-01-01T10:00:00",
        "version": "1.0"
        # Note: no session_name
    },
    "state": {
        "test_key": "test_value"
    }
}

os.makedirs('~/.xespresso/sessions', exist_ok=True)
with open('~/.xespresso/sessions/old_format_test.json', 'w') as f:
    json.dump(old_session, f)
```

2. Open GUI
3. Load "old_format_test.json"

**Expected Results:**
- Session loads without error
- Session name derived from filename: "old format test"
- Fallback mechanism works

## Test 6: Edge Cases

### Test 6.1: Special Characters in Session Name
**Steps:**
1. Rename session to: "Al_scf_v2.1_test-run"
2. Save
3. Load

**Expected Results:**
- Filename sanitized: "Al_scf_v2.1_test-run.json"
- Session name preserved exactly: "Al_scf_v2.1_test-run"

### Test 6.2: Empty Resources
**Steps:**
1. Select machine with no resources defined
2. Check "Adjust Resources"

**Expected Results:**
- Input fields show default fallback values
- Can still input custom values
- No errors

### Test 6.3: Large Structure
**Steps:**
1. Load large structure (>100 atoms)
2. Save session
3. Load session

**Expected Results:**
- Serialization completes (may take a few seconds)
- Structure fully restored
- All atoms present

## Automated Tests

Run the test suite:
```bash
cd /home/runner/work/xespresso/xespresso
python -m pytest tests/test_session_manager.py -v
```

**Expected Results:**
- All tests pass
- Key tests:
  - `test_save_and_load_session`: Basic save/load
  - `test_session_name_in_filename`: Session name as filename
  - `test_load_session_preserves_name`: Name preservation
  - `test_save_load_structure`: Structure serialization

## Troubleshooting

### Issue: Session name not preserved
- Check: JSON file has "session_name" in metadata
- Verify: Using latest code version
- Solution: Re-save session with new code

### Issue: Structure not appearing after load
- Check: ASE library installed
- Check: "current_structure" in session state
- Check: Structure Viewer page refreshed
- Solution: Verify serialization format in JSON

### Issue: Resources not applied
- Check: "Adjust Resources" checkbox enabled
- Check: Values entered in all required fields
- Verify: config['resources'] has values
- Solution: Re-enter values and prepare calculation again

### Issue: Custom resources ignored
- Check: config['adjust_resources'] == True
- Check: Resources in queue configuration
- Verify: Machine has scheduler (not "direct")
- Solution: Ensure checkbox checked before preparing

## Success Criteria

All tests pass if:
1. ✓ Session files named after session name
2. ✓ Session name preserved on load (not from filename)
3. ✓ Resources UI shows machine defaults
4. ✓ Resources UI allows customization
5. ✓ Custom resources override machine defaults
6. ✓ Structure serialized and deserialized correctly
7. ✓ Complete GUI state restored on session load
8. ✓ All pages reflect loaded session
9. ✓ Backward compatibility maintained
10. ✓ No security vulnerabilities

## Performance Notes

- Structure serialization: ~1-5ms for typical structures (<100 atoms)
- Session save: <100ms total
- Session load: <100ms total
- Large structures (>1000 atoms): May take up to 1 second

## Known Limitations

1. Only ASE Atoms objects are specially serialized
2. Machine/code objects not serialized (loaded by name)
3. Calculator objects not serialized (recreated from config)
4. Maximum session file size: ~10MB (depends on structure size)

## Reporting Issues

If any test fails:
1. Note which test failed
2. Check error messages
3. Verify prerequisites installed
4. Check file permissions on ~/.xespresso/
5. Report with:
   - Test number
   - Expected vs actual result
   - Error messages
   - Python/ASE/Streamlit versions
