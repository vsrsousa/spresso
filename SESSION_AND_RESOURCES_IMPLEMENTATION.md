# Session Management and Resources Configuration Implementation

## Overview

This document describes the implementation of improved session management and resources configuration features in xespresso GUI.

## Features Implemented

### 1. Session Naming and File Management

#### Problem Addressed
Previously, session JSON files were named with timestamps, and when loading a session, the name would be derived from the filename rather than preserving the original session name.

#### Solution
- **Save with Session Name**: When saving a session, the JSON filename is now based on the session name (e.g., "Al_scf" becomes "Al_scf.json")
- **Preserve Session Name**: The session name is stored in the metadata of the JSON file
- **Load with Original Name**: When loading a session, the original session name is restored from metadata, not from the filename

#### Usage
```python
# When user names a session "Al_scf" and saves it:
# - File created: ~/.xespresso/sessions/Al_scf.json
# - Metadata includes: {"session_name": "Al_scf"}

# When loading Al_scf.json:
# - Session name restored as "Al_scf"
# - Not "Al scf" or any transformation of the filename
```

#### Implementation Details
- Modified `save_session()` to accept `session_name` parameter
- Modified `load_session()` to return tuple: `(state, session_name)`
- Updated `render_session_manager()` to use session name directly
- Added `session_name` to metadata in JSON files

### 2. Resources Configuration UI

#### Problem Addressed
Users needed a way to adjust computational resources (nodes, tasks, memory, time limits) for individual calculations and workflows, while still having sensible defaults from machine configuration.

#### Solution
- **Machine Defaults**: By default, resources are taken from the machine configuration
- **Adjust Resources Checkbox**: Users can enable this to customize resources for specific calculations
- **Override Capability**: Custom resources override machine defaults when enabled

#### Resources Configurable
- **nodes**: Number of compute nodes
- **ntasks-per-node**: Number of MPI tasks per node
- **mem**: Memory allocation (e.g., "32G", "64GB")
- **time**: Wall time limit (format: "HH:MM:SS")
- **partition**: Scheduler partition/queue name
- **account**: Optional billing account/project code

#### Locations
1. **Calculation Setup Page** (`xespresso/gui/pages/calculation_setup.py`)
   - Added resources section after machine/code selection
   - Resources passed to queue configuration when preparing calculations

2. **Workflow Builder Page** (`xespresso/gui/pages/workflow_builder.py`)
   - Same resources interface for workflows
   - Resources applied to all steps in the workflow

#### Usage Example
```python
# In GUI:
# 1. Select machine (e.g., "slurm_cluster")
# 2. Machine defaults shown: nodes=1, ntasks=16, mem=32G, time=02:00:00
# 3. Check "Adjust Resources" to customize
# 4. Change to: nodes=2, ntasks=32, time=04:00:00
# 5. Custom resources override machine defaults for this calculation
```

### 3. Enhanced Session State Management

#### Problem Addressed
When loading a saved session, the complete GUI state (including structures) needed to be restored so all tabs would reflect the saved session.

#### Solution
- **Structure Serialization**: ASE Atoms objects are now properly serialized to JSON format
- **Automatic Restoration**: When loading a session, all state including structures is restored
- **GUI Synchronization**: All GUI pages automatically reflect the loaded session state

#### Implementation Details
- Modified `get_serializable_state()` to detect and serialize ASE Atoms objects
- Atoms objects converted to JSON string format with special marker: `{"__type__": "ase.Atoms", "__data__": "..."}`
- Modified `restore_session()` to deserialize Atoms objects back to proper objects
- Excluded non-serializable objects that are recreated from config (calculators, machine objects)

#### Serializable State Includes
- Structure (`current_structure`): ASE Atoms object
- Configuration (`workflow_config`): Calculation parameters
- Working directory
- All user settings and preferences
- Pseudopotential selections
- k-point settings
- Machine and code selections (by name, not object)

#### State Excluded from Serialization
- Streamlit internal keys (prefixed with `_`)
- Widget state (recreated on page load)
- Machine/code objects (stored separately, loaded by name)
- Calculator objects (recreated from configuration)

## File Changes

### Modified Files
1. `xespresso/gui/utils/session_manager.py`
   - Enhanced `save_session()` with session name support
   - Updated `load_session()` to return session name
   - Added ASE Atoms serialization in `get_serializable_state()`
   - Added deserialization in `restore_session()`

2. `xespresso/gui/pages/calculation_setup.py`
   - Added resources configuration section
   - Checkbox for adjusting resources
   - Resource input fields with machine defaults
   - Resources passed to queue configuration

3. `xespresso/gui/pages/workflow_builder.py`
   - Added resources configuration section (same as calculation_setup)
   - Resources applied to workflow queue configuration

4. `tests/test_session_manager.py`
   - Updated tests for new `load_session()` signature
   - Added test for session name preservation
   - Added test for structure serialization/deserialization

## Usage Examples

### Example 1: Save and Load Session with Custom Name
```python
# In GUI:
# 1. Configure calculation for Al SCF
# 2. Rename session to "Al_scf"
# 3. Click "Save" - creates Al_scf.json
# 4. Later, click "Load" on Al_scf.json
# 5. Session restored with name "Al_scf" and all settings
```

### Example 2: Customize Resources
```python
# In GUI - Calculation Setup:
# 1. Select machine "slurm_cluster"
# 2. See default resources from machine config
# 3. Check "Adjust Resources"
# 4. Set custom values:
#    - nodes: 2 (default was 1)
#    - ntasks-per-node: 32 (default was 16)
#    - time: 04:00:00 (default was 02:00:00)
# 5. These custom resources override machine defaults
```

### Example 3: Complete Workflow
```python
# In GUI:
# 1. Load structure: Al bulk
# 2. Configure calculation:
#    - Type: SCF
#    - Machine: slurm_cluster
#    - Adjust Resources: nodes=2, time=04:00:00
# 3. Name session: "Al_scf"
# 4. Save session
# 5. Close application
# 6. Reopen application
# 7. Load "Al_scf.json"
# 8. All state restored:
#    - Structure visible in structure viewer
#    - Configuration in calculation setup
#    - Resources customization preserved
#    - Ready to submit or modify
```

## Technical Details

### ASE Atoms Serialization Format
```json
{
  "current_structure": {
    "__type__": "ase.Atoms",
    "__data__": "{\"1\": {\"cell\": {...}, \"positions\": {...}, ...}}"
  }
}
```

### Resources in Queue Configuration
```python
config["queue"] = {
    "execution": "remote",
    "scheduler": "slurm",
    "resources": {
        "nodes": 2,
        "ntasks-per-node": 32,
        "mem": "64G",
        "time": "04:00:00",
        "partition": "compute"
    },
    # ... other queue settings
}
```

## Testing

### Manual Testing Checklist
- [ ] Save session with custom name (e.g., "Al_scf")
- [ ] Verify JSON file has correct name (Al_scf.json)
- [ ] Load session and verify name is preserved
- [ ] Load session and verify structure appears in structure viewer
- [ ] Configure resources with "Adjust Resources" checkbox
- [ ] Verify machine defaults shown when checkbox disabled
- [ ] Verify custom resources override defaults when enabled
- [ ] Save session with structure and custom resources
- [ ] Load session and verify all settings restored

### Automated Tests
Run: `python -m pytest tests/test_session_manager.py -v`

Tests include:
- Session save/load with session name
- Session name preservation
- Structure serialization/deserialization
- Metadata preservation

## Future Enhancements

Possible improvements:
1. Add versioning to session format for backward compatibility
2. Add session comparison/diff functionality
3. Add session export/import between users
4. Add session templates for common workflows
5. Add resource estimation based on structure size
6. Add resource usage history and recommendations

## Compatibility

- **Python**: 3.8+
- **ASE**: Any version with JSON I/O support
- **Streamlit**: 1.0+
- **Backward Compatibility**: Old session files without session_name metadata will fall back to using filename as name
