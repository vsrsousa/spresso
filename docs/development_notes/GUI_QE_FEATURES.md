# GUI Implementation of QE Code Detection Features

This document describes the implementation of QE code detection features in the GUI.

## ⚠️ Security Update

**Note**: The module listing feature has been removed due to a security vulnerability (command-line injection). See `SECURITY_FIX_CHANGES.md` for details.

## Overview

Two features are available in the Streamlit GUI:

1. **Explicit Version Specification**: Specify QE version to avoid compiler version confusion
2. **Version Selection**: Select different QE versions for different calculations

## Feature 1: Explicit QE Version Specification

### Purpose
Auto-detection may incorrectly identify compiler versions (e.g., "2021.4" from Intel compiler) as QE versions. This feature allows explicit specification of the actual QE version.

### GUI Implementation
- **Location**: Auto-Detect Codes form
- **Input**: "QE Version (optional but recommended)" text field
  - Placeholder: "e.g., 7.2, 7.1, 6.8"
  - Help text explains the issue with compiler versions
- **Info Box**: Warns that auto-detection may pick up compiler versions

### Backend Integration
```python
codes_config = detect_qe_codes(
    machine_name=selected_machine,
    qe_prefix=qe_prefix,
    search_paths=search_paths,
    modules=modules,
    auto_load_machine=True,
    qe_version=qe_version  # Explicitly set version
)
```

### User Workflow
1. Navigate to "Codes Configuration" page
2. Select a machine
3. In the "Auto-Detect Codes" form, fill in:
   - Optional: QE Installation Prefix
   - **Recommended: QE Version** (e.g., "7.2")
   - Optional: Version Label
   - Optional: Modules to Load
4. Click "🔍 Auto-Detect Codes"
5. Review detected codes with correct version

### Benefits
- Prevents confusion with compiler versions
- Ensures accurate version tracking
- Makes configuration more reliable

## Feature 2: Version Selection

### Purpose
Allows users to select and load different QE versions from the same machine configuration, enabling different versions for different calculations.

### GUI Implementation
- **Location**: "Existing Codes Configuration" section
- **Display**: Shows available QE versions
- **Selector**: Dropdown to choose QE version
- **Action**: "Load QE {version} Configuration" button
- **Output**: 
  - Codes for selected version
  - Modules for selected version
  - Version stored in session state for calculations

### Backend Integration
```python
# List available versions
available_versions = existing_codes.list_versions()

# Load specific version
version_config = load_codes_config(
    machine_name, 
    codes_dir, 
    version=selected_version
)
```

### User Workflow
1. Navigate to "Codes Configuration" page
2. Select a machine with existing configuration
3. View "📦 Available QE versions" info box
4. Select version from "Choose QE Version" dropdown
5. Click "Load QE {version} Configuration" button
6. Review codes and modules for that version
7. Use selected version for calculations

### Benefits
- Easy switching between QE versions
- No need to reconfigure for different versions
- Visual confirmation of available versions
- Maintains connection to machine (no reconnection needed)

## Complete User Journey

### Initial Setup
1. **Configure Machine** (if not already done)
   - Set up SSH connection details
   - Configure scheduler settings

2. **Manually Check Available Modules** (via SSH)
   - SSH to the machine: `ssh user@cluster.edu`
   - Run: `module avail espresso` or `module spider quantum`
   - Note down available module names

3. **Configure QE Codes**
   - Fill in QE Version (e.g., "7.2")
   - Type module name into "Modules to Load"
   - Click "Auto-Detect Codes"
   - Review detected codes
   - Save configuration

### Using Multiple Versions
1. **Add Another Version**
   - Change QE Version to "7.1"
   - Update "Modules to Load" to match
   - Click "Auto-Detect Codes"
   - Save configuration (merges with existing)

2. **Switch Between Versions**
   - In "Existing Codes Configuration" section
   - Select version from dropdown (e.g., "7.2" or "7.1")
   - Click "Load QE {version} Configuration"
   - Use selected version for calculations

## Implementation Details

### Files Modified
- `xespresso/gui/pages/codes_config.py`: Main implementation

### New UI Elements
1. **Text Input**: `st.text_input("QE Version")`
2. **Selectbox**: `st.selectbox("Choose QE Version:")`
3. **Session State**: 
   - `selected_qe_version`: Stores selected version

### Error Handling
- Try-catch blocks around all API calls
- User-friendly error messages with `st.error()`
- Detailed tracebacks for debugging

### Backward Compatibility
- All features are optional
- Existing configurations continue to work
- Default behaviors preserved

## Testing

### Automated Tests
File: `tests/test_gui_qe_features.py`

Tests verify:
1. All three features present in code
2. Backend functions properly called
3. Parameters correctly passed
4. UI elements exist
5. Backward compatibility maintained

Run tests:
```bash
python tests/test_gui_qe_features.py
```

### Manual Testing Checklist

#### Module Listing
- [ ] Expander appears and can be expanded
- [ ] Search pattern input works
- [ ] Environment setup input works
- [ ] Button triggers module listing
- [ ] Results display correctly
- [ ] Error handling works for invalid inputs

#### QE Version Specification
- [ ] Version input field appears
- [ ] Help text explains the feature
- [ ] Info box warns about compiler versions
- [ ] Version is passed to backend correctly
- [ ] Detected codes show correct version

#### Version Selection
- [ ] Available versions display correctly
- [ ] Dropdown shows all versions
- [ ] Button loads correct version
- [ ] Codes for version display correctly
- [ ] Modules for version display correctly
- [ ] Session state updated correctly

## Example Screenshots

### Module Listing
```
🔍 Discover Available Modules
├── Search Pattern: espresso
├── Environment Setup: source /etc/profile
└── 🔎 List Available Modules

✅ Found 3 modules!
Available Modules:
- quantum-espresso/7.2
- quantum-espresso/7.1
- espresso/6.8

💡 Copy a module name and paste it in the 'Modules to Load' field below.
```

### QE Version Specification
```
Auto-Detect Codes
├── QE Installation Prefix: /opt/qe-7.2/bin
├── QE Version: 7.2                          ← NEW!
├── Version Label: production
└── Modules to Load: quantum-espresso/7.2

💡 Tip: Explicit Version
Auto-detection may pick up compiler versions.
It's recommended to specify the QE version explicitly!
```

### Version Selection
```
Existing Codes Configuration
✅ Loaded existing configuration
📦 Available QE versions: 7.2, 7.1, 6.8

Select QE Version for Calculations
Choose QE Version: [7.2 ▼]
[Load QE 7.2 Configuration]

✅ Loaded QE 7.2 configuration!
Codes for QE 7.2:
┌──────┬────────────────────────────┬─────────┐
│ Code │ Path                       │ Version │
├──────┼────────────────────────────┼─────────┤
│ pw   │ /opt/qe-7.2/bin/pw.x      │ 7.2     │
│ hp   │ /opt/qe-7.2/bin/hp.x      │ 7.2     │
│ dos  │ /opt/qe-7.2/bin/dos.x     │ 7.2     │
└──────┴────────────────────────────┴─────────┘
```

## Future Enhancements

Possible improvements:
1. Auto-fill version from selected module name
2. Show version in module list
3. Bulk version configuration
4. Version comparison view
5. Module availability check before detection
6. Cache discovered modules

## Related Documentation
- Backend implementation: `CODE_DETECTION_ENHANCEMENTS.md`
- Example scripts: `examples/module_listing_and_version_example.py`
- API documentation: `xespresso/codes/manager.py`
