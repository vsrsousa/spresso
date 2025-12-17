# Changes to Address Security Issue

## Security Issue Addressed

**Issue**: Command-line injection vulnerability in `CodesManager.list_available_modules()` method
- **Location**: `xespresso/codes/manager.py:304`
- **Type**: `py/command-line-injection`
- **Severity**: Medium
- **Risk**: User-provided inputs were used to construct shell commands executed with `subprocess.run(cmd, shell=True, ...)`

## Changes Made

### 1. Backend Changes

**File**: `xespresso/codes/manager.py`
- **Removed**: `list_available_modules()` static method (lines 250-345)
- This method allowed discovering available modules on remote systems but had a security vulnerability

### 2. GUI Changes

**File**: `qtgui/pages/codes_config.py`
- **Removed**: "🔍 Discover Available Modules" expandable section (Feature 1)
- **Removed**: Module discovery UI with search pattern and environment setup inputs
- **Removed**: Auto-population of discovered modules into the detection form
- **Updated**: Feature numbering (Feature 2 → Feature 1, Feature 3 → Feature 2)

### 3. Test Changes

**File**: `tests/test_gui_qe_features.py`
- **Removed**: `test_module_listing_functionality()` test
- **Updated**: `test_codes_config_module_has_updated_features()` to remove module listing assertions
- **Updated**: `test_gui_page_structure_with_new_features()` to remove module listing checks
- **Updated**: Test count from 9 tests to 8 tests

### 4. Files Deleted

- `tests/test_module_listing.py` - Tests for list_available_modules functionality
- `examples/module_listing_and_version_example.py` - Example code demonstrating module listing

## Impact on Users

### What Users Can NO LONGER Do

1. **Discover Available Modules via GUI**
   - Before: Users could click "🔎 List Available Modules" button to see available QE modules
   - After: This feature is removed from the GUI

2. **Discover Available Modules via API**
   - Before: Users could call `CodesManager.list_available_modules()` programmatically
   - After: This function no longer exists

3. **Auto-populate Module Names**
   - Before: Discovered modules were auto-populated into the "Modules to Load" field
   - After: Users must manually type module names

### What Users Can STILL Do

1. **Explicit QE Version Specification** ✅
   - Users can still specify the QE version explicitly in the GUI
   - Field: "QE Version (optional but recommended)"
   - This prevents compiler version confusion

2. **Version Selection** ✅
   - Users can still select and load different QE versions
   - Dropdown: "Choose QE Version"
   - Button: "Load QE {version} Configuration"

3. **Manual Module Entry** ✅
   - Users can still manually enter module names in "Modules to Load" field
   - Example: `quantum-espresso/7.2`

4. **Auto-detect Codes** ✅
   - Code auto-detection still works
   - All other detection functionality remains intact

### Workarounds for Users

Since module discovery is removed, users need to:

1. **Check Available Modules Manually**
   ```bash
   # SSH to the remote machine
   ssh user@cluster.edu
   
   # List available modules
   module avail espresso
   # or
   module spider quantum
   ```

2. **Keep a List of Module Names**
   - Document available module names for each machine
   - Create a reference file or notes for common modules

3. **Ask System Administrators**
   - Contact cluster admins to get list of available QE modules
   - Request documentation on available software modules

## Summary of Remaining Features

After this change, the GUI still provides:

### Feature 1: Explicit QE Version Specification ✅
- Input field to specify exact QE version
- Prevents auto-detection from picking up compiler versions
- Example: "7.2", "7.1", "6.8"

### Feature 2: Version Selection ✅
- Dropdown to select from available QE versions
- Load button to switch between versions
- Display of version-specific codes and modules

## Security Benefit

By removing the `list_available_modules()` function, we eliminate:
- Command-line injection vulnerability
- Potential for malicious code execution
- Security risks from unsanitized user inputs in shell commands

The trade-off is reduced convenience, but improved security for users.

## Recommendation for Future

If module discovery functionality is needed in the future, it should be reimplemented with:
1. **Input validation**: Sanitize all user inputs
2. **Parameterized commands**: Use list-based subprocess arguments instead of shell strings
3. **Whitelist approach**: Only allow predefined, safe commands
4. **No shell=True**: Avoid shell interpretation of commands

Example secure implementation:
```python
# Instead of:
cmd = f"ssh {host} 'module avail {pattern}'"
subprocess.run(cmd, shell=True)  # ❌ Vulnerable

# Use:
subprocess.run(['ssh', host, 'module', 'avail', pattern])  # ✅ Safer
```

## Files Modified in This Change

1. `xespresso/codes/manager.py` - Removed list_available_modules()
2. `qtgui/pages/codes_config.py` - Removed module discovery UI
3. `tests/test_gui_qe_features.py` - Updated tests
4. `tests/test_module_listing.py` - Deleted
5. `examples/module_listing_and_version_example.py` - Deleted

Total: 2 files modified, 2 files deleted
