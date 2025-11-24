# Code Detection Enhancements

This document describes the enhancements made to the Quantum ESPRESSO code detection system.

## Problem Statement

The original issue (in Portuguese) requested the following improvements:

1. When detecting codes for a machine that has a `.json` file or is listed in `machines.json`, automatically detect:
   - Whether it's local or remote access via SSH
   - The SSH port (defaulting to 22)

2. Handle existing configuration files:
   - Check if a configuration already exists before overwriting
   - Provide options to overwrite or merge

3. Handle systems without the `module` command:
   - Many systems don't have environment modules
   - Need graceful fallback

4. Clarify where modules configuration belongs:
   - Should modules be per-version (in codes config)?
   - Or per-machine (in machine config)?

## Solutions Implemented

### 1. Machine Configuration Auto-Loading

The `detect_qe_codes()` function now automatically loads machine configuration:

**Before:**
```python
# Had to manually specify all SSH details
config = detect_qe_codes(
    machine_name="cluster",
    ssh_connection={'host': 'cluster.edu', 'username': 'user'},
    modules=['quantum-espresso/7.2']
)
```

**After:**
```python
# Just provide the machine name - everything is auto-loaded
config = detect_qe_codes(machine_name="cluster")
```

**Implementation:**
- Added `auto_load_machine` parameter (default: `True`)
- Loads machine config from `~/.xespresso/machines.json` or `~/.xespresso/machines/<machine_name>.json`
- Extracts SSH details (host, username, port)
- Uses configured modules automatically
- Falls back to manual parameters if machine config not found

### 2. SSH Port Support

All SSH connections now support custom ports with a default of 22:

**Changes:**
- Added `port` parameter to SSH connection dictionaries
- Default value: 22 (standard SSH port)
- Updated `_detect_remote()` to use port in SSH commands
- Updated `detect_qe_version()` to use port
- Machine class already had port support, now fully integrated

**Example:**
```python
# Explicit port
ssh_connection = {'host': 'server.com', 'username': 'user', 'port': 2222}

# Default port (22)
ssh_connection = {'host': 'server.com', 'username': 'user'}

# From machine config (port stored in machine.json)
machine = load_machine('my_cluster')
# machine.port is automatically used when detecting codes
```

### 3. Module Command Fallback

Added graceful handling for systems without the `module` command:

**Implementation:**
- New method: `CodesManager._check_module_available()`
- Checks if `module` command exists using `command -v module`
- Works for both local and remote systems
- If unavailable, skips module loading with a warning
- Detection continues with other methods (PATH search, manual paths)

**Behavior:**
```python
# If modules specified but command not available:
config = detect_qe_codes(
    machine_name="cluster",
    modules=['quantum-espresso/7.2']  # ⚠️ Skipped if module cmd not found
)
# Output: "⚠️ 'module' command not available, skipping module loads"
```

### 4. Overwrite Protection and Merge

Enhanced `CodesManager.save_config()` with multiple save modes:

**New Parameters:**
- `overwrite` (bool): If True, overwrites without asking
- `merge` (bool): If True, merges with existing config
- `interactive` (bool): If True, prompts user when file exists

**Modes:**

1. **Interactive (default):**
   ```python
   CodesManager.save_config(config)
   # Prompts: "Choose: [o]verwrite, [m]erge, [c]ancel"
   ```

2. **Overwrite:**
   ```python
   CodesManager.save_config(config, overwrite=True)
   # Overwrites existing file without asking
   ```

3. **Merge:**
   ```python
   CodesManager.save_config(config, merge=True)
   # Intelligently merges codes, versions, modules, etc.
   ```

4. **Non-interactive (for scripts/tests):**
   ```python
   CodesManager.save_config(config, interactive=False)
   # Raises FileExistsError if file exists
   ```

**Merge Logic:**
- Combines codes from both configs
- Merges version-specific configurations
- Updates metadata (qe_prefix, qe_version, modules, environment)
- Preserves existing data when not explicitly overridden

### 5. Modules: Per-Version or Per-Machine?

**Answer: Both, with clear separation of concerns**

**Per-Machine (in machine config):**
- Use for system-level modules needed for execution environment
- Examples: compiler toolchains, MPI implementations
- Loaded once when accessing the machine
- Part of the execution environment setup

**Per-Version (in codes config):**
- Use for QE version-specific modules
- Examples: `quantum-espresso/7.2`, `quantum-espresso/6.8`
- Loaded when using a specific QE version
- Part of the codes configuration

**Example Machine Config:**
```json
{
  "machines": {
    "cluster": {
      "use_modules": true,
      "modules": ["intel-compiler/2021", "intel-mpi/2021"]
    }
  }
}
```

**Example Codes Config:**
```json
{
  "machine_name": "cluster",
  "versions": {
    "7.2": {
      "modules": ["quantum-espresso/7.2"],
      "codes": {...}
    },
    "6.8": {
      "modules": ["quantum-espresso/6.8"],
      "codes": {...}
    }
  }
}
```

## API Changes

### Modified Functions

#### `detect_qe_codes()`
**New parameters:**
- `auto_load_machine` (bool, default=True): Enable/disable automatic machine config loading

**Enhanced behavior:**
- Automatically loads machine configuration
- Extracts SSH details including port
- Uses modules from machine config
- Checks module command availability

#### `CodesManager.detect_codes()`
**New parameters:**
- `use_modules` (bool, default=True): Control module usage
- Enhanced `ssh_connection` dict: Now includes `port` parameter

#### `CodesManager.save_config()`
**New parameters:**
- `overwrite` (bool, default=False): Overwrite without asking
- `merge` (bool, default=False): Merge with existing config
- `interactive` (bool, default=True): Prompt user in interactive mode

**New behavior:**
- Raises `FileExistsError` if file exists and no overwrite/merge specified (non-interactive)
- Prompts user in interactive mode
- Intelligently merges configurations

#### `create_codes_config()`
**New parameters:**
- `overwrite` (bool, default=False)
- `merge` (bool, default=True): Default changed to merge mode
- `auto_load_machine` (bool, default=True)

## Testing

### Unit Tests
Added comprehensive tests in `tests/test_codes.py`:

- `TestNewFeatures.test_save_config_overwrite_protection`
- `TestNewFeatures.test_save_config_merge`
- `TestNewFeatures.test_save_config_merge_versions`
- `TestNewFeatures.test_check_module_available_local`
- `TestNewFeatures.test_ssh_connection_with_port`
- `TestNewFeatures.test_detect_codes_with_use_modules_false`
- `TestNewFeatures.test_create_codes_config_with_merge`

### Integration Tests
Created comprehensive integration tests:
- Machine auto-loading from `machines.json`
- Machine loading from individual files
- Port default value verification
- Module detection and fallback

All tests pass successfully.

## Documentation

Updated documentation files:
- `docs/CODES_CONFIGURATION.md`: Added "New Features" section
- `examples/codes_example.py`: Added examples of new features
- Updated inline documentation and docstrings

## Backward Compatibility

All changes are backward compatible:
- Existing code continues to work without modifications
- New parameters have sensible defaults
- Auto-loading can be disabled with `auto_load_machine=False`
- Old configuration files are fully supported

## Usage Examples

### Simple Auto-Loading
```python
from xespresso.codes import detect_qe_codes

# Just works™ - loads everything from machine config
config = detect_qe_codes(machine_name="my_cluster")
```

### With Custom Port
```python
# In machine config
{
  "machines": {
    "cluster": {
      "host": "server.com",
      "username": "user",
      "port": 2222
    }
  }
}

# Detection automatically uses port 2222
config = detect_qe_codes(machine_name="cluster")
```

### Merge Configurations
```python
from xespresso.codes import CodesManager, Code, CodesConfig

# Add pw code
config1 = CodesConfig(machine_name="cluster")
config1.add_code(Code(name="pw", path="/opt/qe/bin/pw.x"))
CodesManager.save_config(config1, overwrite=True)

# Add hp code (merge with existing)
config2 = CodesConfig(machine_name="cluster")
config2.add_code(Code(name="hp", path="/opt/qe/bin/hp.x"))
CodesManager.save_config(config2, merge=True)

# Now both pw and hp are in the config
```

## Benefits

1. **Reduced Boilerplate**: No need to manually specify SSH details repeatedly
2. **Better Error Handling**: Graceful fallback for missing module command
3. **Data Safety**: Overwrite protection prevents accidental data loss
4. **Flexibility**: Merge mode allows incremental configuration building
5. **Clear Separation**: Modules are organized by purpose (machine vs. version)
6. **Port Support**: Full support for non-standard SSH ports

## Migration Guide

No migration needed! Existing code works as-is. To benefit from new features:

1. **Use auto-loading**: Simply pass machine name
2. **Add port to machine configs**: Include `"port": 22` (or custom value)
3. **Use merge mode**: Set `merge=True` when adding codes incrementally
4. **Organize modules**: Put system modules in machine config, QE modules in codes config

## Future Enhancements

Possible future improvements:
- Support for password authentication (currently only key-based)
- Parallel detection across multiple machines
- Cache detection results to avoid repeated SSH connections
- GUI for configuration management
- Auto-detection of available modules on remote systems
