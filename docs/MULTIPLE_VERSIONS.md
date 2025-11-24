# Multiple QE Versions Support

## Overview

XEspresso now supports managing multiple Quantum ESPRESSO versions on the same machine. This feature addresses a common scenario in HPC environments where multiple QE versions are installed system-wide, and users need to choose which version to use for different calculations.

## Problem Solved

The implementation addresses three key concerns from the original issue:

1. **Multiple QE versions on same machine** - How to configure and manage them
2. **Version selection per calculation** - How users can choose which version to use
3. **Connection persistence** - Ensuring SSH connections don't break when switching versions

## Key Features

### 1. Version-Aware Code Configuration

```python
from xespresso.codes import CodesConfig, Code

config = CodesConfig(machine_name="cluster1", qe_version="7.2")  # default

# Add QE 7.2
config.add_code(Code(name="pw", path="/opt/qe-7.2/bin/pw.x"), version="7.2")

# Add QE 6.8 to same machine
config.add_code(Code(name="pw", path="/opt/qe-6.8/bin/pw.x"), version="6.8")
```

### 2. Easy Version Selection

```python
from xespresso.codes import load_codes_config

codes = load_codes_config("cluster1")

# Use QE 7.2
pw_72 = codes.get_code("pw", version="7.2")

# Use QE 6.8
pw_68 = codes.get_code("pw", version="6.8")
```

### 3. Programmatic Version Management

```python
from xespresso.codes import add_version_to_config

# Add a new version to existing config
config = add_version_to_config(
    machine_name="cluster1",
    version="7.3",
    qe_prefix="/opt/qe-7.3/bin",
    modules=["quantum-espresso/7.3"]
)
```

### 4. Automatic Connection Persistence

**Most Important:** SSH connections are automatically maintained when switching between versions!

```python
# Load once
machine = load_machine("cluster1")
codes = load_codes_config("cluster1")

# First calculation with QE 7.2
pw_72 = codes.get_code("pw", version="7.2")
# ... run calculation ... (establishes SSH connection)

# Second calculation with QE 6.8
pw_68 = codes.get_code("pw", version="6.8")
# ... run calculation ... (REUSES connection!)

# No reconnection needed - connection persists!
```

**Why it works:** Connections are identified by `(host, username)`, not by QE version. The RemoteExecutionMixin caches connections and reuses them automatically.

## Configuration File Format

### Multi-Version Format

```json
{
  "machine_name": "cluster1",
  "qe_version": "7.2",
  "versions": {
    "7.2": {
      "qe_prefix": "/opt/qe-7.2/bin",
      "modules": ["quantum-espresso/7.2"],
      "codes": {
        "pw": {"name": "pw", "path": "/opt/qe-7.2/bin/pw.x", "version": "7.2"},
        "hp": {"name": "hp", "path": "/opt/qe-7.2/bin/hp.x", "version": "7.2"}
      }
    },
    "6.8": {
      "qe_prefix": "/opt/qe-6.8/bin",
      "modules": ["quantum-espresso/6.8"],
      "codes": {
        "pw": {"name": "pw", "path": "/opt/qe-6.8/bin/pw.x", "version": "6.8"},
        "hp": {"name": "hp", "path": "/opt/qe-6.8/bin/hp.x", "version": "6.8"}
      }
    }
  }
}
```

## Use Cases

### 1. Feature Requirements

```python
# Use QE 7.2 for new Hubbard parameters format
pw_72 = codes.get_code("pw", version="7.2")
hp_72 = codes.get_code("hp", version="7.2")

# Use QE 6.8 for legacy compatibility
pw_68 = codes.get_code("pw", version="6.8")
```

### 2. Testing and Validation

```python
# Compare results between versions
for version in ["7.2", "6.8"]:
    pw = codes.get_code("pw", version=version)
    # ... run same calculation with different versions ...
```

### 3. Mixed Workflows

```python
# SCF with newer version
pw_72 = codes.get_code("pw", version="7.2")
# ... run SCF ...

# Post-processing with stable version
dos_68 = codes.get_code("dos", version="6.8")
bands_68 = codes.get_code("bands", version="6.8")
# ... run DOS and bands ...
```

## Benefits

✅ **Flexibility** - Use different versions for different tasks
✅ **Performance** - Connection persistence eliminates SSH overhead
✅ **Reliability** - No connection breaks when switching versions
✅ **Organization** - Clear version specification in code
✅ **Backward Compatible** - Old single-version configs still work

## Examples

- `examples/multi_version_example.py` - Setup and basic usage
- `examples/workflow_multi_version.py` - Practical workflow example

## Documentation

- `docs/CODES_CONFIGURATION.md` - Complete API reference
- `docs/REMOTE_CONNECTION_PERSISTENCE.md` - Connection details

## Testing

All features are tested in `tests/test_codes.py`:
- 21 tests pass, including 6 new multi-version tests
- Backward compatibility verified
- Serialization and deserialization tested

## Implementation Details

The implementation extends the existing `CodesConfig` class with:

1. **`versions` field** - Dictionary storing version-specific configurations
2. **Version parameter** - Added to `add_code()`, `get_code()`, `has_code()`, `list_codes()`
3. **New methods** - `list_versions()`, `get_version_config()`
4. **Helper function** - `add_version_to_config()` for easy version management

**Connection persistence** is handled by the existing `RemoteExecutionMixin` - no changes needed! The mixin already caches connections by `(host, username)`, so switching QE versions automatically reuses the connection.

## Migration

No migration needed! Existing single-version configurations continue to work:

```python
# Old config - still works!
config = CodesConfig(machine_name="cluster1", qe_version="7.2")
config.add_code(Code(name="pw", path="/usr/bin/pw.x"))

# Access without version parameter
pw = config.get_code("pw")  # Works as before
```

## Summary

This implementation provides a clean, efficient solution to the multi-version QE problem:

1. ✅ Multiple versions can coexist on same machine
2. ✅ Users can easily select version per calculation
3. ✅ SSH connections persist across version switches
4. ✅ Backward compatible with existing code
5. ✅ Well-tested with comprehensive test suite
6. ✅ Documented with examples and API reference
