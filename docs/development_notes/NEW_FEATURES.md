# New Features: Codes Module and Hubbard Parameters Update

## Overview

This update addresses important issues with xespresso:

1. **Codes Module**: A new module for managing Quantum ESPRESSO code configurations across different machines
2. **Multiple QE Versions**: Support for managing multiple QE versions on the same machine
3. **Hubbard Parameters**: Updated support for the new HUBBARD card format introduced in Quantum ESPRESSO 7.0+

## 1. Codes Module

### What is it?

The codes module (`xespresso.codes`) provides utilities to:
- Detect available QE executables on a machine
- Store QE code paths and versions
- Configure codes per machine
- Support both local and remote machines
- **NEW: Support multiple QE versions on the same machine**

### Why do we need it?

Different machines may have QE installed in different locations, with different versions, and require different module loading. The codes module allows you to:
- Store machine-specific QE configurations
- Automatically detect available codes
- Easily switch between machines
- Integrate with the existing machines configuration
- **NEW: Use different QE versions for different calculations**
- **NEW: Maintain SSH connection when switching versions**

### Quick Start

#### Single Version

```python
from xespresso.codes import detect_qe_codes, load_codes_config

# Detect codes on a machine
config = detect_qe_codes(
    machine_name="cluster1",
    qe_prefix="/opt/qe-7.2/bin",
    modules=['quantum-espresso/7.2']
)

# Configuration is automatically saved to ~/.xespresso/codes/cluster1.json

# Later, load the configuration
codes = load_codes_config("cluster1")
if codes.has_code('pw'):
    pw_code = codes.get_code('pw')
    print(f"pw.x: {pw_code.path} (version {pw_code.version})")
```

#### Multiple Versions (NEW!)

```python
from xespresso.codes import add_version_to_config, load_codes_config

# Add QE 7.2
add_version_to_config(
    machine_name="cluster1",
    version="7.2",
    qe_prefix="/opt/qe-7.2/bin",
    modules=["quantum-espresso/7.2"]
)

# Add QE 6.8 to the same machine
add_version_to_config(
    machine_name="cluster1",
    version="6.8",
    qe_prefix="/opt/qe-6.8/bin",
    modules=["quantum-espresso/6.8"]
)

# Use different versions in your workflow
codes = load_codes_config("cluster1")

# SCF with QE 7.2
pw_72 = codes.get_code("pw", version="7.2")

# DOS with QE 6.8
dos_68 = codes.get_code("dos", version="6.8")

# SSH connection persists throughout!
```

### Documentation

- [docs/CODES_CONFIGURATION.md](docs/CODES_CONFIGURATION.md) - Complete codes module documentation
- [docs/MULTIPLE_VERSIONS.md](docs/MULTIPLE_VERSIONS.md) - **NEW: Multiple QE versions guide**
- [docs/REMOTE_CONNECTION_PERSISTENCE.md](docs/REMOTE_CONNECTION_PERSISTENCE.md) - Connection management details

### Examples

- [examples/codes_example.py](examples/codes_example.py) - Basic codes configuration
- [examples/multi_version_example.py](examples/multi_version_example.py) - **NEW: Multi-version setup**
- [examples/workflow_multi_version.py](examples/workflow_multi_version.py) - **NEW: Practical workflow**

## 2. Hubbard Parameters Update

### What changed?

Quantum ESPRESSO 7.0 introduced a new `HUBBARD` card format that replaces the old `Hubbard_U(i)` parameters in the SYSTEM namelist. The new format is more flexible and allows explicit orbital specification.

### Old Format (QE < 7.0)

```python
input_data = {
    "lda_plus_u": True,
    "input_ntyp": {
        "Hubbard_U": {
            "Fe": 4.3,
            "O": 3.0
        }
    }
}
```

Generates:
```
&SYSTEM
  ...
  Hubbard_U(1) = 4.3
  Hubbard_U(2) = 3.0
  ...
/
```

### New Format (QE >= 7.0)

```python
input_data = {
    "lda_plus_u": True,
    "qe_version": "7.2",
    "hubbard": {
        "projector": "atomic",
        "u": {
            "Fe-3d": 4.3,
            "O-2p": 3.0
        }
    }
}
```

Generates:
```
HUBBARD {atomic}
  U Fe-3d 4.3
  U O-2p 3.0
```

### Benefits of New Format

1. ✅ Explicit orbital specification (e.g., `Fe-3d` instead of just `Fe`)
2. ✅ Easier DFT+U+V specification with inter-site interactions
3. ✅ Better support for multiple Hubbard sites
4. ✅ More flexible projector choices
5. ✅ Recommended for QE 7.0 and later

### Backward Compatibility

**The old format still works!** xespresso maintains full backward compatibility:
- Old format input continues to work
- Automatic format detection based on QE version
- Both formats can coexist in the codebase

### Quick Start

```python
from xespresso import Espresso
from ase.build import bulk

atoms = bulk("Fe", cubic=True)

# Use new format
input_data = {
    "ecutwfc": 30.0,
    "nspin": 2,
    "lda_plus_u": True,
    "qe_version": "7.2",
    "hubbard": {
        "projector": "atomic",
        "u": {"Fe-3d": 4.3}
    }
}

calc = Espresso(
    pseudopotentials={"Fe": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF"},
    label="scf/fe",
    input_data=input_data,
    kpts=(4, 4, 4)
)
atoms.calc = calc
```

### Documentation

See [docs/HUBBARD_PARAMETERS.md](docs/HUBBARD_PARAMETERS.md) for complete documentation.

### Examples

See [examples/hubbard_new_format_example.py](examples/hubbard_new_format_example.py) for examples.

## Integration

The two new features work together:

```python
from xespresso.codes import load_codes_config
from xespresso import Espresso

# Load codes configuration
codes = load_codes_config("cluster1")

# Use QE version from codes config
input_data = {
    "ecutwfc": 30.0,
    "qe_version": codes.qe_version,  # Automatically use correct format
    "hubbard": {
        "u": {"Fe-3d": 4.3}
    }
}

# Rest of calculation setup...
```

## Testing

All features are thoroughly tested:
- **15 tests** for codes module
- **14 tests** for Hubbard parameters
- **Integration tests** verifying both modules work together
- **Backward compatibility tests** ensuring old code still works

Run tests:
```bash
pytest tests/test_codes.py tests/test_hubbard.py -v
```

## File Changes

### New Files
- `xespresso/codes/__init__.py` - Codes module public API
- `xespresso/codes/config.py` - Code and CodesConfig dataclasses
- `xespresso/codes/manager.py` - CodesManager class
- `xespresso/hubbard.py` - Hubbard parameters handling
- `tests/test_codes.py` - Tests for codes module
- `tests/test_hubbard.py` - Tests for Hubbard parameters
- `tests/test_integration.py` - Integration tests
- `examples/codes_example.py` - Codes module examples
- `examples/hubbard_new_format_example.py` - Hubbard examples
- `docs/CODES_CONFIGURATION.md` - Codes documentation
- `docs/HUBBARD_PARAMETERS.md` - Hubbard documentation

### Modified Files
- `xespresso/xio.py` - Added Hubbard card support, updated write_espresso_in

## Migration Guide

### For Existing Code

**No migration required!** Your existing code will continue to work as before.

### To Use New Features

**Codes Module:**
1. Create codes configuration: `detect_qe_codes(machine_name="your_machine")`
2. Configuration is saved to `~/.xespresso/codes/your_machine.json`
3. Load in future runs: `load_codes_config("your_machine")`

**Hubbard Parameters (Optional):**
1. If using QE 7.0+, add `"qe_version": "7.2"` to input_data
2. Replace `input_ntyp` with new `hubbard` dictionary
3. Add orbital specifications (e.g., `"Fe-3d"` instead of `"Fe"`)

## Support

- For codes module: See [docs/CODES_CONFIGURATION.md](docs/CODES_CONFIGURATION.md)
- For Hubbard parameters: See [docs/HUBBARD_PARAMETERS.md](docs/HUBBARD_PARAMETERS.md)
- For issues: Open an issue on GitHub

## References

- [Quantum ESPRESSO 7.0 Release](https://www.quantum-espresso.org)
- Cococcioni & de Gironcoli, PRB 71, 035105 (2005) - DFT+U
- Timrov et al., PRB 98, 085127 (2018) - DFT+U+V
