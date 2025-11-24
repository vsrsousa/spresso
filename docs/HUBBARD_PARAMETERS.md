# Hubbard Parameters: New Format for QE 7.0+

## Overview

Starting with Quantum ESPRESSO 7.0, a new `HUBBARD` card format was introduced that provides more flexibility and clarity for DFT+U and DFT+U+V calculations. The xespresso library now supports both the old and new formats, with automatic detection and backward compatibility.

## Key Changes in QE 7.0+

### Old Format (QE < 7.0)
Hubbard parameters were specified in the `SYSTEM` namelist:
```
&SYSTEM
  ...
  lda_plus_u = .true.
  Hubbard_U(1) = 4.3
  Hubbard_U(2) = 3.0
  Hubbard_V(1,2,1) = 0.5
  ...
/
```

### New Format (QE >= 7.0)
A dedicated `HUBBARD` card with explicit orbital specification:
```
HUBBARD {atomic}
  U Fe-3d 4.3
  U O-2p 3.0
  V Fe-3d O-2p 1 1 0.5
```

## Using the New Format in xespresso

### Basic Example (New Format)

```python
from ase.build import bulk
import numpy as np
from xespresso import Espresso

atoms = bulk("Fe", cubic=True)

# Use the new HUBBARD card format
input_data = {
    "ecutwfc": 30.0,
    "occupations": "smearing",
    "degauss": 0.02,
    "nspin": 2,
    "lda_plus_u": True,
    "qe_version": "7.2",  # Specify QE version
    "hubbard": {
        "projector": "atomic",  # or 'ortho-atomic', 'norm-atomic', 'wf', 'pseudo'
        "u": {
            "Fe-3d": 4.3,  # Explicit orbital specification
        }
    }
}

pseudopotentials = {
    "Fe": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF",
}

calc = Espresso(
    pseudopotentials=pseudopotentials,
    label="scf/fe",
    input_data=input_data,
    kpts=(4, 4, 4),
)
atoms.calc = calc
energy = atoms.get_potential_energy()
```

### DFT+U+V Example (New Format)

```python
input_data = {
    "ecutwfc": 30.0,
    "occupations": "smearing",
    "degauss": 0.02,
    "nspin": 2,
    "lda_plus_u": True,
    "qe_version": "7.2",
    "hubbard": {
        "projector": "atomic",
        "u": {
            "Fe1-3d": 4.3,
            "Fe2-3d": 4.3,
            "O-2p": 3.0,
        },
        "v": [
            {
                "species1": "Fe1",
                "orbital1": "3d",
                "species2": "Fe2",
                "orbital2": "3d",
                "i": 1,
                "j": 1,
                "value": 0.5
            },
            {
                "species1": "Fe1",
                "orbital1": "3d",
                "species2": "O",
                "orbital2": "2p",
                "i": 1,
                "j": 1,
                "value": 0.3
            }
        ]
    }
}
```

## Backward Compatibility (Old Format)

The old format is still supported and works exactly as before:

```python
input_ntyp = {
    "starting_magnetization": {
        "Fe": 0.5,
    },
    "Hubbard_U": {
        "Fe": 4.3,
        "O": 3.0,
    },
}

input_data = {
    "ecutwfc": 30.0,
    "occupations": "smearing",
    "degauss": 0.02,
    "nspin": 2,
    "lda_plus_u": True,
    "input_ntyp": input_ntyp,
    # For V parameters in old format:
    "hubbard_v": {"(1,2,1)": 0.5},
}
```

## Automatic Format Detection

xespresso automatically detects which format to use based on:

1. **Explicit QE version**: If `qe_version >= "7.0"`, use new format
2. **Presence of `hubbard` dictionary**: Indicates new format
3. **Orbital specifications**: Keys like `"Fe-3d"` indicate new format
4. **Default**: Use old format for backward compatibility

### Force a Specific Format

```python
# Force new format
input_data['hubbard_format'] = 'card'

# Force old format
input_data['hubbard_format'] = 'namelist'
```

## Projector Types

The new format supports different projector types:

- `atomic`: Atomic orbitals (default)
- `ortho-atomic`: Orthogonalized atomic orbitals
- `norm-atomic`: Normalized atomic orbitals
- `wf`: Wannier functions
- `pseudo`: Pseudopotential orbitals

```python
input_data["hubbard"]["projector"] = "ortho-atomic"
```

## Benefits of New Format

1. ✅ **Explicit orbital specification**: More precise and clear
2. ✅ **Easier DFT+U+V**: Simpler specification of inter-site interactions
3. ✅ **Multiple Hubbard sites**: Better support for complex systems
4. ✅ **Flexible projectors**: More control over projection method
5. ✅ **Cleaner input files**: More readable and maintainable
6. ✅ **Future-proof**: Recommended for QE 7.0 and later
7. ✅ **Backward compatible**: Old format still works

## Migration Guide

### Converting Old to New Format

**Old Format:**
```python
input_data = {
    "lda_plus_u": True,
    "input_ntyp": {
        "Hubbard_U": {
            "Fe": 4.3,
            "O": 3.0
        }
    },
    "hubbard_v": {"(1,2,1)": 0.5}
}
```

**New Format:**
```python
input_data = {
    "lda_plus_u": True,
    "qe_version": "7.2",
    "hubbard": {
        "projector": "atomic",
        "u": {
            "Fe-3d": 4.3,
            "O-2p": 3.0
        },
        "v": [
            {
                "species1": "Fe",
                "orbital1": "3d",
                "species2": "O",
                "orbital2": "2p",
                "i": 1,
                "j": 2,
                "value": 0.5
            }
        ]
    }
}
```

## Configuration with Codes Module

When you detect QE codes, the version is automatically determined:

```python
from xespresso.codes import detect_qe_codes

config = detect_qe_codes(machine_name="cluster")
if config.qe_version:
    print(f"Detected QE version: {config.qe_version}")
    
    # Use this version in your calculations
    input_data["qe_version"] = config.qe_version
```

## Examples

See `examples/hubbard_new_format_example.py` for complete examples.

## API Reference

### HubbardConfig Class

The `HubbardConfig` class from `xespresso.hubbard` handles both formats:

```python
from xespresso.hubbard import HubbardConfig

# Create config with new format
config = HubbardConfig(use_new_format=True, projector='atomic')
config.add_u("Fe", 4.3, orbital="3d")
config.add_v("Fe", "O", 0.5, orbital1="3d", orbital2="2p")

# Generate HUBBARD card
lines = config.to_new_format_card()

# Or create from input_data
config = HubbardConfig.from_input_data(input_data, qe_version="7.2")
```

## Troubleshooting

### Issue: Old format not working with QE 7.x

**Solution**: The old format is deprecated but still works. However, for new projects, use the new format. If you need to use old format with QE 7.x, explicitly set:
```python
input_data['hubbard_format'] = 'namelist'
```

### Issue: Don't know which orbitals to use

**Solution**: Common orbitals for DFT+U:
- Transition metals: `3d` (3d series), `4d` (4d series), `5d` (5d series), `4f`, `5f`
- Oxygen: `2p`
- Sulfur: `3p`

Check the pseudopotential file or QE documentation for your specific elements.

### Issue: Different U values for same element

**Solution**: Use different species names:
```python
atoms.new_array("species", np.array(atoms.get_chemical_symbols(), dtype="U20"))
atoms.arrays["species"][0] = "Fe1"
atoms.arrays["species"][1] = "Fe2"

input_data["hubbard"]["u"] = {
    "Fe1-3d": 4.3,
    "Fe2-3d": 5.0
}
```

## References

- [QE 7.0 Release Notes](https://www.quantum-espresso.org)
- [INPUT_PW documentation](https://www.quantum-espresso.org/Doc/INPUT_PW.html)
- Cococcioni & de Gironcoli, PRB 71, 035105 (2005)
- Timrov et al., PRB 98, 085127 (2018)
