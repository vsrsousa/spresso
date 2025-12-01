# Pseudopotential Auto-Extraction Enhancement

## Overview

This enhancement simplifies the workflow when using `setup_magnetic_config` by allowing users to pass pseudopotentials inside `input_data`, eliminating the need to manually extract and pass them as a separate parameter.

## Problem Statement

Previously, when using `setup_magnetic_config`, users had to:

1. Call `setup_magnetic_config` to generate configuration
2. Manually extract `pseudopotentials` from the config dict
3. Pass `pseudopotentials` separately to the Espresso calculator or write_espresso_in

Example of the old workaround:
```python
config = setup_magnetic_config(atoms, {'Fe': [1, -1]}, 
                               pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'})

# Had to manually extract and pass separately
calc = Espresso(
    pseudopotentials=config['pseudopotentials'],  # Manual extraction required
    input_data={'input_ntyp': config['input_ntyp']},
    nspin=2,
    ecutwfc=40
)
```

## Solution

The Espresso calculator and `write_espresso_in` now automatically extract pseudopotentials from `input_data` if present. This provides three convenient ways to use `setup_magnetic_config`:

### Method 1: Pass Entire Config Dict (Recommended)

```python
config = setup_magnetic_config(atoms, {'Fe': [1, -1]}, 
                               pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'})

# Pass entire config - most convenient!
calc = Espresso(
    atoms=config['atoms'],
    input_data=config,  # Pass entire dict
    nspin=2,
    ecutwfc=40
)
```

### Method 2: Pass Pseudopotentials Inside input_data

```python
config = setup_magnetic_config(atoms, {'Fe': [1, -1]}, 
                               pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'})

# Explicitly structure input_data
calc = Espresso(
    atoms=config['atoms'],
    input_data={
        'input_ntyp': config['input_ntyp'],
        'pseudopotentials': config['pseudopotentials']  # Inside input_data
    },
    nspin=2,
    ecutwfc=40
)
```

### Method 3: Old Method (Still Supported)

```python
config = setup_magnetic_config(atoms, {'Fe': [1, -1]}, 
                               pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'})

# Old way still works
calc = Espresso(
    atoms=config['atoms'],
    pseudopotentials=config['pseudopotentials'],  # Separate parameter
    input_data={'input_ntyp': config['input_ntyp']},
    nspin=2,
    ecutwfc=40
)
```

## Advanced Features

### Auto-Derivation from Base Elements

When you provide only base element pseudopotentials, the system automatically derives them for species:

```python
# Only provide base element
config = setup_magnetic_config(atoms, {'Fe': [1, -1]}, 
                               pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'})

# Fe1 and Fe2 automatically inherit from Fe
assert config['pseudopotentials']['Fe1'] == 'Fe.pbe-spn.UPF'
assert config['pseudopotentials']['Fe2'] == 'Fe.pbe-spn.UPF'
```

This uses the `species_map` that `setup_magnetic_config` creates:
```python
config['species_map']  # {'Fe1': 'Fe', 'Fe2': 'Fe'}
```

### Multiple Elements

```python
atoms = Atoms('Fe2Mn2', positions=[...])
config = setup_magnetic_config(atoms, {
    'Fe': [1],      # Both Fe equivalent
    'Mn': [1, -1]   # Mn AFM
}, pseudopotentials={
    'Fe': 'Fe.pbe-spn.UPF',
    'Mn': 'Mn.pbe-spn.UPF'
})

# All species get correct pseudopotentials:
# Fe, Mn1, Mn2 all mapped correctly
```

### Override with Top-Level Parameter

If you need to override specific pseudopotentials, you can still pass them at the top level:

```python
calc = Espresso(
    atoms=config['atoms'],
    pseudopotentials={'Fe1': 'Fe_custom.UPF'},  # Override Fe1
    input_data=config,  # Fe2 from config
    nspin=2,
    ecutwfc=40
)
# Fe1 uses Fe_custom.UPF, Fe2 uses Fe.pbe-spn.UPF
```

## Technical Details

### Implementation

The enhancement is implemented in two places:

1. **`xespresso.xio.sort_qe_input`**: Extracts pseudopotentials from `input_data` and uses `species_map` to auto-derive pseudopotentials for species
2. **`xespresso.xio.write_espresso_in`**: Extracts pseudopotentials from `input_data` before writing input file

### Metadata Cleanup

The following metadata fields from `setup_magnetic_config` are automatically removed from `input_data` as they are not QE parameters:
- `expanded`
- `hubbard_format`
- `atoms`
- `species_map`

### Merge Priority

When pseudopotentials appear in multiple places:
1. Top-level `pseudopotentials` parameter takes highest priority
2. Pseudopotentials in `input_data` are added for missing species
3. Auto-derivation from `species_map` fills in remaining gaps

## Benefits

1. **Simpler Code**: Reduce boilerplate when using `setup_magnetic_config`
2. **Less Error-Prone**: No manual extraction means fewer mistakes
3. **Backward Compatible**: Old method still works
4. **Flexible**: Multiple ways to specify pseudopotentials
5. **Smart Defaults**: Auto-derivation from base elements

## Examples

### Example 1: Simple AFM System

```python
from ase.build import bulk
from xespresso import Espresso
from xespresso.tools import setup_magnetic_config

atoms = bulk('Fe', cubic=True)
config = setup_magnetic_config(
    atoms, 
    {'Fe': [1, -1]},
    pseudopotentials={'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'}
)

# One-liner setup!
calc = Espresso(
    atoms=config['atoms'],
    input_data=config,
    nspin=2,
    ecutwfc=40,
    kpts=(4, 4, 4)
)
```

### Example 2: With Hubbard Parameters

```python
config = setup_magnetic_config(
    atoms,
    {'Fe': {'mag': [1, -1], 'U': 4.3}},
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'}
)

calc = Espresso(
    atoms=config['atoms'],
    input_data=config,  # Includes Hubbard parameters
    nspin=2,
    lda_plus_u=True,
    ecutwfc=40,
    kpts=(4, 4, 4)
)
```

### Example 3: Complex Multi-Element System

```python
atoms = Atoms('Fe2Mn2Al4', positions=[...])
config = setup_magnetic_config(atoms, {
    'Fe': [1],        # FM
    'Mn': [1, -1],    # AFM
    'Al': [0]         # Non-magnetic
}, pseudopotentials={
    'Fe': 'Fe.pbe-spn.UPF',
    'Mn': 'Mn.pbe-spn.UPF',
    'Al': 'Al.pbe.UPF'
})

# All pseudopotentials automatically handled
calc = Espresso(
    atoms=config['atoms'],
    input_data=config,
    nspin=2,
    ecutwfc=40
)
```

## Migration Guide

### Before (Old Way)

```python
config = setup_magnetic_config(atoms, {'Fe': [1, -1]}, 
                               pseudopotentials={'Fe': 'Fe.UPF'})
# Manual extraction
calc = Espresso(
    pseudopotentials=config['pseudopotentials'],
    input_data={'input_ntyp': config['input_ntyp']},
    nspin=2
)
```

### After (New Way)

```python
config = setup_magnetic_config(atoms, {'Fe': [1, -1]}, 
                               pseudopotentials={'Fe': 'Fe.UPF'})
# No manual extraction needed
calc = Espresso(
    input_data=config,  # That's it!
    nspin=2
)
```

## Testing

The enhancement includes comprehensive tests:
- 10 unit tests in `test_pseudopotential_extraction.py`
- 6 integration tests in `test_integration_pseudopotentials.py`
- All existing tests continue to pass

## Compatibility

- **Backward Compatible**: All existing code continues to work
- **QE Versions**: Works with all Quantum ESPRESSO versions
- **ASE Version**: Compatible with ASE 3.22.1+

## See Also

- [MAGNETIC_HELPERS.md](MAGNETIC_HELPERS.md) - Documentation for `setup_magnetic_config`
- [IMPLEMENTATION_SUMMARY_MAGNETIC.md](IMPLEMENTATION_SUMMARY_MAGNETIC.md) - Implementation details of magnetic helpers
