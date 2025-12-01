# Magnetic Configuration Helpers

## Overview

This module provides simplified helper functions to configure magnetic moments for spin-polarized and antiferromagnetic calculations in Quantum ESPRESSO. These functions automate the process of creating species labels and setting up `starting_magnetization` parameters.

## New: Element-Based Configuration (Recommended!)

The `setup_magnetic_config()` function provides the most intuitive way to define magnetic configurations. Instead of specifying magnetic moments for individual atom indices, you specify them per element type.

### Hubbard Parameters Support

**NEW**: Full support for both old (QE < 7.0) and new (QE 7.x+) Hubbard parameter formats!

#### Old Format (QE < 7.0)
```python
config = setup_magnetic_config(atoms, {
    'Fe': {'mag': [1, -1], 'U': 4.3}
}, qe_version='6.8')
# Result: Hubbard_U in SYSTEM namelist
```

#### New Format (QE 7.x+) - HUBBARD Card
```python
config = setup_magnetic_config(atoms, {
    'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}
}, qe_version='7.2')
# Result: HUBBARD card with Fe-3d orbital specification
```

#### Inter-site V Parameters
```python
config = setup_magnetic_config(atoms, {
    'Fe': {
        'mag': [1],
        'U': {'3d': 4.3},
        'V': [{'species2': 'O', 'orbital1': '3d', 'orbital2': '2p', 'value': 1.0}]
    },
    'O': [0]
}, qe_version='7.2')
# Result: Both U and V parameters in HUBBARD card
```

### Quick Start

```python
from ase.build import bulk
from xespresso import setup_magnetic_config, Espresso

atoms = bulk('Fe', cubic=True)  # 2 Fe atoms

# Simple: All Fe equivalent with magnetization 1
config = setup_magnetic_config(atoms, {'Fe': [1]})

# AFM: Two non-equivalent Fe atoms
config = setup_magnetic_config(atoms, {'Fe': [1, -1]})

# With Hubbard U
config = setup_magnetic_config(atoms, {
    'Fe': {'mag': [1, -1], 'U': 4.3}
})

# Add pseudopotentials and create calculator
config['pseudopotentials']['Fe'] = 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
config['pseudopotentials']['Fe1'] = 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'

calc = Espresso(
    pseudopotentials=config['pseudopotentials'],
    input_data={'input_ntyp': config['input_ntyp']},
    nspin=2,
    ecutwfc=40,
    kpts=(4, 4, 4)
)
config['atoms'].calc = calc
```

### Understanding the Syntax

**Element-based specification:**
- `{'Fe': [1]}` - All Fe atoms equivalent, magnetization = 1
- `{'Fe': [1, -1]}` - Two Fe atoms, non-equivalent (AFM)
- `{'Fe': [1], 'Mn': [1, -1]}` - Multiple elements

**Automatic pattern replication:**
- If you have 4 Fe atoms but specify `{'Fe': [1, -1]}`, the pattern replicates: [1, -1, 1, -1]
- Useful for periodic magnetic structures

**Supercell expansion:**
- `{'Fe': [1, 1, -1, -1]}` with only 2 Fe atoms → Error (by default)
- Set `expand_cell=True` → Automatically creates supercell with 4 Fe atoms

**With Hubbard parameters:**
- `{'Fe': {'mag': [1, -1], 'U': 4.3}}` - Same U for all Fe species
- `{'Fe': {'mag': [1, -1], 'U': [4.3, 4.5]}}` - Different U for each species

## API Reference

### `setup_magnetic_config(atoms, magnetic_config, pseudopotentials=None, expand_cell=False)`

**Parameters:**
- `atoms`: ASE Atoms object
- `magnetic_config`: Dict mapping elements to magnetic moments
  - Simple: `{'Fe': [1, -1]}`
  - With Hubbard: `{'Fe': {'mag': [1, -1], 'U': 4.3}}`
- `pseudopotentials`: Optional base pseudopotentials dict
- `expand_cell`: If True, automatically expand cell when needed

**Returns:**
- `'atoms'`: Updated atoms object (may be supercell if expanded)
- `'input_ntyp'`: Dict with starting_magnetization and optionally Hubbard_U
- `'pseudopotentials'`: Dict mapping species to pseudopotential files
- `'species_map'`: Dict mapping species labels to base elements
- `'expanded'`: Bool indicating if cell was expanded

### Examples

#### Example 1: FeMnAl₂ System

```python
from ase import Atoms
from xespresso import setup_magnetic_config

# Create structure with 2 Fe, 2 Mn, 4 Al
atoms = Atoms('Fe2Mn2Al4', positions=[...])
atoms.cell = [5, 5, 5]

config = setup_magnetic_config(atoms, {
    'Fe': [1],        # Both Fe equivalent
    'Mn': [1, -1],    # Mn AFM
    'Al': [0]         # Al non-magnetic
})

# Result:
# - Fe: one species with mag=1
# - Mn: two species (Mn with mag=1, Mn1 with mag=-1)
# - Al: non-magnetic, not in starting_magnetization
```

#### Example 2: Complex Magnetic Pattern

```python
# 4 Fe atoms with checkerboard pattern
atoms = bulk('Fe', cubic=True) * (2, 1, 1)

config = setup_magnetic_config(atoms, {
    'Fe': [1, -1, -1, 1]
})

# Creates 4 species: Fe, Fe1, Fe2, Fe3
```

#### Example 3: With DFT+U

```python
config = setup_magnetic_config(atoms, {
    'Fe': {'mag': [1, -1], 'U': 4.3},
    'Mn': {'mag': [1, -1], 'U': [5.7, 5.8]}  # Different U values
})

# Result includes both starting_magnetization and Hubbard_U
```

## Original Functions (Still Supported)

## Motivation

Previously, setting up antiferromagnetic or spin-polarized calculations required:
1. Manually creating a species array
2. Manually renaming atoms to different species (e.g., Fe → Fe, Fe1)
3. Manually defining `input_ntyp` dictionary with `starting_magnetization` for each species
4. Manually setting up pseudopotentials for each species

This was error-prone and tedious, especially for complex magnetic structures.

## New API

Three helper functions simplify this process:

### 1. `set_magnetic_moments(atoms, magnetic_moments, pseudopotentials=None)`

The most flexible function - set arbitrary magnetic moments for individual atoms.

**Parameters:**
- `atoms`: ASE Atoms object
- `magnetic_moments`: List, array, or dict of magnetic moments
  - List/array: magnetic moment for each atom in order
  - Dict: `{atom_index: magnetic_moment}` for specific atoms
- `pseudopotentials`: (optional) Existing pseudopotentials dict

**Returns:** Dict with:
- `'input_ntyp'`: Dict with `starting_magnetization`
- `'pseudopotentials'`: Dict mapping species to pseudopotential files
- `'species_map'`: Dict mapping original symbols to species labels

**Example:**
```python
from ase.build import bulk
from xespresso import Espresso, set_magnetic_moments

atoms = bulk('Fe', cubic=True)

# Simple AFM with list
mag_config = set_magnetic_moments(atoms, [1.0, -1.0])

# Add pseudopotential files
mag_config['pseudopotentials']['Fe'] = 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
mag_config['pseudopotentials']['Fe1'] = 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'

# Create calculator
calc = Espresso(
    pseudopotentials=mag_config['pseudopotentials'],
    input_data={'input_ntyp': mag_config['input_ntyp']},
    nspin=2,
    ecutwfc=40,
    kpts=(4, 4, 4)
)
```

### 2. `set_antiferromagnetic(atoms, sublattice_indices, magnetic_moment=1.0, pseudopotentials=None)`

Simplified function specifically for antiferromagnetic configurations.

**Parameters:**
- `atoms`: ASE Atoms object
- `sublattice_indices`: Two sublattices as list of lists
  - Example: `[[0, 2], [1, 3]]` for checkerboard AFM
- `magnetic_moment`: Magnitude of magnetic moment (default=1.0)
- `pseudopotentials`: (optional) Existing pseudopotentials dict

**Returns:** Same as `set_magnetic_moments`

**Example:**
```python
from ase.build import bulk
from xespresso import set_antiferromagnetic

atoms = bulk('Fe', cubic=True)

# Simple AFM: atom 0 has +1.0, atom 1 has -1.0
afm_config = set_antiferromagnetic(atoms, [[0], [1]])

# Complex AFM with 4 atoms
atoms_large = bulk('Fe', cubic=True) * (2, 1, 1)
afm_config = set_antiferromagnetic(
    atoms_large, 
    [[0, 3], [1, 2]],  # Checkerboard pattern
    magnetic_moment=1.5
)
```

### 3. `set_ferromagnetic(atoms, magnetic_moment=1.0, element=None, pseudopotentials=None)`

Simplified function for ferromagnetic configurations.

**Parameters:**
- `atoms`: ASE Atoms object
- `magnetic_moment`: Magnetic moment for all atoms (default=1.0)
- `element`: (optional) Only set magnetization for this element
- `pseudopotentials`: (optional) Existing pseudopotentials dict

**Returns:** Same as `set_magnetic_moments`

**Example:**
```python
from ase.build import bulk
from xespresso import set_ferromagnetic

atoms = bulk('Fe', cubic=True)

# Ferromagnetic with all atoms at 2.0
fm_config = set_ferromagnetic(atoms, magnetic_moment=2.0)

# Mixed system: only Fe is magnetic
from ase import Atoms
atoms_mixed = Atoms('Fe2O2', positions=[[0,0,0], [1.5,0,0], [0,1.5,0], [1.5,1.5,0]])
atoms_mixed.cell = [5, 5, 5]

fm_config = set_ferromagnetic(atoms_mixed, magnetic_moment=2.0, element='Fe')
```

## Comparison: Old vs New

### Old Method (Manual)

```python
import numpy as np
from ase.build import bulk
from xespresso import Espresso

atoms = bulk('Fe', cubic=True)

# Manually create species array
atoms.new_array('species', np.array(atoms.get_chemical_symbols(), dtype='U20'))

# Manually rename atoms
atoms.arrays['species'][0] = 'Fe'
atoms.arrays['species'][1] = 'Fe1'

# Manually create input_ntyp
input_ntyp = {
    'starting_magnetization': {
        'Fe': 1.0,
        'Fe1': -1.0,
    }
}

# Manually define pseudopotentials for each species
pseudopotentials = {
    'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Fe1': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
}

calc = Espresso(
    pseudopotentials=pseudopotentials,
    input_data={'input_ntyp': input_ntyp},
    nspin=2,
    ecutwfc=40,
    kpts=(4, 4, 4)
)
```

### New Method (Simplified)

```python
from ase.build import bulk
from xespresso import Espresso, set_antiferromagnetic

atoms = bulk('Fe', cubic=True)

# One function call to set up everything
mag_config = set_antiferromagnetic(atoms, [[0], [1]])

# Just add the pseudopotential file names
mag_config['pseudopotentials']['Fe'] = 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
mag_config['pseudopotentials']['Fe1'] = 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'

calc = Espresso(
    pseudopotentials=mag_config['pseudopotentials'],
    input_data={'input_ntyp': mag_config['input_ntyp']},
    nspin=2,
    ecutwfc=40,
    kpts=(4, 4, 4)
)
```

## Advanced Examples

### Complex AFM with DFT+U

```python
from ase.build import bulk
from xespresso import Espresso, set_magnetic_moments

atoms = bulk('Mn', cubic=True) * (2, 2, 1)

# Checkerboard AFM pattern
magnetic_moments = {
    0: 1.0,   # Mn1
    1: -1.0,  # Mn2
    2: -1.0,  # Mn3
    3: 1.0,   # Mn4
}

base_pseudo = {'Mn': 'Mn.pbesol-spn-rrkjus_psl.0.3.1.UPF'}
mag_config = set_magnetic_moments(atoms, magnetic_moments, pseudopotentials=base_pseudo)

# Add Hubbard U for all Mn species
mag_config['input_ntyp']['Hubbard_U'] = {}
for species in mag_config['pseudopotentials'].keys():
    if 'Mn' in species:
        mag_config['input_ntyp']['Hubbard_U'][species] = 5.75

input_data = {
    'ecutwfc': 70.0,
    'nspin': 2,
    'lda_plus_u': True,
    'input_ntyp': mag_config['input_ntyp'],
}

calc = Espresso(
    pseudopotentials=mag_config['pseudopotentials'],
    input_data=input_data,
    kpts=(4, 4, 4)
)
```

### Mixed Elements with Selective Magnetization

```python
from ase import Atoms
from xespresso import set_magnetic_moments

# Create FeO structure
atoms = Atoms('Fe2O2', positions=[
    [0, 0, 0], [1.5, 0, 0],  # Fe atoms
    [0, 1.5, 0], [1.5, 1.5, 0]  # O atoms
])
atoms.cell = [5, 5, 5]

# AFM on Fe, no magnetization on O
mag_config = set_magnetic_moments(atoms, {
    0: 1.0,   # Fe
    1: -1.0,  # Fe
    2: 0.0,   # O (zero moments are automatically excluded)
    3: 0.0    # O
})

# Result: Only Fe atoms appear in starting_magnetization
```

## Benefits

1. **Less Error-Prone**: Automatic species labeling prevents manual mistakes
2. **Less Verbose**: Reduces boilerplate code by ~50%
3. **More Readable**: Intent is clear from function names
4. **Flexible**: Supports simple and complex magnetic structures
5. **Backward Compatible**: Old manual method still works

## Migration Guide

Existing code using the manual method will continue to work. To migrate to the new API:

1. Replace manual species array creation and labeling with helper functions
2. Use the returned `mag_config` dictionary instead of manually creating `input_ntyp`
3. Update pseudopotentials using `mag_config['pseudopotentials']`

Both approaches can coexist in the same project.

## Testing

Comprehensive tests are available in `tests/test_magnetic_helpers.py`. Run with:

```bash
pytest tests/test_magnetic_helpers.py -v
```

## Examples

See working examples in:
- `examples/ex03-spin-simplified.py` - Simple AFM Fe
- `examples/ex02-MnO-afm-simplified.py` - Complex AFM with DFT+U
- `examples/magnetic_helpers_examples.py` - Complete demonstration of all features
