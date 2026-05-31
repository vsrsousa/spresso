# enhance_nbands Feature Guide

## Overview

The `enhance_nbands` parameter provides **exact calculation of band count (nbnd)** based on the structure composition and pseudopotential valence electron counts. This is useful for band structure and wave function analysis calculations where you need nbnd to precisely match the total number of valence electrons.

## Formula

```
nbnd = Σ (N_atoms[element] × z_valence[element])
```

Where:
- `N_atoms[element]` = number of atoms of that element in the structure
- `z_valence[element]` = valence electrons per atom (read from UPF file)

## Example Calculation

For Si₂ (2 Si atoms):
- Si has 4 valence electrons
- Formula: nbnd = 2 × 4 = **8 bands**

For Fe₂O₃:
- Fe has 6 valence electrons (3d⁶4s²)
- O has 6 valence electrons (2p⁶)  
- Formula: nbnd = (2 × 6) + (3 × 6) = 12 + 18 = **30 bands**

## Usage

### With CalculationWorkflow

```python
from xespresso.workflow import CalculationWorkflow
from ase.build import bulk

atoms = bulk('Si', 'diamond', a=5.43)

# Create workflow with enhance_nbands
wf = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={'Si': 'Si.pbe.UPF'},
    protocol='moderate',
    enhance_nbands=True  # Use exact nbnd calculation
)

# Run calculations (will use exact nbnd = 8)
calc = wf.run_scf(label='si_scf')
```

### With EOSWorkflow

```python
from xespresso.workflow import EOSWorkflow

eos = EOSWorkflow(
    atoms=atoms,
    pseudopotentials={'Si': 'Si.pbe.UPF'},
    protocol='moderate',
    enhance_nbands=True  # All volume points will use exact nbnd
)

# Run EOS study
eos.run_eos_study()
```

### With ConvergenceWorkflow

```python
from xespresso.workflow import ConvergenceWorkflow

conv = ConvergenceWorkflow(
    atoms=atoms,
    pseudopotentials={'Si': 'Si.pbe.UPF'},
    protocol='moderate',
    enhance_nbands=True  # Convergence calculations will use exact nbnd
)

# Run convergence study
results = conv.run_convergence_study()
```

## Default Behavior (enhance_nbands=False)

When `enhance_nbands=False` (default), the system uses the traditional estimation:

```python
nbnd = suggest_nbnd_from_pseudos(pseudopotentials, buffer=20)
```

This adds a safety buffer but may overestimate band count for band structure calculations.

## When to Use

### ✓ Use enhance_nbands=True when:
- Performing band structure calculations (need exact occupied bands)
- Computing optical properties where nbnd must match electron count
- Analyzing band structure properties that depend on exact band count
- Working with structures where occupancy is known

### ✗ Use enhance_nbands=False (default) when:
- Running SCF calculations where overestimating nbnd is safer
- Working with mixed valency or partially occupied states
- Using spin-polarized calculations (where nbnd is effectively doubled)
- You want a conservative buffer for safety margin

## How It Works

1. **Parameter Propagation**: The `enhance_nbands` parameter is stored in each workflow
2. **Pass-through**: When creating child CalculationWorkflow instances, the parameter is passed along
3. **Calculation Trigger**: During `_estimate_nbnd()` call, if `enhance_nbands=True`:
   - Calls `calculate_nbnd_from_structure()` function
   - Reads z_valence from each pseudopotential's UPF file
   - Counts atoms of each element in structure
   - Computes nbnd = Σ(count × valence)
4. **Default Fallback**: Falls back to buffer-based estimation if pseudopotential files aren't found

## Implementation Details

### File Modifications

1. **xespresso/utils/pseudo_utils.py**
   - Added `calculate_nbnd_from_structure()` function
   - Reads z_valence from UPF files (XML and text formats)
   - Handles missing files gracefully with default z_valence=8

2. **xespresso/workflow/calculation_workflow.py**
   - Added `enhance_nbands` parameter to `__init__()`
   - Modified `_estimate_nbnd()` to check flag and use exact calculation

3. **xespresso/workflow/eos_workflow.py**
   - Added `enhance_nbands` parameter to `__init__()`
   - Passes to all 3 CalculationWorkflow instantiation points

4. **xespresso/workflow/convergence_workflow.py**
   - Added `enhance_nbands` parameter to `__init__()`
   - Passes to both CalculationWorkflow instantiation points via wf_kwargs

### XML Format Support

The implementation supports multiple z_valence formats in UPF files:

```xml
<!-- Format 1: Standard attribute -->
<z_valence>4</z_valence>

<!-- Format 2: Key-value with equals -->
z_valence = "4"

<!-- Format 3: Key-value with colon -->
z_valence: 4

<!-- Format 4: Uppercase -->
Z_valence = "4"
```

## Testing

Two test files demonstrate the feature:

1. **test_enhance_nbands.py**: Basic functionality test
   - Shows nbnd calculation with different structures
   - Compares with default behavior
   - Verifies formula correctness

2. **test_integration_enhance_nbands.py**: Full integration test
   - Tests all three workflows (CalculationWorkflow, EOSWorkflow, ConvergenceWorkflow)
   - Verifies parameter propagation
   - Confirms nbnd values match expected calculations

Run tests:
```bash
python test_enhance_nbands.py
python test_integration_enhance_nbands.py
```

## Troubleshooting

### Issue: nbnd not changing with enhance_nbands=True

**Check:**
1. Ensure UPF files are accessible (not missing or unreadable)
2. Verify z_valence is present in UPF file headers
3. Check that pseudopotentials dict has correct element mapping

### Issue: Unexpected nbnd values

**Common causes:**
- Spin-polarized calculations: nbnd is for each spin channel separately
- Magnetic moments: Each atom may have different valence treatment
- Hubbard U: May affect effective electron count

## Related Functions

- `calculate_nbnd_from_structure()`: Core calculation (xespresso/utils/pseudo_utils.py)
- `suggest_nbnd_from_pseudos()`: Traditional buffer-based estimation (xespresso/workflow/wannier_workflow.py)
- `_estimate_nbnd()`: Dispatcher method (xespresso/workflow/calculation_workflow.py)

## References

- UPF Format: https://www.quantum-espresso.org/fileformats/
- Quantum ESPRESSO nbnd Documentation: https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm45922584993648
