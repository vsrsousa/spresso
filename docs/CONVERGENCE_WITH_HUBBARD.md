## Hubbard Parameters in ConvergenceWorkflow

### What Was Fixed

Previously, the `ConvergenceWorkflow` had a `hubbard_config` parameter but **wasn't actually forwarding it** to the `CalculationWorkflow` instances during convergence studies. This meant users couldn't use DFT+U in convergence studies.

**Status: FIXED** ✅
- Both Phase 1 (ecutwfc convergence) and Phase 2 (kspacing convergence) now properly forward `hubbard_config`
- Works alongside magnetic configuration for realistic strong-correlation systems

### Hubbard Support in xespresso

xespresso supports **two formats** for Hubbard parameters:

#### Format 1: Old Format (QE < 7.0)
```python
hubbard_config = {
    'Fe': 4.3,      # U value in eV (applied to all Fe)
    'Mn': 5.7,
}
```

#### Format 2: New Format (QE >= 7.0)
```python
hubbard_config = {
    'qe_version': '7.2',
    'projector': 'atomic',  # or 'ortho-atomic', 'norm-atomic', 'wf', 'pseudo'
    'u': {
        'Fe-3d': 4.3,      # Explicit orbital specification
        'Mn-3d': 5.7,
        'O-2p': 3.0,
    },
    'v': [                 # Optional: inter-site V parameters
        {
            'species1': 'Fe',
            'orbital1': '3d',
            'species2': 'O',
            'orbital2': '2p',
            'i': 1,
            'j': 1,
            'value': 0.5
        }
    ]
}
```

### Usage with ConvergenceWorkflow

```python
from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow

atoms = bulk('Fe', cubic=True)

# Setup with AFM magnetism + Hubbard U
workflow = ConvergenceWorkflow(
    atoms=atoms,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
    precision='low',
    
    # Magnetic configuration
    magnetic_config={
        'Fe': {'mag': [1.0, -1.0]}  # AFM
    },
    
    # Hubbard U parameters
    hubbard_config={'Fe': 4.3},    # Old format
    
    # OR new format:
    # hubbard_config={
    #     'qe_version': '7.2',
    #     'u': {'Fe-3d': 4.3}
    # }
)

# Run convergence study - will apply both magnetism AND Hubbard to all calculations
results = workflow.run_convergence_study(
    label='test/fe_afm_dft+u',
    max_iterations=3
)
```

### What Gets Passed Through

When `ConvergenceWorkflow.run_convergence_study()` creates `CalculationWorkflow` instances:

**Phase 1 (ecutwfc convergence):**
```
CalculationWorkflow(
    atoms=...,
    pseudopotentials=...,
    protocol=...,
    kspacing=FIXED,         # Fixed during Phase 1
    ecutwfc=VARYING,        # Varying parameter
    magnetic_config=...,    # ✅ Now forwarded
    hubbard_config=...,     # ✅ Now forwarded (FIXED)
)
```

**Phase 2 (kspacing convergence):**
```
CalculationWorkflow(
    atoms=...,
    pseudopotentials=...,
    protocol=...,
    ecutwfc=OPTIMAL,        # Fixed at best value from Phase 1
    kspacing=VARYING,       # Varying parameter
    magnetic_config=...,    # ✅ Now forwarded
    hubbard_config=...,     # ✅ Now forwarded (FIXED)
)
```

### Complete Example: Fe2O2 with AFM + DFT+U

```python
from ase import Atoms
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow

# Create Fe2O2 structure
atoms = Atoms('Fe2O2',
              positions=[[0, 0, 0], [1.5, 0, 0],
                        [0.75, 0.75, 0], [0.75, 0.75, 1.5]],
              cell=[3, 3, 3])

# Hubbard parameters (new format with orbital specs)
hubbard_config = {
    'qe_version': '7.2',
    'projector': 'atomic',
    'u': {
        'Fe-3d': 4.3,      # Fe 3d gets 4.3 eV correction
        'O-2p': 3.0,       # O 2p gets 3.0 eV correction
    }
}

# Magnetic configuration (ferrimagnetic)
magnetic_config = {
    'Fe': {'mag': [1.0, -1.0]},  # Opposite moments on two Fe
}

workflow = ConvergenceWorkflow(
    atoms=atoms,
    pseudopotentials={
        'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
        'O': 'O.pbe-n-rrkjus_psl.1.0.0.UPF'
    },
    precision='medium',
    magnetic_config=magnetic_config,
    hubbard_config=hubbard_config,
)

# Run convergence: 
# - Phase 1: varies ecutwfc from auto-detected min to 200 Ry
# - Phase 2: varies kspacing from 0.3 to 0.1 Å⁻¹
# All with AFM magnetism + DFT+U applied
results = workflow.run_convergence_study(
    label='fe2o2_afm_dfu',
    max_iterations=4
)
```

### How It Works

The convergence workflow implements independent 2-phase studies:

1. **Phase 1: ecutwfc Convergence**
   - Fixes kspacing at 0.3 Å⁻¹ (coarse k-points)
   - Incrementally increases ecutwfc until energy converges
   - All calculations with your specified AFM + DFT+U
   
2. **Phase 2: kspacing Convergence**
   - Uses optimal ecutwfc from Phase 1
   - Incrementally decreases kspacing (finer k-points) until stress converges
   - All calculations with same AFM + DFT+U

### Automatic Features

When you provide `hubbard_config`:
- ConvergenceWorkflow forwards it to every `CalculationWorkflow` instance
- Each SCF calculation applies both magnetic moments AND Hubbard U correction
- No manual intervention needed across iterations

### Related Auto-Detection Features

ConvergenceWorkflow also auto-detects:

1. **min_ecutwfc**: From UPF headers (max suggested value)
   ```python
   # Example: Gd pseudo has suggested 69 Ry → starts convergence from there
   workflow = ConvergenceWorkflow(
       atoms=gd_structure,
       pseudopotentials={'Gd': 'Gd.upf'}
   )
   # min_ecutwfc auto-detected as 69.0 Ry
   ```

2. **ecutrho_ratio**: From pseudo type (NC=4.0, US/PAW=8.0)
   ```python
   # Automatically uses correct ecutrho for energy stability
   ```

### See Also

- [docs/HUBBARD_PARAMETERS.md](docs/HUBBARD_PARAMETERS.md) - Format specifications
- [examples/convergence_with_hubbard.py](examples/convergence_with_hubbard.py) - Detailed examples
- [xespresso/hubbard.py](xespresso/hubbard.py) - HubbardConfig class
