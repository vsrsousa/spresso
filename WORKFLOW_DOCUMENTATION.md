# Simplified Workflow for Quantum ESPRESSO Calculations

This document describes the new simplified workflow system in xespresso that makes it easy to run common calculations with protocol presets and k-spacing support.

## Features

- **Protocol Presets**: Choose between `fast`, `moderate`, and `accurate` calculations
- **K-spacing Support**: Use k-point spacing instead of explicit k-meshes
- **CIF File Support**: Create workflows directly from CIF files
- **Pseudopotential Management**: Store and manage pseudopotential configurations in JSON format
- **Simple API**: Quick functions for common calculation types

## Installation

The workflow functionality is included in spresso. No additional installation is required.

```bash
pip install spresso
```

**Note:** The PyPI package is named `spresso`, but the Python module is `xespresso` (for backwards compatibility).
Install with `pip install spresso`, but import as `from xespresso import ...`

## Quick Start

### Basic SCF Calculation

```python
from xespresso import quick_scf

# Quick SCF calculation from CIF file
calc = quick_scf(
    'structure.cif',
    {'Fe': 'Fe.pbe-spn.UPF'},
    protocol='moderate',
    label='scf/fe'
)

energy = calc.results['energy']
print(f"Energy: {energy} eV")
```

### Structure Relaxation

```python
from xespresso import quick_relax

# Quick relaxation
calc = quick_relax(
    'structure.cif',
    {'Fe': 'Fe.pbe-spn.UPF'},
    protocol='moderate',
    relax_type='vc-relax',  # relax both cell and ions
    label='relax/fe'
)

relaxed_atoms = calc.results['atoms']
```

## Protocol Presets

Three protocol presets are available:

### Fast
- **ecutwfc**: 30.0 Ry
- **ecutrho**: 240.0 Ry
- **conv_thr**: 1.0e-6
- **kspacing**: 0.5 Å⁻¹
- Best for: Quick tests, structure screening

### Moderate
- **ecutwfc**: 50.0 Ry
- **ecutrho**: 400.0 Ry
- **conv_thr**: 1.0e-8
- **kspacing**: 0.3 Å⁻¹
- Best for: Standard calculations, most production runs

### Accurate
- **ecutwfc**: 80.0 Ry
- **ecutrho**: 640.0 Ry
- **conv_thr**: 1.0e-10
- **kspacing**: 0.15 Å⁻¹
- Best for: High-accuracy results, publication-quality data

## Using the Workflow Class

For more control, use the `CalculationWorkflow` class:

```python
from xespresso import CalculationWorkflow
from ase.io import read

# Load structure
atoms = read('structure.cif')

# Create workflow
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={'Si': 'Si.pbe.UPF'},
    protocol='moderate',
    kspacing=0.3  # Optional: override preset k-spacing
)

# Run SCF
calc = workflow.run_scf(label='scf/silicon')

# Run relaxation
calc = workflow.run_relax(label='relax/silicon', relax_type='relax')
```

## K-spacing Support

Instead of specifying explicit k-points, you can use k-spacing (in Å⁻¹):

```python
from xespresso import CalculationWorkflow
import numpy as np

# Using k-spacing in physical units (Angstrom^-1)
# The workflow automatically handles the 2π normalization internally
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials=pseudopotentials,
    protocol='moderate',
    kspacing=0.20  # Just pass the physical value - no need for /(2*np.pi)!
)

# See what k-points this corresponds to
kpts = workflow._get_kpts()
print(f"K-points: {kpts}")
```

The workflow uses ASE's `kspacing_to_grid` function internally to convert k-spacing to k-points. 
**You don't need to worry about the 2π normalization** - just pass your desired k-spacing in Angstrom^-1.

### Using kpts_from_spacing Outside Workflows

For custom calculation setups outside of workflows, you can use the `kpts_from_spacing` utility function:

```python
from xespresso import Espresso, kpts_from_spacing
from ase.build import bulk

atoms = bulk('Si', cubic=True)

# Convert k-spacing to k-points without manual normalization
kpts = kpts_from_spacing(atoms, 0.20)  # Clean and simple!

# Use in a calculation
calc = Espresso(
    pseudopotentials={'Si': 'Si.pbe.UPF'},
    ecutwfc=50,
    kpts=kpts,
    label='scf/silicon'
)

atoms.calc = calc
energy = atoms.get_potential_energy()
```

This is equivalent to manually calling `kspacing_to_grid(atoms, 0.20/(2*np.pi))` but with a cleaner API.

## Magnetic and Hubbard Support

The workflow seamlessly integrates with xespresso's magnetic configuration and Hubbard parameter functionality.

### Simple Magnetic Configurations

```python
from xespresso import CalculationWorkflow

# Ferromagnetic configuration
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
    protocol='moderate',
    magnetic_config='ferro'  # or 'ferromagnetic'
)

# Antiferromagnetic configuration
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
    protocol='moderate',
    magnetic_config='antiferro'  # or 'antiferromagnetic'
)
```

### Element-based Magnetic Configurations

```python
# Define magnetization per element
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF', 'O': 'O.pbe.UPF'},
    protocol='moderate',
    magnetic_config={
        'Fe': [1, -1],  # Two non-equivalent Fe atoms (AFM)
        'O': [0]        # Non-magnetic oxygen
    }
)

# With automatic cell expansion if needed
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
    protocol='moderate',
    magnetic_config={'Fe': [1, 1, -1, -1]},  # Need 4 Fe but only have 2
    expand_cell=True  # Automatically expands the cell
)
```

### Magnetic Configurations with Hubbard Parameters

```python
# Old format (QE < 7.0) with Hubbard U
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF', 'O': 'O.pbe.UPF'},
    protocol='accurate',
    magnetic_config={
        'Fe': {'mag': [1, -1], 'U': 4.3}  # AFM with Hubbard U
    }
)

# New format (QE >= 7.0) with HUBBARD card
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF', 'O': 'O.pbe.UPF'},
    protocol='accurate',
    magnetic_config={
        'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}},  # U on Fe-3d orbital
        'O': {'mag': [0]}
    },
    input_data={'qe_version': '7.2'}
)

# With inter-site Hubbard V parameters
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF', 'O': 'O.pbe.UPF'},
    protocol='accurate',
    magnetic_config={
        'Fe': {
            'mag': [1],
            'U': {'3d': 4.3},
            'V': [{'species2': 'O', 'orbital2': '2p', 'value': 1.0}]
        },
        'O': {'mag': [0]}
    },
    input_data={'qe_version': '7.2'}
)
```

### Using with Quick Functions

```python
from xespresso import quick_scf, quick_relax

# Quick SCF with antiferromagnetic configuration
calc = quick_scf(
    'structure.cif',
    {'Fe': 'Fe.pbe-spn.UPF'},
    magnetic_config='antiferro',
    protocol='moderate'
)

# Quick relaxation with Hubbard parameters
calc = quick_relax(
    atoms,
    {'Fe': 'Fe.pbe-spn.UPF', 'O': 'O.pbe.UPF'},
    magnetic_config={'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}},
    protocol='accurate',
    relax_type='vc-relax'
)
```

See the [Magnetic Helpers documentation](MAGNETIC_HELPERS.md) for more details on magnetic configurations.

## Pseudopotential Configuration Management

Store and manage pseudopotential configurations in `~/.xespresso`:

### Saving a Configuration

```python
from xespresso.utils import save_pseudo_config

config = {
    "name": "PBE_efficiency",
    "description": "Efficient PBE pseudopotentials",
    "functional": "PBE",
    "pseudopotentials": {
        "H": "H.pbe-rrkjus_psl.1.0.0.UPF",
        "C": "C.pbe-n-kjpaw_psl.1.0.0.UPF",
        "O": "O.pbe-n-kjpaw_psl.0.1.UPF",
        "Fe": "Fe.pbe-spn-kjpaw_psl.0.2.1.UPF",
    }
}

save_pseudo_config("pbe_efficiency", config)
```

### Loading and Using a Configuration

```python
from xespresso.utils import load_pseudo_config
from xespresso import quick_scf

# Load configuration
config = load_pseudo_config("pbe_efficiency")

# Use in calculation
calc = quick_scf(
    'structure.cif',
    config['pseudopotentials'],
    protocol='moderate'
)
```

### Managing Configurations

```python
from xespresso.utils import (
    list_pseudo_configs,
    delete_pseudo_config,
    get_pseudo_info
)

# List all configurations
configs = list_pseudo_configs()
print(f"Available: {configs}")

# Get pseudopotential for a specific element
pseudo = get_pseudo_info("pbe_efficiency", "Fe")
print(f"Fe pseudopotential: {pseudo}")

# Delete a configuration
delete_pseudo_config("old_config")
```

## Remote Job Execution

The workflow system fully supports remote job execution on HPC clusters with SLURM or direct execution. This allows you to run calculations on remote machines with the same simple API as local execution.

### Quick Start with Remote Execution

```python
from xespresso import quick_scf, CalculationWorkflow

# Option 1: Use a machine configuration
calc = quick_scf(
    'structure.cif',
    {'Fe': 'Fe.pbe-spn.UPF'},
    protocol='moderate',
    machine='cluster1',  # Load from ~/.xespresso/machines/cluster1.json
    label='scf/fe-remote'
)

# Option 2: Direct queue configuration
queue = {
    "execution": "remote",
    "scheduler": "slurm",
    "remote_host": "cluster.university.edu",
    "remote_user": "username",
    "remote_dir": "/home/username/calculations",
    "remote_auth": {
        "method": "key",
        "ssh_key": "~/.ssh/id_rsa"
    },
    "nodes": 1,
    "ntasks-per-node": 16,
    "time": "02:00:00",
    "partition": "compute"
}

calc = quick_scf(
    'structure.cif',
    {'Fe': 'Fe.pbe-spn.UPF'},
    protocol='moderate',
    queue=queue,
    label='scf/fe-remote'
)
```

### Machine Configuration Approach (Recommended)

The recommended way to use remote execution is to create machine configuration files:

**1. Create a machine configuration:**

```python
from xespresso.machines import Machine

machine = Machine(
    name="my_cluster",
    execution="remote",
    scheduler="slurm",
    workdir="/home/user/calculations",
    host="cluster.university.edu",
    username="user",
    auth={"method": "key", "ssh_key": "~/.ssh/id_rsa"},
    nprocs=32,
    resources={
        "nodes": 2,
        "ntasks-per-node": 16,
        "time": "04:00:00",
        "partition": "compute"
    }
)

# Save for reuse
machine.to_file("~/.xespresso/machines/my_cluster.json")
```

**2. Use in workflows:**

```python
from xespresso import CalculationWorkflow

workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
    protocol='moderate',
    machine='my_cluster'  # References the saved configuration
)

# Run calculations remotely
calc = workflow.run_scf(label='scf/fe')
calc = workflow.run_relax(label='relax/fe', relax_type='vc-relax')
```

### Direct Queue Configuration

For one-off calculations or testing, you can pass the queue configuration directly:

```python
from xespresso import CalculationWorkflow

queue = {
    "execution": "remote",
    "scheduler": "slurm",
    "remote_host": "cluster.edu",
    "remote_user": "username",
    "remote_dir": "/scratch/calculations",
    "remote_auth": {
        "method": "key",
        "ssh_key": "~/.ssh/id_rsa"
    },
    "nodes": 1,
    "ntasks-per-node": 16,
    "time": "02:00:00"
}

workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={'Si': 'Si.pbe.UPF'},
    protocol='moderate',
    queue=queue
)

calc = workflow.run_scf(label='scf/si')
```

### Remote Execution Features

#### Automatic File Transfer

The workflow automatically handles:
- **Input files**: Transferred to remote machine before calculation
- **Pseudopotentials**: Automatically located and transferred
- **Job scripts**: Generated and transferred based on scheduler
- **Output files**: Retrieved after calculation completes

#### Connection Persistence

SSH connections are automatically cached and reused:
- First calculation to a machine creates a new connection
- Subsequent calculations reuse the existing connection
- Connections are identified by `(hostname, username)`
- No manual connection management required

```python
# All three calculations use the same SSH connection
calc1 = quick_scf(atoms1, pseudos, machine='cluster1', label='scf/1')
calc2 = quick_scf(atoms2, pseudos, machine='cluster1', label='scf/2')
calc3 = quick_scf(atoms3, pseudos, machine='cluster1', label='scf/3')
```

#### SLURM Integration

For SLURM clusters, the workflow:
- Generates appropriate SLURM job scripts
- Submits jobs with `sbatch`
- Monitors job status with `squeue`
- Waits for completion before retrieving results
- Handles job failures gracefully

```python
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials=pseudos,
    protocol='accurate',
    machine='slurm_cluster'
)

# Job is submitted, monitored, and results retrieved automatically
calc = workflow.run_relax(label='relax/structure', relax_type='vc-relax')
```

### Combining Remote Execution with Other Features

All workflow features work seamlessly with remote execution:

**Remote + Magnetic Configurations:**
```python
calc = quick_scf(
    atoms,
    {'Fe': 'Fe.pbe-spn.UPF'},
    protocol='moderate',
    magnetic_config='antiferro',
    machine='cluster1'
)
```

**Remote + Hubbard Parameters:**
```python
calc = quick_scf(
    atoms,
    {'Fe': 'Fe.pbe-spn.UPF', 'O': 'O.pbe.UPF'},
    protocol='accurate',
    magnetic_config={'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}},
    machine='cluster1'
)
```

**Remote + Custom K-spacing:**
```python
calc = quick_relax(
    'structure.cif',
    {'Si': 'Si.pbe.UPF'},
    protocol='moderate',
    kspacing=0.2,
    machine='cluster1',
    relax_type='vc-relax'
)
```

### Local SLURM Execution

You can also use SLURM on your local machine (if available):

```python
local_queue = {
    "execution": "local",
    "scheduler": "slurm",
    "nodes": 1,
    "ntasks-per-node": 4,
    "time": "00:30:00"
}

workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials=pseudos,
    protocol='fast',
    queue=local_queue
)

calc = workflow.run_scf(label='scf/local-slurm')
```

### Configuration Validation

The workflow validates configurations to prevent common errors:

```python
# ERROR: Cannot specify both 'queue' and 'machine'
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials=pseudos,
    protocol='moderate',
    queue=queue_config,      # ❌ 
    machine='cluster1'       # ❌
)  # Raises ValueError

# CORRECT: Use one or the other
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials=pseudos,
    protocol='moderate',
    machine='cluster1'       # ✓
)
```

### Error Handling

The workflow provides clear error messages for common issues:

- **Missing pseudopotentials**: Lists which files are missing
- **Connection failures**: Reports authentication or network issues  
- **Job failures**: Reports SLURM job status and errors
- **File transfer errors**: Retries with exponential backoff

### Best Practices

1. **Use machine configurations** for regular clusters
2. **Test locally first** with 'fast' protocol preset
3. **Check pseudopotential paths** before remote submission
4. **Monitor first job** to ensure configuration is correct
5. **Clean up connections** at program exit (optional):
   ```python
   from xespresso.schedulers.remote_mixin import RemoteExecutionMixin
   RemoteExecutionMixin.close_all_connections()
   ```

### See Also

- [Machine Configuration Guide](docs/MACHINE_CONFIGURATION.md) - Creating machine configs
- [Remote Connection Persistence](docs/REMOTE_CONNECTION_PERSISTENCE.md) - Connection management
- [Examples](examples/workflow_remote_execution.py) - Complete remote execution examples

## Advanced Usage

### Custom Input Parameters

You can override any preset parameter or add custom ones:

```python
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials=pseudopotentials,
    protocol='moderate',
    input_data={
        'mixing_beta': 0.9,  # Override preset
        'nspin': 2,          # Add magnetic calculation
    }
)
```

### Working with Different Calculation Types

```python
# SCF calculation
calc = workflow.run_scf(label='scf/fe')

# Ion relaxation only
calc = workflow.run_relax(label='relax/fe', relax_type='relax')

# Full cell relaxation
calc = workflow.run_relax(label='relax/fe-vc', relax_type='vc-relax')
```

## Examples

Complete working examples are available in the `examples/` directory:

- `workflow_simple_example.py`: Basic workflow usage
- `pseudo_config_example.py`: Pseudopotential configuration management
- `complete_workflow_example.py`: Complete integration example from CIF to results

## API Reference

### CalculationWorkflow

Main class for managing calculations.

**Constructor Parameters:**
- `atoms`: ASE Atoms object
- `pseudopotentials`: Dict mapping elements to pseudopotential files
- `quality`: Protocol preset ('fast', 'moderate', 'accurate')
- `kspacing`: K-point spacing in Å⁻¹ (optional)
- `input_data`: Additional input parameters (optional)
- `magnetic_config`: Magnetic configuration (optional)
- `expand_cell`: Expand cell for magnetic config (optional)
- `queue`: Queue configuration dict for job submission (optional)
- `machine`: Name of machine configuration to load (optional)
- `**kwargs`: Additional parameters for Espresso calculator

**Methods:**
- `run_scf(label, **kwargs)`: Run SCF calculation
- `run_relax(label, relax_type='relax', **kwargs)`: Run relaxation
- `get_atoms()`: Get the atoms object
- `get_preset_info()`: Get preset information

**Class Methods:**
- `from_cif(cif_file, ...)`: Create workflow from CIF file

### Quick Functions

- `quick_scf(structure, pseudopotentials, protocol='moderate', ...)`: Quick SCF calculation
  - Additional parameters: `kspacing`, `magnetic_config`, `expand_cell`, `queue`, `machine`, `label`
- `quick_relax(structure, pseudopotentials, protocol='moderate', ...)`: Quick relaxation
  - Additional parameters: `kspacing`, `magnetic_config`, `expand_cell`, `queue`, `machine`, `label`, `relax_type`

### K-point Utility Functions

- `kpts_from_spacing(atoms, kspacing)`: Convert k-spacing to k-point grid without manual 2π normalization
  - **Parameters:**
    - `atoms`: ASE Atoms object
    - `kspacing`: K-point spacing in Å⁻¹ (physical units)
  - **Returns:** Tuple of k-points (kx, ky, kz)
  - **Example:** `kpts = kpts_from_spacing(atoms, 0.20)`

### Pseudo Configuration Functions

- `save_pseudo_config(name, config, overwrite=False)`: Save configuration
- `load_pseudo_config(name)`: Load configuration
- `list_pseudo_configs()`: List all configurations
- `delete_pseudo_config(name)`: Delete configuration
- `get_pseudo_info(config_name, element)`: Get pseudopotential for element

## Notes

- K-spacing is converted to k-points using `ase.io.espresso.kspacing_to_grid`
- Pseudopotential configurations are stored as JSON in `~/.xespresso/`
- All presets use sensible defaults that can be overridden
- The workflow integrates seamlessly with existing xespresso functionality

## See Also

- [ASE Espresso Calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/espresso.html)
- [Quantum ESPRESSO Documentation](https://www.quantum-espresso.org/)
