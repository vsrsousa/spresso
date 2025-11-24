# Machine Configuration Examples

This directory contains example machine configurations for xespresso.

## Files

### Single File Configuration

- **machines.json** - Traditional format with all machines in one file

### Modular Configuration (machines/ directory)

- **machines/local_desktop.json** - Local execution configuration
- **machines/slurm_cluster.json** - Remote SLURM cluster configuration
- **machines/gpu_cluster.json** - GPU-enabled SLURM cluster configuration
- **machines/default.json** - Default machine specification

## Usage

### Option 1: Copy to Home Directory

Copy the desired format to your home directory:

```bash
# For single file configuration
cp machines.json ~/.xespresso/

# For modular configuration
mkdir -p ~/.xespresso/machines
cp machines/*.json ~/.xespresso/machines/
```

### Option 2: Use Directly

```python
from xespresso.machines import load_machine

# Load from custom location
queue = load_machine(
    config_path="/path/to/machines.json",
    machine_name="slurm_cluster"
)

# Or for modular config
queue = load_machine(
    machines_dir="/path/to/machines",
    machine_name="slurm_cluster"
)
```

## Customization

Before using these examples, customize the following fields:

### For Local Machines
- `workdir`: Set to your preferred calculation directory

### For Remote Machines
- `host`: Your cluster hostname
- `username`: Your SSH username
- `workdir`: Your remote working directory
- `auth.ssh_key`: Path to your SSH private key
- `modules`: Available modules on your system
- `resources`: Adjust for your cluster's requirements

## Examples by Use Case

### Local Testing
Use `local_desktop.json` for:
- Development and testing
- Small calculations
- Systems without job schedulers

### Standard HPC Cluster
Use `slurm_cluster.json` for:
- CPU-based calculations
- Standard compute nodes
- Typical production runs

### GPU Acceleration
Use `gpu_cluster.json` for:
- GPU-enabled Quantum ESPRESSO
- Accelerated DFT calculations
- Specialized GPU partitions

## Testing Your Configuration

```python
from xespresso.machines import Machine

# Test loading
machine = Machine.from_file("machines/slurm_cluster.json")
print(f"✓ Loaded: {machine}")

# Validate
try:
    queue = machine.to_queue()
    print("✓ Configuration is valid")
except Exception as e:
    print(f"✗ Configuration error: {e}")
```

## More Information

See the documentation:
- [Machine Configuration Guide](../docs/MACHINE_CONFIGURATION.md)
- [Remote Connection Persistence](../docs/REMOTE_CONNECTION_PERSISTENCE.md)
