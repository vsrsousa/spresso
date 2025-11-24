# Machine Configuration Guide

## Overview

xespresso supports two ways to organize machine configurations:

1. **Traditional**: Single `machines.json` file with all machines
2. **Modular**: Individual JSON files for each machine in `machines/` directory

Both approaches can be used simultaneously, and the system will look for machines in both places.

**Need to migrate?** See the [Machine Migration Guide](MACHINE_MIGRATION.md) for step-by-step instructions on converting from the traditional format to the modular format.

## Directory Structure

```
~/.xespresso/
├── machines.json          # Optional: Traditional single-file config
└── machines/              # Optional: Modular per-machine configs
    ├── default.json       # Optional: Specifies default machine
    ├── local_desktop.json
    ├── cluster1.json
    └── cluster2.json
```

## Configuration Formats

### Traditional Single File (machines.json)

```json
{
  "default": "cluster1",
  "machines": {
    "local_desktop": {
      "execution": "local",
      "scheduler": "direct",
      "workdir": "./calculations",
      "nprocs": 4,
      "launcher": "mpirun -np {nprocs}"
    },
    "cluster1": {
      "execution": "remote",
      "scheduler": "slurm",
      "workdir": "/home/user/calculations",
      "host": "cluster1.university.edu",
      "username": "user",
      "port": 22,
      "auth": {
        "method": "key",
        "ssh_key": "~/.ssh/id_rsa"
      },
      "nprocs": 16,
      "launcher": "mpirun -np {nprocs}",
      "use_modules": true,
      "modules": ["quantum-espresso/7.0"],
      "resources": {
        "nodes": 1,
        "ntasks-per-node": 16,
        "time": "02:00:00",
        "partition": "compute"
      }
    }
  }
}
```

### Modular Individual Files

**~/.xespresso/machines/local_desktop.json:**
```json
{
  "name": "local_desktop",
  "execution": "local",
  "scheduler": "direct",
  "workdir": "./calculations",
  "nprocs": 4,
  "launcher": "mpirun -np {nprocs}"
}
```

**~/.xespresso/machines/cluster1.json:**
```json
{
  "name": "cluster1",
  "execution": "remote",
  "scheduler": "slurm",
  "workdir": "/home/user/calculations",
  "host": "cluster1.university.edu",
  "username": "user",
  "port": 22,
  "auth": {
    "method": "key",
    "ssh_key": "~/.ssh/id_rsa"
  },
  "nprocs": 16,
  "launcher": "mpirun -np {nprocs}",
  "use_modules": true,
  "modules": ["quantum-espresso/7.0"],
  "resources": {
    "nodes": 1,
    "ntasks-per-node": 16,
    "time": "02:00:00",
    "partition": "compute"
  }
}
```

**~/.xespresso/machines/default.json:**
```json
{
  "default": "cluster1"
}
```

## Configuration Fields

### Required Fields (All Machines)

- `execution` (string): `"local"` or `"remote"`
- `scheduler` (string): `"direct"` or `"slurm"`
- `workdir` (string): Working directory path

### Optional Fields (All Machines)

- `nprocs` (int): Number of processors (default: 1)
- `launcher` (string): MPI launcher command template (default: `"mpirun -np {nprocs}"`)
- `use_modules` (bool): Enable environment modules (default: false)
- `modules` (list): Modules to load (default: [])
- `prepend` (string or list): Commands to run before job
- `postpend` (string or list): Commands to run after job
- `resources` (object): Scheduler-specific resource requirements

### Required Fields (Remote Machines)

- `host` (string): Remote hostname or IP
- `username` (string): SSH username
- `auth` (object): Authentication configuration
  - `method` (string): Must be `"key"` (password auth not supported)
  - `ssh_key` (string): Path to SSH private key (default: `"~/.ssh/id_rsa"`)
  - `port` (int): SSH port (default: 22)

### Optional Fields (Remote Machines)

- `port` (int): SSH port (default: 22)

### SLURM-specific Fields (resources)

When using `"scheduler": "slurm"`, you can specify these in `resources`:

- `nodes` (int): Number of nodes
- `ntasks-per-node` (int): Tasks per node
- `time` (string): Wall time limit (format: HH:MM:SS)
- `partition` (string): SLURM partition name
- `mem` (string): Memory per node (e.g., "32G")
- `gres` (string): Generic resources (e.g., "gpu:2")

## Usage Examples

### Using Traditional Config

```python
from xespresso.machines import load_machine

# Load default machine
queue = load_machine()

# Load specific machine
queue = load_machine(machine_name="cluster1")

# Load as Machine object
machine = load_machine(machine_name="cluster1", return_object=True)
print(machine.name, machine.host, machine.is_remote)
```

### Using Modular Config

```python
from xespresso.machines import load_machine

# Automatically finds cluster1.json in machines/ directory
queue = load_machine(machine_name="cluster1")

# Or specify custom paths
queue = load_machine(
    config_path="~/.xespresso/machines.json",
    machines_dir="~/.xespresso/machines",
    machine_name="cluster1"
)
```

### Creating Machines Programmatically

```python
from xespresso.machines import Machine

# Create a local machine
local = Machine(
    name="my_local",
    execution="local",
    scheduler="direct",
    workdir="/tmp/calculations",
    nprocs=8
)

# Save to file
local.to_file("~/.xespresso/machines/my_local.json")

# Create a remote machine
remote = Machine(
    name="my_cluster",
    execution="remote",
    scheduler="slurm",
    workdir="/home/user/calc",
    host="cluster.edu",
    username="user",
    auth={"method": "key", "ssh_key": "~/.ssh/id_rsa"},
    nprocs=32,
    resources={
        "nodes": 2,
        "ntasks-per-node": 16,
        "time": "04:00:00"
    }
)

# Save and use
remote.to_file("~/.xespresso/machines/my_cluster.json")
queue = remote.to_queue()
```

## Advantages of Modular Configuration

### 1. Better Organization
- Each machine has its own file
- Easier to find and edit specific machines
- Cleaner version control (one file per machine)

### 2. Easier Sharing
- Share individual machine configs between users
- Template files for common setups
- No need to merge entire config files

### 3. Flexible Defaults
- Set default machine per project
- Override in specific scripts
- Environment-specific defaults

### 4. Reduced Conflicts
- Multiple people can add machines without conflicts
- Git merge conflicts are minimized
- Each file is independent

## Best Practices

### 1. Use Modular Config for Teams
```
~/.xespresso/machines/
├── default.json          # Team default
├── local_debug.json      # For local testing
├── dev_cluster.json      # Development cluster
├── prod_cluster.json     # Production cluster
└── gpu_node.json         # Special GPU configuration
```

### 2. Use Single File for Personal Use
```json
{
  "default": "local_desktop",
  "machines": {
    "local_desktop": {...},
    "university_cluster": {...}
  }
}
```

### 3. Mix Both Approaches
- Use `machines.json` for stable/core machines
- Use individual files for experimental/temporary machines
- System looks in both places automatically

### 4. Set Appropriate Defaults
```json
// In machines.json
{
  "default": "local_desktop",  // Safe default for development
  "machines": {...}
}

// In machines/default.json (for production environment)
{
  "default": "prod_cluster"
}
```

### 5. Use Machine Class for Validation
```python
from xespresso.machines import Machine

# Load and validate
machine = Machine.from_file("machines/cluster1.json")

# Machine class validates on creation
# Will raise ValueError if configuration is invalid
```

## Migration from Old Format

### Old Format (Deprecated)
```python
queue = {
    "execution": "remote",
    "remote_host": "cluster.edu",
    "remote_user": "user",
    # ... many fields ...
}
```

### New Format (Recommended)
```python
# Option 1: Load from config
queue = load_machine("cluster1")

# Option 2: Use Machine class
machine = Machine(
    name="cluster1",
    execution="remote",
    host="cluster.edu",
    username="user",
    ...
)
queue = machine.to_queue()
```

### Conversion Script
```python
from xespresso.machines import Machine

# Convert old queue dict to Machine and save
old_queue = {...}  # Your old config

machine = Machine(
    name="my_machine",
    execution=old_queue.get("execution"),
    scheduler=old_queue.get("scheduler", "direct"),
    workdir=old_queue.get("local_dir" if old_queue["execution"] == "local" else "remote_dir"),
    host=old_queue.get("remote_host"),
    username=old_queue.get("remote_user"),
    auth=old_queue.get("remote_auth", {}),
    # ... copy other fields ...
)

machine.to_file("~/.xespresso/machines/my_machine.json")
```

## Troubleshooting

### Machine Not Found
```python
# List all available machines
from xespresso.machines.config import list_machines
machines = list_machines()
print("Available machines:", machines)
```

### Check Configuration Paths
```python
import os
config_path = os.path.expanduser("~/.xespresso/machines.json")
machines_dir = os.path.expanduser("~/.xespresso/machines")

print("Config file exists:", os.path.exists(config_path))
print("Machines dir exists:", os.path.exists(machines_dir))
```

### Validate Machine Config
```python
from xespresso.machines import Machine

# This will raise errors if config is invalid
try:
    machine = Machine.from_file("machines/my_machine.json")
    print("✓ Configuration is valid")
except ValueError as e:
    print(f"✗ Configuration error: {e}")
```

## Summary

- ✅ Two configuration formats supported (single file or modular)
- ✅ Both can be used simultaneously
- ✅ Default machine can be specified
- ✅ Machine class provides validation and type safety
- ✅ Backward compatible with existing code
- ✅ Connection persistence works with both formats
