# Environment Setup (env_setup) for Module Detection

## Problem

When using SSH to detect Quantum ESPRESSO codes on remote systems, the `module` command is often not available in non-interactive SSH sessions. This is because:

1. Non-interactive SSH sessions don't automatically source user profile files (`.bashrc`, `.bash_profile`, etc.)
2. The `module` command is typically defined in these profile files or in system-wide files like `/etc/profile.d/modules.sh`
3. Without the module command, code detection fails or cannot load the necessary QE modules

## Solution

The `env_setup` parameter has been added throughout the codebase to allow sourcing environment setup files before running commands that depend on the module system.

## Usage

### 1. Using env_setup with detect_qe_codes()

```python
from xespresso.codes import detect_qe_codes

# Detect codes with environment setup
config = detect_qe_codes(
    machine_name="my_cluster",
    modules=["quantum-espresso/7.2"],
    env_setup="source /etc/profile",
    ssh_connection={
        'host': 'cluster.edu',
        'username': 'myuser',
        'port': 22
    }
)
```

### 2. Using env_setup in Machine Configuration

```python
from xespresso.machines.machine import Machine

# Create a machine with env_setup
machine = Machine(
    name="remote_cluster",
    execution="remote",
    host="cluster.edu",
    username="myuser",
    auth={"method": "key", "ssh_key": "~/.ssh/id_rsa"},
    use_modules=True,
    modules=["quantum-espresso/7.2"],
    env_setup="source /etc/profile && source ~/.bashrc"
)

# Convert to queue configuration
queue = machine.to_queue()
# queue now contains 'env_setup' key
```

### 3. Machine Configuration File with env_setup

Create or update a machine configuration file (e.g., `~/.xespresso/machines/my_cluster.json`):

```json
{
  "name": "my_cluster",
  "execution": "remote",
  "host": "cluster.edu",
  "username": "myuser",
  "port": 22,
  "auth": {
    "method": "key",
    "ssh_key": "~/.ssh/id_rsa"
  },
  "use_modules": true,
  "modules": ["quantum-espresso/7.2"],
  "env_setup": "source /etc/profile.d/modules.sh",
  "scheduler": "slurm",
  "workdir": "/scratch/myuser"
}
```

When you use this machine for code detection, the `env_setup` will be automatically used:

```python
from xespresso.codes import detect_qe_codes

# Automatically loads env_setup from machine config
config = detect_qe_codes(machine_name="my_cluster")
```

### 4. Using env_setup for Job Submission

The `env_setup` is also used when submitting remote jobs:

```python
from ase.build import bulk
from xespresso import Espresso
from xespresso.machines.config.loader import load_machine

# Load machine with env_setup
machine = load_machine("my_cluster", return_object=True)
queue = machine.to_queue()

# Create calculator
atoms = bulk('Si')
calc = Espresso(
    input_data={
        'calculation': 'scf',
        'ecutwfc': 30,
    },
    pseudopotentials={'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'},
    queue=queue
)

atoms.calc = calc
energy = atoms.get_potential_energy()
# The job submission will use env_setup from queue
```

## Common env_setup Values

### For systems with environment modules:
```bash
# Source the modules initialization file
"source /etc/profile.d/modules.sh"

# Or source the full profile
"source /etc/profile"

# Multiple sources
"source /etc/profile && source ~/.bashrc"
```

### For systems without modules:
```bash
# Just source user profile
"source ~/.bashrc"

# Source a custom QE environment
"source /opt/qe-7.2/env.sh"

# Multiple environment files
"source /opt/intel/bin/compilervars.sh intel64 && source /opt/qe/env.sh"
```

## Automatic Extraction from Machine Config

When using `detect_qe_codes()` with `auto_load_machine=True` (default), the function will:

1. First check if `env_setup` is explicitly provided as a parameter
2. If not, try to load it from the machine's `env_setup` field
3. If still not found, check if the machine's `prepend` field contains environment setup commands
4. Use that for code detection

This means you can define `env_setup` once in your machine configuration and it will be used automatically.

## Default Behavior

- **For code detection**: If no `env_setup` is provided, the system will try to detect codes without it
- **For remote job submission**: If no `env_setup` is provided in the queue, it defaults to `"source /etc/profile"`

## API Reference

### CodesManager.detect_codes()

```python
@classmethod
def detect_codes(cls, 
                 search_paths: Optional[List[str]] = None,
                 qe_prefix: Optional[str] = None,
                 modules: Optional[List[str]] = None,
                 ssh_connection: Optional[Dict] = None,
                 use_modules: bool = True,
                 env_setup: Optional[str] = None) -> Dict[str, str]:
    """
    Args:
        env_setup: Shell commands to set up environment before detection.
                  Useful for sourcing profile files to make module command available.
                  Example: "source /etc/profile" or "source ~/.bashrc"
    """
```

### detect_qe_codes()

```python
def detect_qe_codes(machine_name: str = "local",
                   qe_prefix: Optional[str] = None,
                   search_paths: Optional[List[str]] = None,
                   modules: Optional[List[str]] = None,
                   ssh_connection: Optional[Dict] = None,
                   env_setup: Optional[str] = None,
                   auto_load_machine: bool = True) -> CodesConfig:
    """
    Args:
        env_setup: Shell commands to set up environment before detection.
                  Example: "source /etc/profile" or "source ~/.bashrc"
    """
```

### Machine.__init__()

```python
def __init__(self,
            name: str,
            # ... other parameters ...
            env_setup: Optional[str] = None,
            **kwargs):
    """
    Args:
        env_setup: Environment setup commands for SSH sessions.
                  Example: "source /etc/profile" or "source ~/.bashrc"
                  Useful for making module command available in non-interactive SSH.
    """
```

## Benefits

1. **Flexible**: Works with any environment setup commands
2. **Automatic**: When using machine configs, env_setup is applied automatically
3. **Consistent**: Same env_setup is used for both code detection and job submission
4. **Optional**: If you don't need it, you don't have to use it
5. **Compatible**: Maintains backward compatibility with existing configurations

## Troubleshooting

### Module command not found during detection

**Problem**: Getting "module: command not found" when detecting codes remotely

**Solution**: Add `env_setup` to source the modules initialization:

```python
config = detect_qe_codes(
    machine_name="cluster",
    modules=["quantum-espresso"],
    env_setup="source /etc/profile.d/modules.sh"
)
```

### Different environment for detection vs. execution

**Problem**: Codes are detected but jobs fail with "module: command not found"

**Solution**: Make sure `env_setup` is included in your machine configuration so it's used for both detection and execution:

```json
{
  "name": "cluster",
  "env_setup": "source /etc/profile.d/modules.sh",
  ...
}
```

### Finding the right env_setup command

To find what to use for `env_setup`, SSH into your cluster and check:

```bash
# Check what's in your profile
cat ~/.bash_profile
cat ~/.bashrc
cat /etc/profile

# Check for modules initialization
ls /etc/profile.d/modules.sh
ls /usr/share/Modules/init/bash

# Test if module works after sourcing
source /etc/profile
module avail
```

Then use the working command in your `env_setup`.
