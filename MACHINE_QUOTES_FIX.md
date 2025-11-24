# Machine Configuration Quotes Sanitization Fix

## Problem

Job files were sometimes created with unwanted quotes around shell commands, like:

```bash
#!/bin/bash
#SBATCH --job-name=espresso
...

"export OMP_NUM_THREADS=1"

module load quantum-espresso

"srun --mpi=pmi2" pw.x -in espresso.pwi > espresso.pwo
```

Instead of the correct format:

```bash
#!/bin/bash
#SBATCH --job-name=espresso
...

export OMP_NUM_THREADS=1

module load quantum-espresso

srun --mpi=pmi2 pw.x -in espresso.pwi > espresso.pwo
```

## Root Cause

The issue occurred when machine configuration JSON files contained embedded quotes in string values:

```json
{
  "machines": {
    "my_cluster": {
      "launcher": "\"srun --mpi=pmi2\"",
      "prepend": ["\"export OMP_NUM_THREADS=1\""]
    }
  }
}
```

This could happen when:
1. Users manually edit JSON configuration files and add quotes
2. Configuration values are improperly serialized/deserialized

## Solution

Added automatic sanitization of machine configuration values when loading from JSON files:

### New Functions in `xespresso/machines/config/loader.py`:

1. **`sanitize_string_value(value)`** - Removes leading and trailing quotes from string values
2. **`sanitize_list_values(items)`** - Sanitizes all strings in a list
3. **`sanitize_machine_config(config)`** - Sanitizes all relevant fields in a machine config dictionary

### Modified Functions:

1. **`_load_from_machines_json()`** - Now sanitizes config before creating Machine object
2. **`Machine.from_file()`** - Now sanitizes config loaded from individual JSON files

## Usage

No changes are required in user code. The sanitization happens automatically when loading machine configurations:

```python
from xespresso.machines.config.loader import load_machine

# This will automatically sanitize any embedded quotes
machine = load_machine("~/.xespresso/machines.json", "my_cluster", return_object=True)

# The machine configuration will have clean values without embedded quotes
print(machine.launcher)  # "srun --mpi=pmi2" (not "\"srun --mpi=pmi2\"")
```

## Testing

Comprehensive test suite added in `tests/test_machine_quotes_sanitization.py`:

- `test_sanitize_string_value()` - Tests string sanitization
- `test_sanitize_list_values()` - Tests list sanitization
- `test_machine_config_sanitization()` - Tests full machine config sanitization
- `test_job_file_no_quotes()` - Tests that generated job files don't contain quotes

All tests pass successfully.

## Affected Fields

The following machine configuration fields are sanitized:

- `launcher` - MPI launcher command
- `prepend` - Commands to run before the job
- `postpend` - Commands to run after the job
- `modules` - Environment modules to load
- `workdir` - Working directory path
- `host` - Remote host (for SSH)
- `username` - SSH username
- `env_setup` - Environment setup commands

## Backward Compatibility

This fix is fully backward compatible:
- Configurations without embedded quotes work exactly as before
- Configurations with embedded quotes are now fixed automatically
- No API changes - sanitization is transparent to users

## Example

Before this fix, a config like this would generate a broken job file:

```json
{
  "machines": {
    "cluster": {
      "launcher": "\"srun --mpi=pmi2\"",
      "prepend": ["\"export OMP_NUM_THREADS=1\""]
    }
  }
}
```

After this fix, the same config is automatically sanitized and generates correct job files:

```bash
export OMP_NUM_THREADS=1
srun --mpi=pmi2 pw.x -in espresso.pwi > espresso.pwo
```
