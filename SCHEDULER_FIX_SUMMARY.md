# Scheduler Configuration Fix Summary

## Problem Statement

User encountered error when running calculations through the Qt GUI:
```
Calculator "espresso" failed with command "sbatch job_file" failed in 
/scratch/projects/spresso/AlGdNi4/scf with error code 127
```

Error code 127 means "command not found" - `sbatch` was not available on the system.

## Root Cause

The issue had multiple contributing factors:

### 1. Architecture Issue - Calculator Creation
The original code was creating Espresso calculator objects **twice**:
- Once in `_generate_files()` (dry run)
- Again in `_run_calculation()` (job submission)

This violated xespresso/ASE patterns which expect:
```python
# Create calculator ONCE
calc = Espresso(parameters)

# Attach to atoms
atoms.calc = calc

# Use for dry run
calc.write_input(atoms)

# Use for job execution (reuses same calc)
energy = atoms.get_potential_energy()
```

### 2. Machine Configuration Issue
The user's machine configuration likely had:
- `scheduler: 'slurm'` (expecting SLURM)
- `execution: 'local'` (trying to run locally)
- But `sbatch` was not installed on the local system

The error path `/scratch/projects/spresso/` suggests this should have been configured for **remote** execution instead.

## Solution Implemented

### 1. Proper Calculator Management

Created `_prepare_calculator()` method that:
- Creates Espresso calculator once following xespresso patterns
- Stores calculator and prepared_atoms in session_state
- Reuses the same calculator for both dry run and job submission

```python
def _prepare_calculator(self, atoms, config, label):
    """Create Espresso calculator once and store for reuse"""
    # Convert kspacing to kpts if needed
    # Build queue from machine.to_queue()
    # Create calculator via prepare_calculation_from_gui()
    # Store in session_state for reuse
    return prepared_atoms, calc

def _generate_files(self):
    """Dry run - uses stored calculator"""
    prepared_atoms, calc = self._prepare_calculator(atoms, config, label)
    calc.write_input(prepared_atoms)

def _run_calculation(self):
    """Job run - reuses stored calculator or creates new one"""
    calc = self.session_state.get('espresso_calculator')
    if calc is None or calc.label != label:
        prepared_atoms, calc = self._prepare_calculator(atoms, config, label)
    energy = prepared_atoms.get_potential_energy()
```

### 2. Scheduler Validation

Added `validate_scheduler_availability()` function that:
- Checks if `sbatch` is available when `scheduler: 'slurm'` and `execution: 'local'`
- Does NOT check for remote execution (sbatch runs on remote machine)
- Does NOT modify user's configuration
- Provides clear, actionable error messages with 3 solutions

```python
def validate_scheduler_availability(queue_config):
    """Validate scheduler without modifying user's config"""
    scheduler = queue_config.get('scheduler', '').lower()
    execution = queue_config.get('execution', 'local').lower()
    
    # Only validate SLURM for local execution
    if scheduler == 'slurm' and execution == 'local':
        if shutil.which("sbatch") is None:
            return False, detailed_error_message
    
    return True, None
```

### 3. User-Friendly Error Messages

When validation fails, users see:
```
❌ SLURM Scheduler Not Available

Your machine configuration uses the SLURM scheduler, but the 'sbatch'
command is not found on this system.

To fix this issue, choose ONE of the following:

1. Install SLURM on your system:
   - Linux: sudo apt-get install slurm-wlm (Ubuntu/Debian)
   - Or check your distribution's package manager

2. Change your machine configuration:
   - Go to Machine Configuration page
   - Select or create a machine with 'direct' scheduler
   - Save and select that machine for this calculation

3. For remote execution:
   - Set execution to 'remote' in your machine configuration
   - Configure SSH access to a cluster with SLURM
```

### 4. Logging for Debugging

Added logging to show machine queue configuration:
```python
logging.info(f"Machine '{machine.name}' queue config: {calc_config['queue']}")
```

This helps identify configuration issues during development and debugging.

## How xespresso Remote Execution Works

When `execution: 'remote'` is set in machine configuration:

1. Machine's `to_queue()` returns queue dict with `execution: 'remote'`
2. Calculator created with this queue configuration
3. When `write_input()` is called, scheduler is set up via `set_queue(calc)`
4. Scheduler writes job script locally
5. When `atoms.get_potential_energy()` is called:
   - ASE calls calculator's `execute()` method
   - `execute()` calls `scheduler.run()`
   - Since `execution: 'remote'`, RemoteExecutionMixin's `run()` method:
     - Establishes SSH connection to remote host
     - Transfers input files and job script to remote directory
     - Executes `sbatch job_file` ON THE REMOTE MACHINE
     - Waits for job completion (polls squeue)
     - Retrieves output files back to local machine
     - ASE parses results

## User Action Required

Users experiencing this error should check their machine configuration:

**For Remote Clusters (like the user's case):**
```yaml
name: my_cluster
execution: remote      # ← Must be 'remote'
scheduler: slurm
host: cluster.university.edu
username: myuser
workdir: /scratch/projects/spresso
```

**For Local Execution:**
```yaml
name: my_local_machine
execution: local
scheduler: direct      # ← Use 'direct', not 'slurm'
workdir: ~/calculations
```

**For Local with SLURM:**
```yaml
name: my_workstation
execution: local
scheduler: slurm       # ← OK if sbatch is installed
workdir: ~/calculations
```

## Testing

Created unit tests to verify validation logic:
- SLURM + local + no sbatch → Invalid (shows error)
- SLURM + remote → Valid (no local validation)
- direct + local → Valid (no sbatch needed)

All tests passed successfully.

## Benefits

1. **Correct Architecture**: Calculator created once, reused properly
2. **Respects User Config**: No silent overrides of machine configuration
3. **Clear Error Messages**: Users know exactly how to fix the issue
4. **Remote Execution Works**: Properly uses SSH and remote sbatch
5. **Logging for Debug**: Easy to diagnose configuration issues
6. **Follows xespresso Patterns**: Uses proper ASE/xespresso workflow

## Files Modified

- `qtgui/pages/job_submission.py`:
  - Added `validate_scheduler_availability()` function
  - Added `_prepare_calculator()` method
  - Refactored `_generate_files()` to use `_prepare_calculator()`
  - Refactored `_run_calculation()` to reuse calculator from session state
  - Added logging for machine configuration debugging
  - Improved error messages with actionable steps
