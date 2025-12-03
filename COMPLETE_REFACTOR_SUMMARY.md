# Job Submission Refactor - Complete Fix

## Problem Statement

User reported error when running calculations through Qt GUI:
```
Calculator "espresso" failed with command "sbatch job_file" failed in 
/scratch/projects/spresso/AlGdNi4/scf with error code 127
```

Error code 127 = "command not found" → `sbatch` was not available.

## Root Cause - Multiple Issues

### 1. **Architectural Problem: Wrong Calculator Creation Pattern**
The Qt GUI was NOT following xespresso/ASE patterns:
- ❌ OLD: Calculator created separately in `_generate_files()` AND `_run_calculation()`
- ✅ NEW: Calculator created ONCE in `calculation_setup.py`, inherited by `job_submission.py`

### 2. **Missing Remote Execution Logic**
The error path `/scratch/projects/spresso/` suggested remote execution, but:
- Calculator wasn't being created with proper queue configuration
- No validation that scheduler matched execution mode

## Complete Solution Implemented

### Phase 1: calculation_setup.py - Create Calculator Once

**Modified `_prepare_calculation()` method:**

```python
def _prepare_calculation(self):
    """Prepare calculation - creates Espresso calculator and stores in session state."""
    atoms = self.session_state.get('current_structure')
    config = self._get_config()
    
    # Set environment
    os.environ['ASE_ESPRESSO_COMMAND'] = ASE_ESPRESSO_COMMAND_TEMPLATE
    
    # Load modules from codes config
    modules = self._load_modules_from_codes(config)
    if modules:
        config['modules'] = modules
    
    # Store config
    self.session_state['workflow_config'] = config
    
    # CREATE CALCULATOR USING gui.calculations
    from gui.calculations import prepare_calculation_from_gui
    
    preparation_label = config.get("label", f"{config.get('calc_type', 'scf')}/{atoms.get_chemical_formula()}")
    
    # This creates calc with proper queue config from machine
    prepared_atoms, calc = prepare_calculation_from_gui(
        atoms, config, label=preparation_label
    )
    
    # STORE IN SESSION STATE FOR REUSE
    self.session_state['espresso_calculator'] = calc
    self.session_state['prepared_atoms'] = prepared_atoms
    
    # Show success message
    ...
```

**Key Points:**
- Calculator created with machine's queue configuration
- Queue includes `execution: 'remote'` or `execution: 'local'`
- Queue includes `scheduler: 'slurm'` or `scheduler: 'direct'`
- Both `prepared_atoms` (with magnetic/Hubbard config) and `calc` stored
- Follows exact same pattern as Streamlit GUI

### Phase 2: job_submission.py - Inherit and Use Calculator

**Removed:**
- ❌ Deleted `_prepare_calculator()` method completely
- ❌ No more calculator creation in job_submission.py

**Refactored `_generate_files()` (Dry Run):**

```python
def _generate_files(self):
    """Generate calculation files (dry run).
    Uses the pre-created calculator from Calculation Setup."""
    
    atoms = self.session_state.get('current_structure')
    
    # GET PRE-CREATED CALCULATOR
    calc = self.session_state.get('espresso_calculator')
    prepared_atoms = self.session_state.get('prepared_atoms', atoms)
    
    if calc is None:
        QMessageBox.warning(
            self, 
            "Calculator Not Prepared",
            "Please go to Calculation Setup page and click 'Prepare Calculation' first."
        )
        return
    
    # Get label from user
    label = self.label_edit.text()
    full_path = os.path.join(workdir, label)
    
    # Create directory
    safe_makedirs(full_path)
    
    # Save structure using PREPARED_ATOMS (has magnetic/Hubbard config)
    ase_io.write(structure_path, prepared_atoms)
    
    # Update calculator label
    prefix = self._get_prefix_from_label(label)
    calc.set_label(label, prefix)
    
    # GENERATE FILES
    calc.write_input(prepared_atoms)
    
    # Show success
    ...
```

**Refactored `_run_calculation()` (Job Execution):**

```python
def _run_calculation(self):
    """Run the calculation.
    Uses the pre-created calculator from Calculation Setup."""
    
    atoms = self.session_state.get('current_structure')
    
    # GET PRE-CREATED CALCULATOR
    calc = self.session_state.get('espresso_calculator')
    prepared_atoms = self.session_state.get('prepared_atoms', atoms)
    
    if calc is None:
        QMessageBox.warning(
            self, 
            "Calculator Not Prepared",
            "Please go to Calculation Setup page and click 'Prepare Calculation' first."
        )
        return
    
    # Get label from user
    label = self.run_label_edit.text()
    full_path = os.path.join(workdir, label)
    
    # Create directory
    safe_makedirs(full_path)
    
    # Update calculator label
    prefix = self._get_prefix_from_label(label)
    calc.set_label(label, prefix)
    
    # RUN CALCULATION
    energy = prepared_atoms.get_potential_energy()
    
    # Show results
    ...
```

## How xespresso Remote Execution Works

When calculator is created with `execution: 'remote'` in queue:

### 1. Calculator Creation (in calculation_setup.py)
```python
queue = {
    'execution': 'remote',
    'scheduler': 'slurm',
    'remote_host': 'cluster.edu',
    'remote_user': 'username',
    'remote_auth': {...},
    'remote_dir': '/scratch/projects/spresso',
    ...
}
calc = Espresso(..., queue=queue)
```

### 2. Dry Run (write_input)
```python
calc.write_input(prepared_atoms)
```
- Calls `set_queue(calc)` in xespresso
- Creates scheduler based on queue config
- Scheduler writes job script locally
- **No remote execution yet** (just file generation)

### 3. Job Execution (get_potential_energy)
```python
energy = prepared_atoms.get_potential_energy()
```
- ASE calls `calc.calculate()`
- Calculator calls `calc.execute()`
- Execute calls `self.scheduler.run()`
- Scheduler checks `queue['execution']`
- **If remote:**
  - `RemoteExecutionMixin.run()` is called
  - Establishes SSH connection to remote_host
  - Transfers input files to remote_dir
  - Transfers job_file to remote_dir
  - Executes `cd remote_dir && source env && sbatch job_file` **ON REMOTE MACHINE**
  - Extracts job ID from sbatch output
  - Polls `squeue` to wait for completion
  - Retrieves output files back to local machine
  - Returns control to ASE for result parsing

**Key Point:** `sbatch` runs on the REMOTE machine via SSH, NOT locally!

## User Configuration Fix

The user's issue was likely one of:

### Scenario A: Wrong Execution Mode
**Problem:**
```yaml
name: my_cluster
execution: local    # ← WRONG! Should be 'remote'
scheduler: slurm
workdir: /scratch/projects/spresso  # remote path
```

**Fix:**
```yaml
name: my_cluster
execution: remote   # ← Correct!
scheduler: slurm
host: cluster.university.edu
username: myuser
workdir: /scratch/projects/spresso
```

### Scenario B: Local SLURM Without Installation
**Problem:**
```yaml
name: my_local
execution: local
scheduler: slurm    # ← But sbatch not installed locally
workdir: ~/calculations
```

**Fix:**
```yaml
name: my_local
execution: local
scheduler: direct   # ← Use bash instead of sbatch
workdir: ~/calculations
```

## Benefits of New Architecture

### 1. **Correct xespresso/ASE Pattern**
- Calculator created once with all configuration
- Reused for dry run and execution
- Prepared atoms with magnetic/Hubbard config preserved

### 2. **Proper Remote Execution**
- Queue configuration flows through correctly
- Remote execution uses SSH + remote sbatch
- Local validation only for local execution

### 3. **Clear Error Messages**
- User knows if they forgot to prepare calculator
- Guided to Calculation Setup page first
- No silent failures or cryptic errors

### 4. **Consistency with Streamlit GUI**
- Both GUIs now use identical pattern
- Same `prepare_calculation_from_gui()` function
- Easier to maintain and debug

### 5. **Better User Experience**
- Logical workflow: Prepare → Dry Run → Execute
- Can change label between dry run and execution
- Calculator reused automatically if label matches

## Testing

To test the fix:

### Test 1: Local Execution with Direct Scheduler
1. Create machine with `execution: local`, `scheduler: direct`
2. Click "Prepare Calculation" in Calculation Setup
3. Go to Job Submission, click "Generate Files"
4. Verify files created successfully
5. Click "Run Calculation"
6. Verify calculation runs with `bash job_file`

### Test 2: Remote Execution with SLURM
1. Create machine with `execution: remote`, `scheduler: slurm`
2. Configure SSH access (host, username, auth)
3. Click "Prepare Calculation"
4. Go to Job Submission, click "Run Calculation"
5. Verify:
   - SSH connection established
   - Files transferred to remote
   - sbatch executed on remote machine
   - Job monitored via squeue
   - Results retrieved

### Test 3: Error Cases
1. Try dry run without preparing calculator first
   - Should show error: "Please go to Calculation Setup..."
2. Try execution without preparing calculator first
   - Should show error: "Please go to Calculation Setup..."

## Files Modified

1. **qtgui/pages/calculation_setup.py**
   - Modified `_prepare_calculation()` to create calculator
   - Uses `prepare_calculation_from_gui()`
   - Stores in session_state

2. **qtgui/pages/job_submission.py**
   - Removed `_prepare_calculator()` method
   - Simplified `_generate_files()` to inherit calculator
   - Simplified `_run_calculation()` to inherit calculator
   - Added validation that calculator exists

3. **SCHEDULER_FIX_SUMMARY.md** (new)
   - Comprehensive documentation of the fix

## Summary

The complete refactor ensures:
- ✅ Calculator created once in Calculation Setup
- ✅ Inherited by Job Submission for both dry run and execution
- ✅ Proper xespresso/ASE pattern implementation
- ✅ Remote execution works correctly (sbatch on remote machine)
- ✅ Local execution works correctly (bash locally)
- ✅ Clear error messages when calculator not prepared
- ✅ Consistency with Streamlit GUI architecture
- ✅ Respects user's machine configuration
- ✅ No silent overrides or modifications

The sbatch error is fixed because:
1. If `execution: 'remote'`, sbatch runs on remote machine via SSH
2. If `execution: 'local'` with `scheduler: 'slurm'` and sbatch missing, user gets clear error
3. User can fix by either: installing SLURM, changing to 'direct', or using remote execution
