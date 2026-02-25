## Remote Non-Blocking Execution - Complete Solution

### Summary

This document describes the complete solution for remote non-blocking execution in Spresso, which addresses:

1. **`.pwx` extension bug** - Fixed by adding fallback in `scheduler.py`
2. **File transfer repetition** - Eliminated by intercepting remote non-blocking in workflow
3. **Transparent execution** - RemoteJobMonitor handles monitoring automatically
4. **Non-blocking default** - Users get results transparently without explicit monitor calls

---

## Problems Solved

### Problem 1: Wrong Job File Extension (`.pwx` instead of `.pwi`)

**Root Cause:**
- `ASE_ESPRESSO_COMMAND` environment variable in `~/.bashrc` set with incorrect placeholder
- When env var was empty or unset, code had no fallback to correct template

**Solution:**
- **File**: `xespresso/scheduler.py` (line 54)
- **Change**: Added fallback chain:
  ```python
  command = command or os.environ.get("ASE_ESPRESSO_COMMAND", "") or calc.command
  ```
- **Effect**: When env var empty, uses `calc.command` with correct `.pwi` extension
- **Result**: ✅ Job files now have `.pwi` (not `.pwx`)

---

### Problem 2: Files Transferred 3x Instead of 1x

**Root Cause:**
- `xespresso.run()` has while loop checking convergence
- For non-blocking jobs, `.pwo` file not available yet (job still running)
- `read_convergence()` fails → loop repeats
- Each iteration calls: send_file(pseudo) + send_file(input) + send_file(job_file)

**Solution:**
- **File**: `xespresso/workflow/simple_workflow.py`
- **Methods**: `run_scf()` (lines 515-543) and `run_relax()` (lines 600-628)
- **Change**: Manual execution flow for remote non-blocking:
  ```python
  if self.queue and self.queue.get('execution') == 'remote' \
     and not self.queue.get('wait_for_completion', False):
      # Remote non-blocking path (avoids calc.run() retry loop)
      calc.write_input(self.atoms)      # Write input once
      calc.execute()                    # Submit job once
      monitor = RemoteJobMonitor(calc)  # Monitor waits
      monitor.wait(timeout=3600)        # Wait for completion
      monitor.retrieve_output()         # Get output
      calc.read_results()               # Read results
  else:
      # Local or blocking path
      calc.run(atoms=self.atoms)        # Normal flow
  ```
- **Effect**: Skips `calc.run()` entirely, avoids retry loop
- **Result**: ✅ Files transferred exactly once

---

### Problem 3: Monitoring Requires Manual User Code

**Root Cause:**
- RemoteJobMonitor existed but wasn't integrated into workflow
- Users had to manually call monitor methods

**Solution:**
- **File**: `xespresso/workflow/simple_workflow.py`
- **Integration**: Monitor called automatically inside `run_scf()` and `run_relax()`
- **Effect**: Transparent to user - just call `workflow.run_scf()` and wait
- **Result**: ✅ No manual monitor calls needed

---

## Architecture

### Execution Flow

```
User calls: workflow.run_scf(label='scf/si-test', machine='remote_machine')
    ↓
workflow creates Espresso calculator
    ↓
DETECTS: remote + non-blocking?
    ↓ YES
    ├─ write_input(atoms)              → Writes input files locally
    ├─ execute()                       → Submits to remote machine (ONCE)
    ├─ RemoteJobMonitor.wait()         → Polls job status
    │  ├─ SLURM: squeue + sacct
    │  └─ Direct: ps + check output file
    ├─ retrieve_output()               → Gets .pwo, .json from remote
    └─ read_results()                  → Parses results (energy, forces, etc)
    ↓
RETURNS: calc with results populated
    ↓
User has: calc.results['energy'], calc.results['forces'], etc

    ↓ NO (local or blocking)
    use calc.run(atoms=self.atoms)     → Normal Espresso flow
```

### File Transfer Flow (Remote Non-Blocking)

**Before fix:**
```
write_input() → pseudo + input + job_file sent    (1st attempt)
read_convergence() fails (job running)
retry → pseudo + input + job_file sent            (2nd attempt)
read_convergence() fails (still running)
retry → pseudo + input + job_file sent            (3rd attempt)
```

**After fix:**
```
write_input() → pseudo + input + job_file sent    (1st time only)
execute()     → Job submitted
monitor.wait()→ Polling only (no file transfers)
read_results()→ Results parsed
```

---

## Implementation Details

### 1. RemoteJobMonitor Class

**File**: `xespresso/schedulers/remote_job_monitor.py`

**Features:**
- Auto-detects job type from `job_id` format:
  - SLURM: numeric job ID (e.g., `12345`)
  - Direct: PID format (e.g., `PID:54321`)
- Methods:
  - `status()`: Check job state
  - `wait()`: Poll until completion
  - `retrieve_output()`: Get output files from remote
  - `info()`: Detailed job information
  - `cancel()`: Stop job

**Usage:**
```python
monitor = RemoteJobMonitor(calc)
monitor.wait(timeout=3600, poll_interval=10)
monitor.retrieve_output()
```

### 2. Workflow Integration

**File**: `xespresso/workflow/simple_workflow.py`

**Key additions:**
- Import RemoteJobMonitor at top
- In `run_scf()` and `run_relax()`: detect remote + non-blocking
- Use manual flow for remote non-blocking
- Normal `calc.run()` for local or blocking jobs

**Attributes:**
- `self.last_calc`: Track last calculator for potential monitor access
- `self.queue`: Configuration from machine (contains `execution`, `wait_for_completion`, `job_timeout`)

### 3. Scheduler Fallback

**File**: `xespresso/scheduler.py` (line 54)

**Logic:**
```python
command = command or os.environ.get("ASE_ESPRESSO_COMMAND", "") or calc.command
```

**Precedence:**
1. Explicit `command` argument
2. `ASE_ESPRESSO_COMMAND` environment variable
3. `calc.command` (default template with correct `.pwi`)

---

## Files Modified

| File | Lines | Change | Impact |
|------|-------|--------|--------|
| `xespresso/scheduler.py` | 54 | Add fallback to `calc.command` | .pwi extension fixed |
| `xespresso/workflow/simple_workflow.py` | 20, 515-543, 600-628 | Add RemoteJobMonitor logic | Remote non-blocking handling |
| `xespresso/schedulers/remote_job_monitor.py` | NEW | New RemoteJobMonitor class | Unified job monitoring |
| `xespresso/schedulers/__init__.py` | 3, 5 | Export RemoteJobMonitor | Available for import |

---

## Testing

### Mock Tests
**File**: `test_remote_nonblocking_mock.py`
- Tests remote non-blocking detection
- Tests RemoteJobMonitor with SLURM/direct jobs
- Tests workflow integration
- Status: ✅ All pass

### Integration Tests
**File**: `test_integration_final.py`
- Tests command fallback (`.pwi` extension)
- Tests workflow configuration detection
- Tests manual execution flow
- Tests RemoteJobMonitor functionality
- Status: ✅ All pass

---

## User Experience

### Before Fix
```python
# Problem: Files repeated 3 times, wrong extension
calc = workflow.run_scf(label='scf/si-test', machine='remote_machine')
# Output: send_file pseudo (3x), send_file input (3x), send_file job (3x)
# File: si-test.pwx (WRONG!)
```

### After Fix
```python
# Works transparently!
calc = workflow.run_scf(label='scf/si-test', machine='remote_machine')
# Output: send_file pseudo (1x), send_file input (1x), send_file job (1x)
# File: si-test.pwi (CORRECT!)
# Internally: monitor waits for job, retrieves output
# Returns: calc with energy, forces, etc ready to use
```

---

## Configuration

### Machine Configuration (queue parameter)
```python
queue = {
    'execution': 'remote',           # Local vs remote
    'scheduler': 'slurm',            # SLURM or 'direct'
    'wait_for_completion': False,    # Non-blocking (default)
    'job_timeout': 3600,             # 1 hour timeout
    'launcher': 'srun --mpi=pmi2',   # SLURM launcher
    'nprocs': 16,                    # Number of processes
    'modules': ['quantum-espresso/7.4.1'],
}
```

---

## Verification

Run tests to verify solution:

```bash
# Mock tests
python test_remote_nonblocking_mock.py

# Integration tests
python test_integration_final.py
```

Both should show: ✅ **ALL TESTS PASSED**

---

## Backward Compatibility

✅ **Fully backward compatible:**
- Local calculations unaffected (use `calc.run()` normally)
- Remote blocking calculations unaffected (wait_for_completion=True)
- All existing code works as before
- Non-blocking default doesn't break existing workflows

---

## Known Limitations & Future Improvements

1. **Monitor polling interval**: Currently hardcoded to 10 seconds
   - Could be made configurable via queue settings

2. **Timeout handling**: Simple timeout, no job restart
   - Could implement automatic retry for failed jobs

3. **Output file staging**: Manual `retrieve_output()` call
   - Could be extended to auto-stage additional files

4. **Job reuse**: Monitor designed for single job tracking
   - Could support job arrays or multiple simultaneous jobs

---

## Summary Table

| Issue | Root Cause | Fix | Status |
|-------|-----------|-----|--------|
| Wrong .pwx extension | No fallback in scheduler | Added calc.command fallback | ✅ Fixed |
| 3x file repetition | calc.run() retry loop | Manual execution in workflow | ✅ Fixed |
| Complex user interface | Manual monitor calls | Integrated monitor in workflow | ✅ Fixed |

---

**Last Updated**: After integration test validation
**Status**: Ready for production use ✅
