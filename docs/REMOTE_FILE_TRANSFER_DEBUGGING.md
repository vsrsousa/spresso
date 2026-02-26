# Remote File Transfer Troubleshooting Guide

## Problem
When using batch mode with remote execution (`machine="medusa"`), input files are generated locally but not transferred to the remote machine, causing job submissions to fail silently.

## Root Cause
The `execute()` method in [xespresso/xespresso.py](../xespresso/xespresso.py) was catching all exceptions during scheduler initialization without logging them. This meant:

1. **Input files ARE generated locally** ✓
2. **execute() is called** ✓  
3. **Scheduler initialization fails silently** ✗ (error is caught, no logging)
4. **Synthetic results returned instead** ✗ (user thinks it worked but files weren't transferred)

## Solution Applied

### 1. **Better Error Logging**
Modified `execute()` to log all scheduler initialization and execution errors:
- Scheduler initialization errors are logged with full stack trace
- Scheduler execution errors are logged with full stack trace
- Errors can be made fatal by setting `XESPRESSO_VERBOSE_ERRORS=1`

### 2. **How to Enable Verbose Error Reporting**
```bash
# Enable verbose error mode to see the actual error instead of silent failure
export XESPRESSO_VERBOSE_ERRORS=1

# Then run your code
python -c "from examples.converge_example import prepare_example; ..."
```

With verbose errors enabled, you'll see the actual error that's preventing file transfer (e.g., SSH connection failed, authentication failed, etc.)

## How to Debug

### Quick Test: Run the Debug Script
```bash
cd /home/vinicius/scratch/projects/spresso

# Enable verbose errors
export XESPRESSO_VERBOSE_ERRORS=1

# Run debug tests (step by step)
python examples/debug_remote_transfer.py
```

This script tests:
1. Machine configuration loading
2. Queue setup
3. Input file generation
4. Remote connection (SSH)
5. Scheduler setup

Each test will show you exactly where the problem is.

### Step-by-Step Debugging

#### Step 1: Verify Machine Configuration
```python
from xespresso.machines import load_machine
machine = load_machine('medusa')  # Should load successfully
print(machine)
```

#### Step 2: Verify SSH Connectivity
```bash
# Check if you can SSH to the remote host
ssh medusa.fis.uerj.br "pwd"

# Check if your SSH key is set up correctly
ssh-keyscan -t rsa medusa.fis.uerj.br >> ~/.ssh/known_hosts
```

#### Step 3: Enable Logging in Code
```python
import logging
logging.basicConfig(level=logging.DEBUG)

os.environ['XESPRESSO_VERBOSE_ERRORS'] = '1'

result = ConvergenceWorkflow.optimize_parameters(
    atoms=atoms,
    pseudopotentials_config="default",
    precision='low',
    machine="medusa",
    code_version="7.4.1",
    use_batch_mode=True,
    verbose=True
)
```

## Common Issues and Fixes

### Issue 1: SSH Connection Failed
**Symptom**: "Could not connect to remote host"

**Fix**:
```bash
# Test SSH connection
ssh -v medusa.fis.uerj.br "echo 'Connected!'"

# Check SSH key configuration
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_rsa  # or your key file

# Verify host key is in known_hosts
ssh-keyscan -t rsa medusa.fis.uerj.br >> ~/.ssh/known_hosts
```

### Issue 2: Authentication Failed  
**Symptom**: "Permission denied (publickey,password)"

**Fix**:
```bash
# Check if SSH key permissions are correct
chmod 600 ~/.ssh/id_rsa
chmod 700 ~/.ssh

# Test with password authentication (if configured)
ssh -o PubkeyAuthentication=no medusa.fis.uerj.br
```

### Issue 3: Remote Directory Access Issues
**Symptom**: "Permission denied" when creating remote directories

**Fix**:
```bash
# Check if remote_dir exists and is writable
ssh medusa.fis.uerj.br "ls -la ~/xespresso_jobs/"

# Or check what your default remote_dir is in machines.json
cat ~/.xespresso/machines/medusa.json | grep remote_dir
```

## Files Modified

1. **xespresso/xespresso.py**
   - Added better error logging in `execute()` method
   - Added `VERBOSE_ERRORS` import from config
   - Errors are now logged with full stack traces instead of being silently caught
   - Can raise errors if `XESPRESSO_VERBOSE_ERRORS=1`

2. **examples/debug_remote_transfer.py** (NEW)
   - Comprehensive debug script that tests each part of the remote transfer pipeline
   - Run to identify exactly where the problem is

## Testing the Fix

### Test 1: With Verbose Errors Disabled (Silent Fallback)
```python
os.environ['XESPRESSO_VERBOSE_ERRORS'] = '0'

# Should see error logged to stderr/stdout
result = ConvergenceWorkflow.optimize_parameters(
    atoms=atoms,
    machine="medusa",
    use_batch_mode=True,
)
# ERROR message would show up in logs
```

### Test 2: With Verbose Errors Enabled (Exception Raised)
```python
os.environ['XESPRESSO_VERBOSE_ERRORS'] = '1'

# Should raise exception with full stack trace
result = ConvergenceWorkflow.optimize_parameters(
    atoms=atoms,
    machine="medusa", 
    use_batch_mode=True,
)
# Exception raised with full traceback
```

## Expected Behavior After Fix

### Success Case (SSH Connected, Files Transferred)
```
INFO [execute]: Executing job via scheduler...
INFO [set_queue]: Scheduler executed with command: sbatch ...
INFO [submit_scf_batch]: Submitting SCF calculation: convergence/ecut30_ksp0.50
⚡ BATCH MODE: Submitting jobs in parallel to SLURM
[1/5] convergence/ecut30_ksp0.50 (job_id: 12345)
[2/5] convergence/ecut30_ksp0.40 (job_id: 12346)
...
```

### Failure Case (Now Shows Error Instead of Silent Failure)
```
ERROR [execute]: Failed to initialize scheduler: ConnectionError: SSH connection failed to medusa.fis.uerj.br
Traceback...

⚠️  Remote connection failed - this would prevent file transfer
Check:
  1. SSH key is properly configured
  2. Remote host is reachable  
  3. Remote user credentials are correct
```

## Next Steps

1. Run the debug script: `python examples/debug_remote_transfer.py`
2. Identify where it fails
3. Fix the specific issue (SSH, credentials, etc.)
4. Re-run the convergence workflow

Once SSH connectivity is verified, the batch file transfer will work automatically in `scheduler.run()` through `RemoteExecutionMixin`.
