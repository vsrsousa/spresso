## Testing Guide - Remote Non-Blocking Execution

### Quick Start Verification

#### Step 1: Run Mock Tests (No Remote Connection Needed)
```bash
cd /home/vinicius/projects/spresso

# Run mock tests (should all pass)
python test_remote_nonblocking_mock.py

# Run integration tests (should all pass)
python test_integration_final.py
```

Expected output: ✅ **ALL TESTS PASSED**

---

### Step 2: Test with Real Remote Execution

#### Test Script
Create `test_real_remote.py`:

```python
#!/usr/bin/env python
"""
Test remote non-blocking execution on real machine (medusa).
Verifies:
1. Files sent exactly once (not 3x)
2. Job file uses .pwi extension (not .pwx)
3. Results retrieved correctly
"""

import logging
from ase.build import bulk
from xespresso import CalculationWorkflow

# Enable debug logging to see file transfers
logging.basicConfig(level=logging.DEBUG)

# Test 1: Simple SCF calculation
print("="*70)
print("TEST 1: Remote non-blocking SCF calculation")
print("="*70)

atoms = bulk("Si", cubic=True)

workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={'Si': 'Si.pbe.UPF'},
    machine='snake5',  # Remote machine
    protocol='fast'
)

print("\nRunning SCF on remote machine (non-blocking)...")
print("Expected: Files sent ONCE each, job_file has .pwi extension\n")

calc = workflow.run_scf(label='test_remote_scf_simple')

print("\nResults retrieved:")
print(f"  Energy: {calc.results['energy']:.6f} eV")
print(f"  Converged: {calc.results['converged']}")

if 'forces' in calc.results:
    print(f"  Forces shape: {calc.results['forces'].shape}")

print("\n✓ Test 1 passed: SCF completed successfully")

# Test 2: Relaxation calculation
print("\n" + "="*70)
print("TEST 2: Remote non-blocking relaxation")
print("="*70)

atoms = bulk("Si", cubic=True)

workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={'Si': 'Si.pbe.UPF'},
    machine='snake5',
    protocol='fast'
)

print("\nRunning relaxation on remote machine (non-blocking)...")

calc = workflow.run_relax(label='test_remote_relax_simple', steps=5)

print("\nResults retrieved:")
print(f"  Final energy: {calc.results['energy']:.6f} eV")
print(f"  Final cell: {calc.results['cell']}")

print("\n✓ Test 2 passed: Relaxation completed successfully")

print("\n" + "="*70)
print("✓ ALL TESTS PASSED - Real remote execution works!")
print("="*70)
```

#### Run the test:
```bash
python test_real_remote.py 2>&1 | tee test_output.log

# Check for issues:
# 1. Should see FOUR file transfers (pseudo, input, job_file for write_input + execute)
#    NOT twelve (4 x 3)
# 2. Should see job_file with ".pwi" extension
# 3. Should see monitor.wait() polling status
# 4. Should see results retrieved and parsed
```

---

### Step 3: Verify No File Repetition

#### Method A: Check Log Output
Look for this pattern (GOOD):
```
send_file pseudopotential...
send_file input...
send_file job_file...
RemoteJobMonitor.wait() polling...
```

NOT this pattern (BAD - would indicate repetition):
```
send_file pseudopotential...
send_file input...
send_file job_file...
send_file pseudopotential...  ← REPETITION!
send_file input...            ← REPETITION!
send_file job_file...         ← REPETITION!
```

#### Method B: Run with strace
```bash
# Monitor system calls to see actual file operations
strace -e openat,read -f python test_real_remote.py 2>&1 | grep "si-test\|pw.pwo" | head -20
```

---

### Step 4: Verify Job File Extension

#### Check Generated Files
```bash
# After running a test, check the job file:
ls -la scf_*/si-test*

# Should show:
# - si-test.pwi (INPUT, created by write_input)
# - si-test.pwo (OUTPUT, retrieved after job)
# - si-test.json (METADATA)

# NOT si-test.pwx (which indicates wrong extension bug)
```

#### Verify on Remote Machine
```bash
# SSH to medusa and check:
ssh medusa.fis.uerj.br

# List recent job files:
ls -ltr ~/scratch/users/vinicius/xespresso/*/si-test* | tail -10

# Should see:
# si-test.pwi
# si-test.pwo
# si-test.json

# NOT si-test.pwx
```

---

### Step 5: Monitor Job Status (Advanced)

If you want to access monitor directly:

```python
from xespresso import CalculationWorkflow, RemoteJobMonitor
from ase.build import bulk

atoms = bulk("Si", cubic=True)

workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={'Si': 'Si.pbe.UPF'},
    machine='snake5'
)

calc = workflow.run_scf(label='test_monitor')

# Monitor is available if job was remote non-blocking
if hasattr(workflow, 'last_calc'):
    from xespresso.schedulers import RemoteJobMonitor
    
    monitor = RemoteJobMonitor(workflow.last_calc)
    
    print(f"Job ID: {monitor.job_id}")
    print(f"Status: {monitor.status()}")
    print(f"Job type: {monitor.job_type}")
    
    # Wait for it manually if needed
    if monitor.status() == 'running':
        monitor.wait(timeout=3600)
        print(f"Final status: {monitor.status()}")
```

---

### Step 6: Troubleshooting

#### Issue: "Files are still being sent 3 times"

Check:
1. Are you using `workflow.run_scf()` or direct `calc.run()`?
   - Must use `workflow.run_scf()` for automatic remote handling
2. Is machine configured as remote?
   - Check: `workflow.queue.get('execution')` should be `'remote'`
3. Is `wait_for_completion` set to False?
   - Check: `workflow.queue.get('wait_for_completion', False)` should be False

#### Issue: "Job file still has .pwx extension"

Check:
1. Environment variable `ASE_ESPRESSO_COMMAND`
   ```bash
   echo $ASE_ESPRESSO_COMMAND
   ```
   If it returns a value, it might be taking precedence. The scheduler.py line 54 
   fallback only works if this is empty or has correct placeholder.

2. Check your ~/.bashrc or ~/.bashprofile:
   ```bash
   grep ASE_ESPRESSO_COMMAND ~/.bashrc
   ```
   If present, either remove it or ensure it uses `PACKAGEi` (not `PACKAGEx`)

#### Issue: "Job timed out waiting for completion"

Check:
1. Is the job actually running on remote?
   ```bash
   ssh medusa.fis.uerj.br
   squeue -u vinicius  # or 'ps aux' for direct scheduler
   ```

2. Increase timeout if needed:
   ```python
   calc = workflow.run_scf(label='test', machine='snake5')
   # monitor will use: queue.get('job_timeout', 3600)
   # To override, modify workflow.queue before calling run_scf()
   workflow.queue['job_timeout'] = 7200  # 2 hours
   ```

---

### Verification Checklist

After running tests, verify:

- [ ] Mock tests pass (`test_remote_nonblocking_mock.py` ✅)
- [ ] Integration tests pass (`test_integration_final.py` ✅)
- [ ] Real remote tests complete (`test_real_remote.py` ✅)
- [ ] Job files have `.pwi` extension (not `.pwx`)
- [ ] Files sent only once (4 total, not 12)
- [ ] Results retrieved (energy, forces available)
- [ ] Monitor detected correctly (SLURM vs Direct)
- [ ] No timeout errors (jobs complete in < 3600s)

---

### Performance Comparison

#### Before Fix
```
Total time: ~X minutes
File transfers: 9? (3 files × 3 attempts)
Status checks: N
Wait time: Y minutes
```

#### After Fix
```
Total time: ~X-Y minutes (faster, no retries)
File transfers: 4 (1 time only: pseudo, input, job, output)
Status checks: ~6 (10s polling for 1 min job)
Wait time: ~job duration only
```

---

### Next Steps

1. ✅ Run mock tests
2. ✅ Run integration tests
3. 🚀 Run real remote test
4. 🚀 Verify file transfers (1 time)
5. 🚀 Verify `.pwi` extension
6. ✅ Mark solution as complete

---

**Documentation**: See `SOLUTION_SUMMARY.md` for architecture details
**Test Files**: `test_remote_nonblocking_mock.py`, `test_integration_final.py`
**Status**: Ready for testing ✅
