## Remote Job Monitoring

The `RemoteJobMonitor` class provides a unified interface to monitor remote job execution for both SLURM and direct (bash) schedulers.

### Overview

When running calculations remotely in **non-blocking mode** (`wait_for_completion=False`), the job is submitted to the remote cluster and execution returns immediately. The `RemoteJobMonitor` allows you to:

- Check job status (running, completed, failed)
- Wait for job completion with timeout
- Retrieve output files automatically
- Get detailed job information
- Cancel running jobs

### Installation

The `RemoteJobMonitor` is automatically available when imported:

```python
from xespresso.schedulers import RemoteJobMonitor
from xespresso import CalculationWorkflow
```

### Basic Usage

#### Via CalculationWorkflow

The easiest way to monitor is through the workflow's `get_monitor()` method:

```python
from ase.build import bulk
from xespresso import CalculationWorkflow

atoms = bulk("Si", cubic=True)
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials_config='default',
    machine='snake5'  # Remote machine
)

# Run calculation (non-blocking mode)
calc = workflow.run_scf(label='scf/si-test')

# Get monitor
monitor = workflow.get_monitor()

# Check status
status = monitor.status()  # 'running', 'completed', 'failed'

# Wait for completion
if monitor.wait(timeout=3600):  # Wait up to 1 hour
    local_path, output = monitor.retrieve_output()
    print(f"Output: {local_path}")
```

#### Direct Usage

You can also create a monitor directly from a calculator:

```python
from xespresso.schedulers import RemoteJobMonitor

monitor = RemoteJobMonitor(calc)
```

### API Reference

#### Constructor

```python
RemoteJobMonitor(calc, remote_connection=None)
```

- `calc`: Espresso calculator with `last_job_id` and `last_remote_path` (set by non-blocking remote execution)
- `remote_connection`: Optional remote connection object (auto-detected if available)

#### Methods

##### `status() -> str`
Get current job status.
- Returns: `'running'`, `'completed'`, `'failed'`, or `'unknown'`

```python
status = monitor.status()
print(f"Job status: {status}")
```

##### `wait(timeout: int = 3600, poll_interval: int = 10) -> bool`
Wait for job completion.
- `timeout`: Maximum time to wait in seconds (default: 1 hour)
- `poll_interval`: How often to check status in seconds (default: 10)
- Returns: `True` if completed, `False` if timeout

```python
if monitor.wait(timeout=300):  # 5 minute timeout
    print("Job completed!")
else:
    print("Job timed out")
```

##### `retrieve_output() -> Tuple[str, str]`
Retrieve output files from remote.
- Returns: `(local_output_path, output_file_content)`

```python
local_path, content = monitor.retrieve_output()
with open(local_path) as f:
    results = f.read()
```

##### `info() -> Dict`
Get detailed job information.
- Returns dictionary with: `job_id`, `job_type`, `status`, `remote_path`, `time_used` (SLURM), `time_limit` (SLURM)

```python
info = monitor.info()
print(f"Job ID: {info['job_id']}")
print(f"Status: {info['status']}")
if 'time_used' in info:
    print(f"Time used: {info['time_used']}")
```

##### `cancel() -> bool`
Cancel the running job.
- Returns: `True` if cancellation succeeded, `False` otherwise

```python
if monitor.cancel():
    print("Job cancelled successfully")
```

### Advanced Examples

#### Poll Job Status Manually

```python
import time

monitor = workflow.get_monitor()

# Check status every 30 seconds
while True:
    status = monitor.status()
    print(f"Current status: {status}")
    
    if status == 'completed':
        print("Job finished!")
        break
    elif status == 'failed':
        print("Job failed!")
        break
    
    time.sleep(30)
```

#### Get Partial Results While Running

```python
monitor = workflow.get_monitor()

# Check status periodically
count = 0
while count < 10:  # Check up to 10 times
    status = monitor.status()
    info = monitor.info()
    
    print(f"[{count}] Status: {status}")
    if 'time_used' in info:
        print(f"    Time used: {info['time_used']}")
    
    if status != 'running':
        break
    
    count += 1
    time.sleep(30)
```

#### Cancel Job if Taking Too Long

```python
monitor = workflow.get_monitor()

if not monitor.wait(timeout=300):  # 5 minute timeout
    print("Job taking too long, cancelling...")
    monitor.cancel()
```

### Supported Schedulers

| Scheduler | Status Check | Info | Cancel |
|-----------|--------------|------|--------|
| SLURM | `squeue`/`sacct` | Job time, time limit | `scancel` |
| Direct (bash) | `ps` | Process ID | `kill` |

### Common Patterns

#### Pattern 1: Submit and Forget (with later monitoring)

```python
# Submit job
calc = workflow.run_scf(label='scf/si')

# Do other work...

# Later, check if done
monitor = workflow.get_monitor()
if monitor.status() == 'completed':
    output_path, output = monitor.retrieve_output()
```

#### Pattern 2: Wait with Progress Updates

```python
monitor = workflow.get_monitor()

print(f"Job {monitor.info()['job_id']} started")

if monitor.wait(timeout=3600, poll_interval=60):  # Check every minute
    print("Job completed successfully!")
    monitor.retrieve_output()
else:
    print("Job timed out")
    print(f"Last status: {monitor.status()}")
```

#### Pattern 3: Multiple Jobs with Different Monitors

```python
# Run multiple jobs in parallel
jobs = []
for label in ['scf/si', 'scf/ge', 'scf/sn']:
    calc = workflow.run_scf(label=label)
    monitor = RemoteJobMonitor(calc)
    jobs.append((label, monitor))

# Wait for first completion
for label, monitor in jobs:
    if monitor.status() == 'completed':
        print(f"{label} finished!")
        monitor.retrieve_output()
        break
```

### Troubleshooting

#### Monitor Creation Fails

```python
RemoteJobMonitor(calc)  # ValueError: No remote connection available
```

Make sure:
1. Job was run with remote execution enabled (`machine='...'`)
2. Non-blocking mode was active (`wait_for_completion=False`)
3. Remote connection is still available

#### Status Always "Unknown"

Possible causes:
1. Remote connection lost
2. Job ID format not recognized
3. SLURM/bash not available on remote

#### Output Retrieval Fails

```python
monitor.retrieve_output()  # RuntimeError: Output file does not exist
```

Possible causes:
1. Job failed on remote
2. Output file not generated yet (check status first)
3. File system permission issues

### See Also

- [examples/monitor_remote_execution.py](./examples/monitor_remote_execution.py)
- [examples/workflow_monitor_example.py](./examples/workflow_monitor_example.py)
- [Remote Execution Documentation](./REMOTE_EXECUTION.md)
