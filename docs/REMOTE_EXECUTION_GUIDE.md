# Remote Execution Guide for Wannier Workflows

## Overview

This guide explains how to configure xespresso for remote execution, particularly for multi-stage Wannier workflows where NSCF needs to read density from BANDS calculation results.

## The Challenge

When running computationally intensive Quantum ESPRESSO calculations on a remote cluster via job scheduler (SLURM, PBS, etc.), each stage produces intermediate files that subsequent stages need to read:

```
Local Machine          Remote Cluster
    |                       |
    +---> ssh/SFTP -------> SCF job: creates scf/fe_bcc.save/
    |                       |
    +---> ssh/SFTP -------> BANDS job: reads scf/fe_bcc.save/, creates bands/fe_bcc.save/
    |                       |
    +---> ssh/SFTP -------> NSCF job: reads bands/fe_bcc.save/ ???
```

**Problem**: The NSCF job doesn't know where `bands/fe_bcc.save/` is located on the remote cluster.

**Solution**: Tell xespresso the working directory on the remote cluster using `queue['working_directory']`.

## Configuration

### 1. Basic Remote Queue Setup

```python
from xespresso.workflow.wannier_workflow import WannierWorkflow

# Configuration for SLURM cluster
queue = {
    'execution': 'remote',           # Key: tells xespresso this is remote execution
    'scheduler': 'slurm',
    'working_directory': '/scratch/user/wannier_project',  # Absolute path on remote cluster
    'host': 'cluster.university.edu',
    'user': 'your_username',
    'password': 'your_password',     # Or use SSH keys (preferred)
    'ntasks': 16,
    'mem_per_cpu': '4GB',
    'time': '02:00:00',
}

wf = WannierWorkflow(
    cif_file='Fe.cif',
    pseudopotentials={'Fe': 'Fe.pbe.UPF'},
    protocol='moderate',
    queue=queue,  # Pass queue config
)

# Run workflow (stages execute on remote cluster)
results = wf.run(blocking=True)
```

### 2. Advanced Configuration with SSH Keys

```python
queue = {
    'execution': 'remote',
    'scheduler': 'slurm',
    'working_directory': '/scratch/user/wannier_project',
    'host': 'login.cluster.edu',
    'user': 'your_username',
    'ssh_key': '/home/user/.ssh/id_rsa',  # SSH key authentication
    'ntasks': 32,
    'time': '04:00:00',
}
```

## How It Works

### Execution Flow with Remote Working Directory

```
1. User submits workflow with queue['working_directory'] = '/scratch/user/project'

2. SCF stage:
   - Local: 'scf' directory created locally
   - Remote: Transferred to /scratch/user/project/scf/ on cluster
   - Job runs: pwd = /scratch/user/project/

3. BANDS stage:
   - Remote: Job runs in /scratch/user/project/
   - Reads density from ./scf/fe_bcc.save/ (relative to working_directory)
   - Creates ./bands/fe_bcc.save/

4. NSCF stage:
   - Tells QE: parent_remote_dir = /scratch/user/project/bands
   - QE reads density from /scratch/user/project/bands/fe_bcc.save/
   - Symlinks or absolute paths handled automatically
```

### Path Resolution Logic

The `run_nscf()` method in CalculationWorkflow detects execution mode:

```python
is_remote = self.queue and self.queue.get('execution') == 'remote'

if is_remote:
    # Use explicit remote path
    if parent_remote_dir:
        parent_path = parent_remote_dir
    else:
        # Fallback: use parent_label as relative path
        parent_path = parent_label  # e.g., 'bands'
else:
    # Local: verify path exists locally
    parent_path = Path(parent_label)
    if not parent_path.exists():
        raise FileNotFoundError(...)
```

## Debugging Remote Execution Issues

### Issue: "parent directory not found" Error

**Symptom**: NSCF job fails with message like:
```
Error in routine open_xml_file (1):
   xml data file does not exist: bands/fe_bcc.save/data-file.xml
```

**Solutions**:

1. **Verify working_directory exists on remote**:
   ```bash
   ssh user@cluster.edu "ls -la /scratch/user/project/bands/"
   ```

2. **Check if BANDS output was transferred**:
   ```bash
   ssh user@cluster.edu "find /scratch/user/project/ -name 'data-file.xml' -o -name 'charge-density.dat'"
   ```

3. **Update queue configuration**:
   ```python
   # Wrong:
   queue['working_directory'] = 'project/'  # Relative path won't work
   
   # Correct:
   queue['working_directory'] = '/scratch/user/project'  # Absolute path required
   ```

### Issue: Permission Denied on Remote

**Solution**: Ensure directory is writable:
```bash
ssh user@cluster.edu "mkdir -p /scratch/user/project && chmod 755 /scratch/user/project"
```

### Issue: Different Paths on Local vs Remote

If local working directory differs from remote, you can override:

```python
# Local development directory
local_dir = '/home/user/local_project'

# Remote scratch directory
queue = {
    'execution': 'remote',
    'working_directory': '/scratch/user/remote_project',  # Different from local
    ...
}

# WannierWorkflow handles this automatically via run_nscf(parent_remote_dir=...)
```

## Complete Example: Multi-Stage Remote Workflow

```python
from pathlib import Path
from xespresso.workflow.wannier_workflow import WannierWorkflow

# Configure remote execution
queue = {
    'execution': 'remote',
    'scheduler': 'slurm',
    'working_directory': '/scratch/user/iron_wannier',
    'host': 'hpc.university.edu',
    'user': 'jsmith',
    'ntasks': 24,
    'cpus_per_task': 2,
    'mem_per_cpu': '8GB',
    'time': '06:00:00',
    'partition': 'gpu',  # If using GPU nodes
}

# Create workflow
wf = WannierWorkflow(
    cif_file='Fe_bcc.cif',
    pseudopotentials={
        'Fe': 'Fe.pbe.UPF',
    },
    pseudopotentials_config=None,
    protocol='accurate',
    num_wann=8,
    projections='Fe: dz2, dxz, dyz, dx2-y2, dxy',
    kpts_scf=(4, 4, 4),
    kpts_nscf=(8, 8, 8),
    spinors=False,
    queue=queue,
)

# Run remotely (all stages execute on /scratch/user/iron_wannier/)
try:
    results = wf.run(
        blocking=True,              # Wait for completion
        run_bands_validation=True,  # Include BANDS stage
        run_projwfc_analysis=True,  # Include PROJWFC
    )
    
    print("Wannier workflow completed successfully!")
    print(f"Results directory (remote): {queue['working_directory']}")
    
    # Results are available in results dict
    wannier_output = results['wannier90']
    
except Exception as e:
    print(f"Workflow failed: {e}")
    # Debug: Check remote files
    import subprocess
    result = subprocess.run(
        ['ssh', queue['user']+'@'+queue['host'], 
         'ls', '-la', queue['working_directory']],
        capture_output=True,
        text=True
    )
    print("Remote directory contents:")
    print(result.stdout)
```

## Queue Configuration Reference

| Key | Required | Type | Description | Example |
|-----|----------|------|-------------|---------|
| `execution` | Yes | str | Set to 'remote' for cluster execution | `'remote'` |
| `scheduler` | Yes | str | Job scheduler type | `'slurm'`, `'pbs'`, `'sge'` |
| `working_directory` | **Yes** | str | **Absolute** path on remote where jobs run | `/scratch/user/project` |
| `host` | Yes | str | Cluster login hostname | `cluster.edu` |
| `user` | Yes | str | Username on cluster | `jsmith` |
| `password` | No | str | Password (use ssh_key instead!) | `'secret'` |
| `ssh_key` | No | str | Path to SSH private key | `~/.ssh/id_rsa` |
| `ntasks` | No | int | Number of MPI tasks | `16` |
| `cpus_per_task` | No | int | CPU cores per task | `2` |
| `mem_per_cpu` | No | str | Memory per core | `'4GB'` |
| `time` | No | str | Wall time limit | `'02:00:00'` |
| `partition` | No | str | Cluster partition/queue | `'gpu'`, `'standard'` |

## Security Best Practices

1. **Use SSH keys instead of passwords**:
   ```python
   queue['ssh_key'] = '/home/user/.ssh/id_rsa'
   # Don't use: queue['password'] = '...'
   ```

2. **Protect your SSH key**:
   ```bash
   chmod 600 ~/.ssh/id_rsa
   ```

3. **Use .ssh/config for cluster settings**:
   ```bash
   # ~/.ssh/config
   Host mycluster
       HostName cluster.university.edu
       User jsmith
       IdentityFile ~/.ssh/id_rsa
   ```
   Then reference: `queue['host'] = 'mycluster'`

## Troubleshooting Checklist

- [ ] `working_directory` is an **absolute path** on remote
- [ ] Directory exists on remote cluster: `ssh user@host "ls -d /path/to/dir"`
- [ ] Directory is writable: `ssh user@host "touch /path/to/dir/test && rm /path/to/dir/test"`
- [ ] SSH connectivity works: `ssh user@host "pwd"`
- [ ] No trailing slashes in path: use `/scratch/user/project` not `/scratch/user/project/`
- [ ] Previous stage outputs exist: `ssh user@host "ls -la /path/to/project/bands/"`

## References

- [xespresso CalculationWorkflow](../xespresso/workflow/calculation_workflow.py)
- [xespresso WannierWorkflow](../xespresso/workflow/wannier_workflow.py)
- [Quantum ESPRESSO outdir/prefix](https://www.quantum-espresso.org/input-guide/input-file-description/)
