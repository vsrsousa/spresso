"""
Example: Non-blocking Remote Job Submission

This example demonstrates how to submit remote jobs (both SLURM and direct execution)
without blocking the GUI or main thread. This is useful for interactive applications
where you want the user to maintain control while jobs run on remote clusters.

Key feature: wait_for_completion parameter in queue configuration
- Default (False): Non-blocking mode - job submits and returns immediately
- Set to True: Blocking mode - waits for job completion

Supported Schedulers:
- SLURM: Submits to job queue, returns job ID
- Direct: Runs in background, returns PID
"""

from ase.build import bulk
from xespresso import Espresso

# Example 1: Non-blocking mode with SLURM (default behavior)
# =======================================================
print("Example 1: Non-blocking Remote Execution with SLURM (default)")
print("=" * 60)

atoms = bulk("Si", cubic=True)
pseudopotentials = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}

# Default behavior - job submits and returns immediately (non-blocking)
nonblocking_queue = {
    "execution": "remote",
    "scheduler": "slurm",
    "remote_host": "cluster.university.edu",
    "remote_user": "username",
    "remote_dir": "/home/username/calculations",
    "remote_auth": {
        "method": "key",
        "ssh_key": "~/.ssh/id_rsa"
    },
    "resources": {
        "nodes": 1,
        "ntasks-per-node": 16,
        "time": "02:00:00",
        "partition": "compute"
    },
    # wait_for_completion is False by default (non-blocking mode)
    # The calculation will submit and return immediately
}

"""
calc = Espresso(
    label='si_nonblocking',
    pseudopotentials=pseudopotentials,
    queue=nonblocking_queue,
    input_data={'ecutwfc': 30}
)
atoms.calc = calc

# This will submit the job and return immediately (default behavior)
# The job ID is stored in calc.last_job_id
try:
    # Note: get_potential_energy() will fail because results aren't available yet
    # This is expected behavior in non-blocking mode
    atoms.get_potential_energy()
except Exception as e:
    print(f"Expected: {e}")
    
# Access the job ID
job_id = calc.last_job_id
print(f"Job submitted with ID: {job_id}")
print("GUI remains responsive!")

# Later, you can check job status manually:
# ssh to remote host and run: squeue -j {job_id}
# Or retrieve results when job is done
"""

print("Non-blocking mode (default) submits job and returns immediately.")
print("GUI remains responsive during job execution.")
print("Job ID is stored in calc.last_job_id for later reference.")


# Example 2: Blocking mode (for scripts that need results immediately)
# =====================================================================
print("\n" + "=" * 60)
print("Example 2: Blocking Remote Execution")
print("=" * 60)

# Set wait_for_completion=True to wait for job completion
blocking_queue = {
    "execution": "remote",
    "scheduler": "slurm",
    "remote_host": "cluster.university.edu",
    "remote_user": "username",
    "remote_dir": "/home/username/calculations",
    "remote_auth": {
        "method": "key",
        "ssh_key": "~/.ssh/id_rsa"
    },
    "resources": {
        "nodes": 1,
        "ntasks-per-node": 16,
        "time": "02:00:00",
        "partition": "compute"
    },
    # Enable blocking mode (wait for completion)
    "wait_for_completion": True,
    # Timeout used in blocking mode
    "job_timeout": 3600
}

"""
calc = Espresso(
    label='si_blocking',
    pseudopotentials=pseudopotentials,
    queue=blocking_queue,
    input_data={'ecutwfc': 30}
)
atoms.calc = calc

# This will block until the remote job completes
energy = atoms.get_potential_energy()
print(f"Energy: {energy} eV")
"""

print("Blocking mode waits for job completion.")
print("GUI will freeze during this time (not recommended for interactive use).")


# Example 3: Recommended usage in GUI applications
# =================================================
print("\n" + "=" * 60)
print("Example 3: Recommended Pattern for GUI Applications")
print("=" * 60)

print("""
In a PyQt/PySide GUI application:

1. Submit job in non-blocking mode (default):
   ```python
   queue = {
       ...
       # wait_for_completion defaults to False
   }
   calc = Espresso(label='my_calc', queue=queue, ...)
   atoms.calc = calc
   
   # Submit job (returns immediately)
   try:
       atoms.get_potential_energy()
   except:
       pass  # Expected in non-blocking mode
   
   job_id = calc.last_job_id
   # Display: "Job {job_id} submitted successfully!"
   ```

2. Use the Job Monitor dialog to track jobs:
   - Jobs automatically appear in Job Monitor
   - Check status with "Refresh" button
   - Retrieve results when complete
   - Jobs persist across GUI sessions

Benefits:
- GUI never freezes
- User maintains control
- Can submit multiple jobs
- Jobs tracked in Job Monitor
- Better user experience
""")

# ==================================================================
# Example 4: Non-blocking mode with Direct Scheduler
# ==================================================================
print("\n" + "=" * 60)
print("Example 4: Non-blocking Remote Execution with Direct Scheduler")
print("=" * 60)

atoms = bulk("Al", cubic=True)
pseudopotentials = {"Al": "Al.pbe.UPF"}

# Non-blocking direct execution - runs job in background (default)
direct_nonblocking_queue = {
    "execution": "remote",
    "scheduler": "direct",  # Direct execution (no job queue)
    "remote_host": "cluster.university.edu",
    "remote_user": "username",
    "remote_dir": "/home/username/calculations",
    "remote_auth": {
        "method": "key",
        "ssh_key": "~/.ssh/id_rsa"
    },
    # Non-blocking mode is the default
}

"""
calc = Espresso(
    label='al_direct_nonblocking',
    pseudopotentials=pseudopotentials,
    queue=direct_nonblocking_queue,
    input_data={'ecutwfc': 30}
)
atoms.calc = calc

# This will start the job in background on remote host and return immediately
try:
    atoms.get_potential_energy()
except Exception as e:
    print(f"Expected: {e}")
    
# Access the process ID (stored as "PID:12345")
pid = calc.last_job_id  # e.g., "PID:12345"
print(f"Job started with {pid}")
print("GUI remains responsive!")

# Later, check if process is still running:
# ssh to remote and run: ps -p {pid_number}
# Or kill the process: kill {pid_number}
"""

print("Direct scheduler with non-blocking mode (default):")
print("- Runs bash script in background on remote host")
print("- Returns PID immediately (stored as 'PID:12345' in calc.last_job_id)")
print("- Output redirected to .pwo file")
print("- Check status: ps -p {pid}")
print("- Terminate: kill {pid}")

print("\n" + "=" * 60)
print("Summary")
print("=" * 60)
print("""
Scheduler Support:
==================
- SLURM (non-blocking): Submits to queue, returns job ID immediately [DEFAULT]
- SLURM (blocking):     Submits to queue, polls until complete
- Direct (non-blocking): Runs in background, returns PID immediately [DEFAULT]
- Direct (blocking):     Runs in foreground, blocks until complete  

General:
========
- Non-blocking mode (default): wait_for_completion=False
  Good for: GUI applications, interactive use
  
- Blocking mode: wait_for_completion=True
  Good for: Scripts, batch processing where results are needed immediately
  
- Job/Process ID is stored in calc.last_job_id
  - SLURM: Job ID as string (e.g., "12345")
  - Direct: PID prefixed with "PID:" (e.g., "PID:67890")
  
- Use Job Monitor in GUI to track jobs and retrieve results
- Output files must be retrieved manually in non-blocking mode (or use Job Monitor)
""")
