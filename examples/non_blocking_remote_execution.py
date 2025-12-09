"""
Example: Non-blocking Remote Job Submission

This example demonstrates how to submit remote jobs (both SLURM and direct execution)
without blocking the GUI or main thread. This is useful for interactive applications
where you want the user to maintain control while jobs run on remote clusters.

Key feature: wait_for_completion=False in queue configuration

Supported Schedulers:
- SLURM: Submits to job queue, returns job ID
- Direct: Runs in background, returns PID
"""

from ase.build import bulk
from xespresso import Espresso

# Example 1: Blocking mode with SLURM (default behavior)
# =======================================================
print("Example 1: Blocking Remote Execution with SLURM (default)")
print("=" * 60)

atoms = bulk("Si", cubic=True)
pseudopotentials = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}

# Default behavior - will wait for job to complete
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
    # wait_for_completion is True by default
    # The calculation will block until the job finishes
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

print("Blocking mode will wait for job completion.")
print("GUI will freeze during this time.")


# Example 2: Non-blocking mode (for GUI applications)
# ====================================================
print("\n" + "=" * 60)
print("Example 2: Non-blocking Remote Execution")
print("=" * 60)

# Set wait_for_completion=False to submit job and return immediately
non_blocking_queue = {
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
    # Enable non-blocking mode
    "wait_for_completion": False,
    # Timeout is not used in non-blocking mode
    "job_timeout": 300
}

"""
calc = Espresso(
    label='si_non_blocking',
    pseudopotentials=pseudopotentials,
    queue=non_blocking_queue,
    input_data={'ecutwfc': 30}
)
atoms.calc = calc

# This will submit the job and return immediately
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

print("Non-blocking mode submits job and returns immediately.")
print("GUI remains responsive during job execution.")
print("Job ID is stored in calc.last_job_id for later reference.")


# Example 3: Recommended usage in GUI applications
# =================================================
print("\n" + "=" * 60)
print("Example 3: Recommended Pattern for GUI Applications")
print("=" * 60)

print("""
In a PyQt/PySide GUI application:

1. Submit job in non-blocking mode:
   ```python
   queue = {
       ...
       "wait_for_completion": False
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

2. Provide UI to check status:
   ```python
   # Button or menu item: "Check Job Status"
   def check_status():
       # SSH to remote and check: squeue -j {job_id}
       # Or use RemoteAuth to run command
       pass
   ```

3. Provide UI to retrieve results when done:
   ```python
   # Button or menu item: "Retrieve Results"
   def retrieve_results():
       # Use RemoteAuth to:
       # - Check if output file exists
       # - Download output file
       # - Parse results
       pass
   ```

Benefits:
- GUI never freezes
- User maintains control
- Can submit multiple jobs
- Can cancel operations
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

# Non-blocking direct execution - runs job in background
direct_non_blocking_queue = {
    "execution": "remote",
    "scheduler": "direct",  # Direct execution (no job queue)
    "remote_host": "cluster.university.edu",
    "remote_user": "username",
    "remote_dir": "/home/username/calculations",
    "remote_auth": {
        "method": "key",
        "ssh_key": "~/.ssh/id_rsa"
    },
    # Enable non-blocking mode
    "wait_for_completion": False,
}

"""
calc = Espresso(
    label='al_direct_non_blocking',
    pseudopotentials=pseudopotentials,
    queue=direct_non_blocking_queue,
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

print("Direct scheduler with non-blocking mode:")
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
- SLURM (blocking):    Submits to queue, polls until complete
- SLURM (non-blocking): Submits to queue, returns job ID immediately
- Direct (blocking):    Runs in foreground, blocks until complete  
- Direct (non-blocking): Runs in background, returns PID immediately

General:
========
- Blocking mode (default): wait_for_completion=True
  Good for: Scripts, batch processing
  
- Non-blocking mode: wait_for_completion=False
  Good for: GUI applications, interactive use
  
- Job/Process ID is stored in calc.last_job_id
  - SLURM: Job ID as string (e.g., "12345")
  - Direct: PID prefixed with "PID:" (e.g., "PID:67890")
  
- Output files must be retrieved manually in non-blocking mode
""")
