"""
Example: Monitor remote job execution without blocking.

This example demonstrates how to:
1. Run a calculation remotely in non-blocking mode
2. Monitor job progress
3. Retrieve results when done
"""

from ase.build import bulk
from xespresso import CalculationWorkflow
from xespresso.schedulers import RemoteJobMonitor

# Create workflow
atoms = bulk("Si", cubic=True)

workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials_config='default',
    protocol='moderate',
    machine='snake5'
)

print("=" * 60)
print("Starting non-blocking remote SCF calculation...")
print("=" * 60)

# Run SCF in non-blocking mode (wait_for_completion=False)
# Job will be submitted to remote cluster and return immediately
calc = workflow.run_scf(label='scf/si-test')

print(f"\n✓ Job submitted!")
print(f"  Job ID: {calc.last_job_id}")
print(f"  Remote path: {calc.last_remote_path}")

# Create monitor
print("\n" + "=" * 60)
print("Monitoring job progress...")
print("=" * 60 + "\n")

monitor = RemoteJobMonitor(calc)

# Check job info
print("Job Info:")
for key, value in monitor.info().items():
    print(f"  {key}: {value}")

# Wait for job completion (with 30 second timeout for demo)
print("\nWaiting for job completion...")
if monitor.wait(timeout=300, poll_interval=5):  # 5 min timeout, check every 5s
    print("✓ Job completed!")
    
    # Retrieve output
    local_path, output = monitor.retrieve_output()
    print(f"✓ Output retrieved: {local_path}")
    print(f"✓ Output size: {len(output)} bytes")
    
    # Check convergence
    if "JOB DONE" in output:
        print("✓ Job converged successfully!")
    else:
        print("⚠ Warning: Job completed but may not have converged")
else:
    print("✗ Job timed out or failed")
    status = monitor.status()
    print(f"  Current status: {status}")

print("\n" + "=" * 60)
print("Done!")
print("=" * 60)
