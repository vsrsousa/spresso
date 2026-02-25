"""
Example: Using CalculationWorkflow with remote job monitoring.

Shows how to:
1. Run remote SCF calculation in non-blocking mode
2. Get a monitor via the workflow
3. Track job progress
4. Retrieve results when done
"""

from ase.build import bulk
from xespresso import CalculationWorkflow

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
print("=" * 60 + "\n")

# Run SCF in non-blocking mode
calc = workflow.run_scf(label='scf/si-test')

print(f"✓ Job submitted!")
print(f"  Job ID: {calc.last_job_id}")
print(f"  Remote path: {calc.last_remote_path}\n")

# Get monitor from workflow
print("=" * 60)
print("Monitoring job via workflow.get_monitor()...")
print("=" * 60 + "\n")

monitor = workflow.get_monitor()

# Check info
print("Job Info:")
info = monitor.info()
for key, value in info.items():
    print(f"  {key}: {value}\n")

# Wait for completion
print("Waiting for job completion (timeout: 5 min)...")
if monitor.wait(timeout=300, poll_interval=10):
    print("\n✓ Job completed!\n")
    
    # Retrieve output
    try:
        local_path, output = monitor.retrieve_output()
        print(f"✓ Output retrieved: {local_path}")
        
        if "JOB DONE" in output:
            print("✓ Job converged successfully!")
            
            # Read results
            energy = calc.results.get('energy', None)
            if energy:
                print(f"✓ Final energy: {energy:.6f} eV")
        else:
            print("⚠ Warning: Job finished but may not have converged")
    except Exception as e:
        print(f"✗ Could not retrieve output: {e}")
else:
    status = monitor.status()
    print(f"\n✗ Job timed out or failed")
    print(f"  Current status: {status}")

print("\n" + "=" * 60)
print("Done!")
print("=" * 60)
