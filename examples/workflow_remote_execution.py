"""
Example demonstrating remote job execution with the simplified workflow.

This example shows how to:
1. Run workflows on remote machines using SLURM
2. Use machine configurations from ~/.xespresso/machines/
3. Use direct queue configuration
4. Combine remote execution with quality presets and magnetic configurations
"""

from ase.build import bulk
from xespresso import CalculationWorkflow, quick_scf, quick_relax

# ============================================================================
# Example 1: Remote execution with machine configuration
# ============================================================================
print("=" * 60)
print("Example 1: Remote execution using machine configuration")
print("=" * 60)

atoms = bulk("Fe", cubic=True)
pseudopotentials = {"Fe": "Fe.pbe-spn.UPF"}

# Option 1: Load machine configuration from ~/.xespresso/machines/cluster1.json
# This assumes you have created a machine configuration file
"""
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials=pseudopotentials,
    protocol='moderate',
    machine='cluster1'  # Loads from ~/.xespresso/machines/cluster1.json
)

# Run SCF calculation on remote cluster
calc = workflow.run_scf(label='scf/fe-remote')
print(f"Energy: {calc.results['energy']} eV")
"""

print("To use this example:")
print("1. Create a machine config in ~/.xespresso/machines/cluster1.json")
print("2. Uncomment the code above")
print("3. Run the script")

# ============================================================================
# Example 2: Remote execution with direct queue configuration
# ============================================================================
print("\n" + "=" * 60)
print("Example 2: Remote execution with direct queue configuration")
print("=" * 60)

atoms = bulk("Si", cubic=True)
pseudopotentials = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}

# Define queue configuration directly
queue_config = {
    "execution": "remote",
    "scheduler": "slurm",
    "remote_host": "cluster.university.edu",
    "remote_user": "username",
    "remote_dir": "/home/username/calculations",
    "remote_auth": {
        "method": "key",
        "ssh_key": "~/.ssh/id_rsa"
    },
    "nodes": 1,
    "ntasks-per-node": 16,
    "time": "02:00:00",
    "partition": "compute"
}

"""
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials=pseudopotentials,
    protocol='moderate',
    queue=queue_config
)

# Run relaxation on remote cluster
calc = workflow.run_relax(label='relax/si-remote', relax_type='vc-relax')
"""

print("Direct queue configuration allows full control over job submission")

# ============================================================================
# Example 3: Local execution with SLURM (for testing)
# ============================================================================
print("\n" + "=" * 60)
print("Example 3: Local execution with SLURM scheduler")
print("=" * 60)

atoms = bulk("Al", cubic=True)
pseudopotentials = {"Al": "Al.pbe.UPF"}

# Local SLURM execution (if SLURM is available on local machine)
local_slurm_queue = {
    "execution": "local",
    "scheduler": "slurm",
    "nodes": 1,
    "ntasks-per-node": 4,
    "time": "00:30:00",
    "partition": "debug"
}

"""
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials=pseudopotentials,
    protocol='fast',
    queue=local_slurm_queue
)

calc = workflow.run_scf(label='scf/al-local-slurm')
"""

print("This is useful for testing SLURM scripts locally before remote submission")

# ============================================================================
# Example 4: Quick functions with remote execution
# ============================================================================
print("\n" + "=" * 60)
print("Example 4: Using quick functions with remote execution")
print("=" * 60)

atoms = bulk("Fe", cubic=True)
pseudopotentials = {"Fe": "Fe.pbe-spn.UPF"}

# Quick SCF with remote execution and magnetic configuration
"""
calc = quick_scf(
    atoms,
    pseudopotentials,
    protocol='moderate',
    magnetic_config='ferro',
    machine='cluster1',  # Remote execution
    label='scf/fe-ferro-remote'
)
"""

# Quick relaxation with remote execution
"""
calc = quick_relax(
    'structure.cif',
    {'Fe': 'Fe.pbe-spn.UPF'},
    protocol='accurate',
    relax_type='vc-relax',
    machine='cluster1',  # Remote execution
    label='relax/fe-remote'
)
"""

print("Quick functions make remote execution as simple as local execution")

# ============================================================================
# Example 5: Magnetic calculation with Hubbard U on remote cluster
# ============================================================================
print("\n" + "=" * 60)
print("Example 5: Complex magnetic + Hubbard calculation on remote")
print("=" * 60)

atoms = bulk("Fe", cubic=True)
pseudopotentials = {"Fe": "Fe.pbe-spn.UPF", "O": "O.pbe.UPF"}

# Combine all features: quality presets, magnetic config, Hubbard, remote execution
"""
calc = quick_scf(
    atoms,
    pseudopotentials,
    protocol='accurate',
    magnetic_config={
        'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}  # AFM with Hubbard U
    },
    machine='cluster1',  # Remote execution
    label='scf/fe-hubbard-remote'
)
"""

print("All workflow features work seamlessly with remote execution")

# ============================================================================
# Example 6: Creating and using a machine configuration
# ============================================================================
print("\n" + "=" * 60)
print("Example 6: Creating a machine configuration programmatically")
print("=" * 60)

from xespresso.machines import Machine

# Create a machine configuration
"""
machine = Machine(
    name="my_cluster",
    execution="remote",
    scheduler="slurm",
    workdir="/home/user/calculations",
    host="cluster.university.edu",
    username="user",
    auth={"method": "key", "ssh_key": "~/.ssh/id_rsa"},
    nprocs=32,
    resources={
        "nodes": 2,
        "ntasks-per-node": 16,
        "time": "04:00:00",
        "partition": "compute"
    }
)

# Save for future use
machine.to_file("~/.xespresso/machines/my_cluster.json")

# Use it in workflow
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials=pseudopotentials,
    protocol='moderate',
    machine='my_cluster'
)

calc = workflow.run_scf(label='scf/test')
"""

print("Machine configurations can be created and saved programmatically")

# ============================================================================
# Important Notes
# ============================================================================
print("\n" + "=" * 60)
print("Important Notes about Remote Execution")
print("=" * 60)

print("""
1. **Connection Persistence**: 
   - SSH connections are automatically cached and reused
   - Multiple calculations on the same machine use the same connection
   - No need to manually manage connections

2. **Machine Configuration Files**:
   - Store machine configs in ~/.xespresso/machines/
   - Each machine has its own JSON file
   - See docs/MACHINE_CONFIGURATION.md for details

3. **Queue Configuration**:
   - Use 'machine' parameter to load from config file
   - Use 'queue' parameter for direct configuration
   - Cannot use both 'machine' and 'queue' together

4. **Remote Directories**:
   - Input files are automatically transferred
   - Pseudopotentials are automatically transferred
   - Output files are retrieved after calculation

5. **SLURM Integration**:
   - Job submission is automatic
   - Job status is monitored
   - Output is retrieved when job completes

6. **Error Handling**:
   - File transfer failures are caught and retried
   - Connection errors trigger new connection attempts
   - Missing files raise clear error messages

For complete documentation, see:
- WORKFLOW_DOCUMENTATION.md - Workflow overview
- docs/MACHINE_CONFIGURATION.md - Machine configuration
- docs/REMOTE_CONNECTION_PERSISTENCE.md - Connection management
""")

print("\n" + "=" * 60)
print("Examples complete!")
print("=" * 60)
