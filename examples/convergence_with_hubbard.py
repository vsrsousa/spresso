"""
Example: Convergence Study with DFT+U (Hubbard) Parameters

This example demonstrates how to perform a convergence study (ecutwfc and kspacing)
while applying Hubbard U parameters for strongly correlated systems.

Two approaches are shown:
1. Using the old format (input_ntyp with Hubbard_U dictionary)
2. Using the new format with explicit orbital specifications (QE >= 7.0)
"""

from ase.build import bulk
import numpy as np
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow

print("="*80)
print("Example 1: Convergence Study with Hubbard (Old Format)")
print("="*80)

# Create Fe structure with antiferromagnetic ordering
atoms = bulk("Fe", cubic=True)
atoms.new_array("species", np.array(atoms.get_chemical_symbols(), dtype="U20"))
atoms.arrays["species"][1] = "Fe1"

# Define Hubbard parameters using old format
hubbard_config_old = {
    "Fe": 4.3,      # U value in eV
    "Fe1": 4.3,
}

# Magnetic configuration (AFM)
magnetic_config_afm = {
    'Fe': {'mag': [0.5]},
    'Fe1': {'mag': [-0.5]}
}

# Create convergence workflow with Hubbard
workflow_old = ConvergenceWorkflow(
    atoms=atoms,
    pseudopotentials={
        'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
        'Fe1': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
    },
    precision='low',
    magnetic_config=magnetic_config_afm,  # Magnetism
    hubbard_config=hubbard_config_old,    # Hubbard U
)

print("\n✓ ConvergenceWorkflow created with:")
print(f"  - Magnetic config: AFM (Fe: 0.5, Fe1: -0.5)")
print(f"  - Hubbard config: Fe/Fe1 U = 4.3 eV (old format)")
print(f"  - Will run convergence study on both ecutwfc and kspacing")

# To run the actual convergence study:
# results = workflow_old.run_convergence_study(
#     label='test/fe_afm_dft+u',
#     max_iterations=3
# )

print("\n" + "="*80)
print("Example 2: Convergence Study with Hubbard (New Format QE >= 7.0)")
print("="*80)

# Using the new format with explicit orbital specifications
hubbard_config_new = {
    "qe_version": "7.2",
    "projector": "atomic",
    "u": {
        "Fe-3d": 4.3,   # Explicit: Fe 3d orbital
        "O-2p": 3.0     # Explicit: O 2p orbital
    }
}

# Create an iron oxide structure
from ase import Atoms
atoms_oxide = Atoms('Fe2O2', 
                    positions=[[0, 0, 0], [1.5, 0, 0], 
                              [0.75, 0.75, 0], [0.75, 0.75, 1.5]],
                    cell=[3, 3, 3])

magnetic_config_fm = {
    'Fe': {'mag': [1, -1]}  # Ferrimagnetic
}

workflow_new = ConvergenceWorkflow(
    atoms=atoms_oxide,
    pseudopotentials={
        'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
        'O': 'O.pbe-n-rrkjus_psl.1.0.0.UPF'
    },
    precision='medium',
    magnetic_config=magnetic_config_fm,
    hubbard_config=hubbard_config_new,  # New format with orbital specs
)

print("\n✓ ConvergenceWorkflow created with:")
print(f"  - Magnetic config: FM (Fe: 1, -1)")
print(f"  - Hubbard config (QE 7.0+):")
print(f"    * Projector: atomic")
print(f"    * Fe-3d: 4.3 eV")
print(f"    * O-2p: 3.0 eV")

# To run the actual convergence study:
# results = workflow_new.run_convergence_study(
#     label='test/fe2o2_dft+u_v7.2',
#     max_iterations=3
# )

print("\n" + "="*80)
print("Example 3: Using setup_magnetic_config Helper (Recommended)")
print("="*80)

from xespresso.utils.magnetic_config_helper import setup_magnetic_config

# Combined magnetic + Hubbard config in one call
atoms_combined = bulk('Fe', cubic=True)

config = setup_magnetic_config(
    atoms_combined,
    {
        'Fe': {
            'mag': [1, -1],
            'U': {'3d': 4.3}  # Can include U in magnetic_config directly
        }
    },
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
    qe_version='7.2'
)

print("\n✓ setup_magnetic_config returns dict with:")
print(f"  - atoms: Modified structure")
print(f"  - input_ntyp: With magnetization and Hubbard_U")
print(f"  - pseudopotentials: Resolved paths")
print(f"  - hubbard: New format parameters (if QE 7.2+)")

# Can pass this directly to ConvergenceWorkflow
workflow_combined = ConvergenceWorkflow(
    atoms=config['atoms'],
    pseudopotentials=config['pseudopotentials'],
    precision='low',
    magnetic_config=config['input_ntyp'],  # Already has Hubbard info
)

print("\n✓ Config passed to ConvergenceWorkflow via magnetic_config")

print("\n" + "="*80)
print("Key Parameters for ConvergenceWorkflow with Hubbard")
print("="*80)

print("""
1. magnetic_config (existing)
   - Can be a dict with:
     * 'type': 'ferromagnetic', 'antiferromagnetic', etc. (optional)
     * Starting magnetization per species
     * Can include 'U' key for Hubbard (old format): {'Fe': {'mag': 1.0, 'U': 4.3}}
   
2. hubbard_config (NEW - now properly forwarded)
   - Option A (old format):
     hubbard_config = {'Fe': 4.3, 'Mn': 5.7}
   
   - Option B (new format for QE >= 7.0):
     hubbard_config = {
         'qe_version': '7.2',
         'projector': 'atomic',
         'u': {'Fe-3d': 4.3, 'Mn-3d': 5.7},
         'v': [...]  # Optional V parameters (inter-site)
     }

Example workflow call:

from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow

atoms = bulk('Fe', cubic=True)

workflow = ConvergenceWorkflow(
    atoms=atoms,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
    precision='low',
    
    # Magnetic configuration (AFM)
    magnetic_config={
        'type': 'antiferromagnetic',
        'Fe': {'mag': 1.0, 'opposite': True}
    },
    
    # Hubbard U parameters (DFT+U)
    hubbard_config={
        'Fe': 4.3  # Old format
    },
    
    # Optional: machine/queue settings
    queue={'nodes': 1, 'ntasks-per-node': 12}
)

# Run 2-phase convergence (ecutwfc -> kspacing)
results = workflow.run_convergence_study(
    label='test/fe_dft+u_afm',
    max_iterations=3
)
""")

print("\n✓ All examples show different ways to use Hubbard with ConvergenceWorkflow")
print("✓ Hubbard_config is now properly forwarded to all CalculationWorkflow instances")
