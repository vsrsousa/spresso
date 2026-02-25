#!/usr/bin/env python
"""
Example: CalculationWorkflow with Machine + Pseudopotentials + Codes (Versioned)

This demonstrates the complete integration where the user can specify:
1. A machine (cluster, scheduler, resources)
2. A pseudopotential configuration (elements, filenames)
3. A specific QE version to use on that machine
4. The workflow automatically handles the rest!
"""

from ase.build import bulk
from xespresso import CalculationWorkflow

print("=" * 80)
print("EXAMPLE: Full Integration - Machine + Pseudos + Code Version")
print("=" * 80)

# Example 1: Simple case with explicit queue
print("\n1️⃣  EXAMPLE 1: Local execution with explicit pseudos")
print("-" * 80)

print("""
workflow = CalculationWorkflow(
    atoms=bulk("Fe", cubic=True),
    protocol='moderate',
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'}
)
""")

try:
    workflow1 = CalculationWorkflow(
        atoms=bulk("Fe", cubic=True),
        protocol='moderate',
        pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'}
    )
    print("✓ Created workflow with explicit pseudotypes") 
except Exception as e:
    print(f"✗ Error: {e}")

# Example 2: Using pseudopotentials_config
print("\n2️⃣  EXAMPLE 2: Using pseudopotentials_config (auto-extract)")
print("-" * 80)

print("""
workflow = CalculationWorkflow(
    atoms=bulk("Fe", cubic=True),
    protocol='moderate',
    pseudopotentials_config='SSSP_efficiency'
    # Automatically extracts Fe pseudopotential from config!
)
""")

try:
    workflow2 = CalculationWorkflow(
        atoms=bulk("Fe", cubic=True),
        protocol='moderate',
        pseudopotentials_config='SSSP_efficiency'
    )
    print(f"✓ Created workflow with config 'SSSP_efficiency'")
    print(f"  Status: Config loading works (actual values depend on ~/.xespresso/)")
except ValueError as e:
    print(f"ℹ Config not found (expected): {str(e)[:70]}...")

# Example 3: Complete setup: machine + pseudos + code version
print("\n3️⃣  EXAMPLE 3: Remote machine + pseudos + QE version (THE FULL POWER!)")
print("-" * 80)

print("""
workflow = CalculationWorkflow(
    atoms=bulk("Fe", cubic=True),
    protocol='moderate',
    machine='cluster1',                    # ← Load machine config
    pseudopotentials_config='SSSP_efficiency',  # ← Auto-extract pseudos
    code_version='7.2'                     # ← Auto-load modules for v7.2!
)

This workflow automatically:
1. Loads machine='cluster1' → queue with scheduler, resources, host, SSH settings
2. Loads pseudos='SSSP_efficiency' → extracts Fe.pbe-spn.UPF automatically
3. Loads code_version='7.2' → extracts modules like 'quantum-espresso/7.2'
4. Merges everything into a single queue dictionary

Result:
  queue = {
    'execution': 'remote',
    'scheduler': 'slurm',
    'host': 'cluster.edu',
    'modules': ['quantum-espresso/7.2'],  # ← From code version 7.2
    'resources': {'nodes': 2, 'time': '02:00:00'},
    ...
  }

and pseudopotentials = {'Fe': 'Fe.pbe-spn.UPF'}
""")

# Example 4: Different QE versions on same machine
print("\n4️⃣  EXAMPLE 4: Switching QE versions on same machine")
print("-" * 80)

print("""
# Machine 'snake5' has TWO versions of QE: 7.2 and 6.8
# Each with different modules:
#   7.2: ['quantum-espresso/7.2']
#   6.8: ['quantum-espresso/6.8']

# Use with version 7.2
workflow_v72 = CalculationWorkflow(
    atoms=my_atoms,
    machine='snake5',
    pseudopotentials_config='SSSP_efficiency',
    code_version='7.2'  # ← Load 7.2 specific modules
)

# Use with version 6.8
workflow_v68 = CalculationWorkflow(
    atoms=my_atoms,
    machine='snake5',
    pseudopotentials_config='SSSP_efficiency',
    code_version='6.8'  # ← Load 6.8 specific modules
)

Both use same machine,  but different QE versions!
""")

# Example 5: Configuration separation diagram
print("\n5️⃣  CONFIGURATION FILES (Separate, Independent)")
print("-" * 80)

print("""
~/.xespresso/
├── machines/
│   └── cluster1.json          ← Scheduler, host, resources, SSH
├── pseudopotentials/
│   └── SSSP_efficiency.json   ← Elements, filenames, functional
└── codes/
    └── cluster1.json          ← Versions, modules per version

Example: cluster1 codes config with 2 versions:
{
  "machine_name": "cluster1",
  "versions": {
    "7.2": {
      "qe_prefix": "/opt/qe-7.2/bin",
      "modules": ["quantum-espresso/7.2"],
      "codes": {"pw": {...}, "ph": {...}}
    },
    "6.8": {
      "qe_prefix": "/opt/qe-6.8/bin",
      "modules": ["quantum-espresso/6.8"],
      "codes": {"pw": {...}, "ph": {...}}
    }
  }
}
""")

# Example 6: Benefits
print("\n6️⃣  KEY BENEFITS")
print("-" * 80)

print("""
✅ SEPARATION OF CONCERNS
   - machine_config.json ≠ pseudopotential_config.json ≠ codes_config.json
   - Can use machine without changing pseudo/codes

✅ AUTOMATIC MODULE LOADING
   - User just specifies code_version
   - Workflow extracts and adds modules to queue
   - No manual SSH to check module names!

✅ SINGLE MACHINE, MULTIPLE QE VERSIONS
   - Same cluster with QE 7.2 and 6.8
   - Just change code_version parameter
   - No need to reconfigure machine

✅ AUTO-EXTRACTION OF PSEUDOPOTENTIALS
   - Pass config name, get only needed elements
   - Fail early if element not in config
   - No manual mapping for each structure

✅ USER-FRIENDLY API
   - Before: need to know everything (pseudos, modules, scheduler)
   - After: just pass configuration NAMES
   -All details extracted automatically!

Example of power:
   
   BEFORE (tedious):
   workflow = CalculationWorkflow(
       atoms=my_structure,
       protocol='moderate',
       pseudopotentials={
           'Fe': 'Fe.pbe-spn.UPF',
           'O': 'O.pbe-n.UPF',
       },
       queue={
           'execution': 'remote',
           'scheduler': 'slurm',
           'host': 'snake5.df.ufscar.br',
           'modules': ['quantum-espresso/7.2'],  # Had to know this!
           'resources': {'nodes': 2, 'time': '02:00:00'},
       }
   )
   
   AFTER (clean):
   workflow = CalculationWorkflow(
       atoms=my_structure,
       protocol='moderate',
       machine='snake5',
       pseudopotentials_config='SSSP_efficiency',
       code_version='7.2'
   )
   # Everything else is automatically handled!
""")

print("\n" + "=" * 80)
print("✅ EXAMPLE COMPLETE")
print("=" * 80)
