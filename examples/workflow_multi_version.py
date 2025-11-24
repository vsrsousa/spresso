#!/usr/bin/env python
"""
Practical Example: Using multiple QE versions in a workflow

This example demonstrates a realistic scenario where:
- You need QE 7.2 for SCF calculations (newer features)
- You need QE 6.8 for DOS/bands (stability for post-processing)
- All running on the same remote cluster
- Connection persists throughout the workflow
"""

from ase.build import bulk
from xespresso.codes import load_codes_config
from xespresso.machines import load_machine

# This example assumes you have:
# 1. A machine configured named "cluster1"
# 2. Codes configured for QE 7.2 and 6.8 on that machine
# 
# To set up:
#   from xespresso.codes import add_version_to_config
#   
#   add_version_to_config(
#       machine_name="cluster1",
#       version="7.2",
#       qe_prefix="/opt/qe-7.2/bin",
#       modules=["quantum-espresso/7.2"]
#   )
#   
#   add_version_to_config(
#       machine_name="cluster1",
#       version="6.8",
#       qe_prefix="/opt/qe-6.8/bin",
#       modules=["quantum-espresso/6.8"]
#   )

print("="*70)
print("Practical Workflow: Multiple QE Versions")
print("="*70)

# Load machine and codes configuration ONCE
print("\n1. Loading machine and codes configuration...")
try:
    machine = load_machine('cluster1')
    codes = load_codes_config('cluster1')
    
    print(f"   ✓ Machine loaded: {machine.get('execution', 'N/A')} execution")
    print(f"   ✓ Available QE versions: {', '.join(codes.list_versions())}")
except Exception as e:
    print(f"\n❌ Error loading configuration: {e}")
    print("\nNote: This example requires:")
    print("1. A machine configuration named 'cluster1'")
    print("2. QE codes configured for versions 7.2 and 6.8")
    print("\nSee the header of this file for setup instructions.")
    exit(1)

# Build structure
atoms = bulk('Fe', 'bcc', a=2.87)

# Task 1: SCF calculation with QE 7.2
print("\n2. SCF Calculation (using QE 7.2)...")
print("   Reason: Using newer version for improved SCF convergence")

pw_code_72 = codes.get_code('pw', version='7.2')
if pw_code_72:
    print(f"   ✓ Selected: {pw_code_72.path}")
    
    # Get version-specific modules
    ver_config = codes.get_version_config('7.2')
    modules = ver_config.get('modules', []) if ver_config else []
    print(f"   ✓ Modules to load: {', '.join(modules)}")
    
    # In a real scenario, you would run:
    # calc = Espresso(
    #     atoms=atoms,
    #     command=pw_code_72.path,
    #     queue=machine,
    #     ...
    # )
    # atoms.calc = calc
    # energy = atoms.get_potential_energy()
    
    print("   ✓ SCF would run here...")
    print("   📡 SSH connection to cluster1 established")
else:
    print("   ⚠️  QE 7.2 pw.x not found")

# Task 2: DOS calculation with QE 6.8
print("\n3. DOS Calculation (using QE 6.8)...")
print("   Reason: Using stable version for reliable post-processing")

dos_code_68 = codes.get_code('dos', version='6.8')
if dos_code_68:
    print(f"   ✓ Selected: {dos_code_68.path}")
    
    # Get version-specific modules
    ver_config = codes.get_version_config('6.8')
    modules = ver_config.get('modules', []) if ver_config else []
    print(f"   ✓ Modules to load: {', '.join(modules)}")
    
    print("   ✓ DOS would run here...")
    print("   🔄 REUSING SSH connection to cluster1 (no reconnect!)")
else:
    print("   ⚠️  QE 6.8 dos.x not found")

# Task 3: Bands calculation with QE 6.8
print("\n4. Bands Calculation (using QE 6.8)...")
print("   Reason: Same version for consistency in post-processing")

bands_code_68 = codes.get_code('bands', version='6.8')
if bands_code_68:
    print(f"   ✓ Selected: {bands_code_68.path}")
    print("   ✓ Bands would run here...")
    print("   🔄 STILL using same SSH connection!")
else:
    print("   ⚠️  QE 6.8 bands.x not found")

# Task 4: Another SCF with QE 7.2 (different system)
print("\n5. Another SCF Calculation (back to QE 7.2)...")
print("   Reason: Different structure, using newer version again")

pw_code_72_again = codes.get_code('pw', version='7.2')
if pw_code_72_again:
    print(f"   ✓ Selected: {pw_code_72_again.path}")
    print("   ✓ SCF would run here...")
    print("   🔄 Connection persistence maintained!")
else:
    print("   ⚠️  QE 7.2 pw.x not found")

print("\n" + "="*70)
print("Key Benefits Demonstrated:")
print("="*70)
print("""
✅ Version Flexibility:
   - Use QE 7.2 for SCF (newer features)
   - Use QE 6.8 for post-processing (stability)
   - Switch seamlessly between versions

✅ Connection Efficiency:
   - Single SSH connection throughout
   - No overhead when switching versions
   - Connection identified by (host, username)

✅ Module Management:
   - Each version has its own modules
   - Automatically loaded per calculation
   - No manual module switching needed

✅ Code Organization:
   - Clear version selection per task
   - Easy to maintain and understand
   - Prevents version conflicts
""")

print("\n" + "="*70)
print("Real-World Use Cases:")
print("="*70)
print("""
1. Testing/Validation:
   - Compare results between QE versions
   - Validate new features in QE 7.x
   - Ensure backward compatibility

2. Feature Requirements:
   - Use QE 7.x for Hubbard parameters (new format)
   - Use QE 6.x for legacy input compatibility
   - Mix as needed in same workflow

3. Stability Considerations:
   - Use stable version for production runs
   - Use newer version for exploratory calculations
   - Transition gradually between versions

4. HPC Environment:
   - Multiple QE versions installed system-wide
   - Different users need different versions
   - Easy version selection per job
""")

print("\n✅ Example completed!")
print("\nFor more information:")
print("  - docs/CODES_CONFIGURATION.md - Full documentation")
print("  - examples/multi_version_example.py - Setup examples")
print("  - docs/REMOTE_CONNECTION_PERSISTENCE.md - Connection details")
