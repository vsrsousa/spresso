#!/usr/bin/env python
"""
Example: Managing multiple Quantum ESPRESSO versions on the same machine

This example demonstrates:
1. Setting up multiple QE versions on a single machine
2. Switching between versions for different calculations
3. Connection persistence when using different versions
"""

import tempfile
import os
from xespresso.codes import (
    CodesConfig, Code, CodesManager,
    add_version_to_config, load_codes_config
)

print("="*70)
print("Example: Multiple QE Versions on Same Machine")
print("="*70)

# Create a temporary directory for this example
tmpdir = tempfile.mkdtemp()
print(f"\n📁 Using temporary directory: {tmpdir}")

# Scenario: A machine has both QE 7.2 and QE 6.8 installed
# This is common in HPC clusters where multiple versions coexist

print("\n" + "="*70)
print("Step 1: Create configuration with QE 7.2")
print("="*70)

config = CodesConfig(
    machine_name="cluster1",
    qe_version="7.2"  # Set as default version
)

# Add QE 7.2 codes
print("\nAdding QE 7.2 codes...")
for code_name, executable in [
    ("pw", "/opt/qe-7.2/bin/pw.x"),
    ("hp", "/opt/qe-7.2/bin/hp.x"),
    ("dos", "/opt/qe-7.2/bin/dos.x"),
]:
    code = Code(name=code_name, path=executable, version="7.2")
    config.add_code(code, version="7.2")
    print(f"  ✓ Added {code_name}.x")

# Add version-specific settings for QE 7.2
if not config.versions:
    config.versions = {}
config.versions["7.2"] = {
    "qe_prefix": "/opt/qe-7.2/bin",
    "modules": ["quantum-espresso/7.2"],
    "codes": config.versions.get("7.2", {}).get("codes", {})
}

print("\n" + "="*70)
print("Step 2: Add QE 6.8 to the same machine")
print("="*70)

# Add QE 6.8 codes
print("\nAdding QE 6.8 codes...")
for code_name, executable in [
    ("pw", "/opt/qe-6.8/bin/pw.x"),
    ("hp", "/opt/qe-6.8/bin/hp.x"),
    ("dos", "/opt/qe-6.8/bin/dos.x"),
]:
    code = Code(name=code_name, path=executable, version="6.8")
    config.add_code(code, version="6.8")
    print(f"  ✓ Added {code_name}.x")

# Add version-specific settings for QE 6.8
config.versions["6.8"] = {
    "qe_prefix": "/opt/qe-6.8/bin",
    "modules": ["quantum-espresso/6.8"],
    "codes": config.versions.get("6.8", {}).get("codes", {})
}

# Save configuration
filepath = CodesManager.save_config(config, output_dir=tmpdir)
print(f"\n💾 Configuration saved to: {filepath}")

print("\n" + "="*70)
print("Step 3: Load and inspect configuration")
print("="*70)

loaded_config = load_codes_config("cluster1", codes_dir=tmpdir)

if loaded_config:
    print("\n📋 Configuration details:")
    print(f"   Machine: {loaded_config.machine_name}")
    print(f"   Default version: {loaded_config.qe_version}")
    print(f"   Available versions: {', '.join(loaded_config.list_versions())}")
    
    # Show codes for each version
    for version in loaded_config.list_versions():
        print(f"\n   QE {version} codes:")
        codes = loaded_config.list_codes(version=version)
        for code_name in codes:
            code = loaded_config.get_code(code_name, version=version)
            print(f"     • {code_name}.x → {code.path}")

print("\n" + "="*70)
print("Step 4: Using different versions in calculations")
print("="*70)

print("""
In your calculation script, you can now choose which version to use:

```python
from xespresso.codes import load_codes_config
from xespresso.machines import load_machine

# Load machine and codes configuration
machine = load_machine('cluster1')
codes = load_codes_config('cluster1')

# Calculation 1: Use QE 7.2 (default)
pw_code_72 = codes.get_code('pw', version='7.2')
print(f"Using pw.x from QE 7.2: {pw_code_72.path}")

# Calculation 2: Use QE 6.8
pw_code_68 = codes.get_code('pw', version='6.8')
print(f"Using pw.x from QE 6.8: {pw_code_68.path}")

# The SSH connection to cluster1 is maintained throughout!
# No need to reconnect when switching versions.
```
""")

print("\n" + "="*70)
print("Step 5: Connection Persistence")
print("="*70)

print("""
✅ Key Benefit: Connection Persistence

When you run multiple calculations on the same machine using different
QE versions, the SSH connection is AUTOMATICALLY maintained:

1. First calculation with QE 7.2:
   - Connects to cluster1
   - Loads quantum-espresso/7.2 module
   - Runs pw.x from /opt/qe-7.2/bin/pw.x

2. Second calculation with QE 6.8:
   - REUSES existing connection to cluster1
   - Loads quantum-espresso/6.8 module
   - Runs pw.x from /opt/qe-6.8/bin/pw.x

3. Third calculation back to QE 7.2:
   - STILL using same connection
   - Just switches modules and executable path

Benefits:
  ⚡ Faster job submission (no SSH handshake)
  🎯 More reliable (established connection)
  💪 Less server load (fewer connections)
  🔄 Seamless version switching

The connection is identified by (host, username), not by QE version!
""")

print("\n" + "="*70)
print("Step 6: Programmatic Version Addition")
print("="*70)

print("""
You can also add versions programmatically:

```python
from xespresso.codes import add_version_to_config

# Add QE 7.3 to existing configuration
config = add_version_to_config(
    machine_name="cluster1",
    version="7.3",
    qe_prefix="/opt/qe-7.3/bin",
    modules=["quantum-espresso/7.3"]
)

# For remote machines, provide SSH connection
config = add_version_to_config(
    machine_name="remote_cluster",
    version="7.2",
    qe_prefix="/software/qe-7.2/bin",
    modules=["qe/7.2"],
    ssh_connection={
        'host': 'cluster.edu',
        'username': 'user'
    }
)
```
""")

# Cleanup
import shutil
shutil.rmtree(tmpdir, ignore_errors=True)

print("\n" + "="*70)
print("✅ Example completed successfully!")
print("="*70)
print("\nNext steps:")
print("1. See docs/CODES_CONFIGURATION.md for detailed documentation")
print("2. See examples/codes_example.py for more examples")
print("3. Use xespresso.codes API in your calculation scripts")
print()
