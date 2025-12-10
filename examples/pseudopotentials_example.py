"""
Example demonstrating the pseudopotentials module for xespresso.

This example shows how to:
1. Auto-detect pseudopotentials from a directory
2. Create and save configurations
3. Load and use configurations in calculations
4. Manage multiple pseudopotential libraries
"""

from xespresso.pseudopotentials import (
    create_pseudopotentials_config,
    load_pseudopotentials_config,
    PseudopotentialsManager
)

print("=" * 70)
print("Pseudopotentials Module - Example Usage")
print("=" * 70)

# ============================================================================
# Example 1: Configure local pseudopotentials
# ============================================================================
print("\n" + "=" * 70)
print("Example 1: Configure local pseudopotentials (SSSP)")
print("=" * 70)

# After downloading SSSP pseudopotentials to your local computer:
# (Replace with your actual path)
sssp_path = "/path/to/SSSP_1.1.2_PBE_efficiency"

print(f"""
To create a configuration for SSSP pseudopotentials:

config = create_pseudopotentials_config(
    name="SSSP_efficiency",
    base_path="{sssp_path}",
    description="SSSP 1.1.2 PBE efficiency pseudopotentials",
    functional="PBE",
    library="SSSP",
    version="1.1.2",
    save=True
)
""")

# ============================================================================
# Example 2: Configure pseudopotentials from PSLibrary
# ============================================================================
print("\n" + "=" * 70)
print("Example 2: Configure PSLibrary pseudopotentials")
print("=" * 70)

pslibrary_path = "/path/to/pslibrary/pbe"

print(f"""
For PSLibrary pseudopotentials:

config = create_pseudopotentials_config(
    name="PSLibrary_pbe",
    base_path="{pslibrary_path}",
    description="PSLibrary PBE pseudopotentials",
    functional="PBE",
    library="PSLibrary",
    version="1.0.0",
    save=True
)
""")

# ============================================================================
# Example 3: Configure pseudopotentials on a remote machine
# ============================================================================
print("\n" + "=" * 70)
print("Example 3: Configure pseudopotentials on a remote machine")
print("=" * 70)

print("""
If pseudopotentials are stored on a remote machine (less common):

config = create_pseudopotentials_config(
    name="remote_sssp",
    base_path="/opt/pseudopotentials/SSSP",
    machine_name="my_cluster",  # Reference to a configured machine
    description="SSSP pseudopotentials on remote cluster",
    functional="PBE",
    library="SSSP",
    version="1.1.2",
    save=True
)

Note: The machine must be configured in xespresso first.
""")

# ============================================================================
# Example 4: Load and use configuration in calculations
# ============================================================================
print("\n" + "=" * 70)
print("Example 4: Load and use pseudopotentials in calculations")
print("=" * 70)

print("""
# Load a saved configuration
config = load_pseudopotentials_config("SSSP_efficiency")

# Get pseudopotentials as a dictionary for use in calculations
pseudopotentials = config.get_pseudopotentials_dict()
# Returns: {'Fe': 'Fe.pbe-spn-kjpaw_psl.0.2.1.UPF', 'O': '...', ...}

# Use in xespresso calculation
from xespresso import Espresso
from ase.build import bulk
import os

atoms = bulk('Fe', 'bcc', a=2.87)

# Set ESPRESSO_PSEUDO environment variable to config's base_path
# This allows xespresso to find files without explicitly setting pseudo_dir
os.environ['ESPRESSO_PSEUDO'] = config.base_path

calc = Espresso(
    atoms=atoms,
    pseudopotentials=pseudopotentials,
    # No need to set pseudo_dir - xespresso uses ESPRESSO_PSEUDO
    input_data={'ecutwfc': 50.0},
    # ... other calculation parameters
)

# This approach:
# - Works for both LOCAL and REMOTE execution
# - Respects xespresso's original implementation
# - For remote: xespresso transfers files to ./pseudo automatically
""")

# ============================================================================
# Example 5: List and manage configurations
# ============================================================================
print("\n" + "=" * 70)
print("Example 5: List and manage configurations")
print("=" * 70)

print("""
# List all saved configurations
configs = PseudopotentialsManager.list_configs()
print(f"Available configurations: {configs}")

# Load a specific configuration
config = load_pseudopotentials_config("SSSP_efficiency")

# Check available elements
elements = config.list_elements()
print(f"Elements: {elements}")

# Get info for a specific element
fe_pseudo = config.get_pseudopotential("Fe")
print(f"Fe pseudopotential: {fe_pseudo.filename}")
print(f"Type: {fe_pseudo.type}")
print(f"Functional: {fe_pseudo.functional}")

# Delete a configuration
PseudopotentialsManager.delete_config("old_config")
""")

# ============================================================================
# Example 6: Multiple library versions
# ============================================================================
print("\n" + "=" * 70)
print("Example 6: Managing multiple library versions")
print("=" * 70)

print("""
You can configure multiple versions of the same library:

# SSSP 1.1.2 efficiency
create_pseudopotentials_config(
    name="SSSP_1.1.2_efficiency",
    base_path="/path/to/SSSP_1.1.2_PBE_efficiency",
    library="SSSP",
    version="1.1.2",
    functional="PBE"
)

# SSSP 1.1.2 precision
create_pseudopotentials_config(
    name="SSSP_1.1.2_precision",
    base_path="/path/to/SSSP_1.1.2_PBE_precision",
    library="SSSP",
    version="1.1.2",
    functional="PBE"
)

# SSSP 1.2.0 (future version)
create_pseudopotentials_config(
    name="SSSP_1.2.0_efficiency",
    base_path="/path/to/SSSP_1.2.0_PBE_efficiency",
    library="SSSP",
    version="1.2.0",
    functional="PBE"
)
""")

# ============================================================================
# Example 7: GUI Configuration
# ============================================================================
print("\n" + "=" * 70)
print("Example 7: Using the GUI for configuration")
print("=" * 70)

print("""
For a more user-friendly approach, use the Qt GUI:

1. Start the xespresso Qt GUI:
   $ python -m qtgui

2. Navigate to the Pseudopotentials Configuration page

3. Fill in the form:
   - Configuration Name: e.g., "SSSP_efficiency"
   - Pseudopotentials Directory: path to .UPF files
   - Functional: e.g., "PBE"
   - Library Name: e.g., "SSSP"
   - Library Version: e.g., "1.1.2"
   - Description: Brief description

5. Click "🔍 Auto-Detect Pseudopotentials"

6. Review detected pseudopotentials and click "💾 Save"

The configuration will be saved to ~/.xespresso/pseudopotentials/
""")

# ============================================================================
# Example 8: Best practices
# ============================================================================
print("\n" + "=" * 70)
print("Example 8: Best practices")
print("=" * 70)

print("""
Best Practices:
---------------
1. **Download locally**: Keep pseudopotentials on your local computer
   - Easier to manage and update
   - xespresso copies them to remote machines automatically

2. **Use descriptive names**: Include library, version, and type in config name
   - Good: "SSSP_1.1.2_PBE_efficiency"
   - Bad: "pseudos"

3. **Track versions**: Always specify library version
   - Makes it easier to reproduce calculations
   - Important for scientific reproducibility

4. **Organize by functional**: Keep different functionals separate
   - Create separate configs for PBE, LDA, PBEsol, etc.

5. **Document your choices**: Use the description field
   - Explain why you chose this particular set
   - Note any special considerations

6. **Set ESPRESSO_PSEUDO**: When using configs programmatically
   - os.environ['ESPRESSO_PSEUDO'] = config.base_path
   - Allows xespresso to find files without setting pseudo_dir
   - Works with xespresso's original implementation

7. **Use GUI**: For interactive work, the GUI automatically sets ESPRESSO_PSEUDO
   - Simplifies pseudopotential selection
   - Handles environment variable management

8. **Test first**: Before production calculations
   - Test with a small calculation
   - Verify elements are correctly detected
   - Check that paths are accessible
""")


print("\n" + "=" * 70)
print("For more information:")
print("  - Documentation: https://github.com/vsrsousa/xespresso")
print("  - SSSP Library: https://www.materialscloud.org/sssp")
print("  - PSLibrary: https://dalcorso.github.io/pslibrary/")
print("=" * 70)
