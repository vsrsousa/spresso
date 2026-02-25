#!/usr/bin/env python
"""
Example: Using CalculationWorkflow with pseudopotentials_config

This demonstrates the new simplified workflow where pseudopotentials are
automatically extracted from a configuration based on the elements in your structure.

Key advantages:
1. Just pass config name - workflow extracts needed pseudopotentials
2. No need to manually specify every element
3. Separation of concerns: machine config ≠ pseudopotential config
4. User only passes relevant names, everything else is automatic
"""

from ase.build import bulk
from ase.io import write
import tempfile
from xespresso import CalculationWorkflow
from xespresso.pseudopotentials import create_pseudopotentials_config, load_pseudopotentials_config

print("=" * 80)
print("EXAMPLE: CalculationWorkflow with Automatic Pseudopotential Extraction")
print("=" * 80)

# =============================================================================
# 1. SETUP: Create a pseudopotentials configuration
# =============================================================================
print("\n1️⃣  SETUP: Create a pseudopotentials configuration")
print("-" * 80)

print("\nIn a real workflow, you would:")
print("  a) Download pseudopotentials from SSSP, PSLibrary, etc.")
print("  b) Save them in a directory (e.g., ~/pseudopotentials/SSSP_v1.1.2_PBE/)")
print("  c) Run: create_pseudopotentials_config()")
print("     to scan the directory and save a .json config")
print("\nFor this example, we'll create a mock config:")

# Create temporary directory with mock pseudopotential files
with tempfile.TemporaryDirectory() as tmpdir:
    # Simulate downloading pseudopotentials by creating mock files
    import os
    for element, filename in [("Fe", "Fe.pbe-spn.UPF"), 
                              ("O", "O.pbe-n.UPF"),
                              ("H", "H.pbe-rrkjus.UPF")]:
        filepath = os.path.join(tmpdir, filename)
        with open(filepath, 'w') as f:
            f.write(f"# Mock {element} pseudopotential\n")
    
    print(f"\n✓ Created mock pseudopotentials in: {tmpdir}")
    
    # Create the config by scanning the directory
    config = create_pseudopotentials_config(
        name="example_SSSP_PBE",
        base_path=tmpdir,
        description="Example SSSP PBE pseudopotentials for Fe, O, H",
        functional="PBE",
        library="SSSP",
        version="1.1.2",
        save=True,
        overwrite=True
    )
    
    print(f"✓ Created pseudopotentials config: example_SSSP_PBE")
    print(f"  - Location: ~/.xespresso/pseudopotentials/example_SSSP_PBE.json")
    print(f"  - Available elements: {config.list_elements()}")
    
    # Verify the config was saved
    loaded_config = load_pseudopotentials_config("example_SSSP_PBE")
    print(f"✓ Verified config can be loaded: {loaded_config.name}")
    
    # ==========================================================================
    # 2. EXAMPLE 1: Simple Fe structure
    # ==========================================================================
    print("\n2️⃣  EXAMPLE 1: Simple Fe Structure (Single Element)")
    print("-" * 80)
    
    fe_atoms = bulk("Fe", cubic=True)
    print(f"\nStructure: {len(fe_atoms)} Fe atoms")
    print(f"Elements: {set(fe_atoms.get_chemical_symbols())}")
    
    print("\n✨ NEW WAY - Pass just the config name:")
    print("""
    workflow = CalculationWorkflow(
        atoms=fe_atoms,
        protocol='moderate',
        pseudopotentials_config='example_SSSP_PBE'
    )
    """)
    
    workflow1 = CalculationWorkflow(
        atoms=fe_atoms,
        protocol='moderate',
        pseudopotentials_config='example_SSSP_PBE'
    )
    
    print(f"✓ Workflow created successfully!")
    print(f"✓ Extracted pseudopotentials: {workflow1.original_pseudopotentials}")
    
    # ==========================================================================
    # 3. EXAMPLE 2: Multi-element structure (FeO)
    # ==========================================================================
    print("\n3️⃣  EXAMPLE 2: Multi-element Structure (FeO - Iron Oxide)")
    print("-" * 80)
    
    # Create FeO structure
    from ase import Atoms
    feo_atoms = Atoms(
        symbols=['Fe', 'O'],
        positions=[[0, 0, 0], [2.0, 0, 0]],
        cell=[[4.0, 0, 0], [0, 4.0, 0], [0, 0, 4.0]],
        pbc=True
    )
    
    print(f"\nStructure: FeO compound")
    print(f"Elements: {set(feo_atoms.get_chemical_symbols())}")
    
    print("\n✨ Workflow with multi-element structure:")
    print("""
    workflow = CalculationWorkflow(
        atoms=feo_atoms,
        protocol='fast',
        pseudopotentials_config='example_SSSP_PBE'
    )
    """)
    
    workflow2 = CalculationWorkflow(
        atoms=feo_atoms,
        protocol='fast',
        pseudopotentials_config='example_SSSP_PBE'
    )
    
    print(f"✓ Workflow created successfully!")
    print(f"✓ Extracted pseudopotentials:")
    for element, pseudo_file in workflow2.original_pseudopotentials.items():
        print(f"  - {element}: {pseudo_file}")
    
    # ==========================================================================
    # 4. EXAMPLE 3: Test with from_cif() classmethod
    # ==========================================================================
    print("\n4️⃣  EXAMPLE 3: Using from_cif() with pseudopotentials_config")
    print("-" * 80)
    
    # Save Fe structure as CIF
    cif_file = os.path.join(tmpdir, "fe_structure.cif")
    write(cif_file, fe_atoms)
    
    print(f"\n✓ Created CIF file: {cif_file}")
    
    print("\n✨ Load from CIF and extract pseudos automatically:")
    print("""
    workflow = CalculationWorkflow.from_cif(
        'fe_structure.cif',
        protocol='moderate',
        pseudopotentials_config='example_SSSP_PBE'
    )
    """)
    
    workflow3 = CalculationWorkflow.from_cif(
        cif_file,
        protocol='moderate',
        pseudopotentials_config='example_SSSP_PBE'
    )
    
    print(f"✓ Workflow created from CIF!")
    print(f"✓ Structure read: {len(workflow3.atoms)} Fe atoms")
    print(f"✓ Pseudopotentials extracted: {workflow3.original_pseudopotentials}")
    
    # ==========================================================================
    # 5. COMPARISON: Old vs New
    # ==========================================================================
    print("\n5️⃣  COMPARISON: Old Way vs New Way")
    print("-" * 80)
    
    print("\n❌ OLD WAY - Must manually specify all pseudos:")
    print("""
    workflow = CalculationWorkflow(
        atoms=atoms,
        protocol='moderate',
        pseudopotentials={
            'Fe': 'Fe.pbe-spn.UPF',
            'O': 'O.pbe-n.UPF',
            'H': 'H.pbe-rrkjus.UPF',
        }
    )
    """)
    
    print("\n✅ NEW WAY - Config handles it:")
    print("""
    workflow = CalculationWorkflow(
        atoms=atoms,
        protocol='moderate',
        pseudopotentials_config='example_SSSP_PBE'
    )
    
    # Automatically extracts only needed: {'Fe': 'Fe.pbe-spn.UPF', 'O': 'O.pbe-n.UPF'}
    """)
    
    # ==========================================================================
    # 6. INTEGRATION: Machine + Pseudopotentials Config
    # ==========================================================================
    print("\n6️⃣  INTEGRATION: Machine + Pseudopotentials Config (Separate)")
    print("-" * 80)
    
    print("\n📋 User Configuration Files:")
    print("""
    ~/.xespresso/
    ├── machines/
    │   ├── cluster1.json          ← Machine config (scheduler, resources, etc)
    │   └── local_desktop.json
    ├── pseudopotentials/
    │   ├── SSSP_efficiency.json   ← Pseudo config (elements, filenames)
    │   ├── SSSP_precise.json
    │   └── pbe_standard.json
    """)
    
    print("\n📝 User Code:")
    print("""
    workflow = CalculationWorkflow(
        atoms=my_structure,
        protocol='moderate',
        machine='cluster1',                    # ← Independent: cluster config
        pseudopotentials_config='SSSP_efficiency'  # ← Independent: pseudopotential config
    )
    """)
    
    print("\n✨ How it works:")
    print("  1. machine='cluster1' loads scheduler, resources, host, etc")
    print("  2. pseudopotentials_config='SSSP_efficiency' loads pseudo files")
    print("  3. Workflow extracts only elements needed from your structure")
    print("  4. Everything is automatic - user just provides names!")
    
    # ==========================================================================
    # 7. ERROR HANDLING
    # ==========================================================================
    print("\n7️⃣  ERROR HANDLING")
    print("-" * 80)
    
    # Test missing pseudopotentials config
    print("\n❌ User provides wrong config name:")
    try:
        workflow = CalculationWorkflow(
            atoms=bulk("Fe", cubic=True),
            pseudopotentials_config='nonexistent_config'
        )
    except ValueError as e:
        print(f"   Caught error: {str(e)[:70]}...")
    
    # Test missing element in config
    print("\n❌ Config missing an element from structure:")
    from ase.build import molecule
    h2_atoms = molecule('H2')
    h2_atoms.center(vacuum=5.0)
    try:
        workflow = CalculationWorkflow(
            atoms=h2_atoms,
            pseudopotentials_config='example_SSSP_PBE'
        )
    except ValueError as e:
        print(f"   Caught error: {str(e)[:70]}...")

print("\n" + "=" * 80)
print("✅ EXAMPLE COMPLETE")
print("=" * 80)
print("""
Summary of improvements:
1. ✓ Cleaner API: Just pass config name
2. ✓ Automatic extraction: Only needed elements imported
3. ✓ Better separation: machine config ≠ pseudo config
4. ✓ Backward compatible: Old dict style still works
5. ✓ Better errors: Clear messages when config missing

User workflow:
  1. Create pseudopotential config once: create_pseudopotentials_config()
  2. Create machine config once: create_machine()
  3. Then use both with simple names in calculations!
""")
print("=" * 80)
