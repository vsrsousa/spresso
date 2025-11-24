"""
Example demonstrating the modular GUI design pattern.

This example shows how the GUI calculation and workflow modules work together
following xespresso's design patterns:
- Calculation modules create Espresso and atoms objects
- Job submission receives prepared objects and executes them
"""

from ase.build import bulk
from ase import io

# Example 1: Direct use of calculation module (simulating GUI workflow)
print("="*70)
print("Example 1: Using Calculation Module")
print("="*70)

# In GUI, atoms would come from structure viewer
atoms = bulk('Fe', cubic=True)
print(f"Structure: {atoms.get_chemical_formula()} ({len(atoms)} atoms)")

# In GUI, this config would come from calculation_setup page
config = {
    'pseudopotentials': {'Fe': 'Fe.pbe-spn-kjpaw_psl.0.2.1.UPF'},
    'calc_type': 'scf',
    'ecutwfc': 50.0,
    'ecutrho': 400.0,
    'kpts': (4, 4, 4),
    'occupations': 'smearing',
    'degauss': 0.02,
    'conv_thr': 1.0e-8,
}

# Use calculation module to prepare atoms and Espresso calculator
# This is what calculation_setup.py does
print("\n1. Preparing calculation using calculation module...")
try:
    from gui.calculations import prepare_calculation_from_gui
    
    label = "/tmp/example_scf/fe"
    prepared_atoms, calc = prepare_calculation_from_gui(atoms, config, label=label)
    
    print(f"✓ Calculation module created:")
    print(f"  - Prepared atoms: {prepared_atoms.get_chemical_formula()}")
    print(f"  - Espresso calculator with label: {calc.label}")
    print(f"  - Calculation type: {config['calc_type']}")
    
    # In GUI, these would be stored in session_state
    # st.session_state.espresso_calculator = calc
    # st.session_state.prepared_atoms = prepared_atoms
    
except ImportError as e:
    print(f"⚠ Calculation module not available: {e}")
    print("This example requires the GUI calculation modules to be in PYTHONPATH")

# Example 2: Dry run (generate input files without execution)
print("\n" + "="*70)
print("Example 2: Dry Run (Generate Input Files)")
print("="*70)

try:
    from gui.calculations import dry_run_calculation
    
    label = "/tmp/example_dryrun/fe"
    print("\n2. Performing dry run...")
    prepared_atoms, calc = dry_run_calculation(atoms, config, label=label)
    
    print(f"✓ Dry run complete!")
    print(f"  - Input files written to: {calc.directory}")
    print(f"  - Calculator prepared but not executed")
    
except ImportError as e:
    print(f"⚠ Calculation module not available: {e}")

# Example 3: Workflow with multiple steps
print("\n" + "="*70)
print("Example 3: Multi-Step Workflow")
print("="*70)

try:
    from gui.workflows import GUIWorkflow
    
    print("\n3. Creating workflow with SCF + Relaxation...")
    
    # Create workflow
    base_label = "/tmp/example_workflow"
    workflow = GUIWorkflow(atoms, config, base_label=base_label)
    
    # Add SCF step
    scf_config = config.copy()
    scf_config['calc_type'] = 'scf'
    scf_atoms, scf_calc = workflow.add_calculation('scf', scf_config)
    print(f"✓ Added SCF step")
    
    # Add relaxation step
    relax_config = config.copy()
    relax_config['calc_type'] = 'relax'
    relax_config['forc_conv_thr'] = 1.0e-3
    relax_atoms, relax_calc = workflow.add_calculation('relax', relax_config)
    print(f"✓ Added relaxation step")
    
    print(f"\n✓ Workflow created with {len(workflow.calculations)} steps:")
    for name in workflow.calculations.keys():
        print(f"  - {name}")
    
    # In GUI job_submission, you would execute the workflow:
    # workflow.run_step('scf')
    # workflow.run_step('relax')
    
except ImportError as e:
    print(f"⚠ Workflow module not available: {e}")

print("\n" + "="*70)
print("Summary: Modular Design Pattern")
print("="*70)
print("""
Key principles demonstrated:

1. Calculation modules (gui/calculations/) create Espresso and atoms objects
   - prepare_calculation_from_gui(): Creates calculator from GUI config
   - dry_run_calculation(): Prepares and writes input files

2. Workflow modules (gui/workflows/) orchestrate multiple calculations
   - GUIWorkflow: Coordinates multiple calculation steps
   - Each step uses calculation modules to prepare objects

3. Job submission receives prepared objects
   - Gets atoms and calculator from calculation/workflow modules
   - Only responsible for execution, not object creation

4. Uses xespresso's logic and definitions
   - All Espresso creation uses xespresso's Espresso class
   - Follows xespresso's parameter patterns (input_data, kpts, etc.)

This pattern ensures:
✓ Clear separation of concerns
✓ Reusable calculation preparation logic
✓ Consistent use of xespresso patterns
✓ Easy testing and maintenance
""")
