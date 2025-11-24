"""
Examples demonstrating the new element-based magnetic configuration API.

This shows how to use setup_magnetic_config() for intuitive magnetic setup.
"""

from ase.build import bulk
from ase import Atoms
from xespresso import setup_magnetic_config, Espresso

print("="*70)
print("ELEMENT-BASED MAGNETIC CONFIGURATION EXAMPLES")
print("="*70)

# Example 1: Simple AFM with equivalent Fe atoms
print("\n1. Fe=[1] - All Fe equivalent with magnetization 1")
print("-"*70)
atoms = bulk('Fe', cubic=True)  # 2 Fe atoms
config = setup_magnetic_config(atoms, {'Fe': [1]})

print(f"Number of species: {len(config['input_ntyp']['starting_magnetization'])}")
print(f"Magnetizations: {config['input_ntyp']['starting_magnetization']}")
print("→ Both Fe atoms are equivalent, same species")

# Example 2: AFM with non-equivalent Fe atoms
print("\n2. Fe=[1, -1] - Two non-equivalent Fe atoms")
print("-"*70)
atoms = bulk('Fe', cubic=True)  # 2 Fe atoms
config = setup_magnetic_config(atoms, {'Fe': [1, -1]})

print(f"Number of species: {len(config['input_ntyp']['starting_magnetization'])}")
print(f"Magnetizations: {config['input_ntyp']['starting_magnetization']}")
print("→ Two different species: Fe (mag=1) and Fe1 (mag=-1)")

# Example 3: Complex FeMnAl2 system
print("\n3. FeMnAl2 - Multiple elements with different configs")
print("-"*70)
atoms = Atoms('Fe2Mn2Al4', positions=[
    [0, 0, 0], [1, 0, 0],           # Fe
    [0, 1, 0], [1, 1, 0],           # Mn
    [0, 0, 1], [1, 0, 1],           # Al
    [0, 1, 1], [1, 1, 1]            # Al
])
atoms.cell = [5, 5, 5]

config = setup_magnetic_config(atoms, {
    'Fe': [1],        # Both Fe equivalent, mag=1
    'Mn': [1, -1],    # Mn AFM, non-equivalent
    'Al': [0]         # Al non-magnetic
})

print(f"Magnetic species: {list(config['input_ntyp']['starting_magnetization'].keys())}")
print(f"Magnetizations: {config['input_ntyp']['starting_magnetization']}")
print("→ Fe: both equivalent, Mn: AFM, Al: no magnetization")

# Example 4: Supercell expansion
print("\n4. Fe=[1, 1, -1, -1] with only 2 Fe - auto-expand")
print("-"*70)
atoms = bulk('Fe', cubic=True)  # 2 Fe atoms
config = setup_magnetic_config(
    atoms, 
    {'Fe': [1, 1, -1, -1]},
    expand_cell=True
)

print(f"Cell expanded: {config['expanded']}")
print(f"Number of atoms: {len(config['atoms'])}")
print(f"Number of species: {len(config['input_ntyp']['starting_magnetization'])}")
print(f"Magnetizations: {config['input_ntyp']['starting_magnetization']}")
print("→ Cell expanded from 2 to 4 atoms to accommodate configuration")

# Example 5: With Hubbard U - same for all
print("\n5. Fe=[1, -1] with Hubbard U=4.3")
print("-"*70)
atoms = bulk('Fe', cubic=True)
config = setup_magnetic_config(atoms, {
    'Fe': {'mag': [1, -1], 'U': 4.3}
})

print(f"Magnetizations: {config['input_ntyp']['starting_magnetization']}")
print(f"Hubbard U: {config['input_ntyp']['Hubbard_U']}")
print("→ Both species get U=4.3")

# Example 6: Different Hubbard U for each species
print("\n6. Fe=[1, -1] with different U values")
print("-"*70)
atoms = bulk('Fe', cubic=True)
config = setup_magnetic_config(atoms, {
    'Fe': {'mag': [1, -1], 'U': [4.3, 4.5]}
})

print(f"Magnetizations: {config['input_ntyp']['starting_magnetization']}")
print(f"Hubbard U: {config['input_ntyp']['Hubbard_U']}")
print("→ Fe gets U=4.3, Fe1 gets U=4.5")

# Example 7: Pattern replication
print("\n7. Fe=[1, -1] with 4 Fe atoms - pattern replicates")
print("-"*70)
atoms = bulk('Fe', cubic=True) * (2, 1, 1)  # 4 Fe atoms
config = setup_magnetic_config(atoms, {'Fe': [1, -1]})

species_list = list(config['atoms'].arrays['species'])
print(f"Species pattern: {species_list}")
print(f"Magnetizations: {config['input_ntyp']['starting_magnetization']}")
print("→ Pattern [1, -1] replicates to [1, -1, 1, -1]")

# Example 8: Ready for Espresso calculator
print("\n8. Complete Espresso setup")
print("-"*70)
atoms = bulk('Fe', cubic=True)
config = setup_magnetic_config(atoms, {
    'Fe': {'mag': [1, -1], 'U': 4.3}
})

# Add pseudopotential files
config['pseudopotentials']['Fe'] = 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
config['pseudopotentials']['Fe1'] = 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'

print("Configuration ready:")
print(f"  Pseudopotentials: {config['pseudopotentials']}")
print(f"  input_ntyp: {config['input_ntyp']}")
print("\nEspresso calculator setup:")
print("""
calc = Espresso(
    pseudopotentials=config['pseudopotentials'],
    input_data={'input_ntyp': config['input_ntyp']},
    nspin=2,
    lda_plus_u=True,
    ecutwfc=40,
    kpts=(4, 4, 4)
)
config['atoms'].calc = calc
""")

print("\n" + "="*70)
print("All examples completed successfully!")
print("="*70)
