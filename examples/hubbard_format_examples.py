"""
Examples demonstrating Hubbard parameter support in setup_magnetic_config().

This shows how to use both old format (QE < 7.0) and new format (QE >= 7.0)
with orbital specification and V parameters.
"""

from ase.build import bulk
from ase import Atoms
from xespresso import setup_magnetic_config

print("="*80)
print("HUBBARD PARAMETER SUPPORT - OLD AND NEW FORMATS")
print("="*80)

# Example 1: Old format (QE < 7.0) - Simple U parameter
print("\n1. OLD FORMAT (QE < 7.0) - Simple Hubbard U")
print("-"*80)
atoms = bulk('Fe', cubic=True)

config = setup_magnetic_config(atoms, {
    'Fe': {'mag': [1, -1], 'U': 4.3}
}, qe_version='6.8')

print(f"Format: {config['hubbard_format']}")
print(f"Hubbard_U in input_ntyp: {config['input_ntyp'].get('Hubbard_U', {})}")
print("→ Old format: U values in SYSTEM namelist as Hubbard_U(i)")

# Example 2: New format (QE 7.x) - U with orbital specification
print("\n2. NEW FORMAT (QE 7.x) - Hubbard U with orbital")
print("-"*80)
atoms = bulk('Fe', cubic=True)

config = setup_magnetic_config(atoms, {
    'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}
}, qe_version='7.2')

print(f"Format: {config['hubbard_format']}")
print(f"Hubbard card data: {config.get('hubbard', {})}")
print("→ New format: HUBBARD card with Fe-3d and Fe1-3d")

# Example 3: Different U values per species (new format)
print("\n3. NEW FORMAT - Different U for each species")
print("-"*80)
atoms = bulk('Fe', cubic=True)

config = setup_magnetic_config(atoms, {
    'Fe': {'mag': [1, -1], 'U': {'3d': [4.3, 4.5]}}
}, qe_version='7.2')

hubbard = config['hubbard']
print(f"U parameters:")
for species_orbital, value in hubbard['u'].items():
    print(f"  {species_orbital}: {value} eV")
print("→ Fe-3d gets 4.3 eV, Fe1-3d gets 4.5 eV")

# Example 4: V parameter (inter-site interaction)
print("\n4. NEW FORMAT - V parameter (inter-site interaction)")
print("-"*80)
atoms = Atoms('FeO', positions=[[0, 0, 0], [1.5, 0, 0]])
atoms.cell = [5, 5, 5]

config = setup_magnetic_config(atoms, {
    'Fe': {
        'mag': [1],
        'U': {'3d': 4.3},
        'V': [{'species2': 'O', 'orbital1': '3d', 'orbital2': '2p', 'value': 1.0}]
    },
    'O': [0]
}, qe_version='7.2')

hubbard = config['hubbard']
print(f"U parameters: {hubbard['u']}")
print(f"V parameters:")
for v in hubbard['v']:
    print(f"  {v['species1']}-{v['orbital1']} <-> {v['species2']}-{v['orbital2']}: {v['value']} eV")
print("→ Inter-site interaction between Fe-3d and O-2p")

# Example 5: Complex system with multiple elements
print("\n5. NEW FORMAT - Complex FeMnO system")
print("-"*80)
atoms = Atoms('Fe2Mn2O4', positions=[
    [0, 0, 0], [1, 0, 0],           # Fe
    [0, 1, 0], [1, 1, 0],           # Mn
    [0, 0, 1], [1, 0, 1],           # O
    [0, 1, 1], [1, 1, 1]            # O
])
atoms.cell = [5, 5, 5]

config = setup_magnetic_config(atoms, {
    'Fe': {
        'mag': [1],
        'U': {'3d': 4.3}
    },
    'Mn': {
        'mag': [1, -1],
        'U': {'3d': [5.7, 5.8]}
    },
    'O': [0]
}, qe_version='7.2')

hubbard = config['hubbard']
print(f"U parameters:")
for species_orbital, value in hubbard['u'].items():
    print(f"  {species_orbital}: {value} eV")
print("→ Different U for Fe-3d, Mn-3d, Mn1-3d")

# Example 6: Auto-detection based on format
print("\n6. AUTO-DETECTION - Format determined by U syntax")
print("-"*80)
atoms = bulk('Fe', cubic=True)

# Dict format → auto-detects new format
config = setup_magnetic_config(atoms, {
    'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}
})

print(f"Format (no qe_version specified): {config['hubbard_format']}")
print("→ Dict U format {'3d': 4.3} auto-detects as new format")

# Example 7: Using hubbard_format parameter
print("\n7. EXPLICIT FORMAT CONTROL - Using hubbard_format parameter")
print("-"*80)
atoms = bulk('Fe', cubic=True)

config = setup_magnetic_config(atoms, {
    'Fe': {'mag': [1, -1], 'U': 4.3}
}, hubbard_format='old')

print(f"Format (hubbard_format='old'): {config['hubbard_format']}")
print(f"Hubbard_U: {config['input_ntyp'].get('Hubbard_U', {})}")
print("→ Explicit format specification overrides auto-detection")

# Example 8: Usage with Espresso calculator
print("\n8. USAGE WITH ESPRESSO CALCULATOR")
print("-"*80)
atoms = bulk('Fe', cubic=True)

config = setup_magnetic_config(atoms, {
    'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}
}, qe_version='7.2')

# Add pseudopotentials
config['pseudopotentials']['Fe'] = 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
config['pseudopotentials']['Fe1'] = 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'

print("Configuration ready for Espresso:")
print(f"  - Format: {config['hubbard_format']}")
print(f"  - QE version: {config.get('qe_version', 'not specified')}")
print(f"  - Pseudopotentials: {list(config['pseudopotentials'].keys())}")
if 'hubbard' in config:
    print(f"  - Hubbard U: {config['hubbard']['u']}")
    print(f"  - Hubbard V: {config['hubbard']['v']}")

print("\nExample Espresso setup:")
print("""
input_data = {
    'ecutwfc': 70.0,
    'nspin': 2,
    'lda_plus_u': True,
    'input_ntyp': config['input_ntyp'],
    'qe_version': config.get('qe_version'),
}

# For new format, also add:
if 'hubbard' in config:
    input_data['hubbard'] = config['hubbard']

calc = Espresso(
    pseudopotentials=config['pseudopotentials'],
    input_data=input_data,
    kpts=(4, 4, 4)
)
config['atoms'].calc = calc
""")

print("\n" + "="*80)
print("All Hubbard format examples completed!")
print("="*80)
