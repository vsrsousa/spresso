"""
Example: Using new Hubbard parameters format (QE >= 7.0)

This example demonstrates how to use the new HUBBARD card format
introduced in Quantum ESPRESSO 7.0, which provides more flexibility
for DFT+U and DFT+U+V calculations.
"""

from ase.build import bulk
import numpy as np

# Note: This example shows the input format, but won't run actual calculations
# without a proper QE installation

print("="*70)
print("Example 1: New HUBBARD card format (QE >= 7.0)")
print("="*70)

# Define Hubbard parameters using the new format
# This is more flexible and explicit about orbitals

input_data_new = {
    "ecutwfc": 30.0,
    "occupations": "smearing",
    "degauss": 0.02,
    "nspin": 2,
    "lda_plus_u": True,
    "qe_version": "7.2",  # Specify QE version to use new format
    "hubbard": {
        "projector": "atomic",  # or 'ortho-atomic', 'norm-atomic', 'wf', 'pseudo'
        "u": {
            "Fe-3d": 4.3,  # U value for Fe 3d orbital
            "O-2p": 3.0,   # U value for O 2p orbital
        },
        "v": [
            {
                "species1": "Fe",
                "orbital1": "3d",
                "species2": "O",
                "orbital2": "2p",
                "i": 1,
                "j": 1,
                "value": 0.5
            }
        ]
    }
}

print("\nNew format configuration:")
print(f"  Projector: {input_data_new['hubbard']['projector']}")
print(f"  U parameters:")
for species_orb, value in input_data_new['hubbard']['u'].items():
    print(f"    {species_orb}: {value} eV")
print(f"  V parameters: {len(input_data_new['hubbard']['v'])} inter-site interactions")

print("\nThis will generate a HUBBARD card in the input file:")
print("""
HUBBARD {atomic}
  U Fe-3d 4.3
  U O-2p 3.0
  V Fe-3d O-2p 1 1 0.5
""")

print("\n" + "="*70)
print("Example 2: Old format (QE < 7.0) - backward compatible")
print("="*70)

# Old format using INPUT_NTYP
input_ntyp = {
    "starting_magnetization": {
        "Fe": 0.5,
    },
    "Hubbard_U": {
        "Fe": 4.3,
        "O": 3.0,
    },
}

input_data_old = {
    "ecutwfc": 30.0,
    "occupations": "smearing",
    "degauss": 0.02,
    "nspin": 2,
    "lda_plus_u": True,
    "input_ntyp": input_ntyp,
    # For inter-site V parameters in old format:
    "hubbard_v": {"(1,2,1)": 0.5},
}

print("\nOld format configuration:")
print(f"  U parameters:")
for species, value in input_ntyp['Hubbard_U'].items():
    print(f"    {species}: {value} eV")
print(f"  V parameters: hubbard_v dictionary")

print("\nThis generates parameters in SYSTEM namelist:")
print("""
&SYSTEM
  ...
  lda_plus_u = .true.
  Hubbard_U(1) = 4.3
  Hubbard_U(2) = 3.0
  Hubbard_V(1,2,1) = 0.5
  ...
/
""")

print("\n" + "="*70)
print("Example 3: Automatic format detection")
print("="*70)

print("""
The xespresso library will automatically detect which format to use based on:
1. Explicit 'qe_version' parameter (if >= 7.0, use new format)
2. Presence of 'hubbard' dictionary (indicates new format)
3. Presence of orbital specifications like 'Fe-3d' (indicates new format)
4. Default: use old format for backward compatibility

You can force a specific format:
  input_data['hubbard_format'] = 'card'  # Force new format
  input_data['hubbard_format'] = 'namelist'  # Force old format
""")

print("\n" + "="*70)
print("Example 4: Complete calculation setup (new format)")
print("="*70)

atoms = bulk("Fe", cubic=True)
atoms.new_array("species", np.array(atoms.get_chemical_symbols(), dtype="U20"))

input_data_complete = {
    "ecutwfc": 30.0,
    "ecutrho": 240.0,
    "occupations": "smearing",
    "smearing": "gaussian",
    "degauss": 0.02,
    "nspin": 2,
    "lda_plus_u": True,
    "qe_version": "7.2",
    "hubbard": {
        "projector": "atomic",
        "u": {
            "Fe-3d": 4.3,
        }
    }
}

pseudopotentials = {
    "Fe": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF",
}

print("\nTo use this with xespresso.Espresso:")
print("""
from xespresso import Espresso

calc = Espresso(
    pseudopotentials=pseudopotentials,
    label="scf/fe",
    input_data=input_data_complete,
    kpts=(4, 4, 4),
)
atoms.calc = calc
energy = atoms.get_potential_energy()
""")

print("\n" + "="*70)
print("Benefits of new format:")
print("="*70)
print("""
1. ✅ Explicit orbital specification (more precise)
2. ✅ Easier to specify complex DFT+U+V interactions
3. ✅ Better support for multiple Hubbard sites
4. ✅ More flexible projector choices
5. ✅ Clearer and more readable input files
6. ✅ Recommended for QE 7.0 and later
7. ✅ Backward compatible (old format still works)
""")

print("\n✅ Examples completed!")
