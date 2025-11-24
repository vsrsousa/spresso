"""
Simplified example for MnO antiferromagnetic DFT+U calculation.

This example shows how to use the new helper functions to simplify
setting up complex antiferromagnetic configurations with DFT+U.
"""

from ase.io import read
from xespresso import Espresso
from xespresso.tools import set_magnetic_moments
import numpy as np

# Load structure (assuming datas/MnO.cif exists)
# atoms = read("datas/MnO.cif")

# For demonstration, create a simple structure
from ase.build import bulk
atoms = bulk("Mn", cubic=True) * (2, 2, 1)

# OLD METHOD (manual - for comparison):
# atoms.new_array("species", np.array(atoms.get_chemical_symbols(), dtype="U20"))
# for i in range(4):
#     atoms.arrays["species"][2 * i] = atoms.arrays["species"][2 * i] + "%s" % (i + 1)
#     atoms.arrays["species"][2 * i + 1] = atoms.arrays["species"][2 * i + 1] + "%s" % (i + 1)
# input_ntyp = {
#     "starting_magnetization": {
#         "Mn1": 1.0,
#         "Mn2": -1.0,
#         "Mn3": -1.0,
#         "Mn4": 1.0,
#     },
#     "Hubbard_U": {
#         "Mn1": 5.75,
#         "Mn2": 5.75,
#         "Mn3": 5.75,
#         "Mn4": 5.75,
#     },
# }

# NEW METHOD (simplified):
# Define magnetic moments for each atom
# Checkerboard AFM pattern: +1, -1, -1, +1
magnetic_moments = {
    0: 1.0,   # Mn1
    1: -1.0,  # Mn2
    2: -1.0,  # Mn3
    3: 1.0,   # Mn4
}

# Base pseudopotentials
base_pseudopotentials = {
    "Mn": "Mn.pbesol-spn-rrkjus_psl.0.3.1.UPF",
}

# Set up magnetic configuration
mag_config = set_magnetic_moments(
    atoms, 
    magnetic_moments,
    pseudopotentials=base_pseudopotentials
)

# Add Hubbard U values to input_ntyp
# All Mn species get the same U value
mag_config['input_ntyp']['Hubbard_U'] = {}
for species in mag_config['pseudopotentials'].keys():
    if 'Mn' in species:
        mag_config['input_ntyp']['Hubbard_U'][species] = 5.75

# Create input_data
input_data = {
    "ecutwfc": 70.0,
    "ecutrho": 840.0,
    "occupations": "smearing",
    "degauss": 0.01,
    "nspin": 2,
    "lda_plus_u": True,
    "input_ntyp": mag_config['input_ntyp'],
}

# Setup calculator
calc = Espresso(
    pseudopotentials=mag_config['pseudopotentials'],
    label="scf/mno-afm",
    input_data=input_data,
    kpts=(4, 4, 4),
    debug=True,
)

atoms.set_calculator(calc)
e = atoms.get_potential_energy()

print("Energy: {0:1.3f}".format(e))
print("\nMagnetic configuration:")
print("Species:", list(mag_config['input_ntyp']['starting_magnetization'].keys()))
print("Moments:", list(mag_config['input_ntyp']['starting_magnetization'].values()))
