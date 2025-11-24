"""
Simplified example for antiferromagnetic spin calculation using new helper functions.

This example shows how the new helper functions simplify the setup of
antiferromagnetic calculations compared to the manual method.
"""

from ase.build import bulk
from xespresso import Espresso
from xespresso.tools import set_antiferromagnetic

# Create Fe structure
atoms = bulk("Fe", cubic=True)

# OLD METHOD (manual - for comparison):
# atoms.new_array("species", np.array(atoms.get_chemical_symbols(), dtype="U20"))
# atoms.arrays["species"][0] = "Fe"
# atoms.arrays["species"][1] = "Fe1"
# input_ntyp = {
#     "starting_magnetization": {
#         "Fe": 1.0,
#         "Fe1": -1.0,
#     }
# }
# pseudopotentials = {
#     "Fe": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF",
#     "Fe1": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF",
# }

# NEW METHOD (simplified):
# Automatically set up antiferromagnetic configuration
mag_config = set_antiferromagnetic(
    atoms, 
    sublattice_indices=[[0], [1]],  # Atom 0 has +moment, atom 1 has -moment
    magnetic_moment=1.0
)

# Update with actual pseudopotential file names
mag_config['pseudopotentials']['Fe'] = "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF"
mag_config['pseudopotentials']['Fe1'] = "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF"

# Create calculator with simplified configuration
calc = Espresso(
    pseudopotentials=mag_config['pseudopotentials'],
    label="scf/fe-afm",
    ecutwfc=40,
    occupations="smearing",
    degauss=0.02,
    nspin=2,
    input_data={'input_ntyp': mag_config['input_ntyp']},
    kpts=(4, 4, 4),
    debug=True,
)

atoms.calc = calc
e = atoms.get_potential_energy()
print("Energy: {0:1.3f}".format(e))
