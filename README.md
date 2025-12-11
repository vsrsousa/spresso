## spresso
Qt GUI for Quantum Espresso

*This project is a fork of [xespresso](https://github.com/vsrsousa/xespresso) (qe_gui branch)*

Quantum ESPRESSO Calculator for Atomic Simulation Environment (ASE) with an integrated Qt GUI.

**Note on naming:** The repository and PyPI package are named `spresso`, but the Python module is `xespresso` (for backwards compatibility). Install with `pip install spresso`, but import as `from xespresso import ...`

For the introduction of ASE , please visit https://wiki.fysik.dtu.dk/ase/index.html


### Features

* Support all QE packages, including: pw, band, neb, dos, projwfc, pp ...
* Spin-polarized calculation
* LD(S)A+U
* Automatic submit job
* Automatic check a new calculation required or not
* Automatic set up "nscf" calculation
* Read and plot dos, pdos and layer resolved pdos
* Plot NEB
* **NEW: Simplified workflow with quality presets**
* **NEW: K-spacing support for easy k-point specification**
* **NEW: Pseudopotential configuration management**
* **NEW: Qt GUI for easy configuration and job submission** 🚀

### Author
* Xing Wang  <xingwang1991@gmail.com>

### Dependencies

* Python
* ASE
* numpy
* scipy
* matplotlib

### Installation using pip
```bash
pip install --upgrade --user spresso
```

### Installation from source
You can get the source using git:
``` sh
git clone https://github.com/vsrsousa/spresso.git
cd spresso
pip install -e .
```

**Optional environment variables:**
``` sh
export ASE_ESPRESSO_COMMAND="/path/to/PACKAGE.x  PARALLEL  -in  PREFIX.PACKAGEi  >  PREFIX.PACKAGEo"
export ESPRESSO_PSEUDO="/path/to/pseudo"
```

**Note:** When installed via pip (either from PyPI or from source), you don't need to manually add anything to PYTHONPATH. Simply import with `from xespresso import ...`

### GUI Interface (NEW!) 🎨

spresso includes a user-friendly Qt-based GUI for easy configuration and job submission!

**Features:**
- 🖥️ Configure machines (local/remote) with scheduler support
- ⚙️ Auto-detect Quantum ESPRESSO codes
- 🔬 Upload and visualize structures in 3D
- 📊 Configure calculations with quality presets
- 🔄 Build workflows with easy parameter selection
- 🚀 Submit jobs to configured machines


### Examples

#### Simplified Workflow (NEW!)

**🚀 NEW: Easy calculations with quality presets and k-spacing**

Run calculations from CIF files with quality presets (`fast`, `moderate`, `accurate`):

``` python
from xespresso import quick_scf, quick_relax

# Quick SCF calculation from CIF file
calc = quick_scf(
    'structure.cif',
    {'Fe': 'Fe.pbe-spn.UPF'},
    quality='moderate',
    label='scf/fe'
)

# Quick relaxation with k-spacing (instead of k-points)
calc = quick_relax(
    'structure.cif',
    {'Fe': 'Fe.pbe-spn.UPF'},
    quality='moderate',
    kspacing=0.3,  # Angstrom^-1, converted automatically
    relax_type='vc-relax',
    label='relax/fe'
)
```

**Quality Presets:**
- `fast`: Quick calculations for testing (ecutwfc=30 Ry, kspacing=0.5)
- `moderate`: Standard production runs (ecutwfc=50 Ry, kspacing=0.3)
- `accurate`: High-precision results (ecutwfc=80 Ry, kspacing=0.15)

**Pseudopotential Configuration Management:**

Store and reuse pseudopotential configurations:

``` python
from xespresso.utils import save_pseudo_config, load_pseudo_config

# Save configuration to ~/.xespresso/
config = {
    "name": "my_config",
    "pseudopotentials": {
        "Fe": "Fe.pbe-spn.UPF",
        "O": "O.pbe.UPF"
    }
}
save_pseudo_config("my_config", config)

# Load and use
config = load_pseudo_config("my_config")
calc = quick_scf('structure.cif', config['pseudopotentials'], quality='moderate')
```

See [WORKFLOW_DOCUMENTATION.md](WORKFLOW_DOCUMENTATION.md) for complete documentation.

#### Automatic submit job

A example of setting parameters for the queue. See example/queue.py

``` python
queue = {'nodes': 4,
         'ntasks-per-node': 20,
         'partition': 'all',
         'time': '23:10:00'}
calc = Espresso(queue = queue)
```




#### Automatic check a new calculation required or not

Before the calculation, it will first check the working directory. If the same geometry and parameters are used, try to check whether the results are available or not. Automatic check input parameters with Quantum Espresso document.

``` python
calc = Espresso(label = 'scf/fe')
```

#### Show debug information.

``` python
calc = Espresso(debug = True)
```

#### Magnetic configuration (New Simplified API!)

**🎉 NEW: Element-based magnetic configuration** - The easiest way to set up magnetic systems!

``` python
from xespresso import setup_magnetic_config

atoms = bulk('Fe', cubic=True)

# All Fe equivalent with magnetization 1
config = setup_magnetic_config(atoms, {'Fe': [1]})

# AFM: two non-equivalent Fe
config = setup_magnetic_config(atoms, {'Fe': [1, -1]})

# With Hubbard U
config = setup_magnetic_config(atoms, {
    'Fe': {'mag': [1, -1], 'U': 4.3}
})

# Multiple elements (e.g., FeMnAl2)
config = setup_magnetic_config(atoms, {
    'Fe': [1],        # All Fe equivalent
    'Mn': [1, -1],    # Mn antiferromagnetic
    'Al': [0]         # Al non-magnetic
})

# Auto-expand cell if needed
config = setup_magnetic_config(
    atoms, 
    {'Fe': [1, 1, -1, -1]},  # Need 4 Fe but only have 2
    expand_cell=True
)

# Use in calculator - NOW EVEN SIMPLER! 🎉
# Method 1: Pass entire config (recommended)
config = setup_magnetic_config(atoms, {'Fe': [1, -1]}, 
                               pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'})
calc = Espresso(
    atoms=config['atoms'],
    input_data=config,  # Pass entire dict - no manual extraction!
    nspin=2,
    ecutwfc=40
)

# Method 2: Old way still works
config['pseudopotentials']['Fe'] = 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
config['pseudopotentials']['Fe1'] = 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
calc = Espresso(
    pseudopotentials=config['pseudopotentials'],
    input_data={'input_ntyp': config['input_ntyp']},
    nspin=2
)
```

**✨ New Feature**: Pseudopotentials can now be passed inside `input_data`! The calculator automatically extracts them, making it super convenient to use `setup_magnetic_config`. See `PSEUDOPOTENTIAL_AUTO_EXTRACTION.md` for details.

**Also available:** `set_antiferromagnetic()`, `set_ferromagnetic()`, `set_magnetic_moments()`

See `MAGNETIC_HELPERS.md` for complete documentation and examples.

#### Add new species (Manual method)
Some atoms are special:
+ atoms with different starting_magnetization
+ atoms with different U values
+ atoms with special basis set

For example, Fe with spin state AFM. See example/spin.py

``` python
atoms.new_array('species', np.array(atoms.get_chemical_symbols(), dtype = 'U20'))
atoms.arrays['species'][1] = 'Fe1'
```

#### Setting parameters with "(i), i=1,ntyp"
Hubbard, starting_magnetization, starting_charge and so on. See example/dft+u.py

``` python
input_ntyp = {
'starting_magnetization': {'Fe1': 1.0, 'Fe2': -1.0},
'Hubbard_U': {'Fe1': 3.0, 'Fe2': 3.0},
}
input_data['input_ntyp'] = input_ntyp,
```

#### Setting parameters for "Hubbard_V(na,nb,k)"
Hubbard, starting_magnetization, starting_charge and so on. See example/dft+u.py

``` python
input_data = {
'hubbard_v': {'(1,1,1)': 4.0, '(3,3,1)': 1.0},
}
```

#### Control parallelization levels
To control the number of processors in each group: -ni,
-nk, -nb, -nt, -nd) are used.

``` python
calc = Espresso(pseudopotentials = pseudopotentials,
                 package = 'pw',
                 parallel = '-nk 2 -nt 4 -nd 144',  # parallel parameters
                 }
```

#### Non self-consistent calculation

A example of nscf calculation following the above one.

``` python
# start nscf calculation
from xespresso.post.nscf import EspressoNscf
nscf = EspressoNscf(calc.directory, prefix = calc.prefix,
                occupations = 'tetrahedra',
                kpts = (2, 2, 2),
                debug = True,
                )
nscf.run()
```

#### Calculate dos and pdos

A example of calculating and plotting the pdos from the nscf calculation.

``` python
from xespresso.post.dos import EspressoDos
# dos
dos = EspressoDos(parent_directory = 'calculations/scf/co',
            prefix = calc.prefix,
            Emin = fe - 30, Emax = fe + 30, DeltaE = 0.01)
dos.run()
# pdos
from xespresso.post.projwfc import EspressoProjwfc
projwfc = EspressoProjwfc(parent_directory = 'calculations/scf/co',
            prefix = 'co',
            DeltaE = 0.01)
projwfc.run()
```
<img src="docs/source/_static/images/co-pdos.png" width="500"/>

#### Calculate work function
``` python
from xespresso.post.pp import EspressoPp
pp = EspressoPp(calc.directory, prefix = calc.prefix,
                plot_num = 11,
                fileout = 'potential.cube',
                iflag = 3,
                output_format=6,
                debug = True,
                )
pp.get_work_function()
```

#### Restart from previous calculation
``` python
calc.read_results()
atoms = calc.results['atoms']
calc.run(atoms = atoms, restart = 1)
```

#### NEB calculation
See example/neb.py
``` python
from xespresso.neb import NEBEspresso
calc = NEBEspresso(
                 package = 'neb',
                 images = images,
                 climbing_images = [5],
                 path_data = path_data
                 )
calc.calculate()
calc.read_results()
calc.plot()
```
<img src="docs/source/_static/images/neb.png" width="500"/>


## Workflow
### Oxygen evolution reaction (OER) calculation

The workflow includes four modules:
* OER_bulk
* OER_pourbaix
* OER_surface
* OER_site


The workflow can handle:
* Generate surface slab model from bulk structure
* Determine the surface adsorption site
* Determine the surface coverage(*, O*, OH*), Pourbaix diagram
* Calculate the Zero-point energy


```python
oer = OER_site(slab,
               label = 'oer/Pt-001-ontop',
               site_type = 'ontop',
               site = -1,
               height=2.0,
	           calculator = parameters,
               molecule_energies = molecule_energies,
               )
oer.run()
```

### To do lists:
* add `qPointsSpecs` and `Line-of-input` for phonon input file
