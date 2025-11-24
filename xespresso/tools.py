import os
from ase.geometry import get_layers
from ase.constraints import FixAtoms
from ase.io.espresso import construct_namelist, kspacing_to_grid
from xespresso import Espresso
from xespresso.xio import build_atomic_species_str
from ase.dft.bandgap import bandgap
import pickle
import multiprocessing
import numpy as np

# ====================================================
# K-point utilities
# ====================================================


def kpts_from_spacing(atoms, kspacing):
    """
    Convert k-point spacing to k-point grid without manual 2π normalization.
    
    This is a convenience wrapper around ase.io.espresso.kspacing_to_grid that
    handles the 2π normalization internally, so users can pass k-spacing values
    in physical units (Angstrom^-1) directly.
    
    Parameters
    ----------
    atoms : ase.Atoms
        The atomic structure
    kspacing : float
        K-point spacing in Angstrom^-1 (physical units)
        
    Returns
    -------
    tuple
        K-point grid (kx, ky, kz)
        
    Examples
    --------
    >>> from ase.build import bulk
    >>> from xespresso import kpts_from_spacing
    >>> atoms = bulk('Si', cubic=True)
    >>> 
    >>> # Simple usage - no need for /(2*np.pi)
    >>> kpts = kpts_from_spacing(atoms, 0.20)
    >>> print(kpts)
    (6, 6, 6)
    >>> 
    >>> # This is equivalent to:
    >>> # from ase.io.espresso import kspacing_to_grid
    >>> # kpts = kspacing_to_grid(atoms, 0.20/(2*np.pi))
    
    Notes
    -----
    This function automatically handles the 2π normalization required by
    ase.io.espresso.kspacing_to_grid, providing a cleaner API for users.
    """
    # Apply 2π normalization internally
    kpts = kspacing_to_grid(atoms, kspacing / (2 * np.pi))
    return tuple(kpts)


# ====================================================
# Magnetic configuration helpers
# ====================================================


def setup_magnetic_config(atoms, magnetic_config, pseudopotentials=None, expand_cell=False, 
                         qe_version=None, hubbard_format='auto', projector='ortho-atomic'):
    """
    Simplified and intuitive way to set up magnetic configurations by element.
    
    This function allows you to specify magnetic moments per element type, automatically
    handling species creation, supercell expansion if needed, and Hubbard parameters.
    
    Parameters
    ----------
    atoms : ase.Atoms
        The atomic structure
    magnetic_config : dict
        Magnetic configuration per element. Format:
        {'Element': [mag1, mag2, ...], 'Element2': [mag], ...}
        
        Examples:
        - {'Fe': [1]} - All Fe atoms have magnetization 1 (equivalent)
        - {'Fe': [1, -1]} - Two Fe atoms with different magnetizations (non-equivalent)
        - {'Fe': [1, 1], 'Mn': [1, -1]} - Multiple elements with different configs
        
        If more moments are specified than atoms exist in the cell, and expand_cell=True,
        the cell will be expanded to accommodate the configuration.
        
        Special syntax for Hubbard parameters:
        
        OLD FORMAT (QE < 7.0):
        - {'Fe': {'mag': [1, -1], 'U': 4.3}} - Include Hubbard U
        - {'Fe': {'mag': [1, -1], 'U': [4.3, 4.5]}} - Different U for each species
        
        NEW FORMAT (QE >= 7.0 - HUBBARD card):
        - {'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}} - U on Fe-3d orbital
        - {'Fe': {'mag': [1, -1], 'U': {'3d': [4.3, 4.5]}}} - Different U per species
        - {'Fe': {'mag': [1], 'U': {'3d': 4.3}, 'V': [{'species2': 'O', 'orbital2': '2p', 'value': 1.0}]}}
          - Inter-site V interaction between Fe-3d and O-2p
    
    pseudopotentials : dict, optional
        Pseudopotentials dict mapping elements to UPF files.
        Only base element names are needed (e.g., 'Fe'); derived species 
        (e.g., 'Fe1', 'Fe2') automatically inherit the pseudopotential from 
        their base element. Example: {'Fe': 'Fe.pbe-spn.UPF', 'O': 'O.pbe.UPF'}
    
    expand_cell : bool, default=False
        If True and more moments specified than atoms exist, expand the cell.
        If False and mismatch occurs, raise an error.
    
    qe_version : str, optional
        QE version string (e.g., '7.2', '6.8'). Used to determine format.
        If not provided, uses hubbard_format parameter.
    
    hubbard_format : str, default='auto'
        Hubbard parameter format: 'auto', 'old', or 'new'.
        - 'auto': Determine from qe_version or U format
        - 'old': Use old format (Hubbard_U in SYSTEM namelist)
        - 'new': Use new format (HUBBARD card)
    
    projector : str, default='ortho-atomic'
        Projector type for new Hubbard format. One of:
        - 'ortho-atomic': Orthogonalized atomic orbitals (recommended by QE)
        - 'atomic': Atomic orbitals
        - 'norm-atomic': Normalized atomic orbitals
        - 'wf': Wannier functions
        - 'pseudo': Pseudopotential orbitals
    
    Returns
    -------
    dict
        A dictionary containing:
        - 'atoms': Updated atoms object (may be supercell if expanded)
        - 'input_ntyp': dict with starting_magnetization and optionally Hubbard_U
        - 'pseudopotentials': dict mapping species to pseudopotential files
        - 'species_map': dict mapping species labels to their base element symbols
        - 'expanded': bool indicating if cell was expanded
    
    Examples
    --------
    # Example 1: Simple AFM Fe - all equivalent
    >>> from ase.build import bulk
    >>> atoms = bulk('Fe', cubic=True)  # 2 Fe atoms
    >>> config = setup_magnetic_config(atoms, {'Fe': [1, -1]})
    # Creates Fe1 (mag=1) and Fe2 (mag=-1)
    
    # Example 2: FM with all equivalent
    >>> atoms = bulk('Fe', cubic=True)
    >>> config = setup_magnetic_config(atoms, {'Fe': [1]})
    # Both Fe atoms get mag=1, same species
    
    # Example 3: Complex FeMnAl2 system
    >>> # Assume 2 Fe, 2 Mn, 4 Al in unit cell
    >>> config = setup_magnetic_config(atoms, {
    ...     'Fe': [1],        # Both Fe equivalent, mag=1
    ...     'Mn': [1, -1],    # 2 Mn non-equivalent, AFM
    ...     'Al': [0]         # Al non-magnetic
    ... })
    
    # Example 4: With Hubbard U (old format)
    >>> config = setup_magnetic_config(atoms, {
    ...     'Fe': {'mag': [1, -1], 'U': 4.3}
    ... })
    
    # Example 5: With Hubbard U (new format QE 7.x)
    >>> config = setup_magnetic_config(atoms, {
    ...     'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}
    ... }, qe_version='7.2')
    
    # Example 5b: Specify custom projector (default is 'ortho-atomic')
    >>> config = setup_magnetic_config(atoms, {
    ...     'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}
    ... }, qe_version='7.2', projector='atomic')
    
    # Example 5c: With pseudopotentials (only base element needed)
    >>> pseudopotentials = {'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'}
    >>> config = setup_magnetic_config(atoms, {'Fe': [1, -1]},
    ...                                pseudopotentials=pseudopotentials)
    # Both Fe1 and Fe2 species will use 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
    
    # Example 6: Expand cell if needed
    >>> atoms = bulk('Fe', cubic=True)  # 2 Fe
    >>> config = setup_magnetic_config(
    ...     atoms, 
    ...     {'Fe': [1, 1, -1, -1]},  # Need 4 Fe
    ...     expand_cell=True
    ... )
    # Will create 2x1x1 supercell
    """
    from ase import Atoms
    
    # Validate projector parameter
    valid_projectors = ['ortho-atomic', 'atomic', 'norm-atomic', 'wf', 'pseudo']
    if projector not in valid_projectors:
        raise ValueError(
            f"Invalid projector '{projector}'. Must be one of: {', '.join(valid_projectors)}"
        )
    
    # Determine Hubbard format
    use_new_hubbard_format = False
    if hubbard_format == 'new':
        use_new_hubbard_format = True
    elif hubbard_format == 'old':
        use_new_hubbard_format = False
    elif qe_version:
        # Auto-detect from version
        try:
            major, minor = map(int, qe_version.split('.')[:2])
            use_new_hubbard_format = (major >= 7)
        except (ValueError, AttributeError):
            pass
    
    # Parse magnetic configuration
    element_mags = {}
    element_hubbard = {}
    element_hubbard_v = {}
    
    for element, config in magnetic_config.items():
        if isinstance(config, dict):
            # Extended format with Hubbard parameters
            element_mags[element] = config.get('mag', config.get('magnetization', [0]))
            if 'U' in config:
                element_hubbard[element] = config['U']
            if 'V' in config:
                element_hubbard_v[element] = config['V']
            
            # Auto-detect format if not specified and no qe_version override
            if hubbard_format == 'auto' and not qe_version and 'U' in config:
                u_val = config['U']
                if isinstance(u_val, dict):
                    # New format: {'3d': 4.3}
                    use_new_hubbard_format = True
        elif isinstance(config, (list, tuple)):
            # Simple list of magnetic moments
            element_mags[element] = list(config)
        else:
            # Single value
            element_mags[element] = [config]
    
    # Count atoms per element in current cell
    element_counts = {}
    for atom in atoms:
        element = atom.symbol
        element_counts[element] = element_counts.get(element, 0) + 1
    
    # Check if we need to expand the cell
    needs_expansion = False
    expansion_factors = {}
    
    for element, mags in element_mags.items():
        if element not in element_counts:
            raise ValueError(f"Element {element} not found in atoms structure")
        
        current_count = element_counts[element]
        needed_count = len(mags)
        
        if needed_count > current_count:
            if not expand_cell:
                raise ValueError(
                    f"Element {element}: {needed_count} magnetic moments specified "
                    f"but only {current_count} atoms exist. Set expand_cell=True to auto-expand."
                )
            needs_expansion = True
            factor = (needed_count + current_count - 1) // current_count  # Ceiling division
            expansion_factors[element] = factor
        elif needed_count < current_count:
            # If fewer moments than atoms, replicate the pattern
            # e.g., Fe=[1] with 2 Fe atoms -> both get mag=1
            full_mags = []
            for i in range(current_count):
                full_mags.append(mags[i % len(mags)])
            element_mags[element] = full_mags
    
    # Expand cell if needed
    expanded = False
    if needs_expansion:
        max_factor = max(expansion_factors.values())
        # Simple expansion along first direction
        atoms = atoms * (max_factor, 1, 1)
        expanded = True
        
        # Recount atoms
        element_counts = {}
        for atom in atoms:
            element = atom.symbol
            element_counts[element] = element_counts.get(element, 0) + 1
        
        # Verify expansion is sufficient
        for element, mags in element_mags.items():
            if len(mags) > element_counts[element]:
                raise ValueError(
                    f"Element {element}: Even after expansion, not enough atoms. "
                    f"Need {len(mags)}, have {element_counts[element]}"
                )
    
    # Now assign magnetic moments to atoms
    # Build atom index to magnetic moment mapping
    element_atom_indices = {}
    for i, atom in enumerate(atoms):
        element = atom.symbol
        if element not in element_atom_indices:
            element_atom_indices[element] = []
        element_atom_indices[element].append(i)
    
    mag_dict = {}
    for element, mags in element_mags.items():
        if element in element_atom_indices:
            atom_indices = element_atom_indices[element]
            for i, atom_idx in enumerate(atom_indices):
                if i < len(mags):
                    mag_dict[atom_idx] = mags[i]
                else:
                    # Replicate pattern if more atoms than moments
                    mag_dict[atom_idx] = mags[i % len(mags)]
    
    # Use existing set_magnetic_moments to create species
    result = set_magnetic_moments(atoms, mag_dict, pseudopotentials)
    
    # Add Hubbard parameters if specified
    if element_hubbard or element_hubbard_v:
        if use_new_hubbard_format:
            # NEW FORMAT: Use HUBBARD card (QE 7.x+)
            # Build hubbard dict for new format
            hubbard_dict = {
                'projector': projector,  # Use user-provided projector
                'u': {},
                'v': []
            }
            
            # Process U parameters
            for element, u_value in element_hubbard.items():
                # Find all species for this element
                element_species = [sp for sp in result['species_map'].keys() 
                                 if result['species_map'][sp] == element]
                
                if isinstance(u_value, dict):
                    # New format: {'3d': 4.3} or {'3d': [4.3, 4.5]}
                    for orbital, val in u_value.items():
                        if isinstance(val, (list, tuple)):
                            # Different U for each species
                            for i, species in enumerate(element_species):
                                if i < len(val):
                                    hubbard_dict['u'][f"{species}-{orbital}"] = val[i]
                                else:
                                    hubbard_dict['u'][f"{species}-{orbital}"] = val[-1]
                        else:
                            # Same U for all species of this element
                            for species in element_species:
                                hubbard_dict['u'][f"{species}-{orbital}"] = val
                else:
                    # Old-style value provided but new format requested
                    # Need orbital information - raise error
                    raise ValueError(
                        f"Element {element}: New Hubbard format requires orbital specification. "
                        f"Use 'U': {{'3d': {u_value}}} instead of 'U': {u_value}"
                    )
            
            # Process V parameters
            for element, v_list in element_hubbard_v.items():
                element_species = [sp for sp in result['species_map'].keys() 
                                 if result['species_map'][sp] == element]
                
                if not isinstance(v_list, list):
                    v_list = [v_list]
                
                for v_param in v_list:
                    if not isinstance(v_param, dict):
                        raise ValueError(
                            f"V parameters must be dictionaries with 'species2', 'orbital2', 'value' keys"
                        )
                    
                    species2 = v_param.get('species2')
                    orbital1 = v_param.get('orbital1', v_param.get('orbital'))  # Allow 'orbital' as shorthand
                    orbital2 = v_param.get('orbital2')
                    value = v_param.get('value')
                    i = v_param.get('i', 1)
                    j = v_param.get('j', 1)
                    
                    if not all([species2, orbital1, orbital2, value]):
                        raise ValueError(
                            f"V parameter must specify: species2, orbital1 (or orbital), orbital2, value"
                        )
                    
                    # Add V for each species of this element
                    for species1 in element_species:
                        hubbard_dict['v'].append({
                            'species1': species1,
                            'orbital1': orbital1,
                            'species2': species2,
                            'orbital2': orbital2,
                            'i': i,
                            'j': j,
                            'value': value
                        })
            
            # Store in result for later use in input file generation
            result['hubbard'] = hubbard_dict
            result['qe_version'] = qe_version if qe_version else '7.0'
            
        else:
            # OLD FORMAT: Use input_ntyp (QE < 7.0)
            if 'Hubbard_U' not in result['input_ntyp']:
                result['input_ntyp']['Hubbard_U'] = {}
            
            # Map Hubbard U to species
            for element, u_value in element_hubbard.items():
                # Find all species for this element
                element_species = [sp for sp in result['species_map'].keys() 
                                 if result['species_map'][sp] == element]
                
                if isinstance(u_value, dict):
                    # User specified orbital but we're in old format
                    # Extract first value and warn
                    first_orbital, first_val = next(iter(u_value.items()))
                    print(f"Warning: Orbital specification ('{first_orbital}') ignored in old Hubbard format")
                    u_value = first_val
                
                if isinstance(u_value, (list, tuple)):
                    # Different U for each species
                    for i, species in enumerate(element_species):
                        if i < len(u_value):
                            result['input_ntyp']['Hubbard_U'][species] = u_value[i]
                        else:
                            result['input_ntyp']['Hubbard_U'][species] = u_value[-1]
                else:
                    # Same U for all species of this element
                    for species in element_species:
                        result['input_ntyp']['Hubbard_U'][species] = u_value
            
            # V parameters in old format
            if element_hubbard_v:
                print("Warning: V parameters specified but old Hubbard format may not fully support them")
    
    result['atoms'] = atoms
    result['expanded'] = expanded
    result['hubbard_format'] = 'new' if use_new_hubbard_format else 'old'
    
    return result


def set_magnetic_moments(atoms, magnetic_moments, pseudopotentials=None):
    """
    Simplified way to set up magnetic moments for spin-polarized calculations.
    
    Automatically creates species labels and input_ntyp dictionary for 
    starting_magnetization based on the provided magnetic moments.
    
    Parameters
    ----------
    atoms : ase.Atoms
        The atomic structure
    magnetic_moments : list, dict, or array
        Magnetic moments for each atom. Can be:
        - list/array: magnetic moment for each atom in order
        - dict: {atom_index: magnetic_moment} for specific atoms
    pseudopotentials : dict, optional
        Existing pseudopotentials dict. If provided, it will be updated
        with new species labels. Otherwise, returns the species mapping.
    
    Returns
    -------
    dict
        A dictionary containing:
        - 'input_ntyp': dict with starting_magnetization
        - 'pseudopotentials': dict mapping species to pseudopotential files
        - 'species_map': dict mapping species labels to their base element symbols
                         (e.g., {'Fe1': 'Fe', 'Fe2': 'Fe'} or {'Fe': 'Fe'} if single species)
    
    Examples
    --------
    # Example 1: Simple antiferromagnetic Fe
    >>> from ase.build import bulk
    >>> from xespresso import Espresso
    >>> atoms = bulk('Fe', cubic=True)
    >>> mag_config = set_magnetic_moments(atoms, [1.0, -1.0])
    >>> calc = Espresso(
    ...     pseudopotentials=mag_config['pseudopotentials'],
    ...     input_data={'input_ntyp': mag_config['input_ntyp']},
    ...     nspin=2
    ... )
    
    # Example 2: Specific atoms with magnetic moments
    >>> atoms = bulk('Fe', cubic=True) * (2, 2, 1)
    >>> mag_config = set_magnetic_moments(atoms, {0: 1.0, 1: -1.0, 2: -1.0, 3: 1.0})
    """
    # Ensure species array exists
    if 'species' not in atoms.arrays:
        atoms.new_array('species', np.array(atoms.get_chemical_symbols(), dtype='U20'))
    
    # Convert magnetic_moments to dict format
    if isinstance(magnetic_moments, (list, np.ndarray)):
        mag_dict = {i: mag for i, mag in enumerate(magnetic_moments)}
    else:
        mag_dict = magnetic_moments
    
    # Group atoms by (symbol, magnetic_moment)
    species_groups = {}
    species_counter = {}
    
    # First pass: identify all unique species
    for i in range(len(atoms)):
        symbol = atoms[i].symbol
        mag = mag_dict.get(i, 0.0)
        
        key = (symbol, mag)
        if key not in species_groups:
            # Count how many species of this element we already have
            if symbol not in species_counter:
                species_counter[symbol] = 0
            else:
                species_counter[symbol] += 1
            species_groups[key] = species_counter[symbol]
    
    # Now determine if we need to number species
    # If an element has multiple species, number all of them (Fe1, Fe2, etc)
    # If an element has only one species, keep it unnumbered (Fe)
    element_species_count = {}
    for (symbol, mag) in species_groups.keys():
        element_species_count[symbol] = element_species_count.get(symbol, 0) + 1
    
    # Create final species labels
    final_species_labels = {}
    for (symbol, mag), counter in species_groups.items():
        if element_species_count[symbol] > 1:
            # Multiple species of this element - number all (Fe1, Fe2, ...)
            species_label = f"{symbol}{counter + 1}"
        else:
            # Only one species of this element - no number (Fe)
            species_label = symbol
        final_species_labels[(symbol, mag)] = species_label
    
    # Second pass: assign final labels to atoms
    for i in range(len(atoms)):
        symbol = atoms[i].symbol
        mag = mag_dict.get(i, 0.0)
        key = (symbol, mag)
        atoms.arrays['species'][i] = final_species_labels[key]
    
    # Create input_ntyp dictionary for starting_magnetization
    input_ntyp = {'starting_magnetization': {}}
    species_map = {}
    
    for key, species_label in final_species_labels.items():
        symbol, mag = key
        if mag != 0.0:
            input_ntyp['starting_magnetization'][species_label] = mag
        # Map species_label back to base element symbol
        species_map[species_label] = symbol
    
    # Handle pseudopotentials
    if pseudopotentials is None:
        # Return template that user should fill
        pseudo_dict = {}
        for species_label in final_species_labels.values():
            # Get the base element symbol
            base_symbol = ''.join([c for c in species_label if not c.isdigit()])
            pseudo_dict[species_label] = f"{base_symbol}.UPF"
    else:
        # Update existing pseudopotentials
        pseudo_dict = pseudopotentials.copy()
        for key, species_label in final_species_labels.items():
            symbol, mag = key
            if species_label not in pseudo_dict:
                # Try to find pseudopotential for base element
                base_pseudo = pseudopotentials.get(symbol)
                if base_pseudo:
                    pseudo_dict[species_label] = base_pseudo
                else:
                    pseudo_dict[species_label] = f"{symbol}.UPF"
    
    return {
        'input_ntyp': input_ntyp,
        'pseudopotentials': pseudo_dict,
        'species_map': species_map
    }


def set_antiferromagnetic(atoms, sublattice_indices, magnetic_moment=1.0, pseudopotentials=None):
    """
    Simplified way to set up antiferromagnetic configurations.
    
    Divides atoms into two sublattices with opposite magnetic moments.
    
    Parameters
    ----------
    atoms : ase.Atoms
        The atomic structure
    sublattice_indices : list of lists
        Two sublattices: [sublattice_A_indices, sublattice_B_indices]
        Example: [[0, 2], [1, 3]] for 4 atoms in checkerboard AFM
    magnetic_moment : float, default=1.0
        Magnitude of magnetic moment (absolute value)
    pseudopotentials : dict, optional
        Existing pseudopotentials dict
    
    Returns
    -------
    dict
        Same as set_magnetic_moments: input_ntyp, pseudopotentials, species_map
    
    Examples
    --------
    # Simple AFM with alternating spins
    >>> from ase.build import bulk
    >>> from xespresso import Espresso
    >>> atoms = bulk('Fe', cubic=True)
    >>> afm_config = set_antiferromagnetic(atoms, [[0], [1]])
    >>> calc = Espresso(
    ...     pseudopotentials=afm_config['pseudopotentials'],
    ...     input_data={'input_ntyp': afm_config['input_ntyp']},
    ...     nspin=2
    ... )
    """
    if len(sublattice_indices) != 2:
        raise ValueError("sublattice_indices must contain exactly 2 sublattices")
    
    # Create magnetic moments array
    mag_dict = {}
    for idx in sublattice_indices[0]:
        mag_dict[idx] = magnetic_moment
    for idx in sublattice_indices[1]:
        mag_dict[idx] = -magnetic_moment
    
    return set_magnetic_moments(atoms, mag_dict, pseudopotentials)


def set_ferromagnetic(atoms, magnetic_moment=1.0, element=None, pseudopotentials=None):
    """
    Simplified way to set up ferromagnetic configurations.
    
    Sets all atoms (or all atoms of a specific element) to the same magnetic moment.
    
    Parameters
    ----------
    atoms : ase.Atoms
        The atomic structure
    magnetic_moment : float, default=1.0
        Magnetic moment for all atoms
    element : str, optional
        If specified, only set magnetic moment for atoms of this element
    pseudopotentials : dict, optional
        Existing pseudopotentials dict
    
    Returns
    -------
    dict
        Same as set_magnetic_moments: input_ntyp, pseudopotentials, species_map
    
    Examples
    --------
    # Ferromagnetic Fe
    >>> from ase.build import bulk
    >>> from xespresso import Espresso
    >>> atoms = bulk('Fe', cubic=True)
    >>> fm_config = set_ferromagnetic(atoms, magnetic_moment=2.0)
    >>> calc = Espresso(
    ...     pseudopotentials=fm_config['pseudopotentials'],
    ...     input_data={'input_ntyp': fm_config['input_ntyp']},
    ...     nspin=2
    ... )
    """
    mag_dict = {}
    for i, atom in enumerate(atoms):
        if element is None or atom.symbol == element:
            mag_dict[i] = magnetic_moment
    
    return set_magnetic_moments(atoms, mag_dict, pseudopotentials)


# ====================================================


def get_nbnd(atoms=None, scale=1.2, pseudopotentials={}, nspin=1, input_data={}):
    input_parameters = construct_namelist(input_data)
    atomic_species_str, species_info, total_valence = build_atomic_species_str(
        atoms, input_parameters, pseudopotentials
    )
    nbnd = int(total_valence / (2.0 / nspin))
    nbnd_scale = int(nbnd * scale)
    print(
        " total valence: %s\n nbnd: %s\n scaled nbnd: %s"
        % (total_valence, nbnd, nbnd_scale)
    )
    return nbnd


def merge_slab(slab1, slab2, index=2):
    """ """
    slab2.cell[index] = slab1.cell[index]
    slab2.set_cell(slab1.cell, scale_atoms=True)
    slab1 = slab1 + slab2
    slab1.cell[2][2] = max(slab1.positions[:, 2] + 15)
    slab1.wrap()
    return slab1


def qeinp(
    calculation,
    ecutwfc=30,
    mixing_beta=0.5,
    conv_thr=1.0e-8,
    edir=False,
    input_ntyp={},
    atoms=None,
):
    #
    inp = {
        # control
        "calculation": calculation,
        "max_seconds": 78000,
        "verbosity": "high",
        "tprnfor": True,
        # system
        "ecutwfc": ecutwfc,
        "ecutrho": ecutwfc * 8,
        "occupations": "smearing",
        "degauss": 0.01,
        "input_ntyp": input_ntyp,
        # electrons
        "mixing_beta": mixing_beta,
        "conv_thr": conv_thr,
        "electron_maxstep": 400,
    }
    #
    if edir:
        inp.update(dipole_correction(atoms, edir))

    return inp


def dipole_correction(atoms, edir=3):
    inp = {
        "dipfield": True,
        "tefield": True,
        "edir": edir,
        "eamp": 0.001,
        "eopreg": 0.05,
    }
    emaxpos = (
        max(atoms.positions[:, edir - 1] + atoms.cell[edir - 1][edir - 1])
        / 2.0
        / atoms.cell[edir - 1][edir - 1]
    )
    inp["emaxpos"] = round(emaxpos, 2)
    return inp


# tools


def build_oer(atoms):
    """ """
    from ase.atoms import Atoms

    # ----------------------------------
    ooh = Atoms("O2H", positions=[[0, 0, 0], [1.4, 0, 0], [1.4, 0, 1.0]])
    ooh.rotate("z", np.pi / 4)
    mols = {
        "o": Atoms("O"),
        "oh": Atoms("OH", positions=[[0, 0, 0], [0, 0, 1.0]]),
        "ooh": ooh,
    }
    #
    jobs = {}
    maxz = max(atoms.positions[:, 2])
    indm = [atom.index for atom in atoms if atom.z > maxz - 1.0 and atom.symbol != "O"][
        0
    ]
    for job, mol in mols.items():
        # print(job, mol)
        ads = mol.copy()
        natoms = atoms.copy()
        ads.translate(atoms[indm].position - ads[0].position + [0, 0, 1.9])
        natoms = natoms + ads
        jobs[job] = natoms
    return jobs


def fix_layers(atoms, miller=(0, 0, 1), tol=1.0, n=[0, 4]):
    """ """
    layers = get_layers(atoms, miller, tol)[0]
    index = [j for j in range(len(atoms)) if layers[j] in range(n[0], n[1])]
    constraint = FixAtoms(indices=index)
    atoms.set_constraint(constraint)
    return atoms


def mypool(jobs, func, showInfo=False):
    """ """
    from random import random
    from time import sleep

    pool = multiprocessing.Pool(len(jobs))
    results = []
    images = []
    for job, atoms in jobs.items():
        if showInfo:
            print(job, len(atoms), atoms)
        sleep(random() * 2)
        results.append(pool.apply_async(func, (job, atoms)))
    for r in results:
        r.get()
    pool.close()
    pool.join()


def dwubelix(updates=[]):
    file = "pw err out dos pdos projwfc int xyz path a.xml txt png"
    print("Downloading.....")
    cwd = os.getcwd()
    for update in updates:
        os.chdir(update)
        os.system("dwubelix-sc.py %s" % file)
        os.chdir(cwd)
    print("Finished")


def ana(dire, calc):
    atoms = calc.results["atoms"]
    results = [dire, atoms, atoms.cell, atoms.positions]
    for prop in ["energy", "forces", "stress", "magmoms"]:
        if prop in calc.results:
            prop = calc.results[prop]
        else:
            prop = None
        results.append(prop)
    return results


# tools


def summary(updates=[], prefix="datas"):
    import pandas as pd

    columns = ["label", "atoms", "cell", "positions", "energy", "forces", "stress"]
    file = "%s.pickle" % prefix
    db = "%s.db" % prefix
    if os.path.exists(file):
        with open(file, "rb") as f:
            datas, df = pickle.load(f)
    else:
        datas = {}
        df = pd.DataFrame(columns=columns)
    calc = Espresso()
    print("Reading.....")
    for update in updates:
        cwd = os.getcwd()
        for i, j, y in os.walk(update):
            output = is_espresso(i)
            if output:
                os.chdir(i)
                print("=" * 30)
                print("Reading dire:", i)
                calc.directory = cwd + "/" + i
                calc.prefix = output[0:-4]
                try:
                    calc.results = {}
                    calc.read_results()
                    datas[i] = calc.results
                    atoms = calc.results["atoms"]
                    atoms.write(os.path.join(calc.directory, "%s.cif" % calc.prefix))
                except Exception as e:
                    print("error: %s \n" % e)
            os.chdir(cwd)
    with open(file, "wb") as f:
        pickle.dump([datas, df], f)
    print("Finished")


def is_espresso(path):
    """
    check espresso
    """
    dirs = os.listdir(path)
    # print(dirs)
    # flag = True

    for qefile in [".pwi"]:
        flag = False
        for file in dirs:
            if qefile in file:
                return file
        if not flag:
            return False
    # return flag


def grep_valence_configuration(pseudopotential):
    """
    Given a UPF pseudopotential file, find the valence configuration.

    Valence configuration:
    nl pn  l   occ       Rcut    Rcut US       E pseu
    3S  1  0  2.00      0.700      1.200    -6.910117

    """
    orbitals = ["S", "P", "D", "F"]
    valence = {}
    with open(pseudopotential) as psfile:
        lines = psfile.readlines()
        for i in range(len(lines)):
            if "valence configuration:" in lines[i].lower():
                j = i + 2
                ob = lines[j].split()[0]
                while ob[1] in orbitals:
                    valence[ob] = lines[j].split()[3]
                    j += 1
                    ob = lines[j].split()[0]
                return valence
    if not valence:
        raise ValueError("Valence configuration missing in {}".format(pseudopotential))
