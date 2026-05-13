from xespresso.post.base import PostCalculation


class EspressoPw2wannier90(PostCalculation):
    """
    Interface to pw2wannier90.x for converting QE wavefunctions to Wannier90 format.
    
    This class handles:
    - Generation of pw2wannier.in input files
    - Execution of pw2wannier90.x conversion
    - Collection of output files (amn, mmn, eig, unk)
    
    The workflow is typically:
    1. Run NSCF calculation to get wavefunctions
    2. Use EspressoPw2wannier90 to convert them
    3. Pass outputs to EspressoWannier90 for function generation
    
    Example:
        >>> pw2wan = EspressoPw2wannier90(
        ...     parent_directory='/path/to/nscf_output',
        ...     prefix='Fe_bcc',           # Read from Fe_bcc.save
        ...     seedname='wannier_seed',   # Output files: wannier_seed.amn, etc.
        ...     directory='./runs/05-wan'  # Where to generate input file
        ... )
        >>> result = pw2wan.run()
    """
    
    package = "pw2wannier90"
    # QE 7.5 official documentation: all &INPUTPP parameters
    package_parameters = {
        "inputpp": [
            # Basic I/O
            "prefix",
            "outdir",
            "seedname",
            # Spin-related
            "spin_component",
            "wan_mode",
            # UNK file control (wavefunction real-space plotting)
            "write_unk",
            "reduce_unk",
            "reduce_unk_factor",
            "wvfn_formatted",
            # Main output matrices
            "write_amn",
            "write_mmn",
            # SCDM projection method
            "scdm_proj",
            "scdm_entanglement",
            "scdm_mu",
            "scdm_sigma",
            # Atomic projection method
            "atom_proj",
            "atom_proj_exclude",
            "atom_proj_ext",
            "atom_proj_dir",
            "atom_proj_ortho",
            # Additional matrices
            "write_spn",
            "spn_formatted",
            "write_uHu",
            "uHu_formatted",
            "write_uIu",
            "uIu_formatted",
            "write_sHu",
            "sHu_formatted",
            "write_sIu",
            "sIu_formatted",
            "write_unkg",
            # Symmetry
            "irr_bz",
            "write_dmn",
            "read_sym",
        ]
    }
    
    def __init__(
        self,
        parent_directory,
        prefix,
        seedname='wannier',
        spin_component='none',
        wan_mode='standalone',
        write_amn=True,
        write_mmn=True,
        write_unk=False,
        write_spn=False,
        write_unkg=False,
        atom_proj=True,
        atom_proj_ortho=True,
        irr_bz=False,
        write_dmn=False,
        queue=None,
        parallel='',
        debug=False,
        dry_run=False,
        directory=None,
        **kwargs,
    ):
        """
        Initialize EspressoPw2wannier90 with QE 7.5 parameters.
        
        Parameters (QE official documentation):
        -----------
        parent_directory : str
            Directory where QE NSCF outputs are located (prefix.save)
        prefix : str
            Prefix for NSCF calculation (e.g., 'Fe_bcc' reads Fe_bcc.save)
        seedname : str, optional
            Base name for Wannier90 output files (default 'wannier').
            For spin-polarized calculations, '_up' or '_dn' suffixes are 
            automatically appended based on spin_component:
            - spin_component='up' → seedname_up.amn, seedname_up.mmn, etc.
            - spin_component='down' → seedname_dn.amn, seedname_dn.mmn, etc.
            - spin_component='none' → seedname.amn, seedname.mmn, etc.
        spin_component : str, optional
            Spin component for collinear/non-collinear calculations (default 'none').
            Valid values:
            - 'up': spin up component (collinear spin calculation)
            - 'down': spin down component (collinear spin calculation)
            - 'none': no-spin or non-collinear calculation
        wan_mode : str, optional
            'standalone': standalone execution, 'library': library mode (default 'standalone')
        write_amn : bool, optional
            Write A(k) projection matrix (default True). Set False if not required.
        write_mmn : bool, optional
            Write M(k,b) overlap matrix (default True). Set False if not required.
        write_unk : bool, optional
            Write periodic part of Bloch functions (default False). 
            WARNING: Creates large files (~GB), only set True for plotting.
        write_spn : bool, optional
            Write spin operator matrix elements (non-collinear only, default False)
        write_unkg : bool, optional
            Write first few Fourier components of periodic Bloch parts (default False)
        atom_proj : bool, optional
            Use pseudo-atomic wavefunctions as initial projection (default True).
            This is the recommended method for Wannier90 initial projections.
        atom_proj_ortho : bool, optional
            Orthonormalize pseudo-atomic wavefunctions before computing 
            inner product with Bloch states (default True).
            Recommended to keep True. Only set to False if you know what you are doing.
            Only relevant if atom_proj=True.
        irr_bz : bool, optional
            Use irreducible BZ for amn/mmn files (default False). 
            Changes output extensions to iamn/immn/ieig.
        write_dmn : bool, optional
            Construct symmetry-adapted Wannier functions (default False)
        queue : dict, optional
            Scheduler configuration (local/remote job submission)
        parallel : str, optional
            Parallelization flags for pw2wannier90.x (e.g. '-nk 4 -nd 2')
        debug : bool, optional
            Enable debug logging
        dry_run : bool, optional
            If True, only generate input without execution
        directory : str, optional
            Directory for pw2wannier.in file (uses parent_directory if None)
        **kwargs : dict
            Additional parameters (SCDM, atom_proj_*, reduce_unk, etc.)
        """
        # Validate spin_component first
        valid_spin_components = ['none', 'up', 'down']
        if spin_component not in valid_spin_components:
            raise ValueError(
                f"Invalid spin_component '{spin_component}'. "
                f"Valid values: {valid_spin_components}"
            )
        
        # Append spin suffix to seedname if spin-polarized
        if spin_component == 'up':
            seedname_full = f"{seedname}_up"
            prefix_full = f"{prefix}_up"
        elif spin_component == 'down':
            seedname_full = f"{seedname}_dn"
            prefix_full = f"{prefix}_dn"
        else:
            seedname_full = seedname
            prefix_full = prefix
        
        self.seedname = seedname
        self.seedname_full = seedname_full  # With spin suffix if applicable
        self.prefix_base = prefix  # Store original prefix for &inputpp parameter
        self.spin_component = spin_component
        self.parallel = parallel
        self.dry_run = dry_run
        
        # Call parent __init__ with prefix_full (includes _up/_dn suffix)
        # This generates Fe_bcc_up.post_asei, Fe_bcc_dn.post_asei, etc.
        super().__init__(
            parent_directory=parent_directory,
            prefix=prefix_full,
            queue=queue,
            debug=debug,
            directory=directory,
            **kwargs,
        )
        
        # Set main defaults (QE 7.5 documentation)
        if 'seedname' not in self.parameters:
            self.parameters['seedname'] = seedname_full
        # Override prefix to use base (without spin suffix) for &inputpp
        self.parameters['prefix'] = self.prefix_base
        if 'spin_component' not in self.parameters:
            self.parameters['spin_component'] = spin_component
        if 'wan_mode' not in self.parameters:
            self.parameters['wan_mode'] = wan_mode
            
        # Write flags with proper defaults
        if 'write_amn' not in self.parameters:
            self.parameters['write_amn'] = '.true.' if write_amn else '.false.'
        if 'write_mmn' not in self.parameters:
            self.parameters['write_mmn'] = '.true.' if write_mmn else '.false.'
        if 'write_unk' not in self.parameters:
            self.parameters['write_unk'] = '.true.' if write_unk else '.false.'
        if 'write_spn' not in self.parameters:
            self.parameters['write_spn'] = '.true.' if write_spn else '.false.'
        if 'write_unkg' not in self.parameters:
            self.parameters['write_unkg'] = '.true.' if write_unkg else '.false.'
            
        # Atomic projection defaults
        if 'atom_proj' not in self.parameters:
            self.parameters['atom_proj'] = '.true.' if atom_proj else '.false.'
        if 'atom_proj_ortho' not in self.parameters:
            self.parameters['atom_proj_ortho'] = '.true.' if atom_proj_ortho else '.false.'
            
        # Warn if atom_proj_ortho=False with atom_proj=True
        if atom_proj and not atom_proj_ortho:
            import warnings
            warnings.warn(
                "atom_proj_ortho=False with atom_proj=True: "
                "It is recommended to keep atom_proj_ortho=True unless you know what you are doing.",
                UserWarning
            )
            
        # Symmetry defaults
        if 'irr_bz' not in self.parameters:
            self.parameters['irr_bz'] = '.true.' if irr_bz else '.false.'
        if 'write_dmn' not in self.parameters:
            self.parameters['write_dmn'] = '.true.' if write_dmn else '.false.'
    
    def write_package_input(self):
        """Write pw2wannier90 input file with proper formatting.
        
        String parameters get quotes, logical parameters don't.
        Files are named with spin suffix: Fe_bcc_up.pw2wannier90i, Fe_bcc_dn.pw2wannier90i
        """
        import os
        filename = os.path.join(self.directory, "%s.%si" % (self.prefix, self.package))
        defaults = self.get_defaults()
        
        with open(filename, "w") as f:
            for section, parameters in self.package_parameters.items():
                if section != "LINE":
                    f.write("&%s\n" % section)
                    for key, value in self.parameters.items():
                        if key in parameters:
                            # Skip if matches default
                            if key in defaults and defaults[key] == value:
                                continue
                            
                            # Logical parameters: .true., .false.
                            if value in ['.true.', '.false.']:
                                f.write("  {0:15s} = {1}\n".format(key, value))
                            # String parameters (all others, including spin_component)
                            else:
                                f.write('  {0:15s} = "{1}"\n'.format(key, value))
                    f.write("/ \n")
    
    def get_defaults(self) -> dict:
        """Get default values for pw2wannier90 parameters per QE 7.5 documentation.
        
        Only parameters that match defaults will be skipped from the input file.
        """
        return {
            'spin_component': 'none',
            'wan_mode': 'standalone',
            'write_unk': '.false.',
            'reduce_unk': '.false.',
            'reduce_unk_factor': 2,  # Only relevant if write_unk=.true.
            'wvfn_formatted': '.false.',
            'write_amn': '.true.',
            'scdm_proj': '.false.',
            'scdm_entanglement': 'isolated',
            'atom_proj': '.false.',
            'atom_proj_ext': '.false.',
            'atom_proj_ortho': '.true.',
            'write_mmn': '.true.',
            'write_spn': '.false.',
            'spn_formatted': '.false.',
            'write_uHu': '.false.',
            'uHu_formatted': '.false.',
            'write_uIu': '.false.',
            'uIu_formatted': '.false.',
            'write_sHu': '.false.',
            'sHu_formatted': '.false.',
            'write_sIu': '.false.',
            'sIu_formatted': '.false.',
            'write_unkg': '.false.',
            'irr_bz': '.false.',
            'write_dmn': '.false.',
            'read_sym': '.false.',
        }
