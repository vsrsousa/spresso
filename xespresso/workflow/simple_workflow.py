"""
Simplified workflow for common Quantum ESPRESSO calculations.

This module provides an easy-to-use interface for running common calculations
like SCF and structure relaxation from CIF files or ASE Atoms objects.
"""

import numpy as np
from typing import Dict, Optional, Union, Tuple, List
from pathlib import Path
from ase import Atoms
from ase.io import read
from ase.io.espresso import kspacing_to_grid
from xespresso import Espresso
from xespresso.tools import setup_magnetic_config
from xespresso.machines import load_machine


# Preset configurations for different calculation protocols
PRESETS = {
    'fast': {
        'ecutwfc': 30.0,
        'ecutrho': 240.0,
        'conv_thr': 1.0e-6,
        'kspacing': 0.5,  # Angstrom^-1
        'mixing_beta': 0.7,
        'electron_maxstep': 100,
    },
    'moderate': {
        'ecutwfc': 50.0,
        'ecutrho': 400.0,
        'conv_thr': 1.0e-8,
        'kspacing': 0.3,  # Angstrom^-1
        'mixing_beta': 0.5,
        'electron_maxstep': 200,
    },
    'accurate': {
        'ecutwfc': 80.0,
        'ecutrho': 640.0,
        'conv_thr': 1.0e-10,
        'kspacing': 0.15,  # Angstrom^-1
        'mixing_beta': 0.3,
        'electron_maxstep': 300,
    }
}


class CalculationWorkflow:
    """
    Simplified workflow for running Quantum ESPRESSO calculations.
    
    This class provides an easy interface for:
    - Reading structures from CIF files
    - Setting up calculations with protocol presets
    - Running SCF or relaxation calculations
    - Using k-spacing instead of explicit k-points
    
    Examples:
        >>> # SCF calculation from CIF file
        >>> workflow = CalculationWorkflow.from_cif(
        ...     'structure.cif',
        ...     pseudopotentials={'Si': 'Si.pbe.UPF'},
        ...     protocol='moderate'
        ... )
        >>> workflow.run_scf(label='scf/silicon')
        
        >>> # Relax calculation with custom k-spacing
        >>> workflow = CalculationWorkflow.from_cif(
        ...     'structure.cif',
        ...     pseudopotentials={'Si': 'Si.pbe.UPF'},
        ...     protocol='fast',
        ...     kspacing=0.4
        ... )
        >>> workflow.run_relax(label='relax/silicon')
    """
    
    def __init__(
        self,
        atoms: Atoms,
        pseudopotentials: Dict[str, str],
        protocol: str = 'moderate',
        kspacing: Optional[float] = None,
        input_data: Optional[Dict] = None,
        magnetic_config: Optional[Union[str, Dict]] = None,
        expand_cell: bool = False,
        queue: Optional[Dict] = None,
        machine: Optional[str] = None,
        **kwargs
    ):
        """
        Initialize a calculation workflow.
        
        Args:
            atoms: ASE Atoms object representing the structure
            pseudopotentials: Dictionary mapping element symbols to pseudopotential files
            protocol: Protocol preset: 'fast', 'moderate', or 'accurate'
            kspacing: K-point spacing in Angstrom^-1 (physical units). If None, uses preset value.
                     The workflow automatically handles the 2π normalization when converting to k-points.
                     Example: kspacing=0.20 will give the same k-points as
                     ase.io.espresso.kspacing_to_grid(atoms, 0.20/(2*np.pi))
            input_data: Additional input parameters (merged with preset)
            magnetic_config: Magnetic configuration. Can be:
                           - 'ferro' or 'ferromagnetic': All atoms ferromagnetic
                           - 'antiferro' or 'antiferromagnetic': Alternating spin
                           - Dict: Element-based config, e.g. {'Fe': [1, -1], 'O': [0]}
                           Also supports Hubbard parameters in the dict format
            expand_cell: If True, expand cell to accommodate magnetic configuration
            queue: Queue configuration dictionary for job submission (local or remote).
                   This is directly passed to the Espresso calculator.
            machine: Name of a machine configuration to load from ~/.xespresso/machines/.
                    If provided, the machine configuration is loaded and converted to a queue dict.
                    Cannot be used together with 'queue' parameter.
            **kwargs: Additional parameters passed to Espresso calculator
        """
        self.atoms = atoms.copy()  # Work with a copy to avoid modifying original
        self.original_pseudopotentials = pseudopotentials
        self.protocol = protocol
        self.extra_kwargs = kwargs
        self.expand_cell = expand_cell
        
        # Handle machine and queue configuration
        if queue is not None and machine is not None:
            raise ValueError(
                "Cannot specify both 'queue' and 'machine' parameters. "
                "Use 'queue' for direct configuration or 'machine' to load from config."
            )
        
        if machine is not None:
            # Load machine configuration
            self.queue = load_machine(machine_name=machine)
        else:
            self.queue = queue
        
        # Get preset configuration
        if protocol not in PRESETS:
            raise ValueError(
                f"Protocol must be one of {list(PRESETS.keys())}, got '{protocol}'"
            )
        
        self.preset = PRESETS[protocol].copy()
        
        # Override k-spacing if provided
        if kspacing is not None:
            self.preset['kspacing'] = kspacing
        
        # Initialize input_data early so it can be used in magnetic config
        self.input_data = self.preset.copy()
        if input_data:
            self.input_data.update(input_data)
        
        # Handle magnetic configuration if provided
        if magnetic_config is not None:
            self._apply_magnetic_config(magnetic_config)
        else:
            self.pseudopotentials = pseudopotentials
        
        # Remove kspacing from input_data as it will be converted to kpts
        self.kspacing = self.input_data.pop('kspacing', None)
    
    def _apply_magnetic_config(self, magnetic_config: Union[str, Dict]):
        """Apply magnetic configuration using setup_magnetic_config."""
        from xespresso.tools import set_ferromagnetic, set_antiferromagnetic
        
        if isinstance(magnetic_config, str):
            magnetic_config = magnetic_config.lower()
            if magnetic_config in ['ferro', 'ferromagnetic']:
                # Simple ferromagnetic configuration
                config = set_ferromagnetic(
                    self.atoms, 
                    magnetic_moment=1.0, 
                    pseudopotentials=self.original_pseudopotentials
                )
                # atoms modified in-place, just update config
                self.pseudopotentials = config.get('pseudopotentials', self.original_pseudopotentials.copy())
                if 'input_ntyp' in config:
                    if 'input_ntyp' not in self.input_data:
                        self.input_data['input_ntyp'] = {}
                    self.input_data['input_ntyp'].update(config['input_ntyp'])
            elif magnetic_config in ['antiferro', 'antiferromagnetic']:
                # Simple antiferromagnetic configuration
                # For antiferromagnetic, we need to determine sublattices
                # Simple approach: alternate atoms
                n_atoms = len(self.atoms)
                sublattice1 = list(range(0, n_atoms, 2))
                sublattice2 = list(range(1, n_atoms, 2))
                
                config = set_antiferromagnetic(
                    self.atoms,
                    sublattice_indices=[sublattice1, sublattice2],
                    magnetic_moment=1.0,
                    pseudopotentials=self.original_pseudopotentials
                )
                # atoms modified in-place, just update config
                self.pseudopotentials = config.get('pseudopotentials', self.original_pseudopotentials.copy())
                if 'input_ntyp' in config:
                    if 'input_ntyp' not in self.input_data:
                        self.input_data['input_ntyp'] = {}
                    self.input_data['input_ntyp'].update(config['input_ntyp'])
            else:
                raise ValueError(
                    f"Unknown magnetic configuration: '{magnetic_config}'. "
                    "Use 'ferro', 'antiferro', or a dict with element-based config."
                )
        elif isinstance(magnetic_config, dict):
            # Element-based configuration with possible Hubbard parameters
            config = setup_magnetic_config(
                self.atoms,
                magnetic_config,
                pseudopotentials=self.original_pseudopotentials,
                expand_cell=self.expand_cell
            )
            self.atoms = config['atoms']
            self.pseudopotentials = config.get('pseudopotentials', self.original_pseudopotentials)
            
            # Merge special input_data from magnetic config
            if 'input_ntyp' in config:
                if 'input_ntyp' not in self.input_data:
                    self.input_data['input_ntyp'] = {}
                self.input_data['input_ntyp'].update(config['input_ntyp'])
            
            # Handle Hubbard parameters in new format
            if 'hubbard' in config:
                self.input_data['hubbard'] = config['hubbard']
            if 'hubbard_v' in config:
                self.input_data['hubbard_v'] = config['hubbard_v']
            if 'qe_version' in config:
                self.input_data['qe_version'] = config.get('qe_version')
            if 'lda_plus_u' in config:
                self.input_data['lda_plus_u'] = config['lda_plus_u']
        else:
            raise TypeError(
                f"magnetic_config must be str or dict, got {type(magnetic_config)}"
            )
    
    @classmethod
    def from_cif(
        cls,
        cif_file: Union[str, Path],
        pseudopotentials: Dict[str, str],
        protocol: str = 'moderate',
        kspacing: Optional[float] = None,
        input_data: Optional[Dict] = None,
        magnetic_config: Optional[Union[str, Dict]] = None,
        expand_cell: bool = False,
        queue: Optional[Dict] = None,
        machine: Optional[str] = None,
        **kwargs
    ) -> 'CalculationWorkflow':
        """
        Create a workflow from a CIF file.
        
        Args:
            cif_file: Path to CIF file
            pseudopotentials: Dictionary mapping element symbols to pseudopotential files
            protocol: Protocol preset: 'fast', 'moderate', or 'accurate'
            kspacing: K-point spacing in Angstrom^-1 (physical units)
            input_data: Additional input parameters
            magnetic_config: Magnetic configuration ('ferro', 'antiferro', or element dict)
            expand_cell: If True, expand cell to accommodate magnetic configuration
            queue: Queue configuration dictionary for job submission
            machine: Name of a machine configuration to load
            **kwargs: Additional parameters passed to Espresso calculator
            
        Returns:
            CalculationWorkflow: Initialized workflow object
        """
        atoms = read(str(cif_file))
        return cls(atoms, pseudopotentials, protocol, kspacing, input_data, 
                   magnetic_config, expand_cell, queue, machine, **kwargs)
    
    def _get_kpts(self) -> Union[Tuple[int, int, int], str]:
        """
        Calculate k-points from k-spacing using ase.io.espresso.kspacing_to_grid.
        
        Returns:
            Tuple of k-points or 'gamma'
        """
        if self.kspacing is not None:
            # Convert kspacing to k-point grid
            # Note: kspacing_to_grid expects spacing in units of 2*pi/Angstrom
            # So we need to convert from Angstrom^-1
            kpts = kspacing_to_grid(self.atoms, self.kspacing / (2 * np.pi))
            return tuple(kpts)
        else:
            # Default to gamma point if no k-spacing specified
            return (1, 1, 1)
    
    def run_scf(
        self,
        label: str = 'scf',
        **calc_kwargs
    ) -> Espresso:
        """
        Run a self-consistent field (SCF) calculation.
        
        Args:
            label: Directory/label for the calculation
            **calc_kwargs: Additional parameters for the Espresso calculator
            
        Returns:
            Espresso: Calculator object with results
        """
        # Prepare parameters
        params = {
            'pseudopotentials': self.pseudopotentials,
            'label': label,
            'calculation': 'scf',
            'input_data': self.input_data.copy(),
            'kpts': self._get_kpts(),
        }
        
        # Add ecutwfc and ecutrho at top level
        params['ecutwfc'] = self.input_data.get('ecutwfc', 50.0)
        params['ecutrho'] = self.input_data.get('ecutrho', 400.0)
        
        # Add queue configuration if provided
        if self.queue is not None:
            params['queue'] = self.queue
        
        # Merge with extra kwargs
        params.update(self.extra_kwargs)
        params.update(calc_kwargs)
        
        # Create and run calculator
        calc = Espresso(**params)
        self.atoms.calc = calc
        calc.run(atoms=self.atoms)
        
        return calc
    
    def run_relax(
        self,
        label: str = 'relax',
        relax_type: str = 'relax',
        **calc_kwargs
    ) -> Espresso:
        """
        Run a structure relaxation calculation.
        
        Args:
            label: Directory/label for the calculation
            relax_type: Type of relaxation: 'relax' (ions only) or 'vc-relax' (ions + cell)
            **calc_kwargs: Additional parameters for the Espresso calculator
            
        Returns:
            Espresso: Calculator object with results
        """
        if relax_type not in ['relax', 'vc-relax']:
            raise ValueError(
                f"relax_type must be 'relax' or 'vc-relax', got '{relax_type}'"
            )
        
        # Prepare parameters
        params = {
            'pseudopotentials': self.pseudopotentials,
            'label': label,
            'calculation': relax_type,
            'input_data': self.input_data.copy(),
            'kpts': self._get_kpts(),
        }
        
        # Add ecutwfc and ecutrho at top level
        params['ecutwfc'] = self.input_data.get('ecutwfc', 50.0)
        params['ecutrho'] = self.input_data.get('ecutrho', 400.0)
        
        # Add queue configuration if provided
        if self.queue is not None:
            params['queue'] = self.queue
        
        # Merge with extra kwargs
        params.update(self.extra_kwargs)
        params.update(calc_kwargs)
        
        # Create and run calculator
        calc = Espresso(**params)
        self.atoms.calc = calc
        calc.run(atoms=self.atoms)
        
        return calc
    
    def get_atoms(self) -> Atoms:
        """Get the current atoms object."""
        return self.atoms
    
    def get_preset_info(self) -> Dict:
        """Get information about the current protocol preset."""
        return {
            'protocol': self.protocol,
            'preset': self.preset,
            'kpts': self._get_kpts(),
            'kspacing': self.kspacing,
        }


def quick_scf(
    structure: Union[str, Path, Atoms],
    pseudopotentials: Dict[str, str],
    label: str = 'scf',
    protocol: str = 'moderate',
    kspacing: Optional[float] = None,
    magnetic_config: Optional[Union[str, Dict]] = None,
    expand_cell: bool = False,
    queue: Optional[Dict] = None,
    machine: Optional[str] = None,
    **kwargs
) -> Espresso:
    """
    Quick SCF calculation helper function.
    
    Args:
        structure: CIF file path or ASE Atoms object
        pseudopotentials: Dictionary mapping element symbols to pseudopotential files
        label: Directory/label for the calculation
        protocol: Protocol preset: 'fast', 'moderate', or 'accurate'
        kspacing: K-point spacing in Angstrom^-1 (physical units)
        magnetic_config: Magnetic configuration ('ferro', 'antiferro', or element dict)
        expand_cell: If True, expand cell to accommodate magnetic configuration
        queue: Queue configuration dictionary for job submission (local or remote)
        machine: Name of a machine configuration to load from ~/.xespresso/machines/
        **kwargs: Additional parameters for the calculator
        
    Returns:
        Espresso: Calculator object with results
        
    Example:
        >>> calc = quick_scf(
        ...     'structure.cif',
        ...     {'Si': 'Si.pbe.UPF'},
        ...     protocol='fast'
        ... )
        >>> # With magnetic configuration
        >>> calc = quick_scf(
        ...     atoms,
        ...     {'Fe': 'Fe.pbe-spn.UPF'},
        ...     magnetic_config='antiferro',
        ...     protocol='moderate'
        ... )
        >>> # With remote execution
        >>> calc = quick_scf(
        ...     'structure.cif',
        ...     {'Fe': 'Fe.pbe-spn.UPF'},
        ...     protocol='moderate',
        ...     machine='cluster1'  # Load from ~/.xespresso/machines/cluster1.json
        ... )
    """
    if isinstance(structure, (str, Path)):
        workflow = CalculationWorkflow.from_cif(
            structure, pseudopotentials, protocol, kspacing, 
            magnetic_config=magnetic_config, expand_cell=expand_cell,
            queue=queue, machine=machine, **kwargs
        )
    else:
        workflow = CalculationWorkflow(
            structure, pseudopotentials, protocol, kspacing,
            magnetic_config=magnetic_config, expand_cell=expand_cell,
            queue=queue, machine=machine, **kwargs
        )
    
    return workflow.run_scf(label=label)


def quick_relax(
    structure: Union[str, Path, Atoms],
    pseudopotentials: Dict[str, str],
    label: str = 'relax',
    protocol: str = 'moderate',
    kspacing: Optional[float] = None,
    relax_type: str = 'relax',
    magnetic_config: Optional[Union[str, Dict]] = None,
    expand_cell: bool = False,
    queue: Optional[Dict] = None,
    machine: Optional[str] = None,
    **kwargs
) -> Espresso:
    """
    Quick structure relaxation helper function.
    
    Args:
        structure: CIF file path or ASE Atoms object
        pseudopotentials: Dictionary mapping element symbols to pseudopotential files
        label: Directory/label for the calculation
        protocol: Protocol preset: 'fast', 'moderate', or 'accurate'
        kspacing: K-point spacing in Angstrom^-1 (physical units)
        relax_type: Type of relaxation: 'relax' or 'vc-relax'
        magnetic_config: Magnetic configuration ('ferro', 'antiferro', or element dict)
        expand_cell: If True, expand cell to accommodate magnetic configuration
        queue: Queue configuration dictionary for job submission (local or remote)
        machine: Name of a machine configuration to load from ~/.xespresso/machines/
        **kwargs: Additional parameters for the calculator
        
    Returns:
        Espresso: Calculator object with results
        
    Example:
        >>> calc = quick_relax(
        ...     'structure.cif',
        ...     {'Si': 'Si.pbe.UPF'},
        ...     protocol='moderate',
        ...     relax_type='vc-relax'
        ... )
        >>> # With Hubbard parameters
        >>> calc = quick_relax(
        ...     atoms,
        ...     {'Fe': 'Fe.pbe-spn.UPF', 'O': 'O.pbe.UPF'},
        ...     magnetic_config={'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}},
        ...     protocol='accurate'
        ... )
        >>> # With remote execution on SLURM cluster
        >>> calc = quick_relax(
        ...     'structure.cif',
        ...     {'Fe': 'Fe.pbe-spn.UPF'},
        ...     protocol='moderate',
        ...     machine='slurm_cluster'  # Load from config
        ... )
    """
    if isinstance(structure, (str, Path)):
        workflow = CalculationWorkflow.from_cif(
            structure, pseudopotentials, protocol, kspacing,
            magnetic_config=magnetic_config, expand_cell=expand_cell,
            queue=queue, machine=machine, **kwargs
        )
    else:
        workflow = CalculationWorkflow(
            structure, pseudopotentials, protocol, kspacing,
            magnetic_config=magnetic_config, expand_cell=expand_cell,
            queue=queue, machine=machine, **kwargs
        )
    
    return workflow.run_relax(label=label, relax_type=relax_type)
