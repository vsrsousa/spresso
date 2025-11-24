"""
Calculation preparation module for xespresso GUI.

This module handles the creation of Espresso calculator and atoms objects
from GUI configuration, following xespresso's design patterns.

This module properly uses xespresso's setup_magnetic_config() function
to handle magnetic and Hubbard configurations correctly.
"""

from typing import Dict, Tuple, Optional
from ase import Atoms
from xespresso import Espresso
from xespresso.tools import setup_magnetic_config
from xespresso.gui.calculations.base import BaseCalculationPreparation
import logging

logger = logging.getLogger(__name__)


class CalculationPreparation(BaseCalculationPreparation):
    """
    Prepares calculations (SCF, relax, vc-relax) from GUI configuration.
    
    This class encapsulates the logic for creating Espresso calculator
    objects from GUI workflow_config dictionary, ensuring that calculator
    creation happens in a dedicated module rather than in main GUI code
    or job submission code.
    
    Examples:
        >>> # In GUI calculation setup or workflow builder:
        >>> prep = CalculationPreparation(atoms, config, label='scf/fe')
        >>> atoms, calc = prep.prepare()
        >>> 
        >>> # In GUI job submission:
        >>> # Just receive the prepared objects and execute
        >>> calc.run(atoms=atoms)
    """
    
    def prepare(self) -> Tuple[Atoms, Espresso]:
        """
        Prepare atoms and Espresso calculator from GUI configuration.
        
        This method creates the Espresso calculator with parameters from
        the GUI's workflow_config dictionary, following xespresso patterns.
        
        Uses xespresso's setup_magnetic_config() to properly handle magnetic
        and Hubbard configurations.
        
        Returns:
            tuple: (atoms, calculator) ready for execution or dry run
        """
        config = self.config
        atoms = self.atoms.copy()  # Work with a copy
        
        # Validate required configuration
        if 'pseudopotentials' not in config or not config['pseudopotentials']:
            raise ValueError("Configuration must include pseudopotentials")
        
        # Check if magnetic and/or Hubbard configuration is enabled
        enable_magnetism = config.get('enable_magnetism', False)
        enable_hubbard = config.get('enable_hubbard', False)
        
        # Prepare pseudopotentials and atoms with magnetic/Hubbard config if needed
        if enable_magnetism or enable_hubbard:
            # Build magnetic_config dict for setup_magnetic_config()
            magnetic_config = {}
            
            # Get unique elements in the structure
            elements = set(atoms.get_chemical_symbols())
            
            for element in elements:
                element_config = {}
                
                # Add magnetic moments if magnetism is enabled
                if enable_magnetism and 'magnetic_config' in config:
                    if element in config['magnetic_config']:
                        mag_moments = config['magnetic_config'][element]
                        # Ensure it's a list
                        if not isinstance(mag_moments, list):
                            mag_moments = [mag_moments]
                        element_config['mag'] = mag_moments
                    else:
                        # Default to non-magnetic
                        element_config['mag'] = [0]
                else:
                    # If magnetism not enabled, use [0] for all
                    element_config['mag'] = [0]
                
                # Add Hubbard U if enabled
                if enable_hubbard and 'hubbard_u' in config:
                    if element in config['hubbard_u']:
                        u_value = config['hubbard_u'][element]
                        
                        # Check if using new format with orbital specification
                        hubbard_format = config.get('hubbard_format', 'old')
                        if hubbard_format == 'new':
                            # New format requires orbital specification
                            orbital = config.get(f'hubbard_orbital_{element}', '3d')
                            element_config['U'] = {orbital: u_value}
                        else:
                            # Old format - just the value
                            element_config['U'] = u_value
                
                # Only add element if it has configuration
                if 'mag' in element_config or 'U' in element_config:
                    magnetic_config[element] = element_config
            
            # Call setup_magnetic_config if we have any configuration
            if magnetic_config:
                logger.info(f"Setting up magnetic/Hubbard configuration: {magnetic_config}")
                
                # Get QE version and other parameters
                qe_version = config.get('qe_version', None)
                if qe_version == 'auto':
                    qe_version = None
                
                hubbard_format = config.get('hubbard_format', 'auto')
                projector = config.get('hubbard_projector', 'ortho-atomic')
                expand_cell = config.get('expand_cell', False)
                
                # Call xespresso's setup_magnetic_config
                mag_result = setup_magnetic_config(
                    atoms,
                    magnetic_config,
                    pseudopotentials=config['pseudopotentials'],
                    expand_cell=expand_cell,
                    qe_version=qe_version,
                    hubbard_format=hubbard_format,
                    projector=projector
                )
                
                # Extract results from setup_magnetic_config
                atoms = mag_result['atoms']  # May be modified/expanded
                pseudopotentials = mag_result['pseudopotentials']  # Updated with species
                input_ntyp = mag_result.get('input_ntyp', {})  # starting_magnetization, Hubbard_U
                
                logger.info(f"Magnetic/Hubbard setup complete. Species: {list(pseudopotentials.keys())}")
            else:
                # No magnetic/Hubbard config
                pseudopotentials = config['pseudopotentials']
                input_ntyp = {}
        else:
            # Neither magnetism nor Hubbard enabled
            pseudopotentials = config['pseudopotentials']
            input_ntyp = {}
        
        # Build calculator parameters from GUI configuration
        calc_params = {
            'pseudopotentials': pseudopotentials,
            'label': self.label,
        }
        
        # Build input_data dictionary following xespresso patterns
        input_data = {}
        
        # Add basic parameters
        if 'ecutwfc' in config:
            input_data['ecutwfc'] = config['ecutwfc']
        if 'ecutrho' in config:
            input_data['ecutrho'] = config['ecutrho']
        if 'occupations' in config:
            input_data['occupations'] = config['occupations']
        if 'conv_thr' in config:
            input_data['conv_thr'] = config['conv_thr']
        
        # Add smearing parameters if using smearing occupation
        if config.get('occupations') == 'smearing':
            input_data['smearing'] = config.get('smearing', 'gaussian')
            input_data['degauss'] = config.get('degauss', 0.02)
        
        # Add spin polarization if magnetism is enabled
        if enable_magnetism:
            input_data['nspin'] = 2
        elif 'nspin' in config:
            input_data['nspin'] = config['nspin']
        
        # Add magnetic configuration from input_ntyp (set by setup_magnetic_config)
        if input_ntyp:
            input_data['input_ntyp'] = input_ntyp
        
        # Add Hubbard parameters if using new format
        if enable_hubbard and enable_magnetism:
            # Check if new format was used
            if 'hubbard' in mag_result:
                input_data['hubbard'] = mag_result['hubbard']
                if 'hubbard_format' in mag_result:
                    input_data['hubbard_format'] = mag_result['hubbard_format']
                if 'qe_version' in mag_result:
                    input_data['qe_version'] = mag_result['qe_version']
        
        # Add lda_plus_u flag if Hubbard is enabled
        if enable_hubbard:
            input_data['lda_plus_u'] = True
        
        # Add calculation type
        calc_type = config.get('calc_type', 'scf')
        if calc_type in ['relax', 'vc-relax']:
            input_data['calculation'] = calc_type
            # Add relaxation-specific parameters
            if 'forc_conv_thr' in config:
                input_data['forc_conv_thr'] = config['forc_conv_thr']
        else:
            input_data['calculation'] = 'scf'
        
        calc_params['input_data'] = input_data
        
        # Add k-points from config
        # Note: kspacing is converted to kpts in the GUI, so we only handle kpts here
        if 'kpts' in config:
            calc_params['kpts'] = config['kpts']
        else:
            # Default to gamma point
            calc_params['kpts'] = (1, 1, 1)
        
        # Add queue configuration if present (for job submission)
        if 'queue' in config and config['queue']:
            calc_params['queue'] = config['queue']
        
        # Add parallel execution parameters if present
        if 'parallel' in config and config['parallel']:
            calc_params['parallel'] = config['parallel']
        
        # Add any additional parameters from config
        for key in ['package', 'debug']:
            if key in config:
                calc_params[key] = config[key]
        
        # Create Espresso calculator using xespresso
        logger.info(f"Creating Espresso calculator with label={self.label}")
        self.calculator = Espresso(**calc_params)
        
        # Update the atoms object we're returning
        self.atoms = atoms
        
        logger.info(f"Successfully prepared {calc_type} calculation")
        
        return self.atoms, self.calculator


def prepare_calculation_from_gui(
    atoms: Atoms,
    config: Dict,
    label: str = "calculation"
) -> Tuple[Atoms, Espresso]:
    """
    Convenience function to prepare calculation from GUI configuration.
    
    This function encapsulates the calculator preparation logic,
    ensuring that atoms and Espresso objects are created in the
    calculation module rather than in main GUI or job submission code.
    
    Args:
        atoms: ASE Atoms object
        config: GUI workflow configuration dictionary
        label: Calculation label/directory
        
    Returns:
        tuple: (atoms, calculator) ready for execution
        
    Example:
        >>> # In GUI pages (calculation_setup, workflow_builder):
        >>> atoms, calc = prepare_calculation_from_gui(
        ...     st.session_state.current_structure,
        ...     st.session_state.workflow_config,
        ...     label='scf/fe'
        ... )
        >>> st.session_state.espresso_calculator = calc
        >>> st.session_state.prepared_atoms = atoms
    """
    prep = CalculationPreparation(atoms, config, label)
    return prep.prepare()


def dry_run_calculation(
    atoms: Atoms,
    config: Dict,
    label: str = "calculation"
) -> Tuple[Atoms, Espresso]:
    """
    Prepare calculation and generate input files (dry run).
    
    This function creates the calculator and writes input files
    without executing the calculation. Useful for reviewing
    parameters before submission.
    
    Args:
        atoms: ASE Atoms object
        config: GUI workflow configuration dictionary
        label: Calculation label/directory
        
    Returns:
        tuple: (atoms, calculator) with input files written
        
    Example:
        >>> # In GUI dry run tab:
        >>> atoms, calc = dry_run_calculation(
        ...     st.session_state.current_structure,
        ...     st.session_state.workflow_config,
        ...     label='scf/fe'
        ... )
        >>> st.success(f"Input files written to {calc.directory}")
    """
    atoms, calc = prepare_calculation_from_gui(atoms, config, label)
    
    # Write input files using xespresso's method
    calc.write_input(atoms)
    
    logger.info(f"Dry run complete - input files written to {calc.directory}")
    
    return atoms, calc
