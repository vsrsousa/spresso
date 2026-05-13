"""
Band structure calculator (post-processor for SCF density).

Similar to NSCF, BANDS reads the charge density from a parent SCF calculation
and computes the band structure along high-symmetry k-points.

Structure:
    parent_dir/
      {prefix}.save/    ← charge density (from SCF)
      bands/            ← BANDS subfolder
        {prefix}.pwi    ← BANDS input with outdir="../"
"""

import os
import logging
from pathlib import Path
from xespresso.xespresso import Espresso
from xespresso.xio import read_espresso_asei

logger = logging.getLogger(__name__)


class EspressoBands(Espresso):
    """
    Band structure calculator using parent SCF density.
    
    Inherits from Espresso to leverage scheduler integration and execution.
    Handles:
    - Reading SCF parameters from parent calculation
    - Creating bands subfolder automatically
    - Setting outdir="../" to read parent density
    - Same prefix consistency
    
    Parameters:
        label: Directory where BANDS will be created (e.g., 'scf/bands')
        scf_directory: Directory where parent SCF was calculated (e.g., 'scf')
        prefix: Prefix of parent calculation (same as SCF)
        **kwargs: Additional Espresso parameters (kpts, parallel, queue, etc.)
    
    Example:
        >>> bands = EspressoBands(
        ...     scf_directory='scf',
        ...     prefix='si',
        ...     kpts=bandpath  # seekpath object with high-symmetry points
        ... )
        >>> bands.run(atoms)
    """
    
    def __init__(
        self,
        label=None,
        scf_directory=None,
        prefix=None,
        atoms=None,
        parallel="",
        queue=None,
        debug=False,
        **kwargs
    ):
        """
        Initialize BANDS calculator.
        
        Args:
            label: Directory where BANDS will be created.
                   If None and scf_directory provided, will be scf_directory/bands/
            scf_directory: Directory where parent SCF was calculated (e.g., 'scf')
            prefix: Prefix of parent calculation (same as SCF)
            atoms: Atoms object (optional, will be loaded from .asei if not provided)
            parallel: Parallelization options
            queue: Job submission config for remote execution
            debug: Debug logging level
            **kwargs: Additional Espresso parameters
        """
        print("{0:=^60}".format("bands"))
        
        # Handle backward compatibility with old signature:
        # EspressoBands(scf_directory='scf', prefix='si')
        # becomes: Espresso(label='scf/bands', prefix='si')
        if scf_directory is not None and label is None:
            # Old signature detected
            label = os.path.join(scf_directory, "bands")
            logger.info(f"Backward compatibility: EspressoBands(scf_directory='{scf_directory}', prefix='{prefix}')")
            logger.info(f"  → Espresso(label='{label}', prefix='{prefix}')")
        
        if label is None:
            label = "bands"
        
        # Store parent directory for load_scf()
        self.scf_directory = scf_directory
        
        # Load SCF parameters BEFORE calling Espresso.__init__
        # This ensures atoms and parameters are available for Espresso
        if scf_directory and prefix:
            self.load_scf(scf_directory, prefix)
            atoms = self.atoms  # Use loaded atoms
            input_data = kwargs.get('input_data', self.parameters.get('input_data'))
        else:
            # If no scf_directory, proceed like normal Espresso
            # (may fail if atoms not provided or calculations not available)
            input_data = kwargs.get('input_data', {})
        
        # Call parent Espresso.__init__
        # This will handle label, prefix, atoms, etc.
        Espresso.__init__(
            self,
            label=label,
            prefix=prefix,
            atoms=atoms,
            parallel=parallel,
            queue=queue,
            debug=debug,
            **kwargs
        )
        
        # Set package for this calculator
        self.package = "pw"
        
        # BANDS-specific modifications to parameters
        self._configure_bands_parameters()
    
    def load_scf(self, scf_directory, prefix):
        """
        Load SCF calculation parameters and results.
        
        Args:
            scf_directory: Directory where SCF calculation was done
            prefix: Prefix of SCF calculation (same prefix used for BANDS)
        """
        scf_path = Path(scf_directory)
        asei_file = scf_path / f"{prefix}.asei"
        
        if not asei_file.exists():
            raise FileNotFoundError(
                f"SCF results not found: {asei_file}\n"
                f"Please run SCF first: workflow.run_scf(label='{scf_directory}', prefix='{prefix}')"
            )
        
        # Load SCF state
        self.atoms, self.parameters = read_espresso_asei(str(asei_file))
        logger.info(f"Loaded SCF parameters from: {asei_file}")
    
    def _configure_bands_parameters(self):
        """Configure parameters for band structure calculation."""
        # Set calculation to 'bands'
        if 'CONTROL' not in self.parameters['input_data']:
            self.parameters['input_data']['CONTROL'] = {}
        self.parameters['input_data']['CONTROL']['calculation'] = 'bands'
        
        # Read density from parent: outdir="../"
        # This points to scf_directory (one level up from bands/)
        self.parameters['input_data']['CONTROL']['outdir'] = '../'
        
        logger.info(f"BANDS will read density from: ../")
