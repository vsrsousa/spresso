from ase.calculators.calculator import (
    FileIOCalculator,
    CalculationFailed,
    equal,
    compare_atoms,
    PropertyNotPresent,
)
from xespresso.xio import read_espresso_asei
from xespresso.xespresso import Espresso
import os
import logging
import copy

logger = logging.getLogger(__name__)


class EspressoNscf(Espresso):
    """
    Non-Self-Consistent Field (NSCF) calculation using Quantum ESPRESSO.
    
    Inherits from Espresso to leverage scheduler, execute(), and other infrastructure.
    Automatically reads SCF parameters and structures NSCF as a subfolder of the parent.
    
    Key behavior:
    - Creates NSCF as subfolder of parent: parent_dir/label/
    - Uses SAME prefix as parent calculation
    - Reads density from parent via outdir="../"
    - Inherits all Espresso features (scheduler, remote execution, etc)
    
    Backward compatibility:
    - Old signature: EspressoNscf(scf_directory='scf', prefix='si', label='nscf', ...)
    - New signature: Espresso(...) but with automatic parent reading
    """

    package = "pw"

    def __init__(
        self,
        label=None,
        scf_directory=None,
        prefix=None,
        atoms=None,
        parallel="",
        queue=None,
        debug=False,
        kpts=(10, 10, 10),
        **kwargs
    ):
        """
        Initialize NSCF calculator.
        
        Args:
            label: Directory where NSCF will be created. 
                   If None and scf_directory provided, will be scf_directory/nscf/
            scf_directory: Directory where parent SCF was calculated (e.g., 'scf', '01_scf')
                          Used for backward compatibility.
            prefix: Prefix of parent calculation (same as SCF).
                   If None, will try to detect from scf_directory
            atoms: Atoms object (optional, will be loaded from .asei if not provided)
            parallel: Parallelization options (e.g., '-npools 4')
            queue: Job submission config for remote execution
            debug: Debug logging level
            kpts: K-point mesh tuple
            **kwargs: Additional Espresso parameters
        """
        print("{0:=^60}".format("nscf"))
        
        # Handle backward compatibility with old signature:
        # EspressoNscf(scf_directory='scf', prefix='si')
        # becomes: Espresso(label='scf/nscf', prefix='si')
        if scf_directory is not None and label is None:
            # Old signature detected
            label = os.path.join(scf_directory, "nscf")
            logger.info(f"Backward compatibility: EspressoNscf(scf_directory='{scf_directory}', prefix='{prefix}')")
            logger.info(f"  → Espresso(label='{label}', prefix='{prefix}')")
        
        if label is None:
            label = "nscf"
        
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
            package=self.package,
            parallel=parallel,
            queue=queue,
            debug=debug,
            **kwargs
        )
        
        # Override kpts if provided
        self.parameters['kpts'] = kpts
        
        # NSCF-specific modifications to parameters
        self._configure_nscf_parameters()
    
    def load_scf(self, scf_directory, prefix):
        """
        Load SCF calculation parameters and results.
        
        Args:
            scf_directory: Directory where SCF calculation was done
            prefix: Prefix of SCF calculation (same prefix used for NSCF)
        """
        asei_file = os.path.join(scf_directory, f"{prefix}.asei")
        
        if not os.path.exists(asei_file):
            raise FileNotFoundError(
                f"SCF results not found: {asei_file}\n"
                f"Please run SCF first: workflow.run_scf(label='{scf_directory}', prefix='{prefix}')"
            )
        
        self.atoms, scf_parameters = read_espresso_asei(asei_file, "PW")
        self.parameters = copy.deepcopy(scf_parameters)
        self.scf_parameters = scf_parameters
        logger.info(f"Loaded SCF parameters from: {asei_file}")
    
    def _configure_nscf_parameters(self):
        """
        Configure parameters for NSCF calculation.
        Modifies SCF parameters to be suitable for NSCF.
        """
        # Ensure input_data structure exists
        if 'input_data' not in self.parameters:
            self.parameters['input_data'] = {}
        
        # Ensure CONTROL section exists
        if 'CONTROL' not in self.parameters['input_data']:
            self.parameters['input_data']['CONTROL'] = {}
        
        # Ensure SYSTEM section exists
        if 'SYSTEM' not in self.parameters['input_data']:
            self.parameters['input_data']['SYSTEM'] = {}
        
        # Set NSCF-specific parameters
        self.parameters['input_data']['CONTROL']['calculation'] = 'nscf'
        self.parameters['input_data']['CONTROL']['verbosity'] = 'high'
        self.parameters['input_data']['CONTROL']['outdir'] = '../'  # Read from parent
        
        # NSCF uses tetrahedra integration
        self.parameters['input_data']['SYSTEM']['occupations'] = 'tetrahedra'
        
        # Remove smearing parameters (not used with tetrahedra)
        self.parameters['input_data']['SYSTEM'].pop('degauss', None)
        self.parameters['input_data']['SYSTEM'].pop('smearing', None)
        
        logger.debug("Configured NSCF parameters")
    
    def set_label(self, label, prefix):
        """
        Set directory and prefix from label.
        Override Espresso.set_label to ensure NSCF structure.
        """
        # Call parent's set_label
        Espresso.set_label(self, label, prefix)
        
        # NSCF-specific file naming
        # Note: keep .asei filename distinct from SCF to avoid confusion
        self.asei = os.path.join(self.directory, f"{self.prefix}.nscf_asei")
        
        logger.debug(f"NSCF directory: {self.directory}")
        logger.debug(f"NSCF prefix: {self.prefix}")
        logger.debug(f"NSCF label: {self.label}")

    def write_input(self, atoms, properties=None, system_changes=None):
        """
        Write NSCF input files.
        Override Espresso.write_input to handle NSCF-specific scheduling.
        """
        from xespresso.xio import write_espresso_asei, write_espresso_in
        from xespresso.scheduler import set_queue

        # Call parent to handle FileIOCalculator boilerplate
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        
        # Set up queue configuration if provided
        # Handle case where profile may not exist (dry_run scenarios)
        if self.queue:
            try:
                set_queue(self, package=self.package, parallel=self.parallel, queue=self.queue)
            except AttributeError as e:
                logger.debug(f"Could not set queue configuration: {e}")
        
        # Write QE input file
        write_espresso_in(self.label + ".pwi", atoms, **self.parameters)
        logger.debug(f"Wrote NSCF input: {self.label}.pwi")
        
        # Write ASE state file
        write_espresso_asei(self.asei, self.state_info, self.parameters)
        logger.debug(f"Wrote NSCF state: {self.asei}")
        logger.debug(f"Parameters: {self.parameters}")

    @property
    def state_info(self):
        """
        Get state info from parent SCF charge-density file.
        Used to detect if recalculation is needed.
        
        Returns dummy hash if file doesn't exist (dry_run scenario).
        """
        from xespresso.utils import get_hash

        # Parent directory is one level up from NSCF directory
        parent_dir = os.path.dirname(self.directory)
        
        # Try charge-density.hdf5 first (faster), then .dat format
        charge_hdf5 = os.path.join(parent_dir, f"{self.prefix}.save/charge-density.hdf5")
        charge_dat = os.path.join(parent_dir, f"{self.prefix}.save/charge-density.dat")
        
        if os.path.isfile(charge_hdf5):
            state_info = get_hash(charge_hdf5)
            logger.debug(f"Using charge-density.hdf5 for state info")
        elif os.path.isfile(charge_dat):
            state_info = get_hash(charge_dat)
            logger.debug(f"Using charge-density.dat for state info")
        else:
            # During dry_run or before SCF runs, file doesn't exist yet
            # Return a constant hash so we can still write input files
            logger.debug(f"charge-density file not found (dry_run or pending SCF)")
            state_info = "dry_run_state"
        
        return state_info

    def check_state(self):
        """
        Check if NSCF needs to be recalculated.
        Compares current state with previous NSCF output.
        """
        # Read state information from parent charge-density file
        self.state_parameters = self.parameters
        
        # Check if NSCF output exists
        output, message = self.read_convergence_post("pw")
        logger.debug(f"Check state: {message}")
        
        if output:
            # NSCF output exists, check if parameters changed
            if os.path.isfile(self.asei):
                system_changes = self.check_state_post(self.asei, package="PW")
                if not system_changes:
                    logger.debug("Using previous NSCF results (no changes)")
                    return False
            else:
                logger.debug("No NSCF state file. Recalculating...")
        else:
            logger.debug("No NSCF output file. Starting new calculation...")
        
        return True

    def check_state_post(self, asei, package):
        """
        Check if NSCF state file matches current parameters.
        """
        logger.debug(f"Checking state file: {asei}")
        
        if not os.path.exists(asei):
            return True  # No state file, need to recalculate
        
        try:
            old_state_info, old_state_parameters = read_espresso_asei(asei, package)
        except Exception as e:
            logger.debug(f"Could not read old state: {e}")
            return True  # If we can't read old state, recalculate
        
        # Check if parent charge-density changed
        if self.state_info != old_state_info:
            logger.debug("Parent charge-density changed")
            return True
        
        # Check if parameters changed
        if self.state_parameters != old_state_parameters:
            logger.debug("NSCF parameters changed")
            return True
        
        return False

    def read_convergence_post(self, package="pw"):
        """
        Check if NSCF calculation completed successfully.
        """
        output_file = f"{self.label}.{package}o"
        logger.debug(f"Reading output: {output_file}")
        
        if not os.path.exists(output_file):
            return False, "Output file not found"
        
        try:
            with open(output_file, "r") as f:
                lines = f.readlines()
                if not lines:
                    return False, "Output file is empty"
                
                # Check last 100 lines for JOB DONE message
                nlines = len(lines)
                n = min(100, nlines)
                for line in lines[-n:]:
                    if "JOB DONE" in line:
                        return True, line.strip()
                
                return False, "JOB DONE not found in output"
        except Exception as e:
            return False, f"Error reading output: {e}"

    def run(self):
        """
        Run NSCF calculation.
        Legacy method for backward compatibility.
        """
        if self.check_state():
            self.calculate()
        else:
            logger.info("Skipping NSCF calculation (no changes)")

