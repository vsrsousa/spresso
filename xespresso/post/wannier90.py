from xespresso.post.base import PostCalculation
import os
from pathlib import Path
from types import SimpleNamespace


class EspressoWannier90(PostCalculation):
    """
    Interface to wannier90.x for generating Wannier functions from QE outputs.
    
    This class handles:
    - Generation of wannier_seed.win input files  
    - Execution of wannier90.x for Wannier function generation
    - Collection of output files (wout, centres, xsf, etc.)
    
    The workflow is typically:
    1. Use EspressoPw2wannier90 to convert QE wavefunctions
    2. Use EspressoWannier90 to generate Wannier functions
    3. Post-process results (band interpolation, visualization, etc.)
    
    Example:
        >>> w90 = EspressoWannier90(
        ...     parent_directory='/path/to/pw2wannier_output',
        ...     seedname='wannier_seed',
        ...     num_wann=16,
        ...     projections='Fe: d',
        ...     directory='./runs/05-wan'
        ... )
        >>> result = w90.run()
    """
    
    package = "wannier90"
    package_parameters = {}  # Parameters go into .win file, not Namelist sections
    
    def __init__(
        self,
        parent_directory,
        seedname='wannier',
        num_wann=None,
        projections=None,
        spinors=False,
        dis_num_iter=1000,
        spin_component='none',
        queue=None,
        parallel='',
        debug=False,
        dry_run=False,
        directory=None,
        **kwargs,
    ):
        """
        Initialize EspressoWannier90.
        
        Parameters:
        - parent_directory: Directory where pw2wannier90 outputs are (amn, mmn, eig files)
        - seedname: Base name for Wannier seed (e.g., 'wannier_seed')
        - num_wann: Number of Wannier functions to generate
        - projections: Initial projection string (multi-line suitable for begin projections block)
        - spinors: If True, include spinor calculation (spinors = .true.)
        - dis_num_iter: Number of disentanglement iterations (default 1000)
        - spin_component: Spin component ('up', 'down', 'none'). For spin-polarized calculations,
                         seedname gets _up or _dn suffix automatically.
        - queue: Scheduler configuration dict
        - parallel: Parallelization flags for wannier90.x
        - debug: Enable debug logging
        - dry_run: If True, only generate input without execution
        - directory: Directory for wannier_seed.win file (uses parent_directory if None)
        - **kwargs: Additional Wannier90 parameters
        """
        self.seedname = seedname
        self.num_wann = num_wann
        self.projections = projections or 'auto'
        self.spinors = spinors
        self.dis_num_iter = dis_num_iter
        self.spin_component = spin_component
        self.parallel = parallel
        self.dry_run = dry_run
        
        # Apply spin suffix to seedname if spin-polarized
        if spin_component == 'up':
            seedname_full = f"{seedname}_up"
        elif spin_component == 'down':
            seedname_full = f"{seedname}_dn"
        else:
            seedname_full = seedname
        
        # Prefix for wannier90 is the full seedname (with spin suffix if applicable)
        super().__init__(
            parent_directory=parent_directory,
            prefix=seedname_full,
            queue=queue,
            debug=debug,
            directory=directory,
            **kwargs,
        )
        
        # wannier90 doesn't use Namelist style input, so parameters dict is not used
        # Configuration goes into the .win file instead
    
    def _generate_win_content(self) -> str:
        """Generate the content of the .win (input) file for wannier90."""
        import textwrap
        
        lines = []
        
        # Basic block
        lines.append(f"num_wann = {self.num_wann if self.num_wann else 4}")
        
        # Projections block
        if self.projections and self.projections != 'auto':
            lines.append("")
            lines.append("begin projections")
            if isinstance(self.projections, str):
                # Handle multiline projections
                for proj_line in self.projections.strip().split('\n'):
                    lines.append(f"  {proj_line}")
            lines.append("end projections")
        
        # Spinors
        if self.spinors:
            lines.append("")
            lines.append("spinors = .true.")
        
        # Disentanglement
        if self.dis_num_iter:
            lines.append("")
            lines.append(f"dis_num_iter = {self.dis_num_iter}")
            lines.append("dis_froz_max = 8.0")
        
        # Default settings
        lines.append("")
        lines.append("# Default bands and parameters")
        lines.append("exclude_bands = ")
        lines.append("iprint = 2")
        lines.append("restart = default")
        lines.append("wvfn_formatted = .false.")
        lines.append("dis_mix_ratio = 0.5")
        lines.append("kmesh_tol = 0.0000001")
        
        return '\n'.join(lines) + '\n'
    
    def write_input(self):
        """Generate the .win file for wannier90."""
        win_content = self._generate_win_content()
        win_file = os.path.join(self.directory, f"{self.prefix}.win")
        with open(win_file, 'w') as f:
            f.write(win_content)
        return win_file
    
    def run(self, dry_run: bool = None, blocking: bool = True, pp: bool = True, **kwargs):
        """
        Run wannier90 with optional preprocessing, following PostCalculation pattern.
        
        Parameters:
        - dry_run: If True, only generate .win file without execution (default: self.dry_run)
        - blocking: If True, wait for job to complete (default: True)
        - pp: If True, run wannier90.x -pp seedname preprocessing first (default: True)
        - **kwargs: Additional parameters (ignored, for compatibility)
        
        Returns:
        - Dict with keys: status, run_dir, outputs, job_id, message
        """
        from xespresso.schedulers.factory import get_scheduler
        import logging
        
        logger = logging.getLogger(__name__)
        
        if dry_run is None:
            dry_run = self.dry_run
        
        print("{0:=^60}".format(self.package))
        
        # Generate .win input file
        self.write_input()
        logger.info(f"{self.package} input files generated in {self.directory}")
        
        # Dry run: return early after generating input
        if dry_run:
            win_file = os.path.join(self.directory, f"{self.prefix}.win")
            return {
                "status": "input_files_generated",
                "run_dir": self.directory,
                "message": f"{self.package} input files generated (dry_run mode)",
                "input_file": win_file if os.path.exists(win_file) else None,
                "outputs": {},
            }
        
        # Non-dry-run: execute preprocessing and main run
        results = {"outputs": {}, "run_dir": self.directory, "job_id": None}
        
        try:
            queue = self.queue if self.queue else {}
            
            # Step 1: Preprocessing (wannier90.x -pp seedname)
            if pp:
                cmd_pp = f"wannier90.x -pp {self.prefix}"
                calc_stub = SimpleNamespace(
                    directory=self.directory,
                    prefix=self.prefix,
                    queue=queue,
                )
                scheduler_pp = get_scheduler(calc_stub, queue, cmd_pp)
                scheduler_pp.write_script()
                scheduler_pp.run()
            
            # Step 2: Main wannier90 run (wannier90.x seedname)
            cmd = f"wannier90.x {self.prefix}"
            calc_stub = SimpleNamespace(
                directory=self.directory,
                prefix=self.prefix,
                queue=queue,
            )
            self.scheduler = get_scheduler(calc_stub, queue, cmd)
            self.scheduler.write_script()
            self.scheduler.run()
            
            results["job_id"] = getattr(calc_stub, 'last_job_id', None)
            
            # Step 3: Collect outputs
            wout_file = os.path.join(self.directory, f"{self.prefix}.wout")
            if os.path.exists(wout_file):
                results["outputs"]["wout"] = wout_file
                # Check for other output files
                for ext in ['centres', 'spreads', 'xsf']:
                    out_file = os.path.join(self.directory, f"{self.prefix}_{ext}.dat")
                    if os.path.exists(out_file):
                        results["outputs"][ext] = out_file
            
            results["status"] = "finished" if blocking else "submitted"
            results["message"] = f"{self.package} calculation completed"
            
            return results
            
        except Exception as e:
            logger.error(f"{self.package} execution failed: {e}")
            return {
                "status": "error",
                "run_dir": self.directory,
                "message": str(e),
                "outputs": {},
            }
        finally:
            print("Done: %s" % self.package)
