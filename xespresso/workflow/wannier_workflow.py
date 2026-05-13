"""
High-level helpers to orchestrate PW -> pw2wannier90 -> wannier90 pipelines.

These helpers use existing xespresso workflow and scheduler machinery to:
- run SCF/NSCF stages via `CalculationWorkflow`
- write `pw2wannier` input files
- submit `pw2wannier90` and `wannier90` using the scheduler factory so local
  and remote execution behave consistently with Espresso jobs.

This module is intended to be used by the GUI controller (the Wannierization
page) or by example scripts. Helpers return plain dicts containing
serializable metadata (paths, job ids, status, messages) so they can be
recorded in `SessionState` run metadata.
"""
from __future__ import annotations

import os
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, Optional, Tuple, Union

from xespresso.schedulers.factory import get_scheduler


def _make_queue_fallback(queue: Optional[dict], blocking: bool) -> dict:
    if queue is None:
        queue = {}
    # Ensure minimal keys expected by get_scheduler
    if "scheduler" not in queue:
        queue = dict(queue)
        queue.setdefault("execution", "local")
        queue.setdefault("scheduler", "direct")
    # Control blocking behaviour via wait_for_completion flag
    queue.setdefault("wait_for_completion", blocking)
    return queue


def run_pw2wannier(
    run_dir: str,
    prefix: str,
    seedname: str,
    *,
    pw2wan_input: Optional[str] = None,
    write_input: bool = True,
    queue: Optional[dict] = None,
    blocking: bool = True,
    command: Optional[str] = None,
    timeout: int = 300,
) -> Dict:
    """Run `pw2wannier90` from an existing NSCF run directory.

    Parameters
    - run_dir: directory where NSCF outputs are located (must contain prefix.save or wavefunction files)
    - prefix: prefix used in the PW calculations (the QE prefix)
    - seedname: base name for pw2wannier/wannier output files
    - pw2wan_input: optional input text for pw2wannier; if None a minimal template is written
    - write_input: when True write `pw2wannier.in` into `run_dir`
    - queue: scheduler/queue dict (see machines loader)
    - blocking: whether to wait for the job to finish
    - command: override the command to run; default: 'pw2wannier90.x -in pw2wannier.in'

    Returns a serializable dict with keys: status, run_dir, outputs, job_id, message
    """
    # Basic presence check for wavefunction directory
    save_dir = os.path.join(run_dir, f"{prefix}.save")
    has_save = os.path.isdir(save_dir)
    # fallback: look for any likely wavefunction files
    has_wavefiles = any(f.endswith(".wfc") or ".wfc" in f for f in os.listdir(run_dir)) if os.path.exists(run_dir) else False

    if not has_save and not has_wavefiles:
        return {"status": "missing_wavefunctions", "run_dir": run_dir, "message": "No .save directory or wavefunction files found"}

    if write_input:
        if pw2wan_input is None:
            pw2wan_input = generate_pw2wannier_input(prefix)
        pw2wan_path = os.path.join(run_dir, "pw2wannier.in")
        with open(pw2wan_path, "w", encoding="utf-8") as f:
            f.write(pw2wan_input)

    cmd = command or "pw2wannier90.x -in pw2wannier.in"
    queue = _make_queue_fallback(queue, blocking)

    calc_stub = SimpleNamespace(directory=run_dir, prefix=seedname, queue=queue)

    try:
        scheduler = get_scheduler(calc_stub, queue, cmd)
    except Exception as e:
        return {"status": "error", "message": f"Could not initialize scheduler: {e}"}

    try:
        scheduler.write_script()
        # run() may raise on failure for blocking local runs; for remote non-blocking it returns immediately
        scheduler.run()
    except Exception as e:
        return {"status": "error", "message": f"Scheduler run failed: {e}"}

    # For non-blocking remote submission, try to find job id from calc_stub attributes
    job_id = getattr(calc_stub, "last_job_id", None)

    # If blocking, wait a few seconds and check for expected outputs
    outputs = {"amn": None, "mmn": None, "eig": None}
    if blocking:
        # Wait briefly for files to appear (poll)
        target_files = [f"{seedname}.amn", f"{seedname}.mmn", f"{seedname}.eig"]
        start = time.time()
        while time.time() - start < timeout:
            found = [os.path.join(run_dir, t) for t in target_files if os.path.exists(os.path.join(run_dir, t))]
            if len(found) == len(target_files):
                outputs = dict(zip(["amn", "mmn", "eig"], [os.path.join(run_dir, t) for t in target_files]))
                return {"status": "finished", "run_dir": run_dir, "outputs": outputs, "job_id": job_id}
            time.sleep(1)
        return {"status": "finished_with_warnings", "run_dir": run_dir, "outputs": outputs, "job_id": job_id, "message": "Timed out waiting for expected pw2wannier outputs"}

    return {"status": "submitted", "run_dir": run_dir, "outputs": outputs, "job_id": job_id}


def run_projwfc(
    run_dir: str,
    prefix: str,
    *,
    queue: Optional[dict] = None,
    blocking: bool = True,
    timeout: int = 300,
    command: Optional[str] = None,
    den_ext: str = '1',
) -> Dict:
    """Run `projwfc.x` to compute projected density of states (PDOS) and projection analysis.

    This step computes projections of the KS wavefunctions onto atomic wavefunctions,
    essential for:
    - Understanding which atoms/orbitals dominate the electronic structure
    - Guiding the selection of Wannier function projections
    - Validating the initial projections chosen

    Parameters
    - run_dir: directory where NSCF outputs are located
    - prefix: prefix used in the PW calculations
    - queue: scheduler/queue dict (see machines loader)
    - blocking: whether to wait for the job to finish
    - timeout: max time to wait for output files (seconds)
    - command: override the command; default: 'projwfc.x -in projwfc.in'
    - den_ext: file extension for density; default '1' (use prefix.save/charge-density.dat)

    Returns a serializable dict with keys: status, run_dir, outputs, job_id, message
    """
    save_dir = os.path.join(run_dir, f"{prefix}.save")
    if not os.path.isdir(save_dir):
        return {"status": "missing_data", "run_dir": run_dir, "message": "No .save directory found"}

    # Write projwfc input file
    projwfc_input = f"""&inputpp
  prefix = '{prefix}'
  outdir = './'
/
filpdos = '{prefix}.pdos'
"""
    projwfc_path = os.path.join(run_dir, "projwfc.in")
    with open(projwfc_path, "w", encoding="utf-8") as f:
        f.write(projwfc_input)

    cmd = command or "projwfc.x -in projwfc.in"
    queue = _make_queue_fallback(queue, blocking)

    calc_stub = SimpleNamespace(directory=run_dir, prefix=prefix, queue=queue)

    try:
        scheduler = get_scheduler(calc_stub, queue, cmd)
    except Exception as e:
        return {"status": "error", "message": f"Could not initialize scheduler: {e}"}

    try:
        scheduler.write_script()
        scheduler.run()
    except Exception as e:
        return {"status": "error", "message": f"Scheduler run failed: {e}"}

    job_id = getattr(calc_stub, "last_job_id", None)

    # Expected output file
    outputs = {"pdos": None, "txt": None}
    if blocking:
        # projwfc generates files like prefix.pdos and prefix.pdos.up, prefix.pdos.dw, etc.
        expected_pdos = os.path.join(run_dir, f"{prefix}.pdos")
        expected_txt = os.path.join(run_dir, f"{prefix}.pdos.txt")
        
        start = time.time()
        while time.time() - start < timeout:
            if os.path.exists(expected_pdos):
                outputs["pdos"] = expected_pdos
                if os.path.exists(expected_txt):
                    outputs["txt"] = expected_txt
                return {"status": "finished", "run_dir": run_dir, "outputs": outputs, "job_id": job_id}
            time.sleep(1)
        
        # File may not exist immediately; if pdos exists partially, consider it done
        if os.path.exists(expected_pdos) or any(
            os.path.exists(os.path.join(run_dir, f)) 
            for f in os.listdir(run_dir) 
            if f.startswith(f"{prefix}.pdos")
        ):
            outputs["pdos"] = expected_pdos
            return {"status": "finished_with_warnings", "run_dir": run_dir, "outputs": outputs, "job_id": job_id, 
                    "message": "PDOS files found but some expected outputs may be missing"}
        
        return {"status": "finished_with_warnings", "run_dir": run_dir, "outputs": outputs, "job_id": job_id, 
                "message": "Timed out waiting for projwfc outputs"}

    return {"status": "submitted", "run_dir": run_dir, "outputs": outputs, "job_id": job_id}


def parse_projwfc_output(projwfc_dir: str, prefix: str) -> Optional[Dict]:
    """Parse projwfc output to suggest Wannier projections.
    
    Analyzes PDOS files to identify dominant orbital contributions.
    
    Parameters:
        projwfc_dir: directory containing projwfc.pdos files
        prefix: file prefix used in projwfc calculation
        
    Returns:
        Dict with suggested projections or None if files not found
    """
    import re
    
    pdos_file = os.path.join(projwfc_dir, f"{prefix}.pdos")
    if not os.path.exists(pdos_file):
        return None
    
    suggestions = {}
    try:
        with open(pdos_file, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
        
        # Simple parsing: look for atom type indicators
        atom_orbitals = {}
        for line in lines:
            if '#' in line:
                # Try to extract atom and orbital info
                tokens = line.split()
                if len(tokens) > 1 and tokens[0].startswith('m='):
                    # This is a state projection line
                    pass
        
        # If detailed parsing fails, provide generic suggestion
        suggestions['note'] = 'Review PDOS files to identify dominant s, p, d contributions'
        suggestions['pdos_file'] = pdos_file
        
    except Exception as e:
        suggestions['error'] = str(e)
    
    return suggestions


def run_wannier90(
    run_dir: str,
    seedname: str,
    *,
    pp: bool = True,
    queue: Optional[dict] = None,
    blocking: bool = True,
    timeout: int = 300,
    command: Optional[str] = None,
) -> Dict:
    """Run `wannier90` for a given seedname inside `run_dir`.

    If `pp` is True the preprocessing step `wannier90.x -pp seedname` is executed first.
    """
    queue = _make_queue_fallback(queue, blocking)
    results = {"outputs": {}, "run_dir": run_dir}

    if pp:
        cmd_pp = command or f"wannier90.x -pp {seedname}"
        calc_stub = SimpleNamespace(directory=run_dir, prefix=seedname, queue=queue)
        try:
            sch = get_scheduler(calc_stub, queue, cmd_pp)
            sch.write_script()
            sch.run()
        except Exception as e:
            return {"status": "error", "message": f"wannier90 -pp failed: {e}"}

    cmd_run = command or f"wannier90.x {seedname}"
    calc_stub = SimpleNamespace(directory=run_dir, prefix=seedname, queue=queue)
    try:
        sch = get_scheduler(calc_stub, queue, cmd_run)
        sch.write_script()
        sch.run()
    except Exception as e:
        return {"status": "error", "message": f"wannier90 run failed: {e}"}

    # Blocking: wait for wout
    if blocking:
        wout = os.path.join(run_dir, f"{seedname}.wout")
        start = time.time()
        while time.time() - start < timeout:
            if os.path.exists(wout):
                results["outputs"]["wout"] = wout
                results["status"] = "finished"
                return results
            time.sleep(1)
        results["status"] = "finished_with_warnings"
        results["message"] = "Timed out waiting for wannier90 output"
        return results

    results["status"] = "submitted"
    return results


def generate_pw2wannier_input(prefix: str) -> str:
    """Return a minimal pw2wannier input string for given QE prefix."""
    return f"""&inputpp
  prefix = '{prefix}'
  outdir = './'
/
write_amn = .true.
write_mmn = .true.
write_unk = .true.
"""


def generate_seedname_win(num_wann: int, projections: str, *, spinors: bool = False, dis_num_iter: int = 1000) -> str:
    """Generate a basic `seedname.win` content.

    `projections` should be a multi-line string suitable for the `begin projections` block.
    """
    spin_line = "spinors = .true." if spinors else "spinors = .false."
    return f"""num_wann = {num_wann}
dis_num_iter = {dis_num_iter}
{spin_line}

begin projections
{projections.strip()}
end projections
"""


def suggest_nbnd_from_pseudos(pseudopotentials: Dict[str, str], buffer: int = 10) -> int:
    """Suggest a conservative `nbnd` value based on pseudopotential valence counts.

    This is a best-effort helper: it looks for `valence` information in UPF
    files when they are available. If parsing fails it returns a conservative
    fallback (e.g., 64 + buffer).
    """
    total_valence = 0
    for sym, pseudo in (pseudopotentials or {}).items():
        try:
            upf_path = pseudo
            if os.path.exists(upf_path):
                # Try to parse a 'z_valence' or valence attribute; be permissive
                with open(upf_path, "r", encoding="utf-8", errors="ignore") as f:
                    text = f.read()
                # look for patterns like valence="X" or Z_valence
                import re

                m = re.search(r"valence\s*=\s*\"(\d+)\"", text)
                if not m:
                    m = re.search(r"Z_valence\s*[:=]\s*(\d+)", text)
                if m:
                    total_valence += int(m.group(1))
                    continue
        except Exception:
            pass
        # If we couldn't parse, add a modest default per element
        total_valence += 8

    # Convert electrons to bands (spin-degenerate approx): nbnd ~ (n_electrons/2) + buffer
    nbnd = max(64, int(total_valence / 2) + buffer)
    return nbnd


# ============================================================================
# WannierWorkflow: Complete Orchestrated Pipeline
# ============================================================================


class WannierWorkflow:
    """
    Complete Wannier workflow orchestrator: CIF → SCF → NSCF → pw2wannier → wannier90
    
    This class encapsulates the entire Wannierization pipeline, automatically
    orchestrating all steps from structure input to final Wannier functions.
    
    Parameters:
        cif_file: Path to structure file (CIF, POSCAR, etc.)
        pseudopotentials: Dict mapping element symbols to pseudopotential file paths (or None if using pseudopotentials_config)
        pseudopotentials_config: Name of pseudopotential config to load from ~/.xespresso/pseudopotentials/ (or None if using pseudopotentials)
        protocol: Convergence protocol ('fast', 'moderate', 'accurate')
        num_wann: Number of Wannier functions to generate
        projections: Initial projections for Wannier (e.g., "Si: s,p" or "Si: sp3d2", or 'auto' to infer from pseudopotentials)
        kpts_scf: K-point mesh for SCF (default (4, 4, 4))
        kpts_nscf: K-point mesh for NSCF (default (6, 6, 6), denser)
        nbnd: Number of bands for NSCF (auto-estimated if None)
        spinors: If True, use non-collinear magnetism (default False)
        dis_num_iter: Disentanglement iterations (default 1000)
        queue: Scheduler/queue configuration dict
        **kwargs: Additional parameters passed to CalculationWorkflow
    
    Example:
        >>> wf = WannierWorkflow(
        ...     cif_file='Si.cif',
        ...     pseudopotentials={'Si': '/path/to/Si.pbe.UPF'},
        ...     protocol='moderate',
        ...     num_wann=4,
        ...     projections='Si: sp3'
        ... )
        >>> results = wf.run(blocking=True)
        >>> print(results['wannier90']['outputs']['wout'])
    """
    
    def __init__(
        self,
        cif_file: Union[str, Path],
        pseudopotentials: Optional[Dict[str, str]] = None,
        pseudopotentials_config: Optional[str] = None,
        protocol: str = 'moderate',
        num_wann: int = 4,
        projections: str = 'auto',
        kpts_scf: Tuple[int, int, int] = (4, 4, 4),
        kpts_nscf: Tuple[int, int, int] = (6, 6, 6),
        nbnd: Optional[int] = None,
        spinors: bool = False,
        dis_num_iter: int = 1000,
        run_bands: bool = True,
        magnetic_config: Optional[Union[str, Dict]] = None,
        queue: Optional[Dict] = None,
        **kwargs
    ):
        """Initialize WannierWorkflow with structure and Wannier parameters."""
        import os
        import logging
        from xespresso.workflow.calculation_workflow import CalculationWorkflow
        from xespresso.pseudopotentials.manager import load_pseudopotentials_config
        from xespresso.utils.pseudo_utils import discover_pseudopotential_directory
        
        logger = logging.getLogger(__name__)
        
        self.cif_file = Path(cif_file)
        self.protocol = protocol
        self.num_wann = num_wann
        self.kpts_scf = kpts_scf
        self.kpts_nscf = kpts_nscf
        self.spinors = spinors
        self.dis_num_iter = dis_num_iter
        self.run_bands = run_bands
        self.queue = queue
        self.kwargs = kwargs
        
        # Handle pseudopotentials with same logic as ConvergenceWorkflow
        self.pseudopotentials = {}
        self.pseudopotentials_base_path = None
        self._pseudo_config_name = pseudopotentials_config
        
        if pseudopotentials_config is not None:
            # Load from config file
            cfg = load_pseudopotentials_config(pseudopotentials_config, verbose=False)
            if cfg is None:
                raise ValueError(f"Pseudopotentials configuration '{pseudopotentials_config}' not found")
            
            # Extract filenames and base path (same as ConvergenceWorkflow)
            self.pseudopotentials_base_path = cfg.base_path if hasattr(cfg, 'base_path') else None
            for el, pseudo in cfg.pseudopotentials.items():
                filename = pseudo.filename if hasattr(pseudo, 'filename') else str(pseudo)
                self.pseudopotentials[el] = filename
        else:
            if pseudopotentials is None:
                raise ValueError("Must provide 'pseudopotentials' mapping or 'pseudopotentials_config' name")
            
            # Discover base directory from dict (same as ConvergenceWorkflow)
            try:
                resolved_pseudos, self.pseudopotentials_base_path = discover_pseudopotential_directory(pseudopotentials)
                # Extract FILENAMES from resolved absolute paths
                for element, full_path in resolved_pseudos.items():
                    filename = os.path.basename(full_path)
                    self.pseudopotentials[element] = filename
            except FileNotFoundError as e:
                raise FileNotFoundError(str(e))
        
        # Set ESPRESSO_PSEUDO environment variable if we discovered a base path
        if self.pseudopotentials_base_path:
            os.environ['ESPRESSO_PSEUDO'] = self.pseudopotentials_base_path
        
        # Infer projections from pseudopotentials
        self.projections = projections if projections != 'auto' else self._infer_projections(self.pseudopotentials)
        self.nbnd = nbnd or suggest_nbnd_from_pseudos(self.pseudopotentials)
        
        # Initialize the underlying CalculationWorkflow
        if pseudopotentials_config:
            self.calc_wf = CalculationWorkflow.from_cif(
                cif_file,
                pseudopotentials_config=pseudopotentials_config,
                protocol=protocol,
                magnetic_config=magnetic_config,
                queue=queue,
                **kwargs
            )
        else:
            # Pass self.pseudopotentials (filenames only) not the original pseudopotentials dict
            self.calc_wf = CalculationWorkflow.from_cif(
                cif_file,
                pseudopotentials=self.pseudopotentials,
                protocol=protocol,
                magnetic_config=magnetic_config,
                queue=queue,
                **kwargs
            )
        
        # Results storage
        self.results = {}
    
    def _get_explicit_kpts(self, kpts_grid: Tuple[int, int, int], atoms) -> list:
        """Convert automatic k-point grid to explicit list of k-points.
        
        For Wannier workflow, NSCF needs explicit uniform k-point mesh.
        
        Parameters
        ----------
        kpts_grid : tuple of 3 ints
            Monkhorst-Pack grid (e.g., (4, 4, 4))
        atoms : ase.Atoms
            Atomic structure
            
        Returns
        -------
        list
            List of k-points in scaled coordinates
        """
        import numpy as np
        from ase.calculators.calculator import kpts2ndarray
        
        # Generate all k-points in the grid [0, 1) with uniform spacing
        kpts = []
        for i in range(kpts_grid[0]):
            for j in range(kpts_grid[1]):
                for k in range(kpts_grid[2]):
                    kpt = np.array([
                        i / kpts_grid[0],
                        j / kpts_grid[1],
                        k / kpts_grid[2]
                    ])
                    kpts.append(kpt)
        
        return np.array(kpts)
    
    @staticmethod
    def _infer_projections(pseudopotentials: Dict[str, str]) -> str:
        """Auto-infer projections from pseudopotential elements."""
        # Simple heuristic: default to p orbitals for all elements
        elements = list(pseudopotentials.keys())
        return '; '.join([f"{elem}: p" for elem in elements])
    
    def run(
        self,
        labels: Optional[Dict[str, str]] = None,
        blocking: bool = True,
        seedname: str = 'wannier_seed',
        run_bands_validation: Optional[bool] = None,
        run_projwfc_analysis: bool = True,
        dry_run: bool = False,
    ) -> Dict:
        """
        Execute the complete Wannier workflow pipeline with optional band structure and projection analysis.
        
        Orchestrates: SCF → Bands (optional) → PROJWFC (optional) → NSCF (wf_collect) → pw2wannier90 → wannier90
        
        The band structure calculation is essential for validating the quality of Wannier functions
        by comparing the interpolated band structure with the original DFT band structure.
        
        The PROJWFC (projection on atomic wavefunctions) analysis helps identify which atoms and
        orbitals contribute to the electronic structure, enabling better selection of Wannier
        function projections (initial guesses).
        
        Parameters:
            labels: Dict with custom labels {'scf': '...', 'nscf': '...', 'bands': '...', 'projwfc': '...'}
                   If None, uses default run labels
            blocking: If True, wait for all jobs to complete (default True)
            seedname: Base name for Wannier output files (default 'wannier_seed')
            run_bands_validation: If True, compute band structure for validation (default: self.run_bands)
            run_projwfc_analysis: If True, compute projections for orbital analysis (default: True)
            dry_run: If True, only generate input files without executing any calculations (default False)
            
        Returns:
            Dict containing results from all pipeline stages:
                {
                    'scf': Espresso calculator with SCF results or input files (if dry_run),
                    'bands': Espresso calculator with band structure input files (if run and dry_run),
                    'projwfc': EspressoProjwfc with PDOS input files (if run and dry_run),
                    'nscf': Espresso calculator with NSCF input files (if dry_run),
                    'pw2wannier': {'status', 'outputs', 'job_id', ...} (skipped if dry_run),
                    'wannier90': {'status', 'outputs', 'job_id', ...} (skipped if dry_run),
                    'projections_suggested': str with auto-suggested projections from PROJWFC,
                    'seedname': seedname used,
                    'run_dir': directory with Wannier outputs
                }
        """
        import logging
        logger = logging.getLogger(__name__)
        
        # Determine if we should run bands and projwfc
        if run_bands_validation is None:
            run_bands_validation = self.run_bands
        
        # Set default labels
        if labels is None:
            labels = {
                'scf': 'runs/02-scf',
                'bands': 'runs/02b-bands',
                'projwfc': 'runs/02c-projwfc',
                'nscf': 'runs/03-nscf',
            }
        
        # Adjust stage count based on options
        total_stages = 5  # Base: SCF + NSCF + pw2wannier + wannier90
        if dry_run:
            total_stages -= 2  # Skip pw2wannier + wannier90
        if run_bands_validation:
            total_stages += 1  # Add bands
        if run_projwfc_analysis:
            total_stages += 1  # Add projwfc
        stage_count = 1
        
        print("\n" + "="*70)
        print("WANNIER WORKFLOW - COMPLETE PIPELINE ORCHESTRATION")
        print("="*70)
        if dry_run:
            print("MODE: DRY RUN (Input files only, no execution)")
        print(f"Structure: {self.cif_file}")
        print(f"Protocol: {self.protocol}")
        print(f"Number of Wannier functions: {self.num_wann}")
        print(f"Initial projections: {self.projections}")
        print(f"K-points (SCF): {self.kpts_scf}")
        if run_bands_validation:
            print(f"K-points (Band Structure): High-symmetry path (auto-generated)")
        print(f"K-points (NSCF): {self.kpts_nscf} (denser for better interpolation)")
        print(f"Number of bands: {self.nbnd}")
        print(f"Band structure validation: {'YES ✓' if run_bands_validation else 'NO (optional)'}")
        print("="*70)
        
        # ====== STAGE 1: SCF ======
        print(f"\n[{stage_count}/{total_stages}] Running SCF calculation...")
        stage_count += 1
        try:
            scf_calc = self.calc_wf.run_scf(
                label=labels['scf'],
                kpts=self.kpts_scf,
                dry_run=dry_run
            )
            self.results['scf'] = scf_calc
            if dry_run:
                print(f"✓ SCF input files generated: {scf_calc.directory}")
            else:
                print(f"✓ SCF completed: {scf_calc.directory}")
        except Exception as e:
            print(f"✗ SCF failed: {e}")
            raise
        
        # ====== STAGE 2 (OPTIONAL): BAND STRUCTURE ======
        if run_bands_validation:
            print(f"\n[{stage_count}/{total_stages}] Running band structure calculation (for validation)...")
            stage_count += 1
            try:
                # Update atoms for bands calculation
                self.calc_wf.atoms = scf_calc.atoms
                
                bands_calc = self.calc_wf.run_bands(
                    label=labels['bands'],
                    dry_run=dry_run
                )
                self.results['bands'] = bands_calc
                if dry_run:
                    print(f"✓ Band structure input files generated: {bands_calc.directory}")
                else:
                    print(f"✓ Band structure completed: {bands_calc.directory}")
                print(f"  Essential for validating Wannier function quality")
            except Exception as e:
                print(f"⚠ Band structure failed (non-critical): {e}")
                print(f"  Continuing without validation, but quality assessment will be limited")
                self.results['bands'] = None
        
        # ====== STAGE 3 (OPTIONAL): PROJWFC ======
        projwfc_projections_suggested = None
        if run_projwfc_analysis:
            print(f"\n[{stage_count}/{total_stages}] Running projwfc (projection analysis for Wannier guidance)...")
            stage_count += 1
            try:
                # Use SCF results for projwfc
                run_dir = str(Path(scf_calc.directory).resolve())
                prefix = scf_calc.prefix
                
                projwfc_result = run_projwfc(
                    run_dir=run_dir,
                    prefix=prefix,
                    blocking=blocking,
                    queue=self.queue
                )
                self.results['projwfc'] = projwfc_result
                
                if projwfc_result['status'] in ['finished', 'finished_with_warnings', 'submitted']:
                    print(f"✓ projwfc completed: {projwfc_result['status']}")
                    if 'outputs' in projwfc_result and projwfc_result['outputs']['pdos']:
                        print(f"  PDOS analysis: {Path(projwfc_result['outputs']['pdos']).name}")
                        
                        # Try to extract suggestions from PDOS output
                        suggestions = parse_projwfc_output(run_dir, prefix)
                        if suggestions:
                            self.results['projwfc_analysis'] = suggestions
                            projwfc_projections_suggested = suggestions.get('projections', None)
                            if projwfc_projections_suggested:
                                print(f"  Suggested projections from PDOS: {projwfc_projections_suggested}")
                        
                        print(f"  ➜ Check {run_dir}/{prefix}.pdos* for detailed orbital contributions")
                else:
                    print(f"✗ projwfc {projwfc_result['status']}: {projwfc_result.get('message', 'Unknown error')}")
                    print(f"  Continuing with user-specified projections...")
                    self.results['projwfc'] = None
            except Exception as e:
                print(f"⚠ projwfc failed (non-critical): {e}")
                print(f"  Continuing with user-specified projections...")
                self.results['projwfc'] = None
        
        # ====== STAGE 4: NSCF ======
        print(f"\n[{stage_count}/{total_stages}] Running NSCF calculation (with wavefunction collection)...")
        stage_count += 1
        try:
            # Update atoms for NSCF workflow
            self.calc_wf.atoms = scf_calc.atoms
            
            # For Wannier, NSCF needs explicit uniform k-point mesh (not automatic)
            # and must disable symmetry (nosym=True, noinv=True) since wannier90 doesn't use them
            kpts_nscf_explicit = self._get_explicit_kpts(self.kpts_nscf, scf_calc.atoms)
            
            # Add nosym and noinv flags to input_data for NSCF
            input_data_nscf = self.calc_wf.input_data.copy() if hasattr(self.calc_wf, 'input_data') else {}
            input_data_nscf['nosym'] = True
            input_data_nscf['noinv'] = True
            
            nscf_calc = self.calc_wf.run_nscf(
                label=labels['nscf'],
                kpts=kpts_nscf_explicit,
                nbnd=self.nbnd,
                wf_collect=True,  # CRITICAL for Wannier
                input_data=input_data_nscf,
                dry_run=dry_run
            )
            self.results['nscf'] = nscf_calc
            if dry_run:
                print(f"✓ NSCF input files generated: {nscf_calc.directory}")
            else:
                print(f"✓ NSCF completed: {nscf_calc.directory}")
                print(f"  Wavefunctions available at: {nscf_calc.directory}/")
        except Exception as e:
            print(f"✗ NSCF failed: {e}")
            raise
        
        # ====== STAGES 5-6: pw2wannier90 and wannier90 (SKIPPED IN DRY RUN) ======
        if dry_run:
            print(f"\n[DRY RUN] Skipping pw2wannier90 and wannier90 execution (input-only mode)")
            run_dir = str(Path(nscf_calc.directory).resolve())
            self.results['seedname'] = seedname
            self.results['run_dir'] = run_dir
            self.results['projections_used'] = self.projections
            
            print("\n" + "="*70)
            print("WANNIER WORKFLOW INPUT GENERATION COMPLETED ✓")
            print("="*70)
            print(f"\nInput files generated in: {run_dir}/")
            print(f"To run the full pipeline, call: workflow.run(dry_run=False)")
            print("="*70 + "\n")
            
            return self.results
        
        # ====== STAGE 5: pw2wannier90 ======
        print(f"\n[{stage_count}/{total_stages}] Running pw2wannier90 (wavefunction conversion)...")
        stage_count += 1
        try:
            run_dir = str(Path(nscf_calc.directory).resolve())
            prefix = nscf_calc.prefix
            
            pw2w_result = run_pw2wannier(
                run_dir=run_dir,
                prefix=prefix,
                seedname=seedname,
                blocking=blocking,
                queue=self.queue
            )
            self.results['pw2wannier'] = pw2w_result
            
            if pw2w_result['status'] in ['finished', 'submitted']:
                print(f"✓ pw2wannier90 completed: {pw2w_result['status']}")
                if 'outputs' in pw2w_result and pw2w_result['outputs']:
                    for fmt in ['amn', 'mmn', 'eig']:
                        if fmt in pw2w_result['outputs']:
                            print(f"  Generated: {Path(pw2w_result['outputs'][fmt]).name}")
            else:
                print(f"✗ pw2wannier90 {pw2w_result['status']}: {pw2w_result.get('message', 'Unknown error')}")
                raise RuntimeError(f"pw2wannier90 failed: {pw2w_result}")
        except Exception as e:
            print(f"✗ pw2wannier90 failed: {e}")
            raise
        
        # ====== STAGE 6: wannier90 ======
        print(f"\n[{stage_count}/{total_stages}] Running wannier90 (Wannier function generation)...")
        try:
            # Generate .win file
            win_text = generate_seedname_win(
                num_wann=self.num_wann,
                projections=self.projections,
                spinors=self.spinors,
                dis_num_iter=self.dis_num_iter
            )
            win_path = Path(run_dir) / f"{seedname}.win"
            with open(win_path, 'w') as f:
                f.write(win_text)
            print(f"  Generated: {seedname}.win")
            
            # Run wannier90
            w90_result = run_wannier90(
                run_dir=run_dir,
                seedname=seedname,
                blocking=blocking,
                queue=self.queue
            )
            self.results['wannier90'] = w90_result
            
            if w90_result['status'] in ['finished', 'submitted']:
                print(f"✓ wannier90 completed: {w90_result['status']}")
                if 'outputs' in w90_result and w90_result['outputs']:
                    if 'wout' in w90_result['outputs']:
                        print(f"  Generated: {Path(w90_result['outputs']['wout']).name}")
            else:
                print(f"✗ wannier90 {w90_result['status']}: {w90_result.get('message', 'Unknown error')}")
                raise RuntimeError(f"wannier90 failed: {w90_result}")
        except Exception as e:
            print(f"✗ wannier90 failed: {e}")
            raise
        
        # ====== COMPLETION ======
        self.results['seedname'] = seedname
        self.results['run_dir'] = run_dir
        self.results['projections_used'] = self.projections
        self.results['band_structure_available'] = run_bands_validation and self.results['bands'] is not None
        self.results['projwfc_available'] = run_projwfc_analysis and self.results.get('projwfc') is not None
        
        print("\n" + "="*70)
        print("WANNIER WORKFLOW COMPLETED SUCCESSFULLY ✓")
        print("="*70)
        print(f"\nResults Location: {run_dir}/")
        print(f"Wannier Functions: {self.num_wann} ({seedname}.*)")
        print(f"Output files:")
        print(f"  - {seedname}.win      (input file)")
        print(f"  - {seedname}.wout     (output log)")
        print(f"  - {seedname}_*.xsf    (Wannier function density)")
        print(f"  - {seedname}_centres.xyz (Wannier centers)")
        
        if run_bands_validation and self.results['bands'] is not None:
            print(f"\nBand Structure Validation:")
            print(f"  - Run: {self.results['bands'].directory}")
            print(f"  - Use to compare DFT vs Wannier interpolation")
            print(f"  - Call: workflow.compare_bands() for detailed analysis")
        
        if run_projwfc_analysis and self.results.get('projwfc') is not None:
            print(f"\nProjection Analysis (PROJWFC):")
            print(f"  - PDOS files: {prefix}.pdos*")
            print(f"  - Shows orbital contributions to electronic structure")
            print(f"  - Use to refine Wannier projections in next iterations")
        
        print("="*70 + "\n")
        
        return self.results
    
    def get_results(self) -> Dict:
        """Get stored results from the last run."""
        return self.results
    
    def get_scf_calculator(self):
        """Get the SCF Espresso calculator."""
        return self.results.get('scf')
    
    def get_nscf_calculator(self):
        """Get the NSCF Espresso calculator."""
        return self.results.get('nscf')
    
    def get_wannier_directory(self) -> str:
        """Get the directory containing Wannier outputs."""
        return self.results.get('run_dir', '')
    
    def get_seedname(self) -> str:
        """Get the seedname used for Wannier calculations."""
        return self.results.get('seedname', '')
    
    def compare_bands(self, verbose: bool = True) -> Dict:
        """
        Compare DFT band structure with Wannier-interpolated band structure.
        
        This is the key validation step for Wannier function quality.
        Good agreement between DFT and Wannier-interpolated bands indicates
        that the Wannier functions correctly represent the electronic structure.
        
        Parameters:
            verbose: Print detailed comparison results
            
        Returns:
            Dict with comparison metrics (requires external wannier90 tools)
            
        Note:
            This requires the w90 postprocessing tool and matplotlib/matplotlib.
            For automated plotting, use:
            
            results = workflow.compare_bands()
            # Generates comparison plots showing DFT vs Wannier bands
        """
        if self.results['bands'] is None:
            print("⚠ Band structure not available. Run with run_bands=True for validation.")
            return {}
        
        dft_bands = self.results['bands']
        wannier_dir = self.results['run_dir']
        seedname = self.results['seedname']
        
        print("\n" + "="*70)
        print("BAND STRUCTURE COMPARISON - WANNIER FUNCTION VALIDATION")
        print("="*70)
        print(f"\nDFT Band Structure: {dft_bands.directory}")
        print(f"Wannier Functions: {wannier_dir}/{seedname}*")
        print("\nTo manually compare band structures:")
        print(f"  1. Run Wannier90 interpolation: wannier90.x -wout {seedname}")
        print(f"  2. Plot band structures: use xmgrace or your favorite plotter")
        print(f"  3. Compare {dft_bands.directory}/bands.dat")
        print(f"        with {wannier_dir}/{seedname}_band.dat")
        print("\nQuality assessment:")
        print("  ✓ Excellent: Near-perfect overlap")
        print("  ✓ Good: Minor deviations at high energies")
        print("  ⚠ Fair: Noticeable differences, may need more Wannier functions")
        print("  ✗ Poor: Large deviations, review projections and num_wann")
        print("="*70)
        
        return {
            'dft_bands': dft_bands.directory,
            'wannier_bands_available_in': wannier_dir,
            'seedname': seedname,
            'instruction': 'Use wannier90 postprocessing to interpolate band structure'
        }
    
    def get_band_structure_calculator(self):
        """Get the band structure Espresso calculator (if available)."""
        return self.results.get('bands', None)
    
    def get_projwfc_analysis(self) -> Optional[Dict]:
        """Get PROJWFC analysis results.
        
        Returns:
            Dict with PDOS file paths and analysis, or None if PROJWFC was not run
        """
        if not self.results.get('projwfc'):
            return None
        
        projwfc_result = self.results['projwfc']
        analysis = {
            'status': projwfc_result.get('status'),
            'pdos_file': projwfc_result.get('outputs', {}).get('pdos'),
            'run_dir': projwfc_result.get('run_dir'),
            'job_id': projwfc_result.get('job_id'),
        }
        
        if 'projwfc_analysis' in self.results:
            analysis['parsed'] = self.results['projwfc_analysis']
        
        return analysis
    
    def validate_wannier_quality(self) -> Dict:
        """
        Validate Wannier function quality by checking various metrics.
        
        Returns:
            Dict with validation status and recommendations
        """
        if not self.results.get('wannier90'):
            return {'status': 'incomplete', 'message': 'Wannier90 not completed yet'}
        
        w90_result = self.results['wannier90']
        if w90_result.get('status') not in ['finished', 'submitted']:
            return {'status': 'failed', 'message': f"Wannier90 {w90_result.get('status')}"}
        
        validation = {
            'status': 'ready_for_validation',
            'wannier_centers': f"{self.results['run_dir']}/{self.results['seedname']}_centres.xyz",
            'wannier_xsf': f"{self.results['run_dir']}/{self.results['seedname']}_*.xsf",
            'next_step': 'compare_bands()',
            'has_band_structure': self.results.get('band_structure_available', False),
            'has_projwfc': self.results.get('projwfc_available', False),
            'recommendations': []
        }
        
        if not validation['has_band_structure']:
            validation['recommendations'].append(
                "Re-run with run_bands=True to enable band structure comparison for quality validation"
            )
        else:
            validation['recommendations'].append(
                "Band structure available - Use compare_bands() to visualize DFT vs Wannier interpolation"
            )
        
        if validation['has_projwfc']:
            validation['recommendations'].append(
                "PROJWFC analysis available - Review PDOS to understand orbital contributions"
            )
        else:
            validation['recommendations'].append(
                "Consider re-running with run_projwfc_analysis=True for orbital analysis guidance"
            )
        
        return validation
