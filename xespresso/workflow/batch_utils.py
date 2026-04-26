"""
Batch submission utilities for parallel structure calculations.

Provides generic functions for submitting multiple independent structure 
calculations to HPC clusters with job schedulers (SLURM, PBS, etc) in 
parallel mode.

Optimized for workflows where calculations are independent and can run
simultaneously across multiple nodes/cores.

Key advantage over sequential submission:
- All jobs are submitted immediately to the scheduler queue
- Remote scheduler distributes jobs across available resources
- No need to wait for each job to complete before submitting the next
- Much faster overall turnaround for convergence studies

Example:
    >>> from xespresso.workflow.batch_utils import submit_structures_parallel, collect_results
    >>> 
    >>> structures = [
    ...     {'atoms': atoms1, 'label': 'vacuum_5', 'param_key': 5.0, 'num_atoms': 12},
    ...     {'atoms': atoms2, 'label': 'vacuum_7', 'param_key': 7.0, 'num_atoms': 12},
    ... ]
    >>> 
    >>> result = submit_structures_parallel(
    ...     structures=structures,
    ...     kmesh=(6, 6, 1),
    ...     pseudopotentials='default',
    ...     machine='medusa',
    ... )
    >>> 
    >>> # Collect results as jobs complete
    >>> results_dict = collect_results(result, timeout=1800)
"""

import logging
from typing import Dict, List, Tuple, Optional

logger = logging.getLogger(__name__)


def submit_structures_parallel(
    structures: List[Dict],
    kmesh: Tuple[int, int, int],
    pseudopotentials: str,
    pseudopotentials_config: Optional[str] = None,
    protocol: str = 'standard',
    precision: str = 'low',
    machine: Optional[str] = None,
    queue: Optional[Dict] = None,
    code_version: Optional[str] = None,
    input_data: Optional[Dict] = None,
) -> Dict:
    """
    Submit multiple independent structure calculations in parallel batch mode.
    
    This function is optimized for HPC clusters with schedulers (SLURM, PBS, etc).
    All jobs are submitted simultaneously, allowing the scheduler to distribute
    them across available nodes/cores in parallel.
    
    **Key advantage**: Instead of waiting for each calculation to complete before
    submitting the next, all jobs are queued at once. The remote scheduler handles
    parallelization automatically.
    
    Args:
        structures: List of dicts, each containing:
            - 'atoms': ase.Atoms structure
            - 'label': str, unique calculation label (e.g., 'vacuum_10')
            - 'param_key': str/float/int, parameter being varied (for result mapping)
            - 'num_atoms': int, number of atoms in this structure
        kmesh: (nk_x, nk_y, nk_z) k-point mesh tuple
        pseudopotentials: str or dict, pseudopotential configuration
        pseudopotentials_config: str, pseudopotential configuration name (optional)
        protocol: str, calculation protocol ('standard', 'fast', etc)
        precision: str, calculation precision ('low', 'normal', 'high')
        machine: str, remote machine name (e.g., 'medusa', 'localhost')
        queue: Dict, queue configuration for remote submission
        code_version: str, QE version code (e.g., '7.4.1')
        input_data: Dict, optional input parameters (ecutwfc, ecutrho, conv_thr, etc.)
                   These override the preset values. Example:
                   {'ecutwfc': 60.0, 'ecutrho': 480.0}
    
    Returns:
        Dict with keys:
            - 'batch_results': List of submission info dicts from submit_scf_batch()
            - 'structures_map': Dict mapping param_key → structure info
            - 'workflow': Reference to the first CalculationWorkflow
            - 'all_workflows': List of all created CalculationWorkflow instances
            - 'total_submitted': int, number of successfully submitted jobs
    
    Example:
        >>> # Test vacuum convergence with parallel submission
        >>> base_slab = atoms.copy()
        >>> structures = []
        >>> 
        >>> for vacuum in [5, 7, 9, 11]:
        ...     slab = base_slab.copy()
        ...     slab.center(vacuum=vacuum, axis=2)
        ...     structures.append({
        ...         'atoms': slab,
        ...         'label': f'vacuum_{vacuum}',
        ...         'param_key': vacuum,
        ...         'num_atoms': len(slab),
        ...     })
        >>> 
        >>> result = submit_structures_parallel(
        ...     structures=structures,
        ...     kmesh=(6, 6, 1),
        ...     pseudopotentials='default',
        ...     machine='medusa',
        ...     code_version='7.4.1',
        ... )
        >>> 
        >>> # All 4 jobs submitted in parallel, scheduler handles distribution
        >>> print(f"Submitted {result['total_submitted']} jobs")
    """
    from xespresso.workflow.calculation_workflow import CalculationWorkflow
    
    logger.info(f"\n{'='*70}")
    logger.info(f"PARALLEL BATCH SUBMISSION: {len(structures)} independent calculations")
    logger.info(f"{'='*70}")
    logger.info(f"All jobs will be submitted in parallel to remote scheduler")
    logger.info(f"Remote scheduler (e.g., SLURM) will distribute jobs across nodes\n")
    
    if not structures:
        logger.warning("No structures provided")
        return {
            'batch_results': [],
            'structures_map': {},
            'workflow': None,
            'all_workflows': [],
            'total_submitted': 0,
        }
    
    workflows = []
    batch_results_all = []
    structures_map = {}
    
    for i, struct_info in enumerate(structures):
        atoms = struct_info['atoms']
        label = struct_info['label']
        param_key = struct_info['param_key']
        
        structures_map[param_key] = struct_info
        
        try:
            # Create CalculationWorkflow for this structure
            cw = CalculationWorkflow(
                atoms=atoms,
                pseudopotentials=pseudopotentials,
                pseudopotentials_config=pseudopotentials_config,
                protocol=protocol,
                precision=precision,
                machine=machine,
                queue=queue,
                code_version=code_version,
                input_data=input_data,
            )
            
            # Submit non-blocking (job queued on remote scheduler, returns immediately)
            logger.info(f"  [{i+1}/{len(structures)}] Submitting: {label}")
            batch_result = cw.submit_scf_batch(
                label=label,
                kpts=kmesh,
                wait_for_completion=False  # Important: non-blocking
            )
            
            batch_results_all.append(batch_result)
            workflows.append(cw)
            
        except Exception as e:
            logger.warning(f"  [{i+1}/{len(structures)}] Failed to submit {label}: {e}")
            batch_results_all.append({
                'label': label,
                'param_key': param_key,
                'submitted': False,
                'error': str(e),
            })
    
    total_submitted = sum(1 for r in batch_results_all if r.get('submitted', False))
    
    logger.info(f"\n{'='*70}")
    logger.info(f"✓ Batch submission complete")
    logger.info(f"  Submitted: {total_submitted}/{len(structures)} jobs")
    logger.info(f"  Scheduler will execute in parallel as resources available")
    logger.info(f"{'='*70}\n")
    
    return {
        'batch_results': batch_results_all,
        'structures_map': structures_map,
        'workflow': workflows[0] if workflows else None,
        'all_workflows': workflows,
        'total_submitted': total_submitted,
    }


def collect_results(
    batch_result: Dict,
    timeout: int = 3600,
    poll_interval: int = 30,
    verbose: bool = True,
) -> Dict:
    """
    Collect results from a parallel batch submission.
    
    Monitors submitted jobs and extracts results as they complete.
    
    Args:
        batch_result: Dict returned from submit_structures_parallel()
        timeout: Maximum time to wait for all jobs (seconds)
        poll_interval: Time between status checks (seconds)
        verbose: Print progress information
    
    Returns:
        Dict mapping param_key → job completion info with:
            - 'success': bool, whether calculation completed successfully
            - 'energy': float, energy in eV (if available)
            - 'error': str, error message (if failed)
            - Additional properties from the calculation
    
    Example:
        >>> batch_result = submit_structures_parallel(structures, kmesh, ...)
        >>> results = collect_results(batch_result, timeout=1800)
        >>> 
        >>> for param_key, completion_info in results.items():
        ...     if completion_info['success']:
        ...         print(f"{param_key}: E = {completion_info['energy']} eV")
    """
    if not batch_result['workflow']:
        logger.warning("No workflow available in batch results")
        return {}
    
    workflow = batch_result['workflow']
    batch_results = batch_result['batch_results']
    structures_map = batch_result['structures_map']
    
    logger.info(f"\nCollecting results from {len(batch_results)} jobs...")
    logger.info(f"(This may take a while depending on calculation complexity)\n")
    
    # Wait for all jobs to complete
    completion_list = workflow.wait_for_batch_jobs(
        batch_results,
        timeout=timeout,
        poll_interval=poll_interval,
        verbose=verbose
    )
    
    # Map results back to parameter keys
    results_dict = {}
    for i, completion in enumerate(completion_list):
        if i < len(batch_results):
            batch_result_item = batch_results[i]
            label = batch_result_item.get('label', f'job_{i}')
            
            # Find corresponding structure
            param_key = None
            for pk, struct_info in structures_map.items():
                if struct_info.get('label') == label:
                    param_key = pk
                    break
            
            if param_key is not None:
                results_dict[param_key] = completion
    
    logger.info(f"\n✓ Result collection complete\n")
    
    return results_dict


def collect_relax_results(
    batch_result: Dict,
    timeout: int = 3600,
    poll_interval: int = 30,
    verbose: bool = True,
) -> Dict:
    """
    Collect results from a parallel relaxation batch submission.
    
    Monitors submitted relaxation jobs and extracts results as they complete.
    Compatible with submit_relaxations_parallel() output.
    
    Args:
        batch_result: Dict returned from submit_relaxations_parallel()
        timeout: Maximum time to wait for all jobs (seconds)
        poll_interval: Time between status checks (seconds)
        verbose: Print progress information
    
    Returns:
        Dict mapping param_key → job completion info with:
            - 'success': bool, whether calculation completed successfully
            - 'energy': float, energy in eV (if available)
            - 'error': str, error message (if failed)
            - 'converged': bool, whether relaxation converged
            - Additional properties from the calculation
    
    Example:
        >>> batch_result = submit_relaxations_parallel(structures, kmesh, ...)
        >>> results = collect_relax_results(batch_result, timeout=1800)
        >>> 
        >>> for param_key, completion_info in results.items():
        ...     if completion_info['success']:
        ...         print(f"{param_key}: E = {completion_info['energy']} eV")
    """
    batch_results = batch_result.get('batch_results', [])
    structures_map = batch_result.get('structures_map', {})
    all_workflows = batch_result.get('all_workflows', [])
    
    if not all_workflows:
        logger.warning("No workflows available in batch results")
        return {}
    
    logger.info(f"\nCollecting results from {len(batch_results)} relaxations...")
    logger.info(f"(This may take a while depending on calculation complexity)\n")
    
    # Collect results from each relaxation job
    results_dict = {}
    
    for batch_result_item in batch_results:
        param_key = batch_result_item.get('param_key')
        label = batch_result_item.get('label')
        submitted = batch_result_item.get('submitted', False)
        
        if not submitted:
            # Job failed to submit
            results_dict[param_key] = {
                'success': False,
                'energy': None,
                'error': batch_result_item.get('error', 'Submission failed'),
                'converged': False,
            }
            continue
        
        # Job was submitted, monitor and retrieve results
        try:
            calc = batch_result_item.get('calc')
            workflow = batch_result_item.get('workflow')
            
            if calc is None or workflow is None:
                results_dict[param_key] = {
                    'success': False,
                    'energy': None,
                    'error': 'No calculator or workflow available',
                    'converged': False,
                }
                continue
            
            # Monitor job using RemoteJobMonitor
            job_id = getattr(calc, 'last_job_id', None)
            if job_id:
                from xespresso.workflow.remote_job_monitor import RemoteJobMonitor
                monitor = RemoteJobMonitor(calc)
                
                if monitor.wait(timeout=timeout, poll_interval=poll_interval):
                    monitor.retrieve_output()
                    logger.info(f"  ✓ {label}: Job completed, retrieving results...")
                    
                    # Read results from completed calculation
                    try:
                        calc.read_results()
                    except Exception as e:
                        logger.debug(f"  Could not read results: {e}")
                    
                    # Extract energy using get_potential_energy()
                    try:
                        energy = workflow.atoms.get_potential_energy()
                        converged = True
                        success = True
                        error = None
                    except Exception as e:
                        logger.warning(f"  Could not extract energy for {label}: {e}")
                        energy = None
                        converged = False
                        success = False
                        error = str(e)
                else:
                    logger.warning(f"  ⏱ {label}: Job monitor timeout")
                    energy = None
                    success = False
                    converged = False
                    error = "Job monitoring timeout"
            else:
                # No job ID (possibly already cached/completed)
                try:
                    energy = workflow.atoms.get_potential_energy()
                    success = True
                    converged = True
                    error = None
                except Exception as e:
                    logger.warning(f"  Could not extract energy for {label}: {e}")
                    energy = None
                    success = False
                    converged = False
                    error = str(e)
            
            results_dict[param_key] = {
                'success': success,
                'energy': energy,
                'converged': converged,
                'error': error,
            }
            
        except Exception as e:
            logger.error(f"  ✗ Error collecting results for {label}: {e}")
            results_dict[param_key] = {
                'success': False,
                'energy': None,
                'error': str(e),
                'converged': False,
            }
    
    logger.info(f"\n✓ Relaxation result collection complete\n")
    
    return results_dict


def submit_relaxations_parallel(
    structures: List[Dict],
    kmesh: Tuple[int, int, int],
    pseudopotentials_config: str,
    relax_type: str = 'relax',
    protocol: str = 'standard',
    machine: Optional[str] = None,
    queue: Optional[Dict] = None,
    code_version: Optional[str] = None,
    input_data: Optional[Dict] = None,
    dipole_correction: bool = False,
) -> Dict:
    """
    Submit multiple relaxation calculations in parallel.
    
    Similar to submit_structures_parallel() but for structural relaxations
    instead of SCF calculations. All relaxations are submitted to the remote
    scheduler simultaneously.
    
    Args:
        structures: List of dicts with 'atoms', 'label', 'param_key', 'num_atoms'
        kmesh: (nk_x, nk_y, nk_z) k-point mesh tuple
        pseudopotentials_config: str, pseudopotential config name
        relax_type: str, 'relax' (ions only) or 'vc-relax' (ions + cell)
        protocol: str, calculation protocol
        machine: str, remote machine name
        queue: Dict, queue configuration
        code_version: str, QE version
        input_data: Dict, override input parameters (ecutwfc, dipole, etc)
        dipole_correction: bool, add dipole='z' correction
    
    Returns:
        Dict with keys:
            - 'batch_results': List of submission results
            - 'structures_map': Dict mapping param_key → structure info
            - 'workflow': Reference workflow
            - 'total_submitted': Number of successfully submitted jobs
    """
    from xespresso.workflow.calculation_workflow import CalculationWorkflow
    
    logger.info(f"\n{'='*70}")
    logger.info(f"PARALLEL RELAXATION BATCH: {len(structures)} independent relaxations")
    logger.info(f"{'='*70}")
    logger.info(f"All jobs will be submitted in parallel to remote scheduler\n")
    
    if not structures:
        logger.warning("No structures provided")
        return {
            'batch_results': [],
            'structures_map': {},
            'workflow': None,
            'total_submitted': 0,
        }
    
    workflows = []
    batch_results_all = []
    structures_map = {}
    
    # STEP 1: Submit all jobs NON-BLOCKING (parallel submission to scheduler)
    logger.info("STEP 1: Submitting all jobs in parallel (non-blocking)...\n")
    
    for i, struct_info in enumerate(structures):
        atoms = struct_info['atoms']
        label = struct_info['label']
        param_key = struct_info['param_key']
        
        structures_map[param_key] = struct_info
        
        try:
            logger.info(f"  [{i+1}/{len(structures)}] Preparing: {label}")
            
            # Create CalculationWorkflow
            cw = CalculationWorkflow(
                atoms=atoms,
                pseudopotentials_config=pseudopotentials_config,
                protocol=protocol,
                machine=machine,
                queue=queue,
                code_version=code_version,
            )
            
            # Apply input overrides
            if input_data:
                cw.input_data.update(input_data)
            if dipole_correction and 'dipole' not in cw.input_data:
                cw.input_data['dipole'] = 'z'
            
            # Submit relaxation NON-BLOCKING (job submitted to remote scheduler, returns immediately)
            logger.info(f"    Submitting to remote scheduler: {label}")
            
            calc = cw.run_relax(
                label=label,
                relax_type=relax_type,
                kpts=kmesh,
                wait_for_completion=False,  # KEY: Non-blocking submission
            )
            
            batch_results_all.append({
                'label': label,
                'param_key': param_key,
                'submitted': True,
                'workflow': cw,
                'calc': calc,
            })
            workflows.append(cw)
            
        except Exception as e:
            logger.warning(f"  [{i+1}/{len(structures)}] Failed to submit {label}: {e}")
            batch_results_all.append({
                'label': label,
                'param_key': param_key,
                'submitted': False,
                'error': str(e),
            })
    
    total_submitted = sum(1 for r in batch_results_all if r.get('submitted', False))
    
    logger.info(f"\n{'='*70}")
    logger.info(f"✓ PARALLEL SUBMISSION COMPLETE")
    logger.info(f"  Successfully submitted: {total_submitted}/{len(structures)} relaxations")
    logger.info(f"  All jobs now queued on remote scheduler")
    logger.info(f"  Remote scheduler will parallelize execution across available nodes")
    logger.info(f"{'='*70}\n")
    
    return {
        'batch_results': batch_results_all,
        'structures_map': structures_map,
        'workflow': workflows[0] if workflows else None,
        'all_workflows': workflows,
        'total_submitted': total_submitted,
    }
