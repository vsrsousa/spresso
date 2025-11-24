"""
Dry run utilities for the xespresso GUI.

This module provides functionality to generate input files for Quantum ESPRESSO
calculations without actually submitting them, useful for testing and debugging.
"""

import os
import streamlit as st
from pathlib import Path

try:
    from ase import io as ase_io
    from xespresso import Espresso
    from xespresso.xio import write_espresso_in
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False


def generate_input_files(atoms, calc_params, workdir, structure_filename="structure.cif"):
    """
    Generate Quantum ESPRESSO input files for a calculation using xespresso's scheduler system.
    
    This uses xespresso's built-in file generation capabilities to create
    the necessary input files AND job script without running the calculation.
    
    Args:
        atoms: ASE Atoms object with the structure
        calc_params: Dictionary of calculation parameters for Espresso calculator
        workdir: Working directory where files will be created
        structure_filename: Name for the structure file (default: structure.cif)
    
    Returns:
        Dictionary with paths to created files, or None if failed
    """
    if not ASE_AVAILABLE:
        st.error("ASE and xespresso not available for input file generation")
        return None
    
    try:
        # Create working directory
        os.makedirs(workdir, exist_ok=True)
        
        # Save structure file
        structure_path = os.path.join(workdir, structure_filename)
        ase_io.write(structure_path, atoms)
        
        # Create Espresso calculator with the parameters
        calc = Espresso(**calc_params)
        
        # Set the calculator for the atoms
        atoms.calc = calc
        
        # Set the directory for the calculation
        calc.directory = workdir
        calc.prefix = 'espresso'
        
        # Write input file using xespresso's write_input method
        # This automatically calls set_queue and generates both the input file and job script
        calc.write_input(atoms)
        
        # Determine input file path
        input_file = os.path.join(workdir, f"{calc.prefix}.pwi")
        
        # Job file is created by the scheduler system
        job_file = os.path.join(workdir, "job_file")
        
        # Return paths to generated files
        result = {
            'structure': structure_path,
            'input': input_file,
            'workdir': workdir,
            'prefix': calc.prefix
        }
        
        # Add job file if it was created by scheduler
        if os.path.exists(job_file):
            result['job_file'] = job_file
        
        return result
        
    except Exception as e:
        st.error(f"Error generating input files: {e}")
        import traceback
        st.code(traceback.format_exc())
        return None


def preview_input_file(filepath):
    """
    Display preview of a generated input file.
    
    Args:
        filepath: Path to the input file
    """
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        
        st.subheader("Input File Preview")
        st.code(content, language='fortran')
        
        # Show file size and line count
        lines = content.split('\n')
        st.caption(f"File: {filepath} | {len(lines)} lines | {len(content)} bytes")
        
    except Exception as e:
        st.error(f"Error reading input file: {e}")


def create_job_script(workdir, machine_config, code_path, input_file, nprocs=1):
    """
    DEPRECATED: This function is no longer needed.
    
    xespresso automatically generates job scripts through its scheduler system
    when you call write_input() on the calculator. The scheduler system
    (in xespresso.scheduler and xespresso.schedulers modules) handles:
    
    - Job script generation based on the scheduler type (direct, slurm, pbs, etc.)
    - Environment setup (modules, prepend commands, xespressorc)
    - Parallel execution commands (launcher, mpirun, etc.)
    - Scheduler directives (#SBATCH, #PBS, etc.)
    - Post-execution commands (postpend)
    
    To use xespresso's scheduler system properly:
    
    1. Set up your queue configuration dictionary with scheduler parameters:
       queue = {
           'scheduler': 'slurm',  # or 'direct', 'pbs', etc.
           'execution': 'local',  # or 'remote'
           'launcher': 'mpirun -np {nprocs}',
           'modules': ['quantum-espresso'],
           'use_modules': True,
           'resources': {  # For SLURM
               'nodes': 1,
               'ntasks-per-node': nprocs,
               'time': '24:00:00'
           }
       }
    
    2. Create calculator with queue parameter:
       calc = Espresso(queue=queue, **other_params)
    
    3. Write input (this automatically generates job_file):
       calc.write_input(atoms)
    
    The job file will be created as 'job_file' in the calculator's directory.
    
    Args:
        workdir: Working directory (legacy)
        machine_config: Machine configuration dict or Machine object (legacy)
        code_path: Path to the QE executable (legacy)
        input_file: Name of the input file (legacy)
        nprocs: Number of processors (legacy)
    
    Returns:
        None - displays deprecation warning
    """
    st.warning("""
    ⚠️ **Deprecated Function**
    
    This function manually creates job scripts, which bypasses xespresso's 
    built-in scheduler system. 
    
    Instead, use xespresso's scheduler system by:
    1. Setting up a `queue` dictionary with scheduler configuration
    2. Passing it to the Espresso calculator
    3. Calling `calc.write_input(atoms)` - this automatically creates the job script
    
    See the function documentation for details.
    """)
    return None
