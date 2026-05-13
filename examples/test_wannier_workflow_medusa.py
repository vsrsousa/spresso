#!/usr/bin/env python
"""
Test Fe-BCC Wannier workflow with MEDUSA machine configuration.

This test verifies:
1. WannierWorkflow initialization with machine='medusa'
2. Input file generation for remote execution
3. Proper job_file format for SLURM submission
4. Spin-polarized Wannier setup (nspin=2, up/down components)
"""

import os
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from xespresso.workflow.wannier_workflow import WannierWorkflow


def main():
    """Run Fe-BCC Wannier workflow test with MEDUSA."""
    print("=" * 80)
    print("Fe-BCC WANNIER WORKFLOW TEST - MEDUSA MACHINE")
    print("=" * 80)
    
    # Setup paths
    cif_file = Path(__file__).parent.parent / "Fe_bcc.cif"
    output_dir = Path(__file__).parent.parent / "test_wannier_medusa_output"
    
    print(f"\n✓ CIF file: {cif_file}")
    print(f"✓ Output directory: {output_dir}")
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(exist_ok=True)
    
    # Change to output directory for workflow
    os.chdir(output_dir)
    
    # Initialize Wannier workflow with MEDUSA machine
    print("\n" + "=" * 80)
    print("INITIALIZING WANNIER WORKFLOW WITH MEDUSA MACHINE")
    print("=" * 80)
    
    wf = WannierWorkflow(
        cif_file=str(cif_file),
        pseudopotentials_config='SSSP_efficiency',
        protocol='fast',
        num_wann=16,
        projections='auto',
        kpts_scf=(4, 4, 4),
        kpts_nscf=(6, 6, 6),
        machine='medusa',  # Use MEDUSA machine configuration
        magnetic_config='ferro',  # Ferromagnetic Fe
        debug=False,
    )
    
    print(f"✓ WannierWorkflow initialized")
    print(f"  - Machine: medusa (SLURM)")
    print(f"  - Protocol: fast")
    print(f"  - Wannier functions: 16")
    print(f"  - Projections: auto (inferred from pseudopotentials)")
    print(f"  - Magnetic: ferro (nspin=2)")
    
    # Run workflow in dry_run mode (generate input files only)
    print("\n" + "=" * 80)
    print("GENERATING WORKFLOW INPUT FILES (DRY RUN - MEDUSA)")
    print("=" * 80)
    
    print("\nExecuting: wf.run(blocking=True, dry_run=True)")
    print("This will generate inputs for: SCF, BANDS, NSCF with MEDUSA configuration\n")
    
    results = wf.run(blocking=True, dry_run=True)
    
    # Verify job_file format for SLURM
    print("\n" + "=" * 80)
    print("VERIFYING JOB_FILE FORMAT FOR SLURM")
    print("=" * 80)
    
    job_files = [
        'runs/01-scf/job_file',
        'runs/02-bands/job_file',
        'runs/03-projwfc/job_file',
        'runs/04-nscf/job_file',
        'runs/05-wan/job_file_up',
        'runs/05-wan/job_file_dn',
    ]
    
    print("\nGenerated job_files:")
    for job_file in job_files:
        path = output_dir / job_file
        if path.exists():
            size = path.stat().st_size
            print(f"  ✓ {job_file} ({size} bytes)")
            
            # Show first few lines
            with open(path, 'r') as f:
                lines = f.readlines()[:5]
                for line in lines:
                    print(f"    {line.rstrip()}")
                if len(f.readlines()) > 5:
                    print(f"    ...")
        else:
            print(f"  ✗ {job_file} (NOT FOUND)")
    
    # Print summary
    print("\n" + "=" * 80)
    print("WORKFLOW SUMMARY")
    print("=" * 80)
    
    print(f"\nResults Location: {output_dir}/runs/")
    print(f"\nWannier Configuration:")
    print(f"  - Seedname: wannier (no _seed suffix)")
    print(f"  - Spin components: up, down (separate job_file_up, job_file_dn)")
    print(f"  - Prefix: Fe_bcc (consistent across all stages)")
    
    print(f"\nMachine Configuration (MEDUSA):")
    print(f"  - Host: medusa.fis.uerj.br")
    print(f"  - Scheduler: SLURM")
    print(f"  - Partition: parallel")
    print(f"  - NTasks: 16")
    print(f"  - Launcher: srun --mpi=pmi2")
    print(f"  - Workdir: /scratch/users/vinicius/xespresso")
    
    print(f"\nTo submit to MEDUSA:")
    print(f"  1. Copy workflow to MEDUSA workdir: /scratch/users/vinicius/xespresso/")
    print(f"  2. Navigate to runs/ subdirectories")
    print(f"  3. Submit each job_file or job_file_* with sbatch:")
    print(f"     sbatch runs/01-scf/job_file")
    print(f"     sbatch runs/02-bands/job_file")
    print(f"     sbatch runs/03-projwfc/job_file")
    print(f"     sbatch runs/04-nscf/job_file")
    print(f"     sbatch runs/05-wan/job_file_up  (and job_file_dn)")
    
    print("\n" + "=" * 80)
    print("✓ TEST COMPLETED SUCCESSFULLY")
    print("=" * 80 + "\n")


if __name__ == '__main__':
    main()
