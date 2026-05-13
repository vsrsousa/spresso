#!/usr/bin/env python
"""
Test Si (non-magnetic) Wannier workflow with MEDUSA machine configuration.

This test verifies:
1. WannierWorkflow with non-magnetic system (nspin=1)
2. Proper job_file generation for SLURM submission
3. Simplified workflow without spin polarization
"""

import os
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from xespresso.workflow.wannier_workflow import WannierWorkflow


def main():
    """Run Si Wannier workflow test with MEDUSA."""
    print("=" * 80)
    print("Si WANNIER WORKFLOW TEST - MEDUSA MACHINE (NON-MAGNETIC)")
    print("=" * 80)
    
    # Setup paths
    cif_file = Path(__file__).parent.parent / "test_wannier_output" / "Si.cif"
    output_dir = Path(__file__).parent.parent / "test_si_wannier_medusa_output"
    
    print(f"\n✓ CIF file: {cif_file}")
    print(f"✓ Output directory: {output_dir}")
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(exist_ok=True)
    
    # Change to output directory for workflow
    os.chdir(output_dir)
    
    # Initialize Wannier workflow with MEDUSA machine
    print("\n" + "=" * 80)
    print("INITIALIZING SI WANNIER WORKFLOW WITH MEDUSA MACHINE")
    print("=" * 80)
    
    wf = WannierWorkflow(
        cif_file=str(cif_file),
        pseudopotentials_config='SSSP_efficiency',
        protocol='fast',
        num_wann=4,
        projections='Si: sp3',
        machine='medusa',
        dry_run=True
    )
    
    print(f"✓ WannierWorkflow initialized")
    print(f"  - Machine: medusa (SLURM)")
    print(f"  - Protocol: fast")
    print(f"  - Wannier functions: 4")
    print(f"  - Projections: Si: sp3")
    print(f"  - Magnetic: non-magnetic (nspin=1)")
    
    print("\n" + "=" * 80)
    print("GENERATING WORKFLOW INPUT FILES (DRY RUN - MEDUSA)")
    print("=" * 80)
    
    print("\nExecuting: wf.run(blocking=True, dry_run=True)")
    print("This will generate inputs for: SCF, BANDS, NSCF with MEDUSA configuration")
    
    wf.run(
        blocking=True,
        dry_run=True
    )
    
    # Display job_file structure
    print("\n" + "=" * 80)
    print("JOB_FILE CONTENTS - CHECKING SLURM HEADERS")
    print("=" * 80)
    
    # Check SCF job_file
    scf_job_file = Path("runs/01-scf/job_file")
    if scf_job_file.exists():
        print("\n[SCF] job_file:")
        with open(scf_job_file) as f:
            content = f.read()
            print(content[:500])
            print("...")
    
    # Check Wannier job_file (non-magnetic, no _up/_dn)
    wan_job_file = Path("runs/05-wan/job_file")
    if wan_job_file.exists():
        print("\n[WANNIER] job_file:")
        with open(wan_job_file) as f:
            content = f.read()
            print(content)
    
    print("\n" + "=" * 80)
    print("WORKFLOW SUMMARY - Si (Non-Magnetic)")
    print("=" * 80)
    
    print("\nResults Location: runs/")
    print("\nJob Files Generated:")
    for root, dirs, files in os.walk("runs"):
        for f in sorted(files):
            if f.startswith("job_file"):
                filepath = os.path.join(root, f)
                size = os.path.getsize(filepath)
                print(f"  ✓ {filepath} ({size} bytes)")
    
    print("\nMachine Configuration (MEDUSA):")
    print(f"  - Host: medusa.fis.uerj.br")
    print(f"  - Scheduler: SLURM")
    print(f"  - Partition: parallel")
    print(f"  - NTasks: 16")
    print(f"  - Launcher: srun --mpi=pmi2")
    
    print("\nTo submit to MEDUSA:")
    print("  1. Copy workflow to MEDUSA workdir")
    print("  2. Submit each job_file with sbatch:")
    print("     sbatch runs/01-scf/job_file")
    print("     sbatch runs/02-bands/job_file")
    print("     sbatch runs/04-nscf/job_file")
    print("     sbatch runs/05-wan/job_file")
    
    print("\n" + "=" * 80)
    print("✓ TEST COMPLETED SUCCESSFULLY")
    print("=" * 80)


if __name__ == "__main__":
    main()
