#!/usr/bin/env python
"""
Test WannierWorkflow by actually running it and inspecting real output files.
"""

from pathlib import Path
from ase.build import bulk
from ase.io import write
import tempfile
import os


def main():
    print("\n" + "="*80)
    print("TEST WANNIER WORKFLOW - REAL EXECUTION")
    print("="*80)
    
    # Create test structure in workspace for inspection
    temp_dir = Path("/home/vinicius/projects/spresso/test_wannier_output")
    
    # Clean up before starting
    import shutil
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    
    temp_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(temp_dir)
    
    print(f"\nWorking directory: {temp_dir}\n")
    
    # Create structure
    si = bulk('Si', 'diamond', a=5.431)
    cif_file = temp_dir / "Si.cif"
    write(str(cif_file), si)
    print(f"✓ Created Si structure: {cif_file}")
    
    si_pseudo = Path("/home/vinicius/pseudos/SSSP_1.3.0_PBE_efficiency/Si.pbe-n-rrkjus_psl.1.0.0.UPF")
    print(f"✓ Pseudopotential: {si_pseudo.name}")
    
    # Initialize workflow
    from xespresso.workflow.wannier_workflow import WannierWorkflow
    
    print("\n" + "="*80)
    print("INITIALIZING WANNIER WORKFLOW")
    print("="*80)
    
    wf = WannierWorkflow(
        cif_file=str(cif_file),
        pseudopotentials={'Si': str(si_pseudo)},
        protocol='fast',
        num_wann=4,
        projections='Si: sp3',
        kpts_scf=(2, 2, 2),
        kpts_nscf=(4, 4, 4),
    )
    
    print(f"✓ WannierWorkflow initialized")
    print(f"  - Protocol: {wf.protocol}")
    print(f"  - Wannier functions: {wf.num_wann}")
    print(f"  - Projections: {wf.projections}")
    print(f"  - SCF k-points: {wf.kpts_scf}")
    print(f"  - NSCF k-points: {wf.kpts_nscf}")
    print(f"  - Number of bands: {wf.nbnd}")
    
    # Run ALL stages with dry_run=True to generate all input files
    print("\n" + "="*80)
    print("GENERATING ALL WORKFLOW INPUT FILES (DRY RUN)")
    print("="*80)
    print("\nExecuting: wf.run(blocking=True, run_bands_validation=True, run_projwfc_analysis=True, dry_run=True)\n")
    print("This will generate inputs for: SCF, BANDS, PROJWFC, NSCF, PW2WANNIER, WANNIER90\n")
    
    try:
        wf.run(blocking=True, run_bands_validation=True, run_projwfc_analysis=True, dry_run=True)
        print("\n✓ All input files generated (dry_run mode)!")
    except Exception as e:
        print(f"\n✗ Error during SCF: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # List all generated files (recursive)
    print("\n" + "="*80)
    print("GENERATED FILES (RECURSIVE)")
    print("="*80)
    
    import subprocess
    result = subprocess.run(['find', str(temp_dir), '-type', 'f'], 
                          capture_output=True, text=True)
    files = sorted(result.stdout.strip().split('\n'))
    for f in files:
        if f:
            file_path = Path(f)
            rel_path = file_path.relative_to(temp_dir)
            try:
                size = file_path.stat().st_size
                print(f"  {str(rel_path):<45} ({size:>10,} bytes)")
            except:
                print(f"  {str(rel_path):<45}")
    
    # Show input file contents
    print("\n" + "="*80)
    print("INPUT FILES CONTENT")
    print("="*80)
    
    input_files = sorted(temp_dir.glob("**/[0-9]*-*.pwi"))
    if not input_files:
        input_files = sorted(temp_dir.glob("**/*.pwi"))
    
    for input_file in input_files:
        print(f"\n{'='*80}")
        print(f"📄 {input_file.relative_to(temp_dir)}")
        print(f"{'='*80}\n")
        content = input_file.read_text()
        # Show full content for input files
        if len(content) > 3000:
            print(content[:3000] + f"\n\n... ({len(content)} chars total)")
        else:
            print(content)
    
    # Show output file content if exists
    print("\n" + "="*80)
    print("OUTPUT FILES")
    print("="*80)
    
    for suffix in [".out", ".log", ".stdout"]:
        out_files = list(temp_dir.glob(f"*{suffix}"))
        for out_file in out_files:
            print(f"\n📄 {out_file.name}\n")
            content = out_file.read_text()
            # Show first 2000 chars
            if len(content) > 2000:
                print(content[:2000] + f"\n\n... ({len(content)} chars total)")
            else:
                print(content)
    
    print("\n" + "="*80)
    print("TEST COMPLETE")
    print("="*80)
    print(f"\nAll files are in: {temp_dir}/\n")


if __name__ == '__main__':
    main()
