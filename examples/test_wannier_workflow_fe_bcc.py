#!/usr/bin/env python
"""
Test WannierWorkflow with Fe-bcc magnetic system.
"""

from pathlib import Path
from xespresso.tools import read_structure
import os


def main():
    print("\n" + "="*80)
    print("TEST WANNIER WORKFLOW - Fe-BCC WITH MAGNETISM")
    print("="*80)
    
    # Create test structure in workspace for inspection
    temp_dir = Path("/home/vinicius/projects/spresso/test_wannier_fe_bcc_output")
    
    # Clean up before starting
    import shutil
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    
    temp_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(temp_dir)
    
    print(f"\nWorking directory: {temp_dir}\n")
    
    # Use the pre-made Fe-bcc CIF file
    cif_file = Path("/home/vinicius/projects/spresso/Fe_bcc.cif")
    
    # Read structure using read_structure function
    atoms = read_structure(str(cif_file), primitive=True, verbose=False)
    print(f"✓ Structure loaded with read_structure()")
    print(f"  - File: {cif_file}")
    print(f"  - Atoms object: {atoms}")
    print(f"  - Number of atoms: {len(atoms)}")
    print(f"  - Chemical formula: {atoms.get_chemical_formula()}")
    print(f"  - Atomic numbers: {atoms.get_atomic_numbers()}")
    print(f"  - Scaled positions:\n{atoms.get_scaled_positions()}")
    print(f"  - Cell:\n{atoms.cell}")
    print(f"  - Volume: {atoms.get_volume():.4f} ų")
    
    fe_pseudo = Path("/home/vinicius/pseudos/SSSP_1.3.0_PBE_efficiency/Fe.pbe-n-rrkjus_psl.1.0.0.UPF")
    print(f"✓ Pseudopotential: {fe_pseudo.name}")
    
    # Initialize workflow with magnetic parameters - MUCH SIMPLER WAY!
    from xespresso.workflow.wannier_workflow import WannierWorkflow
    
    print("\n" + "="*80)
    print("INITIALIZING WANNIER WORKFLOW - MAGNETIC SYSTEM")
    print("="*80)
    
    # Use magnetic_config parameter directly (no need to call setup_magnetic_config!)
    wf = WannierWorkflow(
        cif_file=str(cif_file),
        pseudopotentials={'Fe': str(fe_pseudo)},
        protocol='fast',
        num_wann=16,
        projections='Fe: d',
        kpts_scf=(4, 4, 4),
        kpts_nscf=(6, 6, 6),
        magnetic_config={'Fe': [2.5]},  # Direct magnetic config - super clean!
    )
    
    print(f"✓ WannierWorkflow initialized with magnetism")
    print(f"  - Protocol: {wf.protocol}")
    print(f"  - Wannier functions: {wf.num_wann}")
    print(f"  - Projections: {wf.projections}")
    print(f"  - SCF k-points: {wf.kpts_scf}")
    print(f"  - NSCF k-points: {wf.kpts_nscf}")
    print(f"  - Number of bands: {wf.nbnd}")
    print(f"  - Magnetic: Fe magnetization = 2.5 µB (automatic nspin=2)")
    
    # Run with dry_run=True to generate input files
    print("\n" + "="*80)
    print("GENERATING WORKFLOW INPUT FILES (DRY RUN)")
    print("="*80)
    print("\nExecuting: wf.run(blocking=True, dry_run=True)\n")
    print("This will generate inputs for: SCF, BANDS, NSCF with magnetic settings\n")
    
    try:
        wf.run(blocking=True, dry_run=True)
        print("\n✓ All input files generated (dry_run mode)!")
        
        # Show file structure
        print("\n" + "="*80)
        print("GENERATED FILES")
        print("="*80)
        
        runs_dir = temp_dir / "runs"
        if runs_dir.exists():
            import subprocess
            result = subprocess.run(['find', str(runs_dir), '-type', 'f'], 
                                  capture_output=True, text=True)
            print(result.stdout)
            
            # Show SCF input
            scf_file = runs_dir / "02-scf" / "02-scf.pwi"
            if scf_file.exists():
                print("\n" + "="*80)
                print("SCF INPUT FILE (02-scf.pwi)")
                print("="*80)
                with open(scf_file) as f:
                    print(f.read())
                    
            # Show NSCF input
            nscf_file = runs_dir / "03-nscf" / "03-nscf.pwi"
            if nscf_file.exists():
                print("\n" + "="*80)
                print("NSCF INPUT FILE (03-nscf.pwi)")
                print("="*80)
                with open(nscf_file) as f:
                    content = f.read()
                    # Show first 50 lines and last 10 lines
                    lines = content.split('\n')
                    print('\n'.join(lines[:50]))
                    if len(lines) > 60:
                        print("\n... [k-points omitted] ...\n")
                        print('\n'.join(lines[-10:]))
                    else:
                        print('\n'.join(lines[50:]))
        
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return


if __name__ == "__main__":
    main()

