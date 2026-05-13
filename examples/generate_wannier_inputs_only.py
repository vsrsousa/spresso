#!/usr/bin/env python
"""
Generate and inspect all WannierWorkflow input files WITHOUT executing.
"""

from pathlib import Path
from ase.build import bulk
from ase.io import write
import tempfile
import os


def main():
    print("\n" + "="*80)
    print("WANNIER WORKFLOW - GENERATE & INSPECT INPUT FILES ONLY")
    print("="*80)
    
    # Create test structure
    temp_dir = Path(tempfile.gettempdir()) / "wannier_inputs_only"
    temp_dir.mkdir(exist_ok=True)
    os.chdir(temp_dir)
    
    print(f"\nWorking directory: {temp_dir}\n")
    
    # Create structure
    si = bulk('Si', 'diamond', a=5.431)
    cif_file = temp_dir / "Si.cif"
    write(str(cif_file), si)
    
    si_pseudo = Path("/home/vinicius/pseudos/SSSP_1.3.0_PBE_efficiency/Si.pbe-n-rrkjus_psl.1.0.0.UPF")
    
    # Initialize workflow
    from xespresso.workflow.wannier_workflow import (
        WannierWorkflow,
        generate_pw2wannier_input,
        generate_seedname_win,
    )
    
    wf = WannierWorkflow(
        cif_file=str(cif_file),
        pseudopotentials={'Si': str(si_pseudo)},
        protocol='fast',
        num_wann=4,
        projections='Si: sp3',
        kpts_scf=(2, 2, 2),
        kpts_nscf=(4, 4, 4),
    )
    
    print("✓ WannierWorkflow initialized")
    print(f"  - Wannier functions: {wf.num_wann}")
    print(f"  - Number of bands: {wf.nbnd}")
    print(f"  - SCF k-points: {wf.kpts_scf}")
    print(f"  - NSCF k-points: {wf.kpts_nscf}")
    
    # Generate all input files
    print("\n" + "="*80)
    print("GENERATING INPUT FILES (NO EXECUTION)")
    print("="*80)
    
    calc_wf = wf.calc_wf
    atoms = calc_wf.get_atoms()
    
    from espresso import Espresso
    
    # 1. SCF Input
    print("\n" + "─"*80)
    print("1️⃣  SCF Input (pwscf.scf)")
    print("─"*80)
    
    scf_calc = Espresso(
        input_data=calc_wf.input_data.copy(),
        pseudopotentials=calc_wf.original_pseudopotentials,
        kpts=wf.kpts_scf,
        calculation='scf',
    )
    scf_calc.atoms = atoms
    scf_calc.write_input(atoms, directory=temp_dir, prefix='pwscf_scf')
    
    scf_in = (temp_dir / "pwscf_scf.in").read_text()
    print(scf_in)
    (temp_dir / "pwscf_scf.in").rename(temp_dir / "pwscf.scf")
    
    # 2. NSCF Input
    print("\n" + "─"*80)
    print("2️⃣  NSCF Input (pwscf.nscf)")
    print("─"*80)
    
    nscf_input_data = calc_wf.input_data.copy()
    nscf_input_data['nbnd'] = wf.nbnd
    
    nscf_calc = Espresso(
        input_data=nscf_input_data,
        pseudopotentials=calc_wf.original_pseudopotentials,
        kpts=wf.kpts_nscf,
        calculation='nscf',
    )
    nscf_calc.atoms = atoms
    nscf_calc.write_input(atoms, directory=temp_dir, prefix='pwscf_nscf')
    
    nscf_in = (temp_dir / "pwscf_nscf.in").read_text()
    print(nscf_in)
    (temp_dir / "pwscf_nscf.in").rename(temp_dir / "pwscf.nscf")
    
    # 3. PROJWFC Input
    print("\n" + "─"*80)
    print("3️⃣  PROJWFC Input (projwfc.in)")
    print("─"*80)
    
    projwfc_input = """&inputpp
    prefix='pwscf'
    outdir='./'
    what='proj'
    ngauss=0
    degauss=0.02
/
"""
    projwfc_file = temp_dir / "projwfc.in"
    projwfc_file.write_text(projwfc_input)
    print(projwfc_input)
    
    # 4. PW2WANNIER90 Input
    print("\n" + "─"*80)
    print("4️⃣  PW2WANNIER90 Input (pw2wannier.in)")
    print("─"*80)
    
    pw2wan_input = generate_pw2wannier_input('pwscf')
    pw2wan_file = temp_dir / "pw2wannier.in"
    pw2wan_file.write_text(pw2wan_input)
    print(pw2wan_input)
    
    # 5. WANNIER90 Input
    print("\n" + "─"*80)
    print("5️⃣  WANNIER90 Input (wannier_seed.win)")
    print("─"*80)
    
    win_input = generate_seedname_win(
        num_wann=wf.num_wann,
        projections=wf.projections,
        spinors=wf.spinors,
        dis_num_iter=1000
    )
    win_file = temp_dir / "wannier_seed.win"
    win_file.write_text(win_input)
    print(win_input)
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY - INPUT FILES GENERATED")
    print("="*80)
    
    input_files = [
        temp_dir / "pwscf.scf",
        temp_dir / "pwscf.nscf",
        temp_dir / "projwfc.in",
        temp_dir / "pw2wannier.in",
        temp_dir / "wannier_seed.win",
    ]
    
    print(f"\n✓ All {len(input_files)} input files generated in: {temp_dir}/\n")
    for f in input_files:
        if f.exists():
            size = f.stat().st_size
            print(f"  ✓ {f.name:<25} ({size:>6,} bytes)")
        else:
            print(f"  ✗ {f.name:<25} (NOT FOUND)")
    
    print("\n" + "="*80)
    print("Ready to run with:")
    print("  - pw.x < pwscf.scf")
    print("  - pw.x < pwscf.nscf")
    print("  - projwfc.x < projwfc.in")
    print("  - pw2wannier90.x < pw2wannier.in")
    print("  - wannier90.x wannier_seed")
    print("="*80 + "\n")


if __name__ == '__main__':
    main()
