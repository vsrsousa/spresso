#!/usr/bin/env python
"""
Extract and inspect ACTUAL input files from WannierWorkflow.

Writes all input files to disk and shows their content.
"""

from pathlib import Path
from ase.build import bulk
from ase.io import write
import tempfile
import os


def main():
    print("\n" + "="*80)
    print("WANNIER WORKFLOW - EXTRACT & SHOW REAL INPUT FILES")
    print("="*80)
    
    # Create test structure
    temp_dir = Path(tempfile.gettempdir()) / "wannier_extract"
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
    
    print("="*80)
    print("INPUT FILES CREATED BY WANNIER WORKFLOW")
    print("="*80)
    
    # Generate SCF input file
    print("\n" + "─"*80)
    print("📄 pwscf.scf")
    print("─"*80)
    calc_wf = wf.calc_wf
    atoms_scf = calc_wf.get_atoms()
    
    # Create temporary calculator to extract SCF input
    from espresso import Espresso
    scf_calc = Espresso(
        input_data=calc_wf.input_data.copy(),
        pseudopotentials=calc_wf.original_pseudopotentials,
        kpts=wf.kpts_scf,
        calculation='scf',
    )
    scf_calc.atoms = atoms_scf
    
    # Write SCF input
    scf_input_file = temp_dir / "pwscf.scf"
    scf_calc.write_input(scf_calc.atoms, directory=temp_dir, prefix='pwscf_scf')
    
    # Read and display
    if (temp_dir / "pwscf_scf.in").exists():
        scf_content = (temp_dir / "pwscf_scf.in").read_text()
        print(scf_content)
        (temp_dir / "pwscf_scf.in").rename(scf_input_file)
        print(f"✓ Written to: {scf_input_file}")
    
    # 3. PROJWFC Input
    print("\n" + "─"*80)
    print("📄 projwfc.in")
    print("─"*80)
    projwfc_input = """&inputpp
    prefix='pwscf'
    outdir='./'
    what='proj'
    ngauss=0
    degauss=0.02
/
"""
    projwfc_path = temp_dir / "projwfc.in"
    projwfc_path.write_text(projwfc_input)
    print(projwfc_input)
    print(f"✓ Written to: {projwfc_path}")
    
    # Generate NSCF input file
    print("\n" + "─"*80)
    print("📄 pwscf.nscf")
    print("─"*80)
    nscf_input_data = calc_wf.input_data.copy()
    nscf_input_data['nbnd'] = wf.nbnd
    nscf_input_data['wf_collect'] = True
    
    nscf_calc = Espresso(
        input_data=nscf_input_data,
        pseudopotentials=calc_wf.original_pseudopotentials,
        kpts=wf.kpts_nscf,
        calculation='nscf',
    )
    nscf_calc.atoms = atoms_scf
    nscf_calc.write_input(nscf_calc.atoms, directory=temp_dir, prefix='pwscf_nscf')
    
    if (temp_dir / "pwscf_nscf.in").exists():
        nscf_content = (temp_dir / "pwscf_nscf.in").read_text()
        print(nscf_content)
        (temp_dir / "pwscf_nscf.in").rename(temp_dir / "pwscf.nscf")
        print(f"✓ Written to: {temp_dir}/pwscf.nscf")
    
    # 5. PW2WANNIER90 Input
    print("\n" + "─"*80)
    print("📄 pw2wannier.in")
    print("─"*80)
    pw2wan_input = generate_pw2wannier_input('pwscf')
    pw2wan_path = temp_dir / "pw2wannier.in"
    pw2wan_path.write_text(pw2wan_input)
    print(pw2wan_input)
    print(f"✓ Written to: {pw2wan_path}")
    
    # 6. WANNIER90 Input
    print("─"*80)
    print("📄 wannier_seed.win")
    print("─"*80)
    win_input = generate_seedname_win(
        num_wann=wf.num_wann,
        projections=wf.projections,
        spinors=wf.spinors,
        dis_num_iter=1000
    )
    win_path = temp_dir / "wannier_seed.win"
    win_path.write_text(win_input)
    print(win_input)
    print(f"✓ Written to: {win_path}")
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY - INPUT FILES CREATED")
    print("="*80)
    
    files_created = list(temp_dir.glob("*.in")) + list(temp_dir.glob("*.win")) + list(temp_dir.glob("*.scf")) + list(temp_dir.glob("*.nscf"))
    print(f"\n✓ Files written to: {temp_dir}/")
    for f in files_created:
        size = f.stat().st_size
        print(f"  - {f.name} ({size} bytes)")
    
    print("="*80 + "\n")


if __name__ == '__main__':
    main()
