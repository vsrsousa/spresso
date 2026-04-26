"""
Advanced example: Manual orchestration of Wannier workflow using
CalculationWorkflow and wannier_workflow helpers separately.

This approach provides more fine-grained control over each stage but requires
more manual work. For most use cases, prefer WannierWorkflow class instead.

This demonstrates the sequence:
  1. (optional) vc-relax
  2. scf
  3. nscf (wf_collect)
  4. run `pw2wannier90`
  5. run `wannier90`

Adjust paths, pseudopotential mapping and machine configuration as needed.
"""
from pathlib import Path

from xespresso.workflow.calculation_workflow import CalculationWorkflow
from xespresso.workflow.wannier_workflow import (
    run_pw2wannier,
    run_wannier90,
    generate_seedname_win,
    suggest_nbnd_from_pseudos,
)


def main():
    # Example inputs (user must adapt)
    cif_file = "structure.cif"
    pseudos = {"Si": "/path/to/Si.pbe.UPF"}  # adjust
    protocol = "moderate"
    queue = None  # use local direct scheduler by default

    # 1) Build workflow from CIF
    print("[1] Creating workflow from CIF...")
    wf = CalculationWorkflow.from_cif(cif_file, pseudos, protocol=protocol, queue=queue)

    # 2) (optional) relax
    # relax_calc = wf.run_relax(label='runs/01-vc-relax', relax_type='vc-relax')

    # 3) SCF
    print("[2] Running SCF...")
    scf_calc = wf.run_scf(label="runs/02-scf")
    print("SCF finished. dir:", scf_calc.directory)

    # 4) NSCF: create a new workflow from the converged atoms
    print("[3] Running NSCF with wavefunction collection...")
    atoms_after = wf.get_atoms()
    nbnd = suggest_nbnd_from_pseudos(pseudos)  # Uses wannier_workflow helper
    
    # Create new CalculationWorkflow for NSCF
    nscf_wf = CalculationWorkflow(atoms_after, pseudos, protocol=protocol, queue=queue)
    nscf_calc = nscf_wf.run_nscf(
        label="runs/03-nscf",
        kpts=(6, 6, 6),
        wf_collect=True,  # CRITICAL: collect wavefunctions for Wannier
        nbnd=nbnd
    )
    print("NSCF submitted/finished. dir:", nscf_calc.directory)

    # 5) pw2wannier90
    print("[4] Running pw2wannier90...")
    run_dir = nscf_calc.directory
    prefix = nscf_calc.prefix
    seedname = "wannier_seed"
    res_pw2 = run_pw2wannier(run_dir, prefix, seedname, blocking=True, queue=queue)
    print("pw2wannier result:", res_pw2)
    
    # Check for errors
    if res_pw2['status'] not in ['finished', 'submitted']:
        print(f"ERROR: pw2wannier90 failed with status '{res_pw2['status']}'")
        print(f"Message: {res_pw2.get('message', 'Unknown error')}")
        return False

    # 6) Generate .win file and run wannier90
    print("[5] Running wannier90...")
    projections = "Si: sp3"
    win_text = generate_seedname_win(
        num_wann=4,
        projections=projections,
        spinors=False,
        dis_num_iter=1000
    )
    
    # Write .win file
    win_path = Path(run_dir) / f"{seedname}.win"
    with open(win_path, "w") as f:
        f.write(win_text)
    print(f"Generated: {win_path}")

    # Run wannier90
    res_wann = run_wannier90(run_dir, seedname, blocking=True, queue=queue)
    print("wannier90 result:", res_wann)
    
    # Check for errors
    if res_wann['status'] not in ['finished', 'submitted']:
        print(f"ERROR: wannier90 failed with status '{res_wann['status']}'")
        print(f"Message: {res_wann.get('message', 'Unknown error')}")
        return False

    print("\n" + "="*70)
    print("WANNIER WORKFLOW COMPLETED ✓")
    print("="*70)
    print(f"Results in: {run_dir}/")
    print(f"Seedname: {seedname}")
    print("="*70)
    
    return True


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
