"""
Example script showing how to orchestrate a Wannierization pipeline using
`CalculationWorkflow` and the `wannier_workflow` helpers.

This script demonstrates the sequence:
  1. (optional) vc-relax
  2. scf
  3. nscf (wf_collect)
  4. run `pw2wannier90`
  5. run `wannier90`

Adjust paths, pseudopotential mapping and machine configuration as needed.
"""
from pathlib import Path

from xespresso.workflow.simple_workflow import CalculationWorkflow
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
    wf = CalculationWorkflow.from_cif(cif_file, pseudos, protocol=protocol, queue=queue)

    # 2) (optional) relax
    # relax_calc = wf.run_relax(label='runs/01-vc-relax', relax_type='vc-relax')

    # 3) scf
    scf_calc = wf.run_scf(label="runs/02-scf")
    print("SCF finished. dir:", scf_calc.directory)

    # 4) nscf: create a new workflow from the converged atoms
    atoms_after = wf.get_atoms()
    nbnd = suggest_nbnd_from_pseudos(pseudos)
    nscf_wf = CalculationWorkflow(atoms_after, pseudos, protocol=protocol, queue=queue)
    nscf_calc = nscf_wf.run_scf(label="runs/03-nscf", calculation="nscf", kpts=(6, 6, 6), wf_collect=True, nbnd=nbnd)
    print("NSCF submitted/finished. dir:", nscf_calc.directory)

    # 5) pw2wannier
    run_dir = nscf_calc.directory
    prefix = nscf_calc.prefix
    seedname = "wannier_seed"
    res_pw2 = run_pw2wannier(run_dir, prefix, seedname, blocking=True, queue=queue)
    print("pw2wannier result:", res_pw2)

    # 6) write a simple seedname.win and run wannier90
    projections = "Si: p"
    win_text = generate_seedname_win(num_wann=4, projections=projections)
    with open(Path(run_dir) / f"{seedname}.win", "w") as f:
        f.write(win_text)

    res_wann = run_wannier90(run_dir, seedname, blocking=True, queue=queue)
    print("wannier90 result:", res_wann)


if __name__ == "__main__":
    main()
