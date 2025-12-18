"""Minimal CLI runner for the new workflow tasks.

Usage:
    python -m examples.workflow_cli

This creates a simple SCF task and runs it in debug mode.
"""
from ase.build import molecule
from xespresso.workflow.tasks import ScfTask, WorkflowRunner


def main():
    atoms = molecule("H2")
    atoms.center(6)
    atoms.pbc = [True, True, True]

    task = ScfTask(
        name="h2_scf",
        params={
            "label": "calculations/scf/h2_cli",
            "pseudopotentials": {"H": "H.pbe-rrkjus_psl.1.0.0.UPF"},
            "protocol": "fast",
        },
        inputs={"atoms": atoms},
    )

    runner = WorkflowRunner(tasks=[task], context={"debug": True})
    results = runner.run_all()
    print(results)
    # Query provenance DB
    from xespresso.provenance import ProvenanceDB

    pdb = ProvenanceDB()
    recs = pdb.get_records(task.outputs['atoms'])
    print('Provenance records for H2:', len(recs))


if __name__ == "__main__":
    main()
