"""
Example: show how to configure and instantiate ConvergenceWorkflow using
an existing pseudopotentials configuration (no heavy calculations on import).

Usage (from shell):
    export PYTHONPATH=/home/vinicius/scratch/projects/spresso:$PYTHONPATH
    python -c "import examples.converge_example; print('example imported')"

Run the example (executes a dry-run setup only):
    python examples/converge_example.py
"""

import os
from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow
from xespresso.pseudopotentials.manager import load_pseudopotentials_config


def prepare_example(config_name: str = 'SSSP_efficiency'):
    """Prepare a ConvergenceWorkflow instance using a saved pseudopotential config.

    This function only prepares the workflow object and returns it. It does
    not run calculations automatically; call `run_convergence_study()` to run.
    """
    # Build a test structure
    atoms = bulk('Si', 'diamond', a=5.43)

    # Create the workflow by passing the configuration name directly.
    # The workflow will load the config and extract only the pseudopotentials
    # required for the structure, building absolute paths internally.
    workflow = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials_config=config_name,
        precision='medium'
    )
    return workflow


if __name__ == '__main__':
    wf = prepare_example()
    print('Prepared ConvergenceWorkflow for:', wf.atoms.get_chemical_formula())
    print('Pseudopotentials:', list(wf.pseudopotentials.keys()))
