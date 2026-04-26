"""
RENAMING SUMMARY: min_ecutwfc → ecut_min, max_ecutwfc → ecut_max

All three methods now use more intuitive parameter names:
- optimize_parameters()
- run_convergence_study()
- run_convergence()

BEFORE (old names):
    wf = ConvergenceWorkflow.optimize_parameters(
        atoms=atoms,
        pseudopotentials_config='SSSP_efficiency',
        min_ecutwfc=30,
        max_ecutwfc=150,
        ecut_step=5,
    )

AFTER (new names - more intuitive!):
    wf = ConvergenceWorkflow.optimize_parameters(
        atoms=atoms,
        pseudopotentials_config='SSSP_efficiency',
        ecut_min=30,
        ecut_max=150,
        ecut_step=5,
    )

The new names are clearer and more consistent with the parameter naming convention!
"""

from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow

atoms = bulk('Au', cubic=True)

# ==========================================
# EXAMPLE 1: Simple (ecut_min, ecut_max, ecut_step)
# ==========================================
wf = ConvergenceWorkflow.optimize_parameters(
    atoms=atoms,
    pseudopotentials_config='SSSP_efficiency',
    precision='low',
    phases='both',
    ecut_min=30,         # ✨ New name!
    ecut_max=150,        # ✨ New name!
    ecut_step=5,
)


# ==========================================
# EXAMPLE 2: With n_points
# ==========================================
wf = ConvergenceWorkflow.optimize_parameters(
    atoms=atoms,
    pseudopotentials_config='SSSP_efficiency',
    precision='low',
    phases='both',
    ecut_min=40,
    ecut_max=160,
    n_ecut=10,
)


# ==========================================
# EXAMPLE 3: With dict-based range
# ==========================================
wf = ConvergenceWorkflow.optimize_parameters(
    atoms=atoms,
    pseudopotentials_config='SSSP_efficiency',
    precision='low',
    phases='both',
    ecut_range={
        'min': 30,
        'max': 150,
        'step': 7.5
    }
)


# ==========================================
# EXAMPLE 4: run_convergence_study() also updated
# ==========================================
wf = ConvergenceWorkflow(
    atoms=atoms,
    pseudopotentials_config='SSSP_efficiency',
    precision='low',
)

wf.run_convergence_study(
    ecut_min=30,
    ecut_max=150,
    ecut_step=5,
    phases='both'
)


# ==========================================
# EXAMPLE 5: run_convergence() also updated
# ==========================================
wf = ConvergenceWorkflow(
    atoms=atoms,
    pseudopotentials_config='SSSP_efficiency',
    precision='low',
)

wf.run_convergence(
    ecut_min=30,
    ecut_max=150,
    ecut_step=5,
    phases='ecut'  # Only Phase 1
)

print("✅ All parameter names updated to ecut_min and ecut_max!")
