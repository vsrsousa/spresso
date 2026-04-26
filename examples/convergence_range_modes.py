"""
Examples of different ecutwfc range specification modes.

Shows all 5 ways to specify ecutwfc values:
1. Simple step-based (default ARANGE)
2. Simple n_points (shorthand for linspace)
3. Explicit list
4. Dict with step (ARANGE mode)
5. Dict with n_points (LINSPACE mode)
"""

from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow

# Create test structure
atoms = bulk('Au', cubic=True)

# ==========================================
# MODE 1: Simple step-based (default ARANGE)
# ==========================================
wf1 = ConvergenceWorkflow.optimize_parameters(
    atoms=atoms,
    pseudopotentials_config='SSSP_efficiency',
    precision='low',
    phases='both',
    ecut_max=150,     # max value
    ecut_step=5,      # step size → [30, 35, 40, 45, ..., 150]
)
print("MODE 1 (Simple ARANGE):", wf1._build_ecut_range(ecut_step=5))


# ==========================================
# MODE 2: Simple n_points (shorthand linspace)
# ==========================================
wf2 = ConvergenceWorkflow.optimize_parameters(
    atoms=atoms,
    pseudopotentials_config='SSSP_efficiency',
    precision='low',
    phases='both',
    ecut_max=150,
    n_ecut=10,  # 10 uniform points → [30, 43.3, 56.7, 70, 83.3, ..., 150]
)
print("MODE 2 (Simple linspace):", wf2._build_ecut_range())


# ==========================================
# MODE 3: Explicit list
# ==========================================
wf3 = ConvergenceWorkflow.optimize_parameters(
    atoms=atoms,
    pseudopotentials_config='SSSP_efficiency',
    precision='low',
    phases='both',
    ecut_vals=[30, 40, 50, 60, 70, 80, 100, 130],  # Exact values
)
print("MODE 3 (Explicit list):", wf3._build_ecut_range())


# ==========================================
# MODE 4: Dict with step (ARANGE)
# ==========================================
wf4 = ConvergenceWorkflow.optimize_parameters(
    atoms=atoms,
    pseudopotentials_config='SSSP_efficiency',
    precision='low',
    phases='both',
    ecut_range={
        'min': 30,
        'max': 150,
        'step': 7.5  # [30, 37.5, 45, 52.5, ..., 150]
    }
)
print("MODE 4 (Dict ARANGE):", wf4._build_ecut_range())


# ==========================================
# MODE 5: Dict with n_points (LINSPACE)
# ==========================================
wf5 = ConvergenceWorkflow.optimize_parameters(
    atoms=atoms,
    pseudopotentials_config='SSSP_efficiency',
    precision='low',
    phases='both',
    ecut_range={
        'min': 40,
        'max': 160,
        'n_points': 8  # 8 uniform points in [40, 160]
    }
)
print("MODE 5 (Dict linspace):", wf5._build_ecut_range())


# ==========================================
# PRIORITY ORDER (if multiple are specified)
# ==========================================
# Highest: ecut_vals (explicit list)
wf_priority = ConvergenceWorkflow.optimize_parameters(
    atoms=atoms,
    pseudopotentials_config='SSSP_efficiency',
    precision='low',
    phases='both',
    ecut_vals=[50, 75, 100],  # ← This wins!
    ecut_range={'min': 30, 'max': 150, 'step': 10},  # Ignored
    n_ecut=20,  # Ignored
    ecut_max=200,  # Ignored
    ecut_step=5,  # Ignored
)
print("PRIORITY (explicit list wins):", wf_priority._build_ecut_range())


# ==========================================
# PHASE CONTROL
# ==========================================
# Run only PHASE 1 (ecutwfc convergence)
wf_phase1 = ConvergenceWorkflow.optimize_parameters(
    atoms=atoms,
    pseudopotentials_config='SSSP_efficiency',
    precision='low',
    phases='ecut',  # Only ecutwfc convergence
    ecut_max=150,
    ecut_step=5,
)

# Run only PHASE 2 (kspacing convergence)
# Requires pre-computed optimal_ecutwfc from PHASE 1
wf_phase2 = ConvergenceWorkflow(
    atoms=atoms,
    pseudopotentials_config='SSSP_efficiency',
    precision='low',
    ecut_max=150,
    ecut_step=5,
)
# First run PHASE 1
wf_phase2.run_convergence(phases='ecut')
# Then run PHASE 2
wf_phase2.run_convergence(phases='kpt')


print("\n✅ All range modes available!")
