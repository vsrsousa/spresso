from xespresso.workflow.calculation_workflow import (
    CalculationWorkflow,
    quick_scf,
    quick_relax,
    PRESETS,
)
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow
from xespresso.workflow.wannier_workflow import WannierWorkflow

# Backward compatibility: alias for old module name
import xespresso.workflow.calculation_workflow as simple_workflow

__all__ = [
    "CalculationWorkflow",
    "ConvergenceWorkflow",
    "WannierWorkflow",
    "quick_scf",
    "quick_relax",
    "PRESETS",
    "simple_workflow",  # Deprecated, use calculation_workflow
]
