from xespresso.workflow.calculation_workflow import (
    CalculationWorkflow,
    quick_scf,
    quick_relax,
    PRESETS,
)

# Backward compatibility: alias for old module name
import xespresso.workflow.calculation_workflow as simple_workflow

__all__ = [
    "CalculationWorkflow",
    "quick_scf",
    "quick_relax",
    "PRESETS",
    "simple_workflow",  # Deprecated, use calculation_workflow
]
