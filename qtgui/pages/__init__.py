"""
GUI pages for the xespresso PyQt application.

Each page module is responsible for rendering a specific section of the GUI.
"""

from .machine_config import MachineConfigPage
from .codes_config import CodesConfigPage
from .pseudopotentials_config import PseudopotentialsConfigPage
from .structure_viewer import StructureViewerPage
from .calculation_setup import CalculationSetupPage
from .workflow_builder import WorkflowBuilderPage
from .job_submission import JobSubmissionPage
from .results_postprocessing import ResultsPostprocessingPage

__all__ = [
    'MachineConfigPage',
    'CodesConfigPage',
    'PseudopotentialsConfigPage',
    'StructureViewerPage',
    'CalculationSetupPage',
    'WorkflowBuilderPage',
    'JobSubmissionPage',
    'ResultsPostprocessingPage',
]
