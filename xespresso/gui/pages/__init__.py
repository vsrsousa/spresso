"""
GUI pages for the xespresso Streamlit application.

Each page module is responsible for rendering a specific section of the GUI.
"""

from .machine_config import render_machine_config_page
from .codes_config import render_codes_config_page
from .pseudopotentials_config import render_pseudopotentials_config_page
from .structure_viewer import render_structure_viewer_page
from .calculation_setup import render_calculation_setup_page
from .workflow_builder import render_workflow_builder_page
from .job_submission import render_job_submission_page
from .results_postprocessing import render_results_postprocessing_page

__all__ = [
    'render_machine_config_page',
    'render_codes_config_page',
    'render_pseudopotentials_config_page',
    'render_structure_viewer_page',
    'render_calculation_setup_page',
    'render_workflow_builder_page',
    'render_job_submission_page',
    'render_results_postprocessing_page',
]
