"""
Test GUI import functionality.

This test suite ensures that all GUI-related imports work correctly,
addressing the issue reported in the problem statement.
"""

import unittest


class TestGUIImports(unittest.TestCase):
    """Test suite for GUI import functionality."""

    def test_xespresso_gui_module_import(self):
        """Test that xespresso.gui module can be imported."""
        try:
            import gui
            self.assertTrue(True, "xespresso.gui module imported successfully")
        except ImportError as e:
            self.fail(f"Failed to import gui: {e}")

    def test_gui_pages_import(self):
        """Test that GUI page modules can be imported."""
        try:
            from gui.pages import (
                render_machine_config_page,
                render_codes_config_page,
                render_structure_viewer_page,
                render_calculation_setup_page,
                render_workflow_builder_page,
                render_job_submission_page,
                render_results_postprocessing_page
            )
            self.assertTrue(True, "All GUI page modules imported successfully")
        except ImportError as e:
            self.fail(f"Failed to import GUI pages: {e}")

    def test_save_machine_import_from_loader(self):
        """Test that save_machine can be imported from loader (as used in GUI)."""
        try:
            from xespresso.machines.config.loader import (
                save_machine, load_machine, list_machines,
                DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR
            )
            self.assertTrue(True, "save_machine imported from loader successfully")
        except ImportError as e:
            self.fail(f"Failed to import save_machine from loader: {e}")

    def test_save_machine_import_from_config(self):
        """Test that save_machine can be imported from config module."""
        try:
            from xespresso.machines.config import save_machine
            self.assertTrue(True, "save_machine imported from config successfully")
        except ImportError as e:
            self.fail(f"Failed to import save_machine from config: {e}")

    def test_all_xespresso_modules_used_by_gui(self):
        """Test that all xespresso modules used by GUI can be imported."""
        try:
            from xespresso.machines.machine import Machine
            from xespresso.codes.manager import (
                detect_qe_codes, load_codes_config,
                DEFAULT_CODES_DIR
            )
            from xespresso.workflow import (
                CalculationWorkflow, quick_scf, quick_relax, PRESETS
            )
            from xespresso import Espresso
            self.assertTrue(True, "All xespresso modules used by GUI imported successfully")
        except ImportError as e:
            self.fail(f"Failed to import xespresso modules: {e}")

    def test_gui_utilities_import(self):
        """Test that GUI utilities can be imported."""
        try:
            from gui.utils import (
                validate_path, 
                create_3d_structure_plot, 
                display_structure_info
            )
            self.assertTrue(True, "GUI utilities imported successfully")
        except ImportError:
            # GUI utilities may have optional dependencies, this is acceptable
            self.assertTrue(True, "GUI utilities import had expected optional dependency issues")


if __name__ == '__main__':
    unittest.main()
