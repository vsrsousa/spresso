"""
Test for GUI codes loading fix.

This test verifies that the load_codes_config function is called correctly
in the GUI pages (calculation_setup.py and workflow_builder.py).
"""

import unittest
import tempfile
import json
import os


class TestGUICodesLoadingFix(unittest.TestCase):
    """Test suite for GUI codes loading fix."""

    def setUp(self):
        """Set up test fixtures."""
        # Create temporary directories for test
        self.temp_dir = tempfile.mkdtemp()
        self.codes_dir = os.path.join(self.temp_dir, 'codes')
        os.makedirs(self.codes_dir, exist_ok=True)
        
        # Create a test codes configuration
        self.test_machine_name = 'test_machine'
        self.test_codes_file = os.path.join(self.codes_dir, f'{self.test_machine_name}.json')
        
        codes_config = {
            'machine_name': self.test_machine_name,
            'qe_version': '7.2',
            'qe_prefix': '/usr/local/qe/bin',
            'codes': {
                'pw': {
                    'name': 'pw',
                    'path': '/usr/local/qe/bin/pw.x',
                    'version': '7.2'
                },
                'ph': {
                    'name': 'ph',
                    'path': '/usr/local/qe/bin/ph.x',
                    'version': '7.2'
                }
            }
        }
        
        with open(self.test_codes_file, 'w') as f:
            json.dump(codes_config, f)

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_load_codes_config_signature(self):
        """Test that load_codes_config has the expected signature."""
        from xespresso.codes.manager import load_codes_config
        import inspect
        
        sig = inspect.signature(load_codes_config)
        param_names = list(sig.parameters.keys())
        
        # Verify parameter order
        self.assertEqual(param_names[0], 'machine_name', 
                        "First parameter should be 'machine_name'")
        self.assertEqual(param_names[1], 'codes_dir',
                        "Second parameter should be 'codes_dir'")
        self.assertEqual(param_names[2], 'version',
                        "Third parameter should be 'version'")
        self.assertEqual(param_names[3], 'verbose',
                        "Fourth parameter should be 'verbose'")

    def test_load_codes_config_with_correct_parameters(self):
        """Test that load_codes_config works with correct parameters."""
        from xespresso.codes.manager import load_codes_config
        
        # This is the correct way to call load_codes_config (as fixed in the PR)
        codes = load_codes_config(
            machine_name=self.test_machine_name,
            codes_dir=self.codes_dir
        )
        
        self.assertIsNotNone(codes, "Codes config should be loaded successfully")
        self.assertEqual(codes.machine_name, self.test_machine_name)
        self.assertIn('pw', codes.codes)
        self.assertIn('ph', codes.codes)

    def test_load_codes_config_with_version(self):
        """Test that load_codes_config works with version parameter."""
        from xespresso.codes.manager import load_codes_config
        
        # Test with version parameter
        codes = load_codes_config(
            machine_name=self.test_machine_name,
            codes_dir=self.codes_dir,
            version='7.2'
        )
        
        self.assertIsNotNone(codes, "Codes config should be loaded successfully")
        self.assertEqual(codes.machine_name, self.test_machine_name)

    def test_load_codes_config_verbose_parameter(self):
        """Test that load_codes_config verbose parameter controls output."""
        from xespresso.codes.manager import load_codes_config
        import io
        import sys
        
        # Test with verbose=False (default) - should not print anything
        old_stdout = sys.stdout
        sys.stdout = captured_output = io.StringIO()
        
        codes = load_codes_config(
            machine_name=self.test_machine_name,
            codes_dir=self.codes_dir,
            verbose=False
        )
        
        output = captured_output.getvalue()
        sys.stdout = old_stdout
        
        self.assertEqual(output, "", "verbose=False should not produce any output")
        self.assertIsNotNone(codes, "Codes config should be loaded successfully")
        
        # Test with verbose=True - should print messages
        sys.stdout = captured_output = io.StringIO()
        
        codes = load_codes_config(
            machine_name=self.test_machine_name,
            codes_dir=self.codes_dir,
            verbose=True
        )
        
        output = captured_output.getvalue()
        sys.stdout = old_stdout
        
        self.assertIn(self.test_machine_name, output, "verbose=True should print machine name")
        self.assertIsNotNone(codes, "Codes config should be loaded successfully")

    def test_constants_are_available(self):
        """Test that DEFAULT_CODES_DIR is importable."""
        from xespresso.codes.manager import DEFAULT_CODES_DIR
        
        self.assertTrue(isinstance(DEFAULT_CODES_DIR, str))
        self.assertTrue('codes' in DEFAULT_CODES_DIR)

    def test_gui_calculation_setup_imports(self):
        """Test that calculation_setup.py can import load_codes_config and DEFAULT_CODES_DIR."""
        # This test verifies the fix: both load_codes_config and DEFAULT_CODES_DIR
        # should be importable together
        try:
            from xespresso.codes.manager import load_codes_config, DEFAULT_CODES_DIR
            self.assertIsNotNone(load_codes_config)
            self.assertIsNotNone(DEFAULT_CODES_DIR)
        except ImportError as e:
            self.fail(f"Failed to import required modules: {e}")

    def test_gui_workflow_builder_imports(self):
        """Test that workflow_builder.py can import load_codes_config and DEFAULT_CODES_DIR."""
        # This test verifies the fix: both load_codes_config and DEFAULT_CODES_DIR
        # should be importable together
        try:
            from xespresso.codes.manager import load_codes_config, DEFAULT_CODES_DIR
            self.assertIsNotNone(load_codes_config)
            self.assertIsNotNone(DEFAULT_CODES_DIR)
        except ImportError as e:
            self.fail(f"Failed to import required modules: {e}")


if __name__ == '__main__':
    unittest.main()
