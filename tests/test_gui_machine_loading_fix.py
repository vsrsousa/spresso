"""
Test for GUI machine loading fix.

This test verifies that the load_machine function is called correctly
in the GUI pages (calculation_setup.py and workflow_builder.py).
"""

import unittest
import tempfile
import json
import os


class TestGUIMachineLoadingFix(unittest.TestCase):
    """Test suite for GUI machine loading fix."""

    def setUp(self):
        """Set up test fixtures."""
        # Create temporary directories for test
        self.temp_dir = tempfile.mkdtemp()
        self.machines_dir = os.path.join(self.temp_dir, 'machines')
        os.makedirs(self.machines_dir, exist_ok=True)
        
        # Create a test machine configuration
        self.test_machine_name = 'test_machine'
        self.test_machine_file = os.path.join(self.machines_dir, f'{self.test_machine_name}.json')
        
        machine_config = {
            'name': self.test_machine_name,
            'execution': 'local',
            'scheduler': None,
            'use_modules': False,
            'modules': [],
            'prepend': [],
            'postpend': [],
            'environment': {}
        }
        
        with open(self.test_machine_file, 'w') as f:
            json.dump(machine_config, f)

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_load_machine_signature(self):
        """Test that load_machine has the expected signature."""
        from xespresso.machines.config.loader import load_machine
        import inspect
        
        sig = inspect.signature(load_machine)
        param_names = list(sig.parameters.keys())
        
        # Verify parameter order
        self.assertEqual(param_names[0], 'config_path', 
                        "First parameter should be 'config_path'")
        self.assertEqual(param_names[1], 'machine_name',
                        "Second parameter should be 'machine_name'")
        self.assertEqual(param_names[2], 'machines_dir',
                        "Third parameter should be 'machines_dir'")
        self.assertEqual(param_names[3], 'return_object',
                        "Fourth parameter should be 'return_object'")

    def test_load_machine_with_correct_parameters(self):
        """Test that load_machine works with correct parameters."""
        from xespresso.machines.config.loader import load_machine
        
        # This is the correct way to call load_machine (as fixed in the PR)
        machine = load_machine(
            config_path=os.path.join(self.temp_dir, 'machines.json'),
            machine_name=self.test_machine_name,
            machines_dir=self.machines_dir,
            return_object=True
        )
        
        self.assertIsNotNone(machine, "Machine should be loaded successfully")
        self.assertEqual(machine.name, self.test_machine_name)

    def test_load_machine_with_keyword_arguments(self):
        """Test that load_machine works with keyword arguments."""
        from xespresso.machines.config.loader import load_machine
        
        # This is another correct way to call load_machine
        machine = load_machine(
            machine_name=self.test_machine_name,
            machines_dir=self.machines_dir,
            return_object=True
        )
        
        self.assertIsNotNone(machine, "Machine should be loaded successfully")
        self.assertEqual(machine.name, self.test_machine_name)

    def test_constants_are_available(self):
        """Test that DEFAULT_CONFIG_PATH and DEFAULT_MACHINES_DIR are importable."""
        from xespresso.machines.config.loader import (
            DEFAULT_CONFIG_PATH,
            DEFAULT_MACHINES_DIR
        )
        
        self.assertTrue(isinstance(DEFAULT_CONFIG_PATH, str))
        self.assertTrue(isinstance(DEFAULT_MACHINES_DIR, str))
        self.assertTrue(DEFAULT_CONFIG_PATH.endswith('machines.json'))
        self.assertTrue(DEFAULT_MACHINES_DIR.endswith('machines'))

    def test_list_machines_function(self):
        """Test that list_machines function works correctly."""
        from xespresso.machines.config.loader import list_machines
        
        machines = list_machines(
            config_path=os.path.join(self.temp_dir, 'machines.json'),
            machines_dir=self.machines_dir
        )
        
        self.assertIn(self.test_machine_name, machines,
                     "Test machine should be listed")


if __name__ == '__main__':
    unittest.main()
