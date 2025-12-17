"""
Integration tests for env_setup functionality.

These tests validate the implementation without requiring full imports.
"""
import unittest
import sys
import os
from pathlib import Path

# Test the structure of the code
class TestEnvSetupImplementation(unittest.TestCase):
    """Test that env_setup is properly implemented in the codebase"""
    
    def test_manager_has_env_setup_in_detect_codes(self):
        """Verify CodesManager.detect_codes has env_setup parameter"""
        import ast
        
        from pathlib import Path
        repo_root = Path(__file__).parent.parent
        file_path = repo_root / 'xespresso' / 'codes' / 'manager.py'
        with open(file_path, 'r') as f:
            tree = ast.parse(f.read())
        
        # Find detect_codes method
        found = False
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef) and node.name == 'detect_codes':
                params = [arg.arg for arg in node.args.args]
                self.assertIn('env_setup', params, 
                             "detect_codes should have env_setup parameter")
                found = True
                break
        
        self.assertTrue(found, "detect_codes method not found")
    
    def test_machine_has_env_setup_in_init(self):
        """Verify Machine.__init__ has env_setup parameter"""
        import ast
        
        repo_root = Path(__file__).parent.parent
        file_path = repo_root / 'xespresso' / 'machines' / 'machine.py'
        with open(file_path, 'r') as f:
            tree = ast.parse(f.read())
        
        # Find Machine class and its __init__
        found = False
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name == 'Machine':
                for item in node.body:
                    if isinstance(item, ast.FunctionDef) and item.name == '__init__':
                        params = [arg.arg for arg in item.args.args]
                        self.assertIn('env_setup', params,
                                     "Machine.__init__ should have env_setup parameter")
                        found = True
                        break
        
        self.assertTrue(found, "Machine.__init__ not found")
    
    def test_remote_mixin_uses_env_setup(self):
        """Verify remote_mixin.py uses queue env_setup"""
        repo_root = Path(__file__).parent.parent
        file_path = repo_root / 'xespresso' / 'schedulers' / 'remote_mixin.py'
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Check that it uses queue.get("env_setup")
        self.assertIn('queue.get("env_setup"', content,
                     "remote_mixin should use queue.get('env_setup')")
        
        # Check default value
        self.assertIn('source /etc/profile', content,
                     "remote_mixin should have default env_setup")
    
    def test_detect_qe_codes_has_env_setup(self):
        """Verify detect_qe_codes function has env_setup parameter"""
        import ast
        
        repo_root = Path(__file__).parent.parent
        file_path = repo_root / 'xespresso' / 'codes' / 'manager.py'
        with open(file_path, 'r') as f:
            tree = ast.parse(f.read())
        
        # Find detect_qe_codes function
        found = False
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef) and node.name == 'detect_qe_codes':
                params = [arg.arg for arg in node.args.args]
                self.assertIn('env_setup', params,
                             "detect_qe_codes should have env_setup parameter")
                found = True
                break
        
        self.assertTrue(found, "detect_qe_codes function not found")
    
    def test_all_helper_functions_have_env_setup(self):
        """Verify all helper functions have env_setup parameter"""
        import ast
        
        repo_root = Path(__file__).parent.parent
        file_path = repo_root / 'xespresso' / 'codes' / 'manager.py'
        with open(file_path, 'r') as f:
            tree = ast.parse(f.read())
        
        functions_to_check = [
            '_check_module_available',
            '_detect_local',
            '_detect_remote',
            'detect_qe_version'
        ]
        
        found_functions = {}
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef) and node.name in functions_to_check:
                params = [arg.arg for arg in node.args.args]
                found_functions[node.name] = 'env_setup' in params
        
        for func_name in functions_to_check:
            self.assertIn(func_name, found_functions, 
                         f"{func_name} not found")
            self.assertTrue(found_functions[func_name],
                           f"{func_name} should have env_setup parameter")
    
    def test_env_setup_used_in_commands(self):
        """Verify env_setup is actually used in command construction"""
        repo_root = Path(__file__).parent.parent
        file_path = repo_root / 'xespresso' / 'codes' / 'manager.py'
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Check that env_prefix is constructed from env_setup
        self.assertIn('env_prefix = f"{env_setup} && "', content,
                     "env_setup should be used to build env_prefix")
        
        # Check that env_prefix is used in commands
        self.assertIn('env_prefix}', content,
                     "env_prefix should be used in command strings")
    
    def test_machine_stores_env_setup(self):
        """Verify Machine class stores env_setup attribute"""
        repo_root = Path(__file__).parent.parent
        file_path = repo_root / 'xespresso' / 'machines' / 'machine.py'
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Check that env_setup is assigned
        self.assertIn('self.env_setup = env_setup', content,
                     "Machine should store env_setup attribute")
    
    def test_machine_includes_env_setup_in_to_dict(self):
        """Verify Machine.to_dict() includes env_setup"""
        repo_root = Path(__file__).parent.parent
        file_path = repo_root / 'xespresso' / 'machines' / 'machine.py'
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Check that env_setup is added to config dict in to_dict
        self.assertIn('config["env_setup"] = self.env_setup', content,
                     "Machine.to_dict() should include env_setup")
    
    def test_machine_includes_env_setup_in_to_queue(self):
        """Verify Machine.to_queue() includes env_setup"""
        repo_root = Path(__file__).parent.parent
        file_path = repo_root / 'xespresso' / 'machines' / 'machine.py'
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Check that env_setup is added to queue dict in to_queue
        self.assertIn('queue["env_setup"] = self.env_setup', content,
                     "Machine.to_queue() should include env_setup")
    
    def test_documentation_updated(self):
        """Verify documentation mentions env_setup"""
        # Check Machine docstring
        repo_root = Path(__file__).parent.parent
        machine_file = repo_root / 'xespresso' / 'machines' / 'machine.py'
        with open(machine_file, 'r') as f:
            content = f.read()
        
        self.assertIn('env_setup', content.lower(),
                     "Machine documentation should mention env_setup")
        
        # Check RemoteExecutionMixin docstring
        repo_root = Path(__file__).parent.parent
        remote_file = repo_root / 'xespresso' / 'schedulers' / 'remote_mixin.py'
        with open(remote_file, 'r') as f:
            content = f.read()
        
        self.assertIn('env_setup', content.lower(),
                     "RemoteExecutionMixin documentation should mention env_setup")


class TestEnvSetupLogic(unittest.TestCase):
    """Test the logic flow of env_setup usage"""
    
    def test_env_setup_extraction_from_machine_config(self):
        """Verify detect_qe_codes extracts env_setup from machine config"""
        repo_root = Path(__file__).parent.parent
        file_path = repo_root / 'xespresso' / 'codes' / 'manager.py'
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Check for logic that extracts env_setup from machine object
        self.assertIn('machine_obj.env_setup', content,
                     "Should check for env_setup in machine object")
        
        # Check for prepend fallback logic
        self.assertIn('machine_obj.prepend', content,
                     "Should have fallback to prepend for env_setup")


if __name__ == '__main__':
    unittest.main(verbosity=2)
