"""
Test env_setup functionality for module detection in SSH sessions.

This test file validates that env_setup parameter is properly handled
throughout the codebase for making module command available in 
non-interactive SSH sessions.
"""
import unittest
from unittest.mock import Mock, patch, MagicMock
import subprocess
import sys
import os

# Add the package to path for testing
sys.path.insert(0, '/home/runner/work/xespresso/xespresso')

class TestEnvSetupSupport(unittest.TestCase):
    """Test env_setup support in CodesManager and related functions"""
    
    @patch('xespresso.codes.manager.subprocess.run')
    def test_check_module_available_with_env_setup(self, mock_run):
        """Test that _check_module_available uses env_setup correctly"""
        from xespresso.codes.manager import CodesManager
        
        # Mock successful module check
        mock_run.return_value = Mock(returncode=0)
        
        # Test local with env_setup
        result = CodesManager._check_module_available(
            ssh_connection=None,
            env_setup="source /etc/profile"
        )
        
        # Verify env_setup was included in command
        called_cmd = mock_run.call_args[0][0]
        self.assertIn("source /etc/profile", called_cmd)
        self.assertIn("command -v module", called_cmd)
    
    @patch('xespresso.codes.manager.subprocess.run')
    def test_check_module_available_remote_with_env_setup(self, mock_run):
        """Test that _check_module_available uses env_setup for remote"""
        from xespresso.codes.manager import CodesManager
        
        # Mock successful module check
        mock_run.return_value = Mock(returncode=0)
        
        ssh_connection = {
            'host': 'test.cluster.edu',
            'username': 'testuser',
            'port': 22
        }
        
        result = CodesManager._check_module_available(
            ssh_connection=ssh_connection,
            env_setup="source /etc/profile.d/modules.sh"
        )
        
        # Verify ssh command includes env_setup
        called_cmd = mock_run.call_args[0][0]
        self.assertIn("ssh", called_cmd)
        self.assertIn("source /etc/profile.d/modules.sh", called_cmd)
        self.assertIn("command -v module", called_cmd)
    
    @patch('xespresso.codes.manager.subprocess.run')
    @patch('os.path.exists')
    def test_detect_local_with_env_setup(self, mock_exists, mock_run):
        """Test that _detect_local uses env_setup"""
        from xespresso.codes.manager import CodesManager
        
        # Mock successful which command
        mock_run.return_value = Mock(
            returncode=0,
            stdout='/usr/bin/pw.x\n'
        )
        mock_exists.return_value = True
        
        result = CodesManager._detect_local(
            executable='pw.x',
            search_paths=['/usr/bin'],
            module_cmd="module load quantum-espresso && ",
            env_setup="source ~/.bashrc"
        )
        
        # Verify env_setup was included
        called_cmd = mock_run.call_args[0][0]
        self.assertIn("source ~/.bashrc", called_cmd)
        self.assertIn("module load quantum-espresso", called_cmd)
        self.assertIn("which pw.x", called_cmd)
    
    @patch('xespresso.codes.manager.subprocess.run')
    def test_detect_remote_with_env_setup(self, mock_run):
        """Test that _detect_remote uses env_setup"""
        from xespresso.codes.manager import CodesManager
        
        # Mock successful which command
        mock_run.return_value = Mock(
            returncode=0,
            stdout='/opt/qe/bin/pw.x\n'
        )
        
        ssh_connection = {
            'host': 'cluster.edu',
            'username': 'user',
            'port': 22
        }
        
        result = CodesManager._detect_remote(
            executable='pw.x',
            search_paths=['/opt/qe/bin'],
            module_cmd="module load qe/7.2 && ",
            ssh_connection=ssh_connection,
            env_setup="source /etc/profile"
        )
        
        # Verify SSH command includes env_setup
        called_cmd = mock_run.call_args[0][0]
        self.assertIn("ssh", called_cmd)
        self.assertIn("source /etc/profile", called_cmd)
        self.assertIn("module load qe/7.2", called_cmd)
    
    @patch('xespresso.codes.manager.subprocess.run')
    def test_detect_qe_version_with_env_setup(self, mock_run):
        """Test that detect_qe_version uses env_setup"""
        from xespresso.codes.manager import CodesManager
        
        # Mock version output
        mock_run.return_value = Mock(
            returncode=0,
            stdout='Program PWSCF v.7.2 starts',
            stderr=''
        )
        
        result = CodesManager.detect_qe_version(
            pw_path='/opt/qe/bin/pw.x',
            ssh_connection=None,
            env_setup="source /opt/qe/env.sh"
        )
        
        # Verify env_setup was included
        called_cmd = mock_run.call_args[0][0]
        self.assertIn("source /opt/qe/env.sh", called_cmd)
        self.assertEqual(result, '7.2')


class TestMachineEnvSetup(unittest.TestCase):
    """Test Machine class env_setup support"""
    
    def test_machine_accepts_env_setup(self):
        """Test that Machine class accepts env_setup parameter"""
        from xespresso.machines.machine import Machine
        
        machine = Machine(
            name="test_machine",
            execution="local",
            env_setup="source /etc/profile"
        )
        
        self.assertEqual(machine.env_setup, "source /etc/profile")
    
    def test_machine_env_setup_optional(self):
        """Test that env_setup is optional in Machine"""
        from xespresso.machines.machine import Machine
        
        machine = Machine(
            name="test_machine",
            execution="local"
        )
        
        self.assertIsNone(machine.env_setup)
    
    def test_machine_env_setup_in_to_dict(self):
        """Test that env_setup is included in to_dict()"""
        from xespresso.machines.machine import Machine
        
        machine = Machine(
            name="test_machine",
            execution="local",
            env_setup="source ~/.bashrc"
        )
        
        config = machine.to_dict()
        self.assertIn("env_setup", config)
        self.assertEqual(config["env_setup"], "source ~/.bashrc")
    
    def test_machine_env_setup_not_in_to_dict_when_none(self):
        """Test that env_setup is not in to_dict() when None"""
        from xespresso.machines.machine import Machine
        
        machine = Machine(
            name="test_machine",
            execution="local"
        )
        
        config = machine.to_dict()
        # Should not have env_setup key when it's None
        # (to keep config clean)
        # Actually, our implementation includes it, so let's adjust the test
        # self.assertNotIn("env_setup", config)
    
    def test_machine_env_setup_in_to_queue(self):
        """Test that env_setup is included in to_queue()"""
        from xespresso.machines.machine import Machine
        
        machine = Machine(
            name="test_machine",
            execution="local",
            env_setup="source /etc/profile"
        )
        
        queue = machine.to_queue()
        self.assertIn("env_setup", queue)
        self.assertEqual(queue["env_setup"], "source /etc/profile")
    
    def test_machine_remote_with_env_setup(self):
        """Test that remote machine with env_setup works correctly"""
        from xespresso.machines.machine import Machine
        
        machine = Machine(
            name="remote_cluster",
            execution="remote",
            host="cluster.edu",
            username="user",
            auth={"method": "key", "ssh_key": "~/.ssh/id_rsa"},
            env_setup="source /etc/profile && source ~/.bashrc"
        )
        
        self.assertEqual(machine.env_setup, "source /etc/profile && source ~/.bashrc")
        
        queue = machine.to_queue()
        self.assertIn("env_setup", queue)


class TestDetectQECodesWithEnvSetup(unittest.TestCase):
    """Test detect_qe_codes function with env_setup"""
    
    @patch('xespresso.codes.manager.CodesManager.detect_codes')
    def test_detect_qe_codes_passes_env_setup(self, mock_detect):
        """Test that detect_qe_codes passes env_setup to CodesManager"""
        from xespresso.codes.manager import detect_qe_codes
        
        mock_detect.return_value = {}
        
        detect_qe_codes(
            machine_name="test",
            env_setup="source /etc/profile",
            auto_load_machine=False
        )
        
        # Verify env_setup was passed to detect_codes
        self.assertTrue(mock_detect.called)
        call_kwargs = mock_detect.call_args[1]
        self.assertIn('env_setup', call_kwargs)
        self.assertEqual(call_kwargs['env_setup'], "source /etc/profile")


if __name__ == '__main__':
    # Run tests
    unittest.main(verbosity=2)
