"""
Tests for automatic SSH key installation feature in RemoteAuth.

This module tests the auto_install_key functionality that allows
RemoteAuth to automatically install SSH keys on the remote server
when authentication fails.
"""

import os
import pytest
from unittest.mock import Mock, patch, MagicMock, call
import paramiko
from xespresso.utils.auth import RemoteAuth, install_ssh_key


class TestAutoInstallKey:
    """Test suite for automatic SSH key installation on connection failure."""

    def test_auto_install_key_disabled_by_default(self):
        """Test that auto_install_key is disabled by default."""
        auth_config = {
            "method": "key",
            "ssh_key": "~/.ssh/id_rsa",
            "port": 22
        }
        remote_auth = RemoteAuth("testuser", "testhost", auth_config)
        assert remote_auth.auto_install_key is False

    def test_auto_install_key_enabled_when_set(self):
        """Test that auto_install_key can be enabled via config."""
        auth_config = {
            "method": "key",
            "ssh_key": "~/.ssh/id_rsa",
            "port": 22,
            "auto_install_key": True
        }
        remote_auth = RemoteAuth("testuser", "testhost", auth_config)
        assert remote_auth.auto_install_key is True

    @patch('xespresso.utils.auth.install_ssh_key')
    @patch('paramiko.SSHClient')
    @patch('os.path.exists')
    def test_connect_with_auto_install_on_auth_failure(self, mock_exists, mock_ssh_client, mock_install):
        """Test that SSH key is installed when authentication fails and auto_install_key is True."""
        # Setup
        auth_config = {
            "method": "key",
            "ssh_key": "~/.ssh/id_rsa",
            "port": 22,
            "auto_install_key": True
        }
        remote_auth = RemoteAuth("testuser", "testhost", auth_config)
        
        # Mock public key file exists
        mock_exists.return_value = True
        
        # Setup mock SSH client
        mock_client_instance = MagicMock()
        mock_ssh_client.return_value = mock_client_instance
        
        # First connect attempt fails with AuthenticationException
        # Second connect attempt (after key install) succeeds
        mock_client_instance.connect.side_effect = [
            paramiko.AuthenticationException("Auth failed"),
            None  # Success on second attempt
        ]
        
        # Mock SFTP
        mock_sftp = MagicMock()
        mock_client_instance.open_sftp.return_value = mock_sftp
        
        # Execute
        remote_auth.connect()
        
        # Verify install_ssh_key was called
        mock_install.assert_called_once_with("testuser", "testhost", os.path.expanduser("~/.ssh/id_rsa.pub"), 22)
        
        # Verify connection was attempted twice
        assert mock_client_instance.connect.call_count == 2
        
        # Verify SFTP session was opened
        mock_client_instance.open_sftp.assert_called_once()
        
        # Verify client is set
        assert remote_auth.client is not None
        assert remote_auth.sftp is not None

    @patch('paramiko.SSHClient')
    def test_connect_without_auto_install_raises_on_auth_failure(self, mock_ssh_client):
        """Test that authentication failure raises error when auto_install_key is False."""
        # Setup
        auth_config = {
            "method": "key",
            "ssh_key": "~/.ssh/id_rsa",
            "port": 22,
            "auto_install_key": False
        }
        remote_auth = RemoteAuth("testuser", "testhost", auth_config)
        
        # Setup mock SSH client
        mock_client_instance = MagicMock()
        mock_ssh_client.return_value = mock_client_instance
        mock_client_instance.connect.side_effect = paramiko.AuthenticationException("Auth failed")
        
        # Execute and verify
        with pytest.raises(RuntimeError) as excinfo:
            remote_auth.connect()
        
        assert "Authentication failed" in str(excinfo.value)

    @patch('xespresso.utils.auth.install_ssh_key')
    @patch('paramiko.SSHClient')
    @patch('os.path.exists')
    def test_connect_fails_when_public_key_missing(self, mock_exists, mock_ssh_client, mock_install):
        """Test that connection fails when public key file doesn't exist."""
        # Setup
        auth_config = {
            "method": "key",
            "ssh_key": "~/.ssh/id_rsa",
            "port": 22,
            "auto_install_key": True
        }
        remote_auth = RemoteAuth("testuser", "testhost", auth_config)
        
        # Mock public key file doesn't exist
        mock_exists.return_value = False
        
        # Setup mock SSH client
        mock_client_instance = MagicMock()
        mock_ssh_client.return_value = mock_client_instance
        mock_client_instance.connect.side_effect = paramiko.AuthenticationException("Auth failed")
        
        # Execute and verify
        with pytest.raises(RuntimeError) as excinfo:
            remote_auth.connect()
        
        assert "Public key file not found" in str(excinfo.value)
        
        # Verify install_ssh_key was not called
        mock_install.assert_not_called()

    @patch('xespresso.utils.auth.install_ssh_key')
    @patch('paramiko.SSHClient')
    @patch('os.path.exists')
    def test_connect_fails_when_install_key_fails(self, mock_exists, mock_ssh_client, mock_install):
        """Test that connection fails gracefully when SSH key installation fails."""
        # Setup
        auth_config = {
            "method": "key",
            "ssh_key": "~/.ssh/id_rsa",
            "port": 22,
            "auto_install_key": True
        }
        remote_auth = RemoteAuth("testuser", "testhost", auth_config)
        
        # Mock public key file exists
        mock_exists.return_value = True
        
        # Setup mock SSH client
        mock_client_instance = MagicMock()
        mock_ssh_client.return_value = mock_client_instance
        mock_client_instance.connect.side_effect = paramiko.AuthenticationException("Auth failed")
        
        # Mock install_ssh_key raises an exception
        mock_install.side_effect = Exception("ssh-copy-id failed")
        
        # Execute and verify
        with pytest.raises(RuntimeError) as excinfo:
            remote_auth.connect()
        
        assert "Failed to install SSH key" in str(excinfo.value)

    @patch('xespresso.utils.auth.install_ssh_key')
    @patch('paramiko.SSHClient')
    @patch('os.path.exists')
    def test_connect_fails_when_retry_after_install_fails(self, mock_exists, mock_ssh_client, mock_install):
        """Test that connection fails when retry fails even after successful key installation."""
        # Setup
        auth_config = {
            "method": "key",
            "ssh_key": "~/.ssh/id_rsa",
            "port": 22,
            "auto_install_key": True
        }
        remote_auth = RemoteAuth("testuser", "testhost", auth_config)
        
        # Mock public key file exists
        mock_exists.return_value = True
        
        # Setup mock SSH client
        mock_client_instance = MagicMock()
        mock_ssh_client.return_value = mock_client_instance
        
        # Both connect attempts fail
        mock_client_instance.connect.side_effect = [
            paramiko.AuthenticationException("Auth failed"),
            paramiko.AuthenticationException("Auth still failed")
        ]
        
        # Execute and verify
        with pytest.raises(RuntimeError) as excinfo:
            remote_auth.connect()
        
        assert "Failed to install SSH key and connect" in str(excinfo.value)
        
        # Verify install_ssh_key was called
        mock_install.assert_called_once()

    @patch('paramiko.SSHClient')
    def test_connect_succeeds_without_key_install_when_auth_works(self, mock_ssh_client):
        """Test that connection succeeds without key installation when authentication works."""
        # Setup
        auth_config = {
            "method": "key",
            "ssh_key": "~/.ssh/id_rsa",
            "port": 22,
            "auto_install_key": True
        }
        remote_auth = RemoteAuth("testuser", "testhost", auth_config)
        
        # Setup mock SSH client - connection succeeds on first try
        mock_client_instance = MagicMock()
        mock_ssh_client.return_value = mock_client_instance
        mock_client_instance.connect.return_value = None
        
        # Mock SFTP
        mock_sftp = MagicMock()
        mock_client_instance.open_sftp.return_value = mock_sftp
        
        # Execute
        remote_auth.connect()
        
        # Verify connection was attempted once
        assert mock_client_instance.connect.call_count == 1
        
        # Verify client is set
        assert remote_auth.client is not None
        assert remote_auth.sftp is not None

    @patch('paramiko.SSHClient')
    def test_connect_non_auth_exception_raises_immediately(self, mock_ssh_client):
        """Test that non-authentication exceptions are raised immediately."""
        # Setup
        auth_config = {
            "method": "key",
            "ssh_key": "~/.ssh/id_rsa",
            "port": 22,
            "auto_install_key": True
        }
        remote_auth = RemoteAuth("testuser", "testhost", auth_config)
        
        # Setup mock SSH client - connection fails with non-auth exception
        mock_client_instance = MagicMock()
        mock_ssh_client.return_value = mock_client_instance
        mock_client_instance.connect.side_effect = Exception("Network error")
        
        # Execute and verify
        with pytest.raises(RuntimeError) as excinfo:
            remote_auth.connect()
        
        assert "Failed to connect" in str(excinfo.value)
        assert "Network error" in str(excinfo.value)

    @patch('paramiko.SSHClient')
    def test_connect_already_connected_does_nothing(self, mock_ssh_client):
        """Test that connect() does nothing when already connected."""
        # Setup
        auth_config = {
            "method": "key",
            "ssh_key": "~/.ssh/id_rsa",
            "port": 22,
            "auto_install_key": True
        }
        remote_auth = RemoteAuth("testuser", "testhost", auth_config)
        
        # Simulate already connected
        remote_auth.client = MagicMock()
        
        # Execute
        remote_auth.connect()
        
        # Verify SSHClient was not instantiated
        mock_ssh_client.assert_not_called()

    def test_auth_config_stores_auto_install_key_value(self):
        """Test that auth_config properly stores and retrieves auto_install_key."""
        test_cases = [
            (True, True),
            (False, False),
            (None, False),  # Default case
        ]
        
        for input_value, expected_value in test_cases:
            auth_config = {
                "method": "key",
                "ssh_key": "~/.ssh/id_rsa",
                "port": 22,
            }
            
            if input_value is not None:
                auth_config["auto_install_key"] = input_value
            
            remote_auth = RemoteAuth("testuser", "testhost", auth_config)
            assert remote_auth.auto_install_key == expected_value, \
                f"Expected {expected_value} for input {input_value}, got {remote_auth.auto_install_key}"


class TestInstallSSHKeyIntegration:
    """Integration tests for install_ssh_key function with auto_install_key feature."""

    @patch('subprocess.run')
    def test_install_ssh_key_calls_ssh_copy_id(self, mock_subprocess):
        """Test that install_ssh_key calls ssh-copy-id with correct arguments."""
        install_ssh_key("testuser", "testhost", "~/.ssh/id_rsa.pub", 2222)
        
        # Verify subprocess.run was called with correct command
        mock_subprocess.assert_called_once()
        call_args = mock_subprocess.call_args
        
        # Check the command structure
        cmd = call_args[0][0]
        assert cmd[0] == "ssh-copy-id"
        assert "-p" in cmd
        assert "2222" in cmd
        assert "-i" in cmd
        assert "testuser@testhost" in cmd


class TestMachineIntegration:
    """Test integration with Machine class."""

    def test_machine_passes_auto_install_key_to_queue(self):
        """Test that Machine properly passes auto_install_key to queue."""
        from xespresso.machines.machine import Machine
        
        machine = Machine(
            name='test_machine',
            execution='remote',
            host='test.host',
            username='testuser',
            auth={
                'method': 'key',
                'ssh_key': '~/.ssh/id_rsa',
                'auto_install_key': True
            },
            scheduler='slurm',
            workdir='/scratch/test'
        )
        
        queue = machine.to_queue()
        assert 'remote_auth' in queue
        assert queue['remote_auth']['auto_install_key'] is True

    def test_machine_omits_auto_install_key_when_not_specified(self):
        """Test that Machine doesn't include auto_install_key when not specified."""
        from xespresso.machines.machine import Machine
        
        machine = Machine(
            name='test_machine',
            execution='remote',
            host='test.host',
            username='testuser',
            auth={
                'method': 'key',
                'ssh_key': '~/.ssh/id_rsa'
            },
            scheduler='slurm',
            workdir='/scratch/test'
        )
        
        queue = machine.to_queue()
        assert 'remote_auth' in queue
        assert 'auto_install_key' not in queue['remote_auth']

    def test_machine_with_auto_install_key_false(self):
        """Test that Machine properly passes auto_install_key=False."""
        from xespresso.machines.machine import Machine
        
        machine = Machine(
            name='test_machine',
            execution='remote',
            host='test.host',
            username='testuser',
            auth={
                'method': 'key',
                'ssh_key': '~/.ssh/id_rsa',
                'auto_install_key': False
            },
            scheduler='slurm',
            workdir='/scratch/test'
        )
        
        queue = machine.to_queue()
        assert 'remote_auth' in queue
        assert queue['remote_auth']['auto_install_key'] is False
