"""
test_connection_persistence.py

Tests to verify remote connection persistence behavior.
"""

import pytest
from unittest.mock import Mock, MagicMock, patch
from xespresso.schedulers.remote_mixin import RemoteExecutionMixin
from xespresso.schedulers.slurm import SlurmScheduler
from xespresso.schedulers.direct import DirectScheduler


class TestConnectionPersistence:
    """Test suite for remote connection persistence."""
    
    def setup_method(self):
        """Set up test fixtures and clear connection cache."""
        # Clear any existing connections before each test
        RemoteExecutionMixin._remote_sessions.clear()
        RemoteExecutionMixin._last_remote_path = None
    
    def teardown_method(self):
        """Clean up after each test."""
        RemoteExecutionMixin._remote_sessions.clear()
        RemoteExecutionMixin._last_remote_path = None
    
    @patch('xespresso.schedulers.remote_mixin.RemoteAuth')
    def test_connection_created_once(self, mock_remote_auth_class):
        """Test that connection is created only once for the same machine."""
        mock_remote = Mock()
        mock_remote_auth_class.return_value = mock_remote
        
        # Create first scheduler
        calc1 = Mock()
        calc1.prefix = "calc1"
        calc1.directory = "/tmp/calc1"
        calc1.parameters = {"pseudopotentials": {}}
        
        queue = {
            "execution": "remote",
            "remote_host": "cluster.edu",
            "remote_user": "user",
            "remote_auth": {"method": "key", "ssh_key": "~/.ssh/id_rsa"},
            "remote_dir": "/home/user"
        }
        
        scheduler1 = DirectScheduler(calc1, queue, "test command")
        scheduler1._setup_remote()
        
        # Verify connection was created
        assert mock_remote_auth_class.call_count == 1
        assert mock_remote.connect.call_count == 1
        
        # Create second scheduler for same machine
        calc2 = Mock()
        calc2.prefix = "calc2"
        calc2.directory = "/tmp/calc2"
        calc2.parameters = {"pseudopotentials": {}}
        
        scheduler2 = DirectScheduler(calc2, queue, "test command")
        scheduler2._setup_remote()
        
        # Verify connection was NOT created again
        assert mock_remote_auth_class.call_count == 1  # Still 1
        assert mock_remote.connect.call_count == 1     # Still 1
        
        # Both schedulers should use the same remote connection
        assert scheduler1.remote is scheduler2.remote
    
    @patch('xespresso.schedulers.remote_mixin.RemoteAuth')
    def test_different_machines_different_connections(self, mock_remote_auth_class):
        """Test that different machines get different connections."""
        # Create two different mock connections
        mock_remote1 = Mock()
        mock_remote2 = Mock()
        mock_remote_auth_class.side_effect = [mock_remote1, mock_remote2]
        
        calc = Mock()
        calc.prefix = "calc"
        calc.directory = "/tmp/calc"
        calc.parameters = {"pseudopotentials": {}}
        
        # First machine
        queue1 = {
            "execution": "remote",
            "remote_host": "cluster1.edu",
            "remote_user": "user",
            "remote_auth": {"method": "key"},
            "remote_dir": "/home/user"
        }
        
        scheduler1 = DirectScheduler(calc, queue1, "test")
        scheduler1._setup_remote()
        
        # Second machine (different host)
        queue2 = {
            "execution": "remote",
            "remote_host": "cluster2.edu",  # Different host
            "remote_user": "user",
            "remote_auth": {"method": "key"},
            "remote_dir": "/home/user"
        }
        
        scheduler2 = DirectScheduler(calc, queue2, "test")
        scheduler2._setup_remote()
        
        # Verify two connections were created
        assert mock_remote_auth_class.call_count == 2
        assert mock_remote1.connect.call_count == 1
        assert mock_remote2.connect.call_count == 1
        
        # Different remote connections
        assert scheduler1.remote is not scheduler2.remote
    
    @patch('xespresso.schedulers.remote_mixin.RemoteAuth')
    def test_different_users_different_connections(self, mock_remote_auth_class):
        """Test that same host but different users get different connections."""
        mock_remote1 = Mock()
        mock_remote2 = Mock()
        mock_remote_auth_class.side_effect = [mock_remote1, mock_remote2]
        
        calc = Mock()
        calc.prefix = "calc"
        calc.directory = "/tmp/calc"
        calc.parameters = {"pseudopotentials": {}}
        
        # User 1
        queue1 = {
            "execution": "remote",
            "remote_host": "cluster.edu",
            "remote_user": "user1",  # User 1
            "remote_auth": {"method": "key"},
            "remote_dir": "/home/user1"
        }
        
        scheduler1 = DirectScheduler(calc, queue1, "test")
        scheduler1._setup_remote()
        
        # User 2 on same host
        queue2 = {
            "execution": "remote",
            "remote_host": "cluster.edu",  # Same host
            "remote_user": "user2",         # Different user
            "remote_auth": {"method": "key"},
            "remote_dir": "/home/user2"
        }
        
        scheduler2 = DirectScheduler(calc, queue2, "test")
        scheduler2._setup_remote()
        
        # Verify two connections were created
        assert mock_remote_auth_class.call_count == 2
        assert scheduler1.remote is not scheduler2.remote
    
    @patch('xespresso.schedulers.remote_mixin.RemoteAuth')
    def test_connection_cache_key(self, mock_remote_auth_class):
        """Test that connection cache uses (host, user) as key."""
        mock_remote = Mock()
        mock_remote_auth_class.return_value = mock_remote
        
        calc = Mock()
        calc.prefix = "calc"
        calc.directory = "/tmp/calc"
        calc.parameters = {"pseudopotentials": {}}
        
        queue = {
            "execution": "remote",
            "remote_host": "cluster.edu",
            "remote_user": "user",
            "remote_auth": {"method": "key"},
            "remote_dir": "/home/user"
        }
        
        scheduler = DirectScheduler(calc, queue, "test")
        scheduler._setup_remote()
        
        # Check that the key is in the cache
        key = ("cluster.edu", "user")
        assert key in RemoteExecutionMixin._remote_sessions
        assert RemoteExecutionMixin._remote_sessions[key] is mock_remote
    
    @patch('xespresso.schedulers.remote_mixin.RemoteAuth')
    def test_path_update_tracking(self, mock_remote_auth_class):
        """Test that remote path is tracked across calls."""
        mock_remote = Mock()
        mock_remote_auth_class.return_value = mock_remote
        
        calc1 = Mock()
        calc1.prefix = "calc1"
        calc1.directory = "calc1"
        calc1.parameters = {"pseudopotentials": {}}
        
        queue = {
            "execution": "remote",
            "remote_host": "cluster.edu",
            "remote_user": "user",
            "remote_auth": {"method": "key"},
            "remote_dir": "/home/user/base"
        }
        
        scheduler1 = DirectScheduler(calc1, queue, "test")
        scheduler1._setup_remote()
        
        # Check first path: allow autogenerated suffix for uniqueness
        assert scheduler1.remote_path.startswith("/home/user/base/calc1")
        # Note: _last_remote_path might be stored as instance variable via self._last_remote_path
        # Check via the scheduler instance
        # allow autogenerated suffix for uniqueness
        assert scheduler1._last_remote_path.startswith("/home/user/base/calc1")
        
        # Second calculation with different directory
        calc2 = Mock()
        calc2.prefix = "calc2"
        calc2.directory = "calc2"
        calc2.parameters = {"pseudopotentials": {}}
        
        scheduler2 = DirectScheduler(calc2, queue, "test")
        scheduler2._setup_remote()
        
        # Path should be updated - accept autogenerated suffix
        assert scheduler2.remote_path.startswith("/home/user/base/calc2")
        assert scheduler2._last_remote_path.startswith("/home/user/base/calc2")
    
    def test_close_all_connections(self):
        """Test that close_all_connections clears the cache."""
        # Manually add a mock connection to cache
        mock_remote = Mock()
        key = ("cluster.edu", "user")
        RemoteExecutionMixin._remote_sessions[key] = mock_remote
        RemoteExecutionMixin._last_remote_path = "/some/path"
        
        # Close all connections
        RemoteExecutionMixin.close_all_connections()
        
        # Verify cleanup
        assert len(RemoteExecutionMixin._remote_sessions) == 0
        assert RemoteExecutionMixin._last_remote_path is None
        assert mock_remote.close.called
    
    @patch('xespresso.schedulers.remote_mixin.RemoteAuth')
    def test_connection_reuse_across_scheduler_types(self, mock_remote_auth_class):
        """Test that connection is reused across different scheduler types."""
        mock_remote = Mock()
        mock_remote_auth_class.return_value = mock_remote
        
        calc1 = Mock()
        calc1.prefix = "calc1"
        calc1.directory = "/tmp/calc1"
        calc1.parameters = {"pseudopotentials": {}}
        
        queue = {
            "execution": "remote",
            "scheduler": "direct",
            "remote_host": "cluster.edu",
            "remote_user": "user",
            "remote_auth": {"method": "key"},
            "remote_dir": "/home/user"
        }
        
        # Create DirectScheduler
        scheduler1 = DirectScheduler(calc1, queue, "test")
        scheduler1._setup_remote()
        
        # Create SlurmScheduler for same machine
        queue["scheduler"] = "slurm"
        calc2 = Mock()
        calc2.prefix = "calc2"
        calc2.directory = "/tmp/calc2"
        calc2.parameters = {"pseudopotentials": {}}
        
        scheduler2 = SlurmScheduler(calc2, queue, "test")
        scheduler2._setup_remote()
        
        # Should reuse the same connection
        assert mock_remote_auth_class.call_count == 1
        assert mock_remote.connect.call_count == 1
        assert scheduler1.remote is scheduler2.remote


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
