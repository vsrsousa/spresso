"""
Remote connection wrapper for job monitoring.

This module provides a simplified interface for remote connections
specifically designed for the job monitor dialog.
"""

from xespresso.utils.auth import RemoteAuth


class RemoteConnection:
    """
    Wrapper class for remote connections in job monitoring context.
    
    This class provides a simplified interface compatible with the job
    monitor dialog, wrapping the RemoteAuth functionality.
    
    Args:
        queue (dict): Queue configuration dictionary containing:
            - remote_host: hostname or IP
            - remote_user: SSH username
            - remote_auth: dict with authentication config
    """
    
    def __init__(self, queue):
        """Initialize remote connection from queue configuration."""
        self.queue = queue
        self.remote_auth = RemoteAuth(
            username=queue.get('remote_user'),
            host=queue.get('remote_host'),
            auth_config=queue.get('remote_auth', {'method': 'key'})
        )
        self.remote_auth.connect()
    
    def run_command(self, command):
        """
        Execute a command on the remote host.
        
        Args:
            command (str): Command to execute
            
        Returns:
            str: Standard output from the command
        """
        stdout, stderr = self.remote_auth.run_command(command)
        return stdout
    
    def get_file(self, remote_path, local_path):
        """
        Download a file from the remote host.
        
        Args:
            remote_path (str): Path to file on remote host
            local_path (str): Local path where file should be saved
        """
        self.remote_auth.retrieve_file(remote_path, local_path)
    
    def disconnect(self):
        """Close the remote connection."""
        self.remote_auth.close()
