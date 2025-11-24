"""
Connection testing utilities for the xespresso GUI.

Integrates with xespresso's connection testing functionality.
"""

import streamlit as st
import subprocess
import os

try:
    from xespresso.utils.auth import test_ssh_connection
    XESPRESSO_AUTH_AVAILABLE = True
except ImportError:
    XESPRESSO_AUTH_AVAILABLE = False
    # Fallback implementation when paramiko is not available
    def test_ssh_connection(username, host, key_path=None, port=22):
        """Fallback SSH connection test using subprocess."""
        key_path = os.path.expanduser(key_path) if key_path else None
        cmd = ["ssh", "-p", str(port), "-o", "PasswordAuthentication=no", "-o", "BatchMode=yes", "-o", "ConnectTimeout=5"]
        if key_path:
            cmd += ["-i", key_path]
        cmd += [f"{username}@{host}", "echo 'Connection successful'"]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False


def test_connection(host, username, port=22, ssh_key=None):
    """
    Test SSH connection to a remote host using xespresso's test_ssh_connection function.
    
    Args:
        host: Remote host address
        username: SSH username
        port: SSH port (default 22)
        ssh_key: Path to SSH private key (optional)
    
    Returns:
        bool: True if connection successful, False otherwise
    """
    try:
        success = test_ssh_connection(username, host, ssh_key, port)
        return success
    except Exception as e:
        st.error(f"Connection test failed: {e}")
        return False
