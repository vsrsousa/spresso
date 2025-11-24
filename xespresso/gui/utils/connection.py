"""
Connection testing utilities for the xespresso GUI.

Integrates with xespresso's connection testing functionality.
"""

import streamlit as st

try:
    from xespresso.utils.auth import test_ssh_connection
    XESPRESSO_AUTH_AVAILABLE = True
except ImportError:
    XESPRESSO_AUTH_AVAILABLE = False


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
    if not XESPRESSO_AUTH_AVAILABLE:
        st.error("xespresso authentication utilities not available")
        return False
    
    try:
        success = test_ssh_connection(username, host, ssh_key, port)
        return success
    except Exception as e:
        st.error(f"Connection test failed: {e}")
        return False
