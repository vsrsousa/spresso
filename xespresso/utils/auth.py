"""
auth.py

Manages SSH authentication and remote file operations for xespresso workflows.

This module supports:
- Key-based SSH authentication only (password-based login is no longer supported)
- Persistent SSH and SFTP sessions via paramiko
- Remote command execution
- File transfer (send/retrieve)
- Remote SHA256 checksum validation
- SSH key generation and installation via ssh-keygen and ssh-copy-id
- Connectivity testing via subprocess

Example usage:
from xespresso.utils.auth import RemoteAuth
auth = RemoteAuth(username="vinicius", host="hpc.example.com", auth_config={...})
auth.connect()
auth.send_file("local.txt", "~/remote.txt")
"""

import os
import subprocess
import paramiko
from xespresso.utils import warnings as warnings
from xespresso.utils.logging import get_logger

logger = get_logger()
warnings.apply_custom_format()

class RemoteAuth:
    """
    Manages persistent SSH authentication and file transfer for remote execution.

    Supports:
    - Key-based login only
    - Custom SSH port (default: 22)
    - Persistent SSH and SFTP sessions
    - Remote command execution
    - File transfer (send/retrieve)
    - Remote SHA256 checksum validation
    - Automatic SSH key installation on authentication failure

    Args:
        username (str): SSH login username.
        host (str): Remote machine hostname or IP.
        auth_config (dict): Authentication configuration with keys:
            - method: must be "key"
            - ssh_key: path to private key
            - port: optional SSH port (default: 22)
            - auto_install_key: optional boolean to auto-install SSH key on connection failure (default: False)
    """
    def __init__(self, username, host, auth_config):
        self.username = username
        self.host = host
        self.port = auth_config.get("port", 22)
        self.method = auth_config.get("method", "key")
        self.ssh_key = os.path.expanduser(auth_config.get("ssh_key", "~/.ssh/id_rsa"))
        self.auto_install_key = auth_config.get("auto_install_key", False)
        self.client = None
        self.sftp = None

        if self.method != "key":
            logger.error(f"Unsupported authentication method: {self.method}")
            raise ValueError(f"Unsupported authentication method: {self.method}")

    def connect(self):
        """Establishes SSH and SFTP sessions if not already connected."""
        if self.client:
            return
        try:
            self.client = paramiko.SSHClient()
            self.client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            self.client.connect(
                hostname=self.host,
                port=self.port,
                username=self.username,
                key_filename=self.ssh_key
            )
            self.sftp = self.client.open_sftp()
            logger.info(f"Connected to {self.username}@{self.host}:{self.port}")
        except paramiko.AuthenticationException as e:
            # Authentication failed - try to install SSH key if auto_install_key is enabled
            if self.auto_install_key:
                logger.warning(f"Authentication failed. Attempting to install SSH key...")
                try:
                    # Derive public key path from private key path
                    public_key_path = self.ssh_key + ".pub"
                    if not os.path.exists(public_key_path):
                        msg = f"Public key file not found at {public_key_path}. Cannot auto-install key."
                        logger.error(msg)
                        raise RuntimeError(msg)
                    
                    # Try to install the SSH key
                    install_ssh_key(self.username, self.host, public_key_path, self.port)
                    logger.info(f"SSH key installed successfully. Retrying connection...")
                    
                    # Retry connection after installing key
                    self.client = paramiko.SSHClient()
                    self.client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
                    self.client.connect(
                        hostname=self.host,
                        port=self.port,
                        username=self.username,
                        key_filename=self.ssh_key
                    )
                    self.sftp = self.client.open_sftp()
                    logger.info(f"Connected to {self.username}@{self.host}:{self.port} after key installation")
                except Exception as install_error:
                    msg = f"Failed to install SSH key and connect to {self.username}@{self.host}:{self.port}: {install_error}"
                    logger.error(msg)
                    raise RuntimeError(msg)
            else:
                msg = f"Authentication failed for {self.username}@{self.host}:{self.port}: {e}"
                logger.error(msg)
                raise RuntimeError(msg)
        except Exception as e:
            msg = f"Failed to connect to {self.username}@{self.host}:{self.port}: {e}"
            logger.error(msg)
            raise RuntimeError(msg)

    def run_command(self, command):
        """Executes a shell command on the remote host."""
        try:
            self.connect()
            stdin, stdout, stderr = self.client.exec_command(command)
            return stdout.read().decode(), stderr.read().decode()
        except Exception as e:
            msg = f"Failed to execute remote command '{command}': {e}"
            logger.error(msg)
            raise RuntimeError(msg)

    def send_file(self, local_path, remote_path):
        """Transfers a file to the remote host."""
        try:
            self.connect()
            self.sftp.put(local_path, remote_path)
            logger.info(f"Sent file '{local_path}' to '{remote_path}'")
        except Exception as e:
            msg = f"Failed to send file '{local_path}' to '{remote_path}': {e}"
            logger.error(msg)
            raise RuntimeError(msg)

    def retrieve_file(self, remote_path, local_path):
        """Retrieves a file from the remote host."""
        try:
            self.connect()
            self.sftp.get(remote_path, local_path)
            logger.info(f"Retrieved file '{remote_path}' to '{local_path}'")
        except Exception as e:
            msg = f"Failed to retrieve file '{remote_path}' to '{local_path}': {e}"
            logger.error(msg)
            raise RuntimeError(msg)

    def sha256(self, remote_path):
        """Computes SHA256 checksum of a file on the remote host."""
        try:
            self.connect()
            cmd = f"sha256sum {remote_path}"
            stdout, stderr = self.run_command(cmd)
            if stderr:
                raise RuntimeError(f"Remote error: {stderr.strip()}")
            return stdout.strip().split()[0]
        except Exception as e:
            msg = f"Failed to compute SHA256 for '{remote_path}': {e}"
            logger.error(msg)
            raise RuntimeError(msg)

    def close(self):
        """Closes SSH and SFTP sessions."""
        try:
            if self.sftp:
                self.sftp.close()
            if self.client:
                self.client.close()
            logger.info(f"Closed session with {self.username}@{self.host}")
        except Exception as e:
            msg = f"Failed to close remote session: {e}"
            logger.error(msg)
            raise RuntimeError(msg)

# 🔧 Auxiliar functions

def generate_ssh_key(private_key_path: str):
    """
    Generates a new RSA SSH key pair at the specified path.

    Args:
        private_key_path (str): Path to the private key file (e.g. ~/.ssh/id_rsa)
    """
    private_key_path = os.path.expanduser(private_key_path)
    subprocess.run(["ssh-keygen", "-t", "rsa", "-b", "4096", "-f", private_key_path], check=True)
    print(f"✅ SSH key pair created at {private_key_path} and {private_key_path}.pub")
    logger.info(f"SSH key pair generated at {private_key_path}")

def install_ssh_key(username: str, host: str, public_key_path: str, port: int = 22):
    """
    Installs the public SSH key on the remote server using ssh-copy-id.

    Args:
        username (str): SSH username.
        host (str): Remote host IP or domain.
        public_key_path (str): Path to the public key file.
        port (int): SSH port number.
    """
    public_key_path = os.path.expanduser(public_key_path)
    subprocess.run(["ssh-copy-id", "-p", str(port), "-i", public_key_path, f"{username}@{host}"], check=True)
    print(f"🔐 SSH key installed on {username}@{host}:{port}")
    logger.info(f"SSH key installed on {username}@{host}:{port}")

def test_ssh_connection(username: str, host: str, key_path: str = None, port: int = 22):
    """
    Tests SSH connectivity to the remote host.

    Args:
        username (str): SSH username.
        host (str): Remote host IP or domain.
        key_path (str, optional): Path to the private key file.
        port (int): SSH port number.
    """
    key_path = os.path.expanduser(key_path) if key_path else None
    cmd = ["ssh", "-p", str(port), "-o", "PasswordAuthentication=no"]
    if key_path:
        cmd += ["-i", key_path]
    cmd += [f"{username}@{host}", "echo 'Connection successful'"]
    try:
        subprocess.run(cmd, check=True)
        logger.info(f"SSH connection to {username}@{host}:{port} successful")
        return True
    except subprocess.CalledProcessError:
        print("❌ SSH connection failed.")
        logger.warning(f"SSH connection to {username}@{host}:{port} failed")
        return False
