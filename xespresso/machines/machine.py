"""
machine.py

Machine class for xespresso workflows.

This module provides an object-oriented interface for machine configurations,
encapsulating all properties and providing validation, serialization, and
convenient access methods.

Example usage:
    from xespresso.machines.machine import Machine
    
    # Create from dict
    config = {...}
    machine = Machine.from_dict("my_cluster", config)
    
    # Access properties
    print(machine.name)
    print(machine.execution)
    print(machine.is_remote)
    
    # Convert to queue dict
    queue = machine.to_queue()
"""

import os
import json
from typing import Dict, List, Optional, Union


class Machine:
    """
    Represents a machine configuration for xespresso calculations.
    
    A Machine encapsulates all configuration needed to execute calculations
    either locally or remotely, including scheduler settings, authentication,
    resource requirements, and environment setup.
    
    Attributes:
        name (str): Machine identifier
        execution (str): Execution mode ("local" or "remote")
        scheduler (str): Scheduler type ("direct", "slurm", etc.)
        workdir (str): Working directory path
        use_modules (bool): Whether to use environment modules
        modules (List[str]): List of modules to load
        resources (Dict): Scheduler resource requirements
        prepend (Union[str, List[str]]): Commands to run before execution
        postpend (Union[str, List[str]]): Commands to run after execution
        launcher (str): MPI launcher command template
        nprocs (int): Number of processors
        host (str): Remote host (for remote execution)
        username (str): SSH username (for remote execution)
        auth (Dict): Authentication configuration (for remote execution)
        port (int): SSH port (for remote execution)
        env_setup (str): Environment setup commands for SSH sessions
                        (e.g., "source /etc/profile"). Required for making
                        module command available in non-interactive SSH sessions.
    """
    
    def __init__(
        self,
        name: str,
        execution: str = "local",
        scheduler: str = "direct",
        workdir: str = "./",
        use_modules: bool = False,
        modules: Optional[List[str]] = None,
        resources: Optional[Dict] = None,
        prepend: Union[str, List[str], None] = None,
        postpend: Union[str, List[str], None] = None,
        launcher: str = "mpirun -np {nprocs}",
        nprocs: int = 1,
        host: Optional[str] = None,
        username: Optional[str] = None,
        auth: Optional[Dict] = None,
        port: int = 22,
        env_setup: Optional[str] = None,
        **kwargs
    ):
        """
        Initialize a Machine instance.
        
        Parameters:
            name (str): Machine identifier
            execution (str): "local" or "remote"
            scheduler (str): Scheduler type
            workdir (str): Working directory
            use_modules (bool): Enable module loading
            modules (List[str]): Modules to load
            resources (Dict): Scheduler resources
            prepend (str|List[str]): Pre-execution commands
            postpend (str|List[str]): Post-execution commands
            launcher (str): MPI launcher template
            nprocs (int): Number of processors
            host (str): Remote hostname
            username (str): SSH username
            auth (Dict): Authentication config
            port (int): SSH port
            env_setup (str): Environment setup commands for SSH sessions.
                           Example: "source /etc/profile" or "source ~/.bashrc"
                           Useful for making module command available in non-interactive SSH.
            **kwargs: Additional configuration options
        """
        self.name = name
        self.execution = execution
        self.scheduler = scheduler
        self.workdir = workdir
        self.use_modules = use_modules
        self.modules = modules or []
        self.resources = resources or {}
        self.prepend = prepend or []
        self.postpend = postpend or []
        self.launcher = launcher
        self.nprocs = nprocs
        
        # Remote execution attributes
        self.host = host
        self.username = username
        self.auth = auth or {}
        self.port = port
        self.env_setup = env_setup
        
        # Store any additional attributes
        self._extra = kwargs
        
        # Validate configuration
        self._validate()
    
    def _validate(self):
        """Validate machine configuration."""
        if self.execution not in ["local", "remote"]:
            raise ValueError(f"Invalid execution mode: {self.execution}")
        
        if self.execution == "remote":
            if not self.host:
                raise ValueError("Remote execution requires 'host' parameter")
            if not self.username:
                raise ValueError("Remote execution requires 'username' parameter")
            
            # Validate authentication
            auth_method = self.auth.get("method", "key")
            if auth_method != "key":
                raise ValueError(f"Only key-based authentication is supported, got: {auth_method}")
    
    @property
    def is_remote(self) -> bool:
        """Check if this is a remote machine."""
        return self.execution == "remote"
    
    @property
    def is_local(self) -> bool:
        """Check if this is a local machine."""
        return self.execution == "local"
    
    @classmethod
    def from_dict(cls, name: str, config: Dict) -> "Machine":
        """
        Create a Machine instance from a configuration dictionary.
        
        Parameters:
            name (str): Machine name (used if 'name' not in config)
            config (Dict): Configuration dictionary
            
        Returns:
            Machine: New Machine instance
        """
        # If config has a name, use it, otherwise use the provided name
        config = config.copy()  # Don't modify the original
        if 'name' not in config:
            config['name'] = name
        return cls(**config)
    
    @classmethod
    def from_file(cls, filepath: str) -> "Machine":
        """
        Load a Machine from a JSON file.
        
        Parameters:
            filepath (str): Path to JSON configuration file
            
        Returns:
            Machine: New Machine instance
            
        Raises:
            FileNotFoundError: If file doesn't exist
            json.JSONDecodeError: If file contains invalid JSON
        """
        filepath = os.path.expanduser(filepath)
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Machine config file not found: {filepath}")
        
        with open(filepath, 'r') as f:
            config = json.load(f)
        
        # Sanitize configuration to remove unwanted quotes
        from xespresso.machines.config.loader import sanitize_machine_config
        config = sanitize_machine_config(config)
        
        # Extract name from filename if not in config
        name = config.get("name", os.path.splitext(os.path.basename(filepath))[0])
        
        return cls.from_dict(name, config)
    
    def to_dict(self) -> Dict:
        """
        Convert Machine to a dictionary (suitable for JSON serialization).
        
        Returns:
            Dict: Machine configuration as dictionary
        """
        config = {
            "name": self.name,  # Include name in serialized format
            "execution": self.execution,
            "scheduler": self.scheduler,
            "workdir": self.workdir,
            "use_modules": self.use_modules,
            "modules": self.modules,
            "resources": self.resources,
            "prepend": self.prepend,
            "postpend": self.postpend,
            "launcher": self.launcher,
            "nprocs": self.nprocs,
        }
        
        if self.is_remote:
            config["host"] = self.host
            config["username"] = self.username
            config["auth"] = self.auth
            config["port"] = self.port
        
        # Add env_setup if set
        if self.env_setup:
            config["env_setup"] = self.env_setup
        
        # Add any extra attributes
        config.update(self._extra)
        
        return config
    
    def to_file(self, filepath: str):
        """
        Save Machine configuration to a JSON file.
        
        Parameters:
            filepath (str): Path to save configuration
        """
        filepath = os.path.expanduser(filepath)
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)
    
    def to_queue(self) -> Dict:
        """
        Convert Machine to a queue dictionary (for scheduler consumption).
        
        This format is what the schedulers expect and is backward compatible
        with the existing code.
        
        The {nprocs} placeholder in the launcher is substituted with the
        actual nprocs value here.
        
        Returns:
            Dict: Queue configuration dictionary
        """
        def normalize_script_block(block):
            """Normalize script blocks to strings."""
            if isinstance(block, list):
                return "\n".join(block)
            return block or ""
        
        # Substitute {nprocs} placeholder in launcher with actual value
        launcher = self.launcher
        if "{nprocs}" in launcher:
            launcher = launcher.replace("{nprocs}", str(self.nprocs))
        
        queue = {
            "execution": self.execution,
            "scheduler": self.scheduler,
            "use_modules": self.use_modules,
            "modules": self.modules,
            "resources": self.resources,
            "prepend": normalize_script_block(self.prepend),
            "postpend": normalize_script_block(self.postpend),
            "launcher": launcher,
            "nprocs": self.nprocs,
        }
        
        # Add env_setup if set
        if self.env_setup:
            queue["env_setup"] = self.env_setup
        
        if self.is_local:
            queue["local_dir"] = self.workdir
        elif self.is_remote:
            queue["remote_host"] = self.host
            queue["remote_user"] = self.username
            queue["remote_auth"] = {
                "method": self.auth.get("method", "key"),
                "ssh_key": self.auth.get("ssh_key", "~/.ssh/id_rsa"),
                "port": self.auth.get("port", self.port)
            }
            # Pass through auto_install_key if specified
            if "auto_install_key" in self.auth:
                queue["remote_auth"]["auto_install_key"] = self.auth["auto_install_key"]
            queue["remote_dir"] = self.workdir
        
        return queue
    
    def __repr__(self) -> str:
        """String representation of Machine."""
        mode = "remote" if self.is_remote else "local"
        location = f"{self.username}@{self.host}" if self.is_remote else self.workdir
        return f"Machine(name='{self.name}', {mode}={location}, scheduler={self.scheduler})"
    
    def __str__(self) -> str:
        """Human-readable string representation."""
        return f"{self.name} ({self.execution}/{self.scheduler})"
