"""
manager.py

Utilities for managing Quantum ESPRESSO code configurations.
"""

import os
import re
import json
import subprocess
from typing import Dict, Optional, List, Union
from pathlib import Path

from .config import Code, CodesConfig

# Common QE executable names
COMMON_QE_CODES = [
    'pw', 'ph', 'pp', 'projwfc', 'dos', 'bands', 
    'neb', 'hp', 'dynmat', 'matdyn', 'q2r',
    'pwcond', 'turbo_lanczos', 'turbo_davidson',
    'pw2wannier90', 'wannier90', 'thermo_pw.x',
]

DEFAULT_CODES_DIR = os.path.expanduser("~/.xespresso/codes")


class CodesManager:
    """
    Manager for Quantum ESPRESSO code configurations.
    
    This class provides utilities to:
    - Detect available QE codes on a machine
    - Create and manage code configurations
    - Load and save code configurations
    """
    
    def __init__(self, config: Optional[CodesConfig] = None):
        """
        Initialize CodesManager.
        
        Args:
            config: Optional CodesConfig object
        """
        self.config = config
    
    @classmethod
    def detect_codes(cls, 
                     search_paths: Optional[List[str]] = None,
                     qe_prefix: Optional[str] = None,
                     modules: Optional[List[str]] = None,
                     ssh_connection: Optional[Dict] = None,
                     use_modules: bool = True,
                     env_setup: Optional[str] = None) -> Dict[str, str]:
        """
        Detect available QE codes on the system.
        
        Args:
            search_paths: List of paths to search for executables
            qe_prefix: QE installation prefix (e.g., '/opt/qe-7.2/bin')
            modules: List of modules to load before detection
            ssh_connection: SSH connection info for remote detection.
                           Dict with keys: 'host', 'username', 'port' (default: 22)
                           Example: {'host': 'cluster.edu', 'username': 'user', 'port': 22}
            use_modules: Whether to use module command (default: True, but will
                        auto-detect if modules are available)
            env_setup: Shell commands to set up environment before detection.
                      Useful for sourcing profile files to make module command available.
                      Example: "source /etc/profile" or "source ~/.bashrc"
        
        Returns:
            Dictionary mapping code name to executable path
        """
        detected_codes = {}
        
        # Build search paths
        if qe_prefix:
            if not search_paths:
                search_paths = []
            search_paths.insert(0, qe_prefix)
        
        if not search_paths:
            # Use PATH environment variable
            path_env = os.environ.get('PATH', '')
            search_paths = path_env.split(':')
        
        # Build module load command if needed
        module_cmd = ""
        if modules and use_modules:
            # Check if module command is available
            has_modules = cls._check_module_available(ssh_connection, env_setup)
            if has_modules:
                module_cmd = " && ".join([f"module load {mod}" for mod in modules])
                module_cmd += " && "
            else:
                print("⚠️  'module' command not available, skipping module loads")
        
        # Check for each common QE code
        for code_name in COMMON_QE_CODES:
            executable = f"{code_name}.x"
            
            if ssh_connection:
                # Remote detection
                found_path = cls._detect_remote(executable, search_paths, 
                                                module_cmd, ssh_connection, env_setup)
            else:
                # Local detection
                found_path = cls._detect_local(executable, search_paths, module_cmd, env_setup)
            
            if found_path:
                detected_codes[code_name] = found_path
        
        return detected_codes
    
    @staticmethod
    def _check_module_available(ssh_connection: Optional[Dict] = None,
                                env_setup: Optional[str] = None) -> bool:
        """
        Check if the 'module' command is available on the system.
        
        Args:
            ssh_connection: SSH connection info for remote check
            env_setup: Shell commands to set up environment (e.g., "source /etc/profile")
        
        Returns:
            True if module command is available, False otherwise
        """
        try:
            # Build environment setup prefix
            env_prefix = f"{env_setup} && " if env_setup else ""
            
            if ssh_connection:
                host = ssh_connection.get('host')
                username = ssh_connection.get('username', os.environ.get('USER'))
                port = ssh_connection.get('port', 22)
                cmd = f"ssh -p {port} {username}@{host} '{env_prefix}command -v module > /dev/null 2>&1'"
            else:
                cmd = f"{env_prefix}command -v module > /dev/null 2>&1"
            
            result = subprocess.run(cmd, shell=True, timeout=5)
            return result.returncode == 0
        except (subprocess.TimeoutExpired, Exception):
            return False
    
    @staticmethod
    def _detect_local(executable: str, search_paths: List[str], 
                     module_cmd: str = "", env_setup: Optional[str] = None) -> Optional[str]:
        """Detect executable on local machine."""
        # Build environment setup prefix
        env_prefix = f"{env_setup} && " if env_setup else ""
        
        # Try using 'which' command
        try:
            cmd = f"{env_prefix}{module_cmd}which {executable}"
            result = subprocess.run(cmd, shell=True, capture_output=True, 
                                  text=True, timeout=10)
            if result.returncode == 0:
                path = result.stdout.strip()
                if path and os.path.exists(path):
                    return path
        except (subprocess.TimeoutExpired, Exception):
            pass
        
        # Manual search in paths
        for search_path in search_paths:
            full_path = os.path.join(search_path, executable)
            if os.path.exists(full_path) and os.access(full_path, os.X_OK):
                return full_path
        
        return None
    
    @staticmethod
    def _detect_remote(executable: str, search_paths: List[str],
                      module_cmd: str, ssh_connection: Dict,
                      env_setup: Optional[str] = None) -> Optional[str]:
        """Detect executable on remote machine via SSH."""
        try:
            host = ssh_connection.get('host')
            username = ssh_connection.get('username', os.environ.get('USER'))
            port = ssh_connection.get('port', 22)
            
            # Build SSH command with port
            ssh_cmd = f"ssh -p {port} {username}@{host}"
            
            # Build environment setup prefix
            env_prefix = f"{env_setup} && " if env_setup else ""
            
            # Try using 'which' command remotely
            cmd = f"{ssh_cmd} '{env_prefix}{module_cmd}which {executable}'"
            result = subprocess.run(cmd, shell=True, capture_output=True,
                                  text=True, timeout=30)
            if result.returncode == 0:
                path = result.stdout.strip()
                if path:
                    return path
            
            # Try manual search in paths
            for search_path in search_paths:
                full_path = os.path.join(search_path, executable)
                cmd = f"{ssh_cmd} 'test -x {full_path} && echo {full_path}'"
                result = subprocess.run(cmd, shell=True, capture_output=True,
                                      text=True, timeout=30)
                if result.returncode == 0 and result.stdout.strip():
                    return result.stdout.strip()
        
        except (subprocess.TimeoutExpired, Exception):
            pass
        
        return None
    
    @staticmethod
    def detect_common_prefix(code_paths: Dict[str, str]) -> Optional[str]:
        """
        Detect common directory prefix for a set of code paths.
        
        Args:
            code_paths: Dictionary mapping code name to full path
        
        Returns:
            Common directory prefix or None if paths don't share a common directory
        
        Example:
            >>> paths = {
            ...     'pw': '/opt/qe/7.5/bin/pw.x',
            ...     'ph': '/opt/qe/7.5/bin/ph.x'
            ... }
            >>> detect_common_prefix(paths)
            '/opt/qe/7.5/bin'
        """
        if not code_paths or len(code_paths) < 2:
            # If only one code or none, return its directory
            if len(code_paths) == 1:
                single_path = next(iter(code_paths.values()))
                return str(Path(single_path).parent)
            return None
        
        paths = list(code_paths.values())
        
        # Get parent directory for each path
        directories = [str(Path(p).parent) for p in paths]
        
        # Check if all directories are the same
        if len(set(directories)) == 1:
            return directories[0]
        
        return None
    
    @staticmethod
    def detect_qe_version(pw_path: str, ssh_connection: Optional[Dict] = None,
                         env_setup: Optional[str] = None) -> Optional[str]:
        """
        Detect Quantum ESPRESSO version from pw.x executable.
        
        Args:
            pw_path: Path to pw.x executable
            ssh_connection: SSH connection info for remote detection
                           {'host': 'hostname', 'username': 'user', 'port': 22}
            env_setup: Shell commands to set up environment (e.g., "source /etc/profile")
        
        Returns:
            Version string (e.g., '7.2') or None
        """
        try:
            # Build environment setup prefix
            env_prefix = f"{env_setup} && " if env_setup else ""
            
            if ssh_connection:
                host = ssh_connection.get('host')
                username = ssh_connection.get('username', os.environ.get('USER'))
                port = ssh_connection.get('port', 22)
                cmd = f"ssh -p {port} {username}@{host} '{env_prefix}{pw_path} --version 2>&1 | head -5'"
            else:
                cmd = f"{env_prefix}{pw_path} --version 2>&1 | head -5"
            
            result = subprocess.run(cmd, shell=True, capture_output=True,
                                  text=True, timeout=10)
            
            if result.returncode == 0 or result.stderr:
                output = result.stdout + result.stderr
                # Look for version pattern like "7.2" or "6.8"
                match = re.search(r'(\d+\.\d+)', output)
                if match:
                    return match.group(1)
        except (subprocess.TimeoutExpired, Exception):
            pass
        
        return None
    
    @classmethod
    def create_config(cls,
                     machine_name: str,
                     detected_codes: Dict[str, str],
                     qe_version: Optional[str] = None,
                     qe_prefix: Optional[str] = None,
                     label: Optional[str] = None,
                     modules: Optional[List[str]] = None,
                     environment: Optional[Dict[str, str]] = None) -> CodesConfig:
        """
        Create a CodesConfig from detected codes.
        
        Args:
            machine_name: Name of the machine
            detected_codes: Dictionary mapping code name to path
            qe_version: QE version string
            qe_prefix: QE installation prefix (will auto-detect if not provided)
            label: Custom label for this version
            modules: List of modules to load
            environment: Environment variables
        
        Returns:
            CodesConfig object
        """
        # Auto-detect common prefix if not provided and codes share common directory
        if not qe_prefix and detected_codes:
            detected_prefix = cls.detect_common_prefix(detected_codes)
            if detected_prefix:
                qe_prefix = detected_prefix
        
        config = CodesConfig(
            machine_name=machine_name,
            qe_prefix=qe_prefix,
            qe_version=qe_version,
            label=label,
            modules=modules,
            environment=environment
        )
        
        # If version or label is specified, store codes in versions structure
        # Otherwise store in main codes dictionary
        use_versions = qe_version is not None or label is not None
        
        for code_name, code_path in detected_codes.items():
            code = Code(
                name=code_name,
                path=code_path,
                version=qe_version
            )
            if use_versions and qe_version:
                # Add to version-specific structure
                config.add_code(code, version=qe_version)
                # Store label, modules, and qe_prefix for this version
                if not config.versions:
                    config.versions = {}
                if qe_version not in config.versions:
                    config.versions[qe_version] = {"codes": {}}
                if label:
                    config.versions[qe_version]["label"] = label
                if modules:
                    config.versions[qe_version]["modules"] = modules
                if qe_prefix:
                    config.versions[qe_version]["qe_prefix"] = qe_prefix
            else:
                # Add to main codes dictionary (backward compatible)
                config.add_code(code)
        
        return config
    
    @staticmethod
    def save_config(config: CodesConfig, 
                   output_dir: str = DEFAULT_CODES_DIR,
                   filename: Optional[str] = None,
                   overwrite: bool = False,
                   merge: bool = False,
                   interactive: bool = True) -> str:
        """
        Save a CodesConfig to JSON file.
        
        Args:
            config: CodesConfig object to save
            output_dir: Directory to save the file
            filename: Optional filename (default: <machine_name>.json)
            overwrite: If True, overwrites existing file without asking
            merge: If True and file exists, merges with existing config
            interactive: If True, prompts user when file exists (default: True)
        
        Returns:
            Path to saved file
        
        Raises:
            FileExistsError: If file exists and overwrite/merge are False
        """
        os.makedirs(output_dir, exist_ok=True)
        
        if not filename:
            # Always use machine name for filename
            # Multiple versions/labels are stored within the file structure
            filename = f"{config.machine_name}.json"
        
        filepath = os.path.join(output_dir, filename)
        
        # Check if file already exists
        if os.path.exists(filepath) and not overwrite and not merge:
            if interactive:
                print(f"⚠️  Configuration file already exists: {filepath}")
                response = input("Choose: [o]verwrite, [m]erge, [c]ancel (default: cancel): ").strip().lower()
                if response == 'o':
                    overwrite = True
                elif response == 'm':
                    merge = True
                else:
                    print("❌ Cancelled. Configuration not saved.")
                    raise FileExistsError(f"Configuration file already exists: {filepath}")
            else:
                # Non-interactive mode: raise error immediately
                raise FileExistsError(f"Configuration file already exists: {filepath}")
        
        # Handle merge
        if merge and os.path.exists(filepath):
            try:
                existing_config = CodesConfig.from_json(filepath)
                print(f"📝 Merging with existing configuration...")
                
                # Merge main codes (only if new config has codes in main dict)
                for code_name, code in config.codes.items():
                    existing_config.codes[code_name] = code
                
                # Merge versions - this is the critical part
                if config.versions:
                    if not existing_config.versions:
                        existing_config.versions = {}
                    for version, version_config in config.versions.items():
                        if version not in existing_config.versions:
                            # New version - add it entirely
                            existing_config.versions[version] = version_config
                            print(f"   ✅ Added new version: {version}")
                        else:
                            # Version exists - merge codes and update metadata
                            print(f"   📝 Merging into existing version: {version}")
                            if 'codes' in version_config:
                                if 'codes' not in existing_config.versions[version]:
                                    existing_config.versions[version]['codes'] = {}
                                existing_config.versions[version]['codes'].update(version_config['codes'])
                            # Update label if provided
                            if 'label' in version_config:
                                existing_config.versions[version]['label'] = version_config['label']
                            # Update modules if provided
                            if 'modules' in version_config:
                                existing_config.versions[version]['modules'] = version_config['modules']
                
                # Update top-level fields (these apply to the whole machine config)
                # Don't overwrite unless explicitly set in new config
                if config.qe_prefix:
                    existing_config.qe_prefix = config.qe_prefix
                # Note: qe_version and label at top level are for default/current version
                # They don't overwrite when merging versions
                if config.environment:
                    if not existing_config.environment:
                        existing_config.environment = {}
                    existing_config.environment.update(config.environment)
                
                config = existing_config
                print(f"✅ Merged configurations")
            except Exception as e:
                print(f"⚠️  Failed to merge configurations: {e}")
                print("   Falling back to overwrite")
                import traceback
                traceback.print_exc()
        
        config.to_json(filepath)
        
        return filepath
    
    @staticmethod
    def load_config(machine_name: str,
                   codes_dir: str = DEFAULT_CODES_DIR) -> Optional[CodesConfig]:
        """
        Load a CodesConfig from JSON file.
        
        Args:
            machine_name: Name of the machine
            codes_dir: Directory containing codes configurations
        
        Returns:
            CodesConfig object or None if not found
        """
        filepath = os.path.join(codes_dir, f"{machine_name}.json")
        
        if not os.path.exists(filepath):
            return None
        
        try:
            return CodesConfig.from_json(filepath)
        except Exception as e:
            print(f"Error loading codes config: {e}")
            return None


def detect_qe_codes(machine_name: str = "local",
                   qe_prefix: Optional[str] = None,
                   search_paths: Optional[List[str]] = None,
                   modules: Optional[List[str]] = None,
                   ssh_connection: Optional[Dict] = None,
                   env_setup: Optional[str] = None,
                   auto_load_machine: bool = True,
                   qe_version: Optional[str] = None,
                   label: Optional[str] = None) -> CodesConfig:
    """
    Convenience function to detect and create a codes configuration.
    
    This function can automatically load machine configuration from machines.json
    or individual machine files to extract SSH connection details and modules.
    
    Args:
        machine_name: Name of the machine
        qe_prefix: QE installation prefix (e.g., '/opt/qe-7.2/bin')
        search_paths: List of paths to search for executables
        modules: List of modules to load before detection
        ssh_connection: SSH connection info for remote detection
                       {'host': 'hostname', 'username': 'user', 'port': 22}
        env_setup: Shell commands to set up environment before detection.
                  Example: "source /etc/profile" or "source ~/.bashrc"
        auto_load_machine: If True, attempts to load machine config automatically
        qe_version: Optional version string to use instead of auto-detection.
                   Useful when auto-detection returns compiler version instead.
                   Example: "7.2", "7.1", "6.8"
        label: Optional custom label for this version.
              Example: "qe-7.2-production", "qe-dev"
    
    Returns:
        CodesConfig object with detected codes
    """
    print(f"🔍 Detecting Quantum ESPRESSO codes on '{machine_name}'...")
    
    # Try to load machine configuration if auto_load_machine is True
    if auto_load_machine and ssh_connection is None:
        try:
            from xespresso.machines.config.loader import load_machine, DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR
            from xespresso.machines.machine import Machine
            
            # Try to load machine config
            machine_obj = load_machine(
                config_path=DEFAULT_CONFIG_PATH,
                machine_name=machine_name,
                machines_dir=DEFAULT_MACHINES_DIR,
                return_object=True
            )
            
            if machine_obj and isinstance(machine_obj, Machine):
                print(f"📋 Loaded machine configuration for '{machine_name}'")
                
                # Extract SSH connection info for remote machines
                if machine_obj.is_remote:
                    ssh_connection = {
                        'host': machine_obj.host,
                        'username': machine_obj.username,
                        'port': machine_obj.port if hasattr(machine_obj, 'port') and machine_obj.port else 22
                    }
                    print(f"🌐 Using remote connection: {machine_obj.username}@{machine_obj.host}:{ssh_connection['port']}")
                
                # Use modules from machine config if not explicitly provided
                if modules is None and machine_obj.use_modules and machine_obj.modules:
                    modules = machine_obj.modules
                    print(f"📦 Using modules from machine config: {', '.join(modules)}")
                
                # Use env_setup from machine config if not explicitly provided
                # Check for env_setup in machine's prepend or as a separate field
                if env_setup is None:
                    if hasattr(machine_obj, 'env_setup') and machine_obj.env_setup:
                        env_setup = machine_obj.env_setup
                        print(f"🔧 Using env_setup from machine config")
                    elif machine_obj.prepend:
                        # Check if prepend contains environment setup commands
                        prepend_str = machine_obj.prepend if isinstance(machine_obj.prepend, str) else '\n'.join(machine_obj.prepend)
                        if 'source' in prepend_str or 'module' in prepend_str:
                            env_setup = prepend_str
                            print(f"🔧 Using prepend as env_setup from machine config")
        except ImportError:
            pass  # Machine loader not available
        except Exception as e:
            # Log but don't fail - continue with manual parameters
            import logging
            logging.debug(f"Could not auto-load machine config: {e}")
    
    detected_codes = CodesManager.detect_codes(
        search_paths=search_paths,
        qe_prefix=qe_prefix,
        modules=modules,
        ssh_connection=ssh_connection,
        env_setup=env_setup
    )
    
    if not detected_codes:
        print("⚠️  No Quantum ESPRESSO codes detected")
        return CodesConfig(
            machine_name=machine_name,
            qe_version=qe_version,
            label=label
        )
    
    print(f"✅ Found {len(detected_codes)} codes: {', '.join(detected_codes.keys())}")
    
    # Try to detect QE version from pw.x if not explicitly provided
    if qe_version is None:
        if 'pw' in detected_codes:
            qe_version = CodesManager.detect_qe_version(
                detected_codes['pw'], 
                ssh_connection=ssh_connection,
                env_setup=env_setup
            )
            if qe_version:
                print(f"📦 Detected Quantum ESPRESSO version: {qe_version}")
    else:
        print(f"📦 Using user-specified Quantum ESPRESSO version: {qe_version}")
    
    config = CodesManager.create_config(
        machine_name=machine_name,
        detected_codes=detected_codes,
        qe_version=qe_version,
        qe_prefix=qe_prefix,
        label=label,
        modules=modules
    )
    
    return config


def create_codes_config(machine_name: str = "local",
                       qe_prefix: Optional[str] = None,
                       search_paths: Optional[List[str]] = None,
                       modules: Optional[List[str]] = None,
                       ssh_connection: Optional[Dict] = None,
                       env_setup: Optional[str] = None,
                       save: bool = True,
                       output_dir: str = DEFAULT_CODES_DIR,
                       overwrite: bool = False,
                       merge: bool = True,
                       auto_load_machine: bool = True,
                       qe_version: Optional[str] = None,
                       label: Optional[str] = None) -> CodesConfig:
    """
    Create a codes configuration (with optional auto-save).
    
    Args:
        machine_name: Name of the machine
        qe_prefix: QE installation prefix
        search_paths: List of paths to search for executables
        modules: List of modules to load before detection
        ssh_connection: SSH connection info for remote detection
        env_setup: Shell commands to set up environment before detection
        save: Whether to save the configuration to file
        output_dir: Directory to save the configuration
        overwrite: If True, overwrites existing file without asking
        merge: If True, merges with existing config (default: True)
        auto_load_machine: If True, attempts to load machine config automatically
        qe_version: Optional version string to use instead of auto-detection.
                   Example: "7.2", "7.1", "6.8"
        label: Optional custom label for this version.
              Example: "qe-7.2-production", "qe-dev"
    
    Returns:
        CodesConfig object
    """
    config = detect_qe_codes(
        machine_name=machine_name,
        qe_prefix=qe_prefix,
        search_paths=search_paths,
        modules=modules,
        ssh_connection=ssh_connection,
        env_setup=env_setup,
        auto_load_machine=auto_load_machine,
        qe_version=qe_version,
        label=label
    )
    
    if save and config.codes:
        try:
            filepath = CodesManager.save_config(
                config, 
                output_dir=output_dir,
                overwrite=overwrite,
                merge=merge
            )
            print(f"💾 Configuration saved to: {filepath}")
        except FileExistsError as e:
            print(f"⚠️  {e}")
    
    return config


def load_codes_config(machine_name: str,
                     codes_dir: str = DEFAULT_CODES_DIR,
                     version: Optional[str] = None,
                     verbose: bool = False) -> Optional[CodesConfig]:
    """
    Load a codes configuration from file.
    
    Args:
        machine_name: Name of the machine
        codes_dir: Directory containing codes configurations
        version: Optional QE version to use. If None, uses default or main codes.
        verbose: If True, prints status messages to stdout. Default is False.
    
    Returns:
        CodesConfig object or None if not found
    """
    config = CodesManager.load_config(machine_name, codes_dir)
    
    if config:
        versions = config.list_versions()
        if versions:
            if verbose:
                print(f"✅ Loaded codes configuration for '{machine_name}'")
                print(f"   Available versions: {', '.join(versions)}")
            
            # Show codes for specified version or default
            target_version = version or config.qe_version or versions[0]
            codes = config.list_codes(version=target_version)
            if codes and verbose:
                print(f"   Codes in version {target_version}: {', '.join(codes)}")
        else:
            if verbose:
                print(f"✅ Loaded codes configuration for '{machine_name}'")
                print(f"   Available codes: {', '.join(config.list_codes())}")
    else:
        if verbose:
            print(f"⚠️  No codes configuration found for '{machine_name}'")
    
    return config


def add_version_to_config(machine_name: str,
                         version: str,
                         qe_prefix: Optional[str] = None,
                         search_paths: Optional[List[str]] = None,
                         modules: Optional[List[str]] = None,
                         ssh_connection: Optional[Dict] = None,
                         env_setup: Optional[str] = None,
                         codes_dir: str = DEFAULT_CODES_DIR) -> Optional[CodesConfig]:
    """
    Add a new QE version to an existing codes configuration.
    
    This allows managing multiple QE versions on the same machine. Each version
    can have its own set of executables, modules, and environment settings.
    
    Args:
        machine_name: Name of the machine
        version: QE version string (e.g., '7.2', '6.8')
        qe_prefix: QE installation prefix for this version
        search_paths: Paths to search for executables
        modules: Modules to load for this version
        ssh_connection: SSH connection info for remote detection
        env_setup: Shell commands to set up environment before detection
        codes_dir: Directory containing codes configurations
    
    Returns:
        Updated CodesConfig object or None if machine config doesn't exist
    
    Example:
        # Add QE 7.2 to existing config
        config = add_version_to_config(
            machine_name="cluster1",
            version="7.2",
            qe_prefix="/opt/qe-7.2/bin",
            modules=["quantum-espresso/7.2"]
        )
        
        # Add QE 6.8
        config = add_version_to_config(
            machine_name="cluster1",
            version="6.8",
            qe_prefix="/opt/qe-6.8/bin",
            modules=["quantum-espresso/6.8"]
        )
    """
    # Load existing config
    config = CodesManager.load_config(machine_name, codes_dir)
    
    if not config:
        print(f"⚠️  No codes configuration found for '{machine_name}'")
        print(f"   Creating new configuration...")
        config = CodesConfig(
            machine_name=machine_name,
            qe_version=version  # Set as default version
        )
    
    print(f"🔍 Detecting Quantum ESPRESSO {version} codes on '{machine_name}'...")
    
    # Detect codes for this version
    detected_codes = CodesManager.detect_codes(
        search_paths=search_paths,
        qe_prefix=qe_prefix,
        modules=modules,
        ssh_connection=ssh_connection,
        env_setup=env_setup
    )
    
    if not detected_codes:
        print(f"⚠️  No codes detected for version {version}")
        return config
    
    print(f"✅ Found {len(detected_codes)} codes: {', '.join(detected_codes.keys())}")
    
    # Initialize versions dict if needed
    if not config.versions:
        config.versions = {}
    
    # Add version configuration
    config.versions[version] = {
        'qe_prefix': qe_prefix,
        'modules': modules or [],
        'codes': {}
    }
    
    # Add detected codes to this version
    for code_name, code_path in detected_codes.items():
        code = Code(
            name=code_name,
            path=code_path,
            version=version
        )
        config.add_code(code, version=version)
    
    # Set as default version if not set
    if not config.qe_version:
        config.qe_version = version
        print(f"📌 Set {version} as default version")
    
    # Save updated configuration (always overwrite when explicitly adding versions)
    filepath = CodesManager.save_config(config, output_dir=codes_dir, overwrite=True)
    print(f"💾 Configuration saved to: {filepath}")
    
    return config
