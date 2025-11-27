"""
Path validation utilities for the xespresso PyQt GUI.
"""

import os


def validate_path(path, allow_creation=False, base_dir=None):
    """
    Validate and sanitize file paths to prevent path injection.
    
    Args:
        path: Path to validate
        allow_creation: If True, allow non-existent paths (for file creation)
        base_dir: If provided, ensure path is under this base directory
    
    Returns:
        Tuple of (is_valid, normalized_path, error_message)
    """
    if not path:
        return False, None, "Path cannot be empty"
    
    try:
        # Normalize and resolve the path
        normalized = os.path.abspath(os.path.expanduser(path))
        
        # Additional security: Check for null bytes
        if '\0' in normalized:
            return False, None, "Path contains null bytes"
        
        # If a base_dir is provided, ensure path is under it (prevents path traversal)
        if base_dir:
            base_normalized = os.path.abspath(os.path.expanduser(base_dir))
            try:
                # Use commonpath to check if normalized is under base_normalized
                common = os.path.commonpath([base_normalized, normalized])
                if common != base_normalized:
                    return False, None, "Path traversal not allowed"
            except ValueError:
                # commonpath raises ValueError if paths are on different drives (Windows)
                return False, None, "Path traversal not allowed"
        
        # Check if path exists (if required)
        if not allow_creation and not os.path.exists(normalized):
            return False, normalized, f"Path does not exist: {normalized}"
        
        return True, normalized, None
        
    except Exception as e:
        return False, None, f"Invalid path: {str(e)}"


def validate_path_under_base(path, base_dir):
    """
    Validate that a path is under a base directory.
    
    This is a simpler function for validating that a combined path
    (e.g., base_dir + label) doesn't escape the base directory.
    
    Args:
        path: The full path to validate
        base_dir: The base directory the path must be under
    
    Returns:
        Tuple of (is_valid, normalized_path, error_message)
    """
    if not path or not base_dir:
        return False, None, "Path and base directory cannot be empty"
    
    try:
        full_path = os.path.realpath(path)
        base_real = os.path.realpath(base_dir)
        
        # Check for null bytes
        if '\0' in full_path or '\0' in base_real:
            return False, None, "Path contains null bytes"
        
        # Use commonpath for robust path comparison
        try:
            common = os.path.commonpath([base_real, full_path])
            if common != base_real:
                return False, None, "Path traversal detected"
        except ValueError:
            return False, None, "Path traversal detected"
        
        return True, full_path, None
        
    except (OSError, ValueError) as e:
        return False, None, f"Invalid path: {str(e)}"


def safe_path_exists(path):
    """
    Safely check if a path exists after validation.
    
    Args:
        path: Path to check (should be pre-validated with validate_path)
    
    Returns:
        bool: True if path exists, False otherwise
    """
    try:
        return os.path.exists(path)
    except (OSError, ValueError):
        return False


def safe_makedirs(path):
    """
    Safely create directories after validation.
    
    Args:
        path: Directory path to create (should be pre-validated with validate_path)
    
    Raises:
        OSError: If directory creation fails
    """
    try:
        os.makedirs(path, exist_ok=True)
    except OSError as e:
        raise OSError(f"Failed to create directory: {e}") from e
