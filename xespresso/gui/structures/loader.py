"""
Structure loading module for xespresso GUI.

This module handles loading structures from various sources (files, uploads, etc.)
following the modular design pattern.
"""

from typing import Optional, List, Tuple, Union
from pathlib import Path
import os
import tempfile
from ase import Atoms
from ase import io as ase_io
from xespresso.gui.structures.base import BaseStructureHandler
import logging

logger = logging.getLogger(__name__)


class StructureLoader(BaseStructureHandler):
    """
    Structure loader for GUI.
    
    This class handles loading structures from files, uploaded data,
    and directory browsing, encapsulating all the loading logic
    in a dedicated module.
    
    Examples:
        >>> # Load from file
        >>> loader = StructureLoader()
        >>> atoms = loader.load_from_file('structure.cif')
        >>> 
        >>> # Load from uploaded bytes
        >>> loader = StructureLoader()
        >>> atoms = loader.load_from_upload(file_bytes, 'structure.cif')
    """
    
    @staticmethod
    def get_supported_extensions() -> List[str]:
        """
        Get list of supported structure file extensions.
        
        Returns:
            List[str]: Supported extensions
        """
        return ['.cif', '.xyz', '.pdb', '.vasp', '.poscar', '.traj', '.json']
    
    def load_from_file(self, filepath: str) -> Atoms:
        """
        Load structure from file.
        
        Args:
            filepath: Path to structure file
            
        Returns:
            Atoms: Loaded ASE Atoms object
            
        Raises:
            FileNotFoundError: If file doesn't exist
            ValueError: If file format is not supported
        """
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
        
        # Validate path (prevent path traversal)
        try:
            filepath = os.path.realpath(filepath)
        except (OSError, ValueError) as e:
            raise ValueError(f"Invalid file path: {e}")
        
        logger.info(f"Loading structure from {filepath}")
        
        try:
            self.atoms = ase_io.read(filepath)
            self.source_path = filepath
            self.format = self._detect_format(filepath)
            
            logger.info(f"Successfully loaded {len(self.atoms)} atoms from {filepath}")
            return self.atoms
            
        except Exception as e:
            logger.error(f"Error loading structure from {filepath}: {e}")
            raise ValueError(f"Error loading structure: {e}")
    
    def load_from_upload(self, file_bytes: bytes, filename: str) -> Atoms:
        """
        Load structure from uploaded file bytes.
        
        Args:
            file_bytes: Uploaded file bytes
            filename: Original filename
            
        Returns:
            Atoms: Loaded ASE Atoms object
        """
        logger.info(f"Loading structure from upload: {filename}")
        
        # Create temporary file
        suffix = f".{filename.split('.')[-1]}" if '.' in filename else '.tmp'
        
        try:
            with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
                tmp.write(file_bytes)
                tmp_path = tmp.name
            
            # Load from temporary file
            self.atoms = ase_io.read(tmp_path)
            self.source_path = filename
            self.format = self._detect_format(filename)
            
            # Clean up temporary file
            os.unlink(tmp_path)
            
            logger.info(f"Successfully loaded {len(self.atoms)} atoms from upload")
            return self.atoms
            
        except Exception as e:
            logger.error(f"Error loading structure from upload: {e}")
            if 'tmp_path' in locals() and os.path.exists(tmp_path):
                os.unlink(tmp_path)
            raise ValueError(f"Error loading structure: {e}")
    
    @staticmethod
    def find_structure_files(
        directory: str,
        max_depth: int = 3,
        validate_safety: bool = True
    ) -> List[str]:
        """
        Find structure files in directory.
        
        Args:
            directory: Directory to search
            max_depth: Maximum recursion depth
            validate_safety: Whether to validate paths for security
            
        Returns:
            List[str]: List of structure file paths
        """
        if not os.path.exists(directory) or not os.path.isdir(directory):
            logger.warning(f"Invalid directory: {directory}")
            return []
        
        # Validate and normalize directory
        if validate_safety:
            try:
                directory = os.path.realpath(directory)
                # Check if directory is under safe base directories
                safe_bases = [
                    os.path.realpath(os.path.expanduser("~")),
                    os.path.realpath("/tmp")
                ]
                is_safe = any(directory.startswith(base) for base in safe_bases)
                
                if not is_safe:
                    logger.warning(f"Directory not in safe location: {directory}")
                    return []
            except (OSError, ValueError) as e:
                logger.error(f"Error validating directory: {e}")
                return []
        
        structure_extensions = StructureLoader.get_supported_extensions()
        structure_files = []
        
        logger.info(f"Searching for structure files in {directory}")
        
        for root, dirs, files in os.walk(directory):
            # Validate path safety
            if validate_safety:
                try:
                    if not os.path.realpath(root).startswith(directory):
                        continue
                except (OSError, ValueError):
                    continue
            
            # Check depth
            depth = root[len(directory):].count(os.sep)
            if depth >= max_depth:
                continue
            
            # Find structure files
            for f in files:
                if any(f.lower().endswith(ext) for ext in structure_extensions):
                    file_path = os.path.join(root, f)
                    
                    # Validate constructed path
                    if validate_safety:
                        try:
                            real_path = os.path.realpath(file_path)
                            if not real_path.startswith(directory):
                                continue
                        except (OSError, ValueError):
                            continue
                    
                    structure_files.append(file_path)
        
        logger.info(f"Found {len(structure_files)} structure files")
        return structure_files


def load_structure_from_file(filepath: str) -> Tuple[Atoms, StructureLoader]:
    """
    Convenience function to load structure from file.
    
    Args:
        filepath: Path to structure file
        
    Returns:
        tuple: (atoms, loader) where loader contains metadata
        
    Example:
        >>> atoms, loader = load_structure_from_file('structure.cif')
        >>> print(loader.get_info())
    """
    loader = StructureLoader()
    atoms = loader.load_from_file(filepath)
    return atoms, loader


def load_structure_from_upload(
    file_bytes: bytes,
    filename: str
) -> Tuple[Atoms, StructureLoader]:
    """
    Convenience function to load structure from upload.
    
    Args:
        file_bytes: Uploaded file bytes
        filename: Original filename
        
    Returns:
        tuple: (atoms, loader) where loader contains metadata
        
    Example:
        >>> # In Streamlit GUI:
        >>> uploaded_file = st.file_uploader("Upload structure")
        >>> if uploaded_file:
        >>>     atoms, loader = load_structure_from_upload(
        >>>         uploaded_file.getvalue(),
        >>>         uploaded_file.name
        >>>     )
    """
    loader = StructureLoader()
    atoms = loader.load_from_upload(file_bytes, filename)
    return atoms, loader
