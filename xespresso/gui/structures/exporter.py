"""
Structure export module for xespresso GUI.

This module handles exporting structures to various formats,
following the modular design pattern.
"""

from typing import Optional
import os
import tempfile
from ase import Atoms
from ase import io as ase_io
import logging

logger = logging.getLogger(__name__)


class StructureExporter:
    """
    Structure exporter for GUI.
    
    This class handles exporting structures to various file formats,
    encapsulating export logic in a dedicated module.
    
    Examples:
        >>> exporter = StructureExporter(atoms)
        >>> file_data = exporter.export_to_bytes('cif')
        >>> # Or save to file
        >>> exporter.export_to_file('structure.cif', 'cif')
    """
    
    def __init__(self, atoms: Atoms):
        """
        Initialize structure exporter.
        
        Args:
            atoms: ASE Atoms object to export
        """
        self.atoms = atoms
        logger.info(f"Initialized structure exporter for {atoms.get_chemical_formula()}")
    
    @staticmethod
    def get_supported_formats() -> list:
        """
        Get list of supported export formats.
        
        Returns:
            list: Supported format names
        """
        return ['cif', 'xyz', 'pdb', 'vasp', 'json']
    
    def export_to_file(self, filepath: str, format: Optional[str] = None) -> None:
        """
        Export structure to file.
        
        Args:
            filepath: Path to output file
            format: Optional format (auto-detected from extension if not provided)
        """
        if format is None:
            # Auto-detect from extension
            ext = os.path.splitext(filepath)[1].lower()
            format = ext[1:] if ext else 'cif'
        
        logger.info(f"Exporting structure to {filepath} (format: {format})")
        
        try:
            ase_io.write(filepath, self.atoms, format=format)
            logger.info(f"Successfully exported structure to {filepath}")
        except Exception as e:
            logger.error(f"Error exporting structure: {e}")
            raise ValueError(f"Error exporting structure: {e}")
    
    def export_to_bytes(self, format: str = 'cif') -> bytes:
        """
        Export structure to bytes.
        
        Args:
            format: Output format
            
        Returns:
            bytes: Exported structure data
        """
        logger.info(f"Exporting structure to bytes (format: {format})")
        
        try:
            # Create temporary file
            with tempfile.NamedTemporaryFile(delete=False, suffix=f".{format}") as tmp:
                tmp_path = tmp.name
            
            # Write to temporary file
            ase_io.write(tmp_path, self.atoms, format=format)
            
            # Read bytes
            with open(tmp_path, 'rb') as f:
                file_data = f.read()
            
            # Clean up
            os.unlink(tmp_path)
            
            logger.info(f"Successfully exported structure to bytes")
            return file_data
            
        except Exception as e:
            logger.error(f"Error exporting structure to bytes: {e}")
            if 'tmp_path' in locals() and os.path.exists(tmp_path):
                os.unlink(tmp_path)
            raise ValueError(f"Error exporting structure: {e}")


def export_structure(atoms: Atoms, format: str = 'cif') -> bytes:
    """
    Convenience function to export structure to bytes.
    
    Args:
        atoms: ASE Atoms object
        format: Output format
        
    Returns:
        bytes: Exported structure data
        
    Example:
        >>> # In Streamlit GUI:
        >>> file_data = export_structure(atoms, 'cif')
        >>> st.download_button(
        >>>     label="Download CIF",
        >>>     data=file_data,
        >>>     file_name="structure.cif"
        >>> )
    """
    exporter = StructureExporter(atoms)
    return exporter.export_to_bytes(format)
