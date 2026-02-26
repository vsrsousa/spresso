"""
Formatting utilities for XESPRESSO Lab.
"""

from datetime import datetime
from typing import Optional


def format_energy(energy: Optional[float]) -> str:
    """Format energy value in eV."""
    if energy is None:
        return 'N/A'
    return f'{energy:.6f} eV'


def format_energy_ry(energy: Optional[float]) -> str:
    """Format energy value in Ry (1 Ry = 13.606 eV)."""
    if energy is None:
        return 'N/A'
    ry = energy / 13.606
    return f'{ry:.6f} Ry'


def format_timestamp(ts: Optional[str]) -> str:
    """Format timestamp."""
    if not ts:
        return 'N/A'
    try:
        dt = datetime.fromisoformat(ts)
        return dt.strftime('%Y-%m-%d %H:%M:%S')
    except:
        return str(ts)


def format_wall_time(seconds: Optional[float]) -> str:
    """Format wall time in human-readable format."""
    if seconds is None:
        return 'N/A'
    
    seconds = float(seconds)
    if seconds < 60:
        return f'{seconds:.1f}s'
    elif seconds < 3600:
        minutes = seconds / 60
        return f'{minutes:.1f}m'
    else:
        hours = seconds / 3600
        return f'{hours:.1f}h'


def format_hash(hash_str: str, length: int = 8) -> str:
    """Format hash to shorter version."""
    if not hash_str:
        return 'N/A'
    return hash_str[:length] + '...'


def format_formula(formula: str) -> str:
    """Format chemical formula nicely."""
    if not formula:
        return 'N/A'
    return formula
