"""
Structure visualization utilities for the Qt GUI.

This module provides reusable visualization components for atomic structures.
"""

import os

try:
    import matplotlib
    if matplotlib.get_backend() != 'QtAgg':
        matplotlib.use('QtAgg')
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure
    import numpy as np
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

# Element color map for visualization
ELEMENT_COLORS = {
    'H': 'white', 'C': 'gray', 'N': 'blue', 'O': 'red',
    'F': 'green', 'P': 'orange', 'S': 'yellow',
    'Cl': 'green', 'Fe': 'brown', 'Cu': 'brown',
    'Al': 'silver', 'Si': 'pink', 'Pt': 'silver'
}


def create_structure_figure(atoms, figure=None):
    """
    Create a 3D visualization of an atomic structure.
    
    Args:
        atoms: ASE Atoms object
        figure: Optional matplotlib Figure to use (creates new if None)
        
    Returns:
        Figure object with the structure visualization
    """
    if not MATPLOTLIB_AVAILABLE:
        return None
    
    if figure is None:
        figure = Figure(figsize=(8, 6))
    else:
        figure.clear()
    
    ax = figure.add_subplot(111, projection='3d')
    
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    
    colors = [ELEMENT_COLORS.get(s, 'purple') for s in symbols]
    
    # Draw bonds between atoms
    _draw_bonds(ax, atoms)
    
    # Draw atoms (after bonds so atoms are on top)
    ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
               c=colors, s=100, edgecolors='black')
    
    # Draw cell if present
    if atoms.cell is not None and atoms.pbc.any():
        _draw_cell(ax, atoms.cell.array)
    
    # Add legend with element labels in corner
    _add_element_legend(ax, symbols, colors)
    
    # Hide axes
    ax.set_axis_off()
    
    figure.tight_layout()
    return figure


def _draw_bonds(ax, atoms):
    """Draw bonds between atoms based on natural cutoff distances."""
    try:
        from ase.neighborlist import natural_cutoffs, NeighborList
    except ImportError:
        # If neighbor list not available, skip bond drawing
        return
    
    # Use natural cutoffs with a multiplier for bond detection
    cutoffs = natural_cutoffs(atoms, mult=1.1)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=False)
    nl.update(atoms)
    
    positions = atoms.get_positions()
    
    # Draw bonds
    for i in range(len(atoms)):
        indices, offsets = nl.get_neighbors(i)
        for j, offset in zip(indices, offsets):
            # Calculate position of neighbor including periodic offset
            pos_j = positions[j] + np.dot(offset, atoms.cell.array if atoms.cell is not None else np.zeros((3, 3)))
            
            # Draw bond as a line
            ax.plot([positions[i][0], pos_j[0]], 
                   [positions[i][1], pos_j[1]], 
                   [positions[i][2], pos_j[2]], 
                   'gray', linewidth=1.5, alpha=0.6)


def _add_element_legend(ax, symbols, colors):
    """Add a legend in the corner showing unique elements and their colors."""
    from matplotlib.patches import Circle
    
    # Get unique elements and their colors
    unique_elements = []
    unique_colors = []
    seen = set()
    
    for sym, color in zip(symbols, colors):
        if sym not in seen:
            unique_elements.append(sym)
            unique_colors.append(color)
            seen.add(sym)
    
    # Create legend handles
    legend_elements = []
    for elem, color in zip(unique_elements, unique_colors):
        legend_elements.append(Circle((0, 0), 1, facecolor=color, 
                                     edgecolor='black', label=elem))
    
    # Add legend in upper right corner
    ax.legend(handles=legend_elements, loc='upper right', 
             framealpha=0.9, fontsize=10)


def _draw_cell(ax, cell):
    """Draw the unit cell edges on a 3D axis."""
    edges = [
        ([0, 0, 0], [1, 0, 0]), ([0, 0, 0], [0, 1, 0]), ([0, 0, 0], [0, 0, 1]),
        ([1, 0, 0], [1, 1, 0]), ([1, 0, 0], [1, 0, 1]),
        ([0, 1, 0], [1, 1, 0]), ([0, 1, 0], [0, 1, 1]),
        ([0, 0, 1], [1, 0, 1]), ([0, 0, 1], [0, 1, 1]),
        ([1, 1, 0], [1, 1, 1]), ([1, 0, 1], [1, 1, 1]), ([0, 1, 1], [1, 1, 1])
    ]
    
    for start, end in edges:
        start_pt = np.dot(start, cell)
        end_pt = np.dot(end, cell)
        ax.plot([start_pt[0], end_pt[0]], [start_pt[1], end_pt[1]], 
                [start_pt[2], end_pt[2]], 'k-', linewidth=0.5)


def get_structure_info_text(atoms):
    """
    Generate a text description of an atomic structure.
    
    Args:
        atoms: ASE Atoms object
        
    Returns:
        String with structure information
    """
    info_lines = []
    info_lines.append(f"Chemical Formula: {atoms.get_chemical_formula()}")
    info_lines.append(f"Number of Atoms: {len(atoms)}")
    
    symbols = atoms.get_chemical_symbols()
    unique_elements = list(set(symbols))
    info_lines.append(f"Elements: {', '.join(sorted(unique_elements))}")
    
    if atoms.cell is not None and atoms.pbc.any():
        info_lines.append(f"\nCell Volume: {atoms.get_volume():.2f} Å³")
        pbc_str = "".join(["T" if p else "F" for p in atoms.pbc])
        info_lines.append(f"PBC: {pbc_str}")
        
        cell_params = atoms.cell.cellpar()
        info_lines.append(f"\nCell Parameters:")
        info_lines.append(f"  a = {cell_params[0]:.3f} Å")
        info_lines.append(f"  b = {cell_params[1]:.3f} Å")
        info_lines.append(f"  c = {cell_params[2]:.3f} Å")
        info_lines.append(f"  α = {cell_params[3]:.2f}°")
        info_lines.append(f"  β = {cell_params[4]:.2f}°")
        info_lines.append(f"  γ = {cell_params[5]:.2f}°")
    
    return "\n".join(info_lines)
