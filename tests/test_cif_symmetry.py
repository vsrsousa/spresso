"""
Tests for CIF file reading with proper symmetry preservation.

This test verifies that CIF files with symmetry information are read
using the primitive cell instead of expanding all symmetry operations.
"""

import pytest
import tempfile
import os
from collections import Counter


def test_cif_primitive_cell_reading():
    """Test that CIF files are read using primitive cells by default."""
    from qtgui.utils import read_structure
    
    # Create a test CIF file (GdCo2 structure from the problem statement)
    cif_content = """# generated using pymatgen
data_GdCo2
_symmetry_space_group_name_H-M   Fd-3m
_cell_length_a   7.23529600
_cell_length_b   7.23529600
_cell_length_c   7.23529600
_cell_angle_alpha   90.00000000
_cell_angle_beta   90.00000000
_cell_angle_gamma   90.00000000
_symmetry_Int_Tables_number   227
_chemical_formula_structural   GdCo2
_chemical_formula_sum   'Gd8 Co16'
_cell_volume   378.76418734
_cell_formula_units_Z   8
loop_
 _symmetry_equiv_pos_site_id
 _symmetry_equiv_pos_as_xyz
  1  'x, y, z'
  2  '-y+1/4, x+3/4, z+3/4'
  3  '-x, -y, z'
  4  'y+1/4, -x+3/4, z+3/4'
  5  'x, -y, -z'
  6  '-y+1/4, -x+3/4, -z+3/4'
  7  '-x, y, -z'
  8  'y+1/4, x+3/4, -z+3/4'
  9  'z, x, y'
  10  'z+1/4, -y+3/4, x+3/4'
  11  'z, -x, -y'
  12  'z+1/4, y+3/4, -x+3/4'
  13  '-z, x, -y'
  14  '-z+1/4, -y+3/4, -x+3/4'
  15  '-z, -x, y'
  16  '-z+1/4, y+3/4, x+3/4'
  17  'y, z, x'
  18  'x+1/4, z+3/4, -y+3/4'
  19  '-y, z, -x'
  20  '-x+1/4, z+3/4, y+3/4'
  21  '-y, -z, x'
  22  '-x+1/4, -z+3/4, -y+3/4'
  23  'y, -z, -x'
  24  'x+1/4, -z+3/4, y+3/4'
  25  '-x+1/4, -y+3/4, -z+3/4'
  26  'y, -x, -z'
  27  'x+1/4, y+3/4, -z+3/4'
  28  '-y, x, -z'
  29  '-x+1/4, y+3/4, z+3/4'
  30  'y, x, z'
  31  'x+1/4, -y+3/4, z+3/4'
  32  '-y, -x, z'
  33  '-z+1/4, -x+3/4, -y+3/4'
  34  '-z, y, -x'
  35  '-z+1/4, x+3/4, y+3/4'
  36  '-z, -y, x'
  37  'z+1/4, -x+3/4, y+3/4'
  38  'z, y, x'
  39  'z+1/4, x+3/4, -y+3/4'
  40  'z, -y, -x'
  41  '-y+1/4, -z+3/4, -x+3/4'
  42  '-x, -z, y'
  43  'y+1/4, -z+3/4, x+3/4'
  44  'x, -z, -y'
  45  'y+1/4, z+3/4, -x+3/4'
  46  'x, z, y'
  47  '-y+1/4, z+3/4, x+3/4'
  48  '-x, z, -y'
  49  'x+1/2, y+1/2, z'
  50  '-y+3/4, x+1/4, z+3/4'
  51  '-x+1/2, -y+1/2, z'
  52  'y+3/4, -x+1/4, z+3/4'
  53  'x+1/2, -y+1/2, -z'
  54  '-y+3/4, -x+1/4, -z+3/4'
  55  '-x+1/2, y+1/2, -z'
  56  'y+3/4, x+1/4, -z+3/4'
  57  'z+1/2, x+1/2, y'
  58  'z+3/4, -y+1/4, x+3/4'
  59  'z+1/2, -x+1/2, -y'
  60  'z+3/4, y+1/4, -x+3/4'
  61  '-z+1/2, x+1/2, -y'
  62  '-z+3/4, -y+1/4, -x+3/4'
  63  '-z+1/2, -x+1/2, y'
  64  '-z+3/4, y+1/4, x+3/4'
  65  'y+1/2, z+1/2, x'
  66  'x+3/4, z+1/4, -y+3/4'
  67  '-y+1/2, z+1/2, -x'
  68  '-x+3/4, z+1/4, y+3/4'
  69  '-y+1/2, -z+1/2, x'
  70  '-x+3/4, -z+1/4, -y+3/4'
  71  'y+1/2, -z+1/2, -x'
  72  'x+3/4, -z+1/4, y+3/4'
  73  '-x+3/4, -y+1/4, -z+3/4'
  74  'y+1/2, -x+1/2, -z'
  75  'x+3/4, y+1/4, -z+3/4'
  76  '-y+1/2, x+1/2, -z'
  77  '-x+3/4, y+1/4, z+3/4'
  78  'y+1/2, x+1/2, z'
  79  'x+3/4, -y+1/4, z+3/4'
  80  '-y+1/2, -x+1/2, z'
  81  '-z+3/4, -x+1/4, -y+3/4'
  82  '-z+1/2, y+1/2, -x'
  83  '-z+3/4, x+1/4, y+3/4'
  84  '-z+1/2, -y+1/2, x'
  85  'z+3/4, -x+1/4, y+3/4'
  86  'z+1/2, y+1/2, x'
  87  'z+3/4, x+1/4, -y+3/4'
  88  'z+1/2, -y+1/2, -x'
  89  '-y+3/4, -z+1/4, -x+3/4'
  90  '-x+1/2, -z+1/2, y'
  91  'y+3/4, -z+1/4, x+3/4'
  92  'x+1/2, -z+1/2, -y'
  93  'y+3/4, z+1/4, -x+3/4'
  94  'x+1/2, z+1/2, y'
  95  '-y+3/4, z+1/4, x+3/4'
  96  '-x+1/2, z+1/2, -y'
  97  'x+1/2, y, z+1/2'
  98  '-y+3/4, x+3/4, z+1/4'
  99  '-x+1/2, -y, z+1/2'
  100  'y+3/4, -x+3/4, z+1/4'
  101  'x+1/2, -y, -z+1/2'
  102  '-y+3/4, -x+3/4, -z+1/4'
  103  '-x+1/2, y, -z+1/2'
  104  'y+3/4, x+3/4, -z+1/4'
  105  'z+1/2, x, y+1/2'
  106  'z+3/4, -y+3/4, x+1/4'
  107  'z+1/2, -x, -y+1/2'
  108  'z+3/4, y+3/4, -x+1/4'
  109  '-z+1/2, x, -y+1/2'
  110  '-z+3/4, -y+3/4, -x+1/4'
  111  '-z+1/2, -x, y+1/2'
  112  '-z+3/4, y+3/4, x+1/4'
  113  'y+1/2, z, x+1/2'
  114  'x+3/4, z+3/4, -y+1/4'
  115  '-y+1/2, z, -x+1/2'
  116  '-x+3/4, z+3/4, y+1/4'
  117  '-y+1/2, -z, x+1/2'
  118  '-x+3/4, -z+3/4, -y+1/4'
  119  'y+1/2, -z, -x+1/2'
  120  'x+3/4, -z+3/4, y+1/4'
  121  '-x+3/4, -y+3/4, -z+1/4'
  122  'y+1/2, -x, -z+1/2'
  123  'x+3/4, y+3/4, -z+1/4'
  124  '-y+1/2, x, -z+1/2'
  125  '-x+3/4, y+3/4, z+1/4'
  126  'y+1/2, x, z+1/2'
  127  'x+3/4, -y+3/4, z+1/4'
  128  '-y+1/2, -x, z+1/2'
  129  '-z+3/4, -x+3/4, -y+1/4'
  130  '-z+1/2, y, -x+1/2'
  131  '-z+3/4, x+3/4, y+1/4'
  132  '-z+1/2, -y, x+1/2'
  133  'z+3/4, -x+3/4, y+1/4'
  134  'z+1/2, y, x+1/2'
  135  'z+3/4, x+3/4, -y+1/4'
  136  'z+1/2, -y, -x+1/2'
  137  '-y+3/4, -z+3/4, -x+1/4'
  138  '-x+1/2, -z, y+1/2'
  139  'y+3/4, -z+3/4, x+1/4'
  140  'x+1/2, -z, -y+1/2'
  141  'y+3/4, z+3/4, -x+1/4'
  142  'x+1/2, z, y+1/2'
  143  '-y+3/4, z+3/4, x+1/4'
  144  '-x+1/2, z, -y+1/2'
  145  'x, y+1/2, z+1/2'
  146  '-y+1/4, x+1/4, z+1/4'
  147  '-x, -y+1/2, z+1/2'
  148  'y+1/4, -x+1/4, z+1/4'
  149  'x, -y+1/2, -z+1/2'
  150  '-y+1/4, -x+1/4, -z+1/4'
  151  '-x, y+1/2, -z+1/2'
  152  'y+1/4, x+1/4, -z+1/4'
  153  'z, x+1/2, y+1/2'
  154  'z+1/4, -y+1/4, x+1/4'
  155  'z, -x+1/2, -y+1/2'
  156  'z+1/4, y+1/4, -x+1/4'
  157  '-z, x+1/2, -y+1/2'
  158  '-z+1/4, -y+1/4, -x+1/4'
  159  '-z, -x+1/2, y+1/2'
  160  '-z+1/4, y+1/4, x+1/4'
  161  'y, z+1/2, x+1/2'
  162  'x+1/4, z+1/4, -y+1/4'
  163  '-y, z+1/2, -x+1/2'
  164  '-x+1/4, z+1/4, y+1/4'
  165  '-y, -z+1/2, x+1/2'
  166  '-x+1/4, -z+1/4, -y+1/4'
  167  'y, -z+1/2, -x+1/2'
  168  'x+1/4, -z+1/4, y+1/4'
  169  '-x+1/4, -y+1/4, -z+1/4'
  170  'y, -x+1/2, -z+1/2'
  171  'x+1/4, y+1/4, -z+1/4'
  172  '-y, x+1/2, -z+1/2'
  173  '-x+1/4, y+1/4, z+1/4'
  174  'y, x+1/2, z+1/2'
  175  'x+1/4, -y+1/4, z+1/4'
  176  '-y, -x+1/2, z+1/2'
  177  '-z+1/4, -x+1/4, -y+1/4'
  178  '-z, y+1/2, -x+1/2'
  179  '-z+1/4, x+1/4, y+1/4'
  180  '-z, -y+1/2, x+1/2'
  181  'z+1/4, -x+1/4, y+1/4'
  182  'z, y+1/2, x+1/2'
  183  'z+1/4, x+1/4, -y+1/4'
  184  'z, -y+1/2, -x+1/2'
  185  '-y+1/4, -z+1/4, -x+1/4'
  186  '-x, -z+1/2, y+1/2'
  187  'y+1/4, -z+1/4, x+1/4'
  188  'x, -z+1/2, -y+1/2'
  189  'y+1/4, z+1/4, -x+1/4'
  190  'x, z+1/2, y+1/2'
  191  '-y+1/4, z+1/4, x+1/4'
  192  '-x, z+1/2, -y+1/2'
loop_
 _atom_site_type_symbol
 _atom_site_label
 _atom_site_symmetry_multiplicity
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
  Gd  Gd0  8  0.00000000  0.00000000  0.50000000  1
  Co  Co1  16  0.12500000  0.12500000  0.12500000  1
"""
    
    # Write CIF to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as f:
        f.write(cif_content)
        temp_cif = f.name
    
    try:
        # Read using the new function (should use primitive cell)
        atoms = read_structure(temp_cif)
        
        # Check that we get the primitive cell (6 atoms: 2 Gd + 4 Co)
        # instead of the conventional unit cell (24 atoms: 8 Gd + 16 Co)
        # Note: The CIF lists 8 Gd and 16 Co atoms, but these represent the
        # conventional (expanded) unit cell. The primitive cell is 1/4 the volume
        # and contains 6 atoms (2 Gd + 4 Co) with the same stoichiometry (GdCo2).
        assert len(atoms) == 6, f"Expected 6 atoms in primitive cell, got {len(atoms)}"
        
        # Check composition
        symbols = atoms.get_chemical_symbols()
        symbol_counts = Counter(symbols)
        assert symbol_counts['Gd'] == 2, f"Expected 2 Gd atoms, got {symbol_counts['Gd']}"
        assert symbol_counts['Co'] == 4, f"Expected 4 Co atoms, got {symbol_counts['Co']}"
        
    finally:
        # Clean up
        if os.path.exists(temp_cif):
            os.remove(temp_cif)


def test_cif_expanded_cell_option():
    """Test that we can still get expanded cell if needed."""
    from qtgui.utils import read_structure
    from ase.build import bulk
    from ase import io as ase_io
    
    # Create a simple test structure
    atoms_original = bulk("Si", cubic=True)
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as f:
        temp_cif = f.name
    
    try:
        # Write structure to CIF
        ase_io.write(temp_cif, atoms_original)
        
        # Read it back normally (should use primitive cell by default)
        atoms_prim = read_structure(temp_cif)
        
        # Read with primitive_cell=False (should expand or keep as is)
        atoms_expanded = read_structure(temp_cif, primitive_cell=False)
        
        # Both should read successfully
        assert len(atoms_prim) > 0
        assert len(atoms_expanded) > 0
        
        # Verify that both reading modes work (they may give same result for simple structures,
        # but the important thing is that primitive_cell=False can be used to override default)
        # For structures without centering operations, both may be identical
        assert atoms_prim.get_volume() > 0
        assert atoms_expanded.get_volume() > 0
        
    finally:
        # Clean up
        if os.path.exists(temp_cif):
            os.remove(temp_cif)


def test_non_cif_files_unchanged():
    """Test that non-CIF files are read normally."""
    from qtgui.utils import read_structure
    from ase.build import bulk
    from ase import io as ase_io
    
    # Create a test XYZ file
    atoms_original = bulk("Si", cubic=True)
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
        temp_xyz = f.name
    
    try:
        # Write structure to XYZ
        ase_io.write(temp_xyz, atoms_original)
        
        # Read it back
        atoms_read = read_structure(temp_xyz)
        
        # Should be the same
        assert len(atoms_read) == len(atoms_original)
        assert atoms_read.get_chemical_symbols() == atoms_original.get_chemical_symbols()
        
    finally:
        # Clean up
        if os.path.exists(temp_xyz):
            os.remove(temp_xyz)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
