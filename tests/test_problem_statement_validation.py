"""
Validation test for the problem statement requirements.

This test validates that the fixes address all issues mentioned in the problem statement:
1. pseudo_dir is set to './pseudo' (not a full path)
2. lda_plus_u is NOT set when using new Hubbard format (QE >= 7.0)
3. kspacing is converted to kpts and passed as K_POINTS automatic format
"""

import pytest
import json
import numpy as np


def test_config_from_problem_statement():
    """Test using the exact configuration from the problem statement."""
    try:
        from gui.calculations.preparation import CalculationPreparation
        from ase import Atoms
    except ImportError:
        pytest.skip("Required modules not available")

    # Create CuGd structure (simplified for test)
    atoms = Atoms(
        "CuGd", positions=[(0, 0, 0), (1.5, 1.5, 1.5)], cell=(3, 3, 3), pbc=True
    )

    # Configuration from problem statement
    config = {
        "calc_type": "scf",
        "label": "CuGd/scf",
        "ecutwfc": 50.0,
        "ecutrho": 400.0,
        "occupations": "smearing",
        "conv_thr": 1e-08,
        "smearing": "marzari-vanderbilt",
        "degauss": 0.02,
        "kspacing": 0.2,  # This should be converted to kpts
        "pseudopotentials": {
            "Cu": "Cu.paw.z_11.ld1.psl.v1.0.0-low.upf",
            "Gd": "Gd.paw.z_18.atompaw.wentzcovitch.v1.2.upf",
        },
        "qe_version": "7.4.1",
        "enable_magnetism": True,
        "magnetic_config": {"Gd": [7.0]},
        "enable_hubbard": True,
        "hubbard_format": "new",  # NEW format
        "hubbard_u": {"Gd": 6.0},
        "hubbard_orbitals": {"Gd": "4f"},
    }

    prep = CalculationPreparation(atoms, config, label="CuGd/scf")
    prepared_atoms, calc = prep.prepare()

    # Validate Fix 1: pseudo_dir should be './pseudo' (not a full path)
    input_data = calc.parameters.get("input_data", {})
    assert "pseudo_dir" in input_data, "pseudo_dir must be set"
    assert (
        input_data["pseudo_dir"] == "./pseudo"
    ), f"pseudo_dir must be './pseudo', not '{input_data.get('pseudo_dir')}'"

    # Validate Fix 2: lda_plus_u should NOT be set for new format
    assert (
        "lda_plus_u" not in input_data
    ), "lda_plus_u must NOT be set when using new Hubbard format (QE >= 7.0)"

    # Validate Fix 3: kspacing should be converted to kpts
    assert "kpts" in calc.parameters, "kpts must be set (converted from kspacing)"
    assert (
        "kspacing" not in calc.parameters
    ), "kspacing should not be in calculator parameters"

    kpts = calc.parameters["kpts"]
    assert (
        isinstance(kpts, tuple) and len(kpts) == 3
    ), f"kpts should be a tuple of 3 integers, got {kpts}"

    # All kpts values should be integers and greater than 0
    assert all(
        isinstance(k, (int, np.integer)) and k > 0 for k in kpts
    ), f"All kpts values should be positive integers, got {kpts}"

    print("\n✅ All validations passed!")
    print(f"   - pseudo_dir = '{input_data['pseudo_dir']}'")
    print(f"   - lda_plus_u not set (correct for new format)")
    print(f"   - kspacing 0.2 converted to kpts {kpts}")
    print(f"   - K_POINTS will be written as: K_POINTS automatic")
    print(f"     {kpts[0]} {kpts[1]} {kpts[2]}  0 0 0")


def test_kpoints_format_automatic():
    """Test that K_POINTS will be written in automatic format."""
    try:
        from xespresso.xio import build_kpts_str
        from ase import Atoms
    except ImportError:
        pytest.skip("xespresso module not available")

    # Create a simple structure
    atoms = Atoms("Fe", positions=[(0, 0, 0)], cell=(3, 3, 3), pbc=True)

    # Test with kpts (as would be generated from kspacing)
    kpts = (5, 5, 5)
    koffset = (0, 0, 0)

    kpts_str = build_kpts_str(atoms, kspacing=None, kpts=kpts, koffset=koffset)

    # Join the list to get the full string
    kpts_output = "".join(kpts_str)

    # Validate the format
    assert "K_POINTS automatic" in kpts_output, "Should use K_POINTS automatic"
    assert "5 5 5  0 0 0" in kpts_output, "Should have correct grid and offsets"

    print("\n✅ K_POINTS format validation passed!")
    print("Generated K_POINTS block:")
    print(kpts_output)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
