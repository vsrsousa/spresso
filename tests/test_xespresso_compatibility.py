"""
Test to verify xespresso compatibility fixes.

This test validates the fixes for:
1. pseudo_dir set to "./pseudo" (handled by remote mixin)
2. lda_plus_u only set for old Hubbard format
3. kspacing converted to kpts before passing to Espresso calculator
"""

import pytest
from ase import Atoms


def test_pseudo_dir_set_to_relative_path():
    """Test that pseudo_dir is set to './pseudo' for remote mixin handling."""
    try:
        from gui.calculations.preparation import CalculationPreparation
    except ImportError:
        pytest.skip("GUI module not available")

    # Create a simple test structure
    atoms = Atoms("Fe", positions=[(0, 0, 0)], cell=(3, 3, 3), pbc=True)

    config = {
        "pseudopotentials": {"Fe": "Fe.paw.z_8.ld1.psl.v1.0.0-high.upf"},
        "calc_type": "scf",
        "ecutwfc": 50.0,
        "ecutrho": 400.0,
        "occupations": "smearing",
        "smearing": "gaussian",
        "degauss": 0.02,
        "conv_thr": 1.0e-8,
    }

    prep = CalculationPreparation(atoms, config, label="test")
    prepared_atoms, calc = prep.prepare()

    # Check that pseudo_dir is set to './pseudo'
    input_data = calc.parameters.get("input_data", {})
    assert "pseudo_dir" in input_data, "pseudo_dir should be set in input_data"
    assert input_data["pseudo_dir"] == "./pseudo", "pseudo_dir should be './pseudo'"


def test_lda_plus_u_not_set_for_new_hubbard_format():
    """Test that lda_plus_u is NOT set when using new HUBBARD card format."""
    try:
        from gui.calculations.preparation import CalculationPreparation
    except ImportError:
        pytest.skip("GUI module not available")

    # Create a simple test structure
    atoms = Atoms("Fe", positions=[(0, 0, 0)], cell=(3, 3, 3), pbc=True)

    config = {
        "pseudopotentials": {"Fe": "Fe.paw.z_8.ld1.psl.v1.0.0-high.upf"},
        "calc_type": "scf",
        "ecutwfc": 50.0,
        "ecutrho": 400.0,
        "enable_hubbard": True,
        "hubbard_format": "new",
        "hubbard_u": {"Fe": 4.0},
        "hubbard_orbitals": {"Fe": "3d"},
        "qe_version": "7.4.1",
    }

    prep = CalculationPreparation(atoms, config, label="test")
    prepared_atoms, calc = prep.prepare()

    # Check that lda_plus_u is not in input_data
    input_data = calc.parameters.get("input_data", {})
    assert "lda_plus_u" not in input_data, "lda_plus_u should not be set for new format"


def test_lda_plus_u_is_set_for_old_hubbard_format():
    """Test that lda_plus_u IS set when using old format."""
    try:
        from gui.calculations.preparation import CalculationPreparation
    except ImportError:
        pytest.skip("GUI module not available")

    # Create a simple test structure
    atoms = Atoms("Fe", positions=[(0, 0, 0)], cell=(3, 3, 3), pbc=True)

    config = {
        "pseudopotentials": {"Fe": "Fe.paw.z_8.ld1.psl.v1.0.0-high.upf"},
        "calc_type": "scf",
        "ecutwfc": 50.0,
        "ecutrho": 400.0,
        "enable_hubbard": True,
        "hubbard_format": "old",
        "hubbard_u": {"Fe": 4.0},
        "qe_version": "6.8",
    }

    prep = CalculationPreparation(atoms, config, label="test")
    prepared_atoms, calc = prep.prepare()

    # Check that lda_plus_u IS in input_data
    input_data = calc.parameters.get("input_data", {})
    assert "lda_plus_u" in input_data, "lda_plus_u should be set for old format"
    assert input_data["lda_plus_u"] == True


def test_kspacing_converted_to_kpts():
    """Test that kspacing is converted to kpts before passing to calculator."""
    try:
        from gui.calculations.preparation import CalculationPreparation
    except ImportError:
        pytest.skip("GUI module not available")

    # Create a simple test structure
    atoms = Atoms("Fe", positions=[(0, 0, 0)], cell=(3, 3, 3), pbc=True)

    config = {
        "pseudopotentials": {"Fe": "Fe.paw.z_8.ld1.psl.v1.0.0-high.upf"},
        "calc_type": "scf",
        "ecutwfc": 50.0,
        "ecutrho": 400.0,
        "kspacing": 0.2,  # Pass kspacing
    }

    prep = CalculationPreparation(atoms, config, label="test")
    prepared_atoms, calc = prep.prepare()

    # Check that kpts is in calculator parameters (not kspacing)
    assert "kpts" in calc.parameters, "kpts should be in calculator parameters"
    assert (
        "kspacing" not in calc.parameters
    ), "kspacing should not be in calculator parameters"
    # kpts should be a tuple of 3 integers
    kpts = calc.parameters["kpts"]
    assert (
        isinstance(kpts, tuple) and len(kpts) == 3
    ), "kpts should be a tuple of 3 values"


def test_kpts_used_when_provided_directly():
    """Test that kpts is used when provided directly (without kspacing)."""
    try:
        from gui.calculations.preparation import CalculationPreparation
    except ImportError:
        pytest.skip("GUI module not available")

    # Create a simple test structure
    atoms = Atoms("Fe", positions=[(0, 0, 0)], cell=(3, 3, 3), pbc=True)

    config = {
        "pseudopotentials": {"Fe": "Fe.paw.z_8.ld1.psl.v1.0.0-high.upf"},
        "calc_type": "scf",
        "ecutwfc": 50.0,
        "ecutrho": 400.0,
        "kpts": (2, 2, 2),  # Pass kpts directly
    }

    prep = CalculationPreparation(atoms, config, label="test")
    prepared_atoms, calc = prep.prepare()

    # Check that kpts is in calculator parameters
    assert "kpts" in calc.parameters, "kpts should be passed to calculator"
    assert calc.parameters["kpts"] == (2, 2, 2)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
