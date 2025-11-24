import pytest


@pytest.fixture
def bulk_h():
    from ase.build import bulk

    h = bulk("Fe", cubic=True)
    h.set_chemical_symbols(["H", "H"])
    return h


@pytest.fixture
def graphene_monolayer():
    """Create a graphene monolayer with vacuum"""
    from ase.build import graphene

    atoms = graphene()
    # Add vacuum in z-direction for 2D material
    atoms.cell[2, 2] = 15.0
    atoms.center(axis=2)
    # Set all directions as periodic
    atoms.pbc = [True, True, True]
    return atoms
