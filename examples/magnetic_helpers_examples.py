"""
Complete examples demonstrating the new magnetic helper functions.

These examples show various use cases for the simplified magnetic configuration API:
1. Simple antiferromagnetic configuration
2. Ferromagnetic configuration
3. Complex magnetic structures
4. Element-specific magnetic moments
"""

from ase.build import bulk
from ase import Atoms
import numpy as np
from xespresso.tools import (
    set_magnetic_moments,
    set_antiferromagnetic,
    set_ferromagnetic
)


def example_1_simple_afm():
    """Example 1: Simple antiferromagnetic Fe."""
    print("\n" + "="*60)
    print("Example 1: Simple Antiferromagnetic Fe")
    print("="*60)
    
    atoms = bulk('Fe', cubic=True)
    
    # Method 1: Using set_magnetic_moments with list
    mag_config = set_magnetic_moments(atoms, [1.0, -1.0])
    
    print("Species created:", list(set(atoms.arrays['species'])))
    print("Magnetizations:", mag_config['input_ntyp']['starting_magnetization'])
    
    # Method 2: Using set_antiferromagnetic (simpler for AFM)
    atoms2 = bulk('Fe', cubic=True)
    afm_config = set_antiferromagnetic(atoms2, [[0], [1]])
    
    print("\nUsing set_antiferromagnetic:")
    print("Species created:", list(set(atoms2.arrays['species'])))
    print("Magnetizations:", afm_config['input_ntyp']['starting_magnetization'])


def example_2_ferromagnetic():
    """Example 2: Ferromagnetic configuration."""
    print("\n" + "="*60)
    print("Example 2: Ferromagnetic Fe")
    print("="*60)
    
    atoms = bulk('Fe', cubic=True)
    
    fm_config = set_ferromagnetic(atoms, magnetic_moment=2.0)
    
    print("Species created:", list(set(atoms.arrays['species'])))
    print("Magnetizations:", fm_config['input_ntyp']['starting_magnetization'])
    print("All atoms have same magnetic moment: 2.0")


def example_3_complex_afm():
    """Example 3: Complex antiferromagnetic structure."""
    print("\n" + "="*60)
    print("Example 3: Complex AFM (2x1x1 supercell, 4 atoms)")
    print("="*60)
    
    atoms = bulk('Fe', cubic=True) * (2, 1, 1)
    
    # Checkerboard AFM: atoms 0,3 spin up, atoms 1,2 spin down
    afm_config = set_antiferromagnetic(
        atoms, 
        [[0, 3], [1, 2]], 
        magnetic_moment=1.5
    )
    
    print("Total atoms:", len(atoms))
    print("Species created:", list(set(atoms.arrays['species'])))
    print("Magnetizations:", afm_config['input_ntyp']['starting_magnetization'])
    
    # Show which atoms have which species
    for i, species in enumerate(atoms.arrays['species']):
        mag = afm_config['input_ntyp']['starting_magnetization'].get(species, 0.0)
        print(f"  Atom {i} ({species}): magnetization = {mag}")


def example_4_mixed_elements():
    """Example 4: Mixed elements with selective magnetization."""
    print("\n" + "="*60)
    print("Example 4: Mixed elements (FeO)")
    print("="*60)
    
    # Create FeO structure
    atoms = Atoms('Fe2O', positions=[[0, 0, 0], [1.5, 0, 0], [0, 1.5, 0]])
    atoms.cell = [5, 5, 5]
    
    # Set magnetic moments: AFM on Fe, none on O
    mag_config = set_magnetic_moments(atoms, {0: 1.0, 1: -1.0, 2: 0.0})
    
    print("Species created:", list(set(atoms.arrays['species'])))
    print("Magnetizations:", mag_config['input_ntyp']['starting_magnetization'])
    print("Note: O has no magnetization in input_ntyp (zero moment)")


def example_5_element_specific_fm():
    """Example 5: Ferromagnetic on specific element only."""
    print("\n" + "="*60)
    print("Example 5: Ferromagnetic only on Fe in FeO")
    print("="*60)
    
    atoms = Atoms('Fe2O2', positions=[
        [0, 0, 0], [1.5, 0, 0],  # Fe atoms
        [0, 1.5, 0], [1.5, 1.5, 0]  # O atoms
    ])
    atoms.cell = [5, 5, 5]
    
    # Set ferromagnetic only on Fe
    fm_config = set_ferromagnetic(atoms, magnetic_moment=2.0, element='Fe')
    
    print("Species created:", list(set(atoms.arrays['species'])))
    print("Magnetizations:", fm_config['input_ntyp']['starting_magnetization'])
    print("Only Fe atoms have magnetic moment")


def example_6_with_pseudopotentials():
    """Example 6: Using with existing pseudopotentials."""
    print("\n" + "="*60)
    print("Example 6: With existing pseudopotentials")
    print("="*60)
    
    atoms = bulk('Fe', cubic=True)
    
    # Start with existing pseudopotentials
    existing_pseudo = {
        'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
    }
    
    mag_config = set_antiferromagnetic(
        atoms, 
        [[0], [1]], 
        pseudopotentials=existing_pseudo
    )
    
    print("Pseudopotentials:", mag_config['pseudopotentials'])
    print("Note: Fe1 automatically gets the same pseudopotential as Fe")
    print("Magnetizations:", mag_config['input_ntyp']['starting_magnetization'])


def example_7_ready_for_espresso():
    """Example 7: Complete setup ready for Espresso calculator."""
    print("\n" + "="*60)
    print("Example 7: Complete Espresso calculator setup")
    print("="*60)
    
    atoms = bulk('Fe', cubic=True)
    
    # Set up AFM configuration
    mag_config = set_antiferromagnetic(atoms, [[0], [1]], magnetic_moment=1.0)
    
    # Add actual pseudopotential files
    mag_config['pseudopotentials']['Fe'] = 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
    mag_config['pseudopotentials']['Fe1'] = 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
    
    print("Configuration ready for Espresso:")
    print("\nPseudopotentials:")
    for species, pseudo in mag_config['pseudopotentials'].items():
        print(f"  {species}: {pseudo}")
    
    print("\nInput data (input_ntyp):")
    print(f"  {mag_config['input_ntyp']}")
    
    print("\nExample Espresso calculator creation:")
    print("""
    from xespresso import Espresso
    
    calc = Espresso(
        pseudopotentials=mag_config['pseudopotentials'],
        label='scf/fe-afm',
        ecutwfc=40,
        occupations='smearing',
        degauss=0.02,
        nspin=2,
        input_data={'input_ntyp': mag_config['input_ntyp']},
        kpts=(4, 4, 4),
    )
    atoms.calc = calc
    energy = atoms.get_potential_energy()
    """)


if __name__ == '__main__':
    print("\n" + "="*60)
    print("MAGNETIC HELPER FUNCTIONS - EXAMPLES")
    print("="*60)
    
    example_1_simple_afm()
    example_2_ferromagnetic()
    example_3_complex_afm()
    example_4_mixed_elements()
    example_5_element_specific_fm()
    example_6_with_pseudopotentials()
    example_7_ready_for_espresso()
    
    print("\n" + "="*60)
    print("All examples completed!")
    print("="*60)
