#!/usr/bin/env python
"""
Example: DFT+U+J magnetic configuration with J and J0 parameters.

This example demonstrates:
1. Adding J (exchange) parameters to Hubbard configuration
2. Using J0 as an alternative to J
3. Both old and new QE format support
4. Flexible orbital specification
"""

from ase.build import bulk
from xespresso import Espresso
from xespresso.tools import setup_magnetic_config


def example1_simple_j_parameter():
    """Example 1: Simple AFM Fe with J parameter (new format)."""
    print("=" * 60)
    print("Example 1: AFM Fe with J parameter (new format)")
    print("=" * 60)
    
    # Create BCC Fe structure (2x1x1 to get 2 Fe atoms)
    atoms = bulk('Fe', 'bcc', a=2.87, cubic=True)
    atoms = atoms * (2, 1, 1)
    atoms.set_pbc([True, True, True])
    
    # Setup magnetic configuration with J parameter
    config = setup_magnetic_config(
        atoms,
        magnetic_config={
            'Fe': {
                'mag': [1, -1],           # Antiferromagnetic
                'U': {'3d': 4.3},         # Hubbard U on 3d orbital
                'J': {'3d': 0.4}          # J parameter on 3d orbital (NEW!)
            }
        },
        qe_version='7.2',
        projector='ortho-atomic'
    )
    
    print("\nGenerated species mapping:")
    for species, element in config['species_map'].items():
        print(f"  {species} → {element}")
    
    if 'hubbard' in config:
        print("\nHubbard configuration (new format):")
        print(f"  Projector: {config['hubbard']['projector']}")
        print(f"  U parameters: {config['hubbard']['u']}")
        print(f"  J parameters: {config['hubbard']['j']}")
    
    return config


def example2_j0_alternative():
    """Example 2: Using J0 as alternative to J (new format)."""
    print("\n" + "=" * 60)
    print("Example 2: Using J0 as alternative to J (new format)")
    print("=" * 60)
    
    atoms = bulk('Fe', 'bcc', a=2.87, cubic=True)
    atoms = atoms * (2, 1, 1)
    atoms.set_pbc([True, True, True])
    
    # Setup with J0 instead of J
    config = setup_magnetic_config(
        atoms,
        magnetic_config={
            'Fe': {
                'mag': [1, -1],
                'U': {'3d': 4.3},
                'J0': {'3d': 0.4}         # J0 instead of J (NEW!)
            }
        },
        qe_version='7.2'
    )
    
    print("\nHubbard configuration with J0:")
    if 'hubbard' in config:
        print(f"  U parameters: {config['hubbard']['u']}")
        print(f"  J0 parameters: {config['hubbard']['j0']}")
    
    return config


def example3_multiple_j_values():
    """Example 3: Different J values for different species."""
    print("\n" + "=" * 60)
    print("Example 3: Different J values per species")
    print("=" * 60)
    
    atoms = bulk('Fe', 'bcc', a=2.87, cubic=True)
    atoms = atoms * (2, 1, 1)
    atoms.set_pbc([True, True, True])
    
    # Different J values for each Fe species (AFM)
    config = setup_magnetic_config(
        atoms,
        magnetic_config={
            'Fe': {
                'mag': [1, -1],
                'U': {'3d': [4.3, 4.5]},   # Different U per species
                'J': {'3d': [0.4, 0.45]}   # Different J per species (NEW!)
            }
        },
        qe_version='7.2'
    )
    
    print("\nHubbard parameters with different values per species:")
    if 'hubbard' in config:
        for key, val in config['hubbard']['u'].items():
            print(f"  {key}: U={val}")
        for key, val in config['hubbard']['j'].items():
            print(f"  {key}: J={val}")
    
    return config


def example4_old_format():
    """Example 4: Old QE format with J parameter."""
    print("\n" + "=" * 60)
    print("Example 4: Old QE format (< 7.0) with J parameter")
    print("=" * 60)
    
    atoms = bulk('Fe', 'bcc', a=2.87, cubic=True)
    atoms = atoms * (2, 1, 1)
    atoms.set_pbc([True, True, True])
    
    # Old format - scalar values without orbitals
    config = setup_magnetic_config(
        atoms,
        magnetic_config={
            'Fe': {
                'mag': [1, -1],
                'U': 4.3,         # Old format: scalar value
                'J': 0.4          # J parameter (NEW!)
            }
        },
        hubbard_format='old'  # Force old format
    )
    
    print("\nOld format Hubbard configuration:")
    print(f"  Hubbard_U: {config['input_ntyp'].get('Hubbard_U', {})}")
    print(f"  Hubbard_J: {config['input_ntyp'].get('Hubbard_J', {})}")
    
    return config


def example5_complex_system():
    """Example 5: Complex system with U, V, J and J0."""
    print("\n" + "=" * 60)
    print("Example 5: Complex system with U, V, J, and J0")
    print("=" * 60)
    
    from ase.build import bulk
    from ase import Atoms
    
    # Create a simple complex structure (rock salt: Fe2O4 would be invalid, use Fe + O)
    atoms = Atoms('Fe2O2', 
                  positions=[[0, 0, 0], [0.5, 0.5, 0.5],
                            [0.5, 0, 0], [0, 0.5, 0]])
    atoms.set_cell([4, 4, 4])
    atoms.set_pbc([True, True, True])
    
    # Complex Hubbard configuration
    config = setup_magnetic_config(
        atoms,
        magnetic_config={
            'Fe': {
                'mag': [1, -1],
                'U': {'3d': 4.3},
                'J': {'3d': 0.4},
                'V': [
                    {
                        'species2': 'O',
                        'orbital1': '3d',
                        'orbital2': '2p',
                        'value': 1.0
                    }
                ]
            },
            'O': {'mag': [0, 0]}
        },
        qe_version='7.2'
    )
    
    print("\nComplex Hubbard configuration:")
    if 'hubbard' in config:
        print(f"  U: {config['hubbard']['u']}")
        print(f"  J: {config['hubbard']['j']}")
        print(f"  V: {len(config['hubbard']['v'])} inter-site interactions")
        for i, v_param in enumerate(config['hubbard']['v']):
            print(f"    V[{i}]: {v_param['species1']}-{v_param['orbital1']} ↔ "
                  f"{v_param['species2']}-{v_param['orbital2']} = {v_param['value']}")
    
    return config


def example6_use_cases():
    """Example 6: Common use cases and patterns."""
    print("\n" + "=" * 60)
    print("Example 6: Common use cases")
    print("=" * 60)
    
    atoms = bulk('Fe', 'bcc', a=2.87, cubic=True)
    atoms = atoms * (2, 1, 1)
    
    # Use Case 1: Typical magnetic Fe with modest J
    print("\nUse Case 1: Fe with U+J (rotationally invariant)")
    config1 = setup_magnetic_config(
        atoms,
        magnetic_config={
            'Fe': {
                'mag': [1, -1],
                'U': {'3d': 4.3},
                'J': {'3d': 0.4}   # Typically J ≈ U/10
            }
        },
        qe_version='7.2'
    )
    print("  ✓ Configured U+J for rotationally invariant DFT+U")
    
    # Use Case 2: Using J0 (simpler alternative)
    print("\nUse Case 2: Fe with U+J0 (simpler formulation)")
    config2 = setup_magnetic_config(
        atoms,
        magnetic_config={
            'Fe': {
                'mag': [1, -1],
                'U': {'3d': 4.3},
                'J0': {'3d': 0.4}  # Alternative to J
            }
        },
        qe_version='7.2'
    )
    print("  ✓ Configured U+J0")
    
    # Use Case 3: Lanthanides with f-orbitals
    print("\nUse Case 3: Lanthanides (Gd) with f-orbitals")
    atoms_gd = bulk('Gd', 'hcp', a=3.6, c=5.8)
    config3 = setup_magnetic_config(
        atoms_gd,
        magnetic_config={
            'Gd': {
                'mag': [5, -5],    # f5 configuration
                'U': {'4f': 5.0},
                'J': {'4f': 0.5}   # Exchange in f-electrons
            }
        },
        qe_version='7.2'
    )
    print("  ✓ Configured Gd with 4f orbitals")
    
    return config1, config2, config3


if __name__ == '__main__':
    # Run all examples
    config1 = example1_simple_j_parameter()
    config2 = example2_j0_alternative()
    config3 = example3_multiple_j_values()
    config4 = example4_old_format()
    config5 = example5_complex_system()
    configs = example6_use_cases()
    
    print("\n" + "=" * 60)
    print("✅ All examples completed successfully!")
    print("=" * 60)
    print("\nKey Features Demonstrated:")
    print("  1. Simple J parameter with orbital specification")
    print("  2. J0 as alternative to J")
    print("  3. Different J values per species (AFM)")
    print("  4. Old QE format compatibility")
    print("  5. Complex systems with U+V+J")
    print("  6. Common use cases (Fe, Gd, etc.)")
