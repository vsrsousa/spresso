#!/usr/bin/env python3
"""
Demo script showing the new pseudopotential auto-extraction feature.

This demonstrates the enhancement to the problem statement where users
can now pass the entire config dict from setup_magnetic_config directly
to Espresso without manually extracting pseudopotentials.
"""

import tempfile
import os
from ase.build import bulk
from ase import Atoms
from xespresso.tools import setup_magnetic_config
from xespresso.xio import write_espresso_in


def demo_simple_afm():
    """Demo 1: Simple AFM Fe system."""
    print("=" * 70)
    print("Demo 1: Simple AFM Fe system")
    print("=" * 70)
    
    atoms = bulk('Fe', cubic=True)
    
    # Setup magnetic configuration with base pseudopotential
    config = setup_magnetic_config(
        atoms, 
        {'Fe': [1, -1]},
        pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'}
    )
    
    print("\n✓ Config generated with setup_magnetic_config")
    print(f"  - Species: {list(config['species_map'].keys())}")
    print(f"  - Pseudopotentials: {list(config['pseudopotentials'].keys())}")
    print(f"  - Magnetization: {config['input_ntyp']['starting_magnetization']}")
    
    # Write input file - EASY WAY!
    with tempfile.TemporaryDirectory() as tmpdir:
        input_file = os.path.join(tmpdir, 'demo1.pwi')
        
        write_espresso_in(
            input_file,
            config['atoms'],
            input_data=config,  # Pass entire dict - that's it!
            ecutwfc=30.0,
            nspin=2,
            kpts=(4, 4, 4)
        )
        
        print("\n✓ Input file written successfully!")
        print(f"  Location: {input_file}")
        
        # Show excerpt
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
        print("\n--- ATOMIC_SPECIES section ---")
        in_species = False
        for line in lines:
            if 'ATOMIC_SPECIES' in line:
                in_species = True
            if in_species:
                print(f"  {line.rstrip()}")
                if line.strip() == '':
                    break
    
    print("\n✓ No manual pseudopotential extraction needed!")


def demo_with_hubbard():
    """Demo 2: AFM with Hubbard U parameters."""
    print("\n" + "=" * 70)
    print("Demo 2: AFM Fe with Hubbard U")
    print("=" * 70)
    
    atoms = bulk('Fe', cubic=True)
    
    # Combine magnetic configuration with Hubbard U parameters
    # This demonstrates how Hubbard U can be specified alongside magnetization
    config = setup_magnetic_config(
        atoms,
        {'Fe': {'mag': [1, -1], 'U': 4.3}},
        pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'}
    )
    
    print("\n✓ Config generated with Hubbard parameters")
    print(f"  - Magnetization: {config['input_ntyp']['starting_magnetization']}")
    print(f"  - Hubbard U: {config['input_ntyp']['Hubbard_U']}")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        input_file = os.path.join(tmpdir, 'demo2.pwi')
        
        write_espresso_in(
            input_file,
            config['atoms'],
            input_data=config,  # Includes Hubbard params
            ecutwfc=30.0,
            nspin=2,
            lda_plus_u=True,
            kpts=(4, 4, 4)
        )
        
        print("\n✓ Input file with Hubbard parameters written!")


def demo_multi_element():
    """Demo 3: Multi-element system."""
    print("\n" + "=" * 70)
    print("Demo 3: Multi-element FeMn system")
    print("=" * 70)
    
    atoms = Atoms('Fe2Mn2', positions=[
        [0, 0, 0], [1.5, 0, 0],
        [0, 1.5, 0], [1.5, 1.5, 0]
    ])
    atoms.cell = [5, 5, 5]
    
    # Setup with multiple elements
    config = setup_magnetic_config(
        atoms,
        {
            'Fe': [1],      # Both Fe equivalent
            'Mn': [1, -1]   # Mn AFM
        },
        pseudopotentials={
            'Fe': 'Fe.pbe-spn.UPF',
            'Mn': 'Mn.pbe-spn.UPF'
        }
    )
    
    print("\n✓ Config generated for multi-element system")
    print(f"  - Species: {list(config['species_map'].keys())}")
    print(f"  - Pseudopotentials: {list(config['pseudopotentials'].keys())}")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        input_file = os.path.join(tmpdir, 'demo3.pwi')
        
        write_espresso_in(
            input_file,
            config['atoms'],
            input_data=config,
            ecutwfc=30.0,
            nspin=2,
            kpts=(4, 4, 4)
        )
        
        print("\n✓ Input file for FeMn system written!")
        
        # Show species section
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
        print("\n--- ATOMIC_SPECIES section ---")
        in_species = False
        for line in lines:
            if 'ATOMIC_SPECIES' in line:
                in_species = True
            if in_species:
                print(f"  {line.rstrip()}")
                if line.strip() == '':
                    break


def demo_comparison():
    """Demo 4: Comparison of old vs new method."""
    print("\n" + "=" * 70)
    print("Demo 4: Old vs New Method Comparison")
    print("=" * 70)
    
    atoms = bulk('Fe', cubic=True)
    config = setup_magnetic_config(
        atoms,
        {'Fe': [1, -1]},
        pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'}
    )
    
    print("\n--- OLD METHOD (still works) ---")
    print("""
# Manual extraction required
calc = Espresso(
    atoms=config['atoms'],
    pseudopotentials=config['pseudopotentials'],  # <-- Extract manually
    input_data={'input_ntyp': config['input_ntyp']},
    nspin=2,
    ecutwfc=40
)
    """)
    
    print("\n--- NEW METHOD (recommended) ---")
    print("""
# No extraction needed!
calc = Espresso(
    atoms=config['atoms'],
    input_data=config,  # <-- Pass entire dict
    nspin=2,
    ecutwfc=40
)
    """)
    
    print("✓ Much cleaner and less error-prone!")


if __name__ == '__main__':
    print("\n")
    print("*" * 70)
    print(" Pseudopotential Auto-Extraction Demo")
    print(" Solution to the problem statement")
    print("*" * 70)
    
    demo_simple_afm()
    demo_with_hubbard()
    demo_multi_element()
    demo_comparison()
    
    print("\n" + "=" * 70)
    print("✓ All demos completed successfully!")
    print("=" * 70)
    print("\nFor more details, see:")
    print("  - PSEUDOPOTENTIAL_AUTO_EXTRACTION.md")
    print("  - MAGNETIC_HELPERS.md")
    print("\n")
