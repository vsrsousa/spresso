"""
Example demonstrating how to use Hubbard parameters in the qtGui.

This example shows how the qtGui now properly supports both old and new
Hubbard parameter formats according to xespresso's design.
"""

from ase.build import bulk
import numpy as np

# Example 1: Configuration for OLD format (QE < 7.0)
print("=" * 80)
print("Example 1: Old Hubbard Format (QE < 7.0)")
print("=" * 80)

config_old = {
    'calc_type': 'scf',
    'label': 'fe/scf',
    'ecutwfc': 50.0,
    'ecutrho': 400.0,
    'occupations': 'smearing',
    'smearing': 'gaussian',
    'degauss': 0.02,
    'conv_thr': 1.0e-8,
    
    # Hubbard configuration - OLD FORMAT
    'enable_hubbard': True,
    'hubbard_format': 'old',  # Explicitly select old format
    'hubbard_u': {
        'Fe': 4.3,
        'O': 3.0,
    },
    
    'pseudopotentials': {
        'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
        'O': 'O.pbe-n-rrkjus_psl.1.0.0.UPF',
    },
    'qe_version': '6.8',
}

print("\nConfiguration for OLD format:")
print(f"  hubbard_format: {config_old['hubbard_format']}")
print(f"  QE version: {config_old['qe_version']}")
print(f"  Hubbard U values: {config_old['hubbard_u']}")

print("\nResulting input_data structure:")
print("""
  input_data = {
    'SYSTEM': {
      'lda_plus_u': True,
      ...
    },
    'input_ntyp': {
      'Hubbard_U': {
        'Fe': 4.3,
        'O': 3.0
      }
    }
  }
""")

# Example 2: Configuration for NEW format (QE >= 7.0) with explicit orbitals
print("\n" + "=" * 80)
print("Example 2: New Hubbard Format (QE >= 7.0) with explicit orbitals")
print("=" * 80)

config_new = {
    'calc_type': 'scf',
    'label': 'fe/scf',
    'ecutwfc': 50.0,
    'ecutrho': 400.0,
    'occupations': 'smearing',
    'smearing': 'gaussian',
    'degauss': 0.02,
    'conv_thr': 1.0e-8,
    
    # Hubbard configuration - NEW FORMAT
    'enable_hubbard': True,
    'hubbard_format': 'new',  # Explicitly select new format
    'hubbard_u': {
        'Fe': 4.3,
        'O': 3.0,
    },
    'hubbard_orbitals': {
        'Fe': '3d',  # Specify orbital for Fe
        'O': '2p',   # Specify orbital for O
    },
    'hubbard_projector': 'atomic',  # Optional: projector type
    
    'pseudopotentials': {
        'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
        'O': 'O.pbe-n-rrkjus_psl.1.0.0.UPF',
    },
    'qe_version': '7.2',
}

print("\nConfiguration for NEW format:")
print(f"  hubbard_format: {config_new['hubbard_format']}")
print(f"  QE version: {config_new['qe_version']}")
print(f"  Hubbard U values: {config_new['hubbard_u']}")
print(f"  Hubbard orbitals: {config_new['hubbard_orbitals']}")
print(f"  Hubbard projector: {config_new['hubbard_projector']}")

print("\nResulting input_data structure:")
print("""
  input_data = {
    'SYSTEM': {
      'lda_plus_u': True,
      ...
    },
    'hubbard': {
      'projector': 'atomic',
      'u': {
        'Fe-3d': 4.3,
        'O-2p': 3.0
      },
      'v': []
    },
    'qe_version': '7.2'
  }
""")

# Example 3: Configuration with default orbitals (NEW format)
print("\n" + "=" * 80)
print("Example 3: New Format with default orbitals")
print("=" * 80)

config_new_default = {
    'calc_type': 'scf',
    'label': 'femno/scf',
    'ecutwfc': 50.0,
    'ecutrho': 400.0,
    'occupations': 'smearing',
    'conv_thr': 1.0e-8,
    
    # Hubbard configuration - NEW FORMAT without orbital specification
    'enable_hubbard': True,
    'hubbard_format': 'new',
    'hubbard_u': {
        'Fe': 4.3,
        'Mn': 5.0,
        'O': 3.0,
    },
    # Note: No hubbard_orbitals specified - will use defaults
    
    'pseudopotentials': {
        'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
        'Mn': 'Mn.pbe-spn-rrkjus_psl.1.0.0.UPF',
        'O': 'O.pbe-n-rrkjus_psl.1.0.0.UPF',
    },
    'qe_version': '7.2',
}

print("\nConfiguration for NEW format with defaults:")
print(f"  hubbard_format: {config_new_default['hubbard_format']}")
print(f"  QE version: {config_new_default['qe_version']}")
print(f"  Hubbard U values: {config_new_default['hubbard_u']}")
print(f"  Hubbard orbitals: Not specified (will use defaults)")

print("\nDefault orbitals used:")
print("  Fe → 3d (transition metal)")
print("  Mn → 3d (transition metal)")
print("  O  → 2p (p-block element)")

print("\nResulting input_data structure:")
print("""
  input_data = {
    'SYSTEM': {
      'lda_plus_u': True,
      ...
    },
    'hubbard': {
      'projector': 'atomic',
      'u': {
        'Fe-3d': 4.3,
        'Mn-3d': 5.0,
        'O-2p': 3.0
      },
      'v': []
    },
    'qe_version': '7.2'
  }
""")

# Example 4: How to use in the GUI
print("\n" + "=" * 80)
print("Example 4: How to use in the qtGui")
print("=" * 80)

print("""
When using the qtGui:

1. Load a structure (e.g., Fe2O3, MnO, etc.)

2. In the "Calculation Setup" page:
   - Enable "Hubbard (DFT+U) Configuration" checkbox
   - Select format:
     * "Old (QE < 7.0)" for QE versions 6.x and below
     * "New (QE >= 7.0)" for QE versions 7.0 and above
   
   - For each element in your structure:
     * Set the U value (in eV)
     * Select the orbital (only for new format)
       - 3d for transition metals (Fe, Mn, Co, Ni, etc.)
       - 4f for rare earths (Ce, Nd, Gd, etc.)
       - 2p for oxygen, nitrogen
       - 3p for sulfur, phosphorus

3. The qtGui will automatically:
   - Generate the correct format based on your selection
   - Use default orbitals if you don't specify them (new format)
   - Pass the parameters correctly to xespresso

4. The generated input file will have:
   - OLD format: Hubbard_U(i) parameters in SYSTEM namelist
   - NEW format: HUBBARD card with explicit orbital specifications
""")

print("\n" + "=" * 80)
print("Benefits of the fix:")
print("=" * 80)
print("""
✅ Respects xespresso's design and documentation
✅ Supports both old (QE < 7.0) and new (QE >= 7.0) formats
✅ Allows explicit orbital specifications for better precision
✅ Provides sensible defaults for common elements
✅ Auto-detects format based on QE version
✅ Fully compatible with xespresso's Espresso calculator
✅ Matches behavior of command-line xespresso usage
""")

print("\n✅ Examples completed!")
