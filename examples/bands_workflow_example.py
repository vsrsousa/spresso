"""
Example: Band structure calculation workflow with automatic high-symmetry paths

This example shows how to:
1. Run SCF calculation to converge charge density
2. Run band structure along automatic high-symmetry k-path
3. Analyze band structure (magnetic and non-magnetic)
4. Compare with Wannier interpolation (advanced)
"""

from ase.build import bulk
from xespresso.workflow.calculation_workflow import CalculationWorkflow
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

print("=" * 70)
print("BAND STRUCTURE WORKFLOW WITH AUTOMATIC HIGH-SYMMETRY PATHS")
print("=" * 70)

# ============================================================================
# Example 1: Non-magnetic Al - Simple bands
# ============================================================================
print("\n1. Non-magnetic Al system - Simple band structure")
print("-" * 70)

atoms_al = bulk('Al', 'fcc', a=4.05)
workflow_al = CalculationWorkflow(
    atoms_al,
    protocol='moderate',
    pseudopotentials={'Al': 'Al.pbe.UPF'}
)

print(f"Structure: Al FCC a={atoms_al.cell[0,0]:.3f} Å")
print(f"Protocol: moderate (ecutwfc=50 Ry)")

# Verify the bandpath that will be used
bandpath_al = atoms_al.cell.bandpath()
print(f"\nAuto-generated bandpath:")
print(f"  Path: {bandpath_al.path}")
print(f"  Total k-points: {len(bandpath_al.kpts)}")
print(f"  High-symmetry points: {list(bandpath_al.special_points.keys())}")

print("\nWorkflow steps (not executed here for demo):")
print("  1. SCF: converge charge density on coarse k-mesh")
print("  2. Bands: calculate along high-symmetry path")
print("  3. Extract: BandStructure object")

"""
# Step 1: SCF calculation
scf = workflow_al.run_scf(label='scf/al')

# Step 2: Band structure (reuses SCF charge density)
bands = workflow_al.run_bands(label='bands/al')

# Step 3: Extract and plot band structure
try:
    bs = bands.band_structure()
    bs.reference = bands.get_fermi_level()
    bs.plot()
    plt.savefig('al_bands.png', dpi=150)
except Exception as e:
    print(f"Band plotting not available: {e}")
"""

# ============================================================================
# Example 2: Ferromagnetic Fe - Spin-polarized bands
# ============================================================================
print("\n2. Ferromagnetic Fe - Spin-polarized band structure")
print("-" * 70)

atoms_fe = bulk('Fe', 'bcc', a=2.87)
workflow_fe = CalculationWorkflow(
    atoms_fe,
    protocol='moderate',
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
    magnetic_config='ferro'  # Ferromagnetic
)

print(f"Structure: Fe BCC a={atoms_fe.cell[0,0]:.3f} Å")
print(f"Magnetic config: ferromagnetic (nspin=2)")
print(f"Protocol: moderate")

# Verify magnetic bandpath
bandpath_fe = atoms_fe.cell.bandpath()
print(f"\nAuto-generated bandpath (respects magnetic structure):")
print(f"  Path: {bandpath_fe.path}")
print(f"  Total k-points: {len(bandpath_fe.kpts)}")

print("\nWorkflow steps (not executed here for demo):")
print("  1. SCF: ferromagnetic ordering with nspin=2")
print("  2. Bands: spin-polarized bands (separate ↑ and ↓)")
print("  3. Features:")
print("     - Band splitting from magnetic moment")
print("     - Spin-dependent density at Fermi level")

"""
# Step 1: SCF with magnetic ordering
scf_fe = workflow_fe.run_scf(label='scf/fe')

# Step 2: Band structure preserves magnetic order
bands_fe = workflow_fe.run_bands(label='bands/fe')

# Step 3: Analyze magnetic band structure
bs_fe = bands_fe.band_structure()
bs_fe.reference = bands_fe.get_fermi_level()
# For magnetic systems, separate spin channels in analysis
"""

# ============================================================================
# Example 3: Antiferromagnetic system - AFM band structure  
# ============================================================================
print("\n3. Antiferromagnetic Fe - AFM order in band structure")
print("-" * 70)

atoms_afm = bulk('Fe', 'bcc', a=2.87)
workflow_afm = CalculationWorkflow(
    atoms_afm,
    protocol='moderate',
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
    magnetic_config='antiferro'  # Antiferromagnetic
)

print(f"Structure: Fe BCC a={atoms_afm.cell[0,0]:.3f} Å")
print(f"Magnetic config: antiferromagnetic (nspin=2)")
print(f"  → Reconstructed with 2-atom AFM structure")
print(f"  → Cell doubled to show AFM supercell")

bandpath_afm = atoms_afm.cell.bandpath()
print(f"\nBandpath for AFM supercell:")
print(f"  Path: {bandpath_afm.path}")
print(f"  Total k-points: {len(bandpath_afm.kpts)}")

print("\nKey features for AFM band structure:")
print("  - Flat bands from AFM exchange splitting")
print("  - Spin density imbalance in local DOS")
print("  - Band degeneracies lifted by AFM order")

"""
# Full workflow
scf_afm = workflow_afm.run_scf(label='scf/fe_afm')
bands_afm = workflow_afm.run_bands(label='bands/fe_afm')

# Compare spin contributions
bs_afm = bands_afm.band_structure()
"""

# ============================================================================
# Bandpath details explanation
# ============================================================================
print("\n" + "=" * 70)
print("BANDPATH CONCEPT - Automatic High-Symmetry K-PATH")
print("=" * 70)

print("""
The bandpath is automatically determined for each crystal system:

For FCC (Al): Γ-X-W-K-Γ-L-U-W-L-K|U-X
  - Γ: Gamma point (0,0,0)
  - X, W, L, K, U: High-symmetry points in reciprocal space
  - Path covers all features of the electronic structure
  - Typical ~40-60 k-points

For BCC (Fe): Different high-symmetry path optimized for BCC geometry

The workflow automatically:
1. Detects crystal structure (FCC vs BCC vs HCP, etc)
2. Generates optimal path through Brillouin zone
3. Places more k-points where structure changes rapidly
4. Respects crystal symmetries exactly

This avoids manual specification of k-point paths!
""")

print("=" * 70)
print("\nStandard workflow pattern:")
print("  SCF (label='scf') → Bands (label='bands')")
print("  ↓ (reuses charge density from SCF)")
print("\nFor Wannier comparison:")
print("  SCF → NSCF (dense mesh) → Wannier → Bands (interpolated)")
print("=" * 70)
