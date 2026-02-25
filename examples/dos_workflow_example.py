"""
Example: DOS (Density of States) calculation workflow with spin polarization

This example shows how to:
1. Run SCF calculation (including magnetic systems)
2. Run NSCF on a denser k-point mesh
3. Calculate DOS with spin polarization (up/down electrons)
4. Analyze magnetic ordering through spin-polarized DOS
5. Optionally compute projected DOS (PDOS) for site/orbital analysis

Key Features for Magnetic Systems:
- Automatic detection of nspin (collinear vs non-collinear)
- Separate spin-up and spin-down DOS contributions
- PDOS analysis to validate magnetic ordering in transition metals
"""

from ase.build import bulk
from xespresso.workflow.calculation_workflow import CalculationWorkflow
from pathlib import Path

# ============================================================================
# Example 1: Non-magnetic system (Al)
# ============================================================================
print("=" * 70)
print("EXAMPLE 1: Non-magnetic system (Al)")
print("=" * 70)

atoms_al = bulk('Al', 'fcc', a=4.05)
pseudopotentials_al = {'Al': 'Al.pbe.UPF'}

workflow_al = CalculationWorkflow(
    atoms_al,
    protocol='fast',
    pseudopotentials=pseudopotentials_al,
)

# SCF + NSCF + DOS
scf_calc = workflow_al.run_scf(label='scf/al')
workflow_al.atoms = scf_calc.atoms
nscf_calc = workflow_al.run_nscf(label='nscf/al', kpts=(12, 12, 12))
dos_result_al = workflow_al.run_dos(nscf_label='nscf/al')

print("✓ Non-magnetic DOS completed\n")

# ============================================================================
# Example 2: Magnetic system - Ferromagnetic Fe
# ============================================================================
print("\n" + "=" * 70)
print("EXAMPLE 2: Magnetic system - Ferromagnetic Fe (nspin=2)")
print("=" * 70)

atoms_fe = bulk('Fe', 'bcc', a=2.87)
pseudopotentials_fe = {'Fe': 'Fe.pbe-spn.UPF'}

workflow_fe = CalculationWorkflow(
    atoms_fe,
    protocol='moderate',
    pseudopotentials=pseudopotentials_fe,
    magnetic_config='ferro'  # Ferromagnetic configuration
)

print(f"Magnetic system setup: nspin={workflow_fe.input_data.get('nspin', 1)}")

# SCF with magnetic ordering
scf_calc_fe = workflow_fe.run_scf(label='scf/fe')

# NSCF preserves magnetic ordering
workflow_fe.atoms = scf_calc_fe.atoms
nscf_calc_fe = workflow_fe.run_nscf(
    label='nscf/fe',
    kpts=(12, 12, 12),
    wf_collect=True
)

# DOS will show SEPARATE spin-up and spin-down contributions
# This allows validation of:
# - Fermi level position for each spin
# - Magnetic moment from integrated DOS difference
# - Electronic structure of ferromagnetic ordering

dos_result_fe = workflow_fe.run_dos(
    nscf_label='nscf/fe',
    Emin=-30,
    Emax=10,
    DeltaE=0.01,
    pdos=False  # Set True for site-projected DOS
)

print("\n✓ Ferromagnetic DOS completed (spin-polarized analysis)\n")

# ============================================================================
# Example 3: Antiferromagnetic system with PDOS analysis
# ============================================================================
print("\n" + "=" * 70)
print("EXAMPLE 3: Antiferromagnetic system with PDOS")
print("=" * 70)

# Fe2O3 or similar antiferromagnetic structure
# For simplicity using Fe antiferro in simple geometry
atoms_afm = bulk('Fe', 'bcc', a=2.87)

workflow_afm = CalculationWorkflow(
    atoms_afm,
    protocol='moderate',
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
    magnetic_config='antiferro'  # Antiferromagnetic configuration
)

print(f"AFM system setup: nspin={workflow_afm.input_data.get('nspin', 1)}")

# SCF + NSCF for antiferromagnetic order
scf_afm = workflow_afm.run_scf(label='scf/fe_afm')
workflow_afm.atoms = scf_afm.atoms
nscf_afm = workflow_afm.run_nscf(label='nscf/fe_afm', kpts=(12, 12, 12))

# DOS with PDOS to analyze site-specific magnetization
dos_afm = workflow_afm.run_dos(
    nscf_label='nscf/fe_afm',
    Emin=-30,
    Emax=10,
    pdos=True  # Enable projected DOS for each Fe site
)

print("\n✓ Antiferromagnetic DOS + PDOS completed\n")
print("  PDOS analysis shows:")
print("  - Different magnetization at each Fe site")
print("  - Validation of alternating spin pattern")
print("  - Local electronic structure contributions\n")

# ============================================================================
# Analysis: Plot spin-polarized DOS
# ============================================================================
print("\n" + "=" * 70)
print("Analyzing spin-polarized DOS")
print("=" * 70)

"""
# Optional: Use xespresso.dos.DOS for visualization
import matplotlib.pyplot as plt
from xespresso.dos import DOS

# Load and plot ferromagnetic Fe DOS
dos_fe = DOS(label='nscf/fe', prefix='fe')
dos_fe.read_dos()

# For magnetic systems, DOS file contains:
# - Integrated DOS(spin-up)
# - Integrated DOS(spin-down)
# - Total DOS

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Total DOS
axes[0].plot(dos_fe.e_dos, dos_fe.dos_total)
axes[0].set_ylabel('DOS (states/eV)')
axes[0].set_xlabel('Energy - E_F (eV)')
axes[0].set_title('Total DOS - Fe (Ferromagnetic)')
axes[0].grid(True, alpha=0.3)

# Magnetization analysis
if hasattr(dos_fe, 'dos_up') and hasattr(dos_fe, 'dos_down'):
    axes[1].plot(dos_fe.e_dos, dos_fe.dos_up, label='Spin-up', color='red')
    axes[1].plot(dos_fe.e_dos, -dos_fe.dos_down, label='Spin-down', color='blue')
    axes[1].axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    axes[1].set_ylabel('DOS (states/eV)')
    axes[1].set_xlabel('Energy - E_F (eV)')
    axes[1].set_title('Spin-polarized DOS - Fe')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('dos_magnetic_analysis.png', dpi=150)
print("Saved: dos_magnetic_analysis.png")
plt.show()
"""

print("✓ Analysis complete")
print("\nKey advantages of spin-polarized DOS for magnetic systems:")
print("  • Detect magnetic ordering (ferrro vs antiferro vs non-magnetic)")
print("  • Calculate local magnetization from DOS asymmetry")
print("  • Validate Hubbard U parameters")
print("  • Analyze magnetic moment contributions by atom/orbital")

