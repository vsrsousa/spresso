"""
Example: PROJWFC (Projected Density of States) Analysis Workflow

This example demonstrates how to use the run_projwfc() method to calculate
orbital-projected properties from a prior NSCF calculation.

PROJWFC is useful for understanding which atoms and orbitals contribute
to the electronic structure, especially important for:
- Validating DFT+U/DFT+DMFT calculations
- Understanding magnetic ordering in transition metal compounds
- Analyzing band character and localization
- Orbital-resolved analysis of d-band metals
"""

from ase.build import bulk
from ase import Atoms
from ase.io import read
from xespresso.workflow.calculation_workflow import CalculationWorkflow


# ============================================================================
# Example 1: Non-magnetic system (Si) - Simple PDOS analysis
# ============================================================================
print("="*70)
print("Example 1: Simple PDOS Analysis (Non-magnetic Si)")
print("="*70)

atoms_si = bulk('Si', 'diamond', a=5.43)

workflow_si = CalculationWorkflow(
    atoms=atoms_si,
    pseudopotentials={'Si': 'Si.pbe.UPF'},
    protocol='moderate',
    kspacing=0.3,
)

# Run SCF calculation
print("\n1. Running SCF calculation...")
scf_calc = workflow_si.run_scf(label='si_scf')
workflow_si.atoms = scf_calc.atoms  # Update atoms after SCF

# Run NSCF calculation with denser k-grid for PDOS
print("\n2. Running NSCF calculation...")
nscf_calc = workflow_si.run_nscf(label='si_nscf', kpts=(10, 10, 10))

# Run PROJWFC for orbital-projected DOS
print("\n3. Running PROJWFC (orbital projections)...")
projwfc_calc = workflow_si.run_projwfc(
    nscf_label='si_nscf',
    projwfc_label='si_projwfc',
    Emin=-30,     # 30 eV below Fermi
    Emax=10,      # 10 eV above Fermi
    DeltaE=0.01,  # 0.01 eV energy grid
)

print("\n✓ PDOS calculated. Output structure:")
print("  si_nscf/projwfc/")
print("    ├── si_nscf.pdos_*  (orbital projections)")
print("    └── si_nscf.projwfco (output file)")
print("\nNext steps:")
print("  from xespresso.dos import DOS")
print("  dos = DOS(label='si_nscf', prefix='si_nscf')")
print("  dos.read_pdos()")
print("  dos.plot_pdos(Emin=-20, Emax=10)")


# ============================================================================
# Example 2: Magnetic system (Fe) - Spin-polarized PDOS for validation
# ============================================================================
print("\n" + "="*70)
print("Example 2: Magnetic PDOS Analysis (Fe)")
print("="*70)

atoms_fe = bulk('Fe', 'bcc', a=2.87)

workflow_fe = CalculationWorkflow(
    atoms=atoms_fe,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
    protocol='moderate',
    kspacing=0.3,
    magnetic_config='ferro',  # Set up ferromagnetic configuration
)

print("\nMagnetic configuration:")
print(f"  nspin={workflow_fe.input_data.get('nspin', 1)}")
print(f"  Initial magnetic moments: {workflow_fe.input_data.get('starting_magnetization', {})}")

# Run SCF with magnetization
print("\n1. Running magnetic SCF calculation...")
scf_calc_fe = workflow_fe.run_scf(label='fe_scf')
workflow_fe.atoms = scf_calc_fe.atoms

# NSCF for magnetic system
print("\n2. Running NSCF calculation (spin-polarized)...")
nscf_calc_fe = workflow_fe.run_nscf(label='fe_nscf', kpts=(12, 12, 12))

# Spin-polarized PDOS
print("\n3. Running spin-polarized PROJWFC...")
projwfc_calc_fe = workflow_fe.run_projwfc(
    nscf_label='fe_nscf',
    projwfc_label='fe_projwfc',
    Emin=-30,
    Emax=10,
    DeltaE=0.01,
    lsym=1,      # Enable symmetrization
    pawproj=0,   # Use Rydberg projectors
)

print("\n✓ Spin-polarized PDOS calculated. Structure:")
print("  fe_nscf/projwfc/")
print("    ├── fe_nscf.pdos_*  (spin-up and spin-down projections)")
print("    └── fe_nscf.projwfco")
print("\nAnalysis:")
print("  • Compare d-orbital projections for spin-up vs spin-down")
print("  • Validate magnetic moment through integrated PDOS")
print("  • Examine band character and localization")


# ============================================================================
# Example 3: Transition metal oxide with DFT+U - Orbitally-resolved analysis
# ============================================================================
print("\n" + "="*70)
print("Example 3: DFT+U with PDOS (Transition Metal Oxide)")
print("="*70)

# Create simple antiferromagnetic Fe2O3 structure (pseudocubic approximation)
from ase.lattice.cubic import Diamond
a = 3.0
atoms_feo = bulk('Fe') * (2, 2, 2)  # For simplicity, using cubic Fe supercell
atoms_feo.set_pbc(True)

workflow_feo = CalculationWorkflow(
    atoms=atoms_feo,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
    protocol='moderate',
    kspacing=0.25,
    magnetic_config={  # Element-based magnetic configuration
        'Fe': {
            'mag': [3.0, -3.0] * 4,  # Alternating up/down for 8 Fe atoms
            'U': {
                '3d': 4.3  # Hubbard U for 3d orbitals
            }
        }
    }
)

print("\nDFT+U Configuration:")
print(f"  Hubbard parameters: {workflow_feo.input_data.get('hubbard_u', {})}")
print(f"  Magnetization: {workflow_feo.input_data.get('starting_magnetization', {})}")

print("\n1. Running DFT+U SCF...")
scf_calc_feo = workflow_feo.run_scf(label='feo_scf')
workflow_feo.atoms = scf_calc_feo.atoms

print("\n2. Running NSCF for PDOS...")
nscf_calc_feo = workflow_feo.run_nscf(label='feo_nscf', kpts=(8, 8, 8))

print("\n3. Running DFT+U PROJWFC...")
projwfc_calc_feo = workflow_feo.run_projwfc(
    nscf_label='feo_nscf',
    projwfc_label='feo_projwfc',
    Emin=-30,
    Emax=10,
    DeltaE=0.01,
    filpdos='feo_projwfc',  # Custom output prefix
)

print("\n✓ DFT+U PDOS calculated. Structure:")
print("  feo_nscf/projwfc/")
print("    └── feo_projwfc.pdos_*  (d-orbital projections)")
print("\nAnalysis for DFT+U:")
print("  • Examine how Hubbard U affects orbital occupations")
print("  • Compare occupied vs unoccupied d-orbital projections")
print("  • Validate magnetic moment through integrated d-PDOS")
print("  • Look for band splitting from magnetic ordering")


# ============================================================================
# Example 4: Custom PDOS parameters for fine energy resolution
# ============================================================================
print("\n" + "="*70)
print("Example 4: High-Resolution PDOS near Band Gap")
print("="*70)

atoms_hse = bulk('GaAs')

workflow_hse = CalculationWorkflow(
    atoms=atoms_hse,
    pseudopotentials={'Ga': 'Ga.pbe.UPF', 'As': 'As.pbe.UPF'},
    protocol='accurate',
    kspacing=0.2,
)

print("\n1. Running high-accuracy SCF...")
scf_calc_hse = workflow_hse.run_scf(label='gaas_scf')
workflow_hse.atoms = scf_calc_hse.atoms

print("\n2. Running NSCF with dense k-grid...")
nscf_calc_hse = workflow_hse.run_nscf(label='gaas_nscf', kpts=(16, 16, 16))

print("\n3. Running fine-resolution PROJWFC...")
projwfc_calc_hse = workflow_hse.run_projwfc(
    nscf_label='gaas_nscf',
    projwfc_label='gaas_projwfc',
    Emin=-25,
    Emax=15,
    DeltaE=0.005,  # Fine 0.005 eV resolution for band gap analysis
    degauss=0.005, # Small broadening for sharp features
    ngauss=0,      # Methfessel-Paxton (cold) smearing
)

print("\n✓ High-resolution PDOS calculated.")
print("\nThis high-resolution PDOS reveals:")
print("  • Sharp band edge features near the band gap")
print("  • Orbital character at gap edges (s vs p vs d)")
print("  • Small features from band anticrossings")


print("\n" + "="*70)
print("All PROJWFC examples completed!")
print("="*70)
print("\nCommon PDOS analysis tasks:")
print("  1. Load and plot PDOS:")
print("     from xespresso.dos import DOS")
print("     dos = DOS(label='label', prefix='prefix')")
print("     dos.read_pdos()")
print("     dos.plot_pdos(Emin=-20, Emax=10)")
print("")
print("  2. Compare spin-up vs spin-down (magnetic systems):")
print("     dos.plot_pdos(ions=[0], Emin=-10, Emax=5)")
print("")
print("  3. Integrated PDOS for orbital occupations:")
print("     dos.dos_dict  # Dictionary with all orbital projections")
print("")
print("  4. Site-projected magnetization:")
print("     # Integrate PDOS(up) - PDOS(down) for each atom")
