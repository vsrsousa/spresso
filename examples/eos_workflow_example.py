"""
Examples: Equation of State (EOS) Workflow

Demonstrates how to use the EOSWorkflow class to study and optimize
structures using Birch-Murnaghan equation of state.

This example covers:
1. Simple EOS study for a bulk metal
2. Extended study with custom parameters
3. Comparison of different protocols
4. Data analysis and visualization
5. Using quick_eos() convenience function
"""

from ase.build import bulk
from ase.io import read
import numpy as np
from xespresso.workflow import EOSWorkflow, quick_eos


# ═══════════════════════════════════════════════════════════════════════════════
# EXAMPLE 1: Simple EOS Study (Si)
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "="*80)
print("EXAMPLE 1: Simple EOS Study - Silicon (Diamond Structure)")
print("="*80)

# Create structure
atoms_si = bulk('Si', 'diamond', a=5.43)
print(f"\nStructure: {atoms_si.get_chemical_formula()}")
print(f"Original volume: {atoms_si.get_volume():.4f} Ų")

# Define pseudopotentials
pseudopotentials = {'Si': 'Si.pbe.UPF'}

# Create EOS workflow
eos_si = EOSWorkflow(
    atoms=atoms_si,
    pseudopotentials=pseudopotentials,
    protocol='moderate',
    debug=False  # Set to True for detailed logging
)

print("\n1. Running EOS study (±5% volume range, 7 points)...")
print("   (This would run SCF calculations on the compute cluster)")

# In a real run, you would execute:
# results = eos_si.run_eos_study(
#     volume_range=(0.95, 1.05),
#     n_points=7,
#     label='eos/si',
#     parallel=True
# )

# For this demo, we'll use mock data to show the workflow
print("   Generating synthetic E-V data for demonstration...")

# Create synthetic data (simulating 7 SCF calculations)
from xespresso.workflow.eos_workflow import birch_murnaghan_eos

# True parameters (what we expect to recover)
E0_true = -10.524       # eV
V0_true = 20.35         # Ų
B0_true = 99.0          # GPa
BP_true = 4.0

volumes_test = np.linspace(19.3, 21.4, 7)
energies_test = birch_murnaghan_eos(volumes_test, E0_true, V0_true, B0_true, BP_true)

# Add tiny noise to simulate real calculations
np.random.seed(42)
energies_test += np.random.normal(0, 1e-4, len(energies_test))

# Manually populate results (in real workflow, this is done by run_eos_study)
import pandas as pd
eos_si.results_df = pd.DataFrame({
    'factor': volumes_test / V0_true,
    'volume': volumes_test,
    'energy': energies_test,
})

print(f"\n   Collected {len(eos_si.results_df)} E-V data points")

print("\n2. Fitting Birch-Murnaghan EOS...")
eos_si.fit_eos()

print("\n3. Extracting equilibrium properties...")
props = eos_si.get_eos_properties()

print(f"\n   ✓ Equilibrium Volume (V₀):  {props['v0']:.6f} Ų")
print(f"   ✓ Equilibrium Energy (E₀): {props['e0']:.6f} eV")
print(f"   ✓ Bulk Modulus (B₀):       {props['bulk_modulus']:.2f} GPa")
print(f"   ✓ Pressure Derivative:     {props['bulk_modulus_prime']:.3f}")
print(f"   ✓ Fit Quality (R²):        {props['r_squared']:.6f}")

print("\n4. Visualizing results...")
print("   (In interactive mode, would display plots)")
# eos_si.plot_eos_curve(save_path='eos_si.png')
# eos_si.plot_residuals(save_path='residuals_si.png')

print("\n5. Exporting data...")
# eos_si.to_csv('eos_si_results.csv')
# eos_si.to_json('eos_si_params.json')
print("   (Would save: eos_si_results.csv and eos_si_params.json)")

print("\n" + eos_si.summary())


# ═══════════════════════════════════════════════════════════════════════════════
# EXAMPLE 2: Extended Study with Custom Parameters (Fe)
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "="*80)
print("EXAMPLE 2: Extended EOS Study - Iron (BCC)")
print("="*80)

# Create structure
atoms_fe = bulk('Fe', 'bcc', a=2.87)
print(f"\nStructure: {atoms_fe.get_chemical_formula()}")
print(f"Original volume: {atoms_fe.get_volume():.4f} Ų")

# Create EOS workflow with custom parameters
eos_fe = EOSWorkflow(
    atoms=atoms_fe,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},  # Spin-polarized pseudopotential
    protocol='moderate',
    magnetic_config='ferro',                     # Ferromagnetic configuration
    kspacing=0.2,                                # Custom k-spacing
)

print("\n1. Configuration:")
print(f"   Protocol:       {eos_fe.protocol}")
print(f"   Magnetic config: ferro (all spins parallel)")
print(f"   K-spacing:      0.2 Ų⁻¹")

print("\n2. Running extended volume range study...")
print("   Volume range: 0.93 to 1.07 (±7%)")
print("   Points:       11")

# Create synthetic E-V data for Fe
E0_fe = -7.832
V0_fe = 11.64
B0_fe = 172.0
BP_fe = 4.5

volumes_fe = np.linspace(10.82, 12.45, 11)
energies_fe = birch_murnaghan_eos(volumes_fe, E0_fe, V0_fe, B0_fe, BP_fe)
energies_fe += np.random.normal(0, 2e-4, len(energies_fe))

eos_fe.results_df = pd.DataFrame({
    'factor': volumes_fe / V0_fe,
    'volume': volumes_fe,
    'energy': energies_fe,
})

print(f"   ✓ Collected {len(eos_fe.results_df)} points")

print("\n3. Fitting EOS...")
eos_fe.fit_eos()

print("\n4. Results for Fe:")
props_fe = eos_fe.get_eos_properties()
print(f"   V₀ = {props_fe['v0']:.6f} Ų")
print(f"   E₀ = {props_fe['e0']:.6f} eV")
print(f"   B₀ = {props_fe['bulk_modulus']:.2f} GPa  (experimental: ~170 GPa)")
print(f"   R² = {props_fe['r_squared']:.6f}")


# ═══════════════════════════════════════════════════════════════════════════════
# EXAMPLE 3: Comparing Different Protocols
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "="*80)
print("EXAMPLE 3: Protocol Comparison")
print("="*80)

# Create Au structure
atoms_au = bulk('Au', 'fcc', a=4.08)
print(f"\nStructure: {atoms_au.get_chemical_formula()}")

protocols = ['fast', 'moderate', 'accurate']
eos_results = {}

for protocol in protocols:
    print(f"\nRunning EOS with protocol: {protocol}")
    
    eos = EOSWorkflow(
        atoms=atoms_au,
        pseudopotentials={'Au': 'Au.pbe.UPF'},
        protocol=protocol,
    )
    
    # Create synthetic data (more accurate with higher protocol)
    E0_au = -3.456
    V0_au = 16.98
    B0_au = 180.0
    
    # Tighter convergence with higher protocol → smaller variation
    noise_level = {'fast': 1e-3, 'moderate': 5e-4, 'accurate': 1e-4}[protocol]
    
    volumes_au = np.linspace(16.0, 17.96, 7)
    energies_au = birch_murnaghan_eos(volumes_au, E0_au, V0_au, B0_au, 4.0)
    energies_au += np.random.normal(0, noise_level, len(energies_au))
    
    eos.results_df = pd.DataFrame({
        'factor': volumes_au / V0_au,
        'volume': volumes_au,
        'energy': energies_au,
    })
    
    eos.fit_eos()
    props = eos.get_eos_properties()
    
    eos_results[protocol] = {
        'v0': props['v0'],
        'b0': props['bulk_modulus'],
        'r2': props['r_squared'],
    }
    
    print(f"   V₀ = {props['v0']:.4f} Ų, B₀ = {props['bulk_modulus']:.1f} GPa, R² = {props['r_squared']:.6f}")

print("\nComparison:")
print("┌─────────┬──────────┬─────────┬────────┐")
print("│Protocol │   V₀     │   B₀    │  R²    │")
print("├─────────┼──────────┼─────────┼────────┤")
for protocol in protocols:
    res = eos_results[protocol]
    print(f"│{protocol:7s} │ {res['v0']:8.4f} │ {res['b0']:7.1f} │ {res['r2']:.4f} │")
print("└─────────┴──────────┴─────────┴────────┘")


# ═══════════════════════════════════════════════════════════════════════════════
# EXAMPLE 4: Property Predictions
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "="*80)
print("EXAMPLE 4: Predicting Properties at Different Volumes")
print("="*80)

print("\nUsing Fe EOS from Example 2...")

# Predict energies at different volumes
test_volumes = [10.5, 11.0, 11.64, 12.0, 12.5]

print("\nEnergy predictions:")
print("┌────────────┬──────────────┐")
print("│  Volume    │    Energy    │")
print("│   (Ų)      │    (eV)      │")
print("├────────────┼──────────────┤")

for v in test_volumes:
    E_pred = eos_fe.predict_energy(v)
    print(f"│ {v:10.2f} │ {E_pred:12.6f} │")

print("└────────────┴──────────────┘")

# Predict pressures
print("\nPressure predictions:")
print("┌────────────┬──────────────┐")
print("│  Volume    │  Pressure    │")
print("│   (Ų)      │   (GPa)      │")
print("├────────────┼──────────────┤")

for v in test_volumes:
    P_pred = eos_fe.calculate_pressure(v)
    marker = "  ← Equilibrium" if np.isclose(v, props_fe['v0']) else ""
    print(f"│ {v:10.2f} │ {P_pred:12.3f} │{marker}")

print("└────────────┴──────────────┘")


# ═══════════════════════════════════════════════════════════════════════════════
# EXAMPLE 5: Using quick_eos() Convenience Function
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "="*80)
print("EXAMPLE 5: One-Call EOS Study with quick_eos()")
print("="*80)

print("\nThe quick_eos() function performs complete EOS study in one call:")
print("  1. Create workflow")
print("  2. Run EOS study")
print("  3. Fit EOS")
print("  4. Return workflow object and properties")

print("\nUsage example:")
print("""
from ase.build import bulk
from xespresso.workflow import quick_eos

atoms = bulk('Si', 'diamond', a=5.43)

eos, props = quick_eos(
    atoms=atoms,
    pseudopotentials={'Si': 'Si.pbe.UPF'},
    volume_range=(0.95, 1.05),
    n_points=7,
    protocol='moderate',
    label='eos/si'
)

print(f"Bulk Modulus: {props['bulk_modulus']:.2f} GPa")
print(f"Equilibrium Volume: {props['v0']:.4f} Ų")

# Optional: visualize and export
eos.plot_eos_curve('eos.png')
eos.to_csv('eos_results.csv')
""")


# ═══════════════════════════════════════════════════════════════════════════════
# EXAMPLE 6: Complete Workflow Integration
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "="*80)
print("EXAMPLE 6: Complete Workflow with Different Machines")
print("="*80)

print("""
The EOSWorkflow fully integrates with xespresso's machine and queue system.

LOCAL EXECUTION:
    eos = EOSWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        protocol='moderate',
        queue={'nodes': 1, 'ntasks-per-node': 8, 'time': '1:00:00'}
    )
    
REMOTE EXECUTION (non-blocking):
    eos = EOSWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        protocol='moderate',
        machine='cluster_name',
        code_version='7.2'  # Optional: specify QE version
    )
    
    # Run automatically monitors jobs
    results = eos.run_eos_study(
        volume_range=(0.95, 1.05),
        n_points=7,
        parallel=True
    )

REMOTE EXECUTION (blocking):
    eos = EOSWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        protocol='accurate',
        machine='cluster_name',
        queue={'wait_for_completion': True}
    )
""")


# ═══════════════════════════════════════════════════════════════════════════════
# EXAMPLE 7: Error Handling and Debugging
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "="*80)
print("EXAMPLE 7: Error Handling and Troubleshooting")
print("="*80)

print("""
Error Handling:

1. Insufficient data points:
   >>> eos.run_eos_study(volume_range=(0.95, 1.05), n_points=2)
   ValueError: Need at least 3 points for EOS fitting

2. Fitting before running study:
   >>> eos.fit_eos()
   ValueError: No EOS data collected. Call run_eos_study() first.

3. Invalid volume range:
   >>> eos.scale_volume_uniformly(atoms, -1.0)
   ValueError: scale_factor must be > 0

Debugging Tips:

1. Enable debug logging:
   eos = EOSWorkflow(atoms, pseudopotentials, debug=True)

2. Check intermediate data:
   print(eos.results_df)        # E-V data collected
   print(eos.eos_structures)    # Scaled structures
   print(eos.eos_calcs)         # Calculator objects
   print(eos.error_log)         # Any errors during runs

3. Validate fit quality:
   print(eos.eos_params['r_squared'])
   eos.plot_residuals()         # Visual inspection of residuals
""")


print("\n" + "="*80)
print("End of Examples")
print("="*80)
