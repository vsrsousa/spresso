"""
Example: Convergence Parameter Optimization Workflow

This example demonstrates how to systematically optimize DFT calculation
parameters (ecutwfc, kspacing) for energy and force convergence.

This is useful for:
- Establishing convergence criteria for a new structure
- Finding optimal balance between accuracy and computational cost
- Comparing different pseudopotentials
- Validating results from production calculations
"""

from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow


# ============================================================================
# Example 1: Basic Convergence Study (Si)
# ============================================================================
print("="*80)
print("Example 1: Basic Convergence Study - Silicon")
print("="*80)

atoms_si = bulk('Si', 'diamond', a=5.43)

conv_si = ConvergenceWorkflow(
    atoms=atoms_si,
    pseudopotentials={'Si': 'Si.pbe.UPF'},
    protocol='moderate',
    # Test modest ranges for quick demo
    ecut_range=[40, 50, 60],
    kspacing_range=[0.3, 0.2, 0.15],
)

print("\n1. Running convergence study...")
# Uncomment to run actual calculations:
# results = conv_si.run_convergence_study()
# conv_si.to_csv('convergence_si.csv')

# For demonstration, we'll use mock results
print("   (Skipped actual calculations for this demo)")
print("\n2. Get recommendations...")
# Uncomment after running:
# recommendations = conv_si.get_recommendations(
#     energy_tolerance=1e-4,    # 0.1 meV/atom
#     force_tolerance=0.1,      # eV/Å
# )

print("""
Expected output would be:
  ⚡ FAST: ecutwfc=40 Ry, kspacing=0.3 Å⁻¹
  ⚖️ BALANCED: ecutwfc=50 Ry, kspacing=0.2 Å⁻¹  
  🎯 ACCURATE: ecutwfc=60 Ry, kspacing=0.15 Å⁻¹
""")


# ============================================================================
# Example 2: Detailed Study with Multiple Parameter Ranges
# ============================================================================
print("\n" + "="*80)
print("Example 2: Extended Convergence Study")
print("="*80)

atoms_fe = bulk('Fe', 'bcc', a=2.87)

conv_fe = ConvergenceWorkflow(
    atoms=atoms_fe,
    pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
    protocol='moderate',
    # More comprehensive ranges
    ecut_range=[30, 40, 50, 60, 70, 80],
    kspacing_range=[0.5, 0.4, 0.3, 0.2, 0.15, 0.1],
)

print("\nParameter ranges:")
print(f"  ecutwfc: {conv_fe.ecut_range}")
print(f"  kspacing: {conv_fe.kspacing_range}")
print(f"  Total tests: {len(conv_fe.ecut_range) * len(conv_fe.kspacing_range)}")

print("""
To run this extended study:

    results = conv_fe.run_convergence_study(verbose=True)
    
This will:
1. Create 36 calculation directories (6 × 6 parameter combinations)
2. Run SCF for each combination
3. Extract total energy and forces
4. Save results in DataFrame

Then analyze with:
    
    recommendations = conv_fe.get_recommendations(
        energy_tolerance=1e-4,
        force_tolerance=0.1
    )
    
    conv_fe.plot_convergence(save_path='fe_convergence.png')
    conv_fe.to_csv('fe_convergence_results.csv')
""")


# ============================================================================
# Example 3: Loading Previous Results and Analysis
# ============================================================================
print("\n" + "="*80)
print("Example 3: Analyzing Previous Convergence Study")
print("="*80)

print("""
If you previously ran a convergence study and saved results:

    conv = ConvergenceWorkflow.from_cif(
        'structure.cif',
        pseudopotentials={'Fe': 'Fe.pbe.UPF'}
    )
    
    # Load previous results
    conv.from_csv('convergence_results.csv')
    
    # Get recommendations with different tolerances
    recs_strict = conv.get_recommendations(
        energy_tolerance=0.1e-3,   # 0.1 meV/atom
        force_tolerance=0.05
    )
    
    recs_loose = conv.get_recommendations(
        energy_tolerance=1e-3,    # 1 meV/atom
        force_tolerance=0.2
    )
    
    # Plot results
    conv.plot_convergence(save_path='convergence.png')
""")


# ============================================================================
# Example 4: Using Convergence Results for Production
# ============================================================================
print("\n" + "="*80)
print("Example 4: Using Recommendations for Production Calculations")
print("="*80)

print("""
After optimization, use recommendations in production workflow:

    from xespresso.workflow import ConvergenceWorkflow, CalculationWorkflow
    
    # Step 1: Run convergence study
    conv = ConvergenceWorkflow.from_cif(...)
    conv.run_convergence_study()
    recs = conv.get_recommendations()
    
    # Step 2: Use balanced recommendation
    params = recs['balanced']
    
    # Step 3: Run production workflow with optimized parameters
    workflow = CalculationWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        kspacing=params['kspacing'],  # Use optimized value
    )
    workflow.input_data['ecutwfc'] = params['ecutwfc']
    
    # Now run full workflow with optimized parameters
    scf_calc = workflow.run_scf(label='production/scf')
    nscf_calc = workflow.run_nscf(label='production/nscf')
    dos_calc = workflow.run_dos(nscf_label='production/nscf')
    bands_calc = workflow.run_bands(nscf_label='production/nscf')
""")


# ============================================================================
# Example 5: Comparing Pseudopotentials
# ============================================================================
print("\n" + "="*80)
print("Example 5: Comparing Pseudopotentials")
print("="*80)

print("""
Use ConvergenceWorkflow to compare different pseudopotentials:

    from xespresso.workflow import ConvergenceWorkflow
    
    # Test different Fe pseudopotentials
    for pp_file in ['Fe.pbe-spn.UPF', 'Fe.pbe-dn-rrkjus_psl.1.0.0.UPF']:
        conv = ConvergenceWorkflow(
            atoms=bulk('Fe'),
            pseudopotentials={'Fe': pp_file},
            ecut_range=[40, 60, 80],
            kspacing_range=[0.3, 0.2, 0.1],
        )
        
        conv.run_convergence_study()
        recs = conv.get_recommendations()
        
        print(f"\\nPseudopotential: {pp_file}")
        print(f"  Balanced: ecutwfc={recs['balanced']['ecutwfc']}, "
              f"kspacing={recs['balanced']['kspacing']}")
        print(f"  Energy/atom: {recs['balanced']['energy_per_atom']:.6f} eV")
        
        conv.to_csv(f'convergence_{pp_file}.csv')
""")


# ============================================================================
# Example 6: Magnetic System Convergence
# ============================================================================
print("\n" + "="*80)
print("Example 6: Convergence of Magnetic Systems")
print("="*80)

atoms_mno = bulk('Mn')  # Ferromagnetic Mn for example

conv_mag = ConvergenceWorkflow(
    atoms=atoms_mno,
    pseudopotentials={'Mn': 'Mn.pbe-spn.UPF'},
    protocol='moderate',
    magnetic_config='ferro',  # Magnetic configuration
    ecut_range=[50, 60, 70],
    kspacing_range=[0.3, 0.2, 0.15],
)

print("""
For magnetic systems, ConvergenceWorkflow automatically:
1. Detects nspin from magnetic_config
2. Runs spin-polarized SCF calculations
3. Tracks energy convergence with magnetization

Example:
    conv = ConvergenceWorkflow(
        atoms=structure,
        pseudopotentials=pseudopotentials,
        magnetic_config={'Mn': [2.0]},  # Initial Mn moment
        ecut_range=[50, 60, 70],
        kspacing_range=[0.3, 0.2]
    )
    
    results = conv.run_convergence_study()
    
The results will show how total energy and magnetic moment
converge with ecutwfc and kspacing.
""")


print("\n" + "="*80)
print("Convergence workflow examples complete!")
print("="*80)
print("""
Key takeaways:

1. Start with ConvergenceWorkflow.run_convergence_study()
   - Tests multiple ecutwfc and kspacing combinations
   - Extracts energy, forces, and k-point info

2. Analyze results with get_recommendations()
   - Suggests 'fast', 'balanced', 'accurate' parameters
   - Based on energy and force convergence

3. Use recommendations in production:
   - Apply optimized parameters to full workflow
   - Balance accuracy needs with computational cost

4. For publication/important results:
   - Use 'accurate' parameters
   - Document convergence tests
   - Include convergence plots
""")
