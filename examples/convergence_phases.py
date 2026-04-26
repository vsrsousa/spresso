"""
Example: Run convergence with independent phases (ecut only, kpt only, or both)

Demonstrates the 3 ways to use the phases parameter:
1. run_convergence(phases='ecut') - PHASE 1 only
2. run_convergence(phases='kpt') - PHASE 2 only
3. run_convergence(phases='both') - Both phases (default)
"""

from ase.build import bulk
from xespresso.workflow import ConvergenceWorkflow


def example_only_ecutwfc():
    """
    Example 1: Run only ecutwfc convergence (PHASE 1)
    
    Use when:
    - You want to study ecutwfc convergence separately
    - You need a convergence plot for a paper
    - You want to establish ecutwfc first, then optimize kspacing later
    """
    print("=" * 80)
    print("EXAMPLE 1: Convergence of ecutwfc only")
    print("=" * 80)
    
    # Create structure
    atoms = bulk('Au', 'fcc', a=4.08)
    
    # Initialize workflow
    wf = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials_config='SSSP_efficiency',
        precision='low'
    )
    
    # Run ONLY PHASE 1 (ecutwfc convergence)
    print("\nRunning PHASE 1 (ecutwfc convergence)...")
    results_ecut = wf.run_convergence(
        max_ecutwfc=100.0,
        ecut_step=10.0,
        verbose=True,
        phases='ecut'  # ← Key: phases parameter
    )
    
    # Check results
    print(f"\nOptimal ecutwfc: {wf.optimal_ecutwfc:.1f} Ry")
    print(f"Results shape: {results_ecut.shape}")
    print(f"Columns: {list(results_ecut.columns)}")
    
    # Now you can plot ecutwfc vs energy
    # plt.plot(results_ecut['ecutwfc'], results_ecut['energy'])
    # plt.xlabel('Ecutwfc (Ry)')
    # plt.ylabel('Energy (eV)')
    # plt.savefig('ecutwfc_convergence.png')
    
    return wf, results_ecut


def example_only_kspacing(optimal_ecutwfc=50.0):
    """
    Example 2: Run only kspacing convergence (PHASE 2)
    
    Use when:
    - ecutwfc is already optimized (from a previous run)
    - You want to optimize only k-point density
    - You want to study kspacing convergence separately
    """
    print("\n" + "=" * 80)
    print("EXAMPLE 2: Convergence of kspacing only")
    print("=" * 80)
    
    # Create structure
    atoms = bulk('Au', 'fcc', a=4.08)
    
    # Initialize workflow
    wf = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials_config='SSSP_efficiency',
        precision='low'
    )
    
    # Set pre-computed ecutwfc (e.g., from PHASE 1)
    print(f"\nUsing fixed ecutwfc: {optimal_ecutwfc:.1f} Ry")
    wf.optimal_ecutwfc = optimal_ecutwfc
    
    # Run ONLY PHASE 2 (kspacing convergence)
    print("Running PHASE 2 (kspacing convergence)...")
    results_kpt = wf.run_convergence(
        min_kspacing_allowed=0.05,
        kspacing_step=0.05,
        verbose=True,
        phases='kpt'  # ← Key: phases parameter
    )
    
    # Check results
    print(f"\nOptimal kspacing: {wf.optimal_kspacing:.3f} Å⁻¹")
    print(f"Results shape: {results_kpt.shape}")
    print(f"Columns: {list(results_kpt.columns)}")
    
    return wf, results_kpt


def example_both_phases():
    """
    Example 3: Run both phases sequentially (PHASE 1 + PHASE 2)
    
    Use when:
    - You want a complete convergence study
    - Both parameters need optimization
    - You want the classic two-phase independent algorithm
    """
    print("\n" + "=" * 80)
    print("EXAMPLE 3: Full convergence (both phases)")
    print("=" * 80)
    
    # Create structure
    atoms = bulk('Au', 'fcc', a=4.08)
    
    # Initialize workflow
    wf = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials_config='SSSP_efficiency',
        precision='low'
    )
    
    # Run both PHASE 1 + PHASE 2 sequentially
    print("\nRunning full convergence (PHASE 1 + PHASE 2)...")
    results = wf.run_convergence(
        max_ecutwfc=100.0,
        ecut_step=10.0,
        kspacing_step=0.05,
        verbose=True,
        phases='both'  # ← Key: phases='both' (or omit for default)
    )
    
    # Check results
    print(f"\nOptimal ecutwfc: {wf.optimal_ecutwfc:.1f} Ry")
    print(f"Optimal kspacing: {wf.optimal_kspacing:.3f} Å⁻¹")
    print(f"Results shape: {results.shape}")
    
    recommendations = wf.get_recommendations(verbose=True)
    print(f"\nRecommendations: {recommendations}")
    
    return wf, results


def example_using_phases_parameter():
    """
    Example 4: Different ways to specify phases parameter
    """
    print("\n" + "=" * 80)
    print("EXAMPLE 4: Different ways to use phases parameter")
    print("=" * 80)
    
    atoms = bulk('Au', 'fcc', a=4.08)
    
    # Way 1: Using phases='ecut' for PHASE 1 only
    print("\nWay 1: phases='ecut' (PHASE 1 only)")
    wf1 = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials_config='SSSP_efficiency',
        precision='low'
    )
    results1 = wf1.run_convergence(phases='ecut')
    print(f"  → Optimal ecutwfc: {wf1.optimal_ecutwfc:.1f} Ry")
    
    # Way 2: Using phases='kpt' for PHASE 2 only
    print("\nWay 2: phases='kpt' (PHASE 2 only)")
    wf2 = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials_config='SSSP_efficiency',
        precision='low'
    )
    wf2.optimal_ecutwfc = wf1.optimal_ecutwfc  # Use ecutwfc from PHASE 1
    results2 = wf2.run_convergence(phases='kpt')
    print(f"  → Optimal kspacing: {wf2.optimal_kspacing:.3f} Å⁻¹")
    
    # Way 3: Using phases='both' (default)
    print("\nWay 3: phases='both' (both phases, default)")
    wf3 = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials_config='SSSP_efficiency',
        precision='low'
    )
    results3 = wf3.run_convergence(phases='both')
    print(f"  → Optimal ecutwfc: {wf3.optimal_ecutwfc:.1f} Ry")
    print(f"  → Optimal kspacing: {wf3.optimal_kspacing:.3f} Å⁻¹")
    
    # Way 4: Omitting phases (uses default 'both')
    print("\nWay 4: Omitting phases (same as phases='both')")
    wf4 = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials_config='SSSP_efficiency',
        precision='low'
    )
    results4 = wf4.run_convergence()  # Default: phases='both'
    print(f"  → Optimal ecutwfc: {wf4.optimal_ecutwfc:.1f} Ry")
    print(f"  → Optimal kspacing: {wf4.optimal_kspacing:.3f} Å⁻¹")


def example_workflow_with_phases():
    """
    Example 5: Typical workflow with two separate phases
    
    This shows the most practical use case:
    1. First optimize ecutwfc with low precision
    2. Then optimize kspacing with higher precision
    """
    print("\n" + "=" * 80)
    print("EXAMPLE 5: Typical workflow (2 phases separately)")
    print("=" * 80)
    
    atoms = bulk('Au', 'fcc', a=4.08)
    
    # PHASE 1: Optimize ecutwfc with low precision
    print("\n--- PHASE 1: Optimize ecutwfc (low precision) ---")
    wf1 = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials_config='SSSP_efficiency',
        precision='low'  # Fast, coarse convergence
    )
    results_ecut = wf1.run_convergence(
        max_ecutwfc=100.0,
        ecut_step=10.0,
        phases='ecut'  # ← Only PHASE 1
    )
    print(f"Found ecutwfc: {wf1.optimal_ecutwfc:.1f} Ry")
    
    # PHASE 2: Optimize kspacing with medium precision
    print("\n--- PHASE 2: Optimize kspacing (medium precision) ---")
    wf2 = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials_config='SSSP_efficiency',
        precision='medium'  # Better convergence criteria
    )
    wf2.optimal_ecutwfc = wf1.optimal_ecutwfc  # Use ecutwfc from PHASE 1
    results_kpt = wf2.run_convergence(
        phases='kpt'  # ← Only PHASE 2
    )
    print(f"Found kspacing: {wf2.optimal_kspacing:.3f} Å⁻¹")
    
    # Summary
    print("\n" + "=" * 80)
    print("CONVERGENCE COMPLETE")
    print("=" * 80)
    print(f"Optimal ecutwfc: {wf1.optimal_ecutwfc:.1f} Ry")
    print(f"Optimal kspacing: {wf2.optimal_kspacing:.3f} Å⁻¹")
    print(f"Expected time: ~2-3 hours total (vs ~6+ hours for nested loops)")


if __name__ == '__main__':
    # Run examples
    print("Convergence Phases Examples")
    print("=" * 80)
    
    # NOTE: These are example structures. To actually run them,
    # you would need:
    # - A valid QE installation
    # - Proper pseudopotentials
    # - A queue configuration
    
    # Uncomment to run:
    # wf1, results_ecut = example_only_ecutwfc()
    # wf2, results_kpt = example_only_kspacing(optimal_ecutwfc=50.0)
    # wf3, results = example_both_phases()
    # example_using_phases_parameter()
    # example_workflow_with_phases()
    
    print("\nTo run these examples:")
    print("1. Uncomment the function calls at the bottom of this file")
    print("2. Make sure you have QE installed and configured")
    print("3. python examples/convergence_phases.py")
