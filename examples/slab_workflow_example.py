"""
Example: Complete slab workflow for Au surfaces.

This example shows how to:
1. Run bulk convergence for Au
2. Generate Au(100), Au(110), Au(111) slabs
3. Relax each slab with optimized parameters
4. Calculate surface energies

Run with:
    python examples/slab_workflow_example.py
"""

from ase.build import bulk
from xespresso.workflow.slab_workflow import SlabWorkflow


def main():
    # =========================================================================
    # Setup
    # =========================================================================
    print("\n" + "="*70)
    print("Au Surface Slab Workflow Example")
    print("="*70)
    
    # Create bulk Au structure (FCC, a=4.0782 Å)
    bulk_au = bulk('Au', 'fcc', a=4.0782)
    print(f"\nBulk structure: {bulk_au.get_chemical_formula()}")
    print(f"Lattice parameter: a = {bulk_au.cell[0,0]:.4f} Å")
    
    # =========================================================================
    # Create workflow
    # =========================================================================
    print("\nInitializing SlabWorkflow...")
    
    wf = SlabWorkflow(
        bulk_atoms=bulk_au,
        surface_indices=[
            (1, 0, 0),  # FCC 100 surface
            (1, 1, 0),  # FCC 110 surface  
            (1, 1, 1),  # FCC 111 surface
        ],
        pseudopotentials_config='default',
        machine='medusa',  # or None for local execution
        code_version='7.4.1',
        protocol='moderate',
        nlayers=4,
        fix_layers=[0, 1],  # Fix first 2 layers
    )
    
    print(f"Slab parameters:")
    print(f"  Number of layers: {wf.nlayers}")
    print(f"  Fixed layers: {wf.fix_layers}")
    print(f"  Min slab size: {wf.min_slab_size} Å")
    print(f"  Min vacuum: {wf.min_vacuum_size} Å")
    
    # =========================================================================
    # Phase 1: Bulk Convergence
    # =========================================================================
    print("\n" + "="*70)
    print("PHASE 1: Bulk Convergence")
    print("="*70)
    
    bulk_rec = wf.run_bulk_convergence(
        label_prefix='au_bulk_conv',
        precision='low',  # Use 'normal' or 'high' for production
        verbose=True
    )
    
    print(f"\nBulk convergence results:")
    print(f"  Optimal ecutwfc: {bulk_rec['optimal_ecutwfc']:.1f} Ry")
    print(f"  Optimal kspacing: {bulk_rec['optimal_kspacing']:.3f} Å⁻¹")
    
    # =========================================================================
    # Phase 2: Generate Slabs
    # =========================================================================
    print("\n" + "="*70)
    print("PHASE 2: Generate Slabs")
    print("="*70)
    
    slabs = wf.generate_slabs()
    
    print(f"\nGenerated {len(slabs)} slabs:")
    for hkl, slab_data in slabs.items():
        atoms = slab_data['atoms']
        print(f"  {hkl}: {len(atoms)} atoms, cell: {atoms.cell.lengths()}")
    
    # =========================================================================
    # Phase 3: Relax Slabs
    # =========================================================================
    print("\n" + "="*70)
    print("PHASE 3: Relax Slabs")
    print("="*70)
    
    relaxed = wf.relax_slabs(
        label_prefix='au_relax',
        relax_type='vc-relax',
        verbose=True
    )
    
    print(f"\nRelaxed {len(relaxed)} slabs")
    
    # =========================================================================
    # Phase 4: Surface Energies
    # =========================================================================
    print("\n" + "="*70)
    print("PHASE 4: Surface Energies")
    print("="*70)
    
    # Get bulk energy from convergence (or use DFT result)
    bulk_energy_per_atom = bulk_rec.get('bulk_energy_per_atom', -3.8)  # eV/atom
    
    surface_energies = wf.calculate_surface_energies(
        bulk_energy_per_atom=bulk_energy_per_atom
    )
    
    print(f"\nSurface Energies (J/m²):")
    for hkl, gamma in sorted(surface_energies.items()):
        print(f"  Au{hkl}: γ = {gamma:.3f} J/m²")
    
    # Experimental reference (at 0 K)
    print(f"\nExperimental Reference (at 0 K):")
    experimental = {
        (1, 0, 0): 1.24,  # Au(100) in J/m²
        (1, 1, 0): 1.38,  # Au(110)
        (1, 1, 1): 1.14,  # Au(111)
    }
    print("  Au(100): 1.24 J/m²")
    print("  Au(110): 1.38 J/m²")
    print("  Au(111): 1.14 J/m²")
    
    # =========================================================================
    # Save Results
    # =========================================================================
    print("\n" + "="*70)
    print("Saving Results")
    print("="*70)
    
    wf.save_results(output_dir='./')
    
    print("\n" + "="*70)
    print("✓ Workflow Complete")
    print("="*70)
    
    return wf


if __name__ == '__main__':
    wf = main()
