"""
Example script showing how to orchestrate a complete Wannierization pipeline using
the WannierWorkflow class with band structure validation.

This example demonstrates the complete workflow including:
  1. SCF
  2. Band Structure (for quality validation)
  3. NSCF (with wavefunction collection)
  4. pw2wannier90
  5. wannier90

Band structure calculation is essential for validating Wannier function quality
by comparing DFT band structure with Wannier-interpolated bands.

Adjust paths, pseudopotential mapping and machine configuration as needed.
"""
from pathlib import Path
from xespresso.workflow import WannierWorkflow


def main():
    """Simple Wannier workflow with band structure validation."""
    
    # Example inputs (user must adapt)
    cif_file = "structure.cif"
    pseudos = {"Si": "/path/to/Si.pbe.UPF"}  # adjust
    protocol = "moderate"
    queue = None  # use local direct scheduler by default
    
    # Create complete Wannier workflow WITH band structure for validation
    wf = WannierWorkflow(
        cif_file=cif_file,
        pseudos=pseudos,
        protocol=protocol,
        num_wann=4,                    # Generate 4 Wannier functions
        projections="Si: sp3",         # Initial projections
        kpts_scf=(4, 4, 4),           # SCF k-point mesh
        kpts_nscf=(6, 6, 6),          # NSCF k-point mesh (denser)
        run_bands=True,               # NEW: Calculate band structure for validation!
        queue=queue,                  # Use local execution
    )
    
    # Execute the complete pipeline
    # This runs: SCF → Band Structure → NSCF (wf_collect) → pw2wannier → wannier90
    results = wf.run(
        blocking=True,
        run_bands_validation=True  # Explicitly enable band structure validation
    )
    
    # Access results
    print("\n" + "="*70)
    print("RESULTS & VALIDATION")
    print("="*70)
    print(f"SCF calculator: {wf.get_scf_calculator()}")
    print(f"Band structure calculator: {wf.get_band_structure_calculator()}")
    print(f"NSCF calculator: {wf.get_nscf_calculator()}")
    print(f"Wannier directory: {wf.get_wannier_directory()}")
    print(f"Seedname: {wf.get_seedname()}")
    
    # Validate Wannier function quality
    validation = wf.validate_wannier_quality()
    print(f"\nWannier Quality Validation:")
    print(f"  Status: {validation['status']}")
    print(f"  Wannier centers: {validation['wannier_centers']}")
    print(f"  Wannier XSF: {validation['wannier_xsf']}")
    if validation.get('recommendations'):
        print(f"  Recommendations:")
        for rec in validation['recommendations']:
            print(f"    - {rec}")
    
    # Compare band structures
    if wf.get_band_structure_calculator() is not None:
        print(f"\nBand Structure Comparison (DFT vs Wannier):")
        comparison = wf.compare_bands(verbose=True)
        print(f"  DFT bands: {comparison.get('dft_bands')}")
        print(f"  Wannier interpolation: {comparison.get('wannier_bands_available_in')}")
    
    print("="*70 + "\n")


if __name__ == "__main__":
    main()


