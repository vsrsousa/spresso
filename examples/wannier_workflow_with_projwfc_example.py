#!/usr/bin/env python
"""
Complete Wannier90 Workflow with PROJWFC Analysis Example

Demonstrates the full pipeline:
- SCF calculation
- Band structure calculation (optional validation)
- NSCF calculation with wavefunction collection
- PROJWFC analysis for orbital contributions
- pw2wannier90 wavefunction conversion
- wannier90 Wannier function generation

The PROJWFC step analyzes projected density of states (PDOS) to identify
which atoms and orbitals dominate the electronic structure, helping guide
the selection of Wannier function projections.
"""

from pathlib import Path
from xespresso.workflow.wannier_workflow import WannierWorkflow

# ============================================================================
# Configuration
# ============================================================================

# Structure and pseudopotentials
CIF_FILE = Path("si.cif")  # Your structure file
PSEUDOS = {
    "Si": "/path/to/Si.pbe.UPF",  # Your pseudopotential files
}

# Wannier parameters
NUM_WANN = 4
PROJECTIONS = "Si: sp3"  # Initial guess for Wannier projections

# K-points
KPTS_SCF = (4, 4, 4)      # Dense for SCF convergence
KPTS_NSCF = (6, 6, 6)     # Even denser for better Wannier interpolation

# Queue/scheduler configuration (optional)
QUEUE = {
    "scheduler": "direct",   # or "slurm", "pbs", etc.
    "execution": "local",    # or "remote"
    "timeout": 3600,
}

# ============================================================================
# Initialize Workflow
# ============================================================================

print("Initializing Wannier Workflow with PROJWFC analysis...")
workflow = WannierWorkflow(
    cif_file=CIF_FILE,
    pseudos=PSEUDOS,
    protocol="moderate",  # fast, moderate, accurate
    num_wann=NUM_WANN,
    projections=PROJECTIONS,
    kpts_scf=KPTS_SCF,
    kpts_nscf=KPTS_NSCF,
    nbnd=None,  # Auto-estimate from pseudopotentials
    spinors=False,
    dis_num_iter=1000,  # Disentanglement iterations
    run_bands=True,  # Enable band structure calculation for validation
    queue=QUEUE,
)

# ============================================================================
# Run Complete Pipeline
# ============================================================================

print("\nStarting complete Wannier workflow pipeline...")
print("This includes: SCF → Bands → NSCF → PROJWFC → pw2wannier90 → wannier90")

results = workflow.run(
    labels={
        'scf': 'runs/01_scf',
        'bands': 'runs/02_bands',
        'nscf': 'runs/03_nscf',
        'projwfc': 'runs/03b_projwfc',
    },
    blocking=True,  # Wait for all jobs to complete
    seedname='si_wannier',
    run_bands_validation=True,  # Compute band structure for comparison
    run_projwfc_analysis=True,   # Compute PDOS for orbital analysis
)

# ============================================================================
# Access Results
# ============================================================================

print("\n" + "="*70)
print("WORKFLOW RESULTS SUMMARY")
print("="*70)

# Get basic info
print(f"\nRun directory: {results['run_dir']}")
print(f"Seedname: {results['seedname']}")
print(f"Projections used: {results['projections_used']}")

# Get calculators
scf_calc = workflow.get_scf_calculator()
nscf_calc = workflow.get_nscf_calculator()
bands_calc = workflow.get_band_structure_calculator()
wannier_dir = workflow.get_wannier_directory()

print(f"\nCalculators:")
print(f"  SCF: {scf_calc.directory}")
print(f"  NSCF: {nscf_calc.directory}")
if bands_calc:
    print(f"  Bands: {bands_calc.directory}")

# ============================================================================
# PROJWFC Analysis Results
# ============================================================================

projwfc_analysis = workflow.get_projwfc_analysis()
if projwfc_analysis:
    print(f"\nPROJWFC Analysis Results:")
    print(f"  Status: {projwfc_analysis['status']}")
    print(f"  PDOS file: {projwfc_analysis['pdos_file']}")
    print(f"  Analysis directory: {projwfc_analysis['run_dir']}")
    
    if projwfc_analysis.get('parsed'):
        print(f"\nParsed PDOS Information:")
        print(f"  {projwfc_analysis['parsed']}")
    
    print(f"\nHow to interpret PDOS results:")
    print(f"  1. Open {Path(projwfc_analysis['pdos_file']).parent}/{Path(projwfc_analysis['pdos_file']).name}*")
    print(f"  2. Identify dominant s, p, d orbital contributions per atom")
    print(f"  3. Use these insights to refine your Wannier projections")
    print(f"  4. Re-run wannier90 with improved projections for better accuracy")
else:
    print("\nNo PROJWFC analysis available. Re-run with run_projwfc_analysis=True")

# ============================================================================
# Band Structure Validation
# ============================================================================

print(f"\nBand Structure Comparison:")
bands_comparison = workflow.compare_bands()
if bands_comparison:
    print(f"  {bands_comparison}")

# ============================================================================
# Wannier Function Quality Validation
# ============================================================================

validation = workflow.validate_wannier_quality()
print(f"\nWannier Function Quality Validation:")
print(f"  Status: {validation['status']}")
print(f"  Wannier centers: {validation['wannier_centers']}")
print(f"  Wannier functions (.xsf): {validation['wannier_xsf']}")
print(f"  Has band structure: {validation['has_band_structure']}")
print(f"  Has PROJWFC analysis: {validation['has_projwfc']}")

print(f"\nRecommendations:")
for i, rec in enumerate(validation['recommendations'], 1):
    print(f"  {i}. {rec}")

# ============================================================================
# Next Steps
# ============================================================================

print("\n" + "="*70)
print("NEXT STEPS")
print("="*70)

print(f"""
1. Review PROJWFC Results:
   - Check {results['run_dir']}/{nscf_calc.prefix}.pdos* files
   - Identify which orbitals dominate at different energies
   - Look for sharp features in PDOS that indicate important bands

2. Refine Wannier Projections (if needed):
   - Based on PDOS analysis, adjust your projections
   - Example: If p orbitals dominate, try projections="Si: sp3d2"
   - Re-run wannier90 with improved projections

3. Validate with Band Structure:
   - Use wannier90 postprocessing to interpolate bands
   - Compare with original DFT band structure
   - Check agreement at important k-points

4. Export Wannier Functions:
   - Wannier centers: {results['run_dir']}/si_wannier_centres.xyz
   - Orbital density: {results['run_dir']}/si_wannier_*.xsf
   - Use for: Electronic transport, optical properties, etc.

5. Further Analysis:
   - Compute Hamiltonian matrix elements with wannier90
   - Interface with other codes (e.g., EPW for e-ph coupling)
   - Generate tight-binding models

For more information on Wannier90, visit: http://www.wannier.org
For PROJWFC documentation: https://www.quantum-espresso.org/
""")

print("="*70)
