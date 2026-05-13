#!/usr/bin/env python
"""
Inspect WannierWorkflow inputs WITHOUT executing calculations.

This script demonstrates what input files would be generated and allows
analysis of the workflow before actually running it.

Usage:
    python examples/wannier_workflow_inspect_inputs.py
"""

from pathlib import Path
import tempfile
from xespresso.workflow.wannier_workflow import (
    WannierWorkflow,
    generate_pw2wannier_input,
    generate_seedname_win,
    suggest_nbnd_from_pseudos,
)


def create_test_cif():
    """Create a minimal test structure (Si diamond cubic)."""
    from ase.build import bulk
    from ase.io import read, write
    
    # Create Si diamond structure
    si = bulk('Si', 'diamond', a=5.431)
    
    # Write to temp CIF
    temp_cif = Path(tempfile.gettempdir()) / "si_test.cif"
    write(str(temp_cif), si)
    print(f"✓ Created test structure: {temp_cif}")
    return str(temp_cif)


def inspect_wannier_inputs(
    num_wann=4,
    projections='Si: sp3',
    protocol='fast',
    spinors=False,
):
    """
    Inspect all input files that would be generated.
    
    Parameters:
        num_wann: Number of Wannier functions
        projections: Initial projections string
        protocol: Convergence protocol
        spinors: Use spinors or not
    """
    print("\n" + "="*80)
    print("WANNIER WORKFLOW - INPUT INSPECTION (NO EXECUTION)")
    print("="*80)
    
    # Create test structure
    cif_file = create_test_cif()
    
    # Define dummy pseudopotentials (paths don't need to exist for inspection)
    pseudos = {
        'Si': '/path/to/Si.pbe-n-kjpaw_psl.1.0.0.UPF'
    }
    
    print(f"\n{'─'*80}")
    print("1. WORKFLOW INITIALIZATION")
    print(f"{'─'*80}")
    
    # Initialize workflow (this doesn't run calculations, just sets up)
    try:
        wf = WannierWorkflow(
            cif_file=cif_file,
            pseudopotentials=pseudopotentials,
            protocol=protocol,
            num_wann=num_wann,
            projections=projections,
            spinors=spinors,
            kpts_scf=(4, 4, 4),
            kpts_nscf=(6, 6, 6),
        )
        print("✓ WannierWorkflow initialized successfully")
    except Exception as e:
        print(f"✗ Error during initialization: {e}")
        return
    
    # Print initialization parameters
    print(f"\nWorkflow Parameters:")
    print(f"  Structure: {wf.cif_file}")
    print(f"  Protocol: {wf.protocol}")
    print(f"  Num Wannier: {wf.num_wann}")
    print(f"  Projections: {wf.projections}")
    print(f"  K-points (SCF): {wf.kpts_scf}")
    print(f"  K-points (NSCF): {wf.kpts_nscf}")
    print(f"  Num bands: {wf.nbnd}")
    print(f"  Spinors: {wf.spinors}")
    
    # =========================================================================
    print(f"\n{'─'*80}")
    print("2. PROJWFC INPUT FILE")
    print(f"{'─'*80}")
    
    # Generate projwfc input
    projwfc_input = f"""&inputpp
  prefix = 'pwscf'
  outdir = './'
/
filpdos = 'pwscf.pdos'
"""
    print("projwfc.in:")
    print("─" * 40)
    print(projwfc_input)
    print("─" * 40)
    
    # =========================================================================
    print(f"\n{'─'*80}")
    print("3. PW2WANNIER90 INPUT FILE")
    print(f"{'─'*80}")
    
    pw2wan_input = generate_pw2wannier_input('pwscf')
    print("pw2wannier.in:")
    print("─" * 40)
    print(pw2wan_input)
    print("─" * 40)
    
    # =========================================================================
    print(f"\n{'─'*80}")
    print("4. WANNIER90 INPUT FILE (.win)")
    print(f"{'─'*80}")
    
    seedname_win = generate_seedname_win(
        num_wann=num_wann,
        projections=projections,
        spinors=spinors,
        dis_num_iter=1000
    )
    print("wannier_seed.win:")
    print("─" * 40)
    print(seedname_win)
    print("─" * 40)
    
    # =========================================================================
    print(f"\n{'─'*80}")
    print("5. PSEUDOPOTENTIAL ANALYSIS")
    print(f"{'─'*80}")
    
    nbnd_suggested = suggest_nbnd_from_pseudos(pseudos, buffer=10)
    print(f"Suggested nbnd (num_bands): {nbnd_suggested}")
    print(f"  Calculation: nbnd = max(64, int(valence_electrons/2) + buffer)")
    print(f"  Buffer: 10 bands")
    print(f"  Valence analysis: Using provided pseudopotential paths")
    
    # =========================================================================
    print(f"\n{'─'*80}")
    print("6. PIPELINE STAGES (WITHOUT EXECUTION)")
    print(f"{'─'*80}")
    
    total_stages = 5  # Base
    if True:  # run_bands_validation
        total_stages += 1
    if True:  # run_projwfc_analysis
        total_stages += 1
    
    stages = [
        f"STAGE 1: SCF Calculation",
        f"STAGE 2: Band Structure Calculation (optional, validation)",
        f"STAGE 3: PROJWFC Analysis (projection analysis for Wannier guidance)",
        f"STAGE 4: NSCF Calculation (with wavefunction collection)",
        f"STAGE 5: pw2wannier90 (wavefunction conversion)",
        f"STAGE 6: wannier90 (Wannier function generation)",
    ]
    
    for stage in stages:
        print(f"  ✓ {stage}")
    
    # =========================================================================
    print(f"\n{'─'*80}")
    print("7. OUTPUT FILES THAT WOULD BE GENERATED")
    print(f"{'─'*80}")
    
    outputs = {
        'PROJWFC outputs': [
            'pwscf.pdos',
            'pwscf.pdos.up',
            'pwscf.pdos.dw',
            'pwscf.pdos.txt',
        ],
        'pw2wannier90 outputs': [
            'wannier_seed.amn',
            'wannier_seed.mmn',
            'wannier_seed.eig',
            'wannier_seed.unk',
        ],
        'wannier90 outputs': [
            'wannier_seed.wout',
            'wannier_seed_centres.xyz',
            'wannier_seed_00001.xsf',
            'wannier_seed_00002.xsf',
            'wannier_seed_00003.xsf',
            'wannier_seed_00004.xsf',
            'wannier_seed.win',
        ],
    }
    
    for category, files in outputs.items():
        print(f"\n{category}:")
        for f in files:
            print(f"  - {f}")
    
    # =========================================================================
    print(f"\n{'─'*80}")
    print("8. VALIDATION & ANALYSIS TOOLS")
    print(f"{'─'*80}")
    
    print("\nAvailable methods after run():")
    print("  ✓ wf.get_scf_calculator()")
    print("  ✓ wf.get_nscf_calculator()")
    print("  ✓ wf.get_band_structure_calculator()")
    print("  ✓ wf.get_projwfc_analysis()")
    print("  ✓ wf.compare_bands()")
    print("  ✓ wf.validate_wannier_quality()")
    
    # =========================================================================
    print(f"\n{'─'*80}")
    print("9. KEY PARAMETERS ANALYSIS")
    print(f"{'─'*80}")
    
    analysis = {
        'Convergence Protocol': {
            'fast': 'Lower cutoffs, coarse k-points (testing)',
            'moderate': 'Balanced convergence (default)',
            'accurate': 'High cutoffs, dense k-points (production)',
        },
        'K-points Impact': {
            'SCF': 'Lower density acceptable (4,4,4)',
            'NSCF': 'Must be denser (6,6,6) for better interpolation',
            'Bands': 'High-symmetry path (automatically generated)',
        },
        'Wannier Functions': {
            'num_wann': f'{num_wann} functions requested',
            'projections': projections,
            'dis_num_iter': '1000 iterations',
        },
    }
    
    for category, items in analysis.items():
        print(f"\n{category}:")
        for key, value in items.items():
            print(f"  {key}: {value}")
    
    print("\n" + "="*80)
    print("✓ INPUT INSPECTION COMPLETE - Ready for analysis")
    print("="*80 + "\n")


if __name__ == '__main__':
    # Example 1: Standard Si with sp3 projections
    print("\n" + "🔍 EXAMPLE 1: Si with sp3 projections")
    inspect_wannier_inputs(
        num_wann=4,
        projections='Si: sp3',
        protocol='fast',
    )
    
    # Example 2: More Wannier functions with p projections
    print("\n" + "🔍 EXAMPLE 2: Si with 8 Wannier functions and p projections")
    inspect_wannier_inputs(
        num_wann=8,
        projections='Si: p',
        protocol='moderate',
    )
