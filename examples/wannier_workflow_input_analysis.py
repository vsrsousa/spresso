#!/usr/bin/env python
"""
Quick visual inspection of WannierWorkflow input generation.

Shows all input files that would be created before running calculations.
No actual DFT or Wannier90 calculations are executed.
"""

from xespresso.workflow.wannier_workflow import (
    generate_pw2wannier_input,
    generate_seedname_win,
)


def show_input_example(title, content):
    """Pretty print an input file example."""
    print(f"\n{'─'*80}")
    print(f"{title}")
    print(f"{'─'*80}")
    print(content)
    print(f"{'─'*80}\n")


def main():
    print("\n" + "="*80)
    print("WANNIER WORKFLOW - INPUT FILES INSPECTION")
    print("="*80)
    
    # =========================================================================
    # PROJWFC Input
    # =========================================================================
    projwfc_input = """&inputpp
  prefix = 'pwscf'
  outdir = './'
/
filpdos = 'pwscf.pdos'
"""
    show_input_example("📄 projwfc.in - Projection Analysis Input", projwfc_input)
    
    # =========================================================================
    # PW2WANNIER90 Input
    # =========================================================================
    pw2wan_input = generate_pw2wannier_input('pwscf')
    show_input_example("📄 pw2wannier.in - Wavefunction Conversion Input", pw2wan_input)
    
    # =========================================================================
    # WANNIER90 Input (Example 1: Simple)
    # =========================================================================
    win_simple = generate_seedname_win(
        num_wann=4,
        projections='Si: sp3',
        spinors=False,
        dis_num_iter=1000
    )
    show_input_example("📄 wannier_seed.win - Wannier90 Input (Si: 4 functions, sp3)", win_simple)
    
    # =========================================================================
    # WANNIER90 Input (Example 2: More complex)
    # =========================================================================
    win_complex = generate_seedname_win(
        num_wann=8,
        projections='Fe: d; O: p',
        spinors=False,
        dis_num_iter=1500
    )
    show_input_example("📄 wannier_seed.win - Wannier90 Input (Fe-O compound, 8 functions)", win_complex)
    
    # =========================================================================
    # Analysis
    # =========================================================================
    print("\n" + "="*80)
    print("ANALYSIS POINTS FOR IMPROVEMENTS")
    print("="*80)
    
    improvements = {
        "projwfc.in": [
            "✓ CURRENT: Minimal input, outputs all projections",
            "? POTENTIAL IMPROVEMENTS:",
            "  1. Add Fermi energy calculation flag",
            "  2. Add dos_proj_type='atomic' for atomic projections",
            "  3. Add energy window specification for PDOS",
            "  4. Consider DeltaE parameter for energy resolution",
        ],
        "pw2wannier.in": [
            "✓ CURRENT: Write amn, mmn, eig, unk",
            "? POTENTIAL IMPROVEMENTS:",
            "  1. Add write_u_matrices = .true. for U matrices",
            "  2. Add scdm_entanglement = 'erfc' for entanglement method",
            "  3. Add scdm_mu and scdm_sigma for entanglement control",
            "  4. Consider adding phonon interpolation flags",
        ],
        "wannier_seed.win": [
            "✓ CURRENT: Basic setup with num_wann, spinors, projections",
            "? POTENTIAL IMPROVEMENTS:",
            "  1. Add 'auto_projections = .true.' for auto-generated projections",
            "  2. Add search_shells directive for automatic band grouping",
            "  3. Add Gaussian smearing parameters",
            "  4. Add plotting keywords (wannier_plot = .true.)",
            "  5. Add band structure interpolation keywords",
            "  6. Add hr_plot keyword for Hamiltonian visualization",
            "  7. Consider timing and convergence control keywords",
        ],
        "General Workflow": [
            "✓ CURRENT: SCF → Bands → PROJWFC → NSCF → pw2wannier → wannier90",
            "? POTENTIAL IMPROVEMENTS:",
            "  1. Add checkpoint/restart mechanism for failed stages",
            "  2. Add automatic backup of important results",
            "  3. Add convergence diagnostics between stages",
            "  4. Add preliminary band structure checks",
            "  5. Add PDOS analysis for projection validation",
            "  6. Add interpolated band structure visualization",
        ],
    }
    
    for section, points in improvements.items():
        print(f"\n{section}:")
        for point in points:
            print(f"  {point}")
    
    print("\n" + "="*80)
    print("EXAMPLE USAGE:")
    print("="*80)
    
    usage = """
# Create workflow instance (no execution)
wf = WannierWorkflow(
    cif_file='structure.cif',
    pseudopotentials={'Si': 'Si.pbe.UPF'},
    protocol='moderate',
    num_wann=4,
    projections='Si: sp3'
)

# Run with all optional stages
results = wf.run(
    blocking=True,
    run_bands_validation=True,
    run_projwfc_analysis=True,
    dry_run=True
)

# Access results
print(results['wannier90']['outputs'])
wf.compare_bands()
wf.validate_wannier_quality()
"""
    print(usage)
    
    print("="*80 + "\n")


if __name__ == '__main__':
    main()
