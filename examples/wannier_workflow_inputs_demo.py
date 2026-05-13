#!/usr/bin/env python
"""
Direct input file generation - NO IMPORTS needed, pure example output.

This shows exactly what inputs would be generated for the Wannier workflow.
Shows all 6 input files: SCF, Bands, PROJWFC, NSCF, pw2wannier, wannier90
"""


def main():
    print("\n" + "="*80)
    print("WANNIER WORKFLOW - ALL INPUT FILES (6 STAGES)")
    print("="*80)
    
    # =========================================================================
    # STAGE 1: SCF Input
    # =========================================================================
    scf_input = """&CONTROL
  calculation = 'scf'
  prefix = 'pwscf'
  outdir = './'
  pseudo_dir = './'
  tstress = .true.
  tprnfor = .true.
/
&SYSTEM
  ibrav = 2
  celldm(1) = 10.26
  nat = 2
  ntyp = 1
  ecutwfc = 60.0
  ecutrho = 240.0
/
&ELECTRONS
  conv_thr = 1.0e-08
  mixing_beta = 0.7
/
ATOMIC_SPECIES
 Si  28.086  Si.pbe-n-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS (alat)
 Si  0.0  0.0  0.0
 Si  0.25  0.25  0.25

K_POINTS (automatic)
 4 4 4 0 0 0
"""
    
    print("\n" + "─"*80)
    print("STAGE 1/6: SCF Input (pwscf.scf)")
    print("─"*80)
    print(scf_input)
    
    # =========================================================================
    # STAGE 2: BANDS Input
    # =========================================================================
    bands_input = """&CONTROL
  calculation = 'bands'
  prefix = 'pwscf'
  outdir = './'
  pseudo_dir = './'
  tstress = .false.
  tprnfor = .false.
/
&SYSTEM
  ibrav = 2
  celldm(1) = 10.26
  nat = 2
  ntyp = 1
  ecutwfc = 60.0
  ecutrho = 240.0
  nbnd = 64
/
&ELECTRONS
  conv_thr = 1.0e-08
  diago_full_acc = .true.
/
ATOMIC_SPECIES
 Si  28.086  Si.pbe-n-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS (alat)
 Si  0.0  0.0  0.0
 Si  0.25  0.25  0.25

K_POINTS (crystal_b)
 5
  0.5000  0.5000  0.5000  20   ! X
  0.0000  0.0000  0.0000  20   ! Gamma
  0.5000  0.5000  0.5000  20   ! X
  0.5000  0.2500  0.7500  20   ! W
  0.0000  0.0000  0.0000  1    ! Gamma
"""
    
    print("─"*80)
    print("STAGE 2/6: BANDS Input (pwscf.bands)")
    print("─"*80)
    print(bands_input)
    
    # =========================================================================
    # STAGE 3: PROJWFC Input
    # =========================================================================
    projwfc_input = """&inputpp
  prefix = 'pwscf'
  outdir = './'
/
filpdos = 'pwscf.pdos'
"""
    
    print("─"*80)
    print("STAGE 3/6: PROJWFC Input (projwfc.in)")
    print("─"*80)
    print(projwfc_input)
    
    # =========================================================================
    # STAGE 4: NSCF Input
    # =========================================================================
    nscf_input = """&CONTROL
  calculation = 'nscf'
  prefix = 'pwscf'
  outdir = './'
  pseudo_dir = './'
  tstress = .false.
  tprnfor = .false.
/
&SYSTEM
  ibrav = 2
  celldm(1) = 10.26
  nat = 2
  ntyp = 1
  ecutwfc = 60.0
  ecutrho = 240.0
  nbnd = 64
/
&ELECTRONS
  conv_thr = 1.0e-08
  diago_full_acc = .true.
/
ATOMIC_SPECIES
 Si  28.086  Si.pbe-n-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS (alat)
 Si  0.0  0.0  0.0
 Si  0.25  0.25  0.25

K_POINTS (automatic)
 6 6 6 0 0 0
"""
    
    print("─"*80)
    print("STAGE 4/6: NSCF Input (pwscf.nscf)")
    print("─"*80)
    print(nscf_input)
    
    # =========================================================================
    # STAGE 5: PW2WANNIER90 Input
    # =========================================================================
    pw2wan_input = """&inputpp
  prefix = 'pwscf'
  outdir = './'
/
write_amn = .true.
write_mmn = .true.
write_unk = .true.
"""
    
    print("─"*80)
    print("STAGE 5/6: PW2WANNIER90 Input (pw2wannier.in)")
    print("─"*80)
    print(pw2wan_input)
    
    # =========================================================================
    # STAGE 6: WANNIER90 Input (Example 1: Simple)
    # =========================================================================
    win_simple = """num_wann = 4
dis_num_iter = 1000
spinors = .false.

begin projections
Si: sp3
end projections
"""
    
    print("─"*80)
    print("STAGE 6/6: WANNIER90 Input (wannier_seed.win)")
    print("─"*80)
    print(win_simple)
    
    
    print("\n" + "="*80)
    print("FILE SUMMARY - ALL 6 INPUTS CREATED BY WORKFLOW")
    print("="*80)
    
    files = {
        "1. pwscf.scf": "SCF Calculation (convergence de densidade)",
        "2. pwscf.bands": "Band Structure (validação da qualidade)",
        "3. projwfc.in": "PROJWFC (análise de projeções)",
        "4. pwscf.nscf": "NSCF Calculation (coleção de funções de onda)",
        "5. pw2wannier.in": "Conversão de funções de onda",
        "6. wannier_seed.win": "Wannier90 (geração de funções Wannier)",
    }
    
    for name, description in files.items():
        print(f"\n{name}")
        print(f"  → {description}")
    
    # =========================================================================
    # Analysis
    # =========================================================================
    print("\n" + "="*80)
    print("IMPROVEMENT ANALYSIS - What Could Be Better")
    print("="*80)
    
    improvements = {
        "1. projwfc.in": {
            "Current": "Minimal, just outputs PDOS",
            "Improvements": [
                "• Add Fermi energy calculation flag",
                "• Add dos_proj_type='atomic' for atomic projections",
                "• Add ngauss parameter for smearing method",
                "• Add energy window specification (emin, emax)",
                "• Add DeltaE for energy grid resolution",
            ]
        },
        "2. pw2wannier.in": {
            "Current": "Writes amn, mmn, eig, unk files",
            "Improvements": [
                "• Add write_u_matrices = .true. for U matrices",
                "• Add scdm_entanglement keywords for entanglement control",
                "• Add kmesh_tol parameter for k-point mapping",
                "• Consider adding spin_axis for magnetic systems",
                "• Add wf_file for reading wavefunctions from disk",
            ]
        },
        "3. wannier_seed.win": {
            "Current": "Basic: num_wann, spinors, projections",
            "Improvements": [
                "• Add 'auto_projections = .true.' for automatic guesses",
                "• Add 'search_shells' for automatic band search",
                "• Add smearing parameters (sigma_T, sigma_G)",
                "• Add 'wannier_plot = .true.' for density visualization",
                "• Add 'bands_plot = .true.' for band structure",
                "• Add 'hr_plot = .true.' for Hamiltonian visualization",
                "• Add 'translate_home_cell = .true.' for physical centers",
                "• Add convergence control (conv_tol, conv_window)",
            ]
        },
        "4. Pipeline & Workflow": {
            "Current": "Sequential: SCF → Bands → PROJWFC → NSCF → pw2wan → w90",
            "Improvements": [
                "• Add checkpoint/restart for failed stages",
                "• Add automatic validation between stages",
                "• Add PDOS analysis for projection sanity checks",
                "• Add band structure pre-check before Wannier",
                "• Add interpolated band structure comparison",
                "• Add wannier centers visualization",
                "• Add effective masses calculation",
            ]
        },
        "5. Input Generation": {
            "Current": "User provides fixed strings for projections",
            "Improvements": [
                "• Add automatic projection inference from structure",
                "• Add per-element projection customization",
                "• Add validation of projection syntax",
                "• Add common projection templates",
                "• Add orbital analysis from pseudopotentials",
            ]
        },
    }
    
    for category, details in improvements.items():
        print(f"\n{category}")
        print(f"  Current: {details['Current']}")
        print(f"  Potential improvements:")
        for improvement in details['Improvements']:
            print(f"    {improvement}")
    
    # =========================================================================
    print(f"\n" + "="*80)
    print("GENERATED OUTPUT FILES (After successful run)")
    print("="*80)
    
    outputs = {
        "PROJWFC outputs": [
            "pwscf.pdos              - Total PDOS",
            "pwscf.pdos.up           - Spin-up PDOS",
            "pwscf.pdos.dw           - Spin-down PDOS",
            "pwscf.pdos.txt          - Human-readable format",
        ],
        "pw2wannier90 outputs": [
            "wannier_seed.amn        - Projection matrix",
            "wannier_seed.mmn        - Overlap matrix",
            "wannier_seed.eig        - Eigenvalues",
            "wannier_seed.unk        - Periodic part of wavefunctions",
        ],
        "wannier90 outputs": [
            "wannier_seed.wout       - Output log (contains spreads, centers)",
            "wannier_seed_centres.xyz - Wannier centers (XYZ format)",
            "wannier_seed_*.xsf      - Wannier function densities (for visualization)",
            "wannier_seed_band.dat   - Interpolated band structure (optional)",
        ],
    }
    
    for category, files in outputs.items():
        print(f"\n{category}:")
        for f in files:
            print(f"  • {f}")
    
    # =========================================================================
    print(f"\n" + "="*80)
    print("USAGE EXAMPLE")
    print("="*80)
    
    usage = """
from xespresso.workflow.wannier_workflow import WannierWorkflow

# Initialize workflow
wf = WannierWorkflow(
    cif_file='Si.cif',
    pseudopotentials={'Si': 'Si.pbe-n-kjpaw_psl.1.0.0.UPF'},
    protocol='moderate',
    num_wann=4,
    projections='Si: sp3',
    kpts_scf=(4, 4, 4),
    kpts_nscf=(6, 6, 6),
)

# Run all stages (SCF, Bands, PROJWFC, NSCF, pw2wannier, wannier90)
results = wf.run(
    blocking=True,
    run_bands_validation=True,
    run_projwfc_analysis=True,
    dry_run=True
)

# Access results
wf.get_projwfc_analysis()
wf.compare_bands()
wf.validate_wannier_quality()
"""
    print(usage)
    
    print("="*80 + "\n")


if __name__ == '__main__':
    main()
