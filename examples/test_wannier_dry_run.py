"""
Test dry_run mode with WannierWorkflow - generates all input files without execution.
"""

from ase.build import bulk
from xespresso.workflow.wannier_workflow import WannierWorkflow
from pathlib import Path

# Create test structure (Si diamond)
si = bulk('Si', 'diamond', a=5.431)

# Define pseudopotentials (using SSSP_1.3.0_PBE_efficiency)
pseudos = {
    'Si': '/pseudo/SSSP_1.3.0_PBE_efficiency/Si.pbe-spn-rrkjus_psl.1.0.0.UPF'
}

print("\n" + "="*70)
print("TEST: WannierWorkflow with dry_run=True")
print("="*70)
print("This will generate all input files without executing calculations\n")

# Initialize workflow
wf = WannierWorkflow(
    atoms=si,
    pseudopotentials=pseudopotentials,
    protocol='fast',  # Fast for testing
    num_wann=4,
    projections='Si: sp3',
    kpts_scf=(2, 2, 2),
    kpts_nscf=(3, 3, 3),
)

try:
    # Run with dry_run=True - generates inputs only
    print("\n[1/1] Calling wf.run(dry_run=True)")
    print("-" * 70)
    
    results = wf.run(
        blocking=True,
        seedname='si_test',
        run_bands_validation=True,  # Also test bands input generation
        run_projwfc_analysis=False,  # Skip PROJWFC for faster test
        dry_run=True  # INPUT FILES ONLY
    )
    
    print("\n" + "="*70)
    print("DRY RUN COMPLETED SUCCESSFULLY ✓")
    print("="*70)
    
    # List generated input files
    run_dir = results['run_dir']
    print(f"\nInput files generated in: {run_dir}/")
    
    if Path(f"{results['scf'].directory}").exists():
        scf_files = list(Path(results['scf'].directory).glob('*'))
        print(f"\nSCF inputs ({results['scf'].directory}/): {len(scf_files)} files")
        for f in scf_files[:5]:
            print(f"  - {f.name}")
    
    if 'bands' in results and results['bands']:
        bands_files = list(Path(results['bands'].directory).glob('*'))
        print(f"\nBands inputs ({results['bands'].directory}/): {len(bands_files)} files")
        for f in bands_files[:5]:
            print(f"  - {f.name}")
    
    if Path(f"{results['nscf'].directory}").exists():
        nscf_files = list(Path(results['nscf'].directory).glob('*'))
        print(f"\nNSCF inputs ({results['nscf'].directory}/): {len(nscf_files)} files")
        for f in nscf_files[:5]:
            print(f"  - {f.name}")
    
    print("\n" + "="*70)
    print("NEXT STEPS:")
    print("- Review generated input files in:", run_dir)
    print("- To execute full pipeline: wf.run(dry_run=False)")
    print("="*70 + "\n")

except Exception as e:
    print(f"\n✗ Test failed: {e}")
    import traceback
    traceback.print_exc()
