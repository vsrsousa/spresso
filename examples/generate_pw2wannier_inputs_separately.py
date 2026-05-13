"""
Example: Generate pw2wannier90 and wannier90 input files independently
from the WannierWorkflow orchestration.

This is useful when you want to:
1. Run SCF and NSCF separately with CalculationWorkflow
2. Generate input files for pw2wannier90 and wannier90 at your own pace
3. Submit them independently or use them in custom scripts
4. Inspect inputs before running

The workflow approach is still available via WannierWorkflow class,
but this shows the standalone helpers for fine-grained control.
"""
from pathlib import Path

from xespresso.workflow.wannier_workflow import (
    generate_pw2wannier_input,
    generate_seedname_win,
)


def main():
    # After running SCF and NSCF separately (see wannier_pipeline_advanced_example.py)
    # You can generate the input files independently
    
    # Paths and parameters
    nscf_run_dir = Path("runs/04-nscf").resolve()  # Where NSCF outputs are
    wan_run_dir = Path("runs/05-wan").resolve()    # Where pw2wannier will run
    system_prefix = "Fe_bcc"                       # System name (from CIF)
    seedname = "wannier_seed"                      # Wannier90 seed name
    
    # Create wan_run_dir if needed
    wan_run_dir.mkdir(parents=True, exist_ok=True)
    
    # ===== 1. Generate pw2wannier90 input =====
    print("[1] Generating pw2wannier90 input file...")
    
    pw2wan_content = generate_pw2wannier_input(
        prefix=system_prefix,              # Read from Fe_bcc.save
        seedname=seedname,                 # Output files: wannier_seed.amn, etc.
        nscf_run_dir=str(nscf_run_dir),   # Where to find Fe_bcc.save
        wan_run_dir=str(wan_run_dir)      # Where pw2wannier will run (calculates outdir)
    )
    
    # Write to file
    pw2wan_path = wan_run_dir / "pw2wannier.in"
    with open(pw2wan_path, "w") as f:
        f.write(pw2wan_content)
    
    print(f"✓ Generated: {pw2wan_path}")
    print("\nContent:")
    print(pw2wan_content)
    
    # ===== 2. Generate wannier90 input (.win file) =====
    print("\n[2] Generating wannier90 input file...")
    
    # System-specific parameters
    num_wann = 16           # For Fe (5d + 3s + 3p + 4s + 4p = 16 bands)
    projections = """Fe: d"""  # Fe d-projections
    spinors = False         # Non-magnetic system (or use True for spinor calc)
    dis_num_iter = 1000     # Disentanglement iterations
    
    win_content = generate_seedname_win(
        num_wann=num_wann,
        projections=projections,
        spinors=spinors,
        dis_num_iter=dis_num_iter
    )
    
    # Write to file
    win_path = wan_run_dir / f"{seedname}.win"
    with open(win_path, "w") as f:
        f.write(win_content)
    
    print(f"✓ Generated: {win_path}")
    print("\nFirst 30 lines:")
    print("\n".join(win_content.split("\n")[:30]))
    
    # ===== 3. Summary =====
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Input files ready in: {wan_run_dir}")
    print(f"  - pw2wannier.in    (generated from NSCF at {nscf_run_dir})")
    print(f"  - {seedname}.win   (Wannier90 input)")
    print("\nNext steps:")
    print(f"1. Run pw2wannier90 from {wan_run_dir}:")
    print(f"   cd {wan_run_dir}")
    print(f"   pw2wannier90.x -in pw2wannier.in")
    print(f"\n2. Run wannier90:")
    print(f"   cd {wan_run_dir}")
    print(f"   wannier90.x {seedname}")


if __name__ == "__main__":
    main()
