# PROJWFC Stage in Wannier90 Workflow

## Overview

The **PROJWFC** (Projection on Atomic Wavefunctions) stage is now integrated into the complete Wannier90 workflow. This stage computes the projected density of states (PDOS) and projects the Kohn-Sham wavefunctions onto atomic basis functions, providing critical insights for selecting optimal Wannier function projections.

## What is PROJWFC?

PROJWFC is a Quantum Espresso tool that:
- Projects KS wavefunctions onto atomic wavefunctions (s, p, d, f)
- Computes partial density of states (PDOS) for each atom and orbital
- Identifies which orbitals contribute to the electronic structure at different energies
- Helps understand the character of bands and their potential grouping

## Pipeline Order

The complete workflow now executes in this order:

```
1. SCF Calculation
   ↓
2. Band Structure (optional, for validation)
   ↓
3. NSCF Calculation (with wavefunction collection)
   ↓
4. PROJWFC Analysis ← NEW STAGE
   ↓
5. pw2wannier90 (Wavefunction conversion)
   ↓
6. wannier90 (Wannier function generation)
```

## Why PROJWFC in the Wannier Workflow?

### 1. **Projection Guidance**
   - PDOS analysis shows which orbitals are "relevant" to your system
   - Helps you understand if initial projections are reasonable
   - Prevents selecting orbitals that don't contribute significantly

### 2. **Orbital Character Analysis**
   - Identify which atoms/orbitals form specific bands
   - Distinguish sp bands from d bands, for example
   - Useful for multi-element systems

### 3. **Quality Diagnostics**
   - If PDOS shows unexpected features, projections might need adjustment
   - Can reveal issues with pseudopotential choice or convergence

### 4. **Iterative Improvement**
   - First run: Use generic projections
   - Check PROJWFC results
   - Refine projections based on PDOS
   - Re-run wannier90 with improved projections

## Using PROJWFC in the Workflow

### Enable PROJWFC Analysis

```python
from xespresso.workflow.wannier_workflow import WannierWorkflow

workflow = WannierWorkflow(
    cif_file='si.cif',
    pseudos={'Si': '/path/to/Si.pbe.UPF'},
    protocol='moderate',
    num_wann=4,
    projections='Si: sp3'
)

results = workflow.run(
    blocking=True,
    run_projwfc_analysis=True,  # ← Enable PROJWFC
    run_bands_validation=True,
)
```

### Disable PROJWFC (Optional)

If you want to skip this stage for quick runs:

```python
results = workflow.run(
    blocking=True,
    run_projwfc_analysis=False,  # Skip PROJWFC
    run_bands_validation=True,
)
```

## PROJWFC Output Files

After PROJWFC completes, the following files are generated:

```
example_directory/
├── si.pdos              (Main PDOS file)
├── si.pdos.up           (Spin-up PDOS, if spinor calculation)
├── si.pdos.dw           (Spin-down PDOS)
├── si.pdos_atm#1(Si)_wfc#1(1s)     (Individual atom/orbital projections)
├── si.pdos_atm#1(Si)_wfc#2(2s)
├── si.pdos_atm#1(Si)_wfc#3(2p)
└── si.pdos_tot         (Total PDOS)
```

## Interpreting PROJWFC Results

### PDOS File Format

Each line in the PDOS file contains:
```
Energy  Total_DOS  s_projection  p_projection  d_projection  f_projection
```

### How to Analyze

1. **Plot PDOS**
   ```bash
   gnuplot -> plot "si.pdos" u 1:2 w l
   ```

2. **Identify Important Bands**
   - Look for sharp peaks (high orbital contributions)
   - Note the energy ranges of different orbitals
   - Check if projections align with your initial guess

3. **Evaluate Coverage**
   - Do selected orbitals cover all important features?
   - Are there unexpected contributions?

### Example Interpretation

For Si (4 valence electrons, sp3 hybridization):

```
s orbital:  Strong contribution ~5-8 eV below Fermi level
p orbital:  Maximum near Fermi level (valence bands)
d orbital:  Very weak (should be unoccupied)
```

**Interpretation**: Projections on s and p orbitals are appropriate.

## Workflow Integration Details

### Automatic PDOS Analysis

The workflow automatically:
1. Writes `projwfc.in` with proper parameters
2. Executes `projwfc.x`
3. Collects PDOS output files
4. Parses results for suggestions (if implemented)

### Access PROJWFC Results Programmatically

```python
# Get PROJWFC analysis results
projwfc_analysis = workflow.get_projwfc_analysis()

if projwfc_analysis:
    print(f"PDOS file: {projwfc_analysis['pdos_file']}")
    print(f"Status: {projwfc_analysis['status']}")
    print(f"Results directory: {projwfc_analysis['run_dir']}")
    
    # Check if suggestions were parsed
    if 'parsed' in projwfc_analysis:
        print(f"Suggestions: {projwfc_analysis['parsed']}")
```

## Best Practices

### 1. Initial Projections
Start with reasonable guesses based on chemistry:
- Transition metals: typically `d s p`
- Main group elements: typically `s p`
- f-block elements: `f d s p`

### 2. Review PDOS
Always review PDOS output to verify projections:
```bash
cd run_directory
plot_pdos.py  # Use your favorite plotting tool
```

### 3. Iterative Refinement
If Wannier functions don't converge well:
1. Check PDOS for unexpected features
2. Adjust projections
3. Re-run wannier90

### 4. Multiple Iterations
For difficult systems, multiple iterations may be needed:

```python
# Run 1: Initial guess
workflow = WannierWorkflow(..., projections='auto')
results = workflow.run(run_projwfc_analysis=True)

# Check PDOS results...

# Run 2: Refined projections
workflow = WannierWorkflow(..., projections='Si: sp3d2')  # Updated
results = workflow.run(run_projwfc_analysis=False)  # Skip PROJWFC
```

## Testing PROJWFC Stage

The `run_projwfc()` function can be used independently:

```python
from xespresso.workflow.wannier_workflow import run_projwfc

result = run_projwfc(
    run_dir='/path/to/nscf/calculation',
    prefix='si',
    blocking=True,
    queue=None,
)

print(f"Status: {result['status']}")
print(f"PDOS file: {result['outputs']['pdos']}")
```

## Troubleshooting

### PROJWFC fails or doesn't generate output

1. **Check NSCF completion**
   - Verify NSCF ran successfully
   - Check for `.save` directory

2. **Check disk space**
   - PDOS files can be large for dense k-point grids
   - Ensure sufficient space in run directory

3. **Review projwfc.in**
   - Check if file was written correctly
   - Verify prefix matches NSCF calculation

### PDOS shows unexpected features

1. **Check pseudopotentials**
   - Verify you're using appropriate pseudopotentials
   - Check for semicore effects

2. **Verify convergence**
   - Check energy cutoff adequacy
   - Verify k-point density

3. **Review projections**
   - May need to include additional orbitals
   - Consider energy window for disentanglement

## Advanced Usage

### Custom PROJWFC Input

For special cases, you can provide custom input:

```python
from xespresso.workflow.wannier_workflow import run_projwfc

custom_input = """&inputpp
  prefix = 'si'
  outdir = './'
/
filpdos = 'si.pdos'
ngauss = 0
degauss = 0.01
"""

result = run_projwfc(
    run_dir='/path/to/nscf',
    prefix='si',
    # Use custom input somehow (needs additional parameter)
)
```

### Parsing PROJWFC Results

Custom parsing for automated analysis:

```python
from xespresso.workflow.wannier_workflow import parse_projwfc_output

suggestions = parse_projwfc_output(
    projwfc_dir='/path/to/results',
    prefix='si'
)

if suggestions:
    print(f"Analysis: {suggestions}")
```

## Integration with Other Tools

### Export to Other Codes

PDOS can be exported to external tools:
- **VASP**: Similar PDOS format
- **AbiPy**: Read and analyze Quantum Espresso PDOS
- **ASE**: Interface with projwfc results
- **PyMatGen**: Materials analysis

### Combine with Band Structure

For complete validation:

```python
# Get band structure comparison
bands_comparison = workflow.compare_bands()

# Get PROJWFC analysis
projwfc_analysis = workflow.get_projwfc_analysis()

# Together, they provide comprehensive understanding of:
# - Which bands are important (PDOS)
# - How well Wannier functions reproduce them (band comparison)
```

## References

1. **Quantum Espresso Documentation**
   - PROJWFC tool: https://www.quantum-espresso.org/
   - Input parameters, output formats

2. **Wannier90 Guide**
   - Initial projections: http://www.wannier.org/
   - Choosing projections chapter

3. **Related Tools**
   - plotband.x: Plot band structures
   - dos.x: Additional DOS calculations

4. **Theory**
   - Marzari & Vanderbilt (1997): Maximally localized generalized Wannier functions
   - Check specific references in wannier.org

## Summary

The PROJWFC stage provides essential insight into your electronic structure, enabling:
- **Better projections** → More accurate Wannier functions
- **Quality diagnostics** → Identify potential issues early
- **Iterative improvement** → Systematically optimize results
- **Scientific understanding** → See what your Wannier functions represent

By including this analysis step in the workflow, you ensure that your Wannier functions are not just numerically converged, but physically meaningful and well-suited to your specific material.
