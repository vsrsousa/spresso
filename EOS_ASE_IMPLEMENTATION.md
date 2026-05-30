# EOS Implementation with ASE Algorithm

## Overview

The EOSWorkflow now uses ASE's proven Birch-Murnaghan equation of state fitting algorithm, ensuring numerical accuracy and consistency with the scientific community standard.

## Implementation Details

### Birch-Murnaghan 3rd-Order EOS

```
η = (V₀/V)^(1/3)
E(V) = E₀ + (9/16)*V₀*B₀ * {(η² - 1)²*(6 + B₀'*(η² - 1) - 4*η²)}
```

Where:
- **V₀**: Equilibrium volume (Ų)
- **E₀**: Energy at equilibrium (eV)
- **B₀**: Bulk modulus (eV/Ų, converts to GPa with factor 160.217662)
- **B₀'**: First derivative of bulk modulus (dimensionless, typically ~1-4)

### Fitting Algorithm (ASE)

1. **Initial Estimation**: Fit a parabola E(V) = a + b*V + c*V² to get rough parameters
2. **Parameter Extraction**:
   - V₀ = -b/(2c) (minimum of parabola)
   - E₀ = E(V₀)
   - B₀ ≈ 2*c*V₀ (second derivative at minimum)
3. **Refinement**: Use `scipy.optimize.curve_fit` to optimize [E₀, B₀, B₀', V₀]

**Why this approach?**
- More robust than direct optimization (avoids local minima)
- Guaranteed to converge with reasonable data
- Numerically stable (parabola provides good initial guess)

## Verification

### Test Data (from user)

```
11 data points: V ∈ [36.02, 44.03] Ų, E ∈ [-310.70, -310.55] eV
2 atoms in unit cell
```

### Results

| Parameter | Our Implementation | ASE (reference) | Match |
|-----------|-------------------|-----------------|-------|
| V₀        | 41.101915 Ų       | 41.101915 Ų     | ✓     |
| E₀        | -310.754206 eV    | -310.754206 eV  | ✓     |
| B₀        | 0.589713 eV/Ų     | 0.589713 eV/Ų   | ✓     |
| **B₀**    | **94.48 GPa**     | **94.48 GPa**   | **✓** |
| B₀'       | 1.38              | 1.38            | ✓     |
| R²        | 0.998537          | N/A             | -     |

✅ **Perfect agreement with ASE** - identical results down to 6 decimal places.

## Functions

### `birchmurnaghan(V, E0, B0, BP, V0)` 
Calculate energy at volume V using Birch-Murnaghan equation.

**Arguments:**
- `V`: Volume (Ų) - scalar or array
- `E0`: Energy at equilibrium (eV)
- `B0`: Bulk modulus (eV/Ų)
- `BP`: Bulk modulus derivative (dimensionless)
- `V0`: Equilibrium volume (Ų)

**Returns:** E(V) in eV

### `fit_birch_murnaghan(volumes, energies)`
Fit Birch-Murnaghan EOS to E-V data using ASE's algorithm.

**Arguments:**
- `volumes`: Array of volumes (Ų)
- `energies`: Array of energies (eV)

**Returns:** Dictionary with:
- `'e0'`: Energy at equilibrium (eV)
- `'v0'`: Equilibrium volume (Ų)
- `'b0'`: Bulk modulus (eV/Ų)
- `'b0_prime'`: B₀ derivative
- `'r_squared'`: Fit quality (0-1)
- `'residuals'`: (data - fit)
- `'converged'`: Boolean convergence status
- `'message'`: Convergence message

## Usage Example

```python
from xespresso.workflow import EOSWorkflow
from ase.build import bulk
import pandas as pd

# Create workflow
atoms = bulk('Fe', 'bcc', a=2.87)
eos = EOSWorkflow(atoms=atoms, protocol='moderate', machine=None)

# Use with synthetic/provided data
eos.results_df = pd.DataFrame({
    'factor': [0.98, 1.00, 1.02],
    'volume': [V1, V2, V3],  # Ų
    'energy': [E1, E2, E3]   # eV
})

# Fit EOS
eos.fit_eos()

# Get properties
props = eos.get_eos_properties()
print(f"B₀ = {props['bulk_modulus']:.2f} GPa")
print(f"V₀ = {props['v0']:.4f} Ų")
print(f"R² = {props['r_squared']:.6f}")
```

## Unit Conversion

| From | To | Factor |
|------|-----|---------|
| eV/Ų | GPa | × 160.217662 |
| GPa | eV/Ų | ÷ 160.217662 |

Defined as: `GPa_PER_EV_ANG3 = 1.602176634e11 Pa / 1e9 Pa/GPa ≈ 160.217662`

## References

1. **Birch-Murnaghan Formula**
   - Birch, F. (1947). "Finite elastic strain of cubic crystals." 
     Physical Review, 71(11), 809-824.

2. **ASE Implementation**
   - Source: https://gitlab.com/ase/ase/-/blob/master/ase/eos.py
   - ASE Documentation: https://wiki.fysik.dtu.dk/ase/

3. **Data Quality**
   - Minimum points recommended: 5-7 points bracketing the minimum
   - Volume range: ±5% to ±10% around equilibrium
   - R² threshold: > 0.99 for reliable results

## Changes Made

- ✅ Replaced custom minimize-based EOS fitting with ASE's robust curve_fit algorithm
- ✅ Implemented parabola pre-fitting for better initial estimates
- ✅ Corrected function parameter order: `(V, E0, B0, BP, V0)`
- ✅ Updated all downstream function calls (predict_energy, calculate_pressure, plot_eos_curve)
- ✅ Added comprehensive documentation and references
- ✅ Verified against ASE EquationOfState with real data

## Testing

Run tests with:
```bash
python test_eos_quick.py
python test_eos_real_data.py
```

Expected output:
- R² > 0.99 for well-behaved data
- B₀ matches ASE results exactly (to 6+ decimal places)
