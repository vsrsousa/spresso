# Convergence Workflow - Implementation Progress

**Last Updated**: 2026-03-03  
**Status**: ✅ Working - Ready for next phase  
**Current Focus**: Multi-property convergence criteria implementation

---

## 📋 Executive Summary

The `ConvergenceWorkflow` class implements a **two-phase independent algorithm** for DFT parameter convergence:
- **Phase 1**: ecutwfc convergence with fixed kspacing (0.3 Å⁻¹)
- **Phase 2**: kspacing convergence with optimal ecutwfc from Phase 1

Currently only `'energy'` criterion is implemented. This document tracks the implementation of additional criteria: `'forces'`, `'stress'`, `'geometry'`, `'magnetic_moments'`.

---

## ✅ COMPLETED WORK

### 1. Core Convergence Framework
- ✅ Two-phase independent algorithm (no nested loops)
- ✅ Reference value handling (ecutwfc=200, kspacing=0.1)
- ✅ Dynamic range expansion with intelligent stepping
- ✅ Result caching to avoid redundant calculations
- ✅ Convergence criteria parameter with default `['energy']`

### 2. Parameter Infrastructure  
- ✅ `convergence_criteria_list: List[str]` parameter added to `__init__()` and `optimize_parameters()`
- ✅ Default: `['energy']` (set in `_get_default_convergence_criteria_list()`)
- ✅ Tolerance handling via `convergence_criteria` dict with keys:
  - `'energy_tolerance'`: meV/atom (default: 3e-3 for 'low' precision)
  - `'force_tolerance'`: eV/Å (currently unused)
  - `'stress_tolerance'`: GPa (currently unused)
  - `'geometry_tolerance'`: Å (currently unused)
  - `'magnetic_tolerance'`: μB (currently unused)

### 3. Energy Convergence (✅ Complete)
**File**: `xespresso/workflow/convergence_workflow.py`

**Extraction** (lines 880-907):
```python
def _extract_property_from_result(completion, num_atoms, 'energy'):
    if 'energy' in completion:
        return completion['energy'] / num_atoms
```

**Convergence Check** (lines 1010-1044):
```python
def _check_convergence_vs_reference(...):
    if criterion == 'energy':
        energies = [props.get('energy', np.nan) for props in results_dict.values()]
        max_energy = max(e for e in energies if not np.isnan(e))
        energy_diff = abs(max_energy - reference_properties.get('energy', 0))
        tolerance = criteria_tolerances.get('energy_tolerance', 1e-3)
        if energy_diff >= tolerance:
            return False  # Not converged
```

---

## 🔄 IN PROGRESS / TODO

### Phase 1: Force Convergence

**Status**: ✅ **COMPLETE** (2026-03-04)  
**Commits**: Implementation + Tests  

**What was implemented**:

1. ✅ **Enable force calculation** (lines 880-920):
   - In `_extract_property_from_result()`, added 'forces' criterion
   - Extract max force magnitude from `completion['forces']` array (shape: N_atoms×3)
   - Compute `np.linalg.norm()` per atom, return max magnitude
   - Unit: eV/Å (QE native)

2. ✅ **QE configuration** (lines 827-857):
   - Updated `_get_calculation_config()` to accept 'forces' criterion
   - Added 'forces' to valid_criteria set
   - Automatically inserts `'tprnfor': True` when 'forces' in criteria_list

3. ✅ **Enable in convergence check** (lines 1010-1044):
   - Added elif for 'forces' criterion in `_check_convergence_vs_reference()`
   - Compares max force in test results vs reference
   - Uses `criteria_tolerances.get('force_tolerance', 0.05)` (eV/Å)

4. ✅ **Unit Tests** (test_convergence_multi_property.py):
   - `test_get_calculation_config_forces_implemented`: Config returns tprnfor=True
   - `test_extract_property_forces_implemented`: Extract max force correctly
   - `test_extract_property_forces_max_magnitude`: Magnitude calculation (3-4-5 triangle)
   - `test_check_convergence_forces_converged`: Force convergence check passes
   - `test_check_convergence_forces_not_converged`: Force convergence check fails properly
   - `test_check_convergence_energy_and_forces_both_required`: ALL criteria must converge
   - `test_check_convergence_energy_and_forces_both_converged`: Both pass → converged

**Test Results**: 22/22 tests PASSING ✅

---

### Phase 2: Stress Convergence

**Status**: ⏳ Ready to implement  
**Priority**: HIGH  
**Effort**: MEDIUM (identical to forces structure)

**What needs to happen**:

1. **Extract hydrostatic pressure** (lines 880-920):
   - In `_extract_property_from_result()`, handle `'stress'` criterion
   - Stress tensor from QE is 3×3 symmetric matrix (in kBar)
   - Hydrostatic pressure = -(trace(σ) / 3)
   - Convert kBar → GPa: `value_kbar * 0.1 = value_GPa`
   - Return: `float` with hydrostatic pressure in GPa

2. **Code location**: `xespresso/workflow/convergence_workflow.py:880-920`
   - Add elif block after forces
   - Parse stress tensor and compute pressure
   - Return `abs(pressure)` or signed depending on convention

3. **Enable in convergence check** (lines 1010-1044):
   - Add elif for `'stress'` criterion
   - Compare hydrostatic pressure in test results vs reference
   - Use `criteria_tolerances.get('stress_tolerance', 1.0)` (GPa)

4. **QE configuration**:
   - Add `tstress=True` when 'stress' criterion requested
   - Location: Update `_get_calculation_config()` (lines 827-857)
   - Logic: If 'stress' in criteria_list → add `'tstress': True` to overrides

**Technical Details**:
- Stress tensor shape: (3, 3) symmetric in kBar
- Definition: Hydrostatic = -trace(σ)/3
- Unit conversion: kBar to GPa (×0.1)

---

### Phase 3: Geometry Convergence

**Status**: ⏳ Complex - requires VC-RELAX  
**Priority**: MEDIUM  
**Effort**: HIGH

**What needs to happen**:

1. **Requires structural relaxation** (not just SCF):
   - Current: Only SCF calculations
   - Need: VC-RELAX (volume + cell relaxation) calculations
   - This affects `_get_calculation_config()` return value
   - Decision needed: Should geometry criteria trigger VC-RELAX instead of SCF?

2. **Extract atomic displacement**:
   - Compare final structure with initial
   - Could be: max displacement of any atom, or RMS displacement
   - Unit: Å (conventional)

3. **Implementation strategy**:
   - Option A: Separate calculation mode (vc-relax vs scf)
   - Option B: Always do VC-RELAX if geometry criterion requested
   - Recommend: **Option B** (simpler, more physically meaningful)

4. **Code changes**:
   - `_get_calculation_config()` (lines 827-857)
   - `_extract_property_from_result()` (lines 880-920)
   - Convergence check in `_check_convergence_vs_reference()` (lines 1010-1044)

---

### Phase 4: Magnetic Moments Convergence

**Status**: ⏳ Requires magnetism setup  
**Priority**: LOW (specialized)  
**Effort**: LOW-MEDIUM

**What needs to happen**:

1. **Requires magnetic calculation**:
   - Set `nspin=2` (or `nspin=4` for noncollinear)
   - Requires initial magnetic moments for structure
   - Location: `_get_calculation_config()` (lines 827-857)

2. **Extract total magnetic moment**:
   - From QE output: total magnetization
   - Unit: μB (Bohr magnetons)
   - Return single float value

3. **Convergence logic**:
   - Check total magnetic moment convergence vs reference
   - Tolerance: `criteria_tolerances.get('magnetic_tolerance', 0.01)` μB

---

## 🎯 NEXT IMMEDIATE STEPS

### Step 1: Implement Stress Extraction & Convergence ← **CURRENT**
**Time**: ~30 minutes  
**Files to edit**: 
- `xespresso/workflow/convergence_workflow.py`: lines 880-920 (extraction)
- `xespresso/workflow/convergence_workflow.py`: lines 827-857 (`_get_calculation_config`)
- `xespresso/workflow/convergence_workflow.py`: lines 1010-1044 (convergence check)

**What needs to happen** (same pattern as forces):
1. Extract hydrostatic pressure from stress tensor
2. Set `tstress=True` when 'stress' in criteria_list
3. Add convergence check for stress vs reference, use `stress_tolerance` (GPa)

### Step 3: Create Unit Tests
**Time**: ~1 hour  
**Location**: Create or update `tests/test_convergence_multi_property.py`

**Test scenarios**:
- Energy only (current baseline)
- Energy + Forces
- Energy + Stress
- Energy + Forces + Stress

---

## 📍 KEY CODE LOCATIONS

### Main Algorithm: `run_convergence_independent()`
- **PHASE 1** (ecutwfc): lines 1149-1315
- **PHASE 2** (kspacing): lines 1323-1538
- Both phases use same convergence logic via `_check_convergence_vs_reference()`

### Convergence Decision Points
1. **Line 1236**: Phase 1 convergence check
   ```python
   converged = self._check_convergence_vs_reference(
       test_ecut_results, reference_properties, criteria_tolerances,
       self.convergence_criteria_list
   )
   ```

2. **Lines 1450-1454**: Phase 2 convergence check (same pattern)
   ```python
   if reference_properties_phase2 is not None:
       test_ksp_results = {k: v for k, v in ksp_results.items() if k != max_kspacing}
       converged = self._check_convergence_vs_reference(
           test_ksp_results, reference_properties_phase2, criteria_tolerances,
           self.convergence_criteria_list
       )
   ```

### Property Extraction Point
- **Lines 1225-1235**: Phase 1 - extracting properties from completed calculation
  ```python
  props = {}
  for prop_name in self.convergence_criteria_list:
      props[prop_name] = self._extract_property_from_result(
          comp, len(self.atoms), prop_name
      )
  ```
  Same pattern in Phase 2 (lines 1413-1423)

### Convergence Tolerances
- **Lines 446-476**: `_get_default_convergence_criteria()` - defines tolerance values per precision level
- **Line 485**: Currently stored in `self.convergence_criteria` dict
- **Access**: `criteria_tolerances.get('energy_tolerance', 1e-3)` in convergence check

---

## 🔧 TECHNICAL DECISIONS MADE

### 1. Multiple Criteria Logic (ALL must converge)
- Current: `_check_convergence_vs_reference()` returns `True` **only if ALL criteria converge**
- Each criterion has independent tolerance
- If ANY criterion doesn't meet tolerance → keep expanding ranges

### 2. Reference Value Handling
- **ecutwfc reference**: 200 Ry (never tested against others, pure reference)
- **kspacing reference**: 0.1 Å⁻¹ (never tested against others, pure reference)
- Convergence check explicitly excludes reference via dictionary filtering:
  ```python
  test_ecut_results = {k: v for k, v in ecut_results.items() if k != max_ecutwfc}
  ```

### 3. Property Storage Format
- Results stored as dict: `parameter_value → {'energy': ..., 'forces': ..., ...}`
- Each property extracted on demand in `_extract_property_from_result()`
- Only properties in `self.convergence_criteria_list` are extracted

### 4. Kspacing Precision
- Using `.3f` format (3 decimal places) to prevent rounding collisions in cache keys
- Avoids issues where 0.272 and 0.303 both round to 0.30

---

## 🧪 TESTING STRATEGY

**Current Test File**: `test_scheduler_detection.py` (exists but limited)

**Recommended New Tests**:
1. Test with `['energy']` only (baseline)
2. Test with `['energy', 'forces']`
3. Test with mixed tolerances (tight energy, loose forces)
4. Test that reference value is properly excluded from convergence check
5. Test that All criteria must pass for convergence

---

## 📝 NOTES FOR CONTINUITY

### If Restarting:
1. Read this file first for context
2. Start with **Step 1: Force Extraction** (lowest hanging fruit)
3. All criteria follow same pattern - code structure is clear
4. Reference value filtering is critical (don't skip that part)
5. Remember: multiple criteria = ALL must converge (AND logic, not OR)

### Common Pitfalls to Avoid:
- ❌ Forgetting to exclude reference value in convergence check
- ❌ Not adding required QE input flags (tprnfor, tstress, nspin)
- ❌ Wrong units in tolerance comparisons
- ❌ Including reference in cache vs excluding during check (be consistent)

### Questions to Ask:
1. Should geometry convergence trigger VC-RELAX instead of SCF?
2. For magnetic moments - should we assume user provides initial moments?
3. Should stress use hydrostatic pressure or max eigenvalue of stress tensor?

---

## 🚀 READY TO START?

You're ready to implement forces convergence. Key files:
- `xespresso/workflow/convergence_workflow.py` (3 locations: extraction, config, check)
- No new files needed - all integrated into existing framework

Good luck! Reference this file anytime you lose context.
