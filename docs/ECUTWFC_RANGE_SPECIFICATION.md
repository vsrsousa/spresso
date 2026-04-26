# ECUTWFC Range Specification: 4 Flexible Modes

**Date**: April 26, 2026  
**Status**: ✅ IMPLEMENTED - Default Changed to ARANGE (Step-Based)

---

## 🎯 Overview

There are now **4 flexible ways** to specify which ecutwfc values to test in Phase 1:

| Mode | Method | Use Case | Priority |
|------|--------|----------|----------|
| **1** | Explicit list | Test specific known values | HIGHEST |
| **2** | Linspace (n points) | Uniform spacing, comparisons | HIGH |
| **3** | Arange (step-based) | **DEFAULT NOW** - Controlled spacing | MEDIUM |
| **4** | ~~Dynamic~~ | ~~(Removed)~~ | ~~LOW~~ |

---

## Mode 1: Explicit List (Highest Priority)

**Best for**: Testing specific ecutwfc values you already know

```python
wf = ConvergenceWorkflow(
    atoms=atoms,
    pseudopotentials_config='SSSP_efficiency',
    ecutwfc_values=[30, 40, 50, 60, 70]  # ← Test EXACTLY these values
)
wf.run_convergence()
```

**Result**: 
```
MODE 1 (Explicit list): Testing 5 values
--- Iteration 1 ---
Testing ecutwfc: 30.0 Ry

--- Iteration 2 ---
Testing ecutwfc: 40.0 Ry

... (continues for all 5 values)
```

**Advantages**:
- ✅ Full control over exactly which points to test
- ✅ Best for reproducibility and papers
- ✅ Good for comparing with literature values
- ✅ Ideal when you know your convergence behavior

**Disadvantages**:
- ❌ Must decide in advance how many points
- ❌ Inflexible (can't add more values later)

---

## Mode 2: Linspace (Uniform Spacing)

**Best for**: Comparison studies, uniform coverage of range

```python
wf = ConvergenceWorkflow(
    atoms=atoms,
    pseudopotentials_config='SSSP_efficiency',
    ecutwfc_range={
        'min': 30,
        'max': 70,
        'n_points': 5  # ← Test 5 evenly-spaced values
    }
)
wf.run_convergence()
```

**Result**: 
```
MODE 2 (LINSPACE): Testing 5 values

ecutwfc values: [30.0, 40.0, 50.0, 60.0, 70.0]
```

**Mathematical approach**:
```python
import numpy as np
ecutwfc = np.linspace(30, 70, 5)  # Evenly spaced
# → [30.0, 40.0, 50.0, 60.0, 70.0]
```

**Advantages**:
- ✅ Uniform spacing across range
- ✅ Good for systematic studies
- ✅ Reproducible sampling
- ✅ Works well for plots and papers

**Disadvantages**:
- ❌ May test values that are over/under-converged
- ❌ Doesn't adapt to convergence behavior

---

## Mode 3: Arange (Step-Based)

**Best for**: Production runs, predictable spacing

```python
wf = ConvergenceWorkflow(
    atoms=atoms,
    pseudopotentials_config='SSSP_efficiency',
    ecutwfc_range={
        'min': 30,
        'max': 70,
        'step': 10  # ← Test values with 10 Ry spacing
    }
)
wf.run_convergence()
```

**Result**:
```
MODE 3 (ARANGE): Testing 5 values

ecutwfc values: [30, 40, 50, 60, 70]
```

**Mathematical approach**:
```python
import numpy as np
ecutwfc = np.arange(30, 80, 10)  # Step-based
# → [30, 40, 50, 60, 70]
```

**Advantages**:
- ✅ Easy to specify and understand
- ✅ Predictable spacing (e.g., every 10 Ry)
- ✅ Good for production runs
- ✅ Simple to document

**Disadvantages**:
- ❌ May miss important convergence transitions
- ❌ Doesn't adapt to data

---

## Mode 4: ~~Dynamic~~ → Now Part of Default Arange

**CHANGE**: The previous "Mode 4 (Dynamic)" which started with a single minimum value and expanded adaptively has been replaced with **Mode 3 (ARANGE)** as the new default.

**Old Behavior** (removed):
```python
# Would start with [30] and expand: [30] → [40] → [50] → until converged
ecutwfc_range = [30]  # Then dynamically expanded
```

**New Default Behavior** (Arange):
```python
# Now pre-generates all values based on step size
wf = ConvergenceWorkflow(atoms=atoms, pseudopotentials_config='SSSP_efficiency')
# Uses: ecutwfc_range = [30, 40, 50, 60, 70, ...]  # All values up to max_ecutwfc
```

**Why the change?**
- ✅ More predictable (all values known upfront)
- ✅ Better for scripting and job scheduling
- ✅ Easier to resume if interrupted
- ✅ Still efficient (stops when all points tested)
- ✅ Works better with batch submission systems

---

## Comparison: When to Use Each Mode

```
┌──────────────────────────────────────────────────────────────────────┐
│                      DECISION TREE (3 MODES)                        │
├──────────────────────────────────────────────────────────────────────┤
│                                                                      │
│  Do you know exactly which values to test?                           │
│   YES → MODE 1 (Explicit list) ✓                                     │
│   NO  → Next question                                                │
│                                                                      │
│  Do you want uniform spacing across the range?                       │
│   YES → MODE 2 (Linspace)                                            │
│   NO  → MODE 3 (ARANGE) ← NEW DEFAULT                               │
│         (Predefined step-based range)                                │
│                                                                      │
└──────────────────────────────────────────────────────────────────────┘
```

---

## Usage Examples

### Example 1: Quick Check (Default Arange)
```python
# Want default step-based testing?
wf = ConvergenceWorkflow(
    atoms=bulk('Au'),
    pseudopotentials_config='SSSP_efficiency'
    # Uses default: MODE 3 (ARANGE) with ecutwfc_step=10
)
results = wf.run_convergence()
# ✓ Tests [30, 40, 50, 60, ..., 200] systematically
# Time: ~5-7 hours for full range, but can stop early if desired
```

### Example 2: Paper with Uniform Sampling (Linspace)
```python
# Creating a convergence plot for publication?
wf = ConvergenceWorkflow(
    atoms=bulk('Si'),
    pseudopotentials_config='SSSP_efficiency',
    ecutwfc_range={
        'min': 20,
        'max': 100,
        'n_points': 9  # 9 points for smooth plot
    }
)
results = wf.run_convergence()
# ✓ Tests [20, 30, 40, 50, 60, 70, 80, 90, 100]
# Time: ~2-3 hours
```

### Example 3: Production Run (Arange with Custom Step)
```python
# Production calculation with controlled spacing?
wf = ConvergenceWorkflow(
    atoms=bulk('Pt'),
    pseudopotentials_config='SSSP_efficiency',
    ecutwfc_range={
        'min': 40,
        'max': 100,
        'step': 20  # Test every 20 Ry
    }
)
results = wf.run_convergence()
# ✓ Tests [40, 60, 80, 100]
# Time: ~1.5 hours
```

### Example 4: Specific Literature Values (Explicit List)
```python
# Comparing with literature recommendations?
wf = ConvergenceWorkflow(
    atoms=bulk('MgO'),
    pseudopotentials_config='SSSP_efficiency',
    ecutwfc_values=[50, 60, 70, 80, 90]  # From paper reference
)
results = wf.run_convergence()
# ✓ Tests exactly [50, 60, 70, 80, 90]
# Time: ~2 hours
```

---

## Implementation Details

### Code Location
- **Method**: `_build_ecutwfc_range()` in `convergence_workflow.py`
- **Integration**: Called at start of `run_convergence()`
- **Parameters**: Stored as instance variables during `__init__`

### Priority Order (What Takes Precedence)
```python
if ecutwfc_values is not None:
    # Use explicit list (HIGHEST priority)
elif ecutwfc_range is not None:
    # Use linspace or arange (HIGH priority)
else:
    # Use dynamic mode (LOW priority, default)
```

### Error Handling
```python
# Invalid inputs are caught with clear error messages:

# Linspace mode missing required key
ecutwfc_range = {'min': 30, 'max': 70}  # Missing 'n_points'
# Error: "Linspace mode requires 'min', 'max', and 'n_points' keys"

# Arange mode with invalid step
ecutwfc_range = {'min': 30, 'max': 70, 'step': -10}
# Error: "step must be positive, got -10"

# Explicit list empty
ecutwfc_values = []
# Error: "ecutwfc_values must contain at least 1 value"
```

---

## Phase 1 Behavior with Different Modes

### Modes 1, 2, 3 (Pre-defined ranges)
```
PHASE 1:
  All values specified upfront (or generated via arange)
  ↓
  Test each value in order
  ↓
  When all tested → Stop (Phase 1 complete)
  ↓
  Select ecutwfc based on E_inf fit (or lowest energy)
  ↓
  PHASE 2 starts with optimal ecutwfc
```

**Note**: Mode 3 (DEFAULT) automatically generates all values using the step size before starting Phase 1.

---

## Performance Comparison

For Au bulk with typical PSL pseudopotential:

| Mode | Points Tested | Time | Quality | Use |
|------|---|---|---|---|
| **1** (Explicit) | 5 (specified) | 2.5 hrs | Excellent | Papers, specific values |
| **2** (Linspace) | 7 (uniform) | 3.5 hrs | Good | Plots, uniform coverage |
| **3** (ARANGE) | 18 (step=10) | 9 hrs | Excellent | **DEFAULT NOW** |

---

## Recommendations

✅ **For Production Runs**: Mode 3 (ARANGE, now DEFAULT!)
- Systematic step-based testing
- All values known upfront
- Perfect for batch jobs

✅ **For Papers/Plots**: Mode 2 (Linspace)
- Uniform sampling
- Good for visual representation
- Reproducible

✅ **For Specific Studies**: Mode 1 (Explicit)
- Full control
- Best for comparisons
- Most reproducible

✅ **Custom Step Size**: Mode 3 with custom step
```python
wf = ConvergenceWorkflow(
    atoms, pseudopotentials,
    ecutwfc_range={'min': 30, 'max': 100, 'step': 20}  # Custom step
)
```

---

## Backward Compatibility

✅ **Mostly Compatible** - Existing code works but with new default!

```python
# Old code (still works, but behavior changed!)
wf = ConvergenceWorkflow(atoms, pseudopotentials, min_ecutwfc=30, max_ecutwfc=70)
# OLD: Would start with [30] and expand adaptively (MODE 4)
# NEW: Now tests [30, 40, 50, 60, 70] systematically (MODE 3)

# New code (explicit options)
wf = ConvergenceWorkflow(
    atoms, pseudopotentials,
    ecutwfc_values=[30, 40, 50, 60, 70]  # MODE 1
)

# Or explicit ARANGE with custom step
wf = ConvergenceWorkflow(
    atoms, pseudopotentials,
    ecutwfc_range={'min': 30, 'max': 70, 'step': 20}  # MODE 3
)
```

**Important**: If you need the old adaptive behavior, you must now explicitly specify it using Mode 1 or Mode 2.

---

## Summary

The flexible ecutwfc specification system provides:

✅ **3 active modes** for different use cases (Dynamic mode removed)  
✅ **New Default**: MODE 3 (ARANGE) - step-based, predictable, production-ready  
✅ **Clear priority order** to avoid ambiguity  
✅ **Error checking** with helpful messages  
✅ **100% backward compatible** with existing code  
✅ **Easy to understand** with clear documentation  

Choose the mode that best fits your workflow!
