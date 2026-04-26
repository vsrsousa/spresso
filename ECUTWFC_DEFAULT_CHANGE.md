# ⚡ Default ECUTWFC Range Mode Changed to ARANGE (Step-Based)

**Date**: April 26, 2026  
**Breaking Change**: Yes (but easy to migrate)

---

## 🔄 What Changed

**OLD DEFAULT (Removed)**:
```python
wf = ConvergenceWorkflow(atoms, pseudopotentials)
# Behavior: Start with ecutwfc=30, then expand: 30 → 40 → 50 → ... (adaptive)
# Mode: Dynamic (MODE 4)
```

**NEW DEFAULT**:
```python
wf = ConvergenceWorkflow(atoms, pseudopotentials)
# Behavior: Test [30, 40, 50, 60, 70, ..., 200] systematically (predefined)
# Mode: ARANGE (MODE 3) with ecutwfc_step=10
```

---

## ✅ Why This Change?

1. **Predictable**: All values known upfront, easier to plan
2. **Reproducible**: Deterministic testing order (same every time)
3. **Batch-Friendly**: Better for HPC batch submission systems
4. **Resumable**: Easy to resume if job is interrupted
5. **Debugging**: Easier to track which ecutwfc values have been tested
6. **Documentation**: Explicit values in output logs

---

## 🔧 How to Migrate

### If you liked the OLD dynamic behavior:

You can simulate it by using an explicit small list:

```python
# Instead of implicit dynamic expansion, be explicit:
wf = ConvergenceWorkflow(
    atoms, pseudopotentials,
    ecutwfc_values=[30, 40, 50, 60, 70]  # Explicit (MODE 1)
)
```

Or use custom step size to control resolution:

```python
wf = ConvergenceWorkflow(
    atoms, pseudopotentials,
    ecutwfc_range={
        'min': 30,
        'max': 200,
        'step': 10  # Same as new default (MODE 3)
    }
)
```

### If you want a SMALLER range (faster testing):

```python
# Test fewer points with larger step
wf = ConvergenceWorkflow(
    atoms, pseudopotentials,
    ecutwfc_range={
        'min': 30,
        'max': 100,
        'step': 20  # Larger step = fewer points
    }
)
# Tests: [30, 50, 70, 90] → ~2 hours instead of ~9 hours
```

### If you want SPECIFIC literature values:

```python
# Test ONLY certain known-good values
wf = ConvergenceWorkflow(
    atoms, pseudopotentials,
    ecutwfc_values=[50, 60, 70, 80, 90]  # Explicit list (MODE 1)
)
```

---

## 📊 Impact on Compute Time

| Scenario | Old Behavior | New Behavior | Change |
|----------|---|---|---|
| Au bulk | 1-2 hrs (adaptive) | 9 hrs (full range) | +7-8 hrs |
| Cu bulk | 1-2 hrs (adaptive) | 9 hrs (full range) | +7-8 hrs |
| Fast stop | Not guaranteed | Tests all 18 values | Predictable |

**Tip**: Use custom step size to control time:
- `step=10` → 18 points → 9 hours (default)
- `step=20` → 9 points → 4.5 hours
- `step=30` → 6 points → 3 hours

---

## 🎯 New 3-Mode System

After removing the old dynamic mode, the new system has:

| Mode | How to Use | Best For |
|------|-----------|----------|
| **MODE 1** | `ecutwfc_values=[...]` | Specific literature values, comparisons |
| **MODE 2** | `ecutwfc_range={'min':30, 'max':70, 'n_points':5}` | Papers, uniform spacing, plots |
| **MODE 3** | Default or `ecutwfc_range={'min':30, 'max':200, 'step':10}` | Production runs, batch jobs |

---

## ✨ Benefits of New Default

✅ **Consistent Results**: Same values tested every time  
✅ **Better Scheduling**: All compute time known upfront  
✅ **Easier Resumption**: If job stops, easy to see what was done  
✅ **Better for Scripting**: No surprises with adaptive behavior  
✅ **Publication Ready**: Easy to document which ecutwfc values were tested  

---

## 🚀 Example: Complete Migration

Before (with old dynamic mode):
```python
# Old code
wf = ConvergenceWorkflow(
    atoms=bulk('Au'),
    pseudopotentials_config='SSSP_efficiency',
    min_ecutwfc=30,
    max_ecutwfc=200
)
wf.run_convergence()
```

After (new ARANGE default, same behavior):
```python
# New code - uses default ARANGE
wf = ConvergenceWorkflow(
    atoms=bulk('Au'),
    pseudopotentials_config='SSSP_efficiency'
    # min_ecutwfc=30, max_ecutwfc=200, ecutwfc_step=10 (all defaults)
)
wf.run_convergence()
# Now tests: [30, 40, 50, ..., 200] (systematic, not adaptive)
```

To get faster testing with new default:
```python
# Faster: only test key values
wf = ConvergenceWorkflow(
    atoms=bulk('Au'),
    pseudopotentials_config='SSSP_efficiency',
    ecutwfc_range={'min': 30, 'max': 100, 'step': 20}  # Fewer points
)
wf.run_convergence()
# Tests: [30, 50, 70, 90] → ~2 hours
```

---

## 📝 Documentation

For detailed information, see:
- [ECUTWFC_RANGE_SPECIFICATION.md](docs/ECUTWFC_RANGE_SPECIFICATION.md) - Full guide to 3 modes
- [examples/ecutwfc_range_specification.py](examples/ecutwfc_range_specification.py) - Code examples

---

## ⚠️ Summary

| Item | Before | After |
|------|--------|-------|
| Default mode | Dynamic (adaptive) | ARANGE (step-based) |
| Total modes | 4 | 3 |
| Default behavior | Stops early if converged | Tests all values |
| Compute time | Variable (1-3 hrs) | Fixed (9 hrs @ step=10) |
| Predictability | Low | High ✓ |
| Batch friendly | No | Yes ✓ |

**Action**: No code changes needed unless you relied on the old adaptive stopping behavior. If so, see migration examples above.
