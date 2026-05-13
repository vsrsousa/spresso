# Simplified Convergence API: phases parameter only

**Date**: April 26, 2026  
**Status**: ✅ IMPLEMENTED

---

## ✨ Simplified Approach

Instead of creating 3 separate methods, we now use a **single method with a parameter**:

```python
# Just use run_convergence() with phases parameter!
wf.run_convergence(phases='ecut')   # Only PHASE 1
wf.run_convergence(phases='kpt')    # Only PHASE 2
wf.run_convergence(phases='both')   # Both (default)
wf.run_convergence()                # Default = both
```

---

## 🗑️ What Was Removed

We removed the extra wrapper methods that were redundant:

- ❌ `run_ecut_convergence()` - Now use `run_convergence(phases='ecut')`
- ❌ `run_kpt_convergence()` - Now use `run_convergence(phases='kpt')`

These were just wrappers calling `run_convergence()` anyway, so having them separate made the API unnecessarily complex.

---

## 🎯 Cleaner API

Instead of:
```python
# 3 different method names
wf.run_ecut_convergence()      # What's a convergence?
wf.run_kpt_convergence()       # Different method
wf.run_convergence()  # And another
```

We now have:
```python
# 1 method with clear parameter
wf.run_convergence(phases='ecut')
wf.run_convergence(phases='kpt')
wf.run_convergence(phases='both')
```

**Benefits**:
- ✅ Simpler to understand
- ✅ Consistent method name
- ✅ Parameter makes intent explicit
- ✅ Easier to document
- ✅ Less code to maintain

---

## 📚 How to Use

### Run only ecutwfc (PHASE 1)
```python
wf = ConvergenceWorkflow(atoms, pseudopotentials)
results = wf.run_convergence(phases='ecut')
print(f"Optimal ecutwfc: {wf.optimal_ecutwfc} Ry")
```

### Run only kspacing (PHASE 2)
```python
wf = ConvergenceWorkflow(atoms, pseudopotentials)
wf.optimal_ecutwfc = 50.0  # Pre-computed ecutwfc
results = wf.run_convergence(phases='kpt')
print(f"Optimal kspacing: {wf.optimal_kspacing} Å⁻¹")
```

### Run both (default)
```python
wf = ConvergenceWorkflow(atoms, pseudopotentials)
results = wf.run_convergence()  # Default: phases='both'
# or explicit:
results = wf.run_convergence(phases='both')
```

---

## 🔄 Migration from old API (if you used private methods)

**Before** (hypothetical):
```python
results = wf.run_ecut_convergence()
results = wf.run_kpt_convergence()
```

**Now**:
```python
results = wf.run_convergence(phases='ecut')
results = wf.run_convergence(phases='kpt')
```

---

## ✅ Implementation Details

**In convergence_workflow.py**:
- ✅ Added `phases` parameter to `run_convergence()`
- ✅ Added validation for `phases` parameter (only 'ecut', 'kpt', 'both')
- ✅ Added early return if `phases='ecut'` (skip PHASE 2)
- ✅ Added check for pre-computed `optimal_ecutwfc` if `phases='kpt'`
- ✅ Updated `run_convergence_study()` to pass through `phases`

**Files Updated**:
- ✅ `xespresso/workflow/convergence_workflow.py` - Core implementation
- ✅ `docs/CONVERGENCE_PHASES.md` - Complete documentation
- ✅ `examples/convergence_phases.py` - 5 working examples

---

## 📖 Documentation

See [CONVERGENCE_PHASES.md](docs/CONVERGENCE_PHASES.md) for complete guide with:
- API reference for `run_convergence(phases=...)`
- Usage examples
- Error handling
- Workflow patterns

---

## ✨ Summary

**One method, three modes**:

```python
wf.run_convergence(phases='ecut')    # PHASE 1 only
wf.run_convergence(phases='kpt')     # PHASE 2 only
wf.run_convergence(phases='both')    # Both (default)
```

Clean, simple, effective! 🎉
