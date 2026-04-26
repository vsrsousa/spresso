# Outlier Detection in Exponential Fit: Practical Guide

## The Problem

When running Phase 1 ecutwfc convergence, first few points can have **massive energy deviations**:

```
Au bulk example:
  ecutwfc = 30 Ry  →  ΔE = 2000 meV  ❌ Huge deviation!
  ecutwfc = 40 Ry  →  ΔE =  100 meV  ⚠ Significant
  ecutwfc = 50 Ry  →  ΔE =   10 meV  ✓ Normal
  ecutwfc = 60 Ry  →  ΔE =    2 meV  ✓ Converged
```

The pseudopotential basis set is incomplete at low ecutwfc, causing non-exponential behavior. Including these in the fit **inflates the recommended ecutwfc** from 50 → 70 Ry.

---

## Two Detection Methods

### Method 1: Statistical (Z-Score) - DEFAULT ✅

**How it works:**
- Fit exponential to ALL data
- Calculate residuals (deviation from fit)
- Find points where |residual| > threshold × σ (standard deviation)
- Exclude and refit

**Threshold options:**
```
2.0σ (default)  → Moderate filtering, keeps most data
1.5σ            → More aggressive
2.5σ            → More conservative
```

**Pros:**
- Automatic, no manual choice needed
- Adapts to your data (different datasets have different noise levels)
- Statistically principled

**Cons:**
- Less intuitive (what is 2.0σ in meV?)
- Can miss patterns if noise is uniform

**Best for:** 
- Automated pipelines
- Don't know pseudopotential behavior

---

### Method 2: Absolute meV Threshold - PRACTICAL

**How it works:**
- Mark any point with |ΔE| > threshold as outlier
- Exclude and refit

**Threshold options for Au-like systems:**

| Threshold | Behavior | Example |
|-----------|----------|---------|
| **5 meV** | ❌ Too aggressive | Removes converged data |
| **10 meV** | ⚠ Aggressive | Removes some valid early points |
| **20 meV** | ✓ Reasonable | Good balance |
| **50 meV** | ✓ Conservative | Keeps more early data |
| **100 meV** | ✓✓ Safe | Only removes truly anomalous |

**Pros:**
- Intuitive (you control exactly what counts as "bad")
- Works well if you know your pseudopotential

**Cons:**
- Requires manual choice
- Different for each pseudopotential/system
- May keep too much or too little data

**Best for:**
- You know pseudopotential behavior
- Want explicit control

---

## Recommendation Strategy

### For Au (and most PBE pseudopotentials)

**Step 1: Inspect data first**
```python
# Run Phase 1
wf = ConvergenceWorkflow(atoms, pseudopotentials, precision='low')
wf.run_convergence_study(max_ecutwfc=80, ecutwfc_step=10)

# Review outliers with automatic detection (Z-score 2.0σ)
fit = wf.phase1_fit_result  # Already fitted with automatic detection
review = wf.review_fit_outliers(fit, verbose=True)
```

You'll see output like:
```
⚠ 2 point(s) excluded as outliers:
  1. ecutwfc = 30 Ry  →  ΔE = 2000.00 meV [MAJOR]
  2. ecutwfc = 40 Ry  →  ΔE =  100.00 meV [MODERATE]

📊 IMPACT ON RECOMMENDATION:
  Without removing outliers:  72.3 Ry  (R² = 0.9821)
  With removing outliers:     50.1 Ry  (R² = 0.9998)
  Difference:                 22.2 Ry  ← BIG IMPACT!
```

**Step 2: Check if automatic detection makes sense**
- If removing outliers improves R² significantly → ✅ Trust it
- If removes too much data → adjust outlier_threshold

**Step 3: Optional - Refit with custom meV threshold**
```python
# If you want to be more aggressive
new_fit = wf.refit_with_absolute_meV_threshold(
    ecut_results,
    criteria_tolerances,
    meV_threshold=50,  # Exclude |ΔE| > 50 meV
    verbose=True
)
```

---

## Practical Values by Pseudopotential Type

### Norm-Conserving (NC) Pseudopotentials
- **Typical behavior:** Smooth exponential, few outliers
- **Suggested Z-score threshold:** 2.0σ (default) ✅
- **Fallback meV threshold:** 50-100 meV

### Ultra-Soft (US) Pseudopotentials  
- **Typical behavior:** More scatter, more outliers at low ecut
- **Suggested Z-score threshold:** 1.5σ (more aggressive)
- **Fallback meV threshold:** 20-50 meV

### PAW Pseudopotentials
- **Typical behavior:** Sharp transition region
- **Suggested Z-score threshold:** 2.5σ (more conservative)
- **Fallback meV threshold:** 50-100 meV

---

## Implementation Details

### Auto-detection with Z-score (DEFAULT)

```python
# Already done automatically in _fit_exponential_decay_phase1()
# No code needed! It's automatic.

# But you can customize:
fit = wf._fit_exponential_decay_phase1(
    ecut_results,
    criteria_tolerances,
    verbose=True,
    auto_exclude_outliers=True,      # Enable detection (default)
    outlier_threshold=2.0              # Z-score threshold (default)
)

# For more aggressive:
fit = wf._fit_exponential_decay_phase1(
    ecut_results,
    criteria_tolerances,
    verbose=True,
    auto_exclude_outliers=True,
    outlier_threshold=1.5              # More aggressive
)

# To disable (use all data):
fit = wf._fit_exponential_decay_phase1(
    ecut_results,
    criteria_tolerances,
    verbose=True,
    auto_exclude_outliers=False        # No detection
)
```

### Manual exclusion with meV threshold

```python
# If you know ecutwfc=30,40 are problematic:
new_fit = wf.refit_with_custom_exclusion(
    ecut_results,
    criteria_tolerances,
    excluded_ecutwfc=[30, 40],
    verbose=True
)
```

---

## Decision Tree

```
Does fitted R² > 0.99?
  ├─ YES (excellent fit)
  │  └─ ✅ Trust the automatic outlier detection
  │     No changes needed
  │
  └─ NO (poor fit)
     │
     └─ Are there obvious anomalies at low ecutwfc?
        ├─ YES (ΔE >> expected for convergence region)
        │  └─ Increase outlier_threshold to 1.5σ
        │     or exclude manually via refit_with_custom_exclusion()
        │
        └─ NO (scattered data everywhere)
           └─ Check calculation quality
              Might need to rerun with stricter tolerances
```

---

## Example: Au Bulk with Different Thresholds

**Data collected:**
```
ecutwfc  |  Energy (eV)  |  ΔE from first (meV)
---------|---------------|--------------------
30       |  -9999.0000   |     0
40       |  -9997.5000   |   1500
50       |  -9995.0000   |   2000  ← OUTLIER!
60       |  -9994.5000   |   2500  ← OUTLIER!
...
```

Wait, let me recalculate. The FIRST point should be the reference (0 meV), not subsequent ones.

```
ecutwfc  |  Energy (eV)      |  ΔE from converged (meV)
---------|-------------------|------------------------
30       |  -9800.0000       |   2000  ← HUGE! (basis incomplete)
40       |    -9995.0000     |    100  ← Still significant
50       |    -9999.5000     |     10  ← Converging
60       |    -9999.9000     |      2  ← Good
70       |   -10000.0500     |      1  ← Very converged
```

### With Z-score 2.0σ (automatic)
```
Residuals std = ~50 meV
2.0σ = 100 meV

Outliers excluded:
  - ecutwfc=30 (deviation 2000 meV >> 100 meV)
  
Remaining: 40, 50, 60, 70
Result: ecut_recommended = 50 Ry ✓ Good!
R² = 0.9999
```

### With Z-score 1.5σ (more aggressive)
```
1.5σ = 75 meV

Outliers excluded:
  - ecutwfc=30 (2000 meV >> 75 meV)
  - ecutwfc=40 (100 meV >> 75 meV)
  
Remaining: 50, 60, 70
Result: ecut_recommended = 45 Ry ✓ Even better!
R² = 0.9999
```

### With meV threshold = 100 meV
```
Outliers excluded: ecutwfc=30

Remaining: 40, 50, 60, 70
Result: ecut_recommended = 50 Ry
R² = 0.9998
```

### With meV threshold = 50 meV
```
Outliers excluded: ecutwfc=30, ecutwfc=40

Remaining: 50, 60, 70
Result: ecut_recommended = 47 Ry
R² = 0.9999
```

---

## My Recommendation for Your Workflow

### Quick Answer:
**Use the default Z-score 2.0σ automatic detection.** 

It's:
- ✅ Automatic (no manual choice)
- ✅ Adaptive to your data
- ✅ Statistically sound
- ✅ Shows impact clearly

### If automatic doesn't work well:
1. **Check R² from automatic detection**
   - R² > 0.99 → Trust it ✅
   - R² < 0.95 → Something's wrong, investigate

2. **Look at outliers reported**
   - If it removed only first 1-2 points → probably correct ✅
   - If it removed many scattered points → maybe too aggressive

3. **If needed, adjust manually:**
   ```python
   # More conservative (keep more data)
   fit = wf._fit_exponential_decay_phase1(
       ecut_results, 
       criteria_tolerances,
       outlier_threshold=2.5  # Higher threshold
   )
   
   # More aggressive (remove more data)
   fit = wf._fit_exponential_decay_phase1(
       ecut_results,
       criteria_tolerances, 
       outlier_threshold=1.5  # Lower threshold
   )
   ```

---

## Summary Table

| Strategy | Threshold | Auto? | Best For | Risk |
|----------|-----------|-------|----------|------|
| Z-score default | 2.0σ | ✅ Yes | Most cases | Low - adaptive |
| Z-score aggressive | 1.5σ | ✅ Yes | Noisy data | Medium - may remove good data |
| meV conservative | 100 meV | ❌ Manual | You know the pseudo | Medium - keeps outliers |
| meV moderate | 50 meV | ❌ Manual | Balance needed | Low-Medium |
| meV aggressive | 20 meV | ❌ Manual | Clean pseudos | Medium-High - too strict |
| Manual exclusion | - | ❌ Manual | You inspected data | Low - most control |

---

## When to Be Worried ⚠

You should manually inspect if:

1. **R² drops below 0.95** after outlier removal
   - Might be removing too much
   - Or data quality is poor

2. **All low ecutwfc removed, but convergence not smooth**
   - Points might be fine, just looks messy
   - Try meV threshold instead

3. **Different runs give very different recommendations**
   - Outlier detection is sensitive
   - Use manual exclusion for consistency

---

## Code Examples

### Scenario 1: Default (recommended)
```python
wf = ConvergenceWorkflow(atoms, pseudopotentials)
wf.run_convergence_study()
fit = wf.phase1_fit_result
# Automatic Z-score 2.0σ detection already applied ✅
```

### Scenario 2: Review and adjust
```python
# Run with custom Z-score threshold
fit = wf._fit_exponential_decay_phase1(
    ecut_results,
    criteria_tolerances,
    outlier_threshold=1.5,  # More aggressive
    verbose=True
)

# Check impact
review = wf.review_fit_outliers(fit, verbose=True)

# If satisfied with fit
rec = wf.get_recommendations()
```

### Scenario 3: Manual control (meV)
```python
# You know ecutwfc=30,40 are problematic
new_fit = wf.refit_with_custom_exclusion(
    ecut_results,
    criteria_tolerances,
    excluded_ecutwfc=[30, 40],
    verbose=True
)
```

---

## Final Answer

**For most cases: Default Z-score 2.0σ ✅**
- Automatic detection
- Good balance
- Clear impact reporting

**Adjust to 1.5σ if:** Fit quality poor with default
**Use manual meV if:** You know your pseudopotential behavior
