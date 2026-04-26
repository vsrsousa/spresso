# 🎯 OUTLIER DETECTION: QUICK DECISION GUIDE

## The Problem You Identified ✅

```
Au bulk ecutwfc convergence:
  ecutwfc = 30 Ry  → ΔE = 2000 meV  ← OUTLIER (basis incomplete)
  ecutwfc = 40 Ry  → ΔE =  100 meV  ← OUTLIER (still poor)
  ecutwfc = 50 Ry  → ΔE =   10 meV  ← Good (convergence region)
  ecutwfc = 60 Ry  → ΔE =    2 meV  ← Excellent

Result if keep all:  ecut_recommended = 70 Ry ❌ WRONG! Too high
Result if remove outliers: ecut_recommended = 50 Ry ✅ CORRECT!
```

---

## Three Practical Thresholds

| Value | Use Case | Example |
|-------|----------|---------|
| **20 meV** | Ultra-aggressive | `Only for clean pseudos like PAW` |
| **50 meV** | Balanced (RECOMMENDED) | `Most norm-conserving pseudos` |
| **100 meV** | Conservative | `If data is very noisy` |

---

## For You (Au bulk)

**DEFAULT**: Use automatic Z-score (2.0σ) 
→ Already active, no configuration needed! ✅

**IF NOT WORKING**: Try meV threshold
```python
fit = wf.refit_with_meV_threshold(
    ecut_results,
    criteria_tolerances,
    meV_threshold=50.0,  # ← Recommended starting point
    verbose=True
)
```

---

## Decision Tree (1 minute)

```
┌─────────────────────────────────────────────┐
│ DID AUTOMATIC FIT LOOK GOOD? (R² > 0.99)    │
└─────────────────────────────────────────────┘
           │
      ┌────┴────┐
      │ NO      │ YES
      │         │
      ▼         └──→ ✅ USE IT!
   Check output    No changes needed
      │
      ├─→ Many outliers removed? (>3 points)
      │   ├─ YES: Try higher threshold (2.5σ)
      │   │       Or use meV_threshold=100
      │   │
      │   └─ NO: Try lower threshold (1.5σ)
      │         Or use meV_threshold=20
      │
      └─→ Too many data points removed?
          └─ Use meV approach (more control)
             Try 50 meV first, adjust up/down
```

---

## Three Methods in Code

### Method 1: AUTOMATIC (Default)
```python
# Already running! No code needed.
fit = wf.phase1_fit_result

# To review what was excluded:
review = wf.review_fit_outliers(fit, verbose=True)
```

### Method 2: ABSOLUTE meV (Practical)
```python
# Try 50 meV first
fit = wf.refit_with_meV_threshold(
    ecut_results,
    criteria_tolerances,
    meV_threshold=50.0,   # ← MODIFY THIS
    verbose=True
)

# Adjust up or down based on output:
# - Too many excluded? Increase to 75 or 100 meV
# - Too few excluded? Decrease to 20 or 30 meV
```

### Method 3: MANUAL (Maximum Control)
```python
# If you saw ecutwfc=30,40 are bad:
fit = wf.refit_with_custom_exclusion(
    ecut_results,
    criteria_tolerances,
    excluded_ecutwfc=[30, 40],
    verbose=True
)
```

---

## The Three Thresholds Explained

### 20 meV (Aggressive)
```
Effect: Removes most low-ecutwfc noise
Good for: Clean PAW pseudopotentials
Risk: Might remove valid data if pseudo is noisy
Use when: You trust your data quality
```

### 50 meV (Balanced) ⭐ RECOMMENDED
```
Effect: Removes obvious outliers, keeps edge cases
Good for: Most norm-conserving pseudos
Risk: Low
Use when: Not sure what to use (start here!)
```

### 100 meV (Conservative)
```
Effect: Only removes severe outliers
Good for: Noisy calculations or soft pseudos
Risk: Might keep bad data
Use when: Data quality is uncertain
```

---

## Real World Scenario

**You collect:**
```
ecutwfc:  30,  40,  50,  60,  70
Energy:  -19.23, -19.24, -19.25, -19.253, -19.2535
ΔE:    2000, 1000,  200,   100,    50 meV
```

**With 50 meV threshold:**
```
ecutwfc=30 (2000 meV > 50) ❌ EXCLUDED
ecutwfc=40 (1000 meV > 50) ❌ EXCLUDED  
ecutwfc=50 (200 meV > 50)  ❌ EXCLUDED
ecutwfc=60 (100 meV > 50)  ❌ EXCLUDED
ecutwfc=70 (50 meV ≤ 50)   ✓ KEPT

Only one point left! Need more data...
```

**Try 100 meV threshold:**
```
ecutwfc=30 (2000 > 100) ❌ EXCLUDED
ecutwfc=40 (1000 > 100) ❌ EXCLUDED
ecutwfc=50 (200 > 100)  ❌ EXCLUDED
ecutwfc=60 (100 ≤ 100)  ✓ KEPT
ecutwfc=70 (50 ≤ 100)   ✓ KEPT

Now have 2 points + 1 more and it works!
```

---

## What Happens Behind the Scenes

### Z-score Automatic (Default)

```python
# Step 1: Fit with ALL data
residuals_all = measured - predicted_from_all_data
std_residuals = standard_deviation(residuals_all)

# Step 2: Find outliers
for each point:
    z_score = |residual| / std_residuals
    if z_score > 2.0:  # threshold
        mark as outlier

# Step 3: Refit without outliers
residuals_clean = measured[not_outliers] - predicted
compute_R_squared and recommendations
```

**Why it works:**
- Adapts to your data's noise level
- Statistical foundation
- No manual tuning

**When it fails:**
- All data is equally noisy (can't find outliers)
- Data has systematic offset (all too high/low)

---

### Absolute meV (Your Control)

```python
# Set reference (most converged, usually highest ecutwfc)
E_ref = energy_at_highest_ecutwfc

# For each ecutwfc:
ΔE = |energy - E_ref|
if ΔE > threshold_meV:
    exclude it

# Refit without excluded points
```

**Why it works:**
- Intuitive (you control the threshold)
- Reproducible (same threshold = same result)
- Easy to explain

**When to use:**
- You know your pseudopotential
- You want explicit control
- Results are consistent across runs

---

## Summary: Which to Use?

| Situation | Strategy | Threshold |
|-----------|----------|-----------|
| **"Just run it"** | Auto Z-score | 2.0σ (default) |
| **"Want to tune"** | meV absolute | 50 meV |
| **"Aggressive"** | meV absolute | 20 meV |
| **"Conservative"** | meV absolute | 100 meV |
| **"I know what's bad"** | Manual | [list ecutwfc values] |

---

## Test the Strategies

Run the provided example to see all three in action:

```bash
python examples/test_outlier_detection_strategies.py
```

Output shows:
- Which points each strategy excludes
- How recommendation changes
- R² improvement
- Visual comparison

---

## Key Insight 💡

You identified something **crucial**: Low ecutwfc outliers can increase the recommended ecutwfc by **20+ Ry** (in your Au case: 50 → 70 Ry).

The automatic detection (Z-score) fixes this **without you needing to do anything**.

But now you have **control** if automatic doesn't work:
- Too aggressive? Use 100 meV threshold
- Too conservative? Use 20 meV threshold
- Custom exclusion? Use manual method

**Status: SOLVED** ✅
- Automatic detection active (default)
- Clear reporting (review_fit_outliers)
- Adjustment methods available (meV, manual)

---

## One More Thing: What If All Points Are Outliers?

If more than 50% of your data gets excluded, something's wrong:

1. **Check your data**
   - Are calculations converged?
   - Are k-points sufficient?
   - Check .out files for warnings

2. **Loosen the threshold**
   ```python
   fit = wf._fit_exponential_decay_phase1(
       ecut_results,
       criteria_tolerances,
       outlier_threshold=3.0  # Much looser
   )
   ```

3. **Use manual approach**
   - Identify ONE clearly bad point
   - Exclude it manually
   - Check R² improvement

---

**Questions? Check:** [OUTLIER_DETECTION_GUIDE.md](OUTLIER_DETECTION_GUIDE.md)
