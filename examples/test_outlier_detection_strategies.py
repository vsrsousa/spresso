#!/usr/bin/env python
"""
Demonstrate three outlier detection strategies for exponential fit.

Strategies:
1. Automatic (Z-score) - DEFAULT
2. Conservative meV threshold (50 meV)
3. Aggressive meV threshold (20 meV)
4. Manual exclusion (explicit control)

This example shows how each affects the exponential fit and recommendations.

Run with:
    python examples/test_outlier_detection_strategies.py
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

print("="*80)
print("OUTLIER DETECTION STRATEGIES - COMPARISON")
print("="*80)

# ============================================================================
# Simulate Au bulk data with outliers (like real case)
# ============================================================================

print("\n1️⃣ SIMULATED DATA (Au bulk with low-ecutwfc outliers)")
print("-"*80)

# Create realistic convergence curve
ecutwfc_all = np.array([30, 40, 50, 60, 70, 80])
E_inf = -19.253789  # Asymptotic energy (eV)
A = -0.0194         # Amplitude (eV)
B = 0.0875          # Decay constant (Ry⁻¹)

# Calculate energies from exponential + noise
energy_clean = E_inf + A * np.exp(-B * ecutwfc_all)

# Add realistic errors (outliers at low ecutwfc)
energy_with_outliers = energy_clean.copy()
energy_with_outliers[0] += 0.02    # ecutwfc=30: +20 meV (pseudopotential basis incomplete)
energy_with_outliers[1] += 0.001   # ecutwfc=40: +1 meV (small residual error)

print(f"\nGenerated energies:")
print(f"{'ecutwfc (Ry)':<15} {'Energy (eV)':<18} {'ΔE from E_inf (meV)':<20}")
print("-"*70)

for ecut, E in zip(ecutwfc_all, energy_with_outliers):
    delta_e = (E - E_inf) * 1000
    print(f"{ecut:<15.0f} {E:<18.8f} {delta_e:<20.2f}")

# ============================================================================
# STRATEGY 1: Automatic (Z-score) - DEFAULT
# ============================================================================

print("\n" + "="*80)
print("STRATEGY 1: AUTOMATIC (Z-SCORE) - DEFAULT ✅")
print("="*80)

def exponential_decay(x, E_inf, A, B):
    return E_inf + A * np.exp(-B * x)

# Initial fit with all data
popt_all, _ = curve_fit(
    exponential_decay,
    ecutwfc_all,
    energy_with_outliers,
    p0=[E_inf, A, B],
    maxfev=10000
)

E_inf_all, A_all, B_all = popt_all
residuals_all = energy_with_outliers - exponential_decay(ecutwfc_all, *popt_all)
residuals_std = np.std(residuals_all)
z_scores = np.abs(residuals_all) / residuals_std

threshold_z = 2.0
outlier_mask = z_scores > threshold_z
excluded_auto = ecutwfc_all[outlier_mask]

print(f"\nZ-score analysis (threshold = {threshold_z}σ):")
print(f"{'ecutwfc':<12} {'Residual (meV)':<18} {'Z-score':<12} {'Status':<15}")
print("-"*70)

for i, ecut in enumerate(ecutwfc_all):
    res_meV = residuals_all[i] * 1000
    z = z_scores[i]
    status = "❌ OUTLIER" if outlier_mask[i] else "✓ Keep"
    print(f"{ecut:<12.0f} {res_meV:<18.2f} {z:<12.2f} {status:<15}")

print(f"\nResidual std: {residuals_std*1000:.2f} meV")
print(f"Threshold: {threshold_z}σ = {threshold_z*residuals_std*1000:.2f} meV")

# Refit without outliers
if np.any(outlier_mask) and np.sum(~outlier_mask) >= 3:
    ecutwfc_fit_auto = ecutwfc_all[~outlier_mask]
    energy_fit_auto = energy_with_outliers[~outlier_mask]
    
    popt_auto, _ = curve_fit(
        exponential_decay,
        ecutwfc_fit_auto,
        energy_fit_auto,
        p0=[E_inf, A, B],
        maxfev=10000
    )
    
    E_inf_auto, A_auto, B_auto = popt_auto
    residuals_auto = energy_fit_auto - exponential_decay(ecutwfc_fit_auto, *popt_auto)
    r2_auto = 1 - np.sum(residuals_auto**2) / np.sum((energy_fit_auto - np.mean(energy_fit_auto))**2)
    
    tolerance = 0.001  # 1 meV tolerance
    ecut_rec_auto = -np.log(tolerance / abs(A_auto)) / B_auto if tolerance < abs(A_auto) else ecutwfc_fit_auto[0]
    
    print(f"\n✓ Refitted with {len(ecutwfc_fit_auto)} points (excluded {len(excluded_auto)})")
    print(f"  E_inf = {E_inf_auto:.8f} eV")
    print(f"  A = {A_auto:.8f} eV")
    print(f"  B = {B_auto:.6f} Ry⁻¹")
    print(f"  R² = {r2_auto:.6f}")
    print(f"  Recommended ecutwfc: {ecut_rec_auto:.1f} Ry")

# ============================================================================
# STRATEGY 2: Conservative meV threshold (50 meV)
# ============================================================================

print("\n" + "="*80)
print("STRATEGY 2: CONSERVATIVE meV THRESHOLD (50 meV)")
print("="*80)

meV_threshold_conservative = 50.0
E_highest = energy_with_outliers[-1]  # Use highest ecutwfc as reference
deviations_meV = (energy_with_outliers - E_highest) * 1000

excluded_meV_conservative = ecutwfc_all[np.abs(deviations_meV) > meV_threshold_conservative]

print(f"\nmeV threshold: {meV_threshold_conservative} meV")
print(f"{'ecutwfc':<12} {'ΔE from high (meV)':<20} {'Status':<15}")
print("-"*70)

for i, ecut in enumerate(ecutwfc_all):
    dev = deviations_meV[i]
    status = "❌ EXCLUDED" if ecut in excluded_meV_conservative else "✓ Kept"
    print(f"{ecut:<12.0f} {dev:<20.2f} {status:<15}")

# Refit
mask_conservative = ~np.isin(ecutwfc_all, excluded_meV_conservative)
if np.sum(mask_conservative) >= 3:
    ecutwfc_fit_cons = ecutwfc_all[mask_conservative]
    energy_fit_cons = energy_with_outliers[mask_conservative]
    
    popt_cons, _ = curve_fit(
        exponential_decay,
        ecutwfc_fit_cons,
        energy_fit_cons,
        p0=[E_inf, A, B],
        maxfev=10000
    )
    
    E_inf_cons, A_cons, B_cons = popt_cons
    residuals_cons = energy_fit_cons - exponential_decay(ecutwfc_fit_cons, *popt_cons)
    r2_cons = 1 - np.sum(residuals_cons**2) / np.sum((energy_fit_cons - np.mean(energy_fit_cons))**2)
    
    ecut_rec_cons = -np.log(tolerance / abs(A_cons)) / B_cons if tolerance < abs(A_cons) else ecutwfc_fit_cons[0]
    
    print(f"\n✓ Refitted with {len(ecutwfc_fit_cons)} points (excluded {len(excluded_meV_conservative)})")
    print(f"  E_inf = {E_inf_cons:.8f} eV")
    print(f"  A = {A_cons:.8f} eV")
    print(f"  B = {B_cons:.6f} Ry⁻¹")
    print(f"  R² = {r2_cons:.6f}")
    print(f"  Recommended ecutwfc: {ecut_rec_cons:.1f} Ry")

# ============================================================================
# STRATEGY 3: Aggressive meV threshold (20 meV)
# ============================================================================

print("\n" + "="*80)
print("STRATEGY 3: AGGRESSIVE meV THRESHOLD (20 meV)")
print("="*80)

meV_threshold_aggressive = 20.0
excluded_meV_aggressive = ecutwfc_all[np.abs(deviations_meV) > meV_threshold_aggressive]

print(f"\nmeV threshold: {meV_threshold_aggressive} meV")
print(f"{'ecutwfc':<12} {'ΔE from high (meV)':<20} {'Status':<15}")
print("-"*70)

for i, ecut in enumerate(ecutwfc_all):
    dev = deviations_meV[i]
    status = "❌ EXCLUDED" if ecut in excluded_meV_aggressive else "✓ Kept"
    print(f"{ecut:<12.0f} {dev:<20.2f} {status:<15}")

# Refit
mask_aggressive = ~np.isin(ecutwfc_all, excluded_meV_aggressive)
if np.sum(mask_aggressive) >= 3:
    ecutwfc_fit_agg = ecutwfc_all[mask_aggressive]
    energy_fit_agg = energy_with_outliers[mask_aggressive]
    
    popt_agg, _ = curve_fit(
        exponential_decay,
        ecutwfc_fit_agg,
        energy_fit_agg,
        p0=[E_inf, A, B],
        maxfev=10000
    )
    
    E_inf_agg, A_agg, B_agg = popt_agg
    residuals_agg = energy_fit_agg - exponential_decay(ecutwfc_fit_agg, *popt_agg)
    r2_agg = 1 - np.sum(residuals_agg**2) / np.sum((energy_fit_agg - np.mean(energy_fit_agg))**2)
    
    ecut_rec_agg = -np.log(tolerance / abs(A_agg)) / B_agg if tolerance < abs(A_agg) else ecutwfc_fit_agg[0]
    
    print(f"\n✓ Refitted with {len(ecutwfc_fit_agg)} points (excluded {len(excluded_meV_aggressive)})")
    print(f"  E_inf = {E_inf_agg:.8f} eV")
    print(f"  A = {A_agg:.8f} eV")
    print(f"  B = {B_agg:.6f} Ry⁻¹")
    print(f"  R² = {r2_agg:.6f}")
    print(f"  Recommended ecutwfc: {ecut_rec_agg:.1f} Ry")

# ============================================================================
# STRATEGY 4: Manual exclusion
# ============================================================================

print("\n" + "="*80)
print("STRATEGY 4: MANUAL EXCLUSION (exclude ecutwfc=30)")
print("="*80)

excluded_manual = [30]
mask_manual = ~np.isin(ecutwfc_all, excluded_manual)

print(f"\nManually excluding: {excluded_manual}")
print(f"Keeping: {ecutwfc_all[mask_manual].tolist()}")

if np.sum(mask_manual) >= 3:
    ecutwfc_fit_manual = ecutwfc_all[mask_manual]
    energy_fit_manual = energy_with_outliers[mask_manual]
    
    popt_manual, _ = curve_fit(
        exponential_decay,
        ecutwfc_fit_manual,
        energy_fit_manual,
        p0=[E_inf, A, B],
        maxfev=10000
    )
    
    E_inf_manual, A_manual, B_manual = popt_manual
    residuals_manual = energy_fit_manual - exponential_decay(ecutwfc_fit_manual, *popt_manual)
    r2_manual = 1 - np.sum(residuals_manual**2) / np.sum((energy_fit_manual - np.mean(energy_fit_manual))**2)
    
    ecut_rec_manual = -np.log(tolerance / abs(A_manual)) / B_manual if tolerance < abs(A_manual) else ecutwfc_fit_manual[0]
    
    print(f"\n✓ Refitted with {len(ecutwfc_fit_manual)} points (excluded {len(excluded_manual)})")
    print(f"  E_inf = {E_inf_manual:.8f} eV")
    print(f"  A = {A_manual:.8f} eV")
    print(f"  B = {B_manual:.6f} Ry⁻¹")
    print(f"  R² = {r2_manual:.6f}")
    print(f"  Recommended ecutwfc: {ecut_rec_manual:.1f} Ry")

# ============================================================================
# COMPARISON TABLE
# ============================================================================

print("\n" + "="*80)
print("COMPARISON TABLE")
print("="*80)

comparison_data = [
    ("No filtering", len(ecutwfc_all), r2_all := 1 - np.sum((energy_with_outliers - exponential_decay(ecutwfc_all, *popt_all))**2) / np.sum((energy_with_outliers - np.mean(energy_with_outliers))**2), 
     ecut_all := -np.log(tolerance / abs(A_all)) / B_all if tolerance < abs(A_all) else ecutwfc_all[-1]),
    ("Z-score 2.0σ (auto)", len(ecutwfc_fit_auto), r2_auto, ecut_rec_auto),
    ("meV: 50 meV", len(ecutwfc_fit_cons), r2_cons, ecut_rec_cons),
    ("meV: 20 meV", len(ecutwfc_fit_agg), r2_agg, ecut_rec_agg),
    ("Manual: exclude 30", len(ecutwfc_fit_manual), r2_manual, ecut_rec_manual),
]

print(f"\n{'Strategy':<25} {'Points Used':<15} {'R²':<12} {'Recommendation':<15}")
print("-"*70)

for name, n_points, r2, ecut_rec in comparison_data:
    print(f"{name:<25} {n_points:<15} {r2:<12.6f} {ecut_rec:<15.1f} Ry")

# ============================================================================
# VISUALIZATION
# ============================================================================

print("\n" + "="*80)
print("Creating visualization...")
print("="*80)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: All data with automatic detection
ax = axes[0, 0]
ecut_smooth = np.linspace(30, 100, 200)

# Data with outliers marked
not_outlier = ~np.isin(ecutwfc_all, excluded_auto)
ax.scatter(ecutwfc_all[not_outlier], energy_with_outliers[not_outlier], 
           s=100, color='blue', label='Kept', zorder=5)
if len(excluded_auto) > 0:
    ax.scatter(ecutwfc_all[np.isin(ecutwfc_all, excluded_auto)], 
               energy_with_outliers[np.isin(ecutwfc_all, excluded_auto)],
               s=100, color='red', marker='x', linewidth=2, label='Outlier (Z-score)', zorder=5)

# Fit lines
E_smooth_auto = exponential_decay(ecut_smooth, E_inf_auto, A_auto, B_auto)
ax.plot(ecut_smooth, E_smooth_auto, 'b-', linewidth=2, label=f'Fit (R²={r2_auto:.4f})')

ax.axvline(ecut_rec_auto, color='green', linestyle='--', linewidth=2, 
           label=f'Recommended: {ecut_rec_auto:.1f} Ry')

ax.set_xlabel('ecutwfc (Ry)', fontsize=11)
ax.set_ylabel('Energy (eV)', fontsize=11)
ax.set_title('Strategy 1: Automatic (Z-score)', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()

# Plot 2: Conservative meV threshold
ax = axes[0, 1]

not_excluded_cons = ~np.isin(ecutwfc_all, excluded_meV_conservative)
ax.scatter(ecutwfc_all[not_excluded_cons], energy_with_outliers[not_excluded_cons],
           s=100, color='blue', label='Kept', zorder=5)
if len(excluded_meV_conservative) > 0:
    ax.scatter(ecutwfc_all[np.isin(ecutwfc_all, excluded_meV_conservative)],
               energy_with_outliers[np.isin(ecutwfc_all, excluded_meV_conservative)],
               s=100, color='orange', marker='x', linewidth=2, label='Outlier (meV)', zorder=5)

E_smooth_cons = exponential_decay(ecut_smooth, E_inf_cons, A_cons, B_cons)
ax.plot(ecut_smooth, E_smooth_cons, 'b-', linewidth=2, label=f'Fit (R²={r2_cons:.4f})')

ax.axvline(ecut_rec_cons, color='green', linestyle='--', linewidth=2,
           label=f'Recommended: {ecut_rec_cons:.1f} Ry')

ax.set_xlabel('ecutwfc (Ry)', fontsize=11)
ax.set_ylabel('Energy (eV)', fontsize=11)
ax.set_title('Strategy 2: Conservative (50 meV)', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()

# Plot 3: Aggressive meV threshold
ax = axes[1, 0]

not_excluded_agg = ~np.isin(ecutwfc_all, excluded_meV_aggressive)
ax.scatter(ecutwfc_all[not_excluded_agg], energy_with_outliers[not_excluded_agg],
           s=100, color='blue', label='Kept', zorder=5)
if len(excluded_meV_aggressive) > 0:
    ax.scatter(ecutwfc_all[np.isin(ecutwfc_all, excluded_meV_aggressive)],
               energy_with_outliers[np.isin(ecutwfc_all, excluded_meV_aggressive)],
               s=100, color='red', marker='x', linewidth=2, label='Outlier (meV)', zorder=5)

E_smooth_agg = exponential_decay(ecut_smooth, E_inf_agg, A_agg, B_agg)
ax.plot(ecut_smooth, E_smooth_agg, 'b-', linewidth=2, label=f'Fit (R²={r2_agg:.4f})')

ax.axvline(ecut_rec_agg, color='green', linestyle='--', linewidth=2,
           label=f'Recommended: {ecut_rec_agg:.1f} Ry')

ax.set_xlabel('ecutwfc (Ry)', fontsize=11)
ax.set_ylabel('Energy (eV)', fontsize=11)
ax.set_title('Strategy 3: Aggressive (20 meV)', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()

# Plot 4: Comparison of all recommendations
ax = axes[1, 1]

strategies = ['No filter', 'Z-score', 'meV 50', 'meV 20', 'Manual']
recommendations = [ecut_all, ecut_rec_auto, ecut_rec_cons, ecut_rec_agg, ecut_rec_manual]
colors = ['gray', 'blue', 'orange', 'red', 'green']

bars = ax.bar(strategies, recommendations, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)

# Add value labels on bars
for bar, val in zip(bars, recommendations):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height,
            f'{val:.1f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

ax.set_ylabel('Recommended ecutwfc (Ry)', fontsize=11)
ax.set_title('Comparison of Recommendations', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')
ax.set_ylim([0, max(recommendations)*1.15])

plt.tight_layout()
plt.savefig('outlier_detection_comparison.png', dpi=150)
print("✓ Saved: outlier_detection_comparison.png")
plt.show()

# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "="*80)
print("SUMMARY & RECOMMENDATIONS")
print("="*80)

print(f"""
For this simulated Au data with outliers:

1. NO FILTERING
   → ecut = {ecut_all:.1f} Ry (INFLATED due to outliers)
   → R² = {r2_all:.6f} (poor fit quality)
   
2. AUTO Z-SCORE (2.0σ) ✅ RECOMMENDED
   → ecut = {ecut_rec_auto:.1f} Ry
   → R² = {r2_auto:.6f}
   → Why: Automatic, adaptive, statistically sound
   
3. CONSERVATIVE meV (50 meV)
   → ecut = {ecut_rec_cons:.1f} Ry
   → R² = {r2_cons:.6f}
   → Excluded: {list(excluded_meV_conservative)}
   
4. AGGRESSIVE meV (20 meV)
   → ecut = {ecut_rec_agg:.1f} Ry
   → R² = {r2_agg:.6f}
   → Excluded: {list(excluded_meV_aggressive)}
   
5. MANUAL EXCLUSION
   → ecut = {ecut_rec_manual:.1f} Ry
   → R² = {r2_manual:.6f}
   → Excluded: {excluded_manual}

CONCLUSION:
- All filtering strategies give similar and more reasonable results
- Z-score (automatic) requires no manual tuning
- meV thresholds offer explicit control if needed
- Default recommendation: Use Z-score 2.0σ (already active in code!)
""")

print("="*80)
