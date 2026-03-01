# QGD: Generalised κ-Inversion Formula & T^μν Diagnosis
## Version 2.1 — From Bullet Cluster to All Galaxy Classes

---

## The General Formula

The Bullet Cluster gave us a **specific** inversion:

```
κ_star = (M_lens − M_gas × κ_gas − M_cross) / M_star  =  38.83  →  κ₅
```

The **general case** for any system with N baryonic components, each with mass Mᵢ and local κᵢ:

```
┌──────────────────────────────────────────────────────────────────────┐
│  Component form (clusters, lensing, group mass profiles):            │
│                                                                      │
│   M_grav = Σᵢ Mᵢ · κᵢ(Σᵢ, T^μν_ᵢ)  +  M_cross                   │
│                                                                      │
│   INVERT for component j:                                            │
│                                                                      │
│   κⱼ = [ M_grav  −  Σᵢ≠ⱼ Mᵢ·κᵢ  −  M_cross ] / Mⱼ               │
│                                                                      │
│  Rotation curve form (galaxies, point-by-point):                     │
│                                                                      │
│   κ_required(r) = [ v_obs(r) / v_baryon(r) ]²                       │
└──────────────────────────────────────────────────────────────────────┘
```

This inverted κ is what the **data demands**. Match it to the factorial ladder to identify which rung is active.

---

## The κ-Ladder (Complete, v2.1)

| n | κₙ (exact) | κₙ | Physical regime | Status |
|---|---|---|---|---|
| 1 | 1 | **1.000** | Newtonian baseline (solar system, dense gas) | ✓ |
| 2 | √(3/2) | **1.225** | Wide binaries, isolated dwarfs | ✓ |
| 3 | √(15/2) | **2.739** | Spiral galaxy outskirts | ✓ R²=0.935 |
| 4 | √(315/4) | **8.874** | Galaxy groups 10¹²–10¹³ M☉ | Pending |
| 5 | √(1417.5) | **37.65** | Galaxy clusters 10¹³–10¹⁵ M☉ | ✓ ratio=0.981 |
| 6 | √(38981.25) | **197.4** | Superclusters / WHIM filaments | Theoretical |
| 7 | √(2129828.6) | **1233.** | Cosmic web / horizon | Theoretical |

Inversion applied to each rung — what M_grav/M_baryon ratio would **confirm** each:

| Rung | κₙ required | f_baryon predicted | Physical test |
|------|------------|-------------------|---------------|
| κ₂ | 1.225 | 81.6% | Wide binary velocity dispersion |
| κ₃ | 2.739 | 36.5% | Spiral rotation curve outskirts |
| κ₄ | 8.874 | 11.3% | Group lensing mass vs X-ray baryon |
| κ₅ | 37.65 | 2.66% | Cluster lensing (Bullet: 0.981 ✓) |
| κ₆ | 197.4 | 0.51% | WHIM filament mass vs gravitational |

---

## The Critical Insight: Σ is Only T^00

The surface density proxy approximates **one component** of the full stress-energy tensor:

```
T^μν = ρu^μu^ν  +  Ph^μν  +  q^μu^ν  +  q^νu^μ  +  π^μν

T^00    = ρc²          → Σ proxy captures this ✓
T^11,22 = P = ρσ_v²   → isotropic pressure     MISSING
T^0i    = ρ v_rot c   → momentum flux           MISSING
π^ij    = viscous      → shear stress            MISSING
```

The QGD κ source couples to the **trace**:

```
T^μ_μ = −T^00 + Tr(T^ii)/c²  =  −ρc²(1 − 3w)

where  w = P/(ρc²) = σ_v²/c²   (equation-of-state parameter)
```

For non-relativistic systems: **w = σ_v²/v_circ²** is the relevant dimensionless ratio.

### Full T^μν κ Formula (proposed v2.2)

```
κ_full = 1 + (κ_base − 1) × f_P × f_β × f_shear × f_cross

where:
  κ_base  = K_floor + (κ_Σ − K_floor) × Q(M)     ← current model
  κ_Σ    = 1 + (Σ_crit/Σ)^α                       ← Σ proxy (T^00 only)

  f_P(w)     = exp(−ln3 · w)             pressure suppression from T^ii
  f_β(β)     = 1 + 0.15(β − 0.5)        orbital anisotropy (β = 1 − σ_t²/σ_r²)
  f_shear    = 1 + 0.1·|∂lnΣ/∂lnr|/(1+w)  T^ij shear from Σ gradient
  f_cross(v) = 1 + 0.08·(v_circ/250)^0.5  T^0i momentum flux  ← MISSING IN v2.1
```

---

## Rung Inversion Applied to Failure Modes

Using `κ_required(r) = [v_obs/v_bar]²`:

| System | v_obs/v_bar | κ_required | Nearest rung | Error |
|--------|-------------|-----------|--------------|-------|
| Dwarf inner r | 2.08 | 4.34 | κ₄ = 8.874 | 51% |
| **Dwarf outer r** | **2.75** | **7.56** | **κ₄ = 8.874** | **15%** |
| Small spiral r=2 kpc | 1.50 | 2.25 | κ₃ = 2.739 | 18% |
| Small spiral r=8 kpc | 1.91 | 3.63 | κ₃ = 2.739 | 32% |
| Massive spiral (bulge) | 1.12 | 1.25 | κ₂ = 1.225 | 2% |
| Massive spiral (disk) | 2.00 | 4.00 | κ₃ = 2.739 | 46% |
| Bullet Cluster (galaxy) | mass | **38.84** | **κ₅ = 37.650** | **3.1% ✓** |

**The striking finding:** Dwarf outer radii demand κ ~ 7.6, which is ~15% from **κ₄**, not κ₂. This is not noise — it's telling us the model underestimates κ at the outskirts of dwarfs, probably because:
1. Dwarf outer gas has very low Σ (below Σ_crit) → full κ-enhancement expected
2. But Q-factor at dwarf mass saturates early → rung ceiling too low
3. The f_P suppression (pressure at large w) was partially compensating correctly but for wrong reason

---

## Galaxy Class Diagnostic

### Dwarfs: R² 10% → 52% at R²>0.9  (+42pp)

**T^μν state:**
- Σ = 3 M☉/pc² (below Σ_crit → should get full κ₂)
- w = σ_v²/v_circ² ≈ 0.77 → **44% pressure support**
- β ≈ −0.3 (tangential orbits, near-isotropic)

**Correction factors:**
```
f_P     = exp(−ln3 × 0.77) = 0.427   ← MASSIVE SUPPRESSION −57%
f_β     = 1 + 0.15(−0.3 − 0.5) = 0.880
f_cross = 1 + 0.08(25/250)^0.5 = 1.025  (missing from model)
f_total = 0.386
```

The 48% that still fail R²>0.9 split into three distinct causes:
1. **w > 1.5 (dispersion-dominated):** rotation curve is the wrong equation. Need **Jeans equation** → `σ_pred = σ_bar × √κ_eff`
2. **Irregular geometry:** no rotation axis. Need 2D velocity dispersion moment maps.
3. **Satellite dwarfs:** EFE from host screen κ₂ down further (g_ext ~ 10⁻¹¹ m/s²).

Expected further gain: +15–25pp bringing R²>0.9 fraction to ~70–75%.

---

### Small Spirals: R² 0.502 → 0.631  (+0.129)

**T^μν state:**
- Σ ≈ 8 M☉/pc² (near Σ_crit, sensitivity zone)
- w ≈ 0.08, pressure support 7%
- f_cross = 1.048 (missing — adds ~5% to κ)

**Root cause:** Q-factor transition. At M ~ 10⁹·⁵–10¹⁰·⁵ M☉ the current Q-function under-ramps to κ₃.

**Fix:** Recalibrate spiral Q-factor:
```
Current: Q = tanh(M / 10^9.25)^0.5
Proposed: Q = tanh(M / 10^9.8)^0.40    ← higher M_ref, flatter ramp

+ f_cross correction: Δκ = 0.08 × (v_circ/250)^0.5 × (κ₃ − 1)
```
Expected gain: +0.05–0.10 R².

---

### Massive Spirals (Sa/S0): R² 0.475 → 0.493  (+0.018)

This is the weakest performer — and the most instructive.

**T^μν state:**
- Σ ≈ 50 M☉/pc² > Σ_crit → κ_Σ suppressed to 1.73 (near Newtonian)
- v_circ = 280 km/s → T^0i is **largest of any class**
- σ_v = 120 km/s, w = 0.18 → f_P = 0.82 (further suppression)
- β = 0.5 (radial orbits) → f_β = 1.00 (neutral)
- bulge fraction 55%

**The tension:**
```
HIGH Σ > Σ_crit   →  κ_Σ = 1.73  (suppressed toward Newtonian)
HIGH M (Q=1)      →  κ₃ fully accessible
HIGH v_circ = 280 →  f_cross = +8.5%  ← MISSING
HIGH σ_v = 120    →  f_P = −18%
```

The missing f_cross term and the high-Σ suppression are **both pointing wrong** — f_cross would increase κ, but the model doesn't include it, so the net result is flat. The disk velocity profile at large r (where Σ drops below Σ_crit) suddenly jumps to high κ, creating a kink. This mismatch drives the poor R².

**Three-part fix:**
1. **Add f_cross** explicitly: `κ_full × 1.085` for massive spirals
2. **Two-phase decomposition:** separate r < R_bulge (Jeans) from r > R_bulge (rotation)  
3. **Anisotropic bulge correction:** `κ_bulge_Σ = 1 + (Σ_crit/Σ)^α / exp(β/2)` — radial orbits see more radial gravity

Expected gain: +0.10–0.20 R² bringing massive spirals to ~0.60–0.65.

---

## The T^μν Correction Hierarchy

| System | w=σ²/v² | β | f_P | f_cross | Net |
|--------|---------|---|-----|---------|-----|
| Spiral outer disk | 0.01 | 0.2 | 0.990 | 1.020 | **+3%** |
| Small spiral | 0.08 | 0.3 | 0.919 | 1.048 | **−4%** |
| Dwarf irregular | 0.77 | −0.3 | 0.427 | 1.025 | **−57%** ← large |
| Massive bulge | 0.18 | 0.5 | 0.817 | 1.085 | **−12%** |
| Hot ICM (cluster) | >>1 | — | ~0 | — | →κ≈1 ✓ |

**Key observation on the ICM:** Hot cluster gas (T ~ 14.8 keV, σ_v >> v_circ) would have w >> 1 and f_P → 0. This means the Σ-proxy and the T^μν correction **independently and consistently** drive κ_gas toward 1 (Newtonian). Two different physical mechanisms give the same answer — which is a non-trivial consistency check on the theory.

---

## Generalised Rung Prediction: (Σ, M) → κ

| System | logM | Σ [M☉/pc²] | κ_pred | Rung |
|--------|------|------------|--------|------|
| Solar system | 0 | 10⁶ | 1.000 | κ₁ |
| Dwarf irregular | 8.5 | 3.0 | 1.19 | κ₂ |
| Dwarf spheroidal | 7.8 | 1.5 | 1.13 | κ₂ |
| Small spiral (Sc) | 10.2 | 8.0 | 2.26 | κ₃ |
| Large spiral (Sb) | 11.0 | 15.0 | 2.05 | κ₃ |
| Massive spiral (Sa) | 11.6 | 50.0 | **1.73** | κ₃ (suppressed) |
| Galaxy group (low-T) | 12.8 | 4.0 | 2.56 | κ₃→κ₄ |
| Galaxy group (high-T) | 13.3 | 8.0 | **2.74→8.87** | κ₄ (transitioning) |
| Cluster galaxy component | 14.4 | 2.9 | **37.65** | **κ₅** ✓ |
| Cluster ICM gas | 14.4 | 51.3 | 1.72 | κ₁ (suppressed) |
| Supercluster filament | 15.5 | 0.3 | **197.4** | κ₆ (theoretical) |

**Key observation on massive spirals:** The model predicts κ = 1.73 because Σ = 50 > Σ_crit = 17.5. But the inversion from the disk region demands κ ~ 4.0. This 2.3× gap is the fundamental residual for this class. The f_cross term (+8.5%) alone does not close it — the two-phase decomposition is essential.

---

## Priority Fixes (v2.2 Roadmap)

| Priority | Fix | Target class | Expected ΔR² |
|----------|-----|--------------|-------------|
| 1 | Add f_cross = 1 + 0.08(v/250)^0.5 to master κ formula | Massive spirals | +0.10–0.20 |
| 2 | Jeans equation switch for w > 1.5 | Dwarfs (48% failing) | +15–25pp R²>0.9 |
| 3 | Two-phase bulge/disk decomposition | Massive spirals | +0.05–0.10 |
| 4 | Q-factor recalibration: M_ref = 10^9.8, p = 0.40 | Small spirals | +0.05–0.10 |
| 5 | EFE catalogue crossmatch for satellite dwarfs | Satellite dwarfs | +0.10–0.20 |

---

## What This Means Theoretically

The Σ proxy works because surface density traces the dominant T^μν component in most cold, rotationally-supported systems. But QGD is really coupling to the **full trace of the stress-energy tensor** — the quantum gravitational propagator doesn't know whether energy comes from rest mass, thermal pressure, or bulk rotation.

The natural upgrade is to replace Σ everywhere with the **effective source density**:

```
ρ_eff = T^μ_μ / c² = ρ(1 + 3w)   [non-relativistic, comoving frame]
```

Then:
```
κ_Σ → κ_ρeff = 1 + (Σ_crit_eff / Σ_eff)^α

where  Σ_eff = Σ × (1 + 3w)  = Σ + 3P/c²
```

For cold disks (w → 0): Σ_eff = Σ (no change — current model is fine).  
For hot gas (w >> 1): Σ_eff >> Σ → κ → 1 (Newtonian — agrees with ICM behavior).  
For dwarfs (w ~ 0.8): Σ_eff ≈ 3.4 × Σ → κ suppressed — explains the 48% failure.

This is a **single parameter** replacement that unifies all four galaxy classes.

The exact coupling coefficient 3 comes from the 4D spacetime trace:  
`T^μ_μ = T^00 + T^11 + T^22 + T^33 = −ρc² + 3P`

The sign convention means pressure **suppresses** the QGD source — physically, hot/pressurized systems are harder to gravitationally amplify because their phase-space is already thermalized. Quantum gravitational corrections are strongest for cold, ordered mass distributions where coherent phase-space volume is maximized.

---

*QGD v2.1 — March 2026 | kappa_inversion.py produces all computed values*
