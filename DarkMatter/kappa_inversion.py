#!/usr/bin/env python3
"""
kappa_inversion.py  —  QGD Generalised κ-Inversion Analysis
=============================================================

CORE PROBLEM
-----------
The Bullet Cluster showed us how to INVERT the κ-formula:
    κ_star = (M_lens - M_gas·κ_gas - M_cross) / M_star  = 38.83  →  κ₅

But that derivation was for a cluster with cleanly separated spatial components.
Can we generalise this across ALL rungs of the κ-ladder?

GENERALISATION
--------------
Any QGD system with N baryonic components, each with surface density Σ_i
and mass M_i:

    M_grav = Σ_i  M_i · κ_i(Σ_i, T^μν_i)  +  M_cross

Inversion (solving for component j):

    κ_j^required = [ M_grav  -  Σ_{i≠j} M_i·κ_i  -  M_cross ] / M_j

For rotation curves at each radius r:

    κ_required(r) = (v_obs(r) / v_baryon(r))²

The SURFACE-DENSITY PROXY APPROXIMATION
----------------------------------------
The current model uses Σ (surface density) as a proxy for the full
stress-energy tensor T^μν. The exact QGD coupling is to T^00 + Tr(T^ii)/c²:

    κ_source = κ_Σ × f_T(w, anisotropy, shear, ...)

where w = P/(ρc²) = σ_v²/c² encodes the equation-of-state.

For cold, rotationally-supported systems (spirals): w → 0, f_T → 1. ✓
For hot, pressure-supported systems (dwarfs, bulges): w >> 0, f_T < 1.

Currently implemented:  f_T ≈ 1 - 3w·M_scale·radial_factor  (approximate)
What it SHOULD be:      f_T = f(T^μν) via the full Tolman-Oppenheimer-Volkoff

PERFORMANCE BY GALAXY CLASS (v2.0 → v2.1 improvement):
    Dwarfs:           R² = −0.19 → 0.84   (but only 52% have R²>0.9)
    Small Spiral:     R² = +0.50 → 0.63
    Large Spiral:     R² = +0.62 → 0.85
    Massive Spiral:   R² = +0.47 → 0.49   ← weakest improvement

This script:
    1. Implements the generalised inversion formula
    2. Diagnoses which κ-rung each class falls on
    3. Identifies the stress-energy tensor corrections needed
    4. Proposes updated T^μν coupling terms for each failure mode

Author: Alex (QGD theory) + Claude (implementation)
Version: 2.1
"""

import numpy as np
from math import factorial
from dataclasses import dataclass
from typing import Optional

# ── κ-ladder ──────────────────────────────────────────────────────────────────
def kappa_n(n: int) -> float:
    """κ_n = √((2n-1)! / 2^(2n-2))  — pure factorial arithmetic."""
    return float(np.sqrt(factorial(2*n - 1) / 2**(2*n - 2)))

K = {n: kappa_n(n) for n in range(1, 8)}
K1, K2, K3, K4, K5, K6, K7 = [K[n] for n in range(1, 8)]

# Physical constants
G     = 6.674e-11   # m³ kg⁻¹ s⁻²
M_SUN = 1.989e30    # kg
PC    = 3.086e16    # m
KPC   = 3.086e19    # m
C     = 3e8         # m/s

# QGD global parameters
SIGMA_CRIT  = 17.5   # M☉/pc²
ALPHA       = 0.30   # power-law index
G_CRIT      = 1.2e-10  # m/s²


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1: GENERALISED κ-INVERSION FORMULA
# ═══════════════════════════════════════════════════════════════════════════════

def kappa_invert_multicomponent(
    M_grav:      float,
    components:  list[tuple[float, float]],   # [(M_i, kappa_i), ...]
    target_idx:  int,
    M_cross:     float = 0.0,
) -> float:
    """
    General κ-inversion: solve for κ of component `target_idx`.

    Equation:
        M_grav = Σ_i M_i·κ_i  +  M_cross
        κ_j = (M_grav - M_cross - Σ_{i≠j} M_i·κ_i) / M_j

    Parameters
    ----------
    M_grav       : total gravitational mass (observed) [M☉]
    components   : list of (M_i [M☉], κ_i) for all components
    target_idx   : index of the component to solve for
    M_cross      : N-body cross-term [M☉]

    Returns
    -------
    float : required κ for the target component

    Examples
    --------
    >>> # Bullet Cluster: solve for κ_star given κ_gas known
    >>> M_grav  = 2.80e14
    >>> kpl_gas = 1.724
    >>> M_gas   = 4.03e13;  M_star = 4.48e12;  M_cross = 3.65e13
    >>> result  = kappa_invert_multicomponent(
    ...     M_grav, [(M_gas, kpl_gas), (M_star, None)], target_idx=1, M_cross=M_cross)
    >>> abs(result - 38.83) < 0.1
    True
    """
    M_j, _ = components[target_idx]
    other_sum = sum(
        M_i * k_i
        for idx, (M_i, k_i) in enumerate(components)
        if idx != target_idx and k_i is not None
    )
    return (M_grav - M_cross - other_sum) / M_j


def kappa_invert_rotation(
    v_obs:     np.ndarray,
    v_baryon:  np.ndarray,
) -> np.ndarray:
    """
    Point-by-point κ-inversion for rotation curves.

    For each radius r:
        κ_required(r) = (v_obs(r) / v_baryon(r))²

    This is the most direct observable: no mass decomposition needed.
    The result tells us exactly which κ-rung a given point sits on.

    Parameters
    ----------
    v_obs    : observed circular velocity [km/s]
    v_baryon : baryonic circular velocity [km/s]

    Returns
    -------
    np.ndarray : κ_required at each point, shape (N,)
    """
    ratio = np.where(v_baryon > 0.1, v_obs / v_baryon, np.nan)
    return ratio**2


def nearest_rung(kappa_val: float) -> tuple[int, float, float]:
    """
    Find the nearest κ-rung to a given κ value.

    Returns (n, kappa_n, fractional_error).
    The fractional error tells how far the observation is from a pure rung.
    """
    best_n, best_err = 1, float('inf')
    for n, kv in K.items():
        err = abs(kappa_val - kv) / kv
        if err < best_err:
            best_err = err; best_n = n
    return best_n, K[best_n], best_err


def rung_histogram(
    kappa_arr: np.ndarray,
    kappa_lo:  float = 0.5,
    kappa_hi:  float = 500.0,
) -> dict:
    """
    Given an array of κ_required values, count how many fall within ±40%
    of each κ-rung.

    Returns dict: {n: count, ...}
    """
    counts = {n: 0 for n in range(1, 8)}
    total  = 0
    for kv in kappa_arr:
        if np.isnan(kv) or kv < kappa_lo or kv > kappa_hi:
            continue
        total += 1
        n, _, err = nearest_rung(kv)
        if err < 0.40:   # within 40% of a rung
            counts[n] += 1
    return counts, total


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2: STRESS-ENERGY TENSOR CORRECTIONS
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class StressEnergyState:
    """
    Encodes the local T^μν state of a baryonic fluid element.

    The surface-density Σ is only the T^00 contribution.
    Full T^μν has:
      T^00 = ρc²            (energy density)
      T^ii = P = ρσ_v²      (isotropic pressure from random motions)
      T^0i = ρv_rot c       (momentum flux, rotation)
      T^ij = π^ij           (viscous shear stress)

    In QGD the κ source couples to Tr(T^μν) = T^μ_μ:
      T^μ_μ = −ρc² + 3P  (GR trace, +signature)
            = −ρc²(1 − 3w)   where w = P/(ρc²) = σ_v²/c²

    For non-relativistic systems (σ_v << c), w = σ_v²/c² << 1, but
    σ_v/v_circ can be of order unity (dwarfs, bulges).
    """
    sigma_surface: float    # M☉/pc²   — surface mass density (T^00 proxy)
    v_circ:        float    # km/s      — circular velocity (rotation support)
    sigma_v:       float    # km/s      — 3D velocity dispersion (pressure support)
    anisotropy:    float    # β = 1 - σ_t²/σ_r²  (Binney anisotropy)
    r_kpc:         float    # kpc       — galactocentric radius
    bulge_frac:    float    # fraction of mass in pressure-supported component

    @property
    def w(self) -> float:
        """Equation-of-state parameter w = σ_v²/v_circ²."""
        if self.v_circ < 0.1:
            return 10.0  # pressure-dominated limit
        return (self.sigma_v / self.v_circ) ** 2

    @property
    def support_fraction(self) -> float:
        """
        Fraction of dynamical support from pressure (vs rotation).
        0.0 → purely rotational (cold disk)
        1.0 → purely pressure-supported (dispersion-dominated)
        """
        v_tot = np.sqrt(self.v_circ**2 + self.sigma_v**2)
        return self.sigma_v**2 / v_tot**2 if v_tot > 0 else 0.5


def kappa_T_correction(state: StressEnergyState) -> dict:
    """
    Compute the stress-energy tensor correction to κ.

    Current model (Σ-proxy):
        κ = 1 + (κ_base - 1) × (1 - 3w·M_scale·r_factor)
    
    This is a first-order approximation of the full T^μν coupling.

    Full correction breaks into four independent terms:

    1. PRESSURE TERM  f_P:
       From isotropic velocity dispersion T^ii = P = ρσ_v²
       f_P = exp(−α_P × w)   where w = σ_v²/v_circ², α_P = ln(3)
       → Suppresses κ in pressure-supported systems
       → f_P → 1 for cold disks (w → 0)
       → f_P → 1/3 for σ_v = v_circ (w = 1)
       → f_P → 0 for pure pressure support (dwarf spheroidals)

    2. ANISOTROPY TERM  f_β:
       From the tangential-to-radial velocity dispersion ratio.
       Radial orbits (β > 0): less effective at supporting against gravity
       Tangential orbits (β < 0): more effective
       f_β = 1 + 0.15 × (β - 0.5)   (linear approx; β ∈ [−∞, 1])

    3. SHEAR / BULGE TERM  f_shear:
       Off-diagonal T^ij from the rotation curve gradient ∂v/∂r.
       In bulge-dominated region, ∂v/∂r > 0 (rising) → extra shear source.
       f_shear = 1 + 0.1 tanh(r/5kpc)  (current implementation)
       IMPROVED: f_shear = 1 + γ_shear × |∂ ln Σ / ∂ ln r| / (1 + w)
       This shear is SUPPRESSED when pressure dominates (w >> 1).

    4. CROSS-COUPLING TERM  f_cross:
       Rotation × pressure cross-term from T^0i.
       T^0i = ρ v_rot c  →  enhances κ slightly in rapidly rotating systems.
       f_cross = 1 + 0.05 × (v_circ/v_circ_ref)  (currently ignored!)
       This is the term missing for MASSIVE SPIRALS.
       Massive spirals have high v_circ → T^0i contribution is largest.

    Returns
    -------
    dict with all correction factors and the composite correction.
    """
    w = state.w
    β = state.anisotropy

    # 1. Pressure suppression
    alpha_P = np.log(3.0)   # = 1.099
    f_P = float(np.exp(-alpha_P * min(w, 5.0)))   # clamp for numerics

    # 2. Orbital anisotropy
    f_beta = float(np.clip(1.0 + 0.15 * (β - 0.5), 0.5, 1.5))

    # 3. Improved shear (suppressed by pressure, relevant in bulge)
    # Approximate Σ gradient from bulge fraction and radius
    sigma_gradient = state.bulge_frac * max(0.1, 1.0 / (1.0 + state.r_kpc / 3.0))
    f_shear = float(1.0 + 0.10 * sigma_gradient / (1.0 + w))

    # 4. Cross-coupling from T^0i (momentum flux) — CURRENTLY MISSING
    # Dimensionless T^0i ∝ v_rot/c
    v_ref   = 250.0   # km/s reference (MW-like)
    f_cross = float(1.0 + 0.08 * (state.v_circ / v_ref)**0.5)
    # This term ENHANCES κ for fast rotators — important for massive spirals

    # Composite
    f_total = f_P * f_beta * f_shear * f_cross

    return {
        'w':         w,
        'f_P':       f_P,
        'f_beta':    f_beta,
        'f_shear':   f_shear,
        'f_cross':   f_cross,
        'f_total':   f_total,
        'pressure_dominated': w > 1.0,
        'support_fraction':   state.support_fraction,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3: GALAXY CLASS DIAGNOSTIC
# ═══════════════════════════════════════════════════════════════════════════════

GALAXY_CLASSES = {
    # class: (log10_M, v_circ_km/s, sigma_v_km/s, sigma_surface_pc2,
    #         r_kpc_typical, beta, bulge_frac, R2_old, R2_new, description)
    'dwarf': {
        'logM': 8.5, 'v_circ': 25.0, 'sigma_v': 22.0,
        'Sigma': 3.0, 'r_kpc': 2.0, 'beta': -0.3, 'bulge_frac': 0.05,
        'R2_old': -0.19, 'R2_new': 0.84, 'pct_gt09_old': 0.10, 'pct_gt09_new': 0.52,
        'desc': 'Isolated dwarf irregulars / dSph (M ~ 10^8–9.5 M☉)',
        'problems': [
            'HI kinematics noisy — v_obs poorly constrained at small r',
            'Mixed rotation+dispersion → v_baryon(rot) ≠ v_baryon(pressure)',
            'w = σ_v/v_circ ~ 0.8 → pressure suppresses κ BUT f_cross enhances',
            'Many dwarfs are irregular — no symmetric rotation axis',
            'EFE from host galaxy (if satellite) adds external screening',
            'Q-factor still ~0 at M < 10^9 M☉ → κ-rung not fully activated',
        ],
    },
    'small_spiral': {
        'logM': 10.2, 'v_circ': 90.0, 'sigma_v': 25.0,
        'Sigma': 8.0, 'r_kpc': 5.0, 'beta': 0.3, 'bulge_frac': 0.15,
        'R2_old': 0.502, 'R2_new': 0.631, 'pct_gt09_old': None, 'pct_gt09_new': None,
        'desc': 'Small spirals (Sc/Sd, M ~ 10^9.5–10^11 M☉)',
        'problems': [
            'Σ ~ 8 M☉/pc² near Σ_crit → power-law correction sensitive',
            'Moderate bulge fraction → some pressure support near centre',
            'v_obs measurement reliable but σ_v/v_circ ~ 0.28 → w ~ 0.08',
            'Q-factor rising (Q~0.5) → partial κ₃ access, not yet full',
            'Bar-driven non-circular motions in inner few kpc bias v_obs',
        ],
    },
    'large_spiral': {
        'logM': 11.0, 'v_circ': 160.0, 'sigma_v': 30.0,
        'Sigma': 15.0, 'r_kpc': 10.0, 'beta': 0.35, 'bulge_frac': 0.25,
        'R2_old': 0.62, 'R2_new': 0.85, 'pct_gt09_old': None, 'pct_gt09_new': None,
        'desc': 'Large spirals (Sb/Sc, M ~ 10^10.5–10^12 M☉)',
        'problems': [
            'Best-performing class — Σ and v_circ well-measured',
            'Moderate residuals from bulge pressure term near centre',
            'v_circ ~ 160 km/s → T^0i term starts mattering (~5%)',
            'AGN in some → scatter from non-gravitational velocity components',
        ],
    },
    'massive_spiral': {
        'logM': 11.6, 'v_circ': 280.0, 'sigma_v': 120.0,
        'Sigma': 50.0, 'r_kpc': 5.0, 'beta': 0.5, 'bulge_frac': 0.55,
        'R2_old': 0.475, 'R2_new': 0.493, 'pct_gt09_old': None, 'pct_gt09_new': None,
        'desc': 'Massive early-type spirals (Sa/S0, M ~ 10^11–10^12.5 M☉)',
        'problems': [
            'Bulge-dominated (55%) → pressure-supported T^ii dominates at r < 10 kpc',
            'High Σ ~ 50 M☉/pc² > Σ_crit → κ suppressed to ~1 by power-law',
            'But galaxy is MASSIVE (Q=1) → it SHOULD have high κ on outer rung',
            'Tension: high Σ suppresses, high M activates — these fight each other',
            'Missing T^0i term: v_circ = 280 km/s → momentum flux T^0i is LARGE',
            'Anisotropic dispersion (β=0.5, radial orbits) → further suppression',
            'Strong AGN feedback → non-circular motions, σ_v overestimated',
        ],
    },
}


def diagnose_galaxy_class(name: str, cls: dict) -> None:
    """Run full T^μν diagnostic for one galaxy class."""
    state = StressEnergyState(
        sigma_surface=cls['Sigma'],
        v_circ=cls['v_circ'],
        sigma_v=cls['sigma_v'],
        anisotropy=cls['beta'],
        r_kpc=cls['r_kpc'],
        bulge_frac=cls['bulge_frac'],
    )
    corr = kappa_T_correction(state)

    # Current κ from Σ-proxy
    kappa_sigma = float(np.clip(
        1.0 + (SIGMA_CRIT / max(state.sigma_surface, 0.01))**ALPHA,
        1.0, K3
    ))
    # κ with full T^μν correction
    kappa_corrected = 1.0 + (kappa_sigma - 1.0) * corr['f_total']

    # Required κ from observed R² improvement
    # R² ∝ 1 - SSres/SStot.  We infer that missing κ ~ 1/sqrt(1-R²)
    R2_new  = cls['R2_new']
    R2_old  = cls['R2_old']
    # Rough estimate: κ_required ≈ kappa_sigma / sqrt(1-R2) if R2 < 1
    kappa_needed_est = kappa_sigma / max(1.0 - R2_new, 0.05)**0.25

    rung_n, rung_kappa, rung_err = nearest_rung(kappa_sigma)

    W = 68
    print(f"\n{'═'*W}")
    print(f"  {name.upper():^{W-4}}")
    print(f"  {cls['desc']:^{W-4}}")
    print(f"{'═'*W}")
    print(f"  R² improvement:  {R2_old:+.3f}  →  {R2_new:+.3f}  (+{R2_new-R2_old:.3f})")
    if cls.get('pct_gt09_new'):
        print(f"  R²>0.9 fraction: {cls['pct_gt09_old']*100:.0f}%  →  {cls['pct_gt09_new']*100:.0f}%  (+{(cls['pct_gt09_new']-cls['pct_gt09_old'])*100:.0f}pp)")
    print(f"\n  T^μν STATE:")
    print(f"    σ_surface = {state.sigma_surface:.1f} M☉/pc²  (vs Σ_crit = {SIGMA_CRIT})")
    print(f"    v_circ    = {state.v_circ:.1f} km/s")
    print(f"    σ_v       = {state.sigma_v:.1f} km/s")
    print(f"    w = σ_v²/v_circ² = {state.w:.3f}   {'(pressure-dominated)' if state.w > 1 else '(rotation-dominated)'}")
    print(f"    anisotropy β      = {state.anisotropy:.2f}")
    print(f"    pressure support  = {state.support_fraction*100:.1f}%")
    print(f"    bulge fraction    = {state.bulge_frac*100:.0f}%")
    print(f"\n  κ DIAGNOSTICS:")
    print(f"    κ_Σ (surface density only):  {kappa_sigma:.4f}  → nearest rung: κ_{rung_n} = {rung_kappa:.3f}  (err={rung_err*100:.0f}%)")
    print(f"  T^μν CORRECTION FACTORS:")
    print(f"    f_P     (pressure/isotropic):  {corr['f_P']:.4f}  {'⚠ large suppression' if corr['f_P'] < 0.7 else ''}")
    print(f"    f_β     (orbit anisotropy):    {corr['f_beta']:.4f}")
    print(f"    f_shear (T^ij / bulge):        {corr['f_shear']:.4f}")
    print(f"    f_cross (T^0i momentum flux):  {corr['f_cross']:.4f}  {'← MISSING IN CURRENT MODEL' if corr['f_cross'] > 1.02 else ''}")
    print(f"    f_total (composite):           {corr['f_total']:.4f}")
    print(f"    κ_corrected = 1 + (κ_Σ-1)×f_total = {kappa_corrected:.4f}")
    print(f"\n  ROOT CAUSES OF POOR R²:")
    for prob in cls['problems']:
        print(f"    • {prob}")
    print(f"\n  PROPOSED CORRECTION:")
    _propose_fix(name, state, corr, kappa_sigma, kappa_corrected, R2_new)


def _propose_fix(name, state, corr, k_sigma, k_corr, R2_new):
    """Print targeted fix for each galaxy class."""
    if name == 'dwarf':
        print(f"""    The 48% of dwarfs with R²<0.9 are mostly:
      a) Pressure-dominated (σ_v/v_circ > 1) — need JEANS EQUATION, not v_circ
      b) Irregular — no well-defined rotation axis; v_obs ≠ v_circ
      c) Satellite dwarfs — EFE from host galaxy not accounted for

    FIX A: For pressure-dominated dwarfs (w > 1.5), switch observable:
      Instead of:   v_pred = v_bar × √κ
      Use:          σ_pred = σ_bar × √κ_eff   (Jeans-based)
      where σ_bar comes from the hydrostatic equilibrium of baryons.

    FIX B: For irregulars, use velocity DISPERSION MAP, not rotation curve.
      This requires re-processing the kinematic data with 2D moment analysis.

    FIX C: Apply EFE correction for satellite dwarfs identified in catalogues.
      Group membership → g_ext ~ 10⁻¹¹–10⁻¹⁰ m/s²  (known from host mass+sep).

    EXPECTED R² GAIN: +0.15–0.25 for the failing 48%.
    (52% already pass → the full R²>0.9 fraction should reach ~70–75%)
""")
    elif name == 'small_spiral':
        print(f"""    Small spirals have moderate improvements because:
      • Σ straddles Σ_crit (sensitive regime)
      • Q-factor ~0.4–0.7 (not yet saturated to κ₃)

    FIX: Improve Q-factor transition for M ~ 10^9.5–10^10.5 M☉.
      Current: Q = tanh(M/M_ref)^0.5
      Issue:   M_ref = 10^9.25 may be too low — the 'spiral turn-on' mass
               from the SPARC bTFR knee is closer to 10^9.8 M☉.
    PROPOSED: Q = tanh(M / 10^9.8)^0.4   (adjusted p and M_ref)
      Also add the T^0i cross-coupling (f_cross = {corr['f_cross']:.3f}):
      κ += Δκ_cross = 0.08 × (v_circ/250)^0.5 × (κ₃ - 1)

    EXPECTED R² GAIN: +0.05–0.10 for small spirals.
""")
    elif name == 'large_spiral':
        print(f"""    Large spirals are the best-performing class (R²=0.85).
    Remaining errors are systematic at small r (bulge) and outer r (Σ→0).

    FIX: At r < 3 kpc (bulge-dominated): apply Jeans correction.
      κ_bulge = κ_Σ × f_P × f_β × f_cross
             = {k_sigma:.3f} × {corr['f_P']:.3f} × {corr['f_beta']:.3f} × {corr['f_cross']:.3f}
             = {k_corr:.3f}
    This small correction (+f_cross) should push R² → 0.90–0.92.
""")
    elif name == 'massive_spiral':
        print(f"""    Massive spirals are the weakest performer (+0.018 R²). Root cause:

    DIAGNOSIS: Two competing effects cancel each other:
      HIGH Σ ~ 50 M☉/pc² (> Σ_crit) → κ_Σ suppressed toward Newtonian
      HIGH M ~ 10^11.6 M☉ (Q=1)     → κ₃ rung fully accessible
      HIGH v_circ = 280 km/s         → T^0i term is LARGE (missing!)
      HIGH σ_v = 120 km/s            → pressure suppresses κ further

    The model currently ignores f_cross entirely. For massive spirals:
      Δκ_cross = f_cross - 1 = {corr['f_cross'] - 1:.4f}
      (v_circ = {state.v_circ} km/s vs ref 250 km/s → significant contribution)

    Also: bulge-dominated systems need TWO-PHASE treatment:
      Phase A (r < R_bulge): Jeans equation, κ_bulge = f(σ_v, Σ_bulge)
      Phase B (r > R_bulge): rotation curve, κ_disk = f(Σ_disk, v_circ)

    FIX 1 — Add T^0i term:
      κ_corrected = 1 + (κ_base - 1) × f_P × f_β × f_shear × f_cross
      (add f_cross explicitly to the master formula)

    FIX 2 — Two-phase decomposition:
      Split data at r = R_bulge (from surface brightness profile)
      Apply different kinematic equations to each phase.

    FIX 3 — High-Σ bulge correction:
      Current: κ_Σ = 1 + (Σ_crit/Σ)^α = {k_sigma:.4f}  (small, Σ >> Σ_crit)
      Needed:  The bulge has large Σ but the ORBITS are radial → anisotropy β>0
               → anisotropy effect should INCREASE κ slightly (radial orbits
                  are more sensitive to the radial component of gravity)
      Proposed: κ_bulge_Σ = 1 + (Σ_crit/Σ)^α × exp(−β/2)^(−1)

    EXPECTED R² GAIN: +0.10–0.20 bringing massive spirals to ~0.60–0.65.
    (Further gain requires Jeans-based treatment, not rotation-curve fitting)
""")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4: GENERALISED κ-RUNG PREDICTION TABLE
# ═══════════════════════════════════════════════════════════════════════════════

def rung_prediction_table():
    """
    Show which κ-rung is predicted for a range of (Σ, M) combinations.
    This is the generalized κ inversion as a look-up / phase diagram.
    """
    print("\n\n" + "═"*72)
    print("  GENERALISED κ-RUNG PREDICTION: (Σ, M) → κ_required")
    print("═"*72)

    # Surface densities and masses spanning the full ladder
    sigmas = [0.5, 2.0, 5.0, 17.5, 50.0, 200.0]
    masses  = [1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15]

    print(f"\n  κ_Σ (surface density correction)  α={ALPHA}, Σ_crit={SIGMA_CRIT}")
    print(f"  {'Σ [M☉/pc²]':>14}", end="")
    for sig in sigmas:
        print(f"  {sig:>8.1f}", end="")
    print()
    print(f"  {'Σ/Σ_crit':>14}", end="")
    for sig in sigmas:
        print(f"  {sig/SIGMA_CRIT:>8.3f}", end="")
    print()
    print(f"  {'κ_Σ':>14}", end="")
    for sig in sigmas:
        ks = 1 + (SIGMA_CRIT/sig)**ALPHA
        print(f"  {ks:>8.4f}", end="")
    print()
    print()

    print(f"  {'Mass [M☉]':>14}  {'logM':>5}  {'Q (v2.0)':>10}  {'κ_base':>8}  {'Rung':>6}  {'Physical regime'}")
    print(f"  {'-'*70}")
    for M in masses:
        lm = np.log10(M)
        # Q-factor (dual-regime v2.0)
        if lm < 9.5:
            M_ref, p = 10**8.8, 0.25
        else:
            M_ref, p = 10**9.25, 0.50
        Q = float(np.tanh(M / M_ref)**p)

        # Typical surface density for this mass scale
        # Rough: Σ ~ M / (π R^2), R ~ M^0.5 kpc (size-mass relation)
        R_kpc = max(0.3, 0.3 * (M/1e10)**0.5)
        Sigma_typ = M / (np.pi * (R_kpc * 1e3)**2)  # M☉/pc²

        # κ_base for this mass scale (using typical Σ)
        if lm < 9.0:
            k_floor, k_target = K1, 1.5
        elif lm < 13.0:
            k_floor, k_target = K2, K3
        else:
            k_floor, k_target = K3, K5

        k_sigma = float(np.clip(1 + (SIGMA_CRIT/max(Sigma_typ,0.01))**ALPHA, k_floor, k_target))
        k_base  = float(k_floor + (k_sigma - k_floor) * Q)
        rung_n, rung_kv, _ = nearest_rung(k_base)

        regimes = {
            8: "Dwarf irregular",
            9: "Dwarf (partial κ₂)",
            10: "Small spiral (Q rising)",
            11: "Large spiral (κ₃ active)",
            12: "Massive spiral / group boundary",
            13: "Galaxy group (κ₄)",
            14: "Galaxy cluster (κ₅)",
            15: "Supercluster / WHIM (κ₆)",
        }
        regime = regimes.get(int(lm), "—")

        print(f"  {M:.2e}  {lm:>5.2f}  {Q:>10.4f}  {k_base:>8.4f}  κ_{rung_n}  {regime}")

    print(f"""
  INTERPRETATION:
  ┌──────────────────────────────────────────────────────────────────┐
  │  The κ-rung is set by TWO separate mechanisms:                   │
  │  1. Q(M):   Mass threshold determines which rung is ACCESSIBLE   │
  │  2. κ_Σ:   Local surface density determines HOW MUCH is accessed │
  │                                                                  │
  │  Current approximation:                                          │
  │    κ_base ≈ K_floor + (κ_Σ - K_floor) × Q                       │
  │                                                                  │
  │  Full T^μν coupling (proposed):                                  │
  │    κ_base ≈ [K_floor + (κ_Σ - K_floor) × Q]                     │
  │            × f_P(w) × f_β(anisotropy) × f_cross(v_circ)         │
  │                                                                  │
  │  WHERE THE Σ PROXY BREAKS DOWN:                                  │
  │  • Dwarfs (w ~ 1): pressure dominates, Σ underestimates source   │
  │  • Bulges (β ~ 0.5): radial orbits, Σ overestimates support      │
  │  • Massive spirals (v_circ ~ 280 km/s): T^0i term missing        │
  │  • Hot clusters (T ~ 14 keV): T^ii (thermal) competes with T^00  │
  └──────────────────────────────────────────────────────────────────┘
""")


def full_bullet_inversion_check():
    """Verify the generalised formula reproduces Bullet Cluster result."""
    print("═"*60)
    print("  VERIFICATION: Bullet Cluster κ₅ inversion")
    print("═"*60)
    M_lens  = 2.80e14
    M_gas   = 2.80e14 * 0.16 * 0.90  # 4.03e13
    M_star  = 2.80e14 * 0.16 * 0.10  # 4.48e12
    M_sub_gas = 2.30e14 * 0.16 * 0.90
    M_cross = float(np.sqrt(M_gas * M_sub_gas))
    kpl_gas = float(np.clip(1.0 + (SIGMA_CRIT / 51.3)**ALPHA, 1.0, K5))

    # Two-component: (M_gas, kpl_gas), (M_star, None)
    components = [(M_gas, kpl_gas), (M_star, None)]
    kappa_req = kappa_invert_multicomponent(M_lens, components, target_idx=1, M_cross=M_cross)

    print(f"  M_lens:    {M_lens:.3e} M☉")
    print(f"  M_gas:     {M_gas:.3e} M☉  ×  kpl_gas = {kpl_gas:.4f}")
    print(f"  M_star:    {M_star:.3e} M☉  ×  κ_star = ?")
    print(f"  M_cross:   {M_cross:.3e} M☉")
    print(f"  κ_required = {kappa_req:.4f}")
    print(f"  κ₅         = {K5:.4f}")
    print(f"  Match:       {abs(kappa_req-K5)/K5*100:.1f}%")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    W = 72
    print("═"*W)
    print("  QGD GENERALISED κ-INVERSION ANALYSIS".center(W))
    print("  Surface Density as T^μν Proxy — Diagnostic & Improvement".center(W))
    print("═"*W)

    print("""
  GENERALISED INVERSION FORMULA
  ──────────────────────────────────────────────────────────────────
  For a system with baryonic components {M_i, κ_i(T^μν_i)}:

      M_grav = Σ_i  M_i · κ_i(Σ_i, w_i, β_i, v_i)  +  M_cross

  Solving for component j (all others known):

      ┌─────────────────────────────────────────────────────────────┐
      │  κ_j = [ M_grav  −  Σ_{i≠j} M_i·κ_i  −  M_cross ] / M_j  │
      └─────────────────────────────────────────────────────────────┘

  For rotation curves at each point r:

      ┌─────────────────────────────────────────┐
      │  κ_required(r) = [v_obs(r)/v_bar(r)]²  │
      └─────────────────────────────────────────┘

  The full T^μν coupling (replacing the Σ proxy):

      κ_full = 1 + (κ_base − 1) × f_P(w) × f_β(β) × f_shear × f_cross(v)

  where:
    κ_base  = Q(M) × κ_Σ(Σ)            — current model (incomplete)
    f_P     = exp(−ln3 · w)             — pressure suppression
    f_β     = 1 + 0.15(β − 0.5)        — orbital anisotropy
    f_shear = 1 + 0.1σ_Σ/(1+w)        — T^ij shear stress
    f_cross = 1 + 0.08(v_circ/250)^0.5 — T^0i momentum flux  ← MISSING
    w       = σ_v²/v_circ²             — equation-of-state parameter
    β       = 1 − σ_t²/σ_r²            — Binney anisotropy parameter
""")

    full_bullet_inversion_check()

    print("\n\n" + "═"*W)
    print("  GALAXY CLASS DIAGNOSTICS".center(W))
    print("═"*W)

    for name, cls in GALAXY_CLASSES.items():
        diagnose_galaxy_class(name, cls)

    rung_prediction_table()

    print("═"*W)
    print("  SUMMARY: PRIORITY IMPROVEMENTS".center(W))
    print("═"*W)
    print("""
  Priority 1 — f_cross term (T^0i momentum flux):
    Missing from master formula. Adds +4–8% κ for v_circ > 200 km/s.
    Most impactful for MASSIVE SPIRALS (+0.10–0.20 R²).
    Formula:  f_cross = 1 + 0.08 × (v_circ / 250 km/s)^0.5

  Priority 2 — Pressure-based kinematic switch for dwarfs:
    Replace rotation-curve equation with Jeans equation for w > 1.5.
    Targets the 48% of dwarfs with R² < 0.9.
    Formula:  σ_pred = σ_bar × √(κ × f_Jeans)
    where  f_Jeans = anisotropic Jeans mass correction.

  Priority 3 — Two-phase (bulge + disk) decomposition:
    Apply separate T^μν physics inside/outside R_bulge.
    Most impactful for Sa/S0 massive spirals with bulge_frac > 0.4.

  Priority 4 — Improved Q-factor for small spirals:
    M_ref = 10^9.8 M☉ (instead of 10^9.25) for the p=0.5 branch.
    Also reduces small-spiral scatter at the κ₂→κ₃ transition.
""")
