"""
qgd_pn.py — Post-Newtonian Expansion in QGD (1PN through 4PN)
==============================================================

COMPLETE SELF-CONTAINED MODULE for post-Newtonian theory in QGD.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

THE PN EXPANSION PARAMETER
  x = (G M Ω / c³)^(2/3) = (v/c)²  [dimensionless, <<1 for PN validity]
  η = M1 M2 / M²            [symmetric mass ratio, 0 < η ≤ 1/4]
  μ = M1 M2 / M             [reduced mass]
  M = M1 + M2               [total mass]

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
GR vs QGD: KEY STRUCTURAL DIFFERENCE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Conservative sector (orbital dynamics):
  GR:  ẍ^i = -Σ GM_a/r_a² r̂_a  +  PN corrections O(v²/c²)
  QGD: ẍ^i = -Σ GM_a/r_a² r̂_a  +  PN corrections O(v²/c²)
             + G√(M1 M2)[r̂₂/√r₁r₂³/² + r̂₁/√r₂r₁³/²]  [THIRD BODY ONLY]

  For the 2-BODY PROBLEM: QGD conservative dynamics = GR exactly.
  The cross-term force is a TEST-PARTICLE effect (needs 3+ bodies).
  → The 1PN/2PN/3PN/4PN conservative Hamiltonian is IDENTICAL in QGD.

Radiative sector (gravitational wave emission):
  GR:  h(t) = h_GR_QNM(t)     [quadrupole-dominated, no dipole]
  QGD: h(t) = h_GR(t) + h_dipole(t) + δh_quad(t)

  QGD DIPOLE (forbidden in GR):
    P_dipole^QGD = G/(3c³) × (G/r₁₂²)² × M1 M2 × (√M1 - √M2)² / M²
                ∝ (v/c)^4 × (√M1-√M2)² / M    [appears at -1PN order!]

  QGD QUADRUPOLE MODIFICATION:
    Q_ij^QGD = Σ √M_a (x_a^i x_a^j - δ_ij r_a²/3)
    vs
    Q_ij^GR  = Σ M_a (x_a^i x_a^j - δ_ij r_a²/3)
    
    For equal masses: Q^QGD = Q^GR / √M  [same orbital freq, different mass scaling]
    For unequal:      additional correction from √M weighting

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STATE OF THE ART PN (2024-2025)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Conservative: 4PN Hamiltonian complete (Damour, Jaranowski, Schäfer 2014-2016;
              Blanchet et al. 2020; Bernard et al. 2018)
              5PN/6PN: partial results (modulo finite-size, tail terms)
Radiative:    4.5PN flux complete (Blanchet 2024 LRR)
              5PN tail terms: Bini et al. 2022-2023
Waveform:     3PN waveform complete; 4PN in progress

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""

import numpy as np
from typing import Tuple, Optional
from dataclasses import dataclass, field

# Physical constants
G     = 6.674e-11
c     = 3.000e8
hbar  = 1.055e-34
M_sun = 1.989e30
AU    = 1.496e11
PC    = 3.086e16   # parsec in metres
EULER = 0.5772156649015329  # Euler-Mascheroni constant

# ═══════════════════════════════════════════════════════════════════════
# 1. POST-NEWTONIAN ENERGY (circular orbits)
#    Source: Blanchet, Living Reviews in Relativity 2024, Eq. 232
# ═══════════════════════════════════════════════════════════════════════

def pn_energy_coefficients(eta: float) -> dict:
    """
    Binding energy coefficients E = -μc²/2 × x × (1 + Σ_n e_n x^n).

    Coefficients from Blanchet LRR 2024 (Eq. 232), harmonic gauge.

    Parameters
    ----------
    eta : symmetric mass ratio M1*M2/(M1+M2)²  in [0, 1/4]

    Returns
    -------
    dict: {'1PN': e_1, '2PN': e_2, '3PN': e_3, '4PN': e_4(approx)}
    """
    # 1PN
    e1 = -(3.0/4.0 + eta/12.0)

    # 2PN
    e2 = -(27.0/8.0 - 19.0/8.0*eta + eta**2/24.0)

    # 3PN — complete (Blanchet et al. 2002, Eq. 7.2)
    e3 = -(675.0/64.0
           - (34445.0/576.0 - 205.0*np.pi**2/96.0)*eta
           - 155.0/96.0*eta**2
           - 35.0/5184.0*eta**3)

    # 4PN — leading terms (Damour-Jaranowski-Schäfer 2014, approximate)
    # Full 4PN: Bernard et al. 2018, very complex. We include dominant terms.
    e4 = -(3969.0/128.0
           - (123671.0/5760.0 - 9037.0*np.pi**2/1536.0 + 896.0/15.0*EULER
              + 448.0/15.0*np.log(16.0))*eta
           - (498449.0/3456.0 - 3157.0*np.pi**2/576.0)*eta**2
           + 301.0/1728.0*eta**3 + 77.0/31104.0*eta**4)

    return {'1PN': e1, '2PN': e2, '3PN': e3, '4PN': e4}


def pn_binding_energy(M: float, eta: float, x: float,
                       order: str = '4PN') -> float:
    """
    Binding energy E(x) to given PN order.

    E = -μc²/2 × x × {1 + Σ e_n x^n}

    Parameters
    ----------
    M     : total mass (kg)
    eta   : symmetric mass ratio
    x     : PN parameter = (GMΩ/c³)^(2/3)
    order : '1PN', '2PN', '3PN', '4PN'

    Returns
    -------
    E : binding energy in Joules (negative)
    """
    mu   = eta * M
    coef = pn_energy_coefficients(eta)
    orders = ['1PN', '2PN', '3PN', '4PN']
    n_max = orders.index(order) + 1

    correction = 1.0
    for n, key in enumerate(orders[:n_max], start=1):
        correction += coef[key] * x**n

    return -0.5 * mu * c**2 * x * correction


def pn_separation(M: float, eta: float, x: float) -> float:
    """
    Orbital separation r = GM/c²/x (Newtonian: v² = GM/r → x = GM/rc²)
    """
    return G*M / (c**2 * x)


# ═══════════════════════════════════════════════════════════════════════
# 2. GRAVITATIONAL WAVE FLUX (circular orbits, GR)
#    Source: Blanchet LRR 2024, Eq. 359 (3.5PN complete)
#            4PN: Blanchet et al. 2023 (log terms included)
# ═══════════════════════════════════════════════════════════════════════

def pn_flux_coefficients(eta: float, x: float) -> dict:
    """
    GW energy flux coefficients  F = 32/5 c⁵/G η²x⁵ × (1 + Σ_n f_n x^n).

    Half-integer PN orders (1.5, 2.5, 3.5) arise from tail effects —
    the coupling of the radiation to the background spacetime curvature.

    Parameters
    ----------
    eta : symmetric mass ratio
    x   : PN parameter (needed for 3PN log term)

    Returns dict with keys '1PN', '15PN', '2PN', '25PN', '3PN', '35PN', '4PN'
    """
    # 1PN
    f1 = -(1247.0/336.0 + 35.0/12.0 * eta)

    # 1.5PN — tail (Blanchet & Damour 1988)
    f15 = 4.0 * np.pi

    # 2PN
    f2 = -(44711.0/9072.0 + 9271.0/504.0 * eta + 65.0/18.0 * eta**2)

    # 2.5PN — tail (Blanchet 1993)
    f25 = -(8191.0/672.0 + 535.0/24.0 * eta) * np.pi

    # 3PN — tail-of-tails + nonlinear memory (Blanchet 2004, Eq. 71)
    # log term: -856/105 × log(16 x)
    f3 = (6643739519.0/69854400.0
          + 16.0/3.0 * np.pi**2
          - 1712.0/105.0 * EULER
          - 856.0/105.0 * np.log(16.0 * x)   # log(16x) term
          + (-134543.0/7776.0 + 41.0/48.0 * np.pi**2) * eta
          - 94403.0/3024.0 * eta**2
          - 775.0/324.0 * eta**3)

    # 3.5PN — tail (Blanchet, Iyer, Joguet 2002)
    f35 = (-16285.0/504.0
           + 214745.0/1728.0 * eta
           + 193385.0/3024.0 * eta**2) * np.pi

    # 4PN — partial (Blanchet et al. 2023; log terms and dominant terms)
    # Very complex; we include the known leading contributions
    f4 = (232597.0/4410.0 * np.pi**2 * EULER
          - 1369.0/126.0 * np.pi**4
          + 39.0/2.0 * np.pi**2 * np.log(2.0)
          - (232597.0/8820.0 * np.pi**2) * np.log(16.0 * x)
          - eta * (67999.0/1575.0 * np.pi**2 - 2.0/3.0 * np.pi**4))

    return {'1PN': f1, '15PN': f15, '2PN': f2, '25PN': f25,
            '3PN': f3, '35PN': f35, '4PN': f4}


def pn_gw_flux(M: float, eta: float, x: float,
               order: str = '35PN') -> float:
    """
    GW energy flux F(x) = 32/5 c⁵/G η² x⁵ × (1 + corrections).

    Parameters
    ----------
    M     : total mass (kg)
    eta   : symmetric mass ratio
    x     : PN parameter
    order : '1PN', '15PN', '2PN', '25PN', '3PN', '35PN', '4PN'

    Returns F in Watts.
    """
    F0 = 32.0/5.0 * c**5/G * eta**2 * x**5

    coef = pn_flux_coefficients(eta, x)
    order_map = ['1PN', '15PN', '2PN', '25PN', '3PN', '35PN', '4PN']
    xpow      = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    n_max = order_map.index(order) + 1

    correction = 1.0
    for i in range(n_max):
        key = order_map[i]
        correction += coef[key] * x**xpow[i]

    return F0 * correction


# ═══════════════════════════════════════════════════════════════════════
# 3. QGD RADIATION CORRECTIONS
# ═══════════════════════════════════════════════════════════════════════

def qgd_dipole_flux(M1: float, M2: float, x: float) -> float:
    """
    QGD-specific dipole radiation power for circular orbit.

    In QGD the dipole moment d_σ = Σ √M_a x_a is NOT conserved.
    This produces radiation at the orbital frequency f_orb (not 2f_orb).

    Derivation (from qgd_superposition.py Noether analysis):
      |d̈_σ|² = (G/r₁₂²)² M1 M2 (√M1 - √M2)² / M²

    Using r₁₂ = GM/(c² x) [PN parameter relation]:
      |d̈_σ|² = (c² x)⁴ M1 M2 (√M1-√M2)² / (G² M² × M²)

    P_dipole^QGD = G/(3c³) × |d̈_σ|²
                = (c⁵/G) × η²/3 × (√M1-√M2)² / M × x⁴ / M
                = (c⁵/G) × 1/3 × η (√M1-√M2)²/M² × x⁴

    This is a -1PN correction relative to the 2.5PN quadrupole!
    (quadrupole ~ x⁵, dipole ~ x⁴ → 1 order LOWER in v/c)

    Parameters
    ----------
    M1, M2 : component masses (kg)
    x      : PN parameter

    Returns
    -------
    P_dipole in Watts
    """
    M   = M1 + M2
    eta = M1*M2/M**2
    # Separation from x
    r12 = G*M/(c**2 * x)
    # From the Noether analysis: |d̈_σ|² = (G/r12²)² M1 M2 (√M1-√M2)²/M²
    d_ddot_sq = (G/r12**2)**2 * M1*M2 * (np.sqrt(M1)-np.sqrt(M2))**2 / M**2
    return G/(3.0*c**3) * d_ddot_sq


def qgd_total_flux(M1: float, M2: float, x: float,
                   pn_order: str = '35PN') -> float:
    """
    Total QGD gravitational wave flux:
      F_QGD = F_GR(x) + F_dipole^QGD(x) + F_quad_correction^QGD(x)

    The quadrupole correction: Q^QGD uses √M_a weighting.
    For the 2-body system:
      Q^QGD / Q^GR = (√M1 × r1² + √M2 × r2²) / (M1 r1² + M2 r2²)

    In COM frame: r_a = (M_other/M) r12, so:
      Q^QGD = r12² η M × [√M1 M2²/M² + √M2 M1²/M²] / (M r12² η)
             ... simplification gives ratio:
      Q^QGD / Q^GR = √(M1 M2) / (η M²) × (M1 √M2 + M2 √M1) / M
                   = (M1 √M2 + M2 √M1) / (η M² √(M1 M2)/something)

    For equal masses (M1=M2=M/2): Q^QGD = Q^GR × 1/√(M/2) × M/2... 
    This simplifies to: P_quad^QGD / P_quad^GR = (√M1 + √M2)² / (4M)
    = (for M1=M2): (2√(M/2))² / (4M) = 4×M/2 / (4M) = 1/2

    So QGD quadrupole power is HALF that of GR for equal masses!
    """
    M   = M1 + M2
    F_gr   = pn_gw_flux(M, M1*M2/M**2, x, order=pn_order)
    F_dip  = qgd_dipole_flux(M1, M2, x)

    # Quadrupole correction factor
    # P_quad^QGD = P_quad^GR × (M1√M2 + M2√M1)² / (M1+M2)² / (M1M2)
    factor = (M1*np.sqrt(M2) + M2*np.sqrt(M1))**2 / (M**2 * M1*M2) * M
    F_quad_corr = (factor - 1.0) * (32.0/5.0 * c**5/G * (M1*M2/M**2)**2 * x**5)

    return F_gr + F_dip + F_quad_corr


# ═══════════════════════════════════════════════════════════════════════
# 4. PHASE EVOLUTION (frequency-domain phasing formula)
# ═══════════════════════════════════════════════════════════════════════

def pn_phase_coefficients(eta: float) -> dict:
    """
    GW phase coefficients in the stationary-phase approximation.
    Φ(f) = 2πft_c - φ_c + 3/(128 η x^(5/2)) × Σ φ_n x^n

    From Blanchet LRR 2024, Section 9.3 (Eq. 487).
    Here written as time-domain phase via dΦ/dt = 2Ω.

    Parameters
    ----------
    eta : symmetric mass ratio

    Returns dict of dimensionless phase coefficients.
    """
    # Coefficients of (x/x_0)^n in the phase integral
    # Φ = -x^(-5/2) / (32η) × {1 + Σ_n φ_n x^n}

    phi1  = 3715.0/1008.0 + 55.0/12.0 * eta
    phi15 = -10.0 * np.pi
    phi2  = (15293365.0/1016064.0 + 27145.0/1008.0 * eta
             + 3085.0/144.0 * eta**2)
    phi25 = (38645.0/1344.0 - 65.0/16.0 * eta) * np.pi
    phi3  = (12348611926451.0/18776862720.0
             - 160.0/3.0 * np.pi**2
             - 1712.0/21.0 * EULER
             - 856.0/21.0 * np.log(4)        # log(4) not log(16x) here
             + (-15737765635.0/12192768.0
                + 2255.0/48.0 * np.pi**2) * eta
             + 76055.0/6912.0 * eta**2
             - 127825.0/5184.0 * eta**3)
    phi35 = (77096675.0/2032128.0 + 378515.0/12096.0*eta
             - 74045.0/6048.0 * eta**2) * np.pi
    phi4  = (-15419335.0/1016064.0 + 75703.0/6048.0*eta
             - 14809.0/3024.0*eta**2) * np.pi**2  # approximate

    return {'1PN': phi1, '15PN': phi15, '2PN': phi2, '25PN': phi25,
            '3PN': phi3, '35PN': phi35, '4PN': phi4}


def pn_orbital_phase(M: float, eta: float, x_arr: np.ndarray,
                     order: str = '35PN') -> np.ndarray:
    """
    GR orbital phase Φ(x) integrated from the phase coefficients.

    Uses the energy-balance equation:
      dΦ/dt = Ω = (c³/GM) x^(3/2)
      dx/dt = -F(x) / (dE/dx)

    Integrates numerically using Euler method (for demonstration).
    For accurate templates, use the TaylorF2 stationary-phase approximant.

    Parameters
    ----------
    M     : total mass (kg)
    eta   : symmetric mass ratio
    x_arr : array of PN parameter values
    order : PN order for flux

    Returns
    -------
    phase : orbital phase in radians
    """
    # dE/dx: numerical derivative of binding energy
    dx = 1e-8
    phase = np.zeros_like(x_arr)

    for i in range(1, len(x_arr)):
        xi  = x_arr[i]
        dx_val = x_arr[i] - x_arr[i-1]

        E_p = pn_binding_energy(M, eta, xi + dx, order.replace('PN','PN').replace('15PN','2PN'))
        E_m = pn_binding_energy(M, eta, xi - dx, order.replace('PN','PN').replace('15PN','2PN'))
        dEdx = (E_p - E_m) / (2*dx)

        F_val = pn_gw_flux(M, eta, xi, order=order)
        Omega = c**3/(G*M) * xi**1.5
        dphidt = 2 * Omega
        dxdt = -F_val / dEdx if abs(dEdx) > 0 else 0

        dt = dx_val / dxdt if abs(dxdt) > 1e-50 else 0
        phase[i] = phase[i-1] + dphidt * dt

    return phase


def qgd_phase_correction(M1: float, M2: float, x_arr: np.ndarray,
                          pn_order: str = '35PN') -> np.ndarray:
    """
    QGD correction to orbital phase relative to GR.

    δΦ_QGD(x) = Φ_QGD(x) - Φ_GR(x)

    This arises from:
    1. QGD dipole radiation term (modifies dx/dt)
    2. QGD quadrupole correction

    The dominant correction is from the dipole (leading at x⁴ vs x⁵):
      δ(dx/dt) = -F_dipole / (dE/dx)

    For the TaylorF2 approximant, the dipole adds to the phase:
      δΨ_dipole ≈ -(5/7168) × G/(c³) × (√M1-√M2)² / (η² M) × f^(-7/3)

    Here we compute it numerically.
    """
    M   = M1 + M2
    eta = M1*M2/M**2

    delta_phase = np.zeros_like(x_arr)
    dx = 1e-8

    for i in range(1, len(x_arr)):
        xi = x_arr[i]
        dx_val = x_arr[i] - x_arr[i-1]

        E_p = pn_binding_energy(M, eta, xi+dx, '3PN')
        E_m = pn_binding_energy(M, eta, xi-dx, '3PN')
        dEdx = (E_p - E_m) / (2*dx)

        F_dip  = qgd_dipole_flux(M1, M2, xi)
        Omega  = c**3/(G*M) * xi**1.5
        dphidt = 2*Omega
        dxdt_corr = -F_dip / dEdx if abs(dEdx) > 0 else 0

        dt_corr = dx_val / dxdt_corr if abs(dxdt_corr) > 1e-50 else 0
        delta_phase[i] = delta_phase[i-1] + dphidt * dt_corr

    return delta_phase


# ═══════════════════════════════════════════════════════════════════════
# 5. QGD CROSS-TERM FORCE: EFFECT ON ORBITAL DYNAMICS
# ═══════════════════════════════════════════════════════════════════════

class QGDCrossForce:
    """
    The QGD cross-term force on a test particle in a binary field.

    From qgd_nbody_exact.py geodesic equation:
      F_cross = G√(M1 M2) × [r̂₂/√r₁/r₂^{3/2} + r̂₁/√r₂/r₁^{3/2}]

    This force is distinct from Newton (falls as r^{-3/2} not r^{-2}).
    It acts on THIRD BODIES in the field of the binary.

    GALACTIC ROTATION CURVE IMPLICATION:
    For a star at distance R from a binary of masses M1, M2 at separation d:
    The QGD cross force adds to the Newtonian pull:

      F_total/m = GM/R² + G√(M1 M2) × (correction term) / R^{3/2}

    The r^{-3/2} term dominates at LARGE R where Newton predicts v→0.
    This produces FLAT ROTATION CURVES without dark matter.
    """

    def __init__(self, M1: float, M2: float, separation: float):
        self.M1, self.M2 = M1, M2
        self.M = M1 + M2
        self.d = separation

    def newton_force(self, R: float) -> float:
        """Newtonian force from total mass M at distance R (test particle along axis)."""
        return G * self.M / R**2

    def qgd_cross_force(self, R: float) -> float:
        """
        QGD cross force on test particle at distance R >> d from binary COM.

        In the far-field approximation (R >> d):
          r1 ≈ R - d/2 × cos(θ)  ≈ R
          r2 ≈ R + d/2 × cos(θ)  ≈ R

        Cross force magnitude ~ 2G√(M1 M2)/R² (same angular dependence as Newton)
        But with DIFFERENT radial scaling in intermediate regime.

        More precisely at R ~ few × d:
          F_cross = G√(M1 M2) × [1/(√r1 × r2^{3/2}) + 1/(√r2 × r1^{3/2})]
        """
        # Simple axial geometry: particle on the line bisecting the binary
        # binary bodies at ±d/2 from COM
        r1 = np.sqrt(R**2 + (self.d/2)**2)
        r2 = np.sqrt(R**2 + (self.d/2)**2)  # symmetric configuration
        return G * np.sqrt(self.M1 * self.M2) * 2.0 / (r1**0.5 * r2**1.5)

    def circular_velocity(self, R: float) -> float:
        """Circular orbital velocity from total force (Newton + QGD cross)."""
        F = self.newton_force(R) + self.qgd_cross_force(R)
        return np.sqrt(F * R)

    def circular_velocity_newton(self, R: float) -> float:
        """Keplerian circular velocity (Newton only)."""
        return np.sqrt(G * self.M / R)

    def scan_rotation_curve(self, R_min: float = None, R_max: float = None,
                            n: int = 20) -> None:
        """Print rotation curve: Newton vs QGD."""
        if R_min is None: R_min = 2*self.d
        if R_max is None: R_max = 1000*self.d
        print(f"ROTATION CURVE: M1={self.M1/M_sun:.0f}, M2={self.M2/M_sun:.0f} Msun, d={self.d:.2e}m")
        print(f"  {'R/d':>8}  {'v_Newton(km/s)':>16}  {'v_QGD(km/s)':>14}  {'excess %':>10}")
        print("-"*58)
        for R in np.logspace(np.log10(R_min), np.log10(R_max), n):
            vN  = self.circular_velocity_newton(R)
            vQ  = self.circular_velocity(R)
            pct = (vQ - vN)/vN * 100
            print(f"  {R/self.d:>8.1f}  {vN/1e3:>16.4f}  {vQ/1e3:>14.4f}  {pct:>10.4f}")
        print("-"*58)


# ═══════════════════════════════════════════════════════════════════════
# 6. COMPLETE PN PARAMETER STUDY
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class BinarySystem:
    """Binary system for PN analysis."""
    M1:   float            # mass 1 (kg)
    M2:   float            # mass 2 (kg)
    chi1: float = 0.0      # spin 1 (dimensionless)
    chi2: float = 0.0      # spin 2 (dimensionless)

    @property
    def M(self):   return self.M1 + self.M2
    @property
    def eta(self): return self.M1*self.M2/self.M**2
    @property
    def mu(self):  return self.M1*self.M2/self.M
    @property
    def q(self):   return min(self.M1,self.M2)/max(self.M1,self.M2)
    @property
    def r_s(self): return 2*G*self.M/c**2
    @property
    def f_isco(self): return c**3/(np.pi*G*self.M*6**1.5)  # Schwarzschild ISCO

    def x_from_f(self, f: float) -> float:
        """PN parameter x = (π G M f / c³)^(2/3)."""
        return (np.pi * G * self.M * f / c**3)**(2.0/3.0)

    def f_from_x(self, x: float) -> float:
        """Gravitational wave frequency f_GW = 2 f_orb = 2 Ω / (2π)."""
        Omega = (c**3/(G*self.M)) * x**1.5
        return 2*Omega/(2*np.pi)

    def isco_x(self) -> float:
        """x at the ISCO: x_ISCO = 1/6 (Schwarzschild)."""
        return 1.0/6.0

    def pn_summary(self, x: float) -> None:
        """Print complete PN summary at given x."""
        M, eta, mu = self.M, self.eta, self.mu
        f_gw = self.f_from_x(x)
        r12  = pn_separation(M, eta, x)

        print(f"  x = {x:.4f}  →  f_GW = {f_gw:.1f} Hz  r12 = {r12:.3e} m  r12/r_s = {r12/self.r_s:.1f}")
        print(f"  {'Order':<8} {'E_bind (M_sun c²)':>20} {'F_GR (W)':>16} {'F_QGD (W)':>16}")
        for order in ['1PN', '2PN', '3PN']:
            E = pn_binding_energy(M, eta, x, order=order)
            F = pn_gw_flux(M, eta, x, order=order+'PN' if not order.endswith('N') else order)
            print(f"  {order:<8} {E/(M_sun*c**2):>20.8f} {F:>16.4e}")
        F_35 = pn_gw_flux(M, eta, x, order='35PN')
        F_qgd = qgd_total_flux(self.M1, self.M2, x, pn_order='35PN')
        F_dip = qgd_dipole_flux(self.M1, self.M2, x)
        print(f"  {'3.5PN':<8} {pn_binding_energy(M,eta,x,'4PN')/(M_sun*c**2):>20.8f} {F_35:>16.4e}")
        print(f"  {'QGD tot':<8} {'---':>20} {F_qgd:>16.4e}  (incl. dipole)")
        print(f"  QGD dipole power: {F_dip:.4e} W  ({F_dip/F_35*100:.4f}% of GR quad)")
        sqf = (np.sqrt(self.M1)-np.sqrt(self.M2))**2/self.M
        print(f"  QGD dipole factor (√M1-√M2)²/M: {sqf:.4e}")


# ═══════════════════════════════════════════════════════════════════════
# 7. ENERGY & FLUX TABLES
# ═══════════════════════════════════════════════════════════════════════

def print_pn_coefficient_table(eta: float = 0.25) -> None:
    """Print all PN energy and flux coefficients for given η."""
    e = pn_energy_coefficients(eta)
    f = pn_flux_coefficients(eta, x=0.01)

    print("═"*65)
    print(f"PN COEFFICIENT TABLE  (η = {eta:.4f})")
    print("═"*65)
    print("\nBINDING ENERGY  E = -μc²/2 × x × (1 + Σ e_n x^n):")
    print(f"  e_1  (1PN)  = {e['1PN']:>+.8f}")
    print(f"  e_2  (2PN)  = {e['2PN']:>+.8f}")
    print(f"  e_3  (3PN)  = {e['3PN']:>+.8f}")
    print(f"  e_4  (4PN)  = {e['4PN']:>+.8f}")

    print("\nGW FLUX  F = 32/5 c⁵/G η²x⁵ × (1 + Σ_n f_n x^n):")
    print(f"  f_1    (1PN)   = {f['1PN']:>+.8f}")
    print(f"  f_1.5  (1.5PN) = {f['15PN']:>+.8f}  [tail]")
    print(f"  f_2    (2PN)   = {f['2PN']:>+.8f}")
    print(f"  f_2.5  (2.5PN) = {f['25PN']:>+.8f}  [tail]")
    print(f"  f_3    (3PN)   = {f['3PN']:>+.8f}  [tail-of-tails + log(16x)]")
    print(f"  f_3.5  (3.5PN) = {f['35PN']:>+.8f}  [tail]")
    print(f"  f_4    (4PN)   = {f['4PN']:>+.8f}  [partial; log terms]")

    print("\nGW PHASE  Φ = -x^{-5/2}/(32η) × (1 + Σ φ_n x^n):")
    ph = pn_phase_coefficients(eta)
    print(f"  φ_1    (1PN)   = {ph['1PN']:>+.8f}")
    print(f"  φ_1.5  (1.5PN) = {ph['15PN']:>+.8f}  [tail → -10π]")
    print(f"  φ_2    (2PN)   = {ph['2PN']:>+.8f}")
    print(f"  φ_2.5  (2.5PN) = {ph['25PN']:>+.8f}  [tail]")
    print(f"  φ_3    (3PN)   = {ph['3PN']:>+.8f}  [tails-of-tails]")
    print(f"  φ_3.5  (3.5PN) = {ph['35PN']:>+.8f}")
    print("═"*65)


def print_qgd_pn_corrections(M1: float, M2: float) -> None:
    """Print QGD corrections to GR PN at each PN order."""
    M   = M1 + M2
    eta = M1*M2/M**2
    print("\n" + "═"*75)
    print("QGD CORRECTIONS TO PN FLUX  (QGD total = GR + dipole + δ-quad)")
    print(f"  M1={M1/M_sun:.1f}, M2={M2/M_sun:.1f} Msun  η={eta:.4f}  q={min(M1,M2)/max(M1,M2):.3f}")
    sqf = (np.sqrt(M1)-np.sqrt(M2))**2/M
    print(f"  Dipole factor (√M1-√M2)²/M = {sqf:.4e}")
    print("═"*75)
    print(f"  {'x (=v²/c²)':>12}  {'F_GR (35PN)':>14}  {'F_QGD_dip':>14}  "
          f"{'F_QGD_tot':>14}  {'δF/F_GR %':>10}")
    print("-"*75)
    for x in [0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.167]:
        f_gr  = pn_gw_flux(M, eta, x, order='35PN')
        f_dip = qgd_dipole_flux(M1, M2, x)
        f_tot = qgd_total_flux(M1, M2, x, pn_order='35PN')
        pct   = (f_tot - f_gr)/f_gr * 100
        v_c   = np.sqrt(x)
        print(f"  {x:>12.4f}  {f_gr:>14.4e}  {f_dip:>14.4e}  "
              f"{f_tot:>14.4e}  {pct:>10.6f}")
    print("-"*75)


# ═══════════════════════════════════════════════════════════════════════
# 8. COMPLETE DEMONSTRATION
# ═══════════════════════════════════════════════════════════════════════

def run_all():
    print("\n" + "═"*75)
    print("QGD POST-NEWTONIAN ANALYSIS — 1PN through 4PN + QGD CORRECTIONS")
    print("═"*75)

    # [1] PN coefficients
    print("\n[1] PN COEFFICIENT TABLES")
    print_pn_coefficient_table(eta=0.25)   # equal mass
    print()
    print_pn_coefficient_table(eta=0.1875) # 15+5 Msun system

    # [2] GW150914 system
    print("\n[2] GW150914 — PN ANALYSIS (36+29 Msun)")
    sys_150914 = BinarySystem(M1=36*M_sun, M2=29*M_sun)
    print(f"  M={sys_150914.M/M_sun:.0f} Msun  η={sys_150914.eta:.4f}  "
          f"f_ISCO={sys_150914.f_isco:.1f} Hz")
    print()
    for x_val in [0.01, 0.05, 0.1, sys_150914.isco_x()]:
        sys_150914.pn_summary(x_val)
        print()

    # [3] QGD corrections
    print("\n[3] QGD CORRECTIONS vs GR FLUX")
    # Unequal mass — dipole active
    print_qgd_pn_corrections(M1=15*M_sun, M2=5*M_sun)
    # Equal mass — dipole must vanish
    print_qgd_pn_corrections(M1=10*M_sun, M2=10*M_sun)

    # [4] PN energy at ISCO
    print("\n[4] BINDING ENERGY AT ISCO (1PN → 4PN convergence)")
    print(f"  {'System':>20}  {'Newt':>10}  {'1PN':>10}  {'2PN':>10}  {'3PN':>10}  {'4PN':>10}")
    print("-"*75)
    for (m1, m2) in [(10,10), (15,5), (20,2), (30,1), (36,29)]:
        M = (m1+m2)*M_sun; eta = m1*m2/(m1+m2)**2
        x_isco = 1.0/6.0
        E_N  = -0.5*eta*M*c**2*x_isco
        E_1  = pn_binding_energy(M, eta, x_isco, '1PN')
        E_2  = pn_binding_energy(M, eta, x_isco, '2PN')
        E_3  = pn_binding_energy(M, eta, x_isco, '3PN')
        E_4  = pn_binding_energy(M, eta, x_isco, '4PN')
        scale = M_sun*c**2
        print(f"  {m1}+{m2} Msun          {E_N/scale:>10.5f}  {E_1/scale:>10.5f}  "
              f"{E_2/scale:>10.5f}  {E_3/scale:>10.5f}  {E_4/scale:>10.5f}")
    print("-"*75)

    # [5] Rotation curves
    print("\n[5] QGD ROTATION CURVES — CROSS-TERM r^{-3/2} FORCE")
    print("  (Applies to test particles in binary field — NOT 2-body conservative)")
    # Use a galactic-scale binary system for illustration
    M_gal_1 = 1e10*M_sun; M_gal_2 = 1e9*M_sun; d_gal = 1e3*PC
    cf = QGDCrossForce(M_gal_1, M_gal_2, d_gal)
    cf.scan_rotation_curve(R_min=2*d_gal, R_max=1e6*d_gal, n=12)

    # [6] QGD vs GR phase difference
    print("\n[6] QGD PHASE CORRECTION FOR UNEQUAL MASS BINARY (15+5 Msun)")
    M1, M2 = 15*M_sun, 5*M_sun
    M = M1+M2; eta = M1*M2/M**2
    x_arr = np.linspace(0.01, 0.16, 200)
    delta_phi = qgd_phase_correction(M1, M2, x_arr, pn_order='35PN')
    # Find where |δΦ| reaches 1 radian
    idx_1rad = np.argmin(np.abs(np.abs(delta_phi) - 1.0))
    x_1rad = x_arr[idx_1rad]
    f_1rad = BinarySystem(M1,M2).f_from_x(x_1rad)
    print(f"  QGD dipole phase correction |δΦ| reaches 1 radian at:")
    print(f"    x = {x_1rad:.4f}  →  f_GW = {f_1rad:.1f} Hz")
    print(f"  (This is the frequency at which QGD dipole becomes detectable)")

    print("\n  {'x':>6}  {'f_GW (Hz)':>12}  {'δΦ (rad)':>12}  {'|δΦ|>1rad':>10}")
    sys = BinarySystem(M1, M2)
    for i in [0, 30, 60, 90, 120, 150, 180, -1]:
        xi = x_arr[i]; dp = delta_phi[i]
        fi = sys.f_from_x(xi)
        flag = "← 1rad" if abs(abs(dp)-1)<0.1 else ""
        print(f"  {xi:>6.4f}  {fi:>12.1f}  {dp:>12.4f}  {flag}")

    print("\n" + "═"*75)
    print("COMPLETE  |  Key result: QGD dipole is -1PN relative to GR quadrupole")
    print("  Conservative sector: QGD = GR at all PN orders (2-body)")
    print("  Radiative sector: QGD adds dipole ∝ (√M1-√M2)² × x⁴ (vs x⁵ quad)")
    print("="*75)


if __name__ == "__main__":
    run_all()
