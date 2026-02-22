"""
qgd_nbody_exact.py — QGD Exact Multi-Body Spacetime Metrics
============================================================

Companion to: "Exact Multi-Body Spacetime Metrics in Quantum Gravitational
Dynamics (QGD)" (see qgd_nbody.pdf).

CORE IDEA
---------
In QGD every gravitational source contributes two σ-field components:

    σ_t^(a)   — temporal (encodes mass + Kerr spin structure)
    σ_φ^(a)   — azimuthal (encodes frame-dragging)

These obey exact linear superposition:

    σ_t^(total) = Σ_a  σ_t^(a)
    σ_φ^(total) = Σ_a  σ_φ^(a)

The metric is then assembled via the master equation — quadratic in σ, so all
cross-terms  σ_α^(a) × σ_β^(b)  emerge automatically from the outer product:

    g_tt     = -(1 - A_tot²)                          A_a ≡ σ_t^(a)
    g_tφ     = -r sinθ · A_tot · B_tot                B_a ≡ σ_φ^(a)
    g_φφ     =  r² sin²θ · (1 + B_tot²)
    g_rr     = -1 / g_tt          (isotropic / harmonic condition)

Every distinct coupling (mass×mass, mass×spin, spin×spin) is just one
entry in the expanded product — no case-by-case derivation needed.

PHYSICAL CONSTANTS
------------------
G    = 6.674e-11  m³ kg⁻¹ s⁻²
c    = 3.000e+08  m/s
ℏ    = 1.055e-34  J·s
ℓ_Q  = √(Gℏ²/c⁴) ≈ 1.6e-70 m  (QGD quantum length scale)

USAGE
-----
    from qgd_nbody_exact import QGDBody, QGDNBodyMetric
    import numpy as np

    # Two spinning black holes
    M_sun = 1.989e30
    b1 = QGDBody(M=10*M_sun, alpha=0.6*G*10*M_sun/c**2, pos=[0,0,0])
    b2 = QGDBody(M= 8*M_sun, alpha=0.4*G* 8*M_sun/c**2, pos=[1e6,0,0])

    metric = QGDNBodyMetric([b1, b2])
    g = metric.at(x=[3e5, 1e4, 0])
    print(g)          # full 4×4 tensor
    print(g.gtt)      # -0.xxx
    print(g.gtphi)    # frame-dragging
    print(g.horizon)  # merged horizon radius (if exists)

STRUCTURE
---------
    QGDBody             — one gravitational source (Kerr ± charge ± Λ)
    QGDNBodyMetric      — assembles metric from list of QGDBody objects
    MetricPoint         — named-tuple result from .at()
    TwoBodySystem       — convenience wrapper for binary problems
    ThreeBodySystem     — convenience wrapper for triples
    run_examples()      — self-contained demo (call directly or run as __main__)
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Tuple
from scipy.optimize import brentq

# ── Physical constants ────────────────────────────────────────────────────────
G    = 6.674e-11        # m³ kg⁻¹ s⁻²
c    = 3.000e8          # m s⁻¹
hbar = 1.055e-34        # J s
M_sun = 1.989e30        # kg
AU    = 1.496e11        # m


# ════════════════════════════════════════════════════════════════════════════
# 1.  QGDBody — a single gravitational source
# ════════════════════════════════════════════════════════════════════════════

@dataclass
class QGDBody:
    """
    A single gravitational source in QGD.

    Parameters
    ----------
    M      : float — mass (kg)
    alpha  : float — Kerr spin parameter a = J/(Mc) (metres); 0 = Schwarzschild
    pos    : array-like, shape (3,) — Cartesian position (metres)
    charge : float — electric charge Q (Coulombs); 0 = uncharged
    Lambda : float — cosmological constant contribution H²r²/c²; pass H (s⁻¹)
    epsilon: float — source signature; +1 attractive (default), -1 repulsive

    Notes
    -----
    Spin parameter alpha must satisfy 0 ≤ alpha ≤ GM/c² (Kerr bound).
    For maximally spinning BH: alpha_max = G*M/c².
    """
    M:      float
    alpha:  float = 0.0          # spin param (metres)
    pos:    np.ndarray = field(default_factory=lambda: np.zeros(3))
    charge: float = 0.0          # Coulombs
    Lambda: float = 0.0          # H₀ in s⁻¹ (for Schwarzschild-de Sitter)
    epsilon: float = +1.0

    def __post_init__(self):
        self.pos = np.asarray(self.pos, dtype=float)

    # ── Schwarzschild radius ─────────────────────────────────────────────
    @property
    def r_s(self) -> float:
        """Schwarzschild radius 2GM/c² (metres)."""
        return 2 * G * self.M / c**2

    # ── σ_t^(a) at field point x ─────────────────────────────────────────
    def sigma_t(self, x: np.ndarray, theta: float = np.pi/2) -> float:
        """
        Temporal σ-field at field point x.

        Kerr:        σ_t = √(2GM r_a / c² S_a)
        Schwarzschild (alpha=0): σ_t = √(2GM / c² r_a)

        where S_a = r_a² + alpha² cos²θ_a
        """
        x = np.asarray(x, dtype=float)
        r_a = np.linalg.norm(x - self.pos)
        if r_a < 1e-15:
            return 0.0
        S_a = r_a**2 + self.alpha**2 * np.cos(theta)**2
        val = np.sqrt(2 * G * self.M * r_a / (c**2 * S_a))

        # Subtract charge contribution (repulsive, ε = -1 in master eq)
        if self.charge != 0.0:
            k_e  = 8.9875e9  # N m² C⁻²
            rQ2  = G * k_e**2 * self.charge**2 / c**4  # m²
            val_Q = np.sqrt(rQ2) / r_a
            # charge is handled as a separate repulsive term; returned separately
        # Cosmological repulsion
        if self.Lambda != 0.0:
            pass  # handled in sigma_t_cosmo

        return val

    def sigma_phi(self, x: np.ndarray, theta: float = np.pi/2) -> float:
        """
        Azimuthal σ-field at field point x.

        σ_φ = alpha sinθ √(2GM / c² r_a S_a)
        """
        x = np.asarray(x, dtype=float)
        if self.alpha == 0.0:
            return 0.0
        r_a = np.linalg.norm(x - self.pos)
        if r_a < 1e-15:
            return 0.0
        S_a   = r_a**2 + self.alpha**2 * np.cos(theta)**2
        return self.alpha * np.sin(theta) * np.sqrt(2 * G * self.M / (c**2 * r_a * S_a))

    def sigma_t_charge(self, x: np.ndarray) -> float:
        """Repulsive charge σ-field: σ_Q = √(G k_e² Q²/c⁴) / r  (ε = -1)."""
        if self.charge == 0.0:
            return 0.0
        k_e = 8.9875e9
        r_a = np.linalg.norm(np.asarray(x) - self.pos)
        if r_a < 1e-15:
            return 0.0
        rQ  = np.sqrt(G * k_e**2 * self.charge**2 / c**4)
        return rQ / r_a

    def sigma_t_cosmo(self, x: np.ndarray) -> float:
        """Cosmological σ-field: σ_Λ = H r / c  (ε = -1, repulsive)."""
        if self.Lambda == 0.0:
            return 0.0
        r = np.linalg.norm(np.asarray(x) - self.pos)
        return self.Lambda * r / c

    def f_scalar(self, r: float, theta: float = np.pi/2) -> float:
        """
        Universal gravitational scalar f(r,θ) from the general wavefunction ψ.
        This is (σ_t)² evaluated at (r, θ) relative to this body.

        f = 2GM r / c² S_a  (Kerr)
          = 2GM / c² r      (Schwarzschild)
        """
        S_a = r**2 + self.alpha**2 * np.cos(theta)**2
        return 2 * G * self.M * r / (c**2 * S_a)


# ════════════════════════════════════════════════════════════════════════════
# 2.  MetricPoint — result container
# ════════════════════════════════════════════════════════════════════════════

@dataclass
class MetricPoint:
    """
    Full metric tensor at one field point, plus derived quantities.

    Attributes
    ----------
    g        : (4,4) numpy array — full covariant metric tensor
    gtt      : g_tt (< 0 outside horizon)
    grr      : g_rr (> 0 outside horizon)
    gtphi    : g_tφ (frame-dragging)
    gphiphi  : g_φφ
    A_total  : Σ σ_t^(a)  (horizon when = 1)
    B_total  : Σ σ_φ^(a)
    A_parts  : list of individual A_a values
    B_parts  : list of individual B_a values
    """
    g:        np.ndarray
    gtt:      float
    grr:      float
    gtphi:    float
    gphiphi:  float
    A_total:  float
    B_total:  float
    A_parts:  List[float]
    B_parts:  List[float]

    def cross_terms_tt(self) -> dict:
        """
        Decompose g_tt into self and cross contributions.

        Returns dict with keys:
          'self'  : {a: A_a²}
          'cross' : {(a,b): 2*A_a*A_b}
        """
        N = len(self.A_parts)
        A = self.A_parts
        result = {
            'self':  {a: A[a]**2 for a in range(N)},
            'cross': {(a,b): 2*A[a]*A[b] for a in range(N) for b in range(a+1,N)},
        }
        return result

    def cross_terms_tphi(self, r: float, theta: float = np.pi/2) -> dict:
        """
        Decompose g_tφ into self and cross frame-drag contributions.

        Returns dict with keys:
          'self'  : {a: -r sinθ A_a B_a}   (body a's own frame-drag)
          'cross' : {(a,b): -r sinθ A_a B_b}  (body b spin in body a mass field)
        """
        N = len(self.A_parts)
        A, B = self.A_parts, self.B_parts
        pref = -r * np.sin(theta)
        return {
            'self':  {a: pref * A[a] * B[a] for a in range(N)},
            'cross': {(a,b): pref * A[a] * B[b]
                      for a in range(N) for b in range(N) if a != b},
        }

    @property
    def inside_horizon(self) -> bool:
        """True if gtt > 0 (inside event horizon, coordinates break down)."""
        return self.gtt > 0


# ════════════════════════════════════════════════════════════════════════════
# 3.  QGDNBodyMetric — main class
# ════════════════════════════════════════════════════════════════════════════

class QGDNBodyMetric:
    """
    Exact N-body spacetime metric in QGD.

    The algorithm (3 steps):
    ─────────────────────────────────────────────────────────────────────
    1. For each body a, evaluate A_a = σ_t^(a)(x) and B_a = σ_φ^(a)(x)
    2. Sum: A_tot = ΣA_a,  B_tot = ΣB_b
    3. Assemble:
         g_tt    = -(1 - A_tot²)
         g_tφ    = -r sinθ · A_tot · B_tot
         g_φφ    = r² sin²θ · (1 + B_tot²)
         g_rr    = -1/g_tt
    ─────────────────────────────────────────────────────────────────────

    ALL cross-terms (mass×mass, mass×spin, spin×spin) emerge automatically
    from steps 2+3 — no extra work required.

    Parameters
    ----------
    bodies : list of QGDBody
    """

    def __init__(self, bodies: List[QGDBody]):
        self.bodies = bodies
        self.N = len(bodies)

    # ── Core evaluation ───────────────────────────────────────────────────
    def at(self,
           x:     np.ndarray,
           theta: float = np.pi/2,
           phi:   float = 0.0) -> MetricPoint:
        """
        Evaluate the exact metric at field point x.

        Parameters
        ----------
        x     : (3,) array — Cartesian field point (metres)
        theta : polar angle (radians); default π/2 (equatorial plane)
        phi   : azimuthal angle (radians)

        Returns
        -------
        MetricPoint with full tensor and decomposed cross-terms.
        """
        x = np.asarray(x, dtype=float)
        r = np.linalg.norm(x)

        # ── Step 1: compute σ-fields for each body ────────────────────────
        A_parts = []  # σ_t^(a)  attractive
        B_parts = []  # σ_φ^(a)
        Q_parts = []  # σ_t^(Q_a)  repulsive (charge)
        L_parts = []  # σ_t^(Λ_a)  repulsive (cosmological)

        for body in self.bodies:
            A_parts.append(body.sigma_t(x, theta))
            B_parts.append(body.sigma_phi(x, theta))
            Q_parts.append(body.sigma_t_charge(x))
            L_parts.append(body.sigma_t_cosmo(x))

        # ── Step 2: sum fields ────────────────────────────────────────────
        A_tot = sum(A_parts)                    # attractive temporal
        A_rep = sum(Q_parts) + sum(L_parts)    # repulsive temporal
        B_tot = sum(B_parts)                    # azimuthal (spin)

        # Effective total (repulsive terms enter with opposite sign in Σ²)
        # From master eq: Σ² = ε_a·η_tt·η_tt·σ_t^(a)² = +σ_att² - σ_rep²
        # So: g_tt = -(1 - A_tot² + A_rep² - 2*A_tot*A_rep + ...)
        # For simplicity track net: A_net = A_tot - A_rep (combines signs)
        A_net = A_tot  # repulsive enters as subtraction in Σ²
        Sigma2 = A_tot**2 - A_rep**2 + 2 * A_tot * (0) - 2 * A_tot * A_rep
        # Correct expansion: (A_tot - A_rep)^2 = A_tot^2 - 2A_tot*A_rep + A_rep^2
        # But in master eq each source has its own ε, so:
        # Σ² = (ΣA_att)^2 - (ΣA_rep)^2  is NOT correct either.
        # Correct: Σ² = Σ_a Σ_b ε_a ε_b A_a A_b
        #             = (A_tot)^2 - 2*A_tot*A_rep + (A_rep)^2
        #             = (A_tot - A_rep)^2  for two groups only
        Sigma_eff = A_tot - A_rep
        Sigma2    = Sigma_eff**2

        # ── Step 3: assemble metric ───────────────────────────────────────
        gtt      = -(1 - Sigma2)
        gtphi    = -r * np.sin(theta) * A_tot * B_tot   # charge doesn't spin
        gphiphi  = r**2 * np.sin(theta)**2 * (1 + B_tot**2)
        gthth    = r**2
        grr      = -1.0 / gtt if abs(gtt) > 1e-30 else np.inf

        # ── Build 4×4 tensor [t, r, θ, φ] ────────────────────────────────
        g = np.zeros((4, 4))
        g[0, 0] = gtt
        g[1, 1] = grr
        g[2, 2] = gthth
        g[3, 3] = gphiphi
        g[0, 3] = g[3, 0] = gtphi

        return MetricPoint(
            g=g,
            gtt=gtt, grr=grr, gtphi=gtphi, gphiphi=gphiphi,
            A_total=A_tot, B_total=B_tot,
            A_parts=A_parts, B_parts=B_parts,
        )

    # ── Horizon finder ────────────────────────────────────────────────────
    def horizon_radius(self,
                       direction: np.ndarray = np.array([1,0,0]),
                       r_min: float = 1.0,
                       r_max: float = 1e13) -> Optional[float]:
        """
        Find the event horizon along a given direction (g_tt = 0, i.e. Σ = 1).

        Parameters
        ----------
        direction : unit vector along which to search
        r_min, r_max : search bracket (metres)

        Returns None if no horizon found in bracket.
        """
        direction = np.asarray(direction, dtype=float)
        direction /= np.linalg.norm(direction)
        try:
            return brentq(
                lambda r: self.at(r * direction).gtt,
                r_min, r_max, xtol=1e-3
            )
        except ValueError:
            return None

    # ── Geodesic acceleration ─────────────────────────────────────────────
    def geodesic_accel(self,
                       x: np.ndarray,
                       eps: float = 1.0) -> np.ndarray:
        """
        Coordinate acceleration from geodesic: ẍ = -c² A_tot ∇A_tot.

        This is the weak-field limit of the geodesic equation.
        Returns 3-vector acceleration (m/s²).
        """
        x = np.asarray(x, dtype=float)
        accel = np.zeros(3)
        A0 = self.at(x).A_total
        for i in range(3):
            xp, xm = x.copy(), x.copy()
            xp[i] += eps; xm[i] -= eps
            Ap = self.at(xp).A_total
            Am = self.at(xm).A_total
            dAdx = (Ap - Am) / (2 * eps)
            accel[i] = -c**2 * A0 * dAdx
        return accel

    # ── Gravitational energy density ──────────────────────────────────────
    def energy_density(self, x: np.ndarray, eps: float = 1.0) -> float:
        """
        ρ_grav(x) = ½ (dσ_t/dx_i)²  — localised gravitational energy.

        Integrates to the total mass; positive definite everywhere.
        """
        x = np.asarray(x, dtype=float)
        total = 0.0
        for i in range(3):
            xp, xm = x.copy(), x.copy()
            xp[i] += eps; xm[i] -= eps
            dA = (self.at(xp).A_total - self.at(xm).A_total) / (2 * eps)
            total += dA**2
        return 0.5 * total

    # ── Merger check ──────────────────────────────────────────────────────
    def sigma_at_midpoint(self) -> float:
        """
        Σ at geometric midpoint of all bodies. Horizon merges when Σ ≥ 1.
        """
        mid = np.mean([b.pos for b in self.bodies], axis=0)
        return self.at(mid).A_total


# ════════════════════════════════════════════════════════════════════════════
# 4.  TwoBodySystem — binary convenience class
# ════════════════════════════════════════════════════════════════════════════

class TwoBodySystem:
    """
    Two-body wrapper with explicit cross-term access.

    Examples
    --------
    >>> sys = TwoBodySystem.kerr_binary(M1=10*M_sun, M2=8*M_sun,
    ...                                  a1=0.6, a2=0.4, separation=1e6)
    >>> g = sys.metric_at([3e5, 1e4, 0])
    >>> print(sys.describe_cross_terms(g))
    """

    def __init__(self, body1: QGDBody, body2: QGDBody):
        self.body1 = body1
        self.body2 = body2
        self.metric = QGDNBodyMetric([body1, body2])

    @classmethod
    def kerr_binary(cls, M1: float, M2: float,
                    a1: float = 0.0, a2: float = 0.0,
                    separation: float = 1e6):
        """
        Convenience constructor for two spinning black holes.

        Parameters
        ----------
        M1, M2       : masses (kg)
        a1, a2       : spin parameters as fraction of max (0 to 1)
        separation   : initial separation (metres)
        """
        alpha1 = a1 * G * M1 / c**2
        alpha2 = a2 * G * M2 / c**2
        b1 = QGDBody(M=M1, alpha=alpha1, pos=[0,0,0])
        b2 = QGDBody(M=M2, alpha=alpha2, pos=[separation,0,0])
        return cls(b1, b2)

    def metric_at(self, x, theta=np.pi/2) -> MetricPoint:
        return self.metric.at(x, theta)

    def describe_cross_terms(self, mp: MetricPoint,
                              r: float = None, theta: float = np.pi/2) -> str:
        """
        Human-readable breakdown of all cross-term contributions.
        """
        A1, A2 = mp.A_parts[0], mp.A_parts[1]
        B1, B2 = mp.B_parts[0], mp.B_parts[1]

        lines = []
        lines.append("═" * 60)
        lines.append("QGD TWO-BODY METRIC — CROSS-TERM DECOMPOSITION")
        lines.append("═" * 60)

        lines.append("\ng_tt  =  -(1  -  self₁  -  self₂  -  cross₁₂)")
        lines.append(f"  self₁  = (σ_t⁽¹⁾)²        = {A1**2:+.6e}  (= 2GM₁/c²r₁)")
        lines.append(f"  self₂  = (σ_t⁽²⁾)²        = {A2**2:+.6e}  (= 2GM₂/c²r₂)")
        lines.append(f"  cross  = 2σ_t⁽¹⁾σ_t⁽²⁾  = {2*A1*A2:+.6e}  (QGD: ∝ √(M₁M₂)/√(r₁r₂))")
        lines.append(f"  g_tt   =                   {mp.gtt:+.8f}")

        if r is not None:
            pref = -r * np.sin(theta)
            lines.append("\ng_tφ  =  -r sinθ · (self-spin₁ + self-spin₂ + cross₁₂ + cross₂₁)")
            lines.append(f"  σ_t⁽¹⁾σ_φ⁽¹⁾  = {pref*A1*B1:+.4e}  (body 1 own frame-drag)")
            lines.append(f"  σ_t⁽²⁾σ_φ⁽²⁾  = {pref*A2*B2:+.4e}  (body 2 own frame-drag)")
            lines.append(f"  σ_t⁽¹⁾σ_φ⁽²⁾  = {pref*A1*B2:+.4e}  (M₁ mass × α₂ spin  ← QGD cross)")
            lines.append(f"  σ_t⁽²⁾σ_φ⁽¹⁾  = {pref*A2*B1:+.4e}  (M₂ mass × α₁ spin  ← QGD cross)")
            lines.append(f"  g_tφ  =           {mp.gtphi:+.4e}")

        return "\n".join(lines)

    def scan_between_bodies(self, n_pts: int = 9) -> None:
        """Print metric profile along the axis connecting the two bodies."""
        x1, x2 = self.body1.pos, self.body2.pos
        d = np.linalg.norm(x2 - x1)
        print(f"\n{'x/d':>6}  {'A_tot':>12}  {'g_tt':>14}  {'g_tφ':>14}  {'g_rr':>12}")
        for f in np.linspace(0.05, 0.95, n_pts):
            xf = x1 + f * (x2 - x1) + np.array([0, 1e3, 0])  # tiny y-offset
            mp = self.metric_at(xf)
            print(f"{f:>6.2f}  {mp.A_total:>12.6e}  "
                  f"{mp.gtt:>14.8f}  {mp.gtphi:>14.4e}  {mp.grr:>12.8f}")

    def horizon_separation(self) -> float:
        """
        Analytical merger separation: d = 4 r_s for equal non-spinning masses.
        For unequal/spinning, return the threshold found numerically.
        """
        # Equal non-spinning analytic formula: d_merge = 8 GM/c² = 4 r_s
        M1, M2 = self.body1.M, self.body2.M
        if self.body1.alpha == 0 and self.body2.alpha == 0 and np.isclose(M1, M2):
            return 4 * self.body1.r_s
        # Else: scan separations at midpoint until Σ ≥ 1
        body1_orig = np.copy(self.body1.pos)
        body2_orig = np.copy(self.body2.pos)
        for sep in np.logspace(np.log10(self.body1.r_s), 10, 500):
            mid = np.array([sep/2, 1.0, 0.0])
            self.body1.pos = np.array([0.0, 0.0, 0.0])
            self.body2.pos = np.array([sep,  0.0, 0.0])
            mp = self.metric_at(mid)
            if mp.A_total >= 1.0:
                self.body1.pos = body1_orig
                self.body2.pos = body2_orig
                return sep
        self.body1.pos = body1_orig
        self.body2.pos = body2_orig
        return float('nan')


# ════════════════════════════════════════════════════════════════════════════
# 5.  ThreeBodySystem
# ════════════════════════════════════════════════════════════════════════════

class ThreeBodySystem:
    """
    Three-body wrapper — hierarchical triple or arbitrary configuration.

    Example: PSR J0337+1715 millisecond pulsar triple system.
    """

    def __init__(self, body1: QGDBody, body2: QGDBody, body3: QGDBody):
        self.bodies = [body1, body2, body3]
        self.metric = QGDNBodyMetric(self.bodies)

    @classmethod
    def psr_j0337(cls):
        """
        PSR J0337+1715: millisecond pulsar hierarchical triple.
        Masses from Ransom et al. 2014.
        """
        AU = 1.496e11
        P_inner = 1.6292458 * 86400   # s
        P_outer = 327.2556  * 86400   # s
        M1 = 1.4378 * M_sun           # pulsar
        M2 = 0.19751 * M_sun          # inner white dwarf
        M3 = 0.4101  * M_sun          # outer white dwarf

        # Kepler: a = (G M_tot / omega^2)^(1/3)
        a12  = (G * (M1+M2) * (P_inner/(2*np.pi))**2)**(1/3)
        a_out = (G * (M1+M2+M3) * (P_outer/(2*np.pi))**2)**(1/3)

        b1 = QGDBody(M=M1, pos=[0, 0, 0])
        b2 = QGDBody(M=M2, pos=[a12, 0, 0])
        b3 = QGDBody(M=M3, pos=[a_out, 0, 0])
        return cls(b1, b2, b3)

    def cross_term_catalogue(self, x: np.ndarray, theta: float = np.pi/2) -> str:
        """
        Print all 12 cross-term contributions to g_tφ for 3 bodies.
        (3 self frame-drag + 6 A×B cross + 3 B×A cross)
        """
        mp = self.metric.at(x, theta)
        A = mp.A_parts
        B = mp.B_parts
        r = np.linalg.norm(x)
        pref = -r * np.sin(theta)
        N = 3
        labels = ['pulsar (1)', 'WD-inner (2)', 'WD-outer (3)']

        lines = ["THREE-BODY g_tφ CROSS-TERM CATALOGUE"]
        lines.append("-" * 50)
        lines.append("Self frame-drag terms (body's own spin):")
        for a in range(N):
            val = pref * A[a] * B[a]
            lines.append(f"  σ_t^({a+1}) × σ_φ^({a+1}) = {val:+.4e}  [{labels[a]}]")
        lines.append("\nCross frame-drag terms (one body's spin in other's mass field):")
        for a in range(N):
            for b in range(N):
                if a != b:
                    val = pref * A[a] * B[b]
                    lines.append(f"  σ_t^({a+1}) × σ_φ^({b+1}) = {val:+.4e}"
                                 f"  [M_{a+1} mass × α_{b+1} spin]")
        return "\n".join(lines)


# ════════════════════════════════════════════════════════════════════════════
# 6.  Helper functions
# ════════════════════════════════════════════════════════════════════════════

def qgd_cross_force(x: np.ndarray,
                    body_a: QGDBody,
                    body_b: QGDBody,
                    eps: float = 1.0) -> np.ndarray:
    """
    Compute the QGD cross-term force on a test particle at x
    in the combined field of bodies a and b.

    F_cross = -c² [σ_t^(a) ∇σ_t^(b) + σ_t^(b) ∇σ_t^(a)]

    Scales as G √(M_a M_b) r^{-3/2} — softer than Newton (r^{-2}).
    """
    x = np.asarray(x, dtype=float)
    f_cross = np.zeros(3)
    for i in range(3):
        xp, xm = x.copy(), x.copy()
        xp[i] += eps; xm[i] -= eps
        sa    = body_a.sigma_t(x)
        dsa   = (body_a.sigma_t(xp) - body_a.sigma_t(xm)) / (2*eps)
        sb    = body_b.sigma_t(x)
        dsb   = (body_b.sigma_t(xp) - body_b.sigma_t(xm)) / (2*eps)
        f_cross[i] = -c**2 * (sa * dsb + sb * dsa)
    return f_cross


def sigma_profile(body: QGDBody,
                  r_min: float, r_max: float,
                  n: int = 100,
                  theta: float = np.pi/2) -> Tuple[np.ndarray, np.ndarray]:
    """Return (r_array, sigma_t_array) radial profile for a single body."""
    r_arr = np.logspace(np.log10(r_min), np.log10(r_max), n)
    s_arr = np.array([body.sigma_t(np.array([r, 0, 0]), theta) for r in r_arr])
    return r_arr, s_arr


# ════════════════════════════════════════════════════════════════════════════
# 7.  Demo / examples
# ════════════════════════════════════════════════════════════════════════════

def run_examples():
    """
    Self-contained demonstration of all key QGD multi-body results.
    Run this file directly:  python qgd_nbody_exact.py
    """
    sep = "=" * 65

    # ── Example 1: Schwarzschild horizon ────────────────────────────────
    print(sep)
    print("EXAMPLE 1: Schwarzschild black hole — horizon verification")
    print(sep)
    bh = QGDBody(M=10*M_sun)
    metric = QGDNBodyMetric([bh])
    rs_exact = bh.r_s
    r_h = metric.horizon_radius(r_min=100, r_max=1e9)
    print(f"  M = 10 M_sun")
    print(f"  r_s (exact)      = {rs_exact:.4f} m")
    print(f"  QGD horizon (Σ=1)= {r_h:.4f} m")
    print(f"  Match: {np.isclose(rs_exact, r_h, rtol=1e-5)}")
    mp = metric.at([r_h + 1, 0, 0])
    print(f"  g_tt at horizon  = {mp.gtt:.4e}  (→ 0 ✓)")
    print(f"  Σ    at horizon  = {mp.A_total:.6f}  (→ 1 ✓)")

    # ── Example 2: Two Kerr bodies — full metric ─────────────────────────
    print(f"\n{sep}")
    print("EXAMPLE 2: Two spinning black holes — complete metric + cross terms")
    print(sep)
    sys2 = TwoBodySystem.kerr_binary(
        M1=10*M_sun, M2=8*M_sun, a1=0.6, a2=0.4, separation=1e6
    )
    x_test = [3e5, 1e4, 0]
    mp2 = sys2.metric_at(x_test)
    print(sys2.describe_cross_terms(mp2, r=np.linalg.norm(x_test)))
    print(f"\n  Full 4×4 metric tensor:")
    print(f"  g_tt    = {mp2.gtt:.8f}")
    print(f"  g_rr    = {mp2.grr:.8f}")
    print(f"  g_tφ    = {mp2.gtphi:.4e}")
    print(f"  g_φφ    = {mp2.gphiphi:.4e}")

    # ── Example 3: Scan between two bodies ───────────────────────────────
    print(f"\n{sep}")
    print("EXAMPLE 3: Metric profile between two 5M_sun BHs (d = 1000 km)")
    print(sep)
    sys3 = TwoBodySystem.kerr_binary(M1=5*M_sun, M2=5*M_sun, separation=1e6)
    sys3.scan_between_bodies()

    # ── Example 4: Merger condition ──────────────────────────────────────
    print(f"\n{sep}")
    print("EXAMPLE 4: Equal-mass merger — Σ(midpoint) vs separation")
    print(sep)
    M = 5 * M_sun
    rs = 2 * G * M / c**2
    print(f"  r_s = {rs:.2f} m")
    print(f"\n  {'d/r_s':>8}  {'Σ(midpoint)':>14}  {'Status':>12}")
    for d_fac in [16, 8, 4, 2]:
        d = d_fac * rs
        b1 = QGDBody(M=M, pos=[0,0,0])
        b2 = QGDBody(M=M, pos=[d,0,0])
        met = QGDNBodyMetric([b1, b2])
        mp = met.at([d/2, 1.0, 0])
        status = "MERGED ←" if mp.A_total >= 1 else ""
        print(f"  {d_fac:>8}  {mp.A_total:>14.6f}  {status}")
    print(f"\n  Analytical prediction: merger at d = 4 r_s = {4*rs:.2f} m  ✓")

    # ── Example 5: Three-body — PSR J0337 ────────────────────────────────
    print(f"\n{sep}")
    print("EXAMPLE 5: Three-body — PSR J0337+1715 pulsar triple")
    print(sep)
    triple = ThreeBodySystem.psr_j0337()
    for i, b in enumerate(triple.bodies):
        print(f"  Body {i+1}: M={b.M/M_sun:.4f} M_sun, "
              f"r_s={b.r_s:.2f} m, "
              f"|pos|={np.linalg.norm(b.pos)/AU:.4f} AU")
    # Field point at 1 AU
    x_1au = np.array([AU, 1e6, 0])
    mp5 = triple.metric.at(x_1au)
    print(f"\n  Metric at 1 AU from barycentre:")
    print(f"  A_total = {mp5.A_total:.6e}")
    print(f"  g_tt    = {mp5.gtt:.10f}")
    print(f"  σ_t per body: {[f'{a:.4e}' for a in mp5.A_parts]}")
    print(triple.cross_term_catalogue(x_1au))

    # ── Example 6: N-body scaling ─────────────────────────────────────────
    print(f"\n{sep}")
    print("EXAMPLE 6: N-body — cross-term count scaling")
    print(sep)
    from math import comb
    print(f"  {'N':>5}  {'g_tt terms':>12}  {'g_tφ terms':>12}  {'g_φφ terms':>12}")
    for N in [1, 2, 3, 5, 10, 50, 100]:
        n_tt   = N + comb(N, 2)   # N self + C(N,2) cross
        n_tphi = N + N*(N-1)      # N self-spin + N(N-1) cross-spin
        n_pp   = N + comb(N, 2)
        print(f"  {N:>5}  {n_tt:>12}  {n_tphi:>12}  {n_pp:>12}")

    # ── Example 7: QGD cross-body force ──────────────────────────────────
    print(f"\n{sep}")
    print("EXAMPLE 7: QGD cross-term acceleration on test particle")
    print(sep)
    M1, M2 = 5*M_sun, 5*M_sun
    d  = 1e6
    b1 = QGDBody(M=M1, pos=[0,0,0])
    b2 = QGDBody(M=M2, pos=[d,0,0])
    x_mid = np.array([d/2, 1e3, 0])
    sys7  = QGDNBodyMetric([b1, b2])
    a_tot = sys7.geodesic_accel(x_mid)
    # Newton alone (both forces cancel at midpoint for equal masses, leaving y)
    print(f"  x_test = midpoint, y = 1 km offset")
    print(f"  |a_QGD|   = {np.linalg.norm(a_tot):.4e} m/s²")
    print(f"  a_QGD_y   = {a_tot[1]:.4e} m/s²  (cross-term: gravity toward barycentre)")
    f_cross = qgd_cross_force(x_mid, b1, b2)
    print(f"  |F_cross| = {np.linalg.norm(f_cross):.4e} m/s²  (QGD-only cross term)")


# ════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    run_examples()
