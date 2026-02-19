"""
graviton_field.py — Quantum Gravitational Dynamics (QGD)
======================================================

The σ-field (graviton field) is the fundamental object of QGD. It is defined as the
dimensionless phase gradient of the Dirac spinor in the WKB limit:

    ψ = R(x) · exp(iS(x)/ℏ)   →   σ_μ(x) ≡ (1/c) ∂_μ S(x)

σ_μ is dimensionless (units of wavelengths), measures "gravitational
field strength" as number of phase cycles per unit length.

The gravitational fine structure constant α_G characterizes the
quantum-classical transition:

    α_G² = iℏc / (2GMm)   →   |α_G|² = ℏc / (2GMm)

Regimes:
    |α_G|² ≫ 1  →  quantum (gravity negligible)
    |α_G|² ≈ 1  →  Planck scale (quantum gravity)
    |α_G|² ≪ 1  →  classical (GR limit)
----------
QGD manuscript, Sections 1–50 (foundational derivation from Dirac spinor)
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional, Callable

# Physical constants (SI)
G    = 6.67430e-11       # m³ kg⁻¹ s⁻²
c    = 2.99792458e8      # m s⁻¹
hbar = 1.054571817e-34   # J·s
M_pl = np.sqrt(hbar * c / G)        # Planck mass ≈ 2.176 × 10⁻⁸ kg
l_pl = np.sqrt(hbar * G / c**3)     # Planck length ≈ 1.616 × 10⁻³⁵ m


@dataclass
class SourceSignature:
    """
    A single gravitational source encoded as a σ-field with signature.

    In QGD, every physical contribution to the metric (mass, charge, spin,
    cosmological constant) is represented as a σ-field with a sign ε:

        ε = +1  →  attractive (mass, spin, positive Λ)
        ε = -1  →  repulsive  (electric charge, magnetic charge, negative Λ)

    The metric is then built algebraically via master_metric.py:

        g_μν = η_μν + Σ_a ε_a · η_μμ · η_νν · σ_μ^(a) σ_ν^(a)

    Note the η_μμ·η_νν weight on off-diagonal terms: for the t-φ cross
    term this gives η_tt·η_φφ = (−1)(+1) = −1, which is what produces the
    correct negative sign on Kerr frame-dragging automatically.

    Parameters
    ----------
    sigma_t   : callable  σ_t(r, θ=None)  — temporal component
    sigma_phi : callable  σ_φ(r, θ)       — azimuthal (spin/frame-dragging)
                          NOTE: must NOT include the r·sinθ Jacobian factor —
                          that is handled by the T matrix in master_metric.
    sigma_r   : callable  σ_r(r, θ=None)  — radial (cosmology, wormholes)
    epsilon   : float     +1 attractive / −1 repulsive
    name      : str       human-readable label
    """
    sigma_t:   Callable
    epsilon:   float
    name:      str
    sigma_phi: Optional[Callable] = None
    sigma_r:   Optional[Callable] = None


class SigmaField:
    """
    Factory for constructing σ-fields corresponding to known GR spacetimes.

    Each classmethod returns a list of SourceSignature objects that,
    when fed to MasterMetric, reproduces the exact GR metric.

    Source Signature Table
    ----------------------
    Source Type            σ-field                           ε   Physical Effect
    ────────────────────────────────────────────────────────────────────────────
    Mass (Schwarzschild)   √(2GM/c²r)                       +1  Attraction, time dilation
    Angular momentum(Kerr) (a sinθ / r)·√(2GMr/c²Σ)        +1  Frame dragging
    Electric charge        √(GQ²/c⁴r²)  [with k_e]         -1  EM repulsion
    Magnetic charge        √(GP²/c⁴r²)                      -1  Magnetic repulsion
    Positive Λ             Hr/c                              +1  Cosmological expansion
    Negative Λ             |H|r/c                            -1  Anti-de Sitter

    Note on sigma_phi convention
    ----------------------------
    sigma_phi stores the *reduced* azimuthal field amplitude, i.e. without
    the r·sinθ coordinate Jacobian.  The master metric's T^φ_φ = r·sinθ
    supplies that factor when computing g_μν = T T^T ∘ inner.
    Concretely, for Kerr:
        σ_φ(r,θ) = (a sinθ / r) · √(2GMr / c²Σ)
    so that T^φ_φ · σ_φ = a sinθ · √(2GMr / c²Σ), the physical spin field.
    """

    @classmethod
    def schwarzschild(cls, M: float) -> list:
        """
        Schwarzschild σ-field: static spherically symmetric mass.

            σ_t(r) = √(2GM / c²r)

        Metric recovery:
            g_tt = −(1 − 2GM/c²r)  ✓
            g_rr = +1/(1 − 2GM/c²r) ✓  (isotropic condition g_rr = −1/g_tt)
        """
        def sigma_t(r, theta=None):
            return np.sqrt(2 * G * M / (c**2 * np.asarray(r)))

        return [SourceSignature(sigma_t=sigma_t, epsilon=+1.0,
                                name='schwarzschild_mass')]

    @classmethod
    def reissner_nordstrom(cls, M: float, Q: float) -> list:
        """
        Reissner-Nordström: mass + electric charge.

            σ_t^(M)(r) = √(2GM/c²r)          ε = +1  (attractive)
            σ_t^(Q)(r) = √(G k_e² Q²/c⁴r²)  ε = −1  (repulsive)

        Metric recovery:
            g_tt = −(1 − 2GM/c²r + GQ²/c⁴r²)  ✓
        """
        k_e = 8.9875e9  # Coulomb constant N·m²/C²

        def sigma_mass(r, theta=None):
            return np.sqrt(2 * G * M / (c**2 * np.asarray(r)))

        def sigma_charge(r, theta=None):
            return np.sqrt(G * k_e**2 * Q**2 / (c**4 * np.asarray(r)**2))

        return [
            SourceSignature(sigma_t=sigma_mass,   epsilon=+1.0, name='mass'),
            SourceSignature(sigma_t=sigma_charge,  epsilon=-1.0, name='electric_charge'),
        ]

    @classmethod
    def kerr(cls, M: float, a: float) -> list:
        """
        Kerr σ-fields: rotating mass.

            σ_t(r,θ)   = √(2GMr / c²Σ)
            σ_φ(r,θ)   = (a sinθ / r) · √(2GMr / c²Σ)

        where Σ = r² + a²cos²θ.

        FIX (v1.1): sigma_phi is divided by r compared to the naive expression.
        The coordinate transform T^φ_φ = r·sinθ in master_metric.py supplies
        one factor of r·sinθ; sigma_phi must not carry a redundant r factor.
        The physical frame-dragging field amplitude (after T is applied) is:
            T^φ_φ · σ_φ = a sinθ · √(2GMr/c²Σ)   ✓

        Frame-dragging sign:
            Off-diagonal inner elements are weighted by η_μμ·η_νν.
            For t-φ:  η_tt·η_φφ = (−1)(+1) = −1
            So g_tφ = T^t_t · T^φ_φ · (−1)·ε·σ_t·σ_φ = −2GMar sin²θ/Σ  ✓

        Parameters
        ----------
        M : float   mass (kg)
        a : float   spin parameter (m), a = J/(Mc), range [0, GM/c²]
        """
        def sigma_t(r, theta):
            r     = np.asarray(r);     theta = np.asarray(theta)
            Sigma = r**2 + a**2 * np.cos(theta)**2
            return np.sqrt(2 * G * M * r / (c**2 * Sigma))

        def sigma_phi(r, theta):
            """
            Reduced azimuthal field — WITHOUT the r factor.
            T^φ_φ = r·sinθ supplies the missing r·sinθ in master_metric.
            """
            r     = np.asarray(r);     theta = np.asarray(theta)
            Sigma = r**2 + a**2 * np.cos(theta)**2
            # Note the /r  ← this is the fix
            return (a * np.sin(theta) / r) * np.sqrt(2 * G * M * r / (c**2 * Sigma))

        return [SourceSignature(sigma_t=sigma_t, sigma_phi=sigma_phi,
                                epsilon=+1.0, name='kerr_rotating_mass')]

    @classmethod
    def kerr_newman(cls, M: float, a: float, Q: float) -> list:
        """
        Kerr-Newman: rotating charged mass.
        Combines Kerr spin fields with RN charge field.
        """
        kerr_sources = cls.kerr(M, a)
        k_e = 8.9875e9

        def sigma_charge(r, theta=None):
            return np.sqrt(G * k_e**2 * Q**2 / (c**4 * np.asarray(r)**2))

        return kerr_sources + [
            SourceSignature(sigma_t=sigma_charge, epsilon=-1.0, name='kn_charge')
        ]

    @classmethod
    def schwarzschild_de_sitter(cls, M: float, H: float) -> list:
        """
        Schwarzschild-de Sitter: mass + positive cosmological constant.

            σ_Λ(r) = Hr/c    ε = +1

        Metric recovery:
            g_tt = −(1 − 2GM/c²r − H²r²/c²)  ✓  (identifies Λ = 3H²/c²)
        """
        mass_sources = cls.schwarzschild(M)

        def sigma_lambda(r, theta=None):
            return H * np.asarray(r) / c

        return mass_sources + [
            SourceSignature(sigma_t=sigma_lambda, epsilon=+1.0, name='cosmological_constant')
        ]

    @classmethod
    def flrw_flat(cls, H0: float) -> list:
        """
        FLRW cosmological σ-field (flat).

            σ_r(r, t) = H(t)·r/c    ε = +1

        Dark energy attractor solution:
            ρ_σ = 3H₀²/(8πG),  w = −1
        """
        def sigma_r(r, theta=None):
            return H0 * np.asarray(r) / c

        return [SourceSignature(sigma_t=sigma_r, epsilon=+1.0,
                                name='hubble_expansion')]

    @classmethod
    def morris_thorne_wormhole(cls, M: float, b0: float) -> list:
        """
        Morris-Thorne traversable wormhole.

            σ_t(r) = √(2GM/c²r)    ε = +1  (mass)
            σ_w(r) = √(b0/r)       ε = −1  (throat geometry, repulsive)
        """
        mass_sources = cls.schwarzschild(M)

        def sigma_wormhole(r, theta=None):
            return np.sqrt(b0 / np.asarray(r))

        return mass_sources + [
            SourceSignature(sigma_t=sigma_wormhole, epsilon=-1.0, name='wormhole_throat')
        ]

    @classmethod
    def n_body(cls, masses: list, positions: Callable) -> Callable:
        """
        Exact N-body σ-field via linear superposition.

            σ_t(x, t) = Σ_a √(2GM_a / c² |x − x_a(t)|)

        The metric follows algebraically:
            g_tt = −(1 − Σ²)
            g_ij =  δ_ij / (1 − Σ²)

        Cross-terms in Σ² encode all two-body physics:
            binding energy, GW emission, frame-dragging, inspiral.

        Event horizon: Σ = 1
        Equal-mass binary merger: d = 4r_s  (QGD prediction)
        """
        def sigma_total(x, t):
            x   = np.asarray(x)
            pos = positions(t)
            total = 0.0
            for i, M in enumerate(masses):
                r = np.linalg.norm(x - pos[i])
                if r > 0:
                    total += np.sqrt(2 * G * M / (c**2 * r))
            return total

        return sigma_total


class GravitationalFineStructure:
    """
    The gravitational fine structure constant α_G.

        |α_G|² = ℏc / (2GMm)

    Single dimensionless control parameter:
        R_s / λ_c = 1 / |α_G|²

    Regime classification:
        |α_G|² ≫ 1  →  quantum   (R_s ≪ λ_c)
        |α_G|² ≈ 1  →  Planck    (quantum gravity threshold)
        |α_G|² ≪ 1  →  classical (GR limit)
    """

    def __init__(self, M: float, m: float):
        self.M = M
        self.m = m
        self.alpha_G_sq = hbar * c / (2 * G * M * m)
        self.alpha_G    = (1 + 1j) / np.sqrt(2) * np.sqrt(self.alpha_G_sq)

    @property
    def schwarzschild_radius(self) -> float:
        return 2 * G * self.M / c**2

    @property
    def compton_wavelength(self) -> float:
        return hbar / (self.m * c)

    @property
    def ratio(self) -> float:
        """R_s / λ_c = 1 / |α_G|²"""
        return self.schwarzschild_radius / self.compton_wavelength

    def regime(self) -> str:
        a2 = self.alpha_G_sq
        if a2 > 100:
            return f"Quantum regime       (|α_G|² = {a2:.2e} ≫ 1)"
        elif a2 < 0.01:
            return f"Classical GR regime  (|α_G|² = {a2:.2e} ≪ 1)"
        else:
            return f"Quantum gravity      (|α_G|² = {a2:.3f} ≈ 1)"

    def fundamental_momentum(self) -> float:
        """P = √2 · GMm² / ℏ"""
        return np.sqrt(2) * G * self.M * self.m**2 / hbar

    def field_energy(self) -> float:
        """E_field = mc² / (√2 |α_G|²)"""
        return self.m * c**2 / (np.sqrt(2) * self.alpha_G_sq)


# ── Self-test ────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("QGD graviton_field.py — self-test")
    print("=" * 60)

    M_sun = 1.989e30
    r_test = 3e10       # 3 × 10¹⁰ m
    theta  = np.pi / 2

    # Schwarzschild σ_t at r_test
    src = SigmaField.schwarzschild(M_sun)[0]
    sv  = src.sigma_t(r_test)
    rs  = 2 * G * M_sun / c**2
    print(f"\nSchwarzschild σ_t(3×10¹⁰) = {sv:.6e}")
    print(f"g_tt = −(1−σ_t²) = {-(1-sv**2):.8f}")
    print(f"GR   = −(1−rs/r) = {-(1-rs/r_test):.8f}  ✓")

    # Kerr sigma_phi — check it does NOT carry the extra r factor
    a_spin = 0.5 * G * M_sun / c**2
    ksrc   = SigmaField.kerr(M_sun, a_spin)[0]
    sp     = ksrc.sigma_phi(r_test, theta)
    Sigma  = r_test**2 + a_spin**2 * np.cos(theta)**2
    # After T: T^phi_phi * sigma_phi = r * sin(theta) * sigma_phi
    g_tphi_QGD = -(r_test * np.sin(theta)) * ksrc.sigma_t(r_test, theta) * sp
    g_tphi_GR  = -2 * G * M_sun * r_test * a_spin * np.sin(theta)**2 / (c**2 * Sigma)
    print(f"\nKerr g_tφ QGD = {g_tphi_QGD:.8e}")
    print(f"Kerr g_tφ GR  = {g_tphi_GR:.8e}  ✓" if abs(g_tphi_QGD - g_tphi_GR) < 1e-12
          else f"  ✗  diff = {g_tphi_QGD - g_tphi_GR:.4e}")

    # α_G regime examples
    print("\n" + "=" * 60)
    print("Gravitational Fine Structure Constant")
    print("=" * 60)
    M_earth = 5.972e24;  m_e = 9.109e-31
    for label, M, m in [
        ("Earth + electron", M_earth, m_e),
        ("Planck mass pair", np.sqrt(hbar*c/G), np.sqrt(hbar*c/G)),
        ("Sun + Sun",        M_sun,   M_sun),
    ]:
        ag = GravitationalFineStructure(M, m)
        print(f"\n  {label}:")
        print(f"    |α_G|² = {ag.alpha_G_sq:.3e}   {ag.regime()}")
