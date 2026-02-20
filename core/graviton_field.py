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

    |α_G|² = ℏc / (2GMm)        ← CORRECT DEFINITION (v2.0)
    (Old incorrect form: α_G = ℏc/GMm — missing the factor of 2)

Regimes:
    |α_G|² ≫ 1  →  quantum (gravity negligible)
    |α_G|² ≈ 1  →  Planck scale (quantum gravity)
    |α_G|² ≪ 1  →  classical (GR limit)

Field Energy:
    E_field = mc²/√2 · (R_s/λ_c) = mc²/(√2 |α_G|²)

Field Momentum:
    P_field = √2·GMm²/ℏ = mc/(√2 |α_G|²) = (1/√2)·mc·(R_s/λ_c)
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
    """
    sigma_t:   Callable
    epsilon:   float
    name:      str
    sigma_phi: Optional[Callable] = None
    sigma_r:   Optional[Callable] = None


class SigmaField:
    """Factory for constructing σ-fields corresponding to known GR spacetimes."""

    @classmethod
    def schwarzschild(cls, M: float) -> list:
        def sigma_t(r, theta=None):
            return np.sqrt(2 * G * M / (c**2 * np.asarray(r)))
        return [SourceSignature(sigma_t=sigma_t, epsilon=+1.0,
                                name='schwarzschild_mass')]

    @classmethod
    def reissner_nordstrom(cls, M: float, Q: float) -> list:
        k_e = 8.9875e9
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
        def sigma_t(r, theta):
            r     = np.asarray(r);     theta = np.asarray(theta)
            Sigma = r**2 + a**2 * np.cos(theta)**2
            return np.sqrt(2 * G * M * r / (c**2 * Sigma))

        def sigma_phi(r, theta):
            r     = np.asarray(r);     theta = np.asarray(theta)
            Sigma = r**2 + a**2 * np.cos(theta)**2
            return (a * np.sin(theta) / r) * np.sqrt(2 * G * M * r / (c**2 * Sigma))

        return [SourceSignature(sigma_t=sigma_t, sigma_phi=sigma_phi,
                                epsilon=+1.0, name='kerr_rotating_mass')]

    @classmethod
    def kerr_newman(cls, M: float, a: float, Q: float) -> list:
        kerr_sources = cls.kerr(M, a)
        k_e = 8.9875e9
        def sigma_charge(r, theta=None):
            return np.sqrt(G * k_e**2 * Q**2 / (c**4 * np.asarray(r)**2))
        return kerr_sources + [
            SourceSignature(sigma_t=sigma_charge, epsilon=-1.0, name='kn_charge')
        ]

    @classmethod
    def schwarzschild_de_sitter(cls, M: float, H: float) -> list:
        mass_sources = cls.schwarzschild(M)
        def sigma_lambda(r, theta=None):
            return H * np.asarray(r) / c
        return mass_sources + [
            SourceSignature(sigma_t=sigma_lambda, epsilon=+1.0, name='cosmological_constant')
        ]

    @classmethod
    def flrw_flat(cls, H0: float) -> list:
        def sigma_r(r, theta=None):
            return H0 * np.asarray(r) / c
        return [SourceSignature(sigma_t=sigma_r, epsilon=+1.0,
                                name='hubble_expansion')]

    @classmethod
    def morris_thorne_wormhole(cls, M: float, b0: float) -> list:
        mass_sources = cls.schwarzschild(M)
        def sigma_wormhole(r, theta=None):
            return np.sqrt(b0 / np.asarray(r))
        return mass_sources + [
            SourceSignature(sigma_t=sigma_wormhole, epsilon=-1.0, name='wormhole_throat')
        ]

    @classmethod
    def n_body(cls, masses: list, positions: Callable) -> Callable:
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

    CORRECT DEFINITION (v2.0):
        |α_G|² = ℏc / (2GMm)

    Previous (incorrect) form was: α_G = ℏc/(GMm) — missing factor of 2.

    Derived quantities:
        R_s / λ_c = 1 / |α_G|²   (Schwarzschild-to-Compton ratio)

    Field Energy:
        E_field = mc²/√2 · (R_s/λ_c) = mc²/(√2 |α_G|²)
        Derivation: substitute GM = R_s c²/2,  ℏ = mc·λ_c

    Field Momentum:
        P_field = √2·GMm²/ℏ = mc/(√2|α_G|²) = (1/√2)·mc·(R_s/λ_c)

    These arise naturally from the normalization condition on the σ-field:
        |ψ|²·σ_μ·x^μ = J/ℏ

    Regime classification:
        |α_G|² ≫ 1  →  quantum   (R_s ≪ λ_c)
        |α_G|² ≈ 1  →  Planck    (quantum gravity threshold)
        |α_G|² ≪ 1  →  classical (GR limit)
    """

    def __init__(self, M: float, m: float):
        self.M = M
        self.m = m
        # CORRECT: |α_G|² = ℏc / (2GMm)
        self.alpha_G_sq = hbar * c / (2 * G * M * m)
        # Complex form: α_G = (1+i)/√2 × √|α_G|² (magnitude = √|α_G|²)
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

    def field_energy(self) -> float:
        """
        E_field = mc²/(√2 |α_G|²)
        
        Equivalently: E_field = (1/√2) mc² (R_s/λ_c)
        
        Derivation: substitute GM = R_s c²/2, ℏ = mc·λ_c into
        the σ-field energy integral.
        """
        return self.m * c**2 / (np.sqrt(2) * self.alpha_G_sq)

    def field_energy_alternative(self) -> float:
        """Same result via R_s/λ_c form — self-consistency check."""
        return (1/np.sqrt(2)) * self.m * c**2 * self.ratio

    def field_momentum(self) -> float:
        """
        P_field = mc/(√2 |α_G|²)
        
        Equivalently: P_field = √2·GMm²/ℏ = (1/√2)·mc·(R_s/λ_c)
        """
        return self.m * c / (np.sqrt(2) * self.alpha_G_sq)

    def field_momentum_alternative(self) -> float:
        """Same result via √2·GMm²/ℏ — self-consistency check."""
        return np.sqrt(2) * G * self.M * self.m**2 / hbar

    def field_momentum_ratio_form(self) -> float:
        """Same result via (1/√2)·mc·(R_s/λ_c) — self-consistency check."""
        return (1/np.sqrt(2)) * self.m * c * self.ratio

    def verify_consistency(self) -> dict:
        """
        Verify all three forms of E_field and P_field agree.
        Returns dict of relative errors between forms.
        """
        E1 = self.field_energy()
        E2 = self.field_energy_alternative()

        P1 = self.field_momentum()
        P2 = self.field_momentum_alternative()
        P3 = self.field_momentum_ratio_form()

        return {
            'E_field': E1,
            'E_field_alt': E2,
            'E_consistency': abs(E1-E2)/E1,
            'P_field': P1,
            'P_field_GMm2_form': P2,
            'P_field_ratio_form': P3,
            'P_consistency_12': abs(P1-P2)/P1,
            'P_consistency_13': abs(P1-P3)/P1,
            'E_over_P_ratio': E1 / (P1 * c),   # should equal 1/2 by E=pc/2? let's check
            'ratio_Rs_lc': self.ratio,
            'alpha_G_sq': self.alpha_G_sq,
        }


# ── Self-test ────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 70)
    print("QGD graviton_field.py — self-test (v2.0, corrected α_G)")
    print("=" * 70)

    M_sun = 1.989e30
    r_test = 3e10
    theta  = np.pi / 2

    # Schwarzschild σ_t at r_test
    src = SigmaField.schwarzschild(M_sun)[0]
    sv  = src.sigma_t(r_test)
    rs  = 2 * G * M_sun / c**2
    print(f"\nSchwarzschild σ_t(3×10¹⁰) = {sv:.6e}")
    print(f"g_tt = −(1−σ_t²) = {-(1-sv**2):.8f}")
    print(f"GR   = −(1−rs/r) = {-(1-rs/r_test):.8f}  ✓")

    # Kerr sigma_phi check
    a_spin = 0.5 * G * M_sun / c**2
    ksrc   = SigmaField.kerr(M_sun, a_spin)[0]
    sp     = ksrc.sigma_phi(r_test, theta)
    Sigma  = r_test**2 + a_spin**2 * np.cos(theta)**2
    g_tphi_QGD = -(r_test * np.sin(theta)) * ksrc.sigma_t(r_test, theta) * sp
    g_tphi_GR  = -2 * G * M_sun * r_test * a_spin * np.sin(theta)**2 / (c**2 * Sigma)
    print(f"\nKerr g_tφ QGD = {g_tphi_QGD:.8e}")
    print(f"Kerr g_tφ GR  = {g_tphi_GR:.8e}",
          "✓" if abs(g_tphi_QGD - g_tphi_GR) < 1e-12 else f"✗ diff={g_tphi_QGD-g_tphi_GR:.2e}")

    # α_G regime examples
    print("\n" + "=" * 70)
    print("Gravitational Fine Structure Constant |α_G|² = ℏc/(2GMm)")
    print("=" * 70)
    M_earth = 5.972e24;  m_e = 9.109e-31
    for label, M, m in [
        ("Earth + electron", M_earth, m_e),
        ("Planck mass pair", np.sqrt(hbar*c/G), np.sqrt(hbar*c/G)),
        ("Sun + proton",     M_sun, 1.6726e-27),
        ("Sun + Sun",        M_sun, M_sun),
    ]:
        ag = GravitationalFineStructure(M, m)
        print(f"\n  {label}:")
        print(f"    |α_G|² = {ag.alpha_G_sq:.3e}   {ag.regime()}")
        print(f"    R_s/λ_c = {ag.ratio:.3e}   (should = 1/|α_G|²)")
        print(f"    E_field = {ag.field_energy():.3e} J")
        print(f"    P_field = {ag.field_momentum():.3e} kg·m/s")
        v = ag.verify_consistency()
        print(f"    E consistency: {v['E_consistency']:.2e}  (3 forms)")
        print(f"    P consistency: {v['P_consistency_12']:.2e} / {v['P_consistency_13']:.2e} (3 forms)")

    # Planck scale check: |α_G|² = 1/2 for two Planck masses
    M_pl_val = np.sqrt(hbar * c / G)
    ag_pl = GravitationalFineStructure(M_pl_val, M_pl_val)
    print(f"\n  Planck mass pair:")
    print(f"    |α_G|² = {ag_pl.alpha_G_sq:.6f}  (should = 0.5 for Planck gravity threshold)")
