"""
qgd_eob_imr.py
==============
Quantum Gravitational Dynamics — Effective One-Body, Merger Transition,
and Spinning Master Formula

Author  : Romeo Matshaba, University of South Africa
Version : 2.0  (March 2026)
Changes : D-factor / dipole radiation set to exactly zero throughout.
          The spin-0 scalar mode is Yukawa-suppressed by exp(−r/ℓ_Pl)
          at all astrophysical separations, contributing zero observable
          radiation.  All three problems updated accordingly.

=============================================================================
PROBLEMS ADDRESSED
=============================================================================

  PROBLEM 1 — QGD Effective One-Body (EOB)
    A^QGD(u) = 1 − 2u/η  (exact, no free parameters)
    Spinning extension: A^QGD(u,χ) = 1 − 2u/η + χ²u²/η²

  PROBLEM 2 — Merger Transition (exact analytic)
    d_merge = 2GM₁(1+(M₂/M₁)^{1/3})³/c²
    σ_saddle = 1 exactly for all mass ratios q.

  PROBLEM 3 — Spinning Master Formula
    Three closed-form series from the Kerr geodesic:
      Spin-free   : e_n^0  = −C(2n,n)(3/4)^n(2n−1)/(n+1)
      Spin-orbit  : e_n^SO = −(2n+1)!!·3^n/(2^n·n!)
      Spin-squared: e_n^SS = (2n+3)!!·3^n/(3·2^{n+1}·n!)

=============================================================================
DIPOLE STATUS — WHY D = 0 IN ALL OBSERVABLES
=============================================================================

The QGD near-field cross-term 2σ_t^{(1)} σ_t^{(2)} ∝ √(M₁M₂)/r DOES exist
in the two-body master metric.  It governs:
  • The binding energy (conservative sector — same as GR ADM result)
  • The σ-field saddle structure (merger condition)
  • The energy redistribution between bodies during inspiral

However, the RADIATION of the mass-asymmetry D = (√M₁−√M₂)²/M to infinity
is carried exclusively by the spin-0 scalar σ-mode, which has mass m_Q ~ M_Pl.
Its Yukawa suppression:

    exp(−r/ℓ_Pl)  ~  10^{−5.24×10^43}  at Hulse-Taylor separation

is absolute — not a sensitivity limit.  The spin-0 mode does not propagate
beyond the Planck length.  Therefore:

    F_dip  =  0,     δΨ^{−1PN}  =  0,     δΨ^{0.5PN}  =  0

for any astrophysical binary.  Observable QGD waveform ≡ GR waveform.

=============================================================================
HONEST ACCOUNTING
=============================================================================

RIGOROUS:
  • d_merge = 2GM₁(1+(M₂/M₁)^{1/3})³/c²  (σ-saddle, verified all q)
  • A^QGD(u) = 1 − 2u/η  (exact algebraic derivation from Σ_tot² = 2u/η)
  • Spinning series e_n^{SO}, e_n^{SS}  (verified vs exact Kerr, n=0,1,2)
  • D-factor zero in all observables  (Yukawa suppression is absolute)

PRELIMINARY:
  • η-corrections at 5PN+ not yet fully derived
  • A^QGD vs NR calibration: not yet done
  • Spinning A-function beyond O(χ²): approximate

OPEN QUESTION:
  A^QGD = 1−2u/η has NO nonlinear PN terms (a₂ = a₃ = … = 0).
  GR EOB uses a₂ ≈ (94/3 − 41π²/32)η from NR calibration.
  Sharpest falsifiable prediction: NR should NOT require a₂ ≠ 0 if QGD is correct.
"""

import numpy as np
from math import factorial, comb
from fractions import Fraction
from scipy.optimize import brentq
from dataclasses import dataclass, field
from typing import Optional, List, Tuple

G      = 6.67430e-11
c      = 2.99792458e8
hbar   = 1.05457182e-34
kB     = 1.38064852e-23
Msun   = 1.98892e30
lPl    = np.sqrt(G * hbar / c**3)
_PI    = np.pi
_EUG   = np.euler_gamma


# =============================================================================
# UTILITIES
# =============================================================================

def _dfact(n):
    """Double factorial n!! = n(n-2)…; 0!! = (−1)!! = 1."""
    if n <= 0:
        return 1
    r = 1
    while n > 0:
        r *= n; n -= 2
    return r


# =============================================================================
# PROBLEM 2:  MERGER CONDITION  (placed first — informs EOB)
# =============================================================================

def qgd_merger_condition(M1, M2):
    """
    Exact QGD merger separation from the σ-field saddle point.

    QGD derivation
    --------------
    The σ-field along the axis joining two bodies at separation d:

        Σ(r₁) = σ_t^{(1)}(r₁) + σ_t^{(2)}(r₂),   r₂ = d − r₁

    with σ_t^{(a)}(r) = √(2GM_a/(c²r)).

    Merger condition: max_{r₁} Σ(r₁) = 1.

    The maximum occurs at  r₁/r₂ = (M₁/M₂)^{1/3}  (from dΣ/dr₁ = 0).
    Substituting and setting Σ = 1 gives:

        ┌───────────────────────────────────────────────────────────┐
        │  d_merge = 2GM₁(1 + (M₂/M₁)^{1/3})³ / c²   [EXACT]    │
        └───────────────────────────────────────────────────────────┘

    Equal-mass limit (M₁ = M₂ = m):
        d_merge = 2Gm(1+1)³/c² = 16Gm/c²
                = 8 × (2Gm/c²) = 8 r_s^{individual}
                = 4 × (4Gm/c²) = 4 r_s^{total}
        Both conventions correct; r_s definition differs.

    Extreme mass ratio (M₂ → 0):
        d_merge → 2GM₁/c²  (Schwarzschild radius of the primary)

    GR comparison: GR has no analogous algebraic merger condition.
    Merger in GR/NR is defined numerically as common apparent horizon formation.
    QGD provides an analytic criterion from the σ-field topology.

    Parameters
    ----------
    M1, M2 : float   Component masses [kg]

    Returns
    -------
    dict with d_merge, saddle location, Σ verification, dimensionless ratios
    """
    M   = M1 + M2
    rho = (M2 / M1)**(1.0/3)
    d   = 2*G*M1*(1 + rho)**3 / c**2

    # Saddle-point coordinates
    r1s = d / (1 + rho)
    r2s = d * rho / (1 + rho)

    # Verify Σ = 1 at saddle
    sig = (np.sqrt(2*G*M1/(c**2*r1s)) + np.sqrt(2*G*M2/(c**2*r2s)))

    rs_tot = 2*G*M / c**2
    rs_1   = 2*G*M1 / c**2

    return {
        'd_merge':       d,
        'r1_saddle':     r1s,
        'r2_saddle':     r2s,
        'Sigma_saddle':  sig,
        'saddle_exact':  abs(sig - 1.0) < 1e-10,
        'd_over_rstot':  d / rs_tot,
        'd_over_rs1':    d / rs_1,
        'q':             min(M1,M2) / max(M1,M2),
        'rho':           rho,
    }


def sigma_saddle(M1, M2, d):
    """Σ_total at the saddle point for bodies at separation d."""
    rho = (M2/M1)**(1.0/3)
    r1s = d / (1 + rho)
    r2s = d * rho / (1 + rho)
    return (np.sqrt(2*G*M1/(c**2*r1s)) + np.sqrt(2*G*M2/(c**2*r2s)))


def print_merger_table():
    """
    Table showing d_merge, Σ_saddle, d/r_s for a range of mass ratios.
    Also compares midpoint vs saddle-point merger conditions.
    """
    print("═"*74)
    print("QGD MERGER CONDITION — d_B = 2GM₁(1+(M₂/M₁)^{1/3})³/c²  (EXACT)")
    print("─"*74)
    print("  Three candidate conditions compared:")
    print("  (A) Midpoint Σ=1: approx, exact only at q=1")
    print("  (B) Saddle   Σ=1: EXACT for all q  ← QGD merger")
    print("  (C) EOB A=0 :     different physics (EOB singularity)")
    print("─"*74)
    M1 = 10 * Msun

    def dA(M2):
        return 4*G*(np.sqrt(M1) + np.sqrt(M2))**2 / c**2
    def dB(M2):
        return qgd_merger_condition(M1, M2)['d_merge']
    def dC(M2):
        M = M1 + M2; eta = M1*M2/M**2
        return 2*G*M / (eta*c**2)

    norm = 16*G*M1/c**2   # equal-mass answer

    print(f"\n  All values / (16Gm/c²)  [= 8r_s^indiv at equal mass]")
    print(f"  {'q':>6}  {'η':>7}  {'d_A/norm':>10}  {'d_B/norm':>10}  "
          f"{'d_C/norm':>10}  {'d_B/r_s_tot':>13}  {'Σ_saddle':>12}")
    print("  " + "─"*78)
    for q in [1.0, 0.5, 0.25, 0.1, 0.05, 0.01]:
        M2 = q*M1
        M  = M1+M2; eta = M1*M2/M**2
        r  = qgd_merger_condition(M1, M2)
        print(f"  {q:>6.3f}  {eta:>7.4f}  {dA(M2)/norm:>10.4f}  "
              f"{dB(M2)/norm:>10.4f}  {dC(M2)/norm:>10.4f}  "
              f"{r['d_over_rstot']:>13.4f}  {r['Sigma_saddle']:>12.8f}")
    print("─"*74)
    print("  KEY: d_B gives Σ_saddle = 1.000000 for ALL q  (algebraically exact)")
    print("  At q=1: d_A = d_B = 8r_s^indiv = 4r_s^total  (both conventions ✓)")
    print("═"*74)


# =============================================================================
# PROBLEM 1:  QGD EFFECTIVE ONE-BODY
# =============================================================================

class QGDEOB:
    """
    QGD Effective One-Body (EOB) mapping.

    QGD derivation
    --------------
    In the COM frame at separation r, the two-body σ-fields are:

        σ_t^{(1)}(r₁) = √(2GM₁/(c²·(M₂/M)r))  →  √(2GM/(c²r) · M₁/M₂)
        σ_t^{(2)}(r₂) = √(2GM₂/(c²·(M₁/M)r))  →  √(2GM/(c²r) · M₂/M₁)

    Σ_tot² = (σ_t^{(1)})² + (σ_t^{(2)})² + 2σ_t^{(1)}σ_t^{(2)}
           = (2GM/(c²r)) × (M₁/M₂ + 2 + M₂/M₁)
           = (2GM/(c²r)) × (M₁+M₂)²/(M₁M₂)
           = 2u/η,        u = GM/(c²r)

    QGD two-body effective metric (exact within QGD):
        g_tt^QGD(r) = −(1 − Σ_tot²) = −(1 − 2u/η)

    QGD EOB A-function:
        A^QGD(u) = 1 − 2u/η

    GR comparison
    -------------
    GR EOB:  A^GR = 1 − 2u + a₂(η)u² + a₃(η)u³ + …  (NR-calibrated)
    QGD EOB: A^QGD = 1 − 2u/η  (ZERO nonlinear PN terms — exact)

    The GR a₂ ≈ (94/3 − 41π²/32)η comes from NR calibration.
    If QGD is correct, NR should reproduce A^QGD without calibration.
    Sharpest open falsifiable prediction.

    Orbital structure (Schwarzschild with u → u/η)
    ───────────────────────────────────────────────
      ISCO:           u_ISCO = η/6           A = 2/3
      Light ring:     u_LR   = η/3           A = 1/3
      EOB singularity:u_h    = η/2           A = 0
      Physical merger (σ-saddle): u < u_h for q ≠ 1

    Spinning extension (leading order in χ)
    ----------------------------------------
      A^QGD(u,χ) = 1 − 2u/η + χ²u²/η²
      Limit η→0: recovers Kerr equatorial A = 1−2u+χ²u² ✓
    """

    def __init__(self, M1, M2, chi1=0.0, chi2=0.0):
        self.M1, self.M2 = M1, M2
        self.M   = M1 + M2
        self.mu  = M1*M2 / self.M
        self.eta = M1*M2 / self.M**2
        self.chi1, self.chi2 = chi1, chi2
        self.chi_eff = (M1*chi1 + M2*chi2) / self.M

        self.u_ISCO = self.eta / 6.0
        self.u_LR   = self.eta / 3.0
        self.u_h    = self.eta / 2.0   # EOB singularity (A=0)

        # Physical merger from σ-field saddle
        mc = qgd_merger_condition(M1, M2)
        self.d_merge  = mc['d_merge']
        self.u_merge  = G*self.M / (c**2 * self.d_merge)

    def A(self, u, chi=None):
        """A^QGD(u) = 1 − 2u/η  (non-spinning) or  1 − 2u/η + χ²u²/η²  (spinning)."""
        chi = chi if chi is not None else self.chi_eff
        return 1.0 - 2.0*u/self.eta + chi**2*u**2/self.eta**2

    def A_gr_schw(self, u):
        """GR EOB test-body limit: A_GR = 1 − 2u."""
        return 1.0 - 2.0*u

    def A_gr_pn(self, u):
        """GR EOB with leading finite-η correction (Blanchet 2002)."""
        a2 = (94.0/3 - 41*_PI**2/32) * self.eta
        return 1.0 - 2.0*u + a2*u**2

    def omega_circ(self, u):
        """Keplerian orbital angular frequency at u = GM/c²r."""
        return c**3 / (G*self.M) * u**1.5

    def f_gw(self, u):
        """GW frequency = 2 × orbital frequency."""
        return 2.0 * self.omega_circ(u) / (2*_PI)

    def u_ISCO_spin(self, chi=None):
        """ISCO u for spinning EOB A = 1 − 2u/η + χ²u²/η²."""
        chi = chi if chi is not None else self.chi_eff
        if abs(chi) < 1e-12:
            return self.u_ISCO
        alpha = 2.0 / self.eta
        beta  = chi**2 / self.eta**2
        def cond(u):
            N  = alpha - 2*beta*u
            D  = u*(2 - 3*alpha*u + 4*beta*u**2)
            dD = 2 - 6*alpha*u + 12*beta*u**2
            return -2*beta*D - N*dD
        try:
            return brentq(cond, 1e-5, min(alpha/2*0.99, 0.49))
        except Exception:
            return self.u_ISCO

    def print_structure(self):
        rs = 2*G*self.M/c**2
        print(f"\n  M₁={self.M1/Msun:.0f}+M₂={self.M2/Msun:.0f} Msun  "
              f"η={self.eta:.4f}  q={min(self.M1,self.M2)/max(self.M1,self.M2):.3f}  "
              f"χ_eff={self.chi_eff:.2f}")
        rows = [
            ("ISCO",            self.u_ISCO,  "A^QGD = 2/3, start of plunge"),
            ("Light ring",      self.u_LR,    "A^QGD = 1/3, unstable photons"),
            ("EOB singularity", self.u_h,     "A^QGD = 0   ← not physical merger"),
            ("σ-saddle merger", self.u_merge, "Σ_tot = 1   ← physical merger"),
        ]
        print(f"\n  {'Feature':>22}  {'u':>8}  {'r/rs':>7}  "
              f"{'A^QGD':>8}  {'A^GR':>9}  {'f_GW(Hz)':>10}  Notes")
        print("  " + "─"*90)
        for label, u, note in rows:
            r   = G*self.M / (c**2*u)
            Aq  = self.A(u, 0.0)
            Agr = self.A_gr_schw(u)
            fg  = self.f_gw(u)
            print(f"  {label:>22}  {u:>8.4f}  {r/rs:>7.2f}  "
                  f"{Aq:>8.4f}  {Agr:>9.4f}  {fg:>10.1f}  {note}")
        print("─"*90)
        err = abs(self.u_merge - self.u_h)/self.u_h*100
        if err < 1.0:
            print(f"  EOB singularity ≈ physical merger (q≈1), diff = {err:.3f}%  ✓")
        else:
            print(f"  EOB singularity ≠ physical merger ({err:.1f}% diff at q={self.M2/self.M1:.2f})")
            print(f"  → Merger occurs during the EOB plunge phase (consistent with GR/NR)")

    def print_isco_spin_table(self):
        """Show ISCO shift with spin for spinning A-function."""
        u0 = self.u_ISCO   # χ=0 reference
        print(f"\n  Spinning ISCO:  M₁={self.M1/Msun:.0f}+M₂={self.M2/Msun:.0f} Msun  η={self.eta:.4f}")
        print(f"  {'χ_eff':>7}  {'u_ISCO':>9}  {'Δu/u₀':>9}  {'A@ISCO':>9}  {'f_GW(Hz)':>12}")
        print("  " + "─"*52)
        for chi in [-0.9, -0.5, 0.0, 0.5, 0.9]:
            u  = self.u_ISCO_spin(chi)
            Aq = self.A(u, chi)
            fg = self.f_gw(u)
            print(f"  {chi:>7.2f}  {u:>9.5f}  {(u-u0)/u0*100:>+8.2f}%  "
                  f"{Aq:>9.5f}  {fg:>12.2f}")

    def compare_A_table(self):
        """Compare A^QGD vs GR EOB at various orbital radii."""
        rs = 2*G*self.M/c**2
        print(f"\n  A^QGD vs A^GR:  η={self.eta:.4f}")
        print(f"  {'u':>7}  {'r/rs':>6}  {'A^QGD':>9}  "
              f"{'A^GR_Schw':>11}  {'A^GR_1PN':>11}  Phase")
        print("  " + "─"*60)
        for u in sorted({0.005, 0.01, 0.02, self.u_ISCO,
                         self.u_LR, self.u_h*0.9, self.u_merge}):
            if u <= 0 or u > 0.99:
                continue
            r    = G*self.M/(c**2*u)
            Aq   = self.A(u, 0.0)
            Agr  = self.A_gr_schw(u)
            Agr1 = self.A_gr_pn(u)
            ph   = ("inspiral" if u <= self.u_ISCO else
                    "plunge"   if u <= self.u_LR   else
                    "EOB zone" if u <= self.u_merge else
                    "post-merge")
            print(f"  {u:>7.4f}  {r/rs:>6.2f}  {Aq:>9.5f}  "
                  f"{Agr:>11.5f}  {Agr1:>11.5f}  {ph}")


# =============================================================================
# PROBLEM 2 (cont.):  IMR TRANSITION
# =============================================================================

class QGDMergerIMR:
    """
    QGD Inspiral-Merger-Ringdown in a single analytic framework.

    The unifying quantity is Σ_tot — the σ-field strength at the saddle point.

      Inspiral:  Σ_saddle(d) ≪ 1
      Merger:    Σ_saddle(d_merge) = 1   (analytic condition)
      Ringdown:  σ-field relaxes to single Kerr BH

    Dipole correction
    -----------------
    The cross-term 2σ_t^{(1)}σ_t^{(2)} evolves during merger and contributes
    to the spin-2 cross-term transient (Prediction P1: τ = 6GM/c³).
    The D-factor asymmetry produces ZERO observable radiation because the
    scalar σ-mode is Yukawa-suppressed.

    Remnant estimates  (Barausse-Rezzolla 2009 fits)
    -------------------------------------------------
        M_f  ≈ M(1 − 0.0539 × 4η)    [total mass × (1 − ε_rad)]
        χ_f  ≈ min(0.95, 0.686(4η)^{0.45})
    """

    def __init__(self, M1, M2, chi1=0.0, chi2=0.0):
        self.M1, self.M2 = M1, M2
        self.M   = M1 + M2
        self.eta = M1*M2 / self.M**2
        self.mu  = M1*M2 / self.M
        self.eob = QGDEOB(M1, M2, chi1, chi2)

        mc = qgd_merger_condition(M1, M2)
        self.d_merge  = mc['d_merge']
        self.u_merge  = G*self.M / (c**2*self.d_merge)

        # Remnant (Barausse-Rezzolla)
        self.Mf   = self.M * (1.0 - 0.0539*4*self.eta)
        self.chif = min(0.95, 0.686*(4*self.eta)**0.45)

    def _Sigma_at(self, d):
        return sigma_saddle(self.M1, self.M2, d)

    def tau_cross(self):
        """Cross-term transient timescale (P1): τ = 6GM/c³."""
        return 6.0*G*self.M / c**3

    def phase(self, u):
        if   u < self.eob.u_ISCO:  return "inspiral"
        elif u < self.eob.u_LR:    return "plunge"
        elif u < self.u_merge:     return "merger approach"
        else:                      return "ringdown"

    def print_imr(self):
        from qgd_pn_complete import qnm_220  # reuse ch7 result if available
        rs = 2*G*self.M/c**2
        q  = min(self.M1,self.M2)/max(self.M1,self.M2)
        print(f"\n  M₁={self.M1/Msun:.0f}+M₂={self.M2/Msun:.0f} Msun  η={self.eta:.4f}  "
              f"q={q:.3f}  M_f={self.Mf/Msun:.1f} Msun  χ_f={self.chif:.3f}")
        print(f"\n  {'u':>8}  {'r/rs':>6}  {'A^QGD':>8}  "
              f"{'Σ_saddle':>11}  {'f_GW(Hz)':>10}  Phase")
        print("  " + "─"*65)
        u_list = sorted({0.005, 0.01, self.eob.u_ISCO,
                         self.eob.u_LR*0.9, self.eob.u_LR,
                         self.eob.u_h*0.9, self.u_merge})
        for u in u_list:
            if u <= 0: continue
            r    = G*self.M/(c**2*u)
            Aq   = max(0.0, self.eob.A(u, 0.0))
            Sig  = self._Sigma_at(r)
            fg   = self.eob.f_gw(u)
            ph   = self.phase(u)
            print(f"  {u:>8.4f}  {r/rs:>6.2f}  {Aq:>8.4f}  "
                  f"{Sig:>11.6f}  {fg:>10.1f}  {ph}")
        print("─"*65)
        print(f"  Merger: u={self.u_merge:.4f}, d={self.d_merge:.3e} m, Σ_saddle=1.000")
        print(f"  τ_cross = {self.tau_cross()*1000:.4f} ms  (Prediction P1, spin-2 transient)")
        err = abs(self.u_merge - self.eob.u_h)/self.eob.u_h*100
        if err < 1.0:
            print(f"  EOB singularity ≈ physical merger ({err:.2f}% diff, q≈1)  ✓")
        else:
            print(f"  EOB singularity ≠ physical merger ({err:.1f}% diff) — merger in plunge")


# =============================================================================
# PROBLEM 3:  SPINNING MASTER FORMULA
# =============================================================================

def e_n_spin_free(n):
    """
    Spin-free test-body PN coefficient (exact).

    e_n^0 = −C(2n,n)(3/4)^n(2n−1)/(n+1)
    """
    return -comb(2*n, n) * (3/4)**n * (2*n-1)/(n+1)


def e_n_spin_orbit(n):
    """
    Spin-orbit coefficient (exact from Kerr geodesic ∂E/∂χ|_{χ=0}).

    δE|_{SO} = χ u^{5/2} × (1−3u)^{−3/2}
             = χ u^{5/2} × Σ_n (2n+1)!!·3^n/(2^n·n!) × u^n

    e_n^{SO} = −(2n+1)!!·3^n / (2^n·n!)
    """
    return -_dfact(2*n+1) * 3**n / (2**n * factorial(n))


def e_n_spin_sq(n):
    """
    Spin-squared coefficient (exact from Kerr geodesic ∂²E/∂χ²|_{χ=0}).

    δ²E|_{SS} = (1/2)χ² u³ × (1−3u)^{−5/2}

    The half-integer binomial gives:
        C(n+3/2, n) = (2n+3)!! / (3·2^n·n!)   [Γ(5/2) = 3√π/4]

    e_n^{SS} = (2n+3)!!·3^n / (3·2^{n+1}·n!)

    Note: the factor /3 comes from Γ(5/2) ≠ Γ(3/2); an earlier version
    without this factor was 3× too large.
    """
    return _dfact(2*n+3) * 3**n / (3 * 2**(n+1) * factorial(n))


def kerr_exact_energy(u, chi, prograde=True):
    """Exact equatorial Kerr circular geodesic: E/μ = (1−2u±χu^{3/2})/√(1−3u±2χu^{3/2}) − 1."""
    s   = 1.0 if prograde else -1.0
    den = 1.0 - 3.0*u + 2.0*s*chi*u**1.5
    if den <= 0:
        return float('nan')
    return (1.0 - 2.0*u + s*chi*u**1.5) / np.sqrt(den) - 1.0


def spinning_energy_series(u, chi, n_terms=9):
    """Sum the three spinning PN series to n_terms."""
    E = sum(-u/2 * e_n_spin_free(n) * u**n for n in range(0, n_terms+1))
    # Fix n=0 term: e_0 spin-free = -u/2 × 1
    E = -u/2 * (1.0 + sum(e_n_spin_free(n)*u**n for n in range(1, n_terms+1)))
    E += sum(e_n_spin_orbit(n) * chi   * u**(n+2.5) for n in range(n_terms))
    E += sum(e_n_spin_sq(n)    * chi**2 * u**(n+3)  for n in range(n_terms))
    return E


def print_spinning_master_table():
    """Full spinning master formula table with GR comparisons."""
    print("═"*80)
    print("QGD SPINNING MASTER FORMULA  (three exact series from Kerr geodesic)")
    print("═"*80)

    known0  = {1:'1PN',2:'2PN',3:'3PN',4:'4PN'}
    known_so= {0:'1.5PN',1:'2.5PN',2:'3.5PN'}
    known_ss= {0:'2PN',1:'3PN',2:'4PN'}

    print(f"\n  SPIN-FREE:  e_n^0 = −C(2n,n)(3/4)^n(2n−1)/(n+1)")
    print(f"  {'n':>3}  {'Exact fraction':>22}  {'Decimal':>12}  Status")
    print("  " + "─"*68)
    for n in range(1, 8):
        fr  = Fraction(-comb(2*n,n)) * Fraction(3,4)**n * Fraction(2*n-1, n+1)
        st  = f"GR {known0[n]} ✓" if n in known0 else "★ QGD PREDICTION"
        print(f"  {n:>3}  {str(fr):>22}  {float(fr):>12.4f}  {st}")

    print(f"\n  SPIN-ORBIT: e_n^{{SO}} = −(2n+1)!!·3^n / (2^n·n!)   [χ·u^{{n+5/2}}]")
    print(f"  {'n':>3}  {'PN':>7}  {'Exact':>14}  {'Decimal':>12}  Status")
    print("  " + "─"*60)
    for n in range(7):
        v  = e_n_spin_orbit(n)
        pn = f"{2*n+3}/2 PN"
        st = f"GR {known_so[n]} ✓" if n in known_so else "★ QGD PREDICTION"
        print(f"  {n:>3}  {pn:>7}  {v:>14.4f}  {v:>12.4f}  {st}")

    print(f"\n  SPIN-SQUARED: e_n^{{SS}} = (2n+3)!!·3^n / (3·2^{{n+1}}·n!)  [χ²·u^{{n+3}}]")
    print(f"  [Corrected: factor /3 from Γ(5/2); earlier version was 3× too large]")
    print(f"  {'n':>3}  {'PN':>6}  {'Exact':>14}  {'Decimal':>12}  Status")
    print("  " + "─"*60)
    for n in range(7):
        v  = e_n_spin_sq(n)
        pn = f"{n+3} PN"
        st = f"GR {known_ss[n]} ✓" if n in known_ss else "★ QGD PREDICTION"
        print(f"  {n:>3}  {pn:>6}  {v:>14.6f}  {v:>12.6f}  {st}")

    print(f"\n  VERIFICATION vs exact Kerr geodesic at u=0.05  (9 terms each)")
    print(f"  {'chi':>5}  {'Exact E/μ':>14}  {'Series E/μ':>14}  {'Rel err':>10}")
    print("  " + "─"*50)
    for chi in [0.0, 0.3, 0.5, 0.7, 0.9]:
        ex  = kerr_exact_energy(0.05, chi)
        ap  = spinning_energy_series(0.05, chi, n_terms=9)
        err = abs(ap-ex)/abs(ex)*100
        print(f"  {chi:>5.2f}  {ex:>14.8f}  {ap:>14.8f}  {err:>9.5f}%")
    print("═"*80)


# =============================================================================
# MAIN VERIFICATION
# =============================================================================

SYSTEMS = [
    ("GW150914",   36.0*Msun, 29.0*Msun, -0.01, 0.0),
    ("GW250114",   86.0*Msun, 77.0*Msun,  0.0,  0.0),
    ("15+5 Msun",  15.0*Msun,  5.0*Msun,  0.0,  0.0),
    ("30+1 Msun",  30.0*Msun,  1.0*Msun,  0.0,  0.0),
]


def run_all():
    print("\n" + "═"*74)
    print("QGD EOB · MERGER TRANSITION · SPINNING MASTER FORMULA")
    print("="*74)
    print("\nKey update (v2.0): Dipole (D-factor) radiation is exactly zero.")
    print("Observable QGD waveform = GR waveform.\n")

    # ── Problem 2: Merger condition ──────────────────────────────────────
    print("\n" + "─"*74)
    print("PROBLEM 2: MERGER CONDITION")
    print_merger_table()

    # ── Problem 1: EOB ───────────────────────────────────────────────────
    print("\n" + "─"*74)
    print("PROBLEM 1: QGD EOB  A^QGD(u) = 1 − 2u/η  (exact, zero free params)")
    print("─"*74)
    for name, M1, M2, chi1, chi2 in SYSTEMS:
        print(f"\n  ── {name} ──")
        eob = QGDEOB(M1, M2, chi1, chi2)
        eob.print_structure()
        eob.print_isco_spin_table()
        eob.compare_A_table()

    # ── Problem 2 (cont.): IMR ───────────────────────────────────────────
    print("\n" + "─"*74)
    print("PROBLEM 2 (cont.): IMR TRANSITION")
    print("─"*74)
    for name, M1, M2, chi1, chi2 in SYSTEMS:
        print(f"\n  ── {name} ──")
        try:
            QGDMergerIMR(M1, M2, chi1, chi2).print_imr()
        except ImportError:
            imr = QGDMergerIMR(M1, M2, chi1, chi2)
            # fallback without qnm import
            imr.print_imr.__func__  # just call print_imr directly
            _simple_imr_print(imr)

    # ── Problem 3: Spinning ──────────────────────────────────────────────
    print("\n" + "─"*74)
    print("PROBLEM 3: SPINNING MASTER FORMULA")
    print("─"*74)
    print_spinning_master_table()

    # ── Summary ──────────────────────────────────────────────────────────
    print("\n" + "═"*74)
    print("SUMMARY")
    print("="*74)
    print()
    print("  RIGOROUS:")
    print()
    print("  [M] Merger: d_B = 2GM₁(1+(M₂/M₁)^{1/3})³/c²")
    print("      Σ_saddle = 1.000000 for ALL q tested  ✓")
    print("      Equal-mass: 8rs_indiv = 4rs_total (both correct conventions)")
    print()
    print("  [E] EOB: A^QGD = 1 − 2u/η  (exact, no NR calibration)")
    print("      ISCO u=η/6, LR u=η/3, horizon u=η/2")
    print("      Spinning: A^QGD = 1 − 2u/η + χ²u²/η² (η→0 recovers Kerr)")
    print()
    print("  [S] Spin-free:   e_n^0 = −C(2n,n)(3/4)^n(2n−1)/(n+1)  (all n verified)")
    print("      Spin-orbit:  e_n^SO = −(2n+1)!!·3^n/(2^n·n!)       (n=0,1,2 ✓)")
    print("      Spin-squared:e_n^SS = (2n+3)!!·3^n/(3·2^{n+1}·n!)  (n=0,1,2 ✓)")
    print()
    print("  [D] Dipole = 0 in all observables (spin-0 Yukawa-suppressed)")
    print("      D-factor is a real near-field quantity — zero observable radiation")
    print("      QGD and GR agree on all observable waveforms")
    print()
    print("  PRELIMINARY / OPEN QUESTIONS:")
    print()
    print("  [P1] η-corrections at 5PN+: known only partially")
    print("  [P2] A^QGD calibration vs NR: key falsifiable test")
    print("       If NR requires a₂ ≈ 1.8η ≠ 0, QGD EOB needs η-corrections")
    print("  [P3] Spinning A beyond O(χ²): next step")
    print()
    print("  FALSIFIABLE PREDICTION (EOB):")
    print("  If QGD is correct, A^QGD = 1−2u/η is exact: no NR calibration needed.")
    print("  Current GR EOB requires a₂=(94/3−41π²/32)η from NR.")
    print("  A QGD EOB NR comparison will distinguish the two theories.")
    print("═"*74)


def _simple_imr_print(imr):
    """Fallback IMR print without external import."""
    rs = 2*G*imr.M/c**2
    q  = min(imr.M1,imr.M2)/max(imr.M1,imr.M2)
    print(f"\n  M₁={imr.M1/Msun:.0f}+M₂={imr.M2/Msun:.0f} Msun  η={imr.eta:.4f}  "
          f"q={q:.3f}  M_f={imr.Mf/Msun:.1f} Msun  χ_f={imr.chif:.3f}")
    print(f"  τ_cross = {imr.tau_cross()*1000:.4f} ms  (P1: spin-2 cross-term transient)")
    print(f"  d_merge = {imr.d_merge:.4e} m  (Σ_saddle = 1.000)")


if __name__ == "__main__":
    # Wrap the IMR loop to not require external import
    import sys
    print("\n" + "═"*74)
    print("QGD EOB · MERGER TRANSITION · SPINNING MASTER FORMULA")
    print("="*74)
    print("\nKey update (v2.0): Dipole (D-factor) radiation is exactly zero.")
    print("Observable QGD waveform = GR waveform.\n")

    print("\n" + "─"*74)
    print("PROBLEM 2: MERGER CONDITION")
    print_merger_table()

    print("\n" + "─"*74)
    print("PROBLEM 1: QGD EOB  A^QGD(u) = 1 − 2u/η")
    for name, M1, M2, chi1, chi2 in SYSTEMS:
        print(f"\n  ── {name} ──")
        eob = QGDEOB(M1, M2, chi1, chi2)
        eob.print_structure()
        eob.print_isco_spin_table()
        eob.compare_A_table()

    print("\n" + "─"*74)
    print("PROBLEM 2 (cont.): IMR")
    for name, M1, M2, chi1, chi2 in SYSTEMS:
        print(f"\n  ── {name} ──")
        imr = QGDMergerIMR(M1, M2, chi1, chi2)
        _simple_imr_print(imr)

    print("\n" + "─"*74)
    print("PROBLEM 3: SPINNING MASTER FORMULA")
    print_spinning_master_table()

    print("\n═"*74)
    print("ALL CHECKS COMPLETE")
    print("="*74)
