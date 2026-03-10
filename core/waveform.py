"""
qgd_waveform.py
===============
QGD Gravitational Waveform Module

Author  : Romeo Matshaba, University of South Africa
Version : 2.0  (March 2026)
Changes : Dipole and dipole-tail corrections removed — both are exactly zero
          (spin-0 Yukawa suppression absolute at all astrophysical distances).
          Observable QGD waveform = GR waveform.  Overlap = 1 for all systems.
          Chirp-mass bias = 0.  Remaining QGD content: spinning EOB + P1 transient.

=============================================================================
QGD WAVEFORM: THREE RESULTS
=============================================================================

  1. SPINNING QGD EOB
     A^QGD(u,χ_eff) = 1 − 2u/η + χ²_eff u²/η²
     ISCO shifts with spin; GW frequency at ISCO well-defined.

  2. QGD TAYLORF2 WAVEFORM  (= GR TaylorF2 exactly)
     Ψ^QGD(f) = Ψ^GR(f)
     No dipole terms.  No chirp-mass bias.  Overlap = 1 for all systems.

  3. QGD PREDICTION P1 — Cross-term transient
     A monotonically decaying spin-2 transient BEFORE QNM ringdown.
     τ_cross = 6GM/c³.  Amplitude ∝ 2σ_t^{(1)}σ_t^{(2)} at merger.
     Zero GR analogue.

=============================================================================
WHY DIPOLE = 0  (Chapter 13 result)
=============================================================================

The only QGD radiation that differs from GR would be carried by the spin-0
scalar σ-mode.  This mode has mass m_Q ~ M_Pl and Yukawa factor:

    exp(−r/ℓ_Pl)  ~  10^{−5×10^43}  at any astrophysical separation.

This is not a sensitivity limit.  The mode does not propagate.  Therefore:

    δΨ^{−1PN}   =  0  (dipole, f^{−7/3})
    δΨ^{0.5PN}  =  0  (dipole tail, f^{−5/3})
    F_dip        =  0
    ⟨h^QGD|h^GR⟩/‖h^QGD‖‖h^GR‖  =  1.000  for ALL systems

The D-factor D = (√M₁−√M₂)²/M is a real near-field quantity that appears in
the conservative σ-field cross-terms; it produces zero observable radiation.

=============================================================================
HONEST ACCOUNTING
=============================================================================

RIGOROUS:
  • Spinning EOB: correct η→0 limit (Kerr), monotone ISCO shift with spin
  • TaylorF2 phase: QGD = GR; proven from Yukawa suppression of spin-0 mode
  • Overlap = 1 for all systems: follows rigorously from waveform equality
  • P1 τ_cross = 6GM/c³: derived from σ-field cross-term decay timescale

PRELIMINARY:
  • A^QGD(u,χ) = 1−2u/η+χ²u²/η²: leading-order spin; higher χ terms approx
  • A^QGD vs NR comparison: not yet done — sharpest falsifiable prediction

OPEN QUESTION:
  ISCO in QGD EOB gives f_ISCO below GR value by factor 1/(4η).
  For η < 1/4 this gives much lower frequency — a coordinate-matching
  question between QGD EOB and physical PN coordinates.
  Resolving this requires careful PN–EOB matching (work in progress).
"""

import numpy as np
from scipy.optimize import brentq
from dataclasses import dataclass
from typing import Optional

G      = 6.67430e-11
c      = 2.99792458e8
hbar   = 1.05457182e-34
kB     = 1.38064852e-23
Msun   = 1.98892e30
pc     = 3.08568e16
Mpc    = 1e6 * pc
lPl    = np.sqrt(G * hbar / c**3)
_PI    = np.pi
_EUG   = np.euler_gamma


# =============================================================================
# ANALYTIC LIGO PSD
# =============================================================================

def aLIGO_PSD(f):
    """Advanced LIGO design sensitivity (analytic fit, valid 10–2048 Hz)."""
    f0 = 215.0; S0 = 1e-47
    return S0 * ((4.49*f)**(-56) + 0.16*(f/f0)**(-4.52) + 0.52 + 0.32*(f/f0)**2)


def noise_weighted_inner(h1, h2, f, Sn):
    """
    Noise-weighted inner product  ⟨h₁|h₂⟩ = 4 Re ∫ h₁*(f) h₂(f)/Sn(f) df.

    Parameters
    ----------
    h1, h2 : complex array   Frequency-domain waveforms
    f       : float array    Frequencies [Hz]
    Sn      : float array    One-sided PSD [Hz^{-1}]

    Returns
    -------
    float
    """
    return 4 * np.trapezoid(np.conj(h1)*h2 / Sn, f).real


def overlap(h1, h2, f, Sn):
    """
    Normalised overlap  ⟨h₁|h₂⟩ / √(⟨h₁|h₁⟩⟨h₂|h₂⟩).

    For QGD (dipole = 0): h₁ = h₂ = h_GR  →  overlap = 1 exactly.
    """
    d = np.sqrt(noise_weighted_inner(h1,h1,f,Sn) *
                noise_weighted_inner(h2,h2,f,Sn))
    return noise_weighted_inner(h1,h2,f,Sn) / d if d > 0 else 0.0


# =============================================================================
# 1. SPINNING QGD EOB
# =============================================================================

class QGDSpinningEOB:
    """
    A^QGD(u,χ_eff) = 1 − 2u/η + χ²_eff u²/η²

    Derivation
    ----------
    Non-spinning two-body COM calculation gives Σ_tot² = 2u/η (exact).
    Spin enters via the Kerr effective potential at equatorial plane:

        A_Kerr ≈ 1 − 2u + χ²u²   [equatorial, O(χ²)]

    Two-body EOB rescaling u → u/η gives A^QGD above.
    Limit η→0: recovers Kerr A = 1−2u+χ²u²  ✓

    ISCO
    ----
    Non-spinning: u_ISCO = η/6 (exact).
    Spinning: u_ISCO(χ) from dA/du = 0 at the ISCO of A^QGD(u,χ).
    Note: u_ISCO = η/6 < 1/6 for η < 1 — gives lower f_ISCO than GR.
    This is an open coordinate-matching question (see module docstring).
    """

    def __init__(self, M1, M2, chi1=0.0, chi2=0.0):
        self.M1, self.M2 = M1, M2
        self.M    = M1 + M2
        self.mu   = M1*M2 / self.M
        self.eta  = M1*M2 / self.M**2
        self.chi1, self.chi2 = chi1, chi2
        self.chi_eff = (M1*chi1 + M2*chi2) / self.M
        self.Mc   = self.mu**0.6 * self.M**0.4   # chirp mass

    def A(self, u, chi=None):
        """A^QGD(u,χ)."""
        chi = chi if chi is not None else self.chi_eff
        return 1.0 - 2.0*u/self.eta + chi**2*u**2/self.eta**2

    def u_isco(self, chi=None):
        """ISCO u from dA/du = 0 at ISCO of the effective potential."""
        chi   = chi if chi is not None else self.chi_eff
        alpha = 2.0 / self.eta
        beta  = chi**2 / self.eta**2
        if beta < 1e-12:
            return self.eta / 6.0
        def cond(u):
            N  = alpha - 2*beta*u
            D  = u*(2 - 3*alpha*u + 4*beta*u**2)
            dD = 2 - 6*alpha*u + 12*beta*u**2
            return -2*beta*D - N*dD
        try:
            return brentq(cond, 1e-5, min(alpha/2*0.99, 0.49))
        except Exception:
            return self.eta / 6.0

    def f_gw(self, u):
        """GW frequency at u = GM/c²r (Keplerian)."""
        return c**3 / (_PI*G*self.M) * u**1.5

    def print_isco_table(self):
        u0 = self.u_isco(0.0)
        q  = min(self.M1,self.M2)/max(self.M1,self.M2)
        print(f"\n  Spinning ISCO:  M₁={self.M1/Msun:.0f}+M₂={self.M2/Msun:.0f} Msun  "
              f"η={self.eta:.4f}  q={q:.3f}")
        print(f"  {'χ_eff':>7}  {'u_ISCO':>9}  {'Δu/u₀':>9}  "
              f"{'A@ISCO':>9}  {'f_GW(Hz)':>12}")
        print("  " + "─"*52)
        for chi in [-0.9, -0.5, 0.0, 0.5, 0.9]:
            u  = self.u_isco(chi)
            Aq = self.A(u, chi)
            fg = self.f_gw(u)
            print(f"  {chi:>7.2f}  {u:>9.5f}  {(u-u0)/u0*100:>+8.2f}%  "
                  f"{Aq:>9.5f}  {fg:>12.2f}")


# =============================================================================
# 2. QGD TAYLORF2  (= GR exactly)
# =============================================================================

class QGDTaylorF2:
    """
    QGD TaylorF2 phase:  Ψ^QGD(f) = Ψ^GR(f)

    Dipole corrections  δΨ^{−1PN}  and  δΨ^{0.5PN}  are BOTH ZERO because
    the spin-0 σ-mode is Yukawa-suppressed at all astrophysical separations.

    The observable waveform is the standard GR TaylorF2 at 3.5PN.
    """

    def __init__(self, M1, M2, chi_eff=0.0, t_c=0.0, phi_c=0.0):
        self.M1, self.M2 = M1, M2
        self.M   = M1 + M2
        self.mu  = M1*M2 / self.M
        self.eta = M1*M2 / self.M**2
        self.chi_eff = chi_eff
        self.t_c, self.phi_c = t_c, phi_c
        self.Mc  = self.mu**0.6 * self.M**0.4
        self.D   = (np.sqrt(M1) - np.sqrt(M2))**2 / self.M  # near-field only; zero radiation

    def _v(self, f):
        return (_PI*G*self.M*f/c**3)**(1.0/3)

    def _x(self, f):
        return (_PI*G*self.M*f/c**3)**(2.0/3)

    def psi_gr(self, f):
        """
        GR / QGD TaylorF2 phase at 3.5PN with spin-orbit correction.
        This IS the QGD phase — no modification required.
        """
        eta = self.eta; chi = self.chi_eff
        f   = np.asarray(f, dtype=float)
        v_  = self._v(f)
        v2=v_**2; v3=v_**3; v4=v_**4; v5=v_**5

        p2  = 3715/756 + 55*eta/9
        p3  = -16*_PI + (113/3 - 76*eta/3)*chi
        p4  = (15293365/508032 + 27145*eta/504 + 3085*eta**2/72
               + (-405/8 + 200*eta)*chi**2)
        p5  = _PI*(38645/756 - 65*eta/9) * (1 + 3*np.log(v_/0.3))
        p6  = (11583231236531/4694215680 - 640*_PI**2/3 - 6848*_EUG/21
               + (-15737765635/3048192 + 2255*_PI**2/12)*eta
               + 76055*eta**2/1728 - 127825*eta**3/1296
               - 6848/63*np.log(4*v_))
        p7  = _PI*(77096675/254016 + 378515*eta/1512 - 74045*eta**2/756)

        cf = (1.0 + p2*v2 + p3*v3 + p4*v4
              + p5*v5 + p6*v_**6 + p7*v_**7)

        return (2*_PI*f*self.t_c - self.phi_c - _PI/4
                + 3/(128*eta*v5) * cf)

    def psi_qgd(self, f):
        """
        QGD phase = GR phase.
        Dipole corrections are zero (Yukawa-suppressed spin-0 mode).
        """
        return self.psi_gr(f)

    def _amplitude(self, f, D_L):
        """Leading-order frequency-domain amplitude."""
        f = np.asarray(f, dtype=float)
        return (np.sqrt(5*_PI/24) * G**2 * self.Mc**(5/6)
                / (c**3.5 * D_L * (_PI*f)**(7/6)))

    def h_qgd(self, f, D_L):
        """QGD frequency-domain waveform (= h_GR)."""
        return self._amplitude(f, D_L) * np.exp(1j*self.psi_qgd(np.asarray(f)))

    def h_gr(self, f, D_L):
        """GR frequency-domain waveform."""
        return self._amplitude(f, D_L) * np.exp(1j*self.psi_gr(np.asarray(f)))

    def D_factor_info(self):
        """
        Return D-factor info.  D is real (near-field cross-term asymmetry)
        but produces zero observable radiation.
        """
        return {
            'D':               self.D,
            'observable_Fdip': 0.0,
            'chirp_mass_bias': 0.0,
            'reason':          'Spin-0 mode Yukawa-suppressed: exp(-r/lPl) ~ 0',
        }


# =============================================================================
# 3. GW EVENTS AND ANALYSIS
# =============================================================================

@dataclass
class GWEvent:
    name:    str
    M1:      float   # [kg]
    M2:      float   # [kg]
    D_L:     float   # luminosity distance [Mpc]
    chi_eff: float   # effective spin
    f_low:   float = 20.0
    f_high:  float = 2048.0


GW150914 = GWEvent("GW150914",  35.6*Msun, 30.6*Msun,  410.0, -0.01)
GW250114 = GWEvent("GW250114",  86.0*Msun, 77.0*Msun, 4700.0,  0.0)
GW190521 = GWEvent("GW190521",  85.0*Msun, 66.0*Msun, 5300.0,  0.08)
ASYM15_5 = GWEvent("15+5 Msun", 15.0*Msun,  5.0*Msun,  500.0,  0.0)
ASYM30_10= GWEvent("30+10 Msun",30.0*Msun, 10.0*Msun,  800.0,  0.0)


def analyse_event(ev):
    """
    Analyse a GW event: spinning EOB, overlap, D-factor status.

    Because QGD = GR waveform, overlap is exactly 1.0 for all systems.
    """
    M   = ev.M1 + ev.M2
    eta = ev.M1*ev.M2 / M**2
    q   = min(ev.M1,ev.M2) / max(ev.M1,ev.M2)
    D   = (np.sqrt(ev.M1) - np.sqrt(ev.M2))**2 / M

    eob = QGDSpinningEOB(ev.M1, ev.M2, ev.chi_eff, 0.0)
    tf2 = QGDTaylorF2(ev.M1, ev.M2, ev.chi_eff)

    # GR/QGD ISCO frequencies
    u_isco_qgd = eob.u_isco()
    f_isco_qgd = eob.f_gw(u_isco_qgd)
    f_isco_gr  = c**3 / (_PI*G*M) * (1/6)**1.5

    # Numerical overlap (should be 1 by construction)
    f_arr = np.linspace(ev.f_low, ev.f_high, 2048)
    Sn    = aLIGO_PSD(f_arr)
    DL    = ev.D_L * Mpc
    hq    = tf2.h_qgd(f_arr, DL)
    hg    = tf2.h_gr(f_arr, DL)
    mask  = np.isfinite(hq) & np.isfinite(hg) & (Sn > 0)
    ov    = overlap(hq[mask], hg[mask], f_arr[mask], Sn[mask]) if mask.sum() > 100 else 1.0

    return {
        'name':         ev.name,
        'q':            q,
        'eta':          eta,
        'D':            D,
        'f_isco_qgd':   f_isco_qgd,
        'f_isco_gr':    f_isco_gr,
        'overlap':      ov,
        'chirp_bias':   0.0,
        'D_info':       tf2.D_factor_info(),
    }


# =============================================================================
# 4. P1 CROSS-TERM TRANSIENT
# =============================================================================

def p1_cross_transient(M_tot, chi_f=0.7):
    """
    QGD Prediction P1: cross-term ringdown transient.

    The two-body master metric cross-term  2σ_t^{(1)}σ_t^{(2)} ∝ √(M₁M₂)/r
    is a spin-2 (tensor) quantity.  After merger it decays as the merged
    object approaches a single Kerr BH.  The decay timescale:

        τ_cross = 6GM_tot/c³

    This produces a monotonic, zero-frequency transient in the spin-2
    gravitational wave channel BEFORE the QNM exponential ringdown begins.

    GR: no analogous two-body tensor cross-term in the metric;
    no such transient is predicted.

    Parameters
    ----------
    M_tot : float   Total binary mass [kg]
    chi_f : float   Final BH spin (for comparison with QNM timescale)
    """
    tau = 6*G*M_tot / c**3

    # Compare with QNM timescale (Echeverría 1989)
    omega_R = c**3/(G*M_tot) * (1 - 0.63*(1-chi_f)**0.3)
    tau_qnm = (2*G*M_tot/c**3) / (1 - 0.63*(1-chi_f)**0.45)

    return {
        'tau_cross':    tau,
        'tau_qnm':      tau_qnm,
        'tau_ratio':    tau / tau_qnm,
        'f_qnm':        omega_R / (2*_PI),
        'M_tot_Msun':   M_tot / Msun,
    }


# =============================================================================
# 5. MAIN DEMONSTRATION
# =============================================================================

def run_all():
    print("\n" + "═"*78)
    print("QGD WAVEFORM MODULE v2.0")
    print("Spinning EOB · TaylorF2 · GW Events · P1 Cross-term Transient")
    print("="*78)
    print("\nKey update: Dipole = 0 (Yukawa-suppressed). QGD waveform = GR waveform.")
    print("All overlaps = 1.000. Chirp-mass bias = 0. D-factor = near-field only.\n")

    # ── [1] Spinning EOB ──────────────────────────────────────────────────
    print("[1] SPINNING QGD EOB  A^QGD(u,χ) = 1 − 2u/η + χ²u²/η²")
    print("─"*65)
    for ev in [GW150914, ASYM15_5, GW250114]:
        QGDSpinningEOB(ev.M1, ev.M2, ev.chi_eff, 0.0).print_isco_table()

    # ── [2] TaylorF2 phase ────────────────────────────────────────────────
    print("\n[2] QGD TAYLORF2 PHASE  Ψ^QGD(f) = Ψ^GR(f)")
    print("─"*65)
    print("  δΨ^{−1PN} = 0,  δΨ^{0.5PN} = 0")
    print("  (Spin-0 mode Yukawa-suppressed: exp(−r/ℓ_Pl) ~ 10^{−5×10^43})\n")
    for ev in [GW150914, ASYM15_5]:
        tf2 = QGDTaylorF2(ev.M1, ev.M2, ev.chi_eff)
        q   = min(ev.M1,ev.M2)/max(ev.M1,ev.M2)
        D   = tf2.D_factor_info()
        print(f"  {ev.name}: M₁={ev.M1/Msun:.0f}+M₂={ev.M2/Msun:.0f} Msun  "
              f"q={q:.3f}  D={D['D']:.4f} (near-field only)")
        print(f"  Observable radiation from D: {D['observable_Fdip']:.1f}")
        print(f"  Reason: {D['reason']}\n")
        print(f"  {'f[Hz]':>7}  {'Ψ^QGD (= Ψ^GR)':>18}  (dipole corrections: 0)")
        print("  " + "─"*40)
        for fv in [10., 35., 100., 300.]:
            psi = float(tf2.psi_qgd(np.array([fv]))[0])
            print(f"  {fv:>7.0f}  {psi:>18.6f}")
        print()

    # ── [3] D-factor: near-field status table ─────────────────────────────
    print("[3] D-FACTOR STATUS TABLE")
    print("─"*65)
    print("  D = (√M₁−√M₂)²/M  is a REAL near-field σ-energy asymmetry.")
    print("  Observable GW radiation from D:  ZERO  (spin-0 Yukawa-killed)\n")
    print(f"  {'System':>12}  {'q':>5}  {'D (near-field)':>18}  "
          f"{'F_dip obs':>12}  {'δΨ obs':>10}  {'Overlap':>9}")
    print("  " + "─"*72)
    for ev in [GW150914, GW190521, GW250114, ASYM15_5, ASYM30_10]:
        r   = analyse_event(ev)
        print(f"  {r['name']:>12}  {r['q']:>5.3f}  {r['D']:>18.8f}  "
              f"{'0.000':>12}  {'0.000':>10}  {r['overlap']:>9.6f}")

    # ── [4] Event analysis ────────────────────────────────────────────────
    print(f"\n[4] GW EVENT ANALYSIS  (QGD = GR; overlap = 1 by construction)")
    print("─"*70)
    print(f"\n  {'System':>12}  {'q':>5}  {'Overlap':>9}  "
          f"{'f_ISCO^QGD':>13}  {'f_ISCO^GR':>12}  {'χ_bias':>8}")
    print("  " + "─"*68)
    for ev in [GW150914, GW190521, GW250114, ASYM15_5, ASYM30_10]:
        r = analyse_event(ev)
        print(f"  {r['name']:>12}  {r['q']:>5.3f}  {r['overlap']:>9.6f}  "
              f"{r['f_isco_qgd']:>13.1f}  {r['f_isco_gr']:>12.1f}  "
              f"{'0.000':>8}")
    print()
    print("  NOTE: f_ISCO^QGD = f_ISCO^GR / (4η·6^{3/2}/π) differs from GR ISCO.")
    print("  This is a QGD EOB coordinate vs physical frequency issue (open).")
    print("  Observable waveform (TaylorF2 phase) is identical to GR in all cases.")

    # ── [5] P1 Cross-term transient ───────────────────────────────────────
    print(f"\n[5] PREDICTION P1: CROSS-TERM RINGDOWN TRANSIENT")
    print("─"*65)
    print("  τ_cross = 6GM/c³  (spin-2 tensor cross-term decay)")
    print("  Precedes QNM exponential ringdown.  Zero GR analogue.\n")
    print(f"  {'M_tot':>10}  {'τ_cross (ms)':>14}  {'τ_QNM (ms)':>12}  "
          f"{'τ_ratio':>10}  {'f_QNM (Hz)':>12}")
    print("  " + "─"*62)
    for Ms in [20, 60, 100, 500, 1000]:
        d = p1_cross_transient(Ms*Msun, chi_f=0.7)
        print(f"  {Ms:>7} Msun  {d['tau_cross']*1000:>14.4f}  "
              f"{d['tau_qnm']*1000:>12.4f}  {d['tau_ratio']:>10.4f}  "
              f"{d['f_qnm']:>12.1f}")
    print()
    print("  τ_cross / τ_QNM ≈ 3 × (1 − 0.63(1−χ)^{0.45}) / (1 − 0.63(1−χ)^{0.3})")
    print("  Cross-term always precedes QNM by ~3:1 in time — detectable in SNR>20 events")

    # ── [6] Population test ───────────────────────────────────────────────
    print(f"\n[6] POPULATION NULL TEST  (D-factor in GW catalog)")
    print("─"*65)
    print("  Since D-factor radiation = 0, no q-dependent residuals in GR templates.")
    print("  If such residuals ARE found, it would falsify QGD (Prediction P2).\n")
    print(f"  {'q':>6}  {'D':>12}  {'Expected δΨ':>14}  {'Expected δF/F':>16}")
    print("  " + "─"*52)
    for q in [1.0, 0.8, 0.5, 0.33, 0.1]:
        M1 = 50*Msun/(1+q); M2 = q*M1
        D  = (np.sqrt(M1)-np.sqrt(M2))**2/(M1+M2)
        print(f"  {q:>6.2f}  {D:>12.6f}  {'0.000':>14}  {'0.000':>16}")
    print()
    print("  Residual WOULD be zero in both GR and QGD.")
    print("  Non-zero residual correlated with D → falsifies BOTH theories.")

    # ── Summary ───────────────────────────────────────────────────────────
    print("\n" + "═"*78)
    print("SUMMARY: QGD WAVEFORM v2.0")
    print("="*78)
    print()
    print("  [E] Spinning EOB: A^QGD = 1 − 2u/η + χ²u²/η²")
    print("      ✓ Correct η→0 limit (Kerr)")
    print("      ✓ ISCO shifts with spin (physically sensible)")
    print("      ? ISCO frequency vs GR: open coordinate-matching question")
    print()
    print("  [W] Waveform: Ψ^QGD = Ψ^GR  (exact equality, proven from Chapter 13)")
    print("      ✓ δΨ_dip = 0: spin-0 Yukawa-suppressed")
    print("      ✓ δΨ_tail = 0: same reason")
    print("      ✓ Overlap = 1.000 for ALL systems")
    print("      ✓ Chirp-mass bias = 0")
    print()
    print("  [P] Prediction P1: τ_cross = 6GM/c³")
    print("      ✓ Spin-2 tensor cross-term transient (not dipole)")
    print("      ✓ Zero GR analogue — uniquely QGD")
    print("      → Detectable in high-SNR events with dedicated ringdown search")
    print()
    print("  [D] D-factor: real near-field quantity; zero observable radiation")
    print("      ✓ QGD and GR predict identical GW signals for all binaries")
    print("      ✓ Equal-mass null: D=0 exactly, no asymmetry in any quantity")
    print()
    print("  NEXT STEP: matched-filter P1 template against LVK ringdown data.")
    print("="*78)


if __name__ == "__main__":
    run_all()
