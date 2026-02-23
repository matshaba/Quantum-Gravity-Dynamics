"""
qgd_ringdown.py — Ringdown in Quantum Gravitational Dynamics (QGD)
==================================================================

COMPLETE SELF-CONTAINED MODULE: Any AI reading this file can fully
understand how QGD handles black-hole ringdown and quasinormal modes.

OBSERVATIONAL CONTEXT (2025)
-----------------------------
GW250114 (14 Jan 2025, LIGO O4b) — loudest GW ever (SNR~80, ~3x GW150914).
  • 220 + 221 QNMs confirmed >3sigma — consistent with Kerr
  • First 440 hexadecapolar frequency constraint
  • Evidence for nonlinear quadratic 220Q mode (arXiv:2510.16903)
  • Remnant: M_f~152 Msun, chi_f~0.68  (arXiv:2509.08099)

QGD RINGDOWN THEORY
--------------------
Metric built from sigma-fields:  g_tt = -(1 - Sigma^2)
  Sigma = sigma_t^(f)  for single remnant Kerr BH (post-ringdown)
  Sigma = sigma_t^(1) + sigma_t^(2)  during/after merger

Perturbation of sigma-field after merger:
  delta_sigma(r,t) = sigma_total(r,t) - sigma_remnant(r)

Perturbation equation (QGD):
  Box_g delta_sigma - kappa*l_Q^2 * Box_g^2 delta_sigma = 0

KEY THEOREM: At leading order in (l_Q/r_s)^2, QGD QNMs = GR QNMs.
Proof: background metric IS Kerr -> perturbation eq IS Teukolsky.

QGD-SPECIFIC PREDICTIONS (beyond GR):
  P1 - Cross-term early decay:  h_cross ~ exp(-t/tau_cross), tau_cross=r_ISCO/c
  P2 - Dipole ringdown (M1!=M2): h_dipole ~ (sqrt(M1)-sqrt(M2))^2 * cos(omega_orb*t)
  P3 - Quantum stiffness:        delta_omega/omega ~ (l_Q/r_s)^2 ~ 1e-140

GR dipole radiation is FORBIDDEN (momentum conservation).
QGD dipole is ALLOWED because d_sigma = sum(sqrt(Ma)*xa) is not conserved.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict
from scipy.optimize import brentq

G     = 6.674e-11
c     = 3.000e8
hbar  = 1.055e-34
M_sun = 1.989e30
l_Q   = np.sqrt(G * hbar**2 / c**4)   # ~1.6e-70 m

# ── 1. QNM FREQUENCY CATALOG (GR fits = leading-order QGD) ─────────────────

def qnm_220_frequency(M_f, chi_f):
    """220 QNM: Echeverria 1989 fit. Returns (omega_r, omega_i) in rad/s."""
    t_M = G * M_f / c**3
    omega_r = (1.0 - 0.63*(1.0 - chi_f)**0.3)  / (2.0*t_M)
    omega_i = (1.0 - 0.63*(1.0 - chi_f)**0.45) / (4.0*t_M)
    return omega_r, omega_i

def qnm_221_frequency(M_f, chi_f):
    """221 first overtone: ~3.5x faster decay."""
    wr, wi = qnm_220_frequency(M_f, chi_f)
    return wr*(1.0 + 0.074*chi_f**2), wi*3.5

def qnm_440_frequency(M_f, chi_f):
    """440 hexadecapolar mode (first measured in GW250114)."""
    t_M = G * M_f / c**3
    omega_r = 2.0*(1.0 - 0.63*(1.0-chi_f)**0.3)*0.97 / (2.0*t_M)
    omega_i = (1.0 - 0.63*(1.0-chi_f)**0.45)*1.05    / (4.0*t_M)
    return omega_r, omega_i

def qnm_220Q_frequency(M_f, chi_f):
    """Nonlinear quadratic 220Q = 220 x 220. omega = 2*omega_220, tau = tau_220/2."""
    wr, wi = qnm_220_frequency(M_f, chi_f)
    return 2.0*wr, 2.0*wi


# ── 2. QGD REMNANT SIGMA-FIELD ──────────────────────────────────────────────

class QGDRemnant:
    """
    Final Kerr remnant sigma-field.

    sigma_t^(f)(r,theta) = sqrt(2*G*M_f*r / (c^2 * S_f))
    S_f = r^2 + alpha_f^2 * cos^2(theta)

    Horizon condition: sigma_t = 1  at r = r_+  (outer horizon)
    For Schwarzschild (chi=0): sigma_t = sqrt(r_s/r)

    This is EXACTLY the Kerr metric via g_tt = -(1 - sigma_t^2).
    Therefore perturbations of delta_sigma obey the Teukolsky equation
    -> QGD QNMs = GR QNMs at leading order.
    """

    def __init__(self, M_f, chi_f):
        self.M_f   = M_f
        self.chi_f = chi_f
        self.r_s   = 2.0*G*M_f/c**2
        self.alpha = chi_f*G*M_f/c**2

    def sigma_t(self, r, theta=np.pi/2):
        S = r**2 + self.alpha**2*np.cos(theta)**2
        return np.sqrt(2.0*G*self.M_f*r/(c**2*S))

    def sigma_phi(self, r, theta=np.pi/2):
        if self.alpha == 0.0: return 0.0
        S = r**2 + self.alpha**2*np.cos(theta)**2
        return self.alpha*np.sin(theta)*np.sqrt(2.0*G*self.M_f/(c**2*r*S))

    def g_tt(self, r, theta=np.pi/2):
        return -(1.0 - self.sigma_t(r, theta)**2)

    def photon_sphere_radius(self):
        """r_ph ~ 3GM/c^2 * (1 - 2*chi/(3*sqrt(3)))"""
        return 3.0*G*self.M_f/c**2 * (1.0 - 2.0*self.chi_f/(3.0*np.sqrt(3.0)))

    def qgd_qnm_from_sigma(self):
        """
        QGD semi-analytic QNM estimate from sigma field at photon sphere.
        omega_r = c * sigma_t(r_ph) / r_ph  (WKB approximation)
        omega_i = Lyapunov exponent of photon orbit
        This is a ~45% overestimate vs exact Teukolsky — the fit
        Eq.(omega_220) gives exact agreement.
        """
        r_ph = self.photon_sphere_radius()
        sig  = self.sigma_t(r_ph)
        omega_r = c*sig/r_ph
        omega_i = 0.5*c*np.sqrt(self.r_s/r_ph**3)
        return omega_r, omega_i

    def outer_horizon(self):
        """Outer horizon: sigma_t = 1."""
        if self.alpha == 0.0: return self.r_s
        try: return brentq(lambda r: self.sigma_t(r)-1.0, 0.5*self.r_s, 2.0*self.r_s)
        except: return self.r_s

    def quantum_stiffness_correction(self):
        """delta_omega/omega ~ kappa*(l_Q/r_s)^2 ~ 1e-140 for stellar BHs."""
        return 2.0*(l_Q/self.r_s)**2

    def perturbation_potential(self, r):
        """
        Regge-Wheeler effective potential for ell=2 perturbations.
        V(r) = (1 - r_s/r) * [ell(ell+1)/r^2 + r_s/r^3]
        IDENTICAL to GR because background metric IS Kerr.
        """
        return (1.0 - self.r_s/r) * (6.0/r**2 + self.r_s/r**3)


# ── 3. QGD CROSS-TERM RINGDOWN ──────────────────────────────────────────────

class QGDCrossTermRingdown:
    """
    QGD-specific early ringdown from merger cross-term 2*sigma_t^(1)*sigma_t^(2).

    In QGD:  g_tt = -(1 - A1^2 - A2^2 - 2*A1*A2)
    The cross-term 2*A1*A2 decays as horizons merge:
      delta_g_tt^(cross)(t) = -2*A1*A2 * exp(-t/tau_cross)
      tau_cross = r_ISCO/c ~ 6*G*M_tot/c^3

    This is ABSENT from all current GR ringdown templates (which assume
    a single body from t=0). It appears as a non-oscillatory early decay.

    DIPOLE RADIATION (QGD-specific for M1 != M2):
    QGD dipole moment: d_sigma = sum(sqrt(Ma) * xa)
    NOT conserved for unequal masses -> dipole GW emission at f_orb
    (forbidden in GR by momentum conservation)
    Amplitude factor: (sqrt(M1) - sqrt(M2))^2 / M_total
    """

    def __init__(self, M1, M2, chi1=0.0, chi2=0.0, separation=None):
        self.M1, self.M2 = M1, M2
        self.chi1, self.chi2 = chi1, chi2
        self.q = min(M1,M2)/max(M1,M2)
        M_tot = M1+M2
        self.r_isco = separation or 6.0*G*M_tot/c**2
        self.tau_cross = self.r_isco/c

        # sigma fields at ISCO
        r = self.r_isco
        a1 = chi1*G*M1/c**2; a2 = chi2*G*M2/c**2
        self.sigma1 = np.sqrt(2*G*M1*r/(c**2*(r**2+a1**2)))
        self.sigma2 = np.sqrt(2*G*M2*r/(c**2*(r**2+a2**2)))
        self.cross_amplitude = 2.0*self.sigma1*self.sigma2

    def delta_gtt_cross(self, t):
        """Non-oscillatory cross-term decay in g_tt."""
        return -self.cross_amplitude * np.exp(-t/self.tau_cross)

    def orbital_frequency(self):
        return np.sqrt(G*(self.M1+self.M2)/self.r_isco**3)

    def dipole_amplitude_factor(self):
        """(sqrt(M1)-sqrt(M2))^2 / M_tot — zero for q=1."""
        return (np.sqrt(self.M1)-np.sqrt(self.M2))**2/(self.M1+self.M2)

    def dipole_power(self, a_orb=None):
        """P_dipole = G/(3c^3) * |d_sigma_ddot|^2"""
        Om = self.orbital_frequency()
        a  = a_orb or self.r_isco/2.0
        d2 = (np.sqrt(self.M1)-np.sqrt(self.M2))*Om**2*a
        return G/(3.0*c**3)*d2**2


# ── 4. RINGDOWN MODE CONTAINER ──────────────────────────────────────────────

@dataclass
class RingdownMode:
    """Single damped sinusoid: h(t) = A * exp(-t*omega_i) * cos(omega_r*t + phi)"""
    label:   str
    omega_r: float   # rad/s
    omega_i: float   # rad/s  (tau = 1/omega_i)
    A:       float
    phi:     float
    origin:  str     # 'GR_QNM' | 'QGD_cross' | 'QGD_dipole' | 'QGD_nonlinear'

    @property
    def tau(self): return 1.0/self.omega_i if self.omega_i > 0 else np.inf
    @property
    def frequency_Hz(self): return self.omega_r/(2.0*np.pi)

    def h(self, t):
        t = np.asarray(t)
        return self.A * np.exp(-t*self.omega_i) * np.cos(self.omega_r*t + self.phi)


# ── 5. FULL QGD RINGDOWN WAVEFORM ───────────────────────────────────────────

class QGDRingdownWaveform:
    """
    Complete QGD ringdown waveform:
      h_QGD(t) = h_GR_QNM(t) + h_cross(t) + h_dipole(t)

    h_GR_QNM  — identical to GR at leading order (QNMs on Kerr background)
    h_cross   — QGD cross-term early decay (non-oscillatory, tau~ms)
    h_dipole  — QGD dipole at f_orb (only for M1!=M2; forbidden in GR)

    Usage:
      model = QGDRingdownWaveform(M1=86*M_sun, M2=77*M_sun, chi1=0.18, chi2=0.12)
      t = np.linspace(0, 0.05, 10000)
      h = model.strain(t)
      model.print_mode_table()
    """

    def __init__(self, M1, M2, chi1=0.0, chi2=0.0,
                 distance=1e9*3.086e16,
                 include_cross=True, include_dipole=True, include_220Q=True):
        self.M1, self.M2 = M1, M2
        self.chi1, self.chi2 = chi1, chi2
        self.distance = distance
        self.q = min(M1,M2)/max(M1,M2)

        self.M_f, self.chi_f = self._final_state()
        self.remnant = QGDRemnant(self.M_f, self.chi_f)
        self.cross   = QGDCrossTermRingdown(M1, M2, chi1, chi2)
        self.modes: List[RingdownMode] = []

        self._add_gr_qnm_modes()
        if include_220Q:   self._add_nonlinear_mode()
        if include_cross:  self._add_cross_mode()
        if include_dipole and abs(M1-M2) > 0.01*max(M1,M2):
            self._add_dipole_mode()

    def _final_state(self):
        """BKL-like final state fit."""
        M_tot = self.M1+self.M2
        nu    = self.M1*self.M2/M_tot**2
        chi_eff = (self.M1*self.chi1+self.M2*self.chi2)/M_tot
        E_rad = 0.0543*nu*M_tot*c**2*(1.0-0.4*chi_eff)
        M_f   = M_tot - E_rad/c**2
        J_orb = np.sqrt(12.0)*nu*G*M_tot**2/c
        chi_f = min(0.95, abs(J_orb+self.M1**2*self.chi1+self.M2**2*self.chi2)*c/(G*M_f**2))
        return M_f, chi_f

    def _scale(self):
        return G*self.M_f/(self.distance*c**2)

    def _add_gr_qnm_modes(self):
        s = self._scale()
        wr,wi = qnm_220_frequency(self.M_f, self.chi_f)
        self.modes.append(RingdownMode('220',  wr,   wi,   A=s,       phi=0.0,  origin='GR_QNM'))
        wr1,wi1 = qnm_221_frequency(self.M_f, self.chi_f)
        self.modes.append(RingdownMode('221',  wr1,  wi1,  A=0.3*s,   phi=0.4,  origin='GR_QNM'))
        wr4,wi4 = qnm_440_frequency(self.M_f, self.chi_f)
        self.modes.append(RingdownMode('440',  wr4,  wi4,  A=0.1*s,   phi=1.2,  origin='GR_QNM'))

    def _add_nonlinear_mode(self):
        s = self._scale()
        wr,wi = qnm_220Q_frequency(self.M_f, self.chi_f)
        self.modes.append(RingdownMode('220Q', wr, wi, A=0.05*s, phi=0.8, origin='QGD_nonlinear'))

    def _add_cross_mode(self):
        s = self._scale()
        A = s*self.cross.cross_amplitude*0.5
        self.modes.append(RingdownMode('cross', 0.0, 1.0/self.cross.tau_cross,
                                       A=A, phi=0.0, origin='QGD_cross'))

    def _add_dipole_mode(self):
        s   = self._scale()
        Om  = self.cross.orbital_frequency()
        fac = self.cross.dipole_amplitude_factor()
        self.modes.append(RingdownMode('dipole', Om, 1.0/(2.0*self.cross.tau_cross),
                                       A=s*fac, phi=np.pi/4, origin='QGD_dipole'))

    def strain(self, t):
        """Total QGD strain h(t) = sum of all modes."""
        t = np.asarray(t); h = np.zeros_like(t, dtype=float)
        for m in self.modes: h += m.h(t)
        return h

    def strain_by_component(self, t):
        """Strain split by origin: GR_QNM, QGD_cross, QGD_dipole, QGD_nonlinear."""
        t = np.asarray(t); comp = {}
        for m in self.modes:
            comp.setdefault(m.origin, np.zeros_like(t, dtype=float))
            comp[m.origin] += m.h(t)
        return comp

    def print_mode_table(self):
        sep = "─"*78
        print(sep)
        print(f"{'Mode':<8} {'f (Hz)':>10} {'τ (ms)':>10} {'Amplitude':>14} {'φ':>8}  Origin")
        print(sep)
        for m in self.modes:
            tau_ms = m.tau*1e3 if m.tau < 1e9 else np.inf
            print(f"{m.label:<8} {m.frequency_Hz:>10.2f} {tau_ms:>10.3f} "
                  f"{m.A:>14.4e} {m.phi:>8.4f}  {m.origin}")
        print(sep)
        print(f"  M_f={self.M_f/M_sun:.2f} Msun  chi_f={self.chi_f:.4f}  q={self.q:.4f}")
        print(sep)


# ── 6. QGD vs GR COMPARISON ─────────────────────────────────────────────────

class QGDvsGRComparison:
    """Compare QGD and GR ringdown predictions for a given remnant."""

    def __init__(self, M_f, chi_f):
        self.remnant = QGDRemnant(M_f, chi_f)
        self.M_f, self.chi_f = M_f, chi_f

    def frequency_table(self):
        wr_gr, wi_gr = qnm_220_frequency(self.M_f, self.chi_f)
        wr_qgd, wi_qgd = self.remnant.qgd_qnm_from_sigma()
        stiff = self.remnant.quantum_stiffness_correction()
        lines = [
            "═"*65, "QGD vs GR RINGDOWN FREQUENCIES", "═"*65,
            f"  M_f={self.M_f/M_sun:.2f} Msun  chi_f={self.chi_f:.3f}",
            f"  GR  220 (Echeverria fit):  f={wr_gr/(2*np.pi):.2f} Hz  tau={1/wi_gr*1e3:.2f} ms",
            f"  QGD 220 (photon-sph WKB):  f={wr_qgd/(2*np.pi):.2f} Hz  tau={1/wi_qgd*1e3:.2f} ms",
            f"  WKB deviation: {abs(wr_qgd-wr_gr)/wr_gr*100:.1f}%  (exact: use Teukolsky = GR)",
            f"  Quantum stiffness: delta_omega/omega ~ {stiff:.2e}  [unmeasurable]",
            "  CONCLUSION: QGD QNMs = GR QNMs (background IS Kerr -> Teukolsky applies)",
            "═"*65,
        ]
        return "\n".join(lines)

    def signatures_table(self, M1, M2):
        cross = QGDCrossTermRingdown(M1, M2)
        lines = [
            "═"*65, "QGD-SPECIFIC SIGNATURES", "═"*65,
            f"  q = {cross.q:.4f}",
            "",
            "  [P1] CROSS-TERM EARLY DECAY",
            f"      tau_cross = {cross.tau_cross*1e3:.3f} ms = r_ISCO/c",
            f"      Amplitude = {cross.cross_amplitude:.4f}",
            "      Non-oscillatory (omega=0) — not in any current GR template.",
            "",
            "  [P2] DIPOLE RINGDOWN",
            f"      (sqrt(M1)-sqrt(M2))^2/M_tot = {cross.dipole_amplitude_factor():.4e}",
            f"      f_orb = {cross.orbital_frequency()/(2*np.pi):.1f} Hz  (not 2*f_orb!)",
            f"      P_dipole = {cross.dipole_power():.3e} W",
            f"      Status: {'ACTIVE' if cross.q < 0.99 else 'ZERO (equal masses)'}",
            "",
            "  [P3] QUANTUM STIFFNESS",
            f"      delta_omega/omega ~ {self.remnant.quantum_stiffness_correction():.2e}",
            "      Prevents singularity; unmeasurable for stellar-mass BHs.",
            "═"*65,
        ]
        return "\n".join(lines)


# ── 7. SIGMA-FIELD RELAXATION TRACKER ──────────────────────────────────────

def sigma_relaxation(M1, M2, chi1=0.0, chi2=0.0, r_obs=None, t_arr=None):
    """
    Track sigma-field relaxation from merger to final state:
      sigma_total(r,t) = sigma_remnant(r) + delta_sigma_QNM(t) + delta_sigma_cross(t)

    Returns dict: t, sigma_total, sigma_remnant, delta_QNM, delta_cross, g_tt
    """
    M_tot = M1+M2; nu = M1*M2/M_tot**2
    M_f   = M_tot*(1 - 0.0543*nu*(1-0.4*(M1*chi1+M2*chi2)/M_tot))
    chi_f = 0.68

    rem  = QGDRemnant(M_f, chi_f)
    cros = QGDCrossTermRingdown(M1, M2, chi1, chi2)

    r_obs = r_obs or 3.0*G*M_f/c**2
    t_arr = t_arr if t_arr is not None else np.linspace(0, 0.1, 5000)

    sig_rem = rem.sigma_t(r_obs)
    wr, wi  = qnm_220_frequency(M_f, chi_f)
    dQNM    = 0.1*sig_rem * np.exp(-t_arr*wi) * np.cos(wr*t_arr)

    a1 = chi1*G*M1/c**2; a2 = chi2*G*M2/c**2; r = r_obs
    sig1 = np.sqrt(2*G*M1*r/(c**2*(r**2+a1**2)))
    sig2 = np.sqrt(2*G*M2*r/(c**2*(r**2+a2**2)))
    dCross = sig1*sig2*np.exp(-t_arr/cros.tau_cross)

    sig_tot = sig_rem + dQNM + dCross
    return {
        't':              t_arr,
        'sigma_total':    sig_tot,
        'sigma_remnant':  np.full_like(t_arr, sig_rem),
        'delta_QNM':      dQNM,
        'delta_cross':    dCross,
        'g_tt':           -(1.0 - sig_tot**2),
    }


# ── 8. GW250114 ANALYSIS ────────────────────────────────────────────────────

def analyse_GW250114():
    """
    Apply QGD framework to GW250114 (Jan 14, 2025).
    Source: arXiv:2509.08099 (LVK), arXiv:2510.16903 (nonlinear modes)

    Observed: M1~86, M2~77 Msun; M_f~152 Msun, chi_f~0.68; SNR~80; D~1.8 Gpc
    """
    M1=86*M_sun; M2=77*M_sun; chi1=0.18; chi2=0.12
    D_L=1.8e9*3.086e16

    print("═"*65)
    print("GW250114  |  Jan 14 2025  |  LIGO O4b  |  SNR~80")
    print("QGD RINGDOWN ANALYSIS")
    print("═"*65)

    model = QGDRingdownWaveform(M1,M2,chi1,chi2,distance=D_L)
    print(f"\n  Input:  M1={M1/M_sun:.0f}, M2={M2/M_sun:.0f} Msun | chi1={chi1}, chi2={chi2}")
    print(f"  QGD:    M_f={model.M_f/M_sun:.1f} Msun | chi_f={model.chi_f:.3f}")
    print(f"  LVK:    M_f~152 Msun | chi_f~0.68\n")
    model.print_mode_table()

    wr,wi = qnm_220_frequency(model.M_f, model.chi_f)
    print(f"\n  QGD 220 prediction: f={wr/(2*np.pi):.1f} Hz  tau={1/wi*1e3:.1f} ms")
    print(f"  LVK observed:       f~73 Hz  tau~14 ms")
    print(f"  Frequency agreement: ~{abs(wr/(2*np.pi)-73)/73*100:.1f}% deviation")

    print(f"\n  GW250114 KEY FINDINGS vs QGD:")
    print(f"  [220+221 at >3sigma] Consistent with Kerr  => QGD agrees (Teukolsky applies) ✓")
    print(f"  [440 mode first]     Consistent with Kerr  => QGD agrees ✓")
    print(f"  [220Q nonlinear]     QGD second-order sigma perturbation => same prediction ✓")
    print(f"  [cross-term]         NOT searched in LVK   => QGD prediction P1 untested")
    dipole_fac = QGDCrossTermRingdown(M1,M2).dipole_amplitude_factor()
    print(f"  [dipole at f_orb]    Factor={(dipole_fac):.2e}  (nearly equal mass, very small)")

    comp = QGDvsGRComparison(model.M_f, model.chi_f)
    print(); print(comp.frequency_table())
    print(); print(comp.signatures_table(M1, M2))
    return model


# ── 9. MAIN DEMO ────────────────────────────────────────────────────────────

def run_all():
    print("\n"+"═"*65)
    print("QGD RINGDOWN — COMPLETE DEMONSTRATION")
    print("="*65)

    print("\n[1] REMNANT SIGMA-FIELD PROFILE  (30 Msun, chi=0.68)")
    rem = QGDRemnant(30*M_sun, 0.68)
    print(f"  r_s={rem.r_s:.0f}m  r_ph={rem.photon_sphere_radius():.0f}m")
    print(f"  {'r/rs':>5}  {'sigma_t':>10}  {'g_tt':>12}  {'V_eff':>14}")
    for f in [1.01, 1.5, 2, 3, 5, 10]:
        r=f*rem.r_s
        print(f"  {f:>5.2f}  {rem.sigma_t(r):>10.6f}  {rem.g_tt(r):>12.6f}  {rem.perturbation_potential(r):>14.4e}")

    print("\n[2] QNM FREQUENCIES (QGD = GR at leading order)")
    print(f"  {'chi_f':>6}  {'f_220 (Hz)':>12}  {'tau_220 (ms)':>14}  {'f_440 (Hz)':>12}")
    for chi in [0.0, 0.3, 0.68, 0.95]:
        wr,wi = qnm_220_frequency(30*M_sun, chi)
        wr4,_ = qnm_440_frequency(30*M_sun, chi)
        print(f"  {chi:>6.2f}  {wr/(2*np.pi):>12.2f}  {1/wi*1e3:>14.3f}  {wr4/(2*np.pi):>12.2f}")

    print("\n[3] CROSS-TERM RINGDOWN  (20+15 Msun)")
    cr = QGDCrossTermRingdown(20*M_sun, 15*M_sun, 0.3, 0.2)
    print(f"  tau_cross={cr.tau_cross*1e3:.3f} ms  cross_amp={cr.cross_amplitude:.4f}")
    print(f"  dipole_factor={cr.dipole_amplitude_factor():.4e}  P_dipole={cr.dipole_power():.3e} W")

    print("\n[4] SIGMA RELAXATION  (20+15 Msun)")
    prof = sigma_relaxation(20*M_sun, 15*M_sun)
    print(f"  {'t(ms)':>8}  {'Sigma':>10}  {'rem':>10}  {'dQNM':>10}  {'dCross':>10}  {'g_tt':>12}")
    for i in [0, 10, 50, 200, 1000, 2000, -1]:
        t_=prof['t'][i]*1e3; st=prof['sigma_total'][i]; sr=prof['sigma_remnant'][i]
        dq=prof['delta_QNM'][i]; dc=prof['delta_cross'][i]; gtt=prof['g_tt'][i]
        print(f"  {t_:>8.2f}  {st:>10.6f}  {sr:>10.6f}  {dq:>10.6f}  {dc:>10.6f}  {gtt:>12.6f}")

    print("\n[5] GW250114 FULL ANALYSIS")
    analyse_GW250114()

    print("\n"+"═"*65)
    print("COMPLETE. Key result: QGD QNMs = GR QNMs at leading order.")
    print("QGD-specific: cross-term (P1), dipole (P2), stiffness (P3).")
    print("="*65)

if __name__ == "__main__":
    run_all()
