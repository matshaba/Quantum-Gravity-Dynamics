"""
qgd_three_problems_fixed.py
============================
QGD: Three Foundational Problems  —  v1.2  (bug-corrected)

  PROBLEM 1 — The QGD Effective One-Body (EOB)
  PROBLEM 2 — The Merger Transition with exact merger condition
  PROBLEM 3 — The Spinning Master Formula

BUG-FIX CHANGELOG (v1.1 → v1.2)
══════════════════════════════════════════════════════════════════════

BUG-FIX 1 — GW250114 masses corrected
  WRONG:   86 + 77 Msun
  CORRECT: 33.6 + 32.2 Msun  (total 65.8, final BH 62.7 Msun, chi_f=0.68)
  Source:  arXiv:2509.08054; PRL 135, 101101 (2025)
  Note:    GW250114 is the loudest GW event ever (SNR=80). It confirmed
           the Kerr QNM spectrum (fundamental + first overtone) to tens
           of percent — directly testing QGD's ringdown = Kerr prediction.

BUG-FIX 2 — Spin-squared PN order labels off by 1
  WRONG:   pn = f'{n+3} PN'   (n=0 labeled '3 PN')
  CORRECT: pn = f'{n+2} PN'   (n=0 labeled '2 PN')
  Reason:  The SS term is chi^2 * u^{n+3}. For circular orbits u ~ (v/c)^2,
           so u^{n+3} counts as (n+3) relative powers of v^2, but the
           LEADING Newtonian energy is already ~u (1 power), so each
           additional u adds one PN order. Leading SS (n=0) is u^3 = 2PN.
           This matches known_ss dict: {0:'2PN', 1:'3PN', 2:'4PN'}.

BUG-FIX 3 — GR EOB A-function comparison formula wrong
  WRONG:   A_gr_pn = 1 - 2u + (94/3 - 41*pi^2/32)*eta * u^2
           [u^2 coefficient is ZERO in GR EOB; eta-correction enters at u^3]
  CORRECT: A_gr_pn = 1 - 2u + 2*eta*u^3 + (94/3 - 41*pi^2/32)*eta*u^4
           [Damour, Jaranowski, Schäfer 2000 — PRD 62, 084011 — 3PN exact]
  Impact:  At u=0.1, eta=0.25 the old code introduced a ~6% error in A(u)
           vs the correct ~0.1% eta-correction; comparison table was misleading.

BUG-FIX 4 — Remnant mass formula lacked EMRI limit
  WRONG:   Mf = M*(1 - 0.0539*4*eta)  → 5.39% radiated for equal mass
           [violates test-body limit; overestimates GW250114 by ~0.7%]
  CORRECT: Mf = M*(1 - eta*(E_ISCO + 0.5712*eta))
           where E_ISCO = 1 - sqrt(8/9) ~ 0.05719 (Schwarzschild ISCO)
  Properties:
    • EMRI limit: E_rad -> eta*E_ISCO as eta->0  [exact test-body limit]
    • Equal mass: E_rad/M ~ 5.0%  [NR SXS catalog: 4.8-5.0%]
    • GW250114:   eta~0.25, measured 4.71% — formula gives 5.0% (~0.3% off)
  Source:  Barausse et al. 2012 (1206.3803); Hofmann et al. 2016 (1605.01938)
"""

import numpy as np
from fractions import Fraction
from math import comb, factorial

G     = 6.674e-11
c     = 3.000e8
M_sun = 1.989e30
EULER = 0.5772156649015329

# Schwarzschild ISCO energy — exact test-body limit for radiated energy
E_ISCO_SCHW = 1.0 - np.sqrt(8.0/9.0)   # ≈ 0.05719 c^2 per unit mu


# ═══════════════════════════════════════════════════════════════════════
# PROBLEM 2: MERGER CONDITION
# ═══════════════════════════════════════════════════════════════════════

def qgd_merger_condition(M1, M2):
    """
    EXACT QGD merger separation from sigma-field saddle point.

    The sigma-field Sigma = sqrt(2GM1/c^2/r1) + sqrt(2GM2/c^2/r2) reaches
    its maximum on the line joining the bodies at:

        r1* = d/(1+rho),   r2* = d*rho/(1+rho),   rho = (M2/M1)^{1/3}

    Setting Sigma(saddle) = 1 gives the EXACT merger condition:

        d_merge = 2*G*M1*(1 + (M2/M1)^{1/3})^3 / c^2

    Equal masses (M1=M2=m):
        d = 16Gm/c^2 = 8 * r_s(individual) = 4 * r_s(total)  [BOTH correct]
    EMRI limit (M2->0):
        d -> 2GM1/c^2 = r_s(primary)  [photon orbit of primary]
    """
    M   = M1 + M2
    eta = M1*M2/M**2
    rho = (M2/M1)**(1.0/3.0)
    d   = 2*G*M1*(1+rho)**3/c**2

    r1_saddle = d/(1+rho)
    r2_saddle = d*rho/(1+rho)

    Sigma_check = (np.sqrt(2*G*M1/c**2/r1_saddle) +
                   np.sqrt(2*G*M2/c**2/r2_saddle))
    Sigma_mid   = (np.sqrt(4*G*M1/c**2/d) + np.sqrt(4*G*M2/c**2/d))

    rs_tot = 2*G*M/c**2
    rs_1   = 2*G*M1/c**2

    return {
        'd_merge':      d,
        'r1_saddle':    r1_saddle,
        'r2_saddle':    r2_saddle,
        'Sigma_saddle': Sigma_check,
        'Sigma_mid':    Sigma_mid,
        'saddle_exact': abs(Sigma_check-1) < 1e-10,
        'd_over_rstot': d/rs_tot,
        'd_over_rs1':   d/rs_1,
        'eta':          eta,
        'q':            min(M1,M2)/max(M1,M2),
        'rho':          rho,
    }


def print_merger_resolution():
    print("="*75)
    print("RESOLVING d=4rs vs d=8rs")
    print("Three candidate merger conditions:")
    print("  (A) Midpoint:    d_A = 4G(sqrt(M1)+sqrt(M2))^2/c^2   [approx, q=1 only]")
    print("  (B) Saddle pt:   d_B = 2GM1*(1+(M2/M1)^{1/3})^3/c^2  [EXACT, all q]")
    print("  (C) EOB A=0:     d_C = 2GM/(eta*c^2)                  [EOB singularity]")
    print("-"*75)

    M1 = 10*M_sun

    def dA(M2): return 4*G*(np.sqrt(M1)+np.sqrt(M2))**2/c**2
    def dB(M2): return qgd_merger_condition(M1,M2)['d_merge']
    def dC(M2):
        M = M1+M2; eta = M1*M2/M**2; return 2*G*M/(eta*c**2)

    norm = 16*G*M1/c**2

    print(f"\n  All normalised by 16Gm/c^2 = 8*rs_indiv (equal-mass answer):")
    print(f"  {'q':>6}  {'eta':>6}  {'dA/8rs':>8}  {'dB/8rs':>8}  "
          f"{'dC/8rs':>8}  {'dB/rs_tot':>10}  {'Sig_saddle':>12}")
    print("-"*75)
    for q in [1.0, 0.5, 0.25, 0.1, 0.05, 0.01]:
        M2 = q*M1
        r  = qgd_merger_condition(M1, M2)
        print(f"  {q:>6.3f}  {r['eta']:>6.4f}  {dA(M2)/norm:>8.4f}  "
              f"{dB(M2)/norm:>8.4f}  {dC(M2)/norm:>8.4f}  "
              f"{r['d_over_rstot']:>10.4f}  {r['Sigma_saddle']:>12.8f}")
    print("-"*75)
    print("  KEY FINDINGS:")
    print("  1. d_B gives Sigma(saddle)=1.000000 for ALL q  [physically exact]")
    print("  2. d_A (midpoint) is only exact at q=1; underestimates for q<1")
    print("  3. d_C (EOB A=0) diverges for q->0  [different physics from merger]")
    print("  4. All three agree at q=1: d = 8rs_indiv = 4rs_total (both correct!)")
    print("="*75)


# ═══════════════════════════════════════════════════════════════════════
# REMNANT FITS
# ═══════════════════════════════════════════════════════════════════════

def remnant_mass_frac(eta):
    """
    Fraction of total mass radiated as GWs (nonspinning, quasi-circular).

    E_rad/M = eta * (E_ISCO_SCHW + 0.5712 * eta)

    where E_ISCO_SCHW = 1 - sqrt(8/9) ~ 0.05719 is the test-body EMRI limit.
    b=0.5712 calibrated so equal mass gives ~5.0% (NR: 4.8-5.0%).

    [BUG-FIX 4: old 0.0539*4*eta violated EMRI limit and overestimated by ~0.7%]
    """
    return eta * (E_ISCO_SCHW + 0.5712*eta)


def remnant_spin(eta):
    """
    Final spin for nonspinning merger: chi_f ~ 0.6864*(4*eta)^0.45
    Equal mass: chi_f = 0.6864 -- matches GW250114 measured 0.68 +/- 0.01.
    Ref: Barausse & Rezzolla 2009 (ApJ 704, L40).
    """
    return min(0.9999, 0.6864*(4.0*eta)**0.45)


# ═══════════════════════════════════════════════════════════════════════
# PROBLEM 1: QGD EFFECTIVE ONE-BODY
# ═══════════════════════════════════════════════════════════════════════

class QGDEOB:
    """
    QGD EOB: A^QGD(u) = 1 - 2u/eta  (exact, no free parameters)

    Derived from Sigma_tot^2 = 2u/eta in the COM frame.
    Compare to GR EOB at 3PN:
      A^GR(u) = 1 - 2u + 2*eta*u^3 + (94/3 - 41*pi^2/32)*eta*u^4
      [Damour-Jaranowski-Schäfer 2000, PRD 62, 084011]

    Note: GR has a^{GR}_2 = 0 (NO u^2 term); first eta-correction at u^3.
    QGD has ALL nonlinear coefficients = 0. This is the key testable difference.

    Orbital structure (analogous to Schwarzschild with u -> u/eta):
      ISCO:         u_ISCO = eta/6
      Light ring:   u_LR   = eta/3
      EOB horizon:  u_h    = eta/2   (A^QGD = 0)
    """

    def __init__(self, M1, M2, chi1=0.0, chi2=0.0):
        self.M1, self.M2 = M1, M2
        self.M   = M1 + M2
        self.mu  = M1*M2/self.M
        self.eta = M1*M2/self.M**2
        self.q   = min(M1,M2)/max(M1,M2)
        self.chi_eff = (M1*chi1 + M2*chi2)/self.M

        self.u_ISCO = self.eta/6.0
        self.u_LR   = self.eta/3.0
        self.u_EOB  = self.eta/2.0

        mc = qgd_merger_condition(M1, M2)
        self.u_merge_phys = G*self.M/(c**2*mc['d_merge'])

    def A(self, u):
        """A^QGD(u) = 1 - 2u/eta."""
        return 1.0 - 2.0*u/self.eta

    def A_gr_schw(self, u):
        """GR EOB test-body limit (eta->0): A = 1 - 2u."""
        return 1.0 - 2.0*u

    def A_gr_pn(self, u):
        """
        GR EOB A-function at 3PN order with eta-corrections.
        [BUG-FIX 3: was 1-2u + a4*eta*u^2, missing u^3 term, wrong power]

        A^GR(u) = 1 - 2u + 2*eta*u^3 + (94/3 - 41*pi^2/32)*eta*u^4
        """
        a4 = 94.0/3.0 - 41.0*np.pi**2/32.0   # ~ 18.688
        return 1.0 - 2.0*u + 2.0*self.eta*u**3 + a4*self.eta*u**4

    def omega_circ(self, u):
        return c**3/(G*self.M) * u**1.5

    def f_gw(self, u):
        return 2.0*self.omega_circ(u)/(2.0*np.pi)

    def print_structure(self):
        print(f"QGD EOB: M1={self.M1/M_sun:.1f}, M2={self.M2/M_sun:.1f} Msun  "
              f"eta={self.eta:.4f}  q={self.q:.3f}")
        rows = [
            ("ISCO",            self.u_ISCO,       "start of plunge"),
            ("Light ring",      self.u_LR,         "unstable photon orbit"),
            ("EOB singularity", self.u_EOB,         "A^QGD=0"),
            ("Phys. merger",    self.u_merge_phys, "Sigma_saddle=1"),
        ]
        rs = 2*G*self.M/c**2
        print(f"  {'Feature':>22}  {'u=GM/c^2 r':>10}  {'r/rs':>8}  "
              f"{'A^QGD':>8}  {'A^GR_Schw':>10}  {'f_GW(Hz)':>12}")
        print("-"*80)
        for label, u, note in rows:
            r   = G*self.M/(c**2*u)
            A_q = self.A(u)
            A_g = self.A_gr_schw(u)
            f   = self.f_gw(u)
            print(f"  {label:>22}  {u:>10.4f}  {r/rs:>8.2f}  "
                  f"{A_q:>8.4f}  {A_g:>10.4f}  {f:>12.1f}   [{note}]")
        print("-"*80)
        err = abs(self.u_merge_phys - self.u_EOB)/self.u_EOB*100
        if err < 0.01:
            print(f"  EOB singularity = physical merger (q~1) CONSISTENT")
        else:
            print(f"  EOB singularity != physical merger: {err:.1f}% diff at q={self.q:.3f}")
            print(f"  (Merger occurs during plunge, consistent with GR EOB)")

    def compare_table(self):
        """A^QGD vs GR EOB (Schwarzschild + 3PN) comparison."""
        rs = 2*G*self.M/c**2
        print(f"\n  {'u':>7}  {'r/rs':>6}  {'A^QGD':>9}  "
              f"{'A^GR_Schw':>11}  {'A^GR_3PN':>11}  {'QGD/Schw':>9}  phase")
        print("-"*80)
        u_list = [0.005, 0.01, 0.02, 0.05, self.u_ISCO, self.u_LR,
                  self.u_EOB,
                  min(self.u_merge_phys, 0.9*self.u_EOB + 0.1*self.u_merge_phys)]
        for u in sorted(set(u_list)):
            if u <= 0 or u >= 0.95: continue
            r    = G*self.M/(c**2*u)
            A_q  = self.A(u)
            A_g  = self.A_gr_schw(u)
            A_3pn= self.A_gr_pn(u)
            ratio = A_q/A_g if abs(A_g)>1e-6 else float('inf')
            if   u <= self.u_ISCO: ph = "inspiral"
            elif u <= self.u_LR:   ph = "plunge"
            elif u <= self.u_EOB:  ph = "EOB zone"
            else:                  ph = "post-EOB"
            print(f"  {u:>7.4f}  {r/rs:>6.2f}  {A_q:>9.5f}  "
                  f"{A_g:>11.5f}  {A_3pn:>11.5f}  {ratio:>9.5f}  {ph}")
        print("-"*80)


# ═══════════════════════════════════════════════════════════════════════
# PROBLEM 2 (cont.): IMR TRANSITION
# ═══════════════════════════════════════════════════════════════════════

class QGDMergerIMR:
    """
    QGD Inspiral-Merger-Ringdown in a single framework.

    Ringdown = Kerr QNM spectrum.
    GW250114 (Sep 2025, SNR=80): confirmed fundamental+overtone of Kerr.
    QGD prediction STATUS: PASSES.
    """

    def __init__(self, M1, M2, chi1=0.0, chi2=0.0):
        self.M1, self.M2 = M1, M2
        self.M   = M1 + M2
        self.eta = M1*M2/self.M**2
        self.eob = QGDEOB(M1, M2, chi1, chi2)

        mc = qgd_merger_condition(M1, M2)
        self.d_merge = mc['d_merge']
        self.u_merge = G*self.M/(c**2*self.d_merge)

        e_rad    = remnant_mass_frac(self.eta)
        self.Mf  = self.M*(1.0 - e_rad)
        self.chif = remnant_spin(self.eta)

    def sigma_saddle(self, d):
        rho = (self.M2/self.M1)**(1.0/3.0)
        r1 = d/(1+rho); r2 = d*rho/(1+rho)
        return (np.sqrt(2*G*self.M1/c**2/r1) +
                np.sqrt(2*G*self.M2/c**2/r2))

    def phase(self, u):
        if   u < self.eob.u_ISCO:  return "inspiral"
        elif u < self.eob.u_LR:    return "plunge"
        elif u < self.u_merge:     return "merger approach"
        else:                      return "ringdown"

    def print_imr(self):
        rs = 2*G*self.M/c**2
        e_pct = remnant_mass_frac(self.eta)*100
        print(f"QGD IMR: M1={self.M1/M_sun:.1f}+M2={self.M2/M_sun:.1f} Msun  "
              f"eta={self.eta:.4f}  Mf={self.Mf/M_sun:.1f} Msun  "
              f"chi_f={self.chif:.3f}  (Erad={e_pct:.2f}%)")
        print()
        u_arr = np.array([0.005, 0.01, self.eob.u_ISCO,
                          self.eob.u_LR*0.8, self.eob.u_LR,
                          self.eob.u_EOB*0.95, self.u_merge])
        print(f"  {'u':>8}  {'r/rs':>6}  {'A^QGD':>8}  "
              f"{'Sig_saddle':>12}  {'f_GW(Hz)':>10}  {'phase':>18}")
        print("-"*75)
        for u in u_arr:
            r   = G*self.M/(c**2*u)
            A_q = max(0.0, self.eob.A(u))
            Sig = self.sigma_saddle(r)
            f   = self.eob.f_gw(u)
            ph  = self.phase(u)
            print(f"  {u:>8.4f}  {r/rs:>6.2f}  {A_q:>8.4f}  "
                  f"{Sig:>12.6f}  {f:>10.1f}  {ph:>18}")
        print("-"*75)
        print(f"  Merger: u={self.u_merge:.4f}  d={self.d_merge:.3e} m  Sigma=1.000")
        d_eob = G*self.M/(c**2*self.eob.u_EOB)
        err   = abs(self.d_merge - d_eob)/d_eob*100
        if err < 1.0:
            print(f"  EOB merger ~ sigma merger ({err:.2f}% diff)")
        else:
            print(f"  EOB merger != sigma merger: {err:.1f}% (merger in plunge)")
        print()


# ═══════════════════════════════════════════════════════════════════════
# PROBLEM 3: SPINNING MASTER FORMULA
# ═══════════════════════════════════════════════════════════════════════

def _dfact(n):
    """Double factorial n!! = n(n-2)...; 0!! = (-1)!! = 1."""
    if n <= 0: return 1
    r = 1
    while n > 0: r *= n; n -= 2
    return r


def spin_free_coeff(n):
    """
    e_n^0 in E/(mu*c^2) = (-u/2) * [1 + sum e_n^0 * u^n]

    EXACT:  e_n^0 = -C(2n,n) * (3/4)^n * (2n-1)/(n+1)
    """
    return (-Fraction(comb(2*n,n), 4**n)
            * Fraction(3**n)
            * Fraction(2*n-1, n+1))


def spin_orbit_coeff(n):
    """
    e_n^{SO} in delta_E_SO/(mu*c^2) = chi * sum e_n^{SO} * u^{n+5/2}

    EXACT:  e_n^{SO} = -(2n+1)!! * 3^n / (2^n * n!)
    PN order of term: (2n+3)/2 PN  (n=0: 1.5PN, n=1: 2.5PN, ...)
    """
    return -Fraction(_dfact(2*n+1)*3**n, 2**n*factorial(n))


def spin_sq_coeff(n):
    """
    e_n^{SS} in delta2_E_SS/(mu*c^2) = chi^2 * sum e_n^{SS} * u^{n+3}

    EXACT:  e_n^{SS} = (2n+3)!! * 3^n / (2^{n+1} * n!)
    PN order of term: n+2 PN  (n=0: 2PN, n=1: 3PN, n=2: 4PN, ...)
    [BUG-FIX 2: was labeled n+3 PN, should be n+2 PN]
    """
    return Fraction(_dfact(2*n+3)*3**n, 2**(n+1)*factorial(n))


def kerr_exact_energy(u, chi, prograde=True):
    """Exact Kerr equatorial circular geodesic: E/mu = (1-2u+s*chi*u^{3/2})/sqrt(den) - 1."""
    s   = 1.0 if prograde else -1.0
    num = 1.0 - 2.0*u + s*chi*u**1.5
    den = 1.0 - 3.0*u + 2.0*s*chi*u**1.5
    return num/den**0.5 - 1.0 if den > 0 else float('nan')


def spinning_energy_series(u, chi, n_max=9):
    """PN series for Kerr binding energy (converges for u < 1/3)."""
    E = -u/2.0
    for n in range(1, n_max+1):
        E += float(spin_free_coeff(n)) * (-u/2.0) * u**n
    for n in range(n_max):
        E += float(spin_orbit_coeff(n)) * chi * u**(n+2.5)
    for n in range(n_max-1):
        E += float(spin_sq_coeff(n)) * chi**2 * u**(n+3)
    return E


def print_spinning_master_table():
    print("="*84)
    print("QGD SPINNING MASTER FORMULA  (Kerr geodesic -> three exact series)")
    print("="*84)

    # -- Spin-free --
    print("\n  SPIN-FREE:  e_n^0 = -C(2n,n)*(3/4)^n*(2n-1)/(n+1)")
    print(f"  {'n':>3}  {'e_n^0 exact':>22}  {'decimal':>12}  {'status':>24}")
    print("-"*67)
    known0 = {1:'1PN known',2:'2PN known',3:'3PN known',4:'4PN known'}
    for n in range(1,8):
        c0 = spin_free_coeff(n)
        st = known0.get(n,'* QGD PREDICTION *')
        print(f"  {n:>3}  {str(c0):>22}  {float(c0):>12.4f}  {st:>24}")

    # -- Spin-orbit --
    print("\n  SPIN-ORBIT: e_n^{SO} = -(2n+1)!!*3^n/(2^n*n!)")
    print(f"  {'n':>3}  {'PN order':>8}  {'e_n^{SO} exact':>22}  "
          f"{'decimal':>12}  {'status':>24}")
    print("-"*75)
    known_so = {0:'1.5PN OK',1:'2.5PN OK',2:'3.5PN OK',3:'4.5PN partial'}
    for n in range(8):
        cSO = spin_orbit_coeff(n)
        pn  = f'{2*n+3}/2 PN'
        st  = known_so.get(n,'* QGD PREDICTION *')
        print(f"  {n:>3}  {pn:>8}  {str(cSO):>22}  {float(cSO):>12.4f}  {st:>24}")

    # -- Spin-squared (BUG-FIX 2: n+2 PN not n+3 PN) --
    print("\n  SPIN-SQUARED: e_n^{SS} = (2n+3)!!*3^n/(2^{n+1}*n!)")
    print("  PN order = n+2  [u^{n+3} with u~(v/c)^2, leading energy ~u^1 = Newt]")
    print(f"  {'n':>3}  {'PN order':>8}  {'e_n^{SS} exact':>22}  "
          f"{'decimal':>12}  {'status':>24}")
    print("-"*75)
    known_ss = {0:'2PN OK',1:'3PN OK',2:'4PN OK'}
    for n in range(7):
        cSS = spin_sq_coeff(n)
        pn  = f'{n+2} PN'    # BUG-FIX 2: was n+3 PN
        st  = known_ss.get(n,'* QGD PREDICTION *')
        print(f"  {n:>3}  {pn:>8}  {str(cSS):>22}  {float(cSS):>12.4f}  {st:>24}")

    # -- Numerical verification --
    print("\n  VERIFICATION: series (n_max=9) vs exact Kerr at u=0.05")
    print(f"  {'chi':>5}  {'exact':>12}  {'series':>14}  {'error %':>10}")
    print("-"*50)
    for chi in [0.0, 0.3, 0.5, 0.7, 0.9]:
        E_ex = kerr_exact_energy(0.05, chi)
        E_s  = spinning_energy_series(0.05, chi, n_max=9)
        err  = abs(E_s-E_ex)/abs(E_ex)*100
        print(f"  {chi:>5.1f}  {E_ex:>12.6f}  {E_s:>14.6f}  {err:>10.4f}%")
    print("="*84)


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def run_all():
    print("\n" + "="*82)
    print("QGD THREE FOUNDATIONAL PROBLEMS  [v1.2 — 4 bugs corrected]")
    print("="*82)

    print("\n" + "-"*82)
    print("PROBLEM 2: MERGER CONDITION")
    print("-"*82)
    print_merger_resolution()

    print("\n" + "-"*82)
    print("PROBLEM 1: QGD EFFECTIVE ONE-BODY")
    print("-"*82)
    print("  A^QGD(u) = 1 - 2u/eta  [exact, zero free parameters]")
    print("  GR 3PN:   A^GR = 1 - 2u + 2*eta*u^3 + 18.69*eta*u^4")
    print("  Key diff: QGD has NO u^2, u^3, u^4 corrections. GR needs them.")
    print()

    # BUG-FIX 1: GW250114 = 33.6+32.2 Msun (not 86+77)
    systems = [
        ("GW150914",  36.0*M_sun,  29.0*M_sun),
        ("15+5 Msun", 15.0*M_sun,   5.0*M_sun),
        ("GW250114",  33.6*M_sun,  32.2*M_sun),   # [BUG-FIX 1] was 86+77
    ]
    for name, M1, M2 in systems:
        print(f"\n  [{name}]")
        eob = QGDEOB(M1, M2)
        eob.print_structure()
        eob.compare_table()

    print("\n" + "-"*82)
    print("PROBLEM 2 (cont.): IMR TRANSITION")
    print("-"*82)
    for name, M1, M2 in systems:
        imr = QGDMergerIMR(M1, M2)
        imr.print_imr()

    print("\n" + "-"*82)
    print("PROBLEM 3: SPINNING MASTER FORMULA")
    print("-"*82)
    print_spinning_master_table()

    print("\n" + "="*82)
    print("SUMMARY")
    print("="*82)
    print()
    print("  RIGOROUS:")
    print("  [M1] d_merge = 2GM1*(1+(M2/M1)^{1/3})^3/c^2  [exact, all q]")
    print("       Sigma(saddle)=1.000 verified; 8rs_indiv = 4rs_total at q=1")
    print("  [E1] A^QGD(u) = 1 - 2u/eta  [exact, no NR calibration]")
    print("       ISCO/LR/horizon structure identical to Schwarzschild scaled by eta")
    print("  [S1-3] Three exact closed-form series; all known GR coefficients matched")
    print()
    print("  OBSERVATIONAL (GW250114, SNR=80, Sep 2025):")
    print("  Kerr QNMs (2 modes) confirmed; no non-GR residual detected.")
    print("  QGD prediction 'ringdown = Kerr' is CONSISTENT with data.")
    print()
    print("  PRELIMINARY:")
    print("  [P1] eta-dependent 5PN+ corrections: not derived")
    print("  [P2] A^QGD calibration vs NR waveforms: not done")
    print("  [P3] LIGO event comparison: pending dipole-flux fix")
    print()
    print("  SHARPEST FALSIFIABLE PREDICTION:")
    print("  QGD: a2=a3=a4=0 in EOB A-function.")
    print("  GR:  a3=2*eta (u^3), a4~18.69*eta (u^4) from 3PN calculation.")
    print("  NR calibration should decide if these are zero or nonzero.")
    print("="*82)


if __name__ == "__main__":
    run_all()
