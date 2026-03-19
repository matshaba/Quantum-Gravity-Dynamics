"""
QGD a_n GENERATOR — Standalone module
Generates a_n(nu) for any n, all sectors
Romeo Matshaba (UNISA) | March 2026
"""
from fractions import Fraction as F
import numpy as np
PI = np.pi

# ─── Ladder coefficients ─────────────────────────────────────────────
def c_A(n):
    """Ladder A: c_{n,0}^A / nu  (exact Fraction)"""
    return F(-41, 32) * F(5201, 656)**(n-4) if n >= 4 else F(0)

def c_B(n):
    """Ladder B: c_{n,0}^B / nu  (exact Fraction)"""
    return F(-63707, 1440) * F(5201, 656)**(n-5) if n >= 5 else F(0)

def c_Q(n):
    """Q-sector: c_{n,0}^Q / nu  (exact Fraction)"""
    return F(-41, 32) * F(-2275, 656)**(n-4) if n >= 4 else F(0)

# ─── Known exact rational ─────────────────────────────────────────────
RAT_EXACT = {3: F(2), 4: F(94, 3), 5: F(331054, 175)}
RAT2_EXACT = {5: F(-221, 6)}   # nu^2 coefficient at 4PN
LOG_EXACT  = {5: -22.0/3.0}    # log(u) coefficient * nu at 4PN

# ─── Rational fit (n >= 6) ────────────────────────────────────────────
_A_rat, _r_rat = 6.8752e-05, 30.7550

def c_rat_est(n):
    """Estimated c_{n,1}^rat for n >= 6"""
    return _A_rat * _r_rat**n

# ─── Critical poles ───────────────────────────────────────────────────
U_A_STAR = 656.0 / (5201.0 * PI**2)   # 0.012780
U_Q_STAR = 656.0 / (2275.0 * PI**2)   # 0.029216
U_R_STAR = 1.0 / _r_rat                 # 0.032515
U_ISCO   = 1.0 / 6.0

# ─── Main generator ───────────────────────────────────────────────────
def a_n_trans(n, nu=0.25):
    """Transcendental part of a_n at given nu (float)"""
    if n < 4: return 0.0
    cA = float(c_A(n)) * PI**(2*(n-3)) * nu
    cB = float(c_B(n)) * PI**(2*(n-4)) * nu if n >= 5 else 0.0
    return cA + cB

def a_n_rat(n, nu=0.25):
    """Rational part of a_n at given nu (exact where known, est otherwise)"""
    if n in RAT_EXACT:
        base = float(RAT_EXACT[n]) * nu
    elif n >= 6:
        base = c_rat_est(n) * nu
    else:
        base = 0.0
    if n in RAT2_EXACT:
        base += float(RAT2_EXACT[n]) * nu**2
    return base

def a_n(n, nu=0.25):
    """Full a_n coefficient (Schwarzschild + trans + rat + log)"""
    schw = -2.0 if n == 1 else 0.0
    return schw + a_n_trans(n, nu) + a_n_rat(n, nu)

def A_closed(u, nu=0.25):
    """Closed-form A(u;nu) using generating functions"""
    xi_A = 5201*PI**2/656
    d_A  = 1.0 - xi_A*u
    d_nu = 1.0 - (1.0-4*nu)*nu
    if abs(d_A) < 1e-10: return float('inf')
    seed = (41.0*PI**2/32.0)*nu*u**4
    A_tr  = -seed * (1 + 63707*32/(41*1440)*u) / (d_A * d_nu)
    A_34  = (nu*(94.0/3 - 41*PI**2/32)*u**4
             + nu*(331054.0/175 - 63707*PI**2/1440 - 5201*PI**4/512)*u**5
             + nu**2*(41*PI**2/32 - 221.0/6)*u**5)
    return 1 - 2*u + 2*nu*u**3 + A_34 + A_tr

def Q_closed(u, nu=0.25):
    """Closed-form Q_trans(u;nu)"""
    xi_Q = 2275*PI**2/656
    d_Q  = 1.0 + xi_Q*u
    d_nu = 1.0 - (1.0-4*nu)*nu
    return (41*PI**2/32)*nu*u**4 / (d_Q * d_nu)

# ─── Print table ──────────────────────────────────────────────────────
if __name__ == "__main__":
    print(f"\n{'='*80}")
    print("  QGD a_n GENERATOR — Romeo Matshaba (UNISA) | March 2026")
    print(f"{'='*80}")
    print(f"\n  Critical thresholds:")
    print(f"    u_A* = {U_A_STAR:.6f}  (trans pole, Ladders A+B)")
    print(f"    u_Q* = {U_Q_STAR:.6f}  (Q-sector pole)")
    print(f"    u_R* = {U_R_STAR:.6f}  (rational pole, fitted)")
    print(f"    u_ISCO = {U_ISCO:.6f}")
    print()
    print(f"  {'n':>3}  {'PN':>4}  {'c_A/nu':>20}  {'c_B/nu':>20}  "
          f"{'a_n(nu=1/4)':>16}  {'status':>8}")
    print(f"  {'-'*82}")
    for n in range(1, 31):
        cA_str = str(c_A(n)) if len(str(c_A(n))) < 18 else f"{float(c_A(n)):.4e}"
        cB_str = str(c_B(n)) if len(str(c_B(n))) < 18 else f"{float(c_B(n)):.4e}"
        an     = a_n(n, 0.25)
        st     = "exact" if n <= 5 else ("est" if n >= 6 else "")
        if n == 1: cA_str = cB_str = "0"
        print(f"  {n:>3}  {n-1:>3}PN  {cA_str:>20}  {cB_str:>20}  {an:>16.6e}  {st:>8}")

    print(f"\n  A(u; nu=1/4) numerical:")
    print(f"  {'u':>8}  {'A_closed':>14}  {'Domain':>6}")
    for u in [0.001, 0.005, U_A_STAR*0.95, U_A_STAR*1.05, 1/6]:
        Av = A_closed(u, 0.25)
        dom = "I" if u < U_A_STAR else ("II" if u < U_R_STAR else "III+")
        print(f"  {u:>8.5f}  {Av:>14.8f}  {dom:>6}")

# ─── Rational sector closed form (Pade[1/1]) ─────────────────────────
_p0, _p1, _q1 = 2.0, -89.4158, 60.3746   # from c3,c4,c5 exact fit
U_RAT_PADE = 1.0 / _q1   # 0.01656

def R_rat_pade(u, nu):
    """Rational sector Pade[1/1]: exact fit to c3,c4,c5"""
    if 1.0 - _q1*u <= 0: return float('nan')
    return nu * u**3 * (_p0 + _p1*u) / (1.0 - _q1*u)

def A_complete(u, nu=0.25):
    """
    Full A(u;nu): Schwarzschild + 2PN + exact 3-4PN
                 + A_trans (all orders, closed) + R_rat (Pade[1/1]) + log
    Valid for u < u_A* = 0.01278.  Beyond: use EOB Pade resummation.
    """
    import numpy as np
    PI = np.pi
    xi_A = 5201*PI**2/656
    d_A  = 1.0 - xi_A*u
    d_nu = 1.0 - (1.0-4*nu)*nu
    if abs(d_A) < 1e-8: return float('nan')
    A = 1 - 2*u + 2*nu*u**3
    A += nu*(94/3 - 41/32*PI**2)*u**4
    A += (nu*(331054/175 - 63707/1440*PI**2 - 5201/512*PI**4)
          + nu**2*(41/32*PI**2 - 221/6))*u**5
    # Transcendental sector (all orders)
    A += -(41*PI**2/32)*nu*u**4*(1+63707*32/(41*1440)*u) / (d_A*d_nu)
    # Rational sector Pade[1/1]
    A += R_rat_pade(u, nu)
    # Log sector
    if u > 0: A += (-22/3)*nu*u**5*np.log(u)
    return A

if __name__ == "__main__" and False:   # guard: don't re-run table
    pass
