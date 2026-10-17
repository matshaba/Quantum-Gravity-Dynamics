#!/usr/bin/env python3
"""
wald_kerr.py -- Wald entropy of the l_Q^2 (nabla nabla sigma)^2 term on Kerr and Kerr-de Sitter horizons.   (G = hbar = c = 1)
================================================================================================================================

Non-rotating case (validation/wald_entropy.py): sigma is null and lies in the binormal plane of the horizon cross-section,
so the term contributes exactly zero.  With spin, sigma has a phi component and the term does not vanish.

KERR (ingoing Kerr-Schild form)
    ds^2 = -dv^2 + 2 dv dr - 2 a sin^2 th dr dphi + Sigma dth^2 + (r^2 + a^2) sin^2 th dphi^2 + (2 M r / Sigma) k k ,
    k = dv - a sin^2 th dphi,   sigma = sqrt(2 M r / Sigma) k,   horizon r_+ = M + sqrt(M^2 - a^2).

KERR-DE SITTER (Carter form, Lambda > 0)
    Delta_r = (r^2 + a^2)(1 - Lambda r^2/3) - 2 M r,  Delta_th = 1 + Lambda a^2 cos^2 th / 3,  Xi = 1 + Lambda a^2 / 3,
    Sigma = r^2 + a^2 cos^2 th, ingoing coordinates  dt = dv - (r^2 + a^2) dr / Delta_r,  dphi = dphi~ - a Xi dr / Delta_r
    (every 1/Delta_r term cancels, so the metric is regular at all roots of Delta_r: black-hole, cosmological and Nariai).
    The regular principal null form is k = dv - (a sin^2 th / Xi) dphi~.  ASSUMPTION: sigma = sqrt(2 M r / Sigma) k, the "mass"
    sigma of the recipe.  For a != 0 the Lambda part of the background is not flat + a rank-one term, so the recipe gives no
    clean sigma_Lambda there; the amplitude of sigma enters the Wald term linearly (it is quadratic in sigma), so this choice
    is a stated assumption.  At Lambda = 0 it is exactly the Kerr case.

The binormal of the cross-section v, r = const is eps = c dv ^ dr (conormals of the horizon cross-section), c fixed by
eps_ab eps^ab = -2.  Metric derivatives are symbolic; Christoffels and Riemann are assembled numerically at each theta, and
S = -2 pi Int E sqrt(h) dtheta dphi, with E the directional derivative of L along Y = eps eps.  Calibration at every horizon:
R/16 pi gives A/4 (A = 4 pi (r_h^2 + a^2)/Xi) and the Gauss-Bonnet term gives 4 pi alpha.
"""
import numpy as np
import sympy as sp

r, th = sp.symbols("r theta", real=True)
Ms, As, Ls = sp.symbols("M a Lambda", positive=True)


class Chart:
    """A (v, r, theta, phi) chart with a stationary axisymmetric metric g(r, theta; params) and a null 1-form k."""

    def __init__(self, g, H, k, params):
        self.params = params
        args = (r, th) + tuple(params)
        idx = [r, th]
        self._g = sp.lambdify(args, g, "numpy")
        self._dg = [sp.lambdify(args, sp.diff(g, x), "numpy") for x in idx]
        self._ddg = [[sp.lambdify(args, sp.diff(g, x, y), "numpy") for y in idx] for x in idx]
        self._H = sp.lambdify(args, H, "numpy")
        self._k = sp.lambdify(args, k, "numpy")

    def geometry(self, rv, tv, *p):
        g = np.array(self._g(rv, tv, *p), dtype=float)
        gi = np.linalg.inv(g)
        dg = np.zeros((4, 4, 4))
        ddg = np.zeros((4, 4, 4, 4))
        for i, c in enumerate([1, 2]):
            dg[c] = np.array(self._dg[i](rv, tv, *p), dtype=float)
            for j, d in enumerate([1, 2]):
                ddg[c, d] = np.array(self._ddg[i][j](rv, tv, *p), dtype=float)
        low = 0.5 * (np.einsum("cfb->fbc", dg) + np.einsum("bfc->fbc", dg) - dg)
        Gam = np.einsum("ef,fbc->ebc", gi, low)
        Riem = 0.5 * (np.einsum("bcad->abcd", ddg) + np.einsum("adbc->abcd", ddg)
                      - np.einsum("bdac->abcd", ddg) - np.einsum("acbd->abcd", ddg))
        Riem += np.einsum("ef,ebc,fad->abcd", g, Gam, Gam) - np.einsum("ef,ebd,fac->abcd", g, Gam, Gam)
        Ric = np.einsum("ac,abcd->bd", gi, Riem)
        Rs = np.einsum("bd,bd->", gi, Ric)
        return g, gi, Riem, Ric, Rs

    def integrands(self, rv, tv, *p, l2=1.0, sign=+1):
        """(E_EH, E_GB, E_sigma, sqrt(h)) at (r, theta)."""
        g, gi, Riem, Ric, Rs = self.geometry(rv, tv, *p)
        eps = np.zeros((4, 4))
        eps[0, 1], eps[1, 0] = 1.0, -1.0
        eps *= np.sqrt(-2.0 / np.einsum("ab,cd,ac,bd->", eps, eps, gi, gi))
        Y = np.einsum("ab,cd->abcd", eps, eps)
        Yu = np.einsum("abcd,ae,bf,cg,dh->efgh", Y, gi, gi, gi, gi)
        RicY = np.einsum("ac,abcd->bd", gi, Y)
        Ys = np.einsum("bd,bd->", gi, RicY)
        E_EH = Ys / (16 * np.pi)
        d_GB = 2 * np.einsum("abcd,abcd->", Riem, Yu) - 8 * np.einsum("ab,cd,ac,bd->", Ric, RicY, gi, gi) + 2 * Rs * Ys
        E_GB = d_GB / (16 * np.pi)
        H = self._H(rv, tv, *p)
        k = np.array(self._k(rv, tv, *p), dtype=float).ravel()
        sig_up = np.sqrt(H) * (gi @ k)

        def Q(A, B):
            a_ = np.einsum("mnab,n->mab", A, sig_up)
            w = np.einsum("lrkt,r->lkt", B, sig_up)
            W = np.einsum("ml,ak,bt,lkt->mab", gi, gi, gi, w)
            return np.einsum("mab,mab->", a_, W)
        E_sig = 2 * sign * l2 / 8.0 * Q(Riem, Y)
        sh = np.sqrt(g[2, 2] * g[3, 3] - g[2, 3] ** 2)
        return E_EH, E_GB, E_sig, sh

    def entropies(self, rh, *p, n=48, l2=1.0, sign=+1):
        """(S_EH, S_GB(alpha=1), S_sigma) on the cross-section r = rh, by Gauss-Legendre in theta."""
        x, w = np.polynomial.legendre.leggauss(n)
        ths, ws = 0.5 * np.pi * (x + 1), 0.5 * np.pi * w
        tot = np.zeros(3)
        for t, wt in zip(ths, ws):
            e = self.integrands(rh, t, *p, l2=l2, sign=sign)
            tot += np.array(e[:3]) * e[3] * wt
        return tuple(-2 * np.pi * 2 * np.pi * tot)


# ------------------------------------------------------------------------------------------------ Kerr
def _kerr_chart():
    Sig = r ** 2 + As ** 2 * sp.cos(th) ** 2
    k = sp.Matrix([1, 0, 0, -As * sp.sin(th) ** 2])
    g0 = sp.zeros(4, 4)
    g0[0, 0] = -1
    g0[0, 1] = g0[1, 0] = 1
    g0[1, 3] = g0[3, 1] = -As * sp.sin(th) ** 2
    g0[2, 2] = Sig
    g0[3, 3] = (r ** 2 + As ** 2) * sp.sin(th) ** 2
    H = 2 * Ms * r / Sig
    return Chart(g0 + H * k * k.T, H, k, (Ms, As))


# ------------------------------------------------------------------------------------------------ Kerr-de Sitter
def _kds_symbolic():
    x = sp.symbols("x", real=True)                               # x = cos(theta); metric depends on x^2 only
    Xi = 1 + Ls * As ** 2 / 3
    Dth = 1 + Ls * As ** 2 * x ** 2 / 3
    rho2 = r ** 2 + As ** 2 * x ** 2
    Dr = (r ** 2 + As ** 2) * (1 - Ls * r ** 2 / 3) - 2 * Ms * r
    s2 = 1 - x ** 2
    e1 = sp.Matrix([1, 0, 0, -As * s2 / Xi])
    e2 = sp.Matrix([As, 0, 0, -(r ** 2 + As ** 2) / Xi])
    gBL = -(Dr / rho2) * e1 * e1.T + (Dth * s2 / rho2) * e2 * e2.T
    gBL[1, 1] += rho2 / Dr
    gBL[2, 2] += rho2 / Dth
    J = sp.Matrix([[1, -(r ** 2 + As ** 2) / Dr, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, -As * Xi / Dr, 0, 1]])
    g = (J.T * gBL * J).applyfunc(sp.cancel)
    return g, x, Xi, rho2, Dr, s2


def _kds_chart():
    g, x, Xi, rho2, Dr, s2 = _kds_symbolic()
    g = g.subs(x, sp.cos(th))
    k = sp.Matrix([1, 0, 0, -As * sp.sin(th) ** 2 / Xi])
    H = 2 * Ms * r / (r ** 2 + As ** 2 * sp.cos(th) ** 2)
    return Chart(g, H, k, (Ms, As, Ls))


KERR = _kerr_chart()
_KDS = None


def kds_chart():
    global _KDS
    if _KDS is None:
        _KDS = _kds_chart()
    return _KDS


def horizon_entropies(M=1.0, a=0.0, n=48, l2=1.0, sign=+1):
    """Kerr: ((S_EH, S_GB, S_sigma), r_+)."""
    rp = M + np.sqrt(M ** 2 - a ** 2)
    return KERR.entropies(rp, M, a, n=n, l2=l2, sign=sign), rp


# ---- Kerr-de Sitter helpers
def kds_horizons(M, a, Lam):
    """Positive real roots of Delta_r (inner, black-hole, cosmological ...) sorted ascending."""
    c = [-Lam / 3.0, 0.0, 1.0 - Lam * a ** 2 / 3.0, -2.0 * M, a ** 2]
    rt = np.roots(c)
    return np.sort(rt[(np.abs(rt.imag) < 1e-9) & (rt.real > 0)].real)


def kds_nariai(a, Lam):
    """(r_N, M_N): double root of Delta_r(r) = (r^2 + a^2)(1 - Lam r^2/3) - 2 M r, where the black-hole and cosmological
    horizons coincide (Nariai limit with spin)."""
    from scipy.optimize import brentq

    def D0(x):
        return (x ** 2 + a ** 2) * (1 - Lam * x ** 2 / 3)

    def dD0(x):
        return 2 * x * (1 - Lam * x ** 2 / 3) - (x ** 2 + a ** 2) * 2 * Lam * x / 3
    F = lambda x: D0(x) - dD0(x) * x                              # Delta_r = 0 with M = dD0/2
    r0 = 1 / np.sqrt(Lam)
    grid = np.linspace(1.2 * r0, 0.05 * r0, 600)                  # first sign change scanning inward from above
    vals = [F(x) for x in grid]
    for i in range(len(grid) - 1):
        if vals[i] * vals[i + 1] < 0:
            rN = brentq(F, grid[i + 1], grid[i], xtol=1e-14)
            M = dD0(rN) / 2
            rts = kds_horizons(M * (1 - 1e-6), a, Lam)
            if len(rts) >= 2 and rts[-1] - rts[-2] < 0.05 * r0:
                return rN, M
            break
    raise ValueError("no Nariai double root with a*sqrt(Lambda) = %g" % (a * np.sqrt(Lam)))


def small_spin_coefficient(M, Lam, which="BH"):
    """c2 in S_sigma = c2 s l_Q^2 a^2 + O(a^4) at the static (a -> 0) horizon, fitted closed form
        c2 = -(4 pi^2/3) sigma_M^2(r_h) F'(r_h) / r_h ,   sigma_M^2 = r_s/r_h,   F = r_s/r + Lam r^2/3,  r_s = 2M.
    (sigma_M^2(r_h) = 1 - Lam r_h^2/3 is the amplitude of the mass sigma at the horizon; with Lam = 0 it is 1 and c2 = pi^2/(3 M^2).)"""
    rts = kds_horizons(M, 1e-9, Lam) if Lam > 0 else np.array([2 * M])
    rh = rts[-2] if which == "BH" and Lam > 0 else (rts[-1] if which == "cosm" else rts[-1])
    rs = 2 * M
    dF = -rs / rh ** 2 + 2 * Lam * rh / 3
    return -(4 * np.pi ** 2 / 3) * (rs / rh) * dF / rh


def main():
    print("KERR:  a/M   S_EH  pi(r+^2+a^2)   S_GB   S_sigma/(s l^2)")
    for a in (0.0, 0.3, 0.6, 0.9, 0.99):
        (seh, sgb, ssig), rp = horizon_entropies(1.0, a)
        print(f"{a:6.2f} {seh:10.5f} {np.pi*(rp**2+a**2):12.5f} {sgb:9.4f} {ssig:12.6f}")


if __name__ == "__main__":
    main()
