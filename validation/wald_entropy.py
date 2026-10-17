#!/usr/bin/env python3
"""
wald_entropy.py -- Wald entropy of the l_Q^2 (grad grad sigma)^2 term of the QGD action.
=========================================================================================

Units G = hbar = c = 1.  Metric (ingoing Eddington-Finkelstein / Kerr-Schild form, the form in which the QGD recipe
is written, g = eta + sigma sigma with sigma_mu = sqrt(F) k_mu, k = dv null):

        ds^2 = -(1 - F(r)) dv^2 + 2 dv dr + rho(r)^2 dOmega^2 ,       sigma_mu dx^mu = sqrt(F) dv .

    Schwarzschild  F = r_s / r            (sigma_t^2 = r_s/r)
    Schwarzschild-de Sitter  F = r_s/r + Lambda r^2 / 3      Nariai: double root of 1 - F.

Wald's formula (Iyer-Wald, derivatives of matter fields treated as independent variables):

        S = -2 pi  Int_H  E^{abcd} eps_ab eps_cd dA ,     E^{abcd} = dL / dR_{abcd} ,
        eps_ab = binormal of the horizon cross-section,  eps_ab eps^ab = -2 .

Because L is at most quadratic in the Riemann tensor, E eps eps is the directional derivative of L along
Y_{abcd} = eps_ab eps_cd (a decomposable 2-form, so Y has the Riemann symmetries: eps ^ eps = 0).  No index
bookkeeping of dL/dR is needed.  Calibration, done below: L = R/16 pi gives A/4; L = (alpha/16 pi) GB gives 4 pi alpha
per spherical horizon (Jacobson-Myers), the same value the Euclidean route gave in tests/test_nariai_entropy.py.

The sigma term.  L_sigma = s (l_Q^2/2) T_{ab m} T^{ab m},  T_{ab m} = nabla_a nabla_b sigma_m.  Split T into its part
symmetric in (a b), which does not involve curvature, and its antisymmetric part
T_[ab]m = (1/2) R_{m n ab} sigma^n, which does.  The cross term vanishes (symmetric x antisymmetric), so

        L_R = (s l_Q^2 / 8)  R_{m n a b} sigma^n  R^m_r^{ab} sigma^r ,

quadratic in Riemann and linear in the curvature after differentiation.  The sign s = +1 follows from the action as
printed (mostly-minus signature) after converting to mostly-plus: nabla nabla sigma nabla nabla sigma has three inverse
metrics and flips sign, R flips sign.  Kinetic term (1/2)(nabla sigma)^2 has no Riemann dependence and adds nothing.

Caveats (stated in the README as well): (i) g is treated as the dynamical field and sigma as matter on it, which is
the reading under which Wald's formula applies; QGD builds g from sigma, which this ignores.  (ii) sigma is regular on
the future horizon in these coordinates but the configuration has no bifurcation surface, so the result is the Iyer-Wald
integrand on a cross-section of the future horizon (stationarity makes it the same on every cross-section).
(iii) The relative normalisation of the (1/2)(nabla sigma)^2 and l_Q^2 terms to R/16 pi G is not fixed dimensionally in
the printed action; Planck units with coefficient exactly (l_Q^2/2) are used.
"""
import sympy as sp

r, v, th, ph = sp.symbols("r v theta phi", positive=True)
X = [v, r, th, ph]
F = sp.Function("F")(r)


def metric(F_expr, rho_expr=None):
    rho = r if rho_expr is None else rho_expr
    return sp.Matrix([[-(1 - F_expr), 1, 0, 0], [1, 0, 0, 0], [0, 0, rho ** 2, 0], [0, 0, 0, rho ** 2 * sp.sin(th) ** 2]])


def riemann_lower(g):
    n = 4
    gi = g.inv()
    Gam = [[[sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                 for d in range(n)) / 2 for c in range(n)] for b in range(n)] for a in range(n)]

    def R_up(a, b, c, d):
        e = sp.diff(Gam[a][b][d], X[c]) - sp.diff(Gam[a][b][c], X[d])
        e += sum(Gam[a][c][k] * Gam[k][b][d] - Gam[a][d][k] * Gam[k][b][c] for k in range(n))
        return e
    Ru = [[[[R_up(a, b, c, d) for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
    Rl = [[[[sp.simplify(sum(g[a, k] * Ru[k][b][c][d] for k in range(n))) for d in range(n)] for c in range(n)]
           for b in range(n)] for a in range(n)]
    return Rl, gi


def binormal():
    """eps = dv ^ dr (eps_{vr} = 1): binormal of the cross-section v, r = r_h, normalised eps_ab eps^ab = -2."""
    e = sp.zeros(4, 4)
    e[0, 1], e[1, 0] = 1, -1
    return e


def _raise_all(T, gi):
    n = 4
    return [[[[sum(gi[a, a2] * gi[b, b2] * gi[c, c2] * gi[d, d2] * T[a2][b2][c2][d2]
                   for a2 in range(n) for b2 in range(n) for c2 in range(n) for d2 in range(n)
                   if T[a2][b2][c2][d2] != 0)
               for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]


class Horizon:
    """Curvature and Wald integrands for  ds^2 = -(1-F) dv^2 + 2 dv dr + rho^2 dOmega^2  at a horizon r = r_h."""

    def __init__(self, F_expr=F, rho_expr=None):
        self.g = metric(F_expr, rho_expr)
        self.Rl, self.gi = riemann_lower(self.g)
        e = binormal()
        self.Y = [[[[e[a, b] * e[c, d] for d in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]
        self.F_expr = F_expr
        self.sigma_up = [0, sp.sqrt(F_expr), 0, 0]                   # sigma^r = sqrt(F), sigma^v = 0
        n = 4
        gi = self.gi
        self.Ric = sp.Matrix(n, n, lambda b, d: sum(gi[a, c] * self.Rl[a][b][c][d] for a in range(n) for c in range(n)))
        self.Rs = sum(gi[b, d] * self.Ric[b, d] for b in range(n) for d in range(n))

    # contractions ----------------------------------------------------------------------------------------------
    def _full(self, A, B):
        gi, n = self.gi, 4
        Bu = _raise_all(B, gi)
        return sum(A[a][b][c][d] * Bu[a][b][c][d] for a in range(n) for b in range(n) for c in range(n)
                   for d in range(n) if A[a][b][c][d] != 0)

    def _ricci_of(self, A):
        gi, n = self.gi, 4
        return sp.Matrix(n, n, lambda b, d: sum(gi[a, c] * A[a][b][c][d] for a in range(n) for c in range(n)))

    # Wald integrands E^{abcd} eps_ab eps_cd -------------------------------------------------------------------------
    def E_einstein(self):
        """L = R / 16 pi."""
        Ys = sum(self.gi[b, d] * self._ricci_of(self.Y)[b, d] for b in range(4) for d in range(4))
        return Ys / (16 * sp.pi)

    def E_gauss_bonnet(self, alpha=1):
        """L = (alpha / 16 pi) (Riem^2 - 4 Ric^2 + R^2): directional derivative along Y."""
        gi = self.gi
        RicY = self._ricci_of(self.Y)
        Ys = sum(gi[b, d] * RicY[b, d] for b in range(4) for d in range(4))
        ric_dot = sum(self.Ric[a, b] * RicY[c, d] * gi[a, c] * gi[b, d] for a in range(4) for b in range(4)
                      for c in range(4) for d in range(4))
        d_GB = 2 * self._full(self.Rl, self.Y) - 8 * ric_dot + 2 * self.Rs * Ys
        return alpha * d_GB / (16 * sp.pi)

    def _Q(self, A, B):
        """Q(A, B) = A_{m n a b} sigma^n  B^{m}_{r}^{ab} sigma^r  (symmetric bilinear)."""
        n, gi, s = 4, self.gi, self.sigma_up
        Bu = _raise_all(B, gi)
        tot = 0
        for m in range(n):
            for a in range(n):
                for b in range(n):
                    Am = sum(A[m][nn][a][b] * s[nn] for nn in range(n))
                    if Am == 0:
                        continue
                    # B^{m}_{r}^{ab} sigma^r = raise m,a,b ; r lowered contracted with sigma^r
                    Bm = sum(sum(self.g[rr, q] * Bu[m][q][a][b] for q in range(n)) * s[rr] for rr in range(n))
                    tot += Am * Bm
        return tot

    def E_sigma(self, l2=1, sign=1):
        """L = sign * (l_Q^2 / 8) Q(R, R)  ->  E = 2 * sign * (l_Q^2 / 8) Q(R, Y)."""
        return 2 * sp.Rational(sign) * l2 / 8 * self._Q(self.Rl, self.Y)


def wald_entropy(E_expr, r_h, rho_h=None):
    """S = -2 pi E * area of the cross-section (area = 4 pi rho_h^2); E evaluated at the horizon."""
    rho_h = r_h if rho_h is None else rho_h
    return sp.simplify(-2 * sp.pi * E_expr * 4 * sp.pi * rho_h ** 2)


def main():
    H = Horizon()
    rs, Lam = sp.symbols("r_s Lambda", positive=True)

    print("calibration at a Schwarzschild horizon (F = r_s/r, r_h = r_s)")
    sub = {F: rs / r}
    def at_h(E):
        return sp.simplify(E.subs(F, rs / r).doit().subs(r, rs))
    S_EH = wald_entropy(at_h(H.E_einstein()), rs)
    S_GB = wald_entropy(at_h(H.E_gauss_bonnet(1)), rs)
    print("  Einstein  S =", S_EH, "   (A/4 =", sp.pi * rs ** 2, ")")
    print("  Gauss-Bonnet alpha=1  S =", S_GB, "   (4 pi alpha expected)")

    print("\nsigma term, S_sigma / (s l_Q^2):")
    Es = H.E_sigma(1, 1)
    S_sig = wald_entropy(at_h(Es), rs)
    print("  Schwarzschild:", S_sig)
    return H


if __name__ == "__main__":
    main()
