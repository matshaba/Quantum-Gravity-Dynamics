#!/usr/bin/env python3
"""
kerr_quantum_term.py -- the literal master-metric quantum term on a rotating (Kerr) source.
==========================================================================================

Setup.  Kerr in ingoing Kerr-Schild form, g = eta + H k k, H = 2Mr/Sigma, k = dv - a sin^2(th) dphi (+ dr),
sigma_mu = sqrt(H) k_mu.  The quantum term of the master metric is  delta g = s * b * C,  b = kappa l_Q^2,
C_ab = eta^{cd} d_a sigma_c d_b sigma_d.

Result (verified symbolically in the tests, by pulling the Cartesian tensor back to spheroidal coordinates):

        delta g = s b (2Mr/Sigma) (d theta^2 + sin^2(th) d phi^2)             (angular only, as for Schwarzschild)

Because Kerr has g^{r phi} = a/Sigma != 0, an angular-only perturbation still moves the horizon.  To first order in b the
horizon r = r_h(theta) is the null surface g^{rr} = 0 (g^{r theta} = 0, so the r_h' term is second order):

        delta r_h(theta) = s b a^2 sin^2(th) (r+^2 + a^2) / (2 Sigma^2 (r+ - M)).

The induced area and the shift  Delta S = Delta A / 4  (G = hbar = 1) follow by integrating over theta.
"""
import numpy as np
import sympy as sp

v, r, th, ph = sp.symbols("v r theta phi", real=True)
M, a, b, s = sp.symbols("M a b s", real=True)

Sigma = r ** 2 + a ** 2 * sp.cos(th) ** 2
Delta = r ** 2 - 2 * M * r + a ** 2
H = 2 * M * r / Sigma


def kerr_metric():
    """Ingoing Kerr-Schild Kerr, coordinates (v, r, theta, phi)."""
    g = sp.zeros(4, 4)
    S2 = sp.sin(th) ** 2
    g[0, 0] = -(1 - H)
    g[0, 1] = g[1, 0] = 1
    g[0, 3] = g[3, 0] = -H * a * S2
    g[1, 3] = g[3, 1] = -a * S2
    g[2, 2] = Sigma
    g[3, 3] = S2 * ((r ** 2 + a ** 2) ** 2 - Delta * a ** 2 * S2) / Sigma
    return g


def delta_metric():
    d = sp.zeros(4, 4)
    d[2, 2] = H
    d[3, 3] = H * sp.sin(th) ** 2
    return d


def perturbed_metric():
    return kerr_metric() + s * b * delta_metric()


def r_plus(Mv, av):
    return Mv + np.sqrt(Mv ** 2 - av ** 2)


def delta_rh_formula():
    rp = M + sp.sqrt(M ** 2 - a ** 2)
    return s * b * a ** 2 * sp.sin(th) ** 2 * (rp ** 2 + a ** 2) / (2 * Sigma.subs(r, rp) ** 2 * (rp - M))


def delta_rh_from_inverse():
    """delta r_h from g^{rr}=0 using the first-order inverse (independent of the closed form above)."""
    g0 = kerr_metric()
    gi = g0.inv()
    dgi = -gi * (s * b * delta_metric()) * gi
    grr0 = sp.simplify(gi[1, 1])
    dgrr = sp.simplify(dgi[1, 1])
    rp = M + sp.sqrt(M ** 2 - a ** 2)
    return sp.simplify((-dgrr / sp.diff(grr0, r)).subs(r, rp))


_INV = {}


def _inverse_funcs():
    """Lambdified g^{rr} and g^{r phi}/g^{r v} of the perturbed metric (symbolic inverse, built once)."""
    if not _INV:
        gi = perturbed_metric().inv()
        args = (r, M, a, b, s, th)
        _INV["grr"] = sp.lambdify(args, gi[1, 1], "numpy")
        _INV["om"] = sp.lambdify(args, gi[3, 1] / gi[0, 1], "numpy")
    return _INV


def exact_horizon_radius(Mv, av, bv, thv, sv=1.0):
    """Root of the full g^{rr}(r, theta) = 0 near r_+ for finite b (numerical, no linearisation)."""
    from scipy.optimize import brentq
    grr = _inverse_funcs()["grr"]
    rp = r_plus(Mv, av)
    lo, hi = rp - 0.5 * (rp - Mv), rp + 0.5 * (rp - Mv) + 0.2
    return brentq(lambda x: float(grr(x, Mv, av, bv, sv, thv)), lo, hi)


def area_shift_integrand():
    """First-order d(sqrt(h)) / d(theta dphi) of the horizon 2-surface, h = induced metric on r = r_h(theta)."""
    g = perturbed_metric()
    drh = delta_rh_formula()
    rp = M + sp.sqrt(M ** 2 - a ** 2)
    hthth = g[2, 2]
    hphph = g[3, 3]
    sq = sp.sqrt(hthth * hphph)                                   # g_rr = g_r theta = 0; cross term h_theta phi is O(b) -> O(b^2) in det
    sq0 = sq.subs(b, 0)
    first = sp.diff(sq, b).subs(b, 0) + sp.diff(sq0, r) * drh / b
    return sp.simplify(first.subs(r, rp) * b)


def area_shift(Mv, av, sv=1.0, bv=1.0):
    """Delta A at first order in b (G = hbar = c = 1), with M, a numeric."""
    from scipy.integrate import quad
    f = sp.lambdify(th, area_shift_integrand().subs({M: Mv, a: av, s: sv, b: bv}), "numpy")
    return 2 * np.pi * quad(lambda x: float(f(x)), 0, np.pi, epsabs=1e-13, epsrel=1e-12)[0]


def entropy_shift_over_pi_b(Mv, av, sv=1.0):
    """Delta S / (pi s b) with Delta S = Delta A / 4; equals 1 for Schwarzschild."""
    return area_shift(Mv, av, sv, 1.0) / (4 * np.pi * sv)


def angular_velocity_profile(Mv, av, bv, thv, sv=1.0):
    """Omega(theta) = g^{r phi} / g^{r v} on the perturbed horizon (generator direction), numeric, finite b."""
    om = _inverse_funcs()["om"]
    return np.array([float(om(exact_horizon_radius(Mv, av, bv, t, sv), Mv, av, bv, sv, t)) for t in thv])


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pathlib import Path
    print("delta r_h closed form == from inverse metric:",
          sp.simplify(delta_rh_formula() - delta_rh_from_inverse()) == 0)
    ratios = np.linspace(0.0, 0.9, 19)
    vals = [entropy_shift_over_pi_b(1.0, x) for x in ratios]
    for x, y in zip(ratios[::3], vals[::3]):
        print("a/M = %.2f   Delta S/(pi s b) = %.6f" % (x, y))
    thv = np.linspace(0.2, np.pi - 0.2, 5)
    om0 = angular_velocity_profile(1.0, 0.6, 0.0, thv)
    om1 = angular_velocity_profile(1.0, 0.6, 1e-3, thv)
    print("Omega(theta), b=0   :", om0)
    print("Omega(theta), b=1e-3:", om1, " spread", np.ptp(om1))
    fig, ax = plt.subplots(1, 2, figsize=(8, 3.3))
    ax[0].plot(ratios, vals, "C3-o", ms=3)
    ax[0].axhline(1, color="gray", lw=0.8)
    ax[0].set_xlabel("a/M"); ax[0].set_ylabel(r"$\Delta S/(\pi s\,\kappa\ell_Q^2)$"); ax[0].set_title("entropy shift of the literal term")
    tt = np.linspace(0, np.pi, 200)
    for x in (0.3, 0.6, 0.9):
        rp = r_plus(1.0, x)
        ax[1].plot(tt, x ** 2 * np.sin(tt) ** 2 * (rp ** 2 + x ** 2) / (2 * (rp ** 2 + x ** 2 * np.cos(tt) ** 2) ** 2 * (rp - 1)),
                   label="a/M=%.1f" % x)
    ax[1].set_xlabel(r"$\theta$"); ax[1].set_ylabel(r"$\delta r_h/(s b)$"); ax[1].legend(fontsize=7)
    ax[1].set_title("horizon displacement")
    fig.tight_layout()
    out = Path(__file__).resolve().parent / "kerr_quantum_term.png"
    fig.savefig(out, dpi=130)
    print("figure ->", out)


if __name__ == "__main__":
    main()
