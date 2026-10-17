#!/usr/bin/env python3
"""
classical_limit_entropy.py -- QGD's "hbar -> 0" against the high-entropy limit of Carcassi, Landini, Aidala (arXiv:2411.00972).
===========================================================================================================================

That paper: classical mechanics is the limit S >> 0 of quantum mechanics; hbar -> 0 is a formal stand-in for it.  For Gaussian
states with uncertainty (action) Sigma, x = Sigma/hbar:

        S_C(x) = ln x + 1                                   (classical Gibbs/Shannon entropy, Sigma/h with h = 2 pi hbar)
        S_Q(x) = (x+1/2) ln(x+1/2) - (x-1/2) ln(x-1/2)      (von Neumann)
        stretch hbar -> hbar/lam:   S_lamQ(x) = S_Q(lam x) - ln lam  ->  S_C.

What this module checks (all exact or numerical, nothing fitted):
  1. S_Q - S_C = -1/(24 x^2) - 1/(320 x^4) - ...: the quantum correction to the entropy is O(hbar^2), the same order as every QGD correction.
  2. The stretch map divides that correction by lam^2, exactly as hbar -> hbar/lam divides l_Q^2 = O(hbar^2).
  3. For a mass M localised within r the action scale is Sigma ~ M c r, so x = r/lambda_C.  At the Schwarzschild radius
     x = r_s/lambda_C = S_BH/(2 pi):  QGD's expansion parameter at the horizon is eps_h = (lambda_C/r_s)^2 = (2 pi/S_BH)^2.
  4. The classical-entropy description breaks down (S_C < 0) for x < 1/e, i.e. r < lambda_C/e.

Not shown: that QGD's field equations have a high-entropy limit, or that Sigma ~ M c r is the right action scale.
"""
import numpy as np
import sympy as sp

x, lam = sp.symbols("x lambda", positive=True)
G, M, hbar, c = sp.symbols("G M hbar c", positive=True)


def S_C(xv):
    return np.log(xv) + 1.0


def S_Q(xv):
    xv = np.asarray(xv, dtype=float)
    up, dn = xv + 0.5, xv - 0.5
    return up * np.log(up) - np.where(dn > 0, dn * np.log(np.where(dn > 0, dn, 1.0)), 0.0)


def S_stretched(xv, lamv):
    return S_Q(lamv * np.asarray(xv, dtype=float)) - np.log(lamv)


def entropy_difference_series(order=7):
    """Series of S_Q - S_C at large x (exact, sympy).  With u = 1/x the ln x terms cancel analytically:
    S_Q - ln x = (1/u + 1/2) ln(1 + u/2) - (1/u - 1/2) ln(1 - u/2)."""
    u = sp.symbols("u", positive=True)
    expr = (1 / u + sp.Rational(1, 2)) * sp.log(1 + u / 2) - (1 / u - sp.Rational(1, 2)) * sp.log(1 - u / 2) - 1
    return sp.series(expr, u, 0, order).removeO().subs(u, 1 / x)


def horizon_ratio():
    """r_s / lambda_C and S_BH as functions of (G, M, hbar, c)."""
    r_s = 2 * G * M / c ** 2
    lamC = hbar / (M * c)
    S_BH = (4 * sp.pi * r_s ** 2) / (4 * G * hbar / c ** 3)
    return sp.simplify(r_s / lamC), sp.simplify(S_BH)


def entropy_zero_crossings():
    """x at which the classical and quantum Gaussian entropies vanish: (1/e, 1/2)."""
    return float(np.exp(-1.0)), 0.5


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pathlib import Path
    print("S_Q - S_C series:", entropy_difference_series())
    ratio, sbh = horizon_ratio()
    print("r_s/lambda_C =", ratio, "   S_BH =", sbh, "   ratio/(S_BH/2pi) =", sp.simplify(ratio / (sbh / (2 * sp.pi))))
    for xv in (0.6, 1, 2, 5, 20, 100):
        print("x=%6.1f  S_C=%.4f  S_Q=%.4f  diff=%.3e  -1/(24x^2)=%.3e" % (xv, S_C(xv), S_Q(xv), S_Q(xv) - S_C(xv), -1 / (24 * xv ** 2)))
    fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.4))
    xs = np.logspace(-0.3, 2, 300)
    ax[0].loglog(xs, np.abs(S_Q(xs) - S_C(xs)), "C0", label=r"$|S_Q-S_C|$")
    for l in (1, 3, 10):
        ax[0].loglog(xs, np.abs(S_stretched(xs, l) - S_C(xs)), "--", label=r"stretched, $\lambda=%d$" % l)
    ax[0].loglog(xs, 1 / (24 * xs ** 2), "k:", label=r"$1/24x^2$")
    ax[0].set_xlabel(r"$x=\Sigma/\hbar$"); ax[0].set_ylabel("entropy difference"); ax[0].legend(fontsize=7)
    Sb = np.logspace(0, 4, 100)
    ax[1].loglog(Sb, (2 * np.pi / Sb) ** 2, "C3")
    ax[1].axvline(2 * np.pi, color="gray", lw=0.8)
    ax[1].set_xlabel(r"$S_{BH}$"); ax[1].set_ylabel(r"$\varepsilon_h=(\lambda_C/r_s)^2$")
    ax[1].set_title("QGD expansion parameter at the horizon", fontsize=9)
    fig.tight_layout()
    out = Path(__file__).resolve().parent / "classical_limit_entropy.png"
    fig.savefig(out, dpi=130)
    print("figure ->", out)


if __name__ == "__main__":
    main()
