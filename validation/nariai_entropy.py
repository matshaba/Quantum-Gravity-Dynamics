#!/usr/bin/env python3
"""
nariai_entropy.py -- test QGD's black-hole entropy correction on the Nariai / de Sitter pair.
=============================================================================================

Units: G = hbar = c = k_B = 1 throughout (Planck units), so S = A/4 and Lambda has units 1/length^2.
The only dimensionful input is Lambda; every other quantity is a pure number times a power of Lambda or of A.

What QGD claims (GR-vs-QGD.md):   S = A/4 + 2 kappa l_Q^2 ,   "constant, mass independent, topological".
What the metric says (README):     f = r_s/r + kappa hbar^2/(M^2 r^2),   g_tt = -(1-f).

PART 1  Euclidean route (no horizon areas used).  For an Einstein space R_ab = Lambda g_ab the on-shell action is
        I = -(1/16 pi) Int (R - 2 Lambda) sqrt(g),  S = -I.
        de Sitter  = S^4,      S_dS = 3 pi / Lambda      (= one horizon of area 12 pi / Lambda  ->  A/4)
        Nariai     = S^2 x S^2, S_N  = 2 pi / Lambda      (= two coincident horizons, each 4 pi/Lambda -> A/4)
        This is compared with the horizon-area route, which is a different calculation.
PART 2  Gauss-Bonnet density integrated over the same manifolds: Int GB = 32 pi^2 chi, with chi(S^4) = 2 and
        chi(S^2 x S^2) = 4.  A TOPOLOGICAL entropy term is proportional to chi, i.e. to the NUMBER OF HORIZONS,
        so it is blind to dS/dM at fixed horizon count but NOT blind to a process that changes it
        (de Sitter -> Nariai pair creation).
PART 3  The quantum term as a metric term, two readings of the README formula:
        (I)  A = kappa * l_Q^2 = const          -> area-law shift  2 pi A  (a constant, but pi x the claimed 2 kappa l_Q^2,
                                                    and it creates a horizon at M = 0)
        (II) A = kappa / M^2   (hbar^2/(M^2 c^2)) -> area-law shift  2 pi kappa / M^2  (mass dependent, NOT constant)
"""
from pathlib import Path

import numpy as np
import sympy as sp

OUT = Path(__file__).resolve().parent / "nariai_entropy.png"


# ---------------------------------------------------------------- symbolic curvature (independent of QGD code)
def curvature(g, X):
    """Christoffel -> Riemann R^a_{bcd} -> Ricci, scalar, Gauss-Bonnet density for metric g in coordinates X."""
    n = len(X)
    gi = g.inv()
    Gam = [[[sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                 for d in range(n)) / 2 for c in range(n)] for b in range(n)] for a in range(n)]

    def Riem(a, b, c, d):
        e = sp.diff(Gam[a][b][d], X[c]) - sp.diff(Gam[a][b][c], X[d])
        e += sum(Gam[a][c][k] * Gam[k][b][d] - Gam[a][d][k] * Gam[k][b][c] for k in range(n))
        return e
    R4 = [[[[sp.simplify(Riem(a, b, c, d)) for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
    Ric = sp.Matrix(n, n, lambda b, d: sp.simplify(sum(R4[a][b][a][d] for a in range(n))))
    Rs = sp.simplify(sum(gi[b, d] * Ric[b, d] for b in range(n) for d in range(n)))
    # all-lower Riemann, then contractions (metric is diagonal in both examples)
    Rl = [[[[sum(g[a, k] * R4[k][b][c][d] for k in range(n)) for d in range(n)] for c in range(n)]
           for b in range(n)] for a in range(n)]
    Riem2 = sp.simplify(sum(Rl[a][b][c][d] * Rl[a2][b2][c2][d2] * gi[a, a2] * gi[b, b2] * gi[c, c2] * gi[d, d2]
                            for a in range(n) for b in range(n) for c in range(n) for d in range(n)
                            if Rl[a][b][c][d] != 0
                            for a2 in [a] for b2 in [b] for c2 in [c] for d2 in [d]))
    Ric2 = sp.simplify(sum(Ric[a, b] * Ric[a2, b2] * gi[a, a2] * gi[b, b2]
                           for a in range(n) for b in range(n) for a2 in [a] for b2 in [b]))
    GB = sp.simplify(Riem2 - 4 * Ric2 + Rs ** 2)
    return Ric, Rs, GB


def euclidean_data():
    """S^4 and S^2 x S^2: Einstein constant, Euclidean entropy -I, Euler characteristic from Int GB / 32 pi^2."""
    a = sp.symbols("a", positive=True)
    th = sp.symbols("t1 t2 t3 t4", positive=True)
    out = {}

    # S^2 x S^2, radius a
    t1, p1, t2, p2 = sp.symbols("t1 p1 t2 p2", positive=True)
    g = sp.diag(a ** 2, a ** 2 * sp.sin(t1) ** 2, a ** 2, a ** 2 * sp.sin(t2) ** 2)
    sg = a ** 4 * sp.sin(t1) * sp.sin(t2)
    ranges = [(t1, 0, sp.pi), (p1, 0, 2 * sp.pi), (t2, 0, sp.pi), (p2, 0, 2 * sp.pi)]
    out["S2xS2"] = _manifold(g, [t1, p1, t2, p2], sg, ranges, a)

    # S^4, radius a:  d al^2 + sin^2 al ( d be^2 + sin^2 be ( d ga^2 + sin^2 ga d ph^2 ) )
    al, be, ga, ph = sp.symbols("al be ga ph", positive=True)
    g4 = sp.diag(a ** 2, a ** 2 * sp.sin(al) ** 2, a ** 2 * sp.sin(al) ** 2 * sp.sin(be) ** 2,
                 a ** 2 * sp.sin(al) ** 2 * sp.sin(be) ** 2 * sp.sin(ga) ** 2)
    sg4 = a ** 4 * sp.sin(al) ** 3 * sp.sin(be) ** 2 * sp.sin(ga)
    r4 = [(al, 0, sp.pi), (be, 0, sp.pi), (ga, 0, sp.pi), (ph, 0, 2 * sp.pi)]
    out["S4"] = _manifold(g4, [al, be, ga, ph], sg4, r4, a)
    return out


def _manifold(g, X, sqrtg, ranges, a):
    Ric, Rs, GB = curvature(g, X)
    Lam = sp.simplify(Ric[0, 0] / g[0, 0])                      # R_ab = Lam g_ab ?
    einstein = all(sp.simplify(Ric[i, j] - Lam * g[i, j]) == 0 for i in range(4) for j in range(4))
    vol = sqrtg
    for lim in ranges:
        vol = sp.integrate(vol, lim)
    vol = sp.simplify(vol)
    I = sp.simplify(-(Rs - 2 * Lam) * vol / (16 * sp.pi))
    chi = sp.simplify(GB * vol / (32 * sp.pi ** 2))
    return dict(einstein=einstein, Lam=Lam, R=Rs, volume=vol, S=sp.simplify(-I), chi=chi,
                S_over_inv_Lam=sp.simplify(-I * Lam), a=a)


# ---------------------------------------------------------------- horizon (area) route with the quantum term
def horizons(M, A, Lam=0.0):
    """Positive roots of r^2 h = r^2 - 2 M r - A - Lam r^4 / 3 = 0  (h = 1 - 2M/r - A/r^2 - Lam r^2/3)."""
    c = [-Lam / 3.0, 0.0, 1.0, -2.0 * M, -A]
    r = np.roots(c)
    r = np.sort(r[(np.abs(r.imag) < 1e-12) & (r.real > 0)].real)
    return r


def area_entropy(M, A):
    """pi r_+^2 for the outer horizon, Lambda = 0."""
    return np.pi * horizons(M, A)[-1] ** 2


def temperature(M, A):
    r = horizons(M, A)[-1]
    return (2 * M / r ** 2 + 2 * A / r ** 3) / (4 * np.pi)


def nariai_with_A(A, Lam=1.0):
    """Double root of r^2 h with the quantum term (Lambda fixed): returns (r_N, M_N)."""
    r, M = sp.symbols("r M", positive=True)
    P = r ** 2 - 2 * M * r - A - sp.Float(Lam) * r ** 4 / 3
    Mf = sp.solve(sp.diff(P, r), M)[0]                          # P' = 0  -> M(r)
    Pr = sp.simplify(P.subs(M, Mf))
    guess = 1.0 / np.sqrt(Lam)
    rN = sp.nsolve(Pr, r, guess)
    return float(rN), float(Mf.subs(r, rN))


def main():
    print("PART 1/2  Euclidean route")
    d = euclidean_data()
    for k in ("S4", "S2xS2"):
        v = d[k]
        print(f"  {k:6s} Einstein={v['einstein']}  Lambda={v['Lam']}  S*Lambda={v['S_over_inv_Lam']}  chi={v['chi']}")

    print("\nPART 3  quantum term as a metric term (G=hbar=c=1)")
    kappa = 1.0
    Ms = np.logspace(0.3, 2.0, 40)
    A_I = 0.05                                                  # reading I: constant A
    dS_I = np.array([area_entropy(M, A_I) - 4 * np.pi * M ** 2 for M in Ms])
    dS_II = np.array([area_entropy(M, kappa / M ** 2) - 4 * np.pi * M ** 2 for M in Ms])
    print(f"  reading I  (A=const={A_I}):  area shift -> {dS_I[-1]:.5f}   2*pi*A = {2*np.pi*A_I:.5f}")
    print(f"  reading II (A=kappa/M^2):    area shift * M^2 -> {dS_II[-1]*Ms[-1]**2:.5f}   2*pi*kappa = {2*np.pi*kappa:.5f}")
    print(f"  reading I at M=0: horizons at r = {horizons(0.0, A_I)} (sqrt(A) = {np.sqrt(A_I):.5f})")

    # first-law check: does dM = T dS_area hold?  (fixed A: exact for RN-like; varying A: not)
    M0, h = 5.0, 1e-5
    for name, Af in (("I ", lambda M: A_I), ("II", lambda M: kappa / M ** 2)):
        dS = (area_entropy(M0 + h, Af(M0 + h)) - area_entropy(M0 - h, Af(M0 - h))) / (2 * h)
        print(f"  first law, reading {name}: (dS_area/dM) * T = {dS * temperature(M0, Af(M0)):.8f}  (1 means dM = T dS)")

    fig, ax = plt_figs(Ms, dS_I, dS_II, A_I, kappa)
    return d


def plt_figs(Ms, dS_I, dS_II, A_I, kappa):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
    ax = axs[0]
    ax.loglog(Ms, dS_II, label=r"reading II: $A=\kappa/M^2$")
    ax.loglog(Ms, 2 * np.pi * kappa / Ms ** 2, "k:", label=r"$2\pi\kappa/M^2$")
    ax.axhline(2 * np.pi * A_I, color="C1", ls="--", label=r"reading I: $A$ const $\to 2\pi A$")
    ax.set_xlabel("M  (Planck units)"); ax.set_ylabel(r"$\pi r_+^2-4\pi M^2$")
    ax.set_title("area-law shift from the quantum term"); ax.legend(fontsize=7)
    ax = axs[1]
    Mgrid = np.linspace(0, 3, 200)
    ax.plot(Mgrid, [horizons(M, A_I)[-1] for M in Mgrid], label=r"$r_+$ (A const)")
    ax.plot(Mgrid, 2 * Mgrid, "k--", label=r"$r_s=2M$ (GR)")
    ax.set_xlabel("M"); ax.set_ylabel("r"); ax.set_title(r"reading I: horizon at $M=0$ ($r=\sqrt{A}$)"); ax.legend(fontsize=7)
    ax = axs[2]
    labels = [r"$S^4$ (dS)", r"$S^2\!\times\!S^2$ (Nariai)"]
    ax.bar(labels, [3 * np.pi, 2 * np.pi], color=["C0", "C2"])
    ax.set_ylabel(r"$S\,\Lambda$"); ax.set_title(r"Euclidean $S=-I$ ($\chi=2$ vs $4$)")
    fig.tight_layout(); fig.savefig(OUT, dpi=130)
    return fig, axs


if __name__ == "__main__":
    main()
    print("figure ->", OUT)
