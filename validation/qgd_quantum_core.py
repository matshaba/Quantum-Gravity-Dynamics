#!/usr/bin/env python3
"""
qgd_quantum_core.py -- what the quantum term of the QGD master metric actually does (G = hbar = c = 1).
=======================================================================================================

The printed master metric is  g = eta -+ sum eps sigma sigma -+ kappa l_Q^2 d_a sigma^c d_b sigma_c,  built in Cartesian
(Kerr-Schild) coordinates and then transformed.  For a radial null sigma_mu = sqrt(F(r)) k_mu  (k = dt + dr) the last term
can be computed exactly, and it is NOT the  kappa hbar^2/(M^2 c^2 r^2)  term that the universal f and
`core/master_metric.py::f9` add to g_tt:

    C_{ab} = eta^{cd} d_a sigma_c d_b sigma_d  =  (F / r^2) P_ij   (spatial transverse projector, no t or r parts)

so only the angular part of the metric changes,

    ds^2 = -(1 - F) dv^2 + 2 dv dr + rho^2 dOmega^2 ,     rho^2 = r^2 + s kappa l_Q^2 F(r) ,   s = +1 or -1,

with s fixed by the signature convention (the formula is written with eta = diag(+,-,-,-) in the text and g_tt = -(1-f)
in the code; the two readings give s = -1 and s = +1).  Consequences, each tested in tests/test_quantum_core.py:

  * g_tt, g_tr unchanged  ->  horizons F = 1 and surface gravity unchanged  ->  T = F'(r_h)/4 pi unchanged.
  * horizon area 4 pi rho_h^2 = 4 pi (r_h^2 + s kappa l_Q^2)  ->  S = A/4 shifts by s pi kappa l_Q^2 on EVERY horizon
    (black hole, cosmological, Nariai pair): mass independent, same on each horizon, consistent with the first law
    (a constant of integration).  The Wald entropy of the l_Q^2 (nabla nabla sigma)^2 term of the action is zero
    (validation/wald_entropy.py), so this is the whole correction in this reading.
  * singularity: for Schwarzschild  K r^6 -> 51 r_s^2/4  as r -> 0 (s = +1; classical value 12 r_s^2), so it is not
    removed; for s = -1 the sphere area vanishes at r_* = (kappa l_Q^2 r_s)^(1/3) and K diverges there as (r - r_*)^-4.
"""
import sympy as sp

import wald_entropy as W

r = W.r


# ---------------------------------------------------------------- the tensor, from the Cartesian definition
def quantum_tensor_cartesian(F_of_r, point=(3, 4, 12), null_sign=+1):
    """C_ab = eta^{cd} d_a sigma_c d_b sigma_d for sigma_mu = sqrt(F) (1, null_sign x_i / r), at a Cartesian point."""
    t, x, y, z = sp.symbols("t x y z", real=True)
    X = [t, x, y, z]
    rr = sp.sqrt(x ** 2 + y ** 2 + z ** 2)
    s = sp.sqrt(F_of_r(rr))
    sig = [s, null_sign * s * x / rr, null_sign * s * y / rr, null_sign * s * z / rr]
    eta_inv = sp.diag(-1, 1, 1, 1)
    d = [[sp.diff(sig[g], X[a]) for g in range(4)] for a in range(4)]
    C = sp.Matrix(4, 4, lambda a, b: sum(eta_inv[g, h] * d[a][g] * d[b][h] for g in range(4) for h in range(4)))
    sub = dict(zip((x, y, z), point))
    Cn = sp.simplify(C.subs(sub))
    r0 = sp.sqrt(sum(p ** 2 for p in point))
    n = sp.Matrix(point) / r0
    P = sp.eye(3) - n * n.T
    pred = sp.Matrix(4, 4, lambda a, b: 0 if a == 0 or b == 0 else F_of_r(r0) / r0 ** 2 * P[a - 1, b - 1])
    return sp.simplify(Cn - pred)


# ---------------------------------------------------------------- the geometry
def core_horizon(rs, b, sign):
    """Horizon object for F = r_s/r, rho^2 = r^2 + sign * b * F   (b = kappa l_Q^2)."""
    return W.Horizon(rs / r, sp.sqrt(r ** 2 + sign * b * rs / r))


def kretschmann_small_r_coefficient(sign=+1):
    rs, b = sp.symbols("r_s b", positive=True)
    H = core_horizon(rs, b, sign)
    K = sp.simplify(H._full(H.Rl, H.Rl))
    return sp.simplify(sp.limit(sp.simplify(K * r ** 6), r, 0)), rs


def kretschmann_function(rs_val, b_val, sign):
    rs, b = sp.symbols("r_s b", positive=True)
    H = core_horizon(rs, b, sign)
    K = sp.simplify(H._full(H.Rl, H.Rl))
    return sp.lambdify(r, K.subs({rs: rs_val, b: b_val}), "mpmath")


# ---------------------------------------------------------------- classical field equation on Schwarzschild
def schwarzschild_field_equation_residual():
    """(box sigma_mu, Q_mu) on Schwarzschild with sigma = sqrt(r_s/r) dv, Q_mu from docs/complete_paper.tex:
    Q_mu = sigma_mu (nabla_a sigma_b nabla^a sigma^b) + (nabla_mu sigma_a)(sigma_b nabla^b sigma^a)."""
    rs = sp.symbols("r_s", positive=True)
    X, n = W.X, 4
    g = W.metric(rs / r)
    gi = g.inv()
    sig = [sp.sqrt(rs / r), 0, 0, 0]
    Gam = [[[sum(gi[a, d] * (sp.diff(g[d, bb], X[c]) + sp.diff(g[d, c], X[bb]) - sp.diff(g[bb, c], X[d]))
                 for d in range(n)) / 2 for c in range(n)] for bb in range(n)] for a in range(n)]
    D = sp.Matrix(n, n, lambda a, bb: sp.diff(sig[bb], X[a]) - sum(Gam[k][a][bb] * sig[k] for k in range(n)))
    DD = [[[sp.diff(D[a, bb], X[c]) - sum(Gam[k][c][a] * D[k, bb] for k in range(n))
            - sum(Gam[k][c][bb] * D[a, k] for k in range(n)) for bb in range(n)] for a in range(n)] for c in range(n)]
    box = [sp.simplify(sum(gi[c, a] * DD[c][a][bb] for c in range(n) for a in range(n))) for bb in range(n)]
    D2 = sum(gi[a, c] * gi[bb, d] * D[a, bb] * D[c, d] for a in range(n) for bb in range(n) for c in range(n)
             for d in range(n))
    sU = [sum(gi[a, bb] * sig[bb] for bb in range(n)) for a in range(n)]
    Q = [sp.simplify(sig[m] * D2 + sum(D[m, a] * sum(sU[bb] * gi[a, c] * D[bb, c] for bb in range(n) for c in range(n))
                                       for a in range(n))) for m in range(n)]
    return box, Q, rs


def main():
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pathlib import Path

    print("C_ab residual vs (F/r^2) P_ij:", quantum_tensor_cartesian(lambda q: 1 / q))
    coef, rs = kretschmann_small_r_coefficient(+1)
    print("K r^6 -> ", coef, "  (classical 12 r_s^2)")
    box, Q, rs = schwarzschild_field_equation_residual()
    print("box sigma_mu =", box, "\nQ_mu =", Q)

    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    rs_v, b_v = 1.0, 0.1
    rr = np.linspace(0.05, 2.0, 400)
    ax[0].plot(rr, rr ** 2, "k--", label=r"$r^2$ (classical)")
    ax[0].plot(rr, rr ** 2 + b_v * rs_v / rr, label=r"$s=+1$: $r^2+\kappa\ell^2 r_s/r$")
    ax[0].plot(rr, np.clip(rr ** 2 - b_v * rs_v / rr, 0, None), label=r"$s=-1$: $r^2-\kappa\ell^2 r_s/r$")
    ax[0].axvline(rs_v, color="gray", lw=0.6)
    ax[0].set_xlabel("r"); ax[0].set_ylabel(r"$\rho^2$ (areal radius squared)"); ax[0].set_ylim(0, 4)
    ax[0].set_title("angular part of the metric (horizon at r = 1)"); ax[0].legend(fontsize=7)
    import mpmath as mp
    for sg, lab in ((+1, "s=+1"), (-1, "s=-1")):
        Kf = kretschmann_function(rs_v, b_v, sg)
        xs = np.logspace(-3, 0, 60)
        ax[1].loglog(xs, [float(Kf(mp.mpf(x))) * x ** 6 for x in xs], label=lab)
    ax[1].axhline(12, color="k", ls="--", lw=0.8, label="Schwarzschild 12")
    ax[1].axhline(51 / 4, color="C0", ls=":", lw=0.8)
    ax[1].set_xlabel("r"); ax[1].set_ylabel(r"$K\,r^6$"); ax[1].set_title("curvature near r = 0"); ax[1].legend(fontsize=7)
    fig.tight_layout()
    out = Path(__file__).resolve().parent / "quantum_core.png"
    fig.savefig(out, dpi=130)
    print("figure ->", out)


if __name__ == "__main__":
    main()
