#!/usr/bin/env python3
"""
sigma_phase_gradient.py -- is the recipe sigma_mu a phase gradient sigma_mu = d_mu S/(m c), and which mass sets the action scale?
===============================================================================================================================

QGD defines sigma_mu = p_mu/(mc) = d_mu S/(mc) from the WKB phase of a mass-m Dirac field (complete_paper.tex, key-equations.md).
Two consequences can be tested against the recipes that are used everywhere else:

  A. A gradient is closed: d(sigma) = 0.  The Schwarzschild recipe in ingoing Kerr-Schild form is sigma = sqrt(F(r)) dv, so
     d(sigma) = (sqrt F)' dr ^ dv  -- nonzero for every non-constant F.  The same holds for the Kerr recipe.  The exact
     Hamilton-Jacobi gradient of a particle falling from rest, u = -dv + dr/(1 + sqrt(r_s/r)), is closed and has norm -1 (the recipe is null).
  B. In sigma = dS/(mc) the quantum scale is (mc/hbar) * (geometric integral of sigma): it contains the *test* mass m.  The radial
     part of the recipe gives  x_sigma(r) = (mc/hbar) int_0^r sqrt(r_s/r') dr' = 2 mc sqrt(r_s r)/hbar = (m/M)(S_BH/pi) sqrt(r/r_s),
     so the earlier Sigma ~ M c r (x = S_BH/2pi at the horizon) is reproduced only up to a factor 2 and only for m = M.
"""
import numpy as np
import sympy as sp

v, r = sp.symbols("v r", positive=True)
rs, G, M, m, hbar, c = sp.symbols("r_s G M m hbar c", positive=True)


def schwarzschild_ks_metric():
    f = rs / r
    return sp.Matrix([[-(1 - f), 1], [1, 0]]), sp.Matrix([[0, 1], [1, 1 - f]])      # g_ab and inverse in (v, r)


def curl_2d(sig):
    """d sigma for a 1-form sig = (sig_v, sig_r) depending on (v, r): the single component omega_{rv}."""
    return sp.diff(sig[0], r) - sp.diff(sig[1], v)


def recipe_one_form(F):
    return (sp.sqrt(F), 0)


def hamilton_jacobi_one_form():
    f = rs / r
    return (-1, 1 / (1 + sp.sqrt(f)))


def norm(sig):
    _, gi = schwarzschild_ks_metric()
    s = sp.Matrix(sig)
    return sp.simplify((s.T * gi * s)[0])


def kerr_curl_at(Mv=1.0, av=0.6, point=(0.7, -0.4, 1.1)):
    """Largest |d sigma| component of the Kerr recipe sigma = sqrt(H) k (Cartesian Kerr-Schild, ingoing) at one point."""
    X = sp.symbols("t x y z", real=True)
    _, x, y, z = X
    a = av
    R2 = x ** 2 + y ** 2 + z ** 2
    r2 = ((R2 - a ** 2) + sp.sqrt((R2 - a ** 2) ** 2 + 4 * a ** 2 * z ** 2)) / 2
    rr = sp.sqrt(r2)
    H = 2 * Mv * rr ** 3 / (rr ** 4 + a ** 2 * z ** 2)
    k = [1, -(rr * x + a * y) / (r2 + a ** 2), -(rr * y - a * x) / (r2 + a ** 2), -z / rr]
    sig = [sp.sqrt(H) * ki for ki in k]
    sub = dict(zip(X, (0.0,) + tuple(point)))
    W = np.array([[float((sp.diff(sig[n], X[mu]) - sp.diff(sig[mu], X[n])).subs(sub)) for n in range(4)] for mu in range(4)])
    return float(np.abs(W).max())


def phase_scale_ratio():
    """x_sigma = 2 m c sqrt(r_s r)/hbar  versus (m/M)(S_BH/pi) sqrt(r/r_s); returns the simplified difference."""
    rsv = 2 * G * M / c ** 2
    S_BH = 4 * sp.pi * G * M ** 2 / (c * hbar)
    xs = (m * c / hbar) * sp.integrate(sp.sqrt(rsv / r), (r, 0, r))
    return sp.simplify(xs - (m / M) * (S_BH / sp.pi) * sp.sqrt(r / rsv))


if __name__ == "__main__":
    F = rs / r
    print("d sigma (recipe)       :", sp.simplify(curl_2d(recipe_one_form(F))))
    print("d u     (HJ gradient)  :", sp.simplify(curl_2d(hamilton_jacobi_one_form())))
    print("|sigma|^2 recipe, HJ   :", norm(recipe_one_form(F)), ",", norm(hamilton_jacobi_one_form()))
    print("Kerr |d sigma|_max     :", kerr_curl_at())
    print("x_sigma identity resid.:", phase_scale_ratio())
