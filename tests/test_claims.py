"""Independent checks of statements the QGD documents make.

Nothing here imports the QGD metric code: the algebra is redone in sympy so
that a pass means the claim is true, not that the code agrees with itself.
Run from the repository root:  python -m pytest -q tests/test_claims.py
"""
import math
import os
import re
import sys

import numpy as np
import pytest
import sympy as sp

ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(ROOT, "validation"))


def kappa(j):
    return math.sqrt(math.factorial(2 * j - 1) / 2 ** (2 * j - 2))


# ---------------------------------------------------------------- kappa_j
def test_kappa_values_in_engine_constants():
    import dark_matter as DM
    C = DM.QGDConstants
    for j, v in enumerate([C.KAPPA_1, C.KAPPA_2, C.KAPPA_3, C.KAPPA_4, C.KAPPA_5], start=1):
        assert v == pytest.approx(kappa(j), rel=1e-3), j


def test_readme_kappa_list_matches_formula():
    text = open(os.path.join(ROOT, "README.md"), encoding="utf-8").read()
    m = re.search(r"\\kappa = \[([^\]]*)\\ldots\]", text)
    assert m, "kappa list not found in README"
    vals = [float(x) for x in re.findall(r"[0-9.]+", m.group(1).replace("\\", " "))]
    assert len(vals) >= 6
    for j, v in enumerate(vals, start=1):
        assert v == pytest.approx(kappa(j), rel=5e-3), "kappa_%d: README %s vs formula %.4g" % (j, v, kappa(j))


def test_kappa_recurrence():
    # kappa_{j+1}^2 / kappa_j^2 = (2j)(2j+1) / 4
    for j in range(1, 8):
        assert kappa(j + 1) ** 2 / kappa(j) ** 2 == pytest.approx(2 * j * (2 * j + 1) / 4)


# ---------------------------------------------------- gravitational alpha_G
HBAR, C_, G_ = 1.054571817e-34, 2.99792458e8, 6.67430e-11


def alpha_G(M, m):
    return np.exp(1j * np.pi / 4) * math.sqrt(C_ * HBAR / (2 * G_ * M * m))


def test_alpha_identities():
    a = alpha_G(5.972e24, 9.109e-31)
    mod2 = C_ * HBAR / (2 * G_ * 5.972e24 * 9.109e-31)
    assert abs(a) ** 2 == pytest.approx(mod2)
    assert a ** 2 == pytest.approx(1j * mod2)
    assert a.real == pytest.approx(math.sqrt(mod2 / 2))
    assert mod2 == pytest.approx(4.36e-11, rel=0.01)    # value given in docs/corrections.md


def test_electron_electron_value():
    # with |alpha_G|^2 = c hbar / (2 G M m) the electron-electron value is ~3e44, not 1e-45
    # (1e-45 is G m_e^2 / (hbar c), the usual gravitational fine-structure constant)
    me = 9.1093837e-31
    assert C_ * HBAR / (2 * G_ * me * me) == pytest.approx(2.85e44, rel=0.01)


def test_docs_do_not_quote_old_earth_electron_value():
    stale = re.compile(r"1\.46\s*(×|\\times)\s*10(\^\{?39\}?|³⁹)")
    for name in ("key-equations.md", "GR-vs-QGD.md", "docs/complete_paper.tex", "docs/Ch1-Foundations.tex",
                 "docs/corrections.md"):
        text = open(os.path.join(ROOT, name), encoding="utf-8").read()
        if name == "docs/corrections.md":
            text = text.replace("`≈ 1.46×10³⁹`", "")        # the correction quotes the old value on purpose
        assert not stale.search(text), name


def test_tex_earth_electron_expression_matches_its_stated_value():
    expr = (3e8 * 1.055e-34) / (2 * 6.67e-11 * 5.97e24 * 9.11e-31)
    for name in ("docs/complete_paper.tex", "docs/Ch1-Foundations.tex"):
        text = open(os.path.join(ROOT, name), encoding="utf-8").read()
        m = re.search(r"\\approx\s*([0-9.]+)\\times10\^\{(-?[0-9]+)\}\.\s*\\end\{align\}", text)
        assert m, name
        assert float(m.group(1)) * 10 ** int(m.group(2)) == pytest.approx(expr, rel=0.01), name


def test_planck_mass_gives_one_half():
    mpl = math.sqrt(HBAR * C_ / G_)
    assert abs(alpha_G(mpl, mpl)) ** 2 == pytest.approx(0.5)


# ------------------------------------- Kerr-Schild form of the "recipes"
t, r, th, ph = sp.symbols("t r theta phi", positive=True)
m, q = sp.symbols("m q", positive=True)
X = (t, r, th, ph)


def ricci_and_einstein(g):
    ginv = g.inv()
    n = 4
    Gam = [[[sum(ginv[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                 for d in range(n)) / 2 for c in range(n)] for b in range(n)] for a in range(n)]

    def R(b, c):
        s = 0
        for a in range(n):
            s += sp.diff(Gam[a][b][c], X[a]) - sp.diff(Gam[a][b][a], X[c])
            for d in range(n):
                s += Gam[a][a][d] * Gam[d][b][c] - Gam[a][c][d] * Gam[d][b][a]
        return sp.simplify(s)
    Ric = sp.Matrix(n, n, lambda b, c: R(b, c))
    Rs = sp.simplify(sum(ginv[a, b] * Ric[a, b] for a in range(n) for b in range(n)))
    return Ric, Rs, sp.simplify(ginv * (Ric - Rs * g / 2))


def kerr_schild(h):
    """g = eta + h(r) k k with k = (1, 1, 0, 0) in (t, r, theta, phi) and eta in spherical form.
    This is eta_{mu nu} + sum_a eps_a sigma_mu^(a) sigma_nu^(a) with sigma^(a) = sqrt(|h_a|) k."""
    eta = sp.diag(-1, 1, r ** 2, r ** 2 * sp.sin(th) ** 2)
    k = sp.Matrix([1, 1, 0, 0])
    return eta + h * k * k.T


def test_schwarzschild_recipe_is_vacuum():
    # sigma_t = sqrt(2 m / r)  ->  h = 2 m / r
    Ric, Rs, _ = ricci_and_einstein(kerr_schild(2 * m / r))
    assert Ric == sp.zeros(4, 4) and Rs == 0


def test_reissner_nordstrom_recipe_has_em_stress():
    # eps = +1 mass term, eps = -1 charge term:  h = 2 m / r - q^2 / r^2
    _, Rs, Gmix = ricci_and_einstein(kerr_schild(2 * m / r - q ** 2 / r ** 2))
    assert Rs == 0                                      # traceless: Maxwell stress
    e = q ** 2 / r ** 4
    assert sp.simplify(Gmix[0, 0] + e) == 0 and sp.simplify(Gmix[1, 1] + e) == 0
    assert sp.simplify(Gmix[2, 2] - e) == 0 and sp.simplify(Gmix[3, 3] - e) == 0


def test_wrong_sign_of_charge_is_not_maxwell():
    # eps = +1 for the charge term gives the opposite-sign stress, i.e. a negative energy density
    _, _, Gmix = ricci_and_einstein(kerr_schild(2 * m / r + q ** 2 / r ** 2))
    assert sp.simplify(Gmix[0, 0] - q ** 2 / r ** 4) == 0
