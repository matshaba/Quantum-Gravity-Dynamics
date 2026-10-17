"""Tests of the places where the QGD claims go beyond what the algebra supports.

Each test states a claim from the notes and either confirms it (passes) or records that it fails
(xfail with the reason).  Geometry is redone independently with sympy/numpy in units G = c = 1.
"""
import os
import sys

import numpy as np
import pytest
import sympy as sp

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "core"))
sys.path.insert(0, os.path.join(ROOT, "validation"))

# ----------------------------------------------------------------- numeric Ricci tensor
X = sp.symbols("t x y z", real=True)
ETA = sp.diag(-1, 1, 1, 1)


def ricci_at(g_sym, point):
    """Ricci tensor R_{mu nu} of a symbolic metric g(x) at a numeric point, from g, dg and ddg."""
    n = 4
    sub = dict(zip(X, point))
    g = np.array(g_sym.subs(sub).evalf(), dtype=float)
    dg = np.zeros((n, n, n))
    ddg = np.zeros((n, n, n, n))
    for l in range(n):
        d1 = g_sym.diff(X[l])
        dg[l] = np.array(d1.subs(sub).evalf(), dtype=float)
        for k in range(n):
            ddg[l, k] = np.array(d1.diff(X[k]).subs(sub).evalf(), dtype=float)
    gi = np.linalg.inv(g)
    # Gamma_{s m n} = 1/2 (d_m g_{s n} + d_n g_{s m} - d_s g_{m n})
    Gl = 0.5 * (np.einsum("msn->smn", dg) + np.einsum("nsm->smn", dg) - dg)
    Gam = np.einsum("rs,smn->rmn", gi, Gl)
    dgi = -np.einsum("ra,lab,bs->lrs", gi, dg, gi)
    dGl = 0.5 * (np.einsum("lmsn->lsmn", ddg) + np.einsum("lnsm->lsmn", ddg) - np.einsum("lsmn->lsmn", ddg))
    dGam = np.einsum("lrs,smn->lrmn", dgi, Gl) + np.einsum("rs,lsmn->lrmn", gi, dGl)
    # explicit loops for the remaining terms (clarity over speed)
    R = np.zeros((n, n))
    for m in range(n):
        for nn in range(n):
            s = 0.0
            for r in range(n):
                s += dGam[r, r, m, nn] - dGam[nn, r, m, r]
                for l in range(n):
                    s += Gam[r, r, l] * Gam[l, m, nn] - Gam[r, nn, l] * Gam[l, m, r]
            R[m, nn] = s
    return R, g


def ks_one_form_term(mass, centre, amp=None):
    """Kerr-Schild piece  h k_mu k_nu  for a static mass at `centre`, k = (1, n_hat)."""
    cx, cy, cz = centre
    dx, dy, dz = X[1] - cx, X[2] - cy, X[3] - cz
    rr = sp.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
    k = sp.Matrix([1, dx / rr, dy / rr, dz / rr])
    h = 2 * mass / rr if amp is None else amp(rr)
    return h, k


P = (0.0, 1.3, 0.7, 0.4)    # field point (t, x, y, z)


def test_ricci_routine_vanishes_for_single_schwarzschild():
    h, k = ks_one_form_term(1.0, (0, 0, 0))
    R, _ = ricci_at(ETA + h * k * k.T, P)
    assert np.max(np.abs(R)) < 1e-8


def test_two_body_superposition_is_not_vacuum():
    """The QGD 'exact N-body' sum of sigma-fields: g = eta + (sum_a sigma_a k_a)(sum_b sigma_b k_b)."""
    parts = [ks_one_form_term(1.0, (-3.0, 0, 0)), ks_one_form_term(1.0, (3.0, 0, 0))]
    s = sum((sp.sqrt(h) * k for h, k in parts), sp.zeros(4, 1))
    R, _ = ricci_at(ETA + s * s.T, P)
    assert np.max(np.abs(R)) > 1e-4          # vacuum would give 0


def test_linear_superposition_of_kerr_schild_terms_is_not_vacuum_either():
    parts = [ks_one_form_term(1.0, (-3.0, 0, 0)), ks_one_form_term(1.0, (3.0, 0, 0))]
    g = ETA + sum((h * k * k.T for h, k in parts), sp.zeros(4, 4))
    R, _ = ricci_at(g, P)
    assert np.max(np.abs(R)) > 1e-4


def test_two_body_residual_is_second_order_in_G():
    """The failure is the cross term, so it scales like (G M)^2: halving the masses quarters it."""
    def resid(m):
        parts = [ks_one_form_term(m, (-3.0, 0, 0)), ks_one_form_term(m, (3.0, 0, 0))]
        g = ETA + sum((h * k * k.T for h, k in parts), sp.zeros(4, 4))
        R, _ = ricci_at(g, P)
        return np.max(np.abs(R))
    assert resid(0.005) / resid(0.01) == pytest.approx(0.25, rel=0.05)


def test_common_null_direction_adds_exactly_in_mixed_ricci():
    """Where the recipe IS exact: terms sharing one geodesic null k (mass + charge on the same centre).
    The mixed Ricci tensor R^mu_nu of a Kerr-Schild metric is linear in h, so Q-terms add to mass terms."""
    h_m, k = ks_one_form_term(1.0, (0, 0, 0))
    rr = sp.sqrt(X[1] ** 2 + X[2] ** 2 + X[3] ** 2)
    h_q = -0.5 / rr ** 2

    def mixed(h):
        R, g = ricci_at(ETA + h * k * k.T, P)
        return np.linalg.inv(g) @ R
    lhs = mixed(h_m + h_q)
    rhs = mixed(h_m) + mixed(h_q)
    assert np.max(np.abs(lhs - rhs)) < 1e-8
    assert np.max(np.abs(mixed(h_q))) > 1e-4      # and the charge term does carry stress-energy


# ----------------------------------------------- static f-metrics: what the universal f implies
r_, G_, c_ = sp.symbols("r G c", positive=True)
f = sp.Function("f")


def einstein_tt(fexpr):
    """G^t_t of ds^2 = -(1-f) dt^2 + dr^2/(1-f) + r^2 dOmega^2 (G=c=1 form)."""
    F = 1 - fexpr
    return sp.simplify((r_ * sp.diff(F, r_) + F - 1) / r_ ** 2)


def implied_density(fexpr):
    """rho = -G^t_t / (8 pi) = (r f)' / (8 pi r^2)."""
    return sp.simplify(-einstein_tt(fexpr) / (8 * sp.pi))


def test_universal_f_is_additive_in_the_einstein_tensor():
    a, b = sp.symbols("a b", positive=True)
    f1, f2 = a / r_, b / r_ ** 2
    assert sp.simplify(einstein_tt(f1 + f2) - einstein_tt(f1) - einstein_tt(f2)) == 0


def test_mass_term_is_vacuum_charge_term_has_positive_density():
    m, q = sp.symbols("m q", positive=True)
    assert implied_density(2 * m / r_) == 0
    assert sp.simplify(implied_density(-q ** 2 / r_ ** 2) - q ** 2 / (8 * sp.pi * r_ ** 4)) == 0


def test_quantum_term_as_written_has_negative_energy_density():
    """README: f contains + kappa hbar^2/(M^2 c^2 r^2).  That is the Reissner-Nordstrom term with Q^2 < 0."""
    A = sp.symbols("A", positive=True)          # A = kappa (hbar/(M c))^2 > 0 in G = c = 1 units
    assert sp.simplify(implied_density(A / r_ ** 2) + A / (8 * sp.pi * r_ ** 4)) == 0


def kretschmann(fexpr):
    F = 1 - fexpr
    return sp.simplify(sp.diff(F, r_, 2) ** 2 + 4 * sp.diff(F, r_) ** 2 / r_ ** 2 + 4 * (1 - F) ** 2 / r_ ** 4)


def test_adding_the_quantum_term_does_not_remove_the_singularity():
    m, A = sp.symbols("m A", positive=True)
    K = kretschmann(2 * m / r_ + A / r_ ** 2)
    lim = sp.limit(K * r_ ** 8, r_, 0)
    assert lim != 0 and lim.is_finite             # K ~ r^-8 as r -> 0: still singular


def test_a_regular_core_needs_f_to_vanish_like_r_squared():
    """For contrast: a Hayward-type f = 2 m r^2/(r^3 + 2 m l^2) has finite curvature at r = 0."""
    m, l = sp.symbols("m l", positive=True)
    K = kretschmann(2 * m * r_ ** 2 / (r_ ** 3 + 2 * m * l ** 2))
    assert sp.limit(K, r_, 0).is_finite


# ------------------------------------------------------------------- CMB acoustic peaks
def test_cmb_claim_less_than_one_percent_is_not_reproduced():
    import check_cmb_peaks as CP
    res = CP.main([])
    rel = np.abs(res["qgd"] / res["observed"] - 1)
    assert rel[0] > 0.2                             # first peak: 279.6 vs ~221
    assert rel.max() > 0.01                         # the "<1%" statement fails
    line_rel = np.abs(res["line"] / res["observed"] - 1)
    assert line_rel.max() < 0.04                    # a two-parameter line does far better


def test_cmb_A_is_not_fitted_to_first_peak():
    import check_cmb_peaks as CP
    res = CP.main([])
    assert abs(res["a_qgd"] - res["observed"][0]) / res["observed"][0] > 0.2


# ------------------------------------------------- reconstructed MasterMetric / N-body
def test_master_metric_reproduces_gr_components():
    import numpy as np
    from graviton_field import SigmaField, G, c
    from master_metric import MasterMetric
    M, r, th = 1.98847e30, 1e10, np.pi / 2
    rs = 2 * G * M / c ** 2
    g = MasterMetric(SigmaField.schwarzschild(M)).at(r, th)
    assert g[0, 0] == pytest.approx(-(1 - rs / r), rel=1e-12)
    assert g[1, 1] == pytest.approx(1 / (1 - rs / r), rel=1e-12)
    a = 0.5 * G * M / c ** 2
    gk = MasterMetric(SigmaField.kerr(M, a)).at(r, th)
    Sig = r ** 2 + a ** 2 * np.cos(th) ** 2
    assert gk[0, 3] == pytest.approx(-2 * G * M * r * a * np.sin(th) ** 2 / (c ** 2 * Sig), rel=1e-9)
    assert MasterMetric(SigmaField.schwarzschild(M)).gravitational_energy_density(r) == pytest.approx(
        G * M / (4 * c ** 2 * r ** 3), rel=1e-9)


def test_binary_merger_prediction_equal_masses_is_four_rs():
    from two_and_three_body_solutions import QGDNBody, Body, G, c
    M = 30 * 1.98847e30
    sysm = QGDNBody([Body(M, [0, 0, 0]), Body(M, [1, 0, 0])])
    assert sysm.binary_merger_prediction() / (2 * G * M / c ** 2) == pytest.approx(4.0)


def test_merger_condition_with_total_sigma_gives_eight_rs():
    """The factor 2 is a choice: setting the TOTAL sigma^2 (self + cross) to 1 gives 8 r_s for equal masses."""
    from two_and_three_body_solutions import G, c
    M = 1.0
    d = 16 * G * M / c ** 2                       # sigma_total^2 = 16 G M / (c^2 d) = 1
    assert d / (2 * G * M / c ** 2) == pytest.approx(8.0)


def test_qgd_vs_gr_script_runs_and_only_fp_check_fails():
    import subprocess
    r = subprocess.run([sys.executable, "QGD_vs_GR.py"], cwd=os.path.join(ROOT, "core"), capture_output=True,
                       text=True, timeout=300, env=dict(os.environ, MPLBACKEND="Agg"))
    assert r.returncode == 0, r.stderr[-400:]
    assert "FINAL RESULTS: 32/33 tests passed" in r.stdout
    assert "f_peak accuracy:  55.1%  ✗" in r.stdout      # GR ISCO frequency (67 Hz) vs the ~150 Hz merger peak


# ------------------------------------------------ Schwarzschild-de Sitter -> Nariai, and the entropy constant
def test_sigma_lambda_recipe_reaches_the_nariai_limit():
    """f = r_s/r + H^2 r^2/c^2 (the sigma_Lambda recipe) has a double horizon at r = 1/sqrt(Lambda) with
    Lambda = 3 H^2/c^2, area 4 pi/Lambda and G M = c^3 / (3 sqrt3 H): the Nariai geometry dS_2 x S^2."""
    x, H, c, G = sp.symbols("x H c G", positive=True)
    rs = sp.symbols("r_s", positive=True)
    P = x - rs - H ** 2 * x ** 3 / c ** 2                       # x (1 - f) = 0
    sol = sp.solve([P, sp.diff(P, x)], [x, rs], dict=True)
    sol = [s for s in sol if s[x].is_positive and s[rs].is_positive]
    assert len(sol) == 1
    xN, rsN = sp.simplify(sol[0][x]), sp.simplify(sol[0][rs])
    Lam = 3 * H ** 2 / c ** 2
    assert sp.simplify(xN - 1 / sp.sqrt(Lam)) == 0
    assert sp.simplify(4 * sp.pi * xN ** 2 - 4 * sp.pi / Lam) == 0
    assert sp.simplify(rsN * c ** 2 / 2 - c ** 3 / (3 * sp.sqrt(3) * H)) == 0     # G M


def test_constant_entropy_correction_drops_out_of_thermodynamic_differences():
    """GR-vs-QGD.md: S = A/4G + const with const mass-independent.  Temperature 1/T = dS/dM and any entropy
    difference are blind to it, so it cannot be observed through them (a statement about observability,
    not about whether the constant is right)."""
    M, k = sp.symbols("M k", positive=True)
    S0 = 4 * sp.pi * M ** 2                                       # A/4 for G = c = hbar = 1, A = 16 pi M^2
    S1 = S0 + k
    assert sp.simplify(sp.diff(S1, M) - sp.diff(S0, M)) == 0
    assert sp.simplify((S1.subs(M, 2) - S1.subs(M, 1)) - (S0.subs(M, 2) - S0.subs(M, 1))) == 0
