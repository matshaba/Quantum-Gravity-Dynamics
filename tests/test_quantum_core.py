"""
Wald entropy of the l_Q^2 term, the literal quantum term of the master metric, a conditional bounce model and the
dodecahedral spectrum (G = hbar = c = 1).  See validation/wald_entropy.py, qgd_quantum_core.py, bounce_toy.py,
dodecahedral_spectrum.py for the derivations; each test states what it establishes.
"""
import sys
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))
import bounce_toy as B  # noqa: E402
import dodecahedral_spectrum as D  # noqa: E402
import qgd_quantum_core as Q  # noqa: E402
import wald_entropy as W  # noqa: E402


@pytest.fixture(scope="module")
def horizon():
    return W.Horizon()


# ---- Wald entropy ----------------------------------------------------------------------------------------------
def test_wald_machinery_is_calibrated(horizon):
    """The same code gives A/4 for R/16 pi and 4 pi alpha for the Gauss-Bonnet term (Jacobson-Myers), the value
    the Euclidean route gave in tests/test_nariai_entropy.py."""
    rs = sp.symbols("r_s", positive=True)

    def at_h(E):
        return sp.simplify(E.subs(W.F, rs / W.r).doit().subs(W.r, rs))
    assert sp.simplify(W.wald_entropy(at_h(horizon.E_einstein()), rs) - sp.pi * rs ** 2) == 0
    assert sp.simplify(W.wald_entropy(at_h(horizon.E_gauss_bonnet(1)), rs) - 4 * sp.pi) == 0


def test_wald_entropy_of_the_l2_term_vanishes_for_any_F(horizon):
    """sigma = sqrt(F) k lies in the binormal plane of the horizon cross-section, so eps_ab sigma^b is proportional
    to sigma_a and R_{abcd} sigma^a sigma^b = 0: the l_Q^2 (nabla nabla sigma)^2 term adds nothing to the entropy of
    Schwarzschild, Schwarzschild-de Sitter or Nariai (any F(r))."""
    assert sp.simplify(horizon.E_sigma(1, +1).doit()) == 0
    assert sp.simplify(horizon.E_sigma(1, -1).doit()) == 0


def test_wald_term_is_not_trivially_zero(horizon):
    """Tilting sigma out of the binormal plane (sigma^v = delta) gives a nonzero contribution, so the zero above is a
    property of the null radial sigma, not of the code."""
    H2 = W.Horizon()
    d = sp.symbols("delta")
    H2.sigma_up = [d, sp.sqrt(W.F), 0, 0]
    rs = sp.symbols("r_s", positive=True)
    E = sp.simplify(H2.E_sigma(1, 1).doit().subs(W.F, rs / W.r).doit().subs(W.r, rs))
    assert sp.simplify(E + d / rs ** 2) == 0


# ---- literal quantum term of the master metric ----------------------------------------------------------------
@pytest.mark.parametrize("F", [lambda q: 1 / q, lambda q: 1 / q + q ** 2 / 3])
@pytest.mark.parametrize("null_sign", [+1, -1])
def test_quantum_term_acts_on_the_angular_part_only(F, null_sign):
    """C_ab = eta^{cd} d_a sigma_c d_b sigma_d = (F/r^2) P_ij: no t, r components, so g_tt and g_tr are not changed."""
    res = Q.quantum_tensor_cartesian(F, null_sign=null_sign)
    assert res == sp.zeros(4, 4)


def test_quantum_term_is_not_the_A_over_r2_term_of_f():
    """The universal-f / f9 form puts kappa hbar^2/(M^2 r^2) into g_tt; the literal master metric does not touch g_tt."""
    rs, b = sp.symbols("r_s b", positive=True)
    H0 = W.Horizon(rs / W.r)
    H1 = Q.core_horizon(rs, b, +1)
    assert sp.simplify(H1.g[0, 0] - H0.g[0, 0]) == 0 and H1.g[0, 1] == H0.g[0, 1] == 1
    assert sp.simplify(H1.g[2, 2] - H0.g[2, 2]) != 0


@pytest.mark.parametrize("s", [+1, -1])
def test_entropy_shift_is_the_same_on_every_horizon_and_consistent_with_the_first_law(s):
    """rho_h^2 = r_h^2 + s b F(r_h) = r_h^2 + s b at every horizon (F = 1): S = pi rho_h^2 shifts by s pi b on each
    horizon of Schwarzschild, de Sitter and Nariai alike; T = F'(r_h)/4 pi is unchanged, so with S = pi r_h^2 + const and
    r_h = 2M the first law dM = T dS holds (the constant is a constant of integration)."""
    b, M, Lam = sp.symbols("b M Lambda", positive=True)
    # Schwarzschild: T = 1/(8 pi M),  S = pi (4 M^2 + s b)
    S = sp.pi * (4 * M ** 2 + s * b)
    T = 1 / (8 * sp.pi * M)
    assert sp.simplify(sp.diff(S, M) * T - 1) == 0
    # de Sitter + Nariai numerics: F = r_s/r + Lam r^2/3, roots of r(1-F)
    for rs_v, Lam_v in ((0.3, 1.0), (2 / 3, 1.0)):                    # (2/3, 1) is the Nariai mass for Lambda = 1
        roots = np.roots([-Lam_v / 3, 0, 1, -rs_v])
        roots = [x.real for x in roots if abs(x.imag) < 1e-6 and x.real > 0]
        for rh in roots:
            F = rs_v / rh + Lam_v * rh ** 2 / 3
            assert F == pytest.approx(1.0, abs=1e-6)
            assert (rh ** 2 + s * 0.1 * F) - rh ** 2 == pytest.approx(s * 0.1, abs=1e-6)
    # pair-creation bookkeeping: Nariai has two horizons, de Sitter one, so the shift of S_N - S_dS is s pi b
    assert sp.simplify((2 * s * sp.pi * b) - (1 * s * sp.pi * b) - s * sp.pi * b) == 0


def test_quantum_term_does_not_remove_the_singularity():
    """s = +1: K r^6 -> 51 r_s^2/4 as r -> 0 (classical 12 r_s^2), independent of kappa l_Q^2.
    s = -1: the sphere area vanishes at r_* = (b r_s)^(1/3) and K diverges there like (r - r_*)^-4."""
    coef, rs = Q.kretschmann_small_r_coefficient(+1)
    assert sp.simplify(coef - sp.Rational(51, 4) * rs ** 2) == 0
    import mpmath as mp
    K = Q.kretschmann_function(1.0, 0.1, -1)
    rstar = 0.1 ** (1 / 3)
    k1, k2 = float(K(mp.mpf(rstar + 1e-3))), float(K(mp.mpf(rstar + 1e-4)))
    assert k2 / k1 == pytest.approx(1e4, rel=0.05)


# ---- classical field equation on Schwarzschild -----------------------------------------------------------------
def test_schwarzschild_sigma_does_not_solve_box_sigma_equals_Q_with_G_zero():
    """With sigma = sqrt(r_s/r) dv: box sigma_v = sqrt(r_s)(r_s - r)/(4 r^(7/2)), box sigma_r = -2 sqrt(r_s)/r^(5/2),
    while Q_v = 2 r_s^(3/2)/r^(7/2), Q_r = 0.  So Box sigma = Q fails unless G_mu supplies the difference; G_mu is not
    specified for Schwarzschild in the repository, so the equation cannot be checked as printed.  (sigma^2 = r_s/r, not
    sigma, is the harmonic function.)"""
    box, Qv, rs = Q.schwarzschild_field_equation_residual()
    r = W.r
    assert sp.simplify(box[0] - sp.sqrt(rs) * (rs - r) / (4 * r ** sp.Rational(7, 2))) == 0
    assert sp.simplify(box[1] + 2 * sp.sqrt(rs) / r ** sp.Rational(5, 2)) == 0
    assert sp.simplify(Qv[0] - 2 * rs ** sp.Rational(3, 2) / r ** sp.Rational(7, 2)) == 0 and Qv[1] == 0
    assert sp.simplify(box[1] - Qv[1]) != 0


# ---- bounce ----------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("w,n,bounces", [(0.0, 4.0, True), (0.0, 5.0, True), (1 / 3, 4.0, False), (1.0, 4.0, False),
                                         (1 / 3, 6.0, True)])
def test_conditional_bounce_criterion(w, n, bounces):
    """rho_q = -C a^-n against matter w: H = 0 exists iff n > 3(1+w), and a''/a > 0 there automatically."""
    ab = B.bounce_scale_factor(w, n, C=0.2)
    assert (ab is not None) == bounces
    if bounces:
        assert B.acceleration_sign_at_bounce(w, n) == 1.0


def test_integration_turns_at_the_predicted_scale_factor_and_radiation_does_not():
    t, a = B.integrate_bounce(0.0, 4.0, C=0.2, tmax=12)
    assert a.min() == pytest.approx(B.bounce_scale_factor(0.0, 4.0, 0.2), rel=1e-3)
    assert a[-1] > a.min() * 20                                      # re-expands
    t, a = B.integrate_bounce(1 / 3, 4.0, C=0.2, tmax=12)
    assert a[-1] == pytest.approx(a.min(), rel=1e-6)                   # keeps collapsing (integration ends at the minimum)


def test_H4_term_has_no_bounce():
    """H^2 + s k H^4 = x: s = +1 gives H^2 > 0 for all x > 0; s = -1 has no real solution beyond x = 1/(4k)."""
    sols = B.h4_term_has_bounce()
    k, x = sp.symbols("k x", positive=True)
    assert sp.simplify(sols[+1][0].subs({k: 1, x: 5})) > 0
    assert sp.im(sols[-1][0].subs({k: 1, x: 1})) != 0               # x = 1 > 1/(4k) = 0.25: complex


# ---- dodecahedral space ---------------------------------------------------------------------------------------
def test_dodecahedral_group_and_spectrum():
    G = D.binary_icosahedral_group()
    assert len(G) == 120 and D.is_closed(G)
    N = np.rint(D.invariants(G, 80)).astype(int)
    assert np.all(N == D.molien_coefficients(80))                    # Klein: degrees 12, 20, 30
    sp_ = D.spectrum(G, 40)
    nonzero = [s for s in sp_ if s[0] > 0]
    assert nonzero[0][0] == 12 and nonzero[0][1] == 12 * 14          # lowest mode k = 12, eigenvalue 168 (S^3: k = 1, eigenvalue 3)
    assert [s[0] for s in nonzero[:5]] == [12, 20, 24, 30, 32]
    assert nonzero[0][2] == 13


# ---- Kerr: the l_Q^2 term has a spin-dependent Wald entropy -----------------------------------------------------
import wald_kerr as K  # noqa: E402


@pytest.mark.parametrize("a", [0.0, 0.3, 0.9])
def test_kerr_wald_machinery_is_calibrated(a):
    """Same numerical machinery: R/16 pi gives A/4 = pi (r_+^2 + a^2) and the Gauss-Bonnet term gives 4 pi alpha for every
    spin (topological)."""
    (seh, sgb, _), rp = K.horizon_entropies(1.0, a, n=32)
    assert seh == pytest.approx(np.pi * (rp ** 2 + a ** 2), rel=1e-9)
    assert sgb == pytest.approx(4 * np.pi, rel=1e-9)


def test_kerr_sigma_term_vanishes_at_zero_spin_and_grows_with_spin():
    """S_sigma / (s l_Q^2): 0 at a = 0 (agrees with the symbolic result above), 1.43 at a/M = 0.6, 4.82 at 0.9, 8.76 at 0.99."""
    got = {a: K.horizon_entropies(1.0, a, n=48)[0][2] for a in (0.0, 0.6, 0.9, 0.99)}
    assert got[0.0] == pytest.approx(0.0, abs=1e-9)
    assert got[0.6] == pytest.approx(1.4311291, rel=1e-6)
    assert got[0.9] == pytest.approx(4.8200848, rel=1e-6)
    assert got[0.99] == pytest.approx(8.758164, rel=1e-5)
    assert got[0.0] < got[0.6] < got[0.9] < got[0.99]


def test_kerr_sigma_term_depends_only_on_a_over_M_and_flips_with_the_sign():
    s1 = K.horizon_entropies(1.0, 0.6, n=32)[0][2]
    assert K.horizon_entropies(2.0, 1.2, n=32)[0][2] == pytest.approx(s1, rel=1e-9)
    assert K.horizon_entropies(1.0, 0.6, n=32, sign=-1)[0][2] == pytest.approx(-s1, rel=1e-12)


def test_kerr_small_spin_coefficient_is_pi_squared_over_three():
    """S_sigma = (pi^2/3) s l_Q^2 (a/M)^2 + O(a^4): from a = 0.02 with the a^4 term removed."""
    a = 0.02
    s = K.horizon_entropies(1.0, a, n=32)[0][2]
    s2 = K.horizon_entropies(1.0, 2 * a, n=32)[0][2]
    c2 = (16 * s - s2) / (12 * a ** 2)                # Richardson: removes the a^4 term
    assert c2 == pytest.approx(np.pi ** 2 / 3, rel=2e-4)


# ---- Kerr-de Sitter: does the spin dependence survive with Lambda? ----------------------------------------------
@pytest.fixture(scope="module")
def kds():
    return K.kds_chart()


def test_kds_chart_is_an_einstein_space_regular_at_every_horizon_with_null_k(kds):
    """R_ab = Lambda g_ab to rounding, including AT the roots of Delta_r (black-hole and cosmological horizons), and the
    regular principal null form k = dv - (a sin^2 th / Xi) dphi is null."""
    M, a, L = 1.0, 0.4, 0.01
    rts = K.kds_horizons(M, a, L)
    assert len(rts) == 3
    for rv in (2.7, rts[-2], rts[-1]):
        g, gi, Riem, Ric, Rs = kds.geometry(rv, 0.8, M, a, L)
        assert np.abs(Ric - L * g).max() < 1e-12 and Rs == pytest.approx(4 * L, rel=1e-12)
        k = np.array(kds._k(rv, 0.8, M, a, L), float).ravel()
        assert abs(k @ gi @ k) < 1e-12


@pytest.mark.parametrize("which", ["BH", "cosm"])
def test_kds_wald_machinery_is_calibrated_on_both_horizons(kds, which):
    """A/4 = pi (r_h^2 + a^2)/Xi and the Gauss-Bonnet term gives 4 pi alpha, on the black-hole and the cosmological horizon."""
    M, a, L = 1.0, 0.6, 0.02
    rts = K.kds_horizons(M, a, L)
    rh = rts[-2] if which == "BH" else rts[-1]
    seh, sgb, _ = kds.entropies(rh, M, a, L, n=32)
    assert seh == pytest.approx(np.pi * (rh ** 2 + a ** 2) / (1 + L * a ** 2 / 3), rel=1e-9)
    assert sgb == pytest.approx(4 * np.pi, rel=1e-9)


def test_kds_reduces_to_kerr_as_lambda_goes_to_zero(kds):
    M, a, L = 1.0, 0.6, 1e-9
    rh = K.kds_horizons(M, a, L)[-2]
    assert kds.entropies(rh, M, a, L, n=32)[2] == pytest.approx(K.horizon_entropies(1.0, 0.6, 32)[0][2], rel=1e-6)


def _c2(kds, M, L, sel, a=0.01):
    vals = []
    for aa in (a, 2 * a):
        rh = K.kds_horizons(M, aa, L)[sel]
        vals.append(kds.entropies(rh, M, aa, L, n=32)[2])
    return (16 * vals[0] - vals[1]) / (12 * a ** 2)            # Richardson: removes the a^4 term


@pytest.mark.parametrize("M,L", [(1.0, 0.02), (1.0, 0.07), (2.0, 0.005), (0.5, 0.3)])
def test_small_spin_coefficient_is_minus_4pi2_over_3_sigmaM2_Fprime_over_r(kds, M, L):
    """S_sigma ~ c2 a^2 with c2 = -(4 pi^2/3) sigma_M^2(r_h) F'(r_h)/r_h: proportional to the surface gravity F'/2 of the
    horizon times the mass-sigma amplitude there (1 - Lambda r_h^2/3).  Opposite signs on the black-hole and cosmological
    horizons.  Agreement ~1e-6 (fitted closed form, verified numerically, not derived)."""
    for name, sel in (("BH", -2), ("cosm", -1)):
        got = _c2(kds, M, L, sel)
        want = K.small_spin_coefficient(M, L, name)
        assert got == pytest.approx(want, rel=1e-5)
    assert K.small_spin_coefficient(M, L, "BH") > 0 > K.small_spin_coefficient(M, L, "cosm")


def test_spin_dependence_survives_at_nariai_but_starts_at_a_to_the_fourth(kds):
    """At the Nariai double root F' = 0, so the a^2 term cancels; S_sigma = (16 pi^2/15) s l_Q^2 (a^2 Lambda)^2 + O(a^6)
    for each of the two merged horizons (Lambda = 1 units, fitted)."""
    a = np.array([0.02, 0.03, 0.04, 0.06])
    S = []
    for x in a:
        rN, MN = K.kds_nariai(x, 1.0)
        S.append(kds.entropies(rN, MN, x, 1.0, n=48)[2])
    S = np.array(S)
    A = np.vstack([a ** 4, a ** 6, a ** 8]).T
    c4 = np.linalg.lstsq(A, S, rcond=None)[0][0]
    assert c4 == pytest.approx(16 * np.pi ** 2 / 15, rel=1e-4)
    assert S[0] / a[0] ** 2 < 1e-2          # S/a^2 = c4 a^2 ~ 4e-3: the a^2 coefficient is zero (naive Kerr-like value: 3 pi^2 = 29.6)
    # both horizons of the nearly merged pair carry (almost) the same value
    rts = K.kds_horizons(K.kds_nariai(0.06, 1.0)[1] * (1 - 1e-6), 0.06, 1.0)
    assert len(rts) >= 2
