"""The literal master-metric quantum term on Kerr: angular-only tensor, horizon displacement, entropy shift."""
import os
import sys

import numpy as np
import pytest
import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "validation"))
import kerr_quantum_term as kq  # noqa: E402


def test_quantum_tensor_pulled_back_to_spheroidal_coordinates_is_angular_only():
    v, r, th, ph = kq.v, kq.r, kq.th, kq.ph
    a, M = kq.a, kq.M
    Y = [v, r, th, ph]
    xx = sp.sin(th) * (r * sp.cos(ph) - a * sp.sin(ph))
    yy = sp.sin(th) * (r * sp.sin(ph) + a * sp.cos(ph))
    zz = r * sp.cos(th)
    X = [v - r, xx, yy, zz]
    J = sp.Matrix(4, 4, lambda A, B: sp.diff(X[A], Y[B]))
    flat = sp.simplify(J.T * sp.diag(-1, 1, 1, 1) * J)
    expect = kq.kerr_metric().subs(M, 0)                       # M = 0 Kerr is flat space in spheroidal coordinates
    assert sp.simplify(flat - expect) == sp.zeros(4, 4)
    k = sp.Matrix([-(r * xx + a * yy) / (r ** 2 + a ** 2), -(r * yy - a * xx) / (r ** 2 + a ** 2), -zz / r])
    C = sp.zeros(4, 4)
    C[1:, 1:] = (sp.eye(3) - k * k.T) * (kq.H / kq.Sigma)
    assert sp.simplify(J.T * C * J - kq.delta_metric()) == sp.zeros(4, 4)


def test_horizon_displacement_closed_form_matches_inverse_metric():
    assert sp.simplify(kq.delta_rh_formula() - kq.delta_rh_from_inverse()) == 0


def test_linearised_horizon_agrees_with_the_exact_root_to_second_order():
    Mv, av, thv = 1.0, 0.6, 1.0
    rp = kq.r_plus(Mv, av)
    f = sp.lambdify((kq.M, kq.a, kq.b, kq.s, kq.th), kq.delta_rh_formula(), "numpy")
    errs = []
    for bv in (4e-3, 2e-3):
        exact = kq.exact_horizon_radius(Mv, av, bv, thv) - rp
        errs.append(abs(exact - float(f(Mv, av, bv, 1.0, thv))))
    assert errs[0] / errs[1] == pytest.approx(4.0, rel=0.1)    # residual is O(b^2)
    assert errs[0] < 1e-4


def test_nonrotating_limit_gives_pi_s_b_and_no_displacement():
    assert kq.entropy_shift_over_pi_b(1.0, 0.0) == pytest.approx(1.0, abs=1e-12)
    assert sp.simplify(kq.delta_rh_formula().subs(kq.a, 0)) == 0


@pytest.mark.parametrize("x", [0.1, 0.3, 0.6, 0.9])
def test_entropy_shift_is_pi_s_b_over_sqrt_one_minus_chi_squared(x):
    assert kq.entropy_shift_over_pi_b(1.0, x) == pytest.approx(1 / np.sqrt(1 - x ** 2), rel=1e-9)


def test_entropy_shift_scales_with_M_only_through_a_over_M_and_flips_with_sign():
    assert kq.entropy_shift_over_pi_b(3.0, 1.8) == pytest.approx(kq.entropy_shift_over_pi_b(1.0, 0.6), rel=1e-9)
    assert kq.area_shift(1.0, 0.6, sv=-1.0) == pytest.approx(-kq.area_shift(1.0, 0.6, sv=1.0), rel=1e-12)


def test_perturbed_horizon_is_not_a_uniform_rotation_at_first_order():
    thv = np.linspace(0.3, np.pi - 0.3, 5)
    assert np.ptp(kq.angular_velocity_profile(1.0, 0.6, 1e-9, thv)) < 1e-10
    sp1 = np.ptp(kq.angular_velocity_profile(1.0, 0.6, 1e-3, thv))
    sp2 = np.ptp(kq.angular_velocity_profile(1.0, 0.6, 2e-3, thv))
    assert sp1 > 1e-6 and sp2 / sp1 == pytest.approx(2.0, rel=0.1)
