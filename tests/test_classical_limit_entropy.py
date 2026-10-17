"""QGD's hbar -> 0 against the high-entropy limit (Carcassi-Landini-Aidala, arXiv:2411.00972)."""
import os
import sys

import numpy as np
import pytest
import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "validation"))
import classical_limit_entropy as cl  # noqa: E402


def test_quantum_entropy_correction_is_minus_one_over_24_x_squared():
    ser = sp.expand(cl.entropy_difference_series())
    assert ser.coeff(cl.x, -2) == sp.Rational(-1, 24)
    assert ser.coeff(cl.x, -1) == 0 and ser.coeff(cl.x, 0) == 0


@pytest.mark.parametrize("xv", [20.0, 100.0])
def test_series_matches_the_exact_gaussian_entropies(xv):
    d = float(cl.S_Q(xv) - cl.S_C(xv))
    assert d == pytest.approx(-1 / (24 * xv ** 2) - 1 / (320 * xv ** 4), rel=1e-4)


def test_entropy_correction_is_second_order_in_hbar():
    # hbar -> hbar/2 at fixed Sigma doubles x; the correction falls by 4 (O(hbar^2)), the same order as every QGD term
    d1, d2 = cl.S_Q(50.0) - cl.S_C(50.0), cl.S_Q(100.0) - cl.S_C(100.0)
    assert d1 / d2 == pytest.approx(4.0, rel=1e-3)


def test_stretch_map_divides_the_correction_by_lambda_squared():
    for lamv in (2.0, 5.0, 10.0):
        d0 = cl.S_stretched(50.0, 1.0) - cl.S_C(50.0)
        dl = cl.S_stretched(50.0, lamv) - cl.S_C(50.0)
        assert d0 / dl == pytest.approx(lamv ** 2, rel=2e-2)


def test_horizon_action_scale_equals_BH_entropy_over_two_pi():
    ratio, sbh = cl.horizon_ratio()
    assert sp.simplify(ratio - sbh / (2 * sp.pi)) == 0
    assert sbh == 4 * sp.pi * cl.G * cl.M ** 2 / (cl.c * cl.hbar)


def test_expansion_parameter_at_horizon_is_two_pi_over_S_BH_squared():
    ratio, sbh = cl.horizon_ratio()
    eps = sp.simplify(1 / ratio ** 2)
    assert sp.simplify(eps - (2 * sp.pi / sbh) ** 2) == 0


def test_stretch_rescales_S_BH_by_lambda_and_eps_by_lambda_inverse_squared():
    ratio, sbh = cl.horizon_ratio()
    lam = cl.lam
    sub = {cl.hbar: cl.hbar / lam}
    assert sp.simplify(sbh.subs(sub) - lam * sbh) == 0
    assert sp.simplify((1 / ratio ** 2).subs(sub) - (1 / ratio ** 2) / lam ** 2) == 0


def test_entropy_zeros_mark_the_end_of_the_classical_description():
    xc, xq = cl.entropy_zero_crossings()
    assert cl.S_C(xc) == pytest.approx(0.0, abs=1e-12)
    assert float(cl.S_Q(xq)) == pytest.approx(0.0, abs=1e-9)
    assert cl.S_C(0.2) < 0                                   # "fictitious" negative-entropy states
