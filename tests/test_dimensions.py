"""SI dimensional analysis of formulas written in the QGD notes.

A formula is a product of constants and variables; its dimension is the vector of exponents of
(mass, length, time).  A passing test means the formula as written is dimensionally what the notes say.
"""
from fractions import Fraction as F

import pytest

# (M, L, T)
DIM = {
    "G": (-1, 3, -2), "c": (0, 1, -1), "hbar": (1, 2, -1), "kB": (1, 2, -2),   # kB per kelvin dropped
    "mass": (1, 0, 0), "length": (0, 1, 0), "time": (0, 0, 1), "one": (0, 0, 0),
    "H": (0, 0, -1),
}


def d(**powers):
    out = [F(0)] * 3
    for name, p in powers.items():
        for i in range(3):
            out[i] += F(p) * DIM[name][i]
    return tuple(out)


LENGTH = d(length=1)
ENERGY = (1, 2, -2)
MOMENTUM = (1, 1, -1)
ACCEL = (0, 1, -2)
DIMLESS = (0, 0, 0)


def eq(a, b):
    return tuple(F(x) for x in a) == tuple(F(x) for x in b)


def test_alpha_G_squared_is_dimensionless():
    assert eq(d(hbar=1, c=1, G=-1, mass=-2), DIMLESS)


def test_field_energy_and_momentum():
    # E = m c^2 / (sqrt2 |alpha|^2) with 1/|alpha|^2 dimensionless
    assert eq(d(mass=1, c=2), ENERGY)
    assert eq(d(mass=1, c=1), MOMENTUM)
    # P = sqrt2 G M m^2 / hbar
    assert eq(d(G=1, mass=3, hbar=-1), MOMENTUM)


def test_sigma_fields_are_dimensionless():
    assert eq(d(G=1, mass=1, c=-2, length=-1), DIMLESS)            # sqrt(2GM/(c^2 r))
    assert eq(d(H=1, length=1, c=-1), DIMLESS)                      # H r / c


def test_quantum_term_in_f_is_dimensionless():
    assert eq(d(hbar=2, mass=-2, c=-2, length=-2), DIMLESS)         # kappa hbar^2/(M^2 c^2 r^2)


def test_quantum_correction_to_g_tt_is_dimensionless():
    assert eq(d(G=1, hbar=2, mass=-1, c=-4, length=-3), DIMLESS)    # G hbar^2/(M c^4 r^3)


def test_maximum_acceleration():
    assert eq(d(mass=1, c=3, hbar=-1), ACCEL)                        # 3 m c^3 / hbar


def test_planck_length_is_a_length():
    assert eq(d(G=F(1, 2), hbar=F(1, 2), c=F(-3, 2)), LENGTH)       # sqrt(G hbar / c^3)


@pytest.mark.xfail(reason="ell_Q = sqrt(G hbar^2 / c^4), as written in README.md and core/master_metric.py, has "
                          "dimension (mass * length^3)^(1/2), not length; the Planck length is sqrt(G hbar / c^3)",
                   strict=True)
def test_ell_Q_as_written_is_a_length():
    assert eq(d(G=F(1, 2), hbar=1, c=-2), LENGTH)


@pytest.mark.xfail(reason="rho_grav = G M / (4 c^2 r^3) = (1/2)(d sigma/dr)^2 has dimension 1/length^2, not an energy "
                          "density; it needs a factor of order c^4/G (or natural units) to be one", strict=True)
def test_rho_grav_is_an_energy_density():
    assert eq(d(G=1, mass=1, c=-2, length=-3), (1, -1, -2))
