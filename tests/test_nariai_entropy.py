"""
Entropy correction vs the Nariai / de Sitter pair (G = hbar = c = 1; see validation/nariai_entropy.py).

Each test states what it establishes.  Where the outcome depends on how a QGD formula is read, both readings
are tested and the dependence is the finding.
"""
import sys
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))
import nariai_entropy as ne  # noqa: E402


@pytest.fixture(scope="module")
def euclid():
    return ne.euclidean_data()


# ---- GR baseline by two independent routes -------------------------------------------------------------------
def test_euclidean_action_entropies(euclid):
    """S = -I: de Sitter (S^4) 3 pi / Lambda, Nariai (S^2 x S^2) 2 pi / Lambda; both are Einstein spaces."""
    assert euclid["S4"]["einstein"] and euclid["S2xS2"]["einstein"]
    assert sp.simplify(euclid["S4"]["S_over_inv_Lam"] - 3 * sp.pi) == 0
    assert sp.simplify(euclid["S2xS2"]["S_over_inv_Lam"] - 2 * sp.pi) == 0


def test_area_route_agrees_with_euclidean_route():
    """Horizon areas of the sigma_Lambda recipe: dS one horizon 12 pi/Lambda, Nariai two coincident 4 pi/Lambda."""
    Lam = 1.0
    r_dS = ne.horizons(0.0, 0.0, Lam)[-1]
    rN, MN = ne.nariai_with_A(0.0, Lam)
    assert np.pi * r_dS ** 2 == pytest.approx(3 * np.pi / Lam, rel=1e-12)          # A/4 = 12 pi/(4 Lam)
    assert 2 * np.pi * rN ** 2 == pytest.approx(2 * np.pi / Lam, rel=1e-9)         # two horizons of 4 pi/Lambda
    assert MN == pytest.approx(1 / (3 * np.sqrt(Lam)), rel=1e-9)                   # G M = c^3 / (3 sqrt3 H)
    assert 2 * np.pi * rN ** 2 - np.pi * r_dS ** 2 == pytest.approx(-np.pi / Lam, rel=1e-9)   # dS -> Nariai: -S_dS/3


# ---- what "topological" would mean ---------------------------------------------------------------------------
def test_euler_characteristics_count_horizons(euclid):
    """Int GB / 32 pi^2 = chi: 2 for S^4 (one horizon), 4 for S^2 x S^2 (two).  A term proportional to chi is
    therefore proportional to the number of horizons."""
    assert euclid["S4"]["chi"] == 2
    assert euclid["S2xS2"]["chi"] == 4


def test_topological_constant_is_invisible_to_dS_dM_but_not_to_pair_creation():
    """CONDITIONAL on the premise that the QGD correction is C * chi/2 (the 'topological' reading).  Then
    d/dM at fixed horizon number is blind to C (as tests/test_gaps.py shows) but the dS -> Nariai entropy
    difference shifts by exactly +C, i.e. the pair-creation rate carries a factor exp(+C)."""
    C, Lam = sp.symbols("C Lambda", positive=True)
    S_dS = 3 * sp.pi / Lam + C * 2 / 2
    S_N = 2 * sp.pi / Lam + C * 4 / 2
    assert sp.simplify((S_N - S_dS) - (-sp.pi / Lam) - C) == 0


# ---- the quantum term as a metric term ------------------------------------------------------------------------
def test_constant_A_area_shift_is_2piA_which_is_pi_times_the_claimed_constant():
    """Reading I: A = kappa l_Q^2 constant.  Outer-horizon area entropy is 4 pi M^2 + 2 pi A + O(A^2/M^2).
    The documented constant is 2 kappa l_Q^2 = 2 A: they differ by a factor pi (normalisation of l_Q in f, or of
    the Gauss-Bonnet coefficient, is not given in the repository)."""
    A = 0.05
    M = 100.0
    shift = ne.area_entropy(M, A) - 4 * np.pi * M ** 2
    assert shift == pytest.approx(2 * np.pi * A, rel=1e-4)
    assert shift / (2 * A) == pytest.approx(np.pi, rel=1e-3)


def test_constant_A_makes_a_horizon_at_zero_mass():
    """With A > 0 constant, h = 1 - A/r^2 vanishes at r = sqrt(A) even for M = 0: the 'quantum horizon' of the
    negative-Q^2 Reissner-Nordstrom form.  Flat space is no longer the M -> 0 limit."""
    for A in (1e-4, 1e-2):
        r = ne.horizons(0.0, A)
        assert len(r) == 1 and r[0] == pytest.approx(np.sqrt(A), rel=1e-9)


def test_mass_dependent_A_is_not_constant_and_breaks_the_first_law():
    """Reading II (README): A = kappa hbar^2 / (M^2 c^2) -> kappa / M^2.  The area shift is 2 pi kappa / M^2 (it
    falls with M, so it is not the documented constant) and dM = T dS_area fails because A varies with M."""
    kappa = 1.0
    shifts = [(M, ne.area_entropy(M, kappa / M ** 2) - 4 * np.pi * M ** 2) for M in (10.0, 20.0, 40.0)]
    for M, s in shifts:
        assert s * M ** 2 == pytest.approx(2 * np.pi * kappa, rel=2e-3)
    assert shifts[0][1] > 3.9 * shifts[1][1] > 15 * shifts[2][1]            # ~ 1/M^2
    M0, h = 5.0, 1e-5
    Af = lambda M: kappa / M ** 2
    dSdM = (ne.area_entropy(M0 + h, Af(M0 + h)) - ne.area_entropy(M0 - h, Af(M0 - h))) / (2 * h)
    assert abs(dSdM * ne.temperature(M0, Af(M0)) - 1) > 1e-4                 # first law violated
    # and for constant A it holds
    dSdM = (ne.area_entropy(M0 + h, 0.05) - ne.area_entropy(M0 - h, 0.05)) / (2 * h)
    assert dSdM * ne.temperature(M0, 0.05) == pytest.approx(1.0, abs=1e-7)


def test_metric_term_shift_is_not_topological_it_changes_sign_on_the_cosmological_horizon():
    """At Lambda = 1, M = 0, constant A: cosmological horizon entropy shifts by -pi A, the small quantum horizon
    carries +pi A.  A term proportional to chi would shift every horizon by the same amount; this does not."""
    A, Lam = 1e-4, 1.0
    rq, rc = ne.horizons(0.0, A, Lam)
    assert np.pi * rc ** 2 - 3 * np.pi / Lam == pytest.approx(-np.pi * A, rel=2e-3)
    assert np.pi * rq ** 2 == pytest.approx(+np.pi * A, rel=2e-3)


def test_nariai_point_moves_with_A():
    """Double root of r^2 h with the quantum term: r_N = 1/sqrt(Lambda) + A/2 and G M_N = 1/(3 sqrt Lambda) - A/2
    to first order (Lambda = 1), so the Nariai mass bound itself is shifted by the term."""
    A = 1e-4
    rN, MN = ne.nariai_with_A(A, 1.0)
    assert rN - 1 == pytest.approx(A / 2, rel=1e-2)
    assert MN - 1 / 3 == pytest.approx(-A / 2, rel=1e-2)


# ---- Gauss-Bonnet reading --------------------------------------------------------------------------------------
def test_gauss_bonnet_entropy_two_routes_and_what_it_fixes(euclid):
    """L = (R + alpha GB)/16 pi.  Euclidean route: S = -I = 2 pi alpha chi.  Jacobson-Myers horizon route:
    S = (1/4) Int_H (1 + 2 alpha R_H) -> alpha * 4 pi per spherical horizon (Int R_H = 8 pi).  They agree on S^4
    (one horizon) and S^2 x S^2 (two).  Matching the documented constant 2 kappa l_Q^2 needs
    alpha = kappa l_Q^2 / (2 pi); and since GB is topological in 4D it adds nothing to the field equations, so it
    cannot generate the kappa l_Q^2 / r^2 metric term: the entropy constant and the f term are independent inputs."""
    alpha, kap = sp.symbols("alpha kappa_lQ2", positive=True)
    for name, n_h in (("S4", 1), ("S2xS2", 2)):
        euclid_route = 2 * sp.pi * alpha * euclid[name]["chi"]
        jm_route = n_h * sp.Rational(1, 4) * 2 * alpha * 8 * sp.pi
        assert sp.simplify(euclid_route - jm_route) == 0
    one_horizon = 4 * sp.pi * alpha
    assert sp.solve(sp.Eq(one_horizon, 2 * kap), alpha)[0] == kap / (2 * sp.pi)
