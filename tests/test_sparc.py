"""Rotation-curve claims, checked with validation/reproduce_sparc.py (no cross-validation: fast)."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "validation"))
import reproduce_sparc as RS  # noqa: E402


@pytest.fixture(scope="module")
def res():
    return RS.main(["--no-cv"])


def test_pooled_r2_reproduces_readme(res):
    assert res["qgd_default"]["R2"] == pytest.approx(0.908, abs=0.01)


def test_engine_beats_baryons_only(res):
    assert res["qgd_default"]["rms"] < 0.6 * res["newtonian_baryons"]["rms"]


@pytest.mark.xfail(reason="README quotes RMS 8.3 km/s; this recipe gives ~26 km/s on the shipped file", strict=False)
def test_readme_rms(res):
    assert res["qgd_default"]["rms"] < 10.0


@pytest.mark.xfail(reason="README quotes chi2_nu ~ 1.2; this recipe gives ~150 on the shipped file", strict=False)
def test_readme_chi2(res):
    assert res["qgd_default"]["chi2_nu"] < 2.0


@pytest.mark.xfail(reason="README says 4,248 measurements / 175 galaxies; the shipped file has 3,827 rows / 225 galaxies",
                   strict=False)
def test_readme_counts(res):
    assert res["rows"] == 4248 and res["galaxies"] == 175


def test_mond_one_constant_baseline_is_reported(res):
    # not a claim about physics: records that the comparison runs and what it returned
    assert 0.5 < res["mond_a0_1.2e-10"]["R2"] < 1.0


def test_engine_has_empirical_parameters(res):
    """The README says 'zero free parameters'; the engine's own constants say otherwise (see QGDConstants)."""
    emp = res["empirical_parameters"]
    assert len([k for k in emp if k in ("SIGMA_CRIT", "ALPHA", "BETA_0", "LOG_M_TRIGGER")]) == 4
