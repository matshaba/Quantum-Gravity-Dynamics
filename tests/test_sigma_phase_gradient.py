"""sigma_mu = d_mu S/(mc): is the recipe a phase gradient, and which mass sets the quantum scale?"""
import os
import sys

import pytest
import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "validation"))
import sigma_phase_gradient as sg  # noqa: E402


def test_schwarzschild_recipe_is_not_closed_so_it_is_not_a_phase_gradient():
    w = sg.curl_2d(sg.recipe_one_form(sg.rs / sg.r))
    assert sp.simplify(w + sp.sqrt(sg.rs) / (2 * sg.r ** sp.Rational(3, 2))) == 0
    assert w != 0


def test_any_nonconstant_profile_is_not_closed_and_a_constant_one_is():
    F = sp.Function("F")(sg.r)
    assert sp.simplify(sg.curl_2d(sg.recipe_one_form(F)) - sp.diff(sp.sqrt(F), sg.r)) == 0
    assert sg.curl_2d(sg.recipe_one_form(sp.Integer(3))) == 0


def test_kerr_recipe_is_not_closed():
    assert sg.kerr_curl_at() > 1e-2


def test_exact_free_fall_gradient_is_closed_and_unit_norm_while_the_recipe_is_null():
    assert sp.simplify(sg.curl_2d(sg.hamilton_jacobi_one_form())) == 0
    assert sg.norm(sg.hamilton_jacobi_one_form()) == -1
    assert sg.norm(sg.recipe_one_form(sg.rs / sg.r)) == 0


def test_phase_scale_from_sigma_is_two_times_test_mass_over_source_mass_times_Mcr():
    assert sg.phase_scale_ratio() == 0
    # at r = r_s and m = M this is S_BH/pi = 2 x (S_BH/2pi): a factor two above Sigma ~ Mcr, and proportional to m/M
    ratio = (sg.m / sg.M)
    assert sp.simplify((2 * sg.m * sg.c * sg.rs / sg.hbar).subs(sg.rs, 2 * sg.G * sg.M / sg.c ** 2) /
                       (4 * sp.pi * sg.G * sg.M ** 2 / (sg.c * sg.hbar) / sp.pi) - ratio) == 0
