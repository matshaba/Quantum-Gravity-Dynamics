#!/usr/bin/env python3
"""
claims_ledger.py -- which QGD claims are verified, which fail, and which cannot be tested here.

The same discipline as the exactness ledger of the CICY / spectre projects: each claim gets one status,
and the status points at the test that decides it.

  VERIFIED     checked independently (sympy / exact arithmetic / numerics) and true as stated
  RECORDED     true but with a stated limitation of what it shows
  FAILS        checked, and false as stated (test is an expected failure)
  UNTESTED     needs input this repository does not contain

`python validation/claims_ledger.py` writes docs/claims_ledger.md.  tests/test_repo.py checks that every test
named here exists and that the written file is current.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "docs", "claims_ledger.md")

CLAIMS = [
    # (claim, where it is made, status, test)
    ("kappa_j = sqrt((2j-1)!/2^(2j-2)) matches the constants used by the engine", "README, validation/dark_matter.py",
     "VERIFIED", "tests/test_claims.py::test_kappa_values_in_engine_constants"),
    ("kappa_7 = 1233 (README previously said 1245)", "README", "VERIFIED", "tests/test_claims.py::test_readme_kappa_list_matches_formula"),
    ("|alpha_G|^2 identities; = 1/2 at the Planck mass", "key-equations.md", "VERIFIED", "tests/test_claims.py::test_planck_mass_gives_one_half"),
    ("Earth-electron |alpha_G|^2 = 1.46e39", "key-equations.md, GR-vs-QGD.md, complete_paper.tex (old)", "FAILS",
     "tests/test_claims.py::test_alpha_identities"),
    ("Schwarzschild recipe (sigma_t = sqrt(2GM/c^2 r)) is vacuum, in Kerr-Schild form", "README", "VERIFIED",
     "tests/test_claims.py::test_schwarzschild_recipe_is_vacuum"),
    ("Reissner-Nordstrom recipe (eps = -1) carries exactly the Maxwell stress", "README", "VERIFIED",
     "tests/test_claims.py::test_reissner_nordstrom_recipe_has_em_stress"),
    ("sigma_Lambda recipe reaches the Nariai double horizon, area 4 pi / Lambda", "README (de Sitter row)", "VERIFIED",
     "tests/test_gaps.py::test_sigma_lambda_recipe_reaches_the_nariai_limit"),
    ("Universal f is additive: terms add in the Einstein tensor (static metrics)", "README", "RECORDED",
     "tests/test_gaps.py::test_universal_f_is_additive_in_the_einstein_tensor"),
    ("Sigma-sum N-body metric is an exact (vacuum) solution", "README, QGD_vs_GR.py", "FAILS",
     "tests/test_gaps.py::test_two_body_superposition_is_not_vacuum"),
    ("Sigma-addition is exact for terms sharing one null direction (mass + charge on one centre)", "implied", "VERIFIED",
     "tests/test_gaps.py::test_common_null_direction_adds_exactly_in_mixed_ricci"),
    ("Quantum term +kappa hbar^2/(M^2 c^2 r^2) in f has positive energy / resolves the singularity", "README, master_metric.py", "FAILS",
     "tests/test_gaps.py::test_adding_the_quantum_term_does_not_remove_the_singularity"),
    ("Merger at d = 4 r_s", "README, QGD_vs_GR.py", "RECORDED", "tests/test_gaps.py::test_merger_condition_with_total_sigma_gives_eight_rs"),
    ("Constant entropy correction 2 kappa l_Q^2/(G hbar) for micro black holes", "GR-vs-QGD.md", "RECORDED",
     "tests/test_gaps.py::test_constant_entropy_correction_drops_out_of_thermodynamic_differences"),
    ("Entropy: de Sitter 3 pi/Lambda and Nariai 2 pi/Lambda agree between the Euclidean action and horizon areas", "implied (Nariai check)",
     "VERIFIED", "tests/test_nariai_entropy.py::test_area_route_agrees_with_euclidean_route"),
    ("A topological (chi-proportional) entropy constant is observable only through horizon-number-changing processes (dS -> Nariai)",
     "refines the row above", "RECORDED",
     "tests/test_nariai_entropy.py::test_topological_constant_is_invisible_to_dS_dM_but_not_to_pair_creation"),
    ("Entropy correction is the constant 2 kappa l_Q^2/(G hbar) from the metric term kappa l_Q^2/r^2", "GR-vs-QGD.md, README", "FAILS",
     "tests/test_nariai_entropy.py::test_constant_A_area_shift_is_2piA_which_is_pi_times_the_claimed_constant"),
    ("Entropy correction is mass independent when the metric term is kappa hbar^2/(M^2 c^2 r^2)", "GR-vs-QGD.md, README", "FAILS",
     "tests/test_nariai_entropy.py::test_mass_dependent_A_is_not_constant_and_breaks_the_first_law"),
    ("Entropy correction is topological (the same on every horizon)", "GR-vs-QGD.md", "FAILS",
     "tests/test_nariai_entropy.py::test_metric_term_shift_is_not_topological_it_changes_sign_on_the_cosmological_horizon"),
    ("Quantum term kappa l_Q^2 d(sigma)d(sigma) of the master metric equals kappa hbar^2/(M^2 c^2 r^2) in g_tt", "README, master_metric.py", "FAILS",
     "tests/test_quantum_core.py::test_quantum_term_acts_on_the_angular_part_only"),
    ("Entropy correction is the Gauss-Bonnet / Wald entropy of the l_Q^2 (nabla nabla sigma)^2 term of the action", "GR-vs-QGD.md", "FAILS",
     "tests/test_quantum_core.py::test_wald_entropy_of_the_l2_term_vanishes_for_any_F"),
    ("Literal master-metric quantum term shifts S by +-pi kappa l_Q^2 on every horizon: constant, first-law consistent", "derived here", "RECORDED",
     "tests/test_quantum_core.py::test_entropy_shift_is_the_same_on_every_horizon_and_consistent_with_the_first_law"),
    ("Literal master-metric quantum term on Kerr: still angular-only, but moves the horizon; Delta S = pi s kappa l_Q^2 / sqrt(1 - a^2/M^2), so the constant is not universal for rotating sources", "derived here", "RECORDED",
     "tests/test_kerr_quantum_term.py::test_entropy_shift_is_pi_s_b_over_sqrt_one_minus_chi_squared"),
    ("Correspondence principle lim hbar->0 g_QGD = g_GR, corrections O(hbar^2) significant at the Compton scale: read as the high-entropy limit, the expansion parameter at the horizon is (2 pi/S_BH)^2 and the Gaussian-state entropy correction is also O(hbar^2)", "key-equations.md, complete_paper.tex", "RECORDED",
     "tests/test_classical_limit_entropy.py::test_expansion_parameter_at_horizon_is_two_pi_over_S_BH_squared"),
    ("The recipe sigma_mu equals d_mu S/(mc) for a single-valued WKB phase S (and carries a geometric phase oint sigma dx)", "key-equations.md, insights-of-QGD.md", "FAILS",
     "tests/test_sigma_phase_gradient.py::test_schwarzschild_recipe_is_not_closed_so_it_is_not_a_phase_gradient"),
    ("Wald entropy of the l_Q^2 term on Kerr is nonzero and spin dependent: (pi^2/3) s l_Q^2 (a/M)^2 at small spin", "derived here", "RECORDED",
     "tests/test_quantum_core.py::test_kerr_small_spin_coefficient_is_pi_squared_over_three"),
    ("Wald term on Kerr-de Sitter: small-spin S = -(4 pi^2/3) s l_Q^2 a^2 sigma_M^2(r_h) F'(r_h)/r_h (opposite signs on the two horizons), a^4 at Nariai", "derived here", "RECORDED",
     "tests/test_quantum_core.py::test_spin_dependence_survives_at_nariai_but_starts_at_a_to_the_fourth"),
    ("Literal master-metric quantum term removes the singularity", "README", "FAILS",
     "tests/test_quantum_core.py::test_quantum_term_does_not_remove_the_singularity"),
    ("Schwarzschild sigma satisfies Box sigma_mu = Q_mu + G_mu + T_mu in vacuum", "README, complete_paper.tex", "UNTESTED",
     "tests/test_quantum_core.py::test_schwarzschild_sigma_does_not_solve_box_sigma_equals_Q_with_G_zero"),
    ("CMB peaks l_n = A kappa_4 n within 1%, A fitted to the first peak", "validation/dark_matter.py", "FAILS",
     "tests/test_gaps.py::test_cmb_claim_less_than_one_percent_is_not_reproduced"),
    ("SPARC: R^2 = 0.908", "README", "VERIFIED", "tests/test_sparc.py::test_pooled_r2_reproduces_readme"),
    ("SPARC: RMS 8.3 km/s, chi2_nu ~ 1.2", "README", "FAILS", "tests/test_sparc.py::test_readme_rms"),
    ("SPARC: zero free parameters", "README", "FAILS", "tests/test_sparc.py::test_engine_has_empirical_parameters"),
    ("ell_Q = sqrt(G hbar^2 / c^4) is a length", "README, master_metric.py", "FAILS", "tests/test_dimensions.py::test_ell_Q_as_written_is_a_length"),
    ("rho_grav = G M / (4 c^2 r^3) is an energy density", "QGD_vs_GR.py", "FAILS", "tests/test_dimensions.py::test_rho_grav_is_an_energy_density"),
    ("Improved master metric is isotropic Schwarzschild at l_Q -> 0", "master_metric.py", "FAILS", "tests/test_repo.py::test_verify_classical_limit_matches"),
    ("Neutron-star mass shift ~0.01 Msun (quantum TOV)", "README", "UNTESTED", ""),
    ("Fourth-order field equation reduces to Einstein's equations as l_Q -> 0", "README", "UNTESTED", ""),
    ("Kerr, Vaidya and interior solutions from the recipe", "README", "UNTESTED", ""),
]

ORDER = ["VERIFIED", "RECORDED", "FAILS", "UNTESTED"]


def render():
    lines = ["# QGD claims ledger", "",
             "Generated by `validation/claims_ledger.py`; each status points at the test that decides it.", "",
             "| status | claim | made in | test |", "|---|---|---|---|"]
    for st in ORDER:
        for claim, where, s, test in CLAIMS:
            if s == st:
                lines.append("| %s | %s | %s | %s |" % (s, claim, where, ("`%s`" % test) if test else "—"))
    counts = {s: sum(1 for c in CLAIMS if c[2] == s) for s in ORDER}
    lines += ["", "Counts: " + ", ".join("%s %d" % (s, counts[s]) for s in ORDER) + ".", "",
              "VERIFIED: independently checked and true as stated. RECORDED: true, with a limitation on what it shows. "
              "FAILS: checked and false as stated (an expected failure in the test suite). UNTESTED: needs input this "
              "repository does not contain.", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    with open(OUT, "w", encoding="utf-8") as f:
        f.write(render())
    print("wrote docs/claims_ledger.md")
