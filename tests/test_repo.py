"""Repository hygiene: links, the Quick Start, scripts that should run, and recorded known issues."""
import hashlib
import os
import py_compile
import re
import subprocess
import sys
import warnings

import numpy as np
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
README = open(os.path.join(ROOT, "README.md"), encoding="utf-8").read()


def md5(path):
    return hashlib.md5(open(os.path.join(ROOT, path), "rb").read()).hexdigest()


def test_readme_relative_links_exist():
    bad = []
    for m in re.finditer(r"\]\(([^)#\s]+)(#[^)]*)?\)", README):
        link = m.group(1)
        if link.startswith(("http", "mailto")):
            continue
        if not os.path.exists(os.path.join(ROOT, link)):
            bad.append(link)
    assert not bad, bad


def test_readme_quoted_paths_exist():
    paths = set(re.findall(r"`((?:[A-Za-z0-9_\-]+/)*[A-Za-z0-9_\-]+\.(?:py|md|tex|txt))`", README))
    missing = sorted(p for p in paths if not os.path.exists(os.path.join(ROOT, p))
                     and not any(os.path.exists(os.path.join(ROOT, d, p)) for d in ("core", "validation", "docs", "DarkMatter/data")))
    assert not missing, missing


def test_quick_start_python_blocks_run():
    quick = README[README.index("## Quick Start"): README.index("## Key References")]
    blocks = re.findall(r"```python\n(.*?)```", quick, flags=re.S)
    assert blocks
    for code in blocks:
        r = subprocess.run([sys.executable, "-c", code], cwd=ROOT, capture_output=True, text=True, timeout=120)
        assert r.returncode == 0, r.stderr[-500:]


@pytest.mark.parametrize("script", ["EFE_solutions_from_QGD_perspective.py",
                                    "comparison/EFE_solutions_from_QGD_perspective.py"])
def test_efe_script_passes_its_checks(script):
    env = dict(os.environ, MPLBACKEND="Agg")
    r = subprocess.run([sys.executable, os.path.join(ROOT, script)], cwd=os.path.dirname(os.path.join(ROOT, script)),
                       capture_output=True, text=True, timeout=300, env=env)
    assert r.returncode == 0, r.stderr[-500:]
    assert "Passed 10/10 checks" in r.stdout


def _python_files():
    for d, _, fs in os.walk(ROOT):
        if ".git" in d or "__pycache__" in d:
            continue
        for f in fs:
            if f.endswith(".py"):
                yield os.path.relpath(os.path.join(d, f), ROOT)


@pytest.mark.parametrize("path", sorted(p for p in _python_files() if p != "DarkMatter/FullStressEnergyTensor.py"))
def test_python_files_compile_without_warnings(path, tmp_path):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        py_compile.compile(os.path.join(ROOT, path), cfile=str(tmp_path / "x.pyc"), doraise=True)


@pytest.mark.xfail(reason="FullStressEnergyTensor.py is a LaTeX document with a .py extension", strict=True)
def test_full_stress_energy_tensor_is_python(tmp_path):
    py_compile.compile(os.path.join(ROOT, "DarkMatter", "FullStressEnergyTensor.py"),
                       cfile=str(tmp_path / "x.pyc"), doraise=True)


# ------------------------------------------------------------ known issues (see README)
def test_master_metrics_have_lorentzian_signature():
    """Fixed: the diagonal and improved metrics had four negative entries; they are now (-,+,+,+) like g_tt = -(1-f),
    and eta is diag(-1,1,1,1)."""
    sys.path.insert(0, os.path.join(ROOT, "core"))
    import master_metric as mm
    f = mm.sp.Float(0.3)
    for g in (mm.improved_master_metric(f), mm.master_metric_diagonal(f), mm.master_metric_diagonal(f, coords="isotropic"),
              mm.master_metric_diagonal(f, isotropic_grr=False), mm.schwarzschild_gauge_metric(f)):
        signs = [mm.sp.N(g[i, i].subs({mm.r: 5, mm.theta: 1.0})) for i in range(4)]
        assert signs[0] < 0 and all(x > 0 for x in signs[1:])
    assert mm.eta == mm.sp.diag(-1, 1, 1, 1)


def _ricci(g, X):
    import sympy as sp
    sys.path.insert(0, os.path.join(ROOT, "validation"))
    import nariai_entropy as ne
    return ne.curvature(g, X)[:2]


def test_areal_gauge_metric_is_GR_and_the_improved_metric_is_not():
    """schwarzschild_gauge_metric(f) is vacuum for f = 2GM/c^2 r, solves the Einstein equations with Lambda for
    f = r_s/r + H^2 r^2, and has the Reissner-Nordstrom Ricci scalar zero; the 'improved' g_rr, which mixes areal g_tt with
    an isotropic-gauge g_rr, is not Ricci flat for the same f."""
    import sympy as sp
    sys.path.insert(0, os.path.join(ROOT, "core"))
    import master_metric as mm
    r, th, ph, t = mm.r, mm.theta, sp.Symbol("phi"), sp.Symbol("t")
    X = [t, r, th, ph]
    rs, H = sp.symbols("r_s H", positive=True)
    Ric, Rs = _ricci(mm.schwarzschild_gauge_metric(rs / r, r, th), X)
    assert all(sp.simplify(Ric[i, j]) == 0 for i in range(4) for j in range(4))
    g = mm.schwarzschild_gauge_metric(rs / r + H ** 2 * r ** 2, r, th)
    Ric, Rs = _ricci(g, X)
    assert all(sp.simplify(Ric[i, j] - 3 * H ** 2 * g[i, j]) == 0 for i in range(4) for j in range(4))
    Q = sp.symbols("Q", positive=True)
    Ric, Rs = _ricci(mm.schwarzschild_gauge_metric(rs / r - Q ** 2 / r ** 2, r, th), X)
    assert sp.simplify(Rs) == 0
    Ric, Rs = _ricci(mm.improved_master_metric(rs / r, r, th), X)
    assert any(sp.simplify(Ric[i, j]) != 0 for i in range(4) for j in range(4))


@pytest.mark.xfail(reason="the 'improved' metric does not equal isotropic Schwarzschild (g_tt differs)", strict=True)
def test_verify_classical_limit_matches():
    sys.path.insert(0, os.path.join(ROOT, "core"))
    import master_metric as mm
    assert mm.verify_classical_limit(mm.f_scalar())["match"]


@pytest.mark.xfail(reason="Ch9-QG_QFT.tex is a byte-for-byte copy of Ch8-GravitonApplications.tex", strict=True)
def test_chapters_8_and_9_differ():
    assert md5("docs/Ch8-GravitonApplications.tex") != md5("docs/Ch9-QG_QFT.tex")


@pytest.mark.xfail(reason="the four notebooks have identical content", strict=True)
def test_notebooks_differ():
    names = os.listdir(os.path.join(ROOT, "notebooks"))
    assert len({md5("notebooks/" + n) for n in names}) == len(names)


def test_claims_ledger_is_current_and_points_at_real_tests():
    sys.path.insert(0, os.path.join(ROOT, "validation"))
    import claims_ledger as CL
    assert open(CL.OUT, encoding="utf-8").read() == CL.render()
    for claim, where, status, test in CL.CLAIMS:
        if not test:
            continue
        path, name = test.split("::")
        assert re.search(r"def %s\(" % re.escape(name), open(os.path.join(ROOT, path), encoding="utf-8").read()), test
