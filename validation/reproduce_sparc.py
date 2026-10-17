#!/usr/bin/env python3
"""
reproduce_sparc.py -- an independent check of the QGD rotation-curve claims.

The README quotes R^2 = 0.908, chi^2_nu ~ 1.2 and "zero free parameters per
galaxy" for DarkMatter/data/SPARC-rotational.txt.  This script re-runs the
engine in validation/dark_matter.py on that file and reports

  1. what the file contains (rows, galaxies, datasets), against the README counts;
  2. how many parameters the engine actually has, read from its own constants;
  3. the engine against two baselines that need no tuning beyond one fixed
     constant: baryons only (Newtonian) and MOND with a0 = 1.2e-10 m/s^2;
  4. an out-of-sample test: refit the engine's empirical parameters
     (SIGMA_CRIT, ALPHA, BETA_0) on a random half of the galaxies and score the
     other half, for several seeded splits, against the shipped defaults.

Choices made here, since the repository does not fix them: the total baryonic
mass of a galaxy is M = V_bar(R_last)^2 R_last / G with V_bar from the fixed
mass-to-light ratios in QGDConstants, and the surface density is the engine's
own estimate from V_bar (sigma_star = None).  External field is 0.

Usage:   python validation/reproduce_sparc.py [--splits 5] [--no-cv] [--json out.json]
"""
import argparse
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
warnings.filterwarnings("ignore")
import dark_matter as DM  # noqa: E402  (validation/dark_matter.py)

DATA = os.path.join(HERE, "..", "DarkMatter", "data", "SPARC-rotational.txt")
KPC, MSUN, G_SI = 3.086e19, 1.989e30, 6.674e-11
A0 = 1.2e-10


def load(path=DATA):
    d = pd.read_csv(path)
    d = d[(d.Radius_kpc > 0) & (d.Velocity_error_km_s > 0) & (d.Velocity_km_s > 0)]
    gal = []
    for name, g in d.groupby("Galaxy", sort=True):
        g = g.sort_values("Radius_kpc")
        r = g.Radius_kpc.values
        vg, vd, vb_ = g.V_gas_km_s.values, g.V_disk_km_s.values, g.V_bulge_km_s.values
        vbar = np.sqrt(np.maximum(vg ** 2 + DM.QGDConstants.UPSILON_DISK * vd ** 2
                                  + DM.QGDConstants.UPSILON_BULGE * vb_ ** 2, 0.0))
        mass = max(vbar[-1] ** 2 * 1e6 * r[-1] * KPC / G_SI / MSUN, 1e6)
        gal.append(dict(name=name, r=r, vobs=g.Velocity_km_s.values, err=g.Velocity_error_km_s.values,
                        vg=vg, vd=vd, vb=vb_, vbar=vbar, mass=mass))
    return d, gal


def predict_qgd(gal, engine):
    """V_pred for one galaxy with a given QGDEngine (the same recipe as predict_rotation_curve)."""
    r, vbar = gal["r"], gal["vbar"]
    sigma = np.clip(vbar ** 2 / (r * DM.QGDConstants.G_KPC * 1e6), 0.1, 1000.0)
    r_m = r * KPC
    g_local = G_SI * ((vbar * 1000) ** 2 * r_m / G_SI) / r_m ** 2
    k = np.array([engine.predict_kappa(sigma_star=sigma[i], total_mass=gal["mass"], g_local=g_local[i],
                                       v_circular=vbar[i], radius_kpc=r[i], g_external=0.0)[0]
                  for i in range(len(r))])
    return vbar * np.sqrt(k)


def predict_mond(gal):
    r_m = gal["r"] * KPC
    gb = (gal["vbar"] * 1000) ** 2 / r_m
    gb = np.maximum(gb, 1e-16)
    g = gb / (1.0 - np.exp(-np.sqrt(gb / A0)))
    return np.sqrt(g * r_m) / 1000.0


def score(gals, preds):
    obs = np.concatenate([g["vobs"] for g in gals])
    err = np.concatenate([g["err"] for g in gals])
    pr = np.concatenate(preds)
    res = obs - pr
    ss_tot = np.sum((obs - obs.mean()) ** 2)
    per = [np.sqrt(np.mean((g["vobs"] - p) ** 2)) for g, p in zip(gals, preds)]
    return dict(n=len(obs), R2=float(1 - np.sum(res ** 2) / ss_tot), rms=float(np.sqrt(np.mean(res ** 2))),
                median_galaxy_rms=float(np.median(per)), chi2_nu=float(np.sum((res / err) ** 2) / len(obs)))


def make_engine(p):
    s, a, b = p
    return DM.QGDEngine(sigma_crit=s, alpha=a, beta_0=b)


def fit(gals, p0, maxiter=60):
    def loss(p):
        if p[0] <= 0 or p[1] <= 0 or not (0 <= p[2] <= 2):
            return 1e9
        e = make_engine(p)
        return float(np.sqrt(np.mean(np.concatenate([(g["vobs"] - predict_qgd(g, e)) ** 2 for g in gals]))))
    r = minimize(loss, p0, method="Nelder-Mead", options=dict(maxiter=maxiter, xatol=1e-2, fatol=1e-2))
    return r.x, r.fun


def figure(gals, eng, path):
    """Observed vs predicted velocity (QGD default, MOND) and per-galaxy RMS."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    obs = np.concatenate([g["vobs"] for g in gals])
    q = np.concatenate([predict_qgd(g, eng) for g in gals])
    mo = np.concatenate([predict_mond(g) for g in gals])
    bar = np.concatenate([g["vbar"] for g in gals])
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
    for a, y, lab in ((ax[0], q, "QGD (shipped defaults)"), (ax[1], mo, "MOND, a0 = 1.2e-10")):
        a.loglog(obs, np.maximum(y, 1e-1), ".", ms=2, alpha=0.4)
        a.plot([1, 500], [1, 500], "k", lw=0.8)
        a.set_xlabel("observed V (km/s)"); a.set_ylabel("predicted V (km/s)"); a.set_title(lab)
    rq = [np.sqrt(np.mean((g["vobs"] - predict_qgd(g, eng)) ** 2)) for g in gals]
    rm = [np.sqrt(np.mean((g["vobs"] - predict_mond(g)) ** 2)) for g in gals]
    rn = [np.sqrt(np.mean((g["vobs"] - g["vbar"]) ** 2)) for g in gals]
    for v, lab in ((rn, "baryons only"), (rq, "QGD"), (rm, "MOND")):
        ax[2].plot(np.sort(v), np.linspace(0, 1, len(v)), label="%s (median %.1f)" % (lab, np.median(v)))
    ax[2].set_xscale("log"); ax[2].set_xlabel("per-galaxy RMS (km/s)"); ax[2].set_ylabel("cumulative fraction")
    ax[2].legend(frameon=False); ax[2].set_title("%d galaxies" % len(gals))
    fig.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--splits", type=int, default=5)
    ap.add_argument("--no-cv", action="store_true")
    ap.add_argument("--json", default=None)
    ap.add_argument("--figure", default=None, help="PNG path")
    args = ap.parse_args(argv)

    d, gals = load()
    C = DM.QGDConstants
    out = dict(rows=int(len(d)), galaxies=len(gals), datasets=d.Dataset.value_counts().to_dict())
    out["readme_claims"] = dict(measurements=4248, galaxies=175, R2=0.908, rms=8.3, chi2_nu=1.2)
    out["empirical_parameters"] = dict(SIGMA_CRIT=C.SIGMA_CRIT, ALPHA=C.ALPHA, BETA_0=C.BETA_0,
                                       LOG_M_TRIGGER=C.LOG_M_TRIGGER, UPSILON_DISK=C.UPSILON_DISK,
                                       UPSILON_BULGE=C.UPSILON_BULGE, G_CRIT=C.G_CRIT)

    eng = make_engine((C.SIGMA_CRIT, C.ALPHA, C.BETA_0))
    out["qgd_default"] = score(gals, [predict_qgd(g, eng) for g in gals])
    out["newtonian_baryons"] = score(gals, [g["vbar"] for g in gals])
    out["mond_a0_1.2e-10"] = score(gals, [predict_mond(g) for g in gals])

    if not args.no_cv:
        rng = np.random.default_rng(0)
        cv = []
        for s in range(args.splits):
            idx = rng.permutation(len(gals))
            tr, te = [gals[i] for i in idx[: len(gals) // 2]], [gals[i] for i in idx[len(gals) // 2:]]
            p, _ = fit(tr, [C.SIGMA_CRIT, C.ALPHA, C.BETA_0])
            e2 = make_engine(p)
            cv.append(dict(split=s, fitted=[float(x) for x in p],
                           test_fitted=score(te, [predict_qgd(g, e2) for g in te]),
                           test_default=score(te, [predict_qgd(g, eng) for g in te]),
                           test_mond=score(te, [predict_mond(g) for g in te])))
        out["cross_validation"] = cv

    def show(k, v):
        print("%-22s R2 %.3f  RMS %6.2f km/s  median per-galaxy RMS %6.2f  chi2_nu %6.2f"
              % (k, v["R2"], v["rms"], v["median_galaxy_rms"], v["chi2_nu"]))
    print("file: %d rows, %d galaxies, datasets %s" % (out["rows"], out["galaxies"], out["datasets"]))
    print("README claims: 4248 measurements, 175 galaxies, R2 0.908, RMS 8.3, chi2_nu ~1.2")
    print("empirical parameters in the engine:", out["empirical_parameters"])
    for k in ("qgd_default", "newtonian_baryons", "mond_a0_1.2e-10"):
        show(k, out[k])
    for c in out.get("cross_validation", []):
        print("split %d  fitted (Sigma_crit, alpha, beta0) = %s" % (c["split"], np.round(c["fitted"], 3)))
        show("  test, refit", c["test_fitted"]); show("  test, default", c["test_default"]); show("  test, MOND", c["test_mond"])
    if args.figure:
        figure(gals, eng, args.figure)
    if args.json:
        with open(args.json, "w") as f:
            json.dump(out, f, indent=1, default=float)
    return out


if __name__ == "__main__":
    main()
