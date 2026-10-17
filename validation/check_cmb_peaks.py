#!/usr/bin/env python3
"""
check_cmb_peaks.py -- test the QGD acoustic-peak claim against the Planck binned TT data.

The claim (validation/dark_matter.py, CMBPredictor): l_n = A kappa_4 n with kappa_4 = 8.8741 and
A = 31.51 "fitted to the first peak", with errors below 1%.

This script takes the 'Planck binned' rows of DarkMatter/data/PLANCK-TT-Spectrum.txt, locates the
first acoustic peaks by a local parabola through the top bins, and compares

  (a) the QGD prediction l_n = 31.51 * 8.8741 * n,
  (b) the same law with A refitted to all measured peaks (one parameter), and
  (c) a two-parameter line l_n = a n + b, the usual phenomenological form.

The location error of a peak is estimated by bootstrapping the bins within their quoted errors.
Usage:  python validation/check_cmb_peaks.py [--figure out.png]
"""
import argparse
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "DarkMatter", "data", "PLANCK-TT-Spectrum.txt")
KAPPA4, A_QGD = 8.8741, 31.51
WINDOWS = [(150, 300), (450, 620), (740, 900), (1050, 1200), (1330, 1500)]   # where peaks 1..5 sit


def load_planck_binned(path=DATA):
    d = pd.read_csv(path, comment="#", skipinitialspace=True)
    d.columns = [c.strip() for c in d.columns]
    d = d[d["Experiment"].astype(str).str.strip() == "Planck binned"]
    for col in ("l_center", "Power", "Sigma_minus", "Sigma_plus"):
        d[col] = pd.to_numeric(d[col], errors="coerce")
    d = d.dropna(subset=["l_center", "Power"]).sort_values("l_center")
    sig = 0.5 * (d["Sigma_minus"].fillna(0) + d["Sigma_plus"].fillna(0)).values
    return d["l_center"].values, d["Power"].values, np.maximum(sig, 1e-6)


def peak_location(l, p, lo, hi, top=5):
    m = (l >= lo) & (l <= hi)
    ll, pp = l[m], p[m]
    idx = np.argsort(pp)[-top:]
    sel = np.sort(idx)
    x, y = ll[sel], pp[sel]
    c = np.polyfit(x, y, 2)
    return -c[1] / (2 * c[0]) if c[0] < 0 else float(ll[np.argmax(pp)])


def peaks(l, p, sig, n_boot=300, seed=0):
    rng = np.random.default_rng(seed)
    loc = np.array([peak_location(l, p, lo, hi) for lo, hi in WINDOWS])
    boots = np.array([[peak_location(l, p + rng.normal(0, sig), lo, hi) for lo, hi in WINDOWS]
                      for _ in range(n_boot)])
    return loc, boots.std(axis=0)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--figure", default=None)
    args = ap.parse_args(argv)
    l, p, sig = load_planck_binned()
    obs, err = peaks(l, p, sig)
    n = np.arange(1, len(obs) + 1)
    qgd = A_QGD * KAPPA4 * n
    a1 = float(np.sum(obs * n) / np.sum(n * n))                      # l = a n
    lin = np.polyfit(n, obs, 1)
    res = dict(n=n, observed=obs, observed_err=err, qgd=qgd, refit_one=a1 * n, line=np.polyval(lin, n),
               a_qgd=A_QGD * KAPPA4, a_refit=a1, line_coeffs=lin)
    print("peak   observed   +-    QGD(279.6 n)  err%   one-par refit   err%   line a n + b   err%")
    for i in range(len(n)):
        print("%2d   %8.1f  %5.1f   %9.1f  %6.1f   %9.1f  %6.1f   %9.1f  %6.1f" % (
            n[i], obs[i], err[i], qgd[i], 100 * (qgd[i] / obs[i] - 1), res["refit_one"][i],
            100 * (res["refit_one"][i] / obs[i] - 1), res["line"][i], 100 * (res["line"][i] / obs[i] - 1)))
    print("line fit: l_n = %.1f n + %.1f   (a = spacing, b = phase offset)" % tuple(lin))
    print("QGD constant A*kappa_4 = %.1f; the first peak is at %.0f, the mean spacing is %.0f"
          % (A_QGD * KAPPA4, obs[0], np.mean(np.diff(obs))))
    if args.figure:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].plot(l, p, ".", ms=3, color="0.5")
        for k in range(len(obs)):
            ax[0].axvline(obs[k], color="C0", lw=0.8)
            ax[0].axvline(qgd[k], color="C3", lw=0.8, ls="--")
        ax[0].set_xlim(0, 1800); ax[0].set_xlabel("l"); ax[0].set_ylabel("l(l+1)C_l/2pi")
        ax[0].set_title("Planck binned TT; blue: measured peaks, red dashed: QGD l_n = 279.6 n")
        ax[1].errorbar(n, obs - qgd, yerr=err, fmt="o", label="measured - QGD")
        ax[1].errorbar(n, obs - res["line"], yerr=err, fmt="s", label="measured - line fit")
        ax[1].axhline(0, color="k", lw=0.6); ax[1].set_xlabel("peak n"); ax[1].set_ylabel("l")
        ax[1].legend(frameon=False)
        fig.tight_layout(); fig.savefig(args.figure, dpi=140)
    return res


if __name__ == "__main__":
    main()
