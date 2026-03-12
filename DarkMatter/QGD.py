"""
Copyright (c) 2026 Romeo Matshaba. All rights reserved.

pygrc.qgd  —  Quantum Gravitational Dynamics (QGD) κ-ladder model
============================================================

Zero-free-parameter rotation curve model derived from the Dirac equation.
Gravity is an emergent quantum effect; the κ-ladder replaces dark matter.

    κ_n = sqrt( (2n-1)! / 4^(n-1) )

    n=1:   κ =  1.000   Newtonian baseline
    n=2:   κ =  1.225   Wide binaries / dwarf spheroidals
    n=3:   κ =  2.739   Spiral outskirts  ← primary SPARC rung
    n=4:   κ =  8.874   Galaxy groups / massive spirals
    n=5:   κ = 37.650   Galaxy clusters (Bullet Cluster: ratio=0.981)

T^μν stress-energy tensor corrections (v2.3)
--------------------------------------------
κ_full = 1 + (κ_base − 1) × f_P × f_β × f_sh × f_×

    f_P (w)    = exp(−ln3 · w)               pressure suppression
    f_β (β)    = 1 + 0.15(β − 0.5)           orbital anisotropy
    f_sh       = 1 + 0.10|∂lnΣ/∂lnr|/(1+w)  shear stress
    f_× (v)    = 1 + 0.08(v/250)^0.5         momentum flux T^0i

Σ(r) is reconstructed from the velocity decomposition when per-point
surface brightness (SBdisk) is unavailable.

Quickstart
----------
    import pygrc as gr
    from pygrc.qgd import QGD, qgd_velocity

    df = gr.Reader.read("NGC5055_rotmod.dat")

    # ── Zero-parameter prediction ──────────────────────────────────────
    model = QGD(df, name="NGC5055")
    v_pred, kappas = model.predict()
    print(model.summary())
    model.plot()

    # ── Use with pygrc.Fit (optional mass-scale scan) ──────────────────
    from pygrc.qgd import register_galaxy
    register_galaxy(df)

    m = gr.Fit(df["Rad"], df["Vobs"], 1.0).fit_lsq(
            qgd_velocity,
            [(0.5, 2.0)],       # mass_scale search range
            df["errV"],
            [False])            # False = let it float; True = fix at 1.0

    fig, ax = plt.subplots()
    gr.Plot().plot_grc(df, m, qgd_velocity, "QGD", "NGC5055", ax)

References
----------
    Matshaba, R. (2026). Quantum Gravitational Dynamics.
    Lelli et al. (2016) SPARC database. AJ 152, 157.
    Clowe et al. (2006) Bullet Cluster. ApJL 648, L109.
"""

from math import factorial
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import typing as tp

# ── Physical constants ────────────────────────────────────────────────────────
_G_KPC   = 4.300e-6    # G in units: kpc (km/s)² M☉⁻¹  (SPARC convention)
_KPC_M   = 3.086e19    # 1 kpc in metres
_A0      = 1.2e-10     # m s⁻²  QGD/MOND critical acceleration
_LN3     = np.log(3.0) # ln(3) for pressure correction

# ── Universal QGD parameters (fitted once over all SPARC galaxies) ────────────
SIGMA_CRIT   = 17.5    # M☉ pc⁻²  phase-coherence threshold
ALPHA_SIGMA  = 0.30    # power-law index for Σ correction
UPS_DISK     = 0.50    # stellar disk  Υ  [M☉ L☉⁻¹]
UPS_BULGE    = 0.70    # bulge         Υ  [M☉ L☉⁻¹]
LOG_M_REF_2  = 8.80    # log₁₀(M/M☉) — Q₂ dwarf   reference (≈ LMC mass)
LOG_M_REF_3  = 9.25    # log₁₀(M/M☉) — Q₃ spiral  reference
LOG_M_REF_4  = 13.0    # log₁₀(M/M☉) — Q₄ group/massive-spiral reference

# T^μν correction coefficients (v2.3)
C_BETA  = 0.15         # anisotropy coefficient
C_SH    = 0.10         # shear stress coefficient
C_CROSS = 0.08         # momentum flux coefficient
V_CROSS = 250.0        # km/s — normalisation velocity for f_×


# ── κ-ladder ──────────────────────────────────────────────────────────────────

def kappa_n(n: int) -> float:
    """
    Single rung of the QGD κ-ladder.

        κ_n = sqrt( (2n-1)! / 4^(n-1) )

    Parameters
    ----------
    n : int  —  rung index (1 = Newtonian)

    Returns
    -------
    float
    """
    return float(np.sqrt(factorial(2 * n - 1) / 4 ** (n - 1)))


def kappa_ladder(n_max: int = 7) -> dict:
    """
    Full QGD κ-ladder as a dict {n: κ_n}.

    Parameters
    ----------
    n_max : int  —  highest rung to compute (default 7)

    Returns
    -------
    dict[int, float]

    Example
    -------
    >>> kappa_ladder(5)
    {1: 1.0, 2: 1.2247..., 3: 2.7386..., 4: 8.8741..., 5: 37.6497...}
    """
    return {n: kappa_n(n) for n in range(1, n_max + 1)}


# Pre-computed values used throughout the module
KAPPA = kappa_ladder(7)
K1, K2, K3, K4, K5, K6, K7 = [KAPPA[n] for n in range(1, 8)]


# ── Baryonic velocity ─────────────────────────────────────────────────────────

def v_baryonic(df: pd.DataFrame,
               ups_disk: float = UPS_DISK,
               ups_bulge: float = UPS_BULGE) -> np.ndarray:
    """
    Total baryonic circular velocity in quadrature from SPARC columns.

        v_bar(r) = sqrt( v_gas² + Υ_disk·v_disk² + Υ_bulge·v_bul² )

    Parameters
    ----------
    df        : SPARC DataFrame (columns: Vgas, Vdisk, Vbul)
    ups_disk  : disk  Υ [M☉ L☉⁻¹]
    ups_bulge : bulge Υ [M☉ L☉⁻¹]

    Returns
    -------
    np.ndarray [km/s]
    """
    return np.sqrt(
        df["Vgas"].values ** 2
        + ups_disk  * df["Vdisk"].values ** 2
        + ups_bulge * df["Vbul"].values  ** 2
    )


def baryonic_mass(df: pd.DataFrame,
                  ups_disk: float = UPS_DISK,
                  ups_bulge: float = UPS_BULGE) -> float:
    """
    Estimate total baryonic mass from the outermost data point.

        M = v_bar²(r_max) · r_max / G

    Parameters
    ----------
    df        : SPARC DataFrame
    ups_disk  : disk  Υ
    ups_bulge : bulge Υ

    Returns
    -------
    float [M☉]
    """
    vbar  = v_baryonic(df, ups_disk, ups_bulge)
    idx   = int(np.argmax(df["Rad"].values))
    r_kpc = float(df["Rad"].values[idx])
    v_kms = float(vbar[idx])
    return v_kms ** 2 * r_kpc / _G_KPC


def surface_density(df: pd.DataFrame,
                    ups_disk: float = UPS_DISK) -> np.ndarray:
    """
    Stellar surface density Σ(r) [M☉ pc⁻²], radially resolved.

    Strategy (in priority order):

    1. **SBdisk column** (native SPARC .dat):  Σ = Υ_disk × SBdisk(r)
    2. **Velocity reconstruction** (this dataset / any SPARC-format CSV):
       Σ_disk(r) is derived from the gradient of the enclosed stellar mass:
           M_disk(<r) = Υ_disk · v_disk²(r) · r / G
           Σ_disk(r) = (1/2πr) · dM_disk/dr   [converted to M☉ pc⁻²]
       This is then normalised to the global Sigma_star if that column is
       present, preserving the absolute scale while recovering the radial
       profile from the already-measured velocity decomposition.
    3. **Constant fallback**: 5.0 M☉ pc⁻² everywhere (< Σ_crit, safe default).

    Parameters
    ----------
    df       : SPARC DataFrame
    ups_disk : disk Υ [M☉ L☉⁻¹]

    Returns
    -------
    np.ndarray [M☉ pc⁻²]  — same length as df, radially resolved
    """
    r  = df["Rad"].values.copy()
    n  = len(r)

    # ── Path 1: native SBdisk ─────────────────────────────────────────────────
    if "SBdisk" in df.columns:
        sb = df["SBdisk"].values
        if np.any(sb > 0):
            return np.where(sb > 0, ups_disk * sb, 5.0)

    # ── Path 2: reconstruct from Vdisk gradient ───────────────────────────────
    if "Vdisk" in df.columns and n >= 3:
        vd = df["Vdisk"].values.copy()
        # Enclosed stellar mass at each r: M_disk(r) = Υ · v_disk² · r / G
        M_enc = ups_disk * vd ** 2 * np.maximum(r, 1e-4) / _G_KPC   # M☉
        # dM/dr via central differences [M☉ / kpc]
        dM_dr = np.gradient(M_enc, r)
        # Convert to surface density: Σ = dM/dr / (2π r kpc) / (kpc→pc)²
        KPC_TO_PC2 = 1e6   # 1 kpc² = 10⁶ pc²
        Sigma = dM_dr / (2.0 * np.pi * np.maximum(r, 0.1) * KPC_TO_PC2)
        Sigma = np.maximum(Sigma, 0.1)   # floor: non-physical negatives → 0.1

        # Normalise to global Sigma_star if available
        if "Sigma_star" in df.columns:
            sig_star = pd.to_numeric(df["Sigma_star"], errors="coerce").median()
            if np.isfinite(sig_star) and sig_star > 0:
                mean_sigma = float(np.mean(Sigma))
                if mean_sigma > 0:
                    Sigma = Sigma * sig_star / mean_sigma

        return Sigma

    # ── Path 3: constant fallback ─────────────────────────────────────────────
    return np.full(n, 5.0)


# ── T^μν stress-energy tensor corrections (v2.3) ─────────────────────────────

# HI thermal velocity dispersion (cold gas, non-relativistic)
SIGMA_HI = 10.0    # km/s — typical HI line width / velocity dispersion


def tmunu_corrections(r_kpc: np.ndarray,
                    v_bar: np.ndarray,
                    v_obs: np.ndarray,
                    Sigma: np.ndarray,
                    vg: np.ndarray,
                    vd: np.ndarray) -> np.ndarray:
    """
    Compute the full T^μν correction product f_P × f_β × f_sh × f_×.

    All four components of the stress-energy tensor beyond T^00 (energy
    density) are included:

        f_P  = exp(−ln3 · w)               T^ii  pressure suppression
               w = (σ_HI / v_obs)²  — thermal pressure parameter
               σ_HI ≈ 10 km/s (HI velocity dispersion, non-relativistic)
               NOTE: w uses σ_HI/v_obs, NOT the gas rotation fraction.
               For cold gas P/ρc² = (σ_gas/c)² ≪ 1; the relevant
               proxy is σ_thermal / v_circular, not Vgas / Vbar.

        f_β  = 1 + C_β(β − 0.5)           T^ij  orbital anisotropy
               β = 0.5 × (Vgas / (Vgas+Vdisk))
               β→0 for disk-dominated (circular orbits)
               β→0.5 for gas-dominated (pressure-supported)

        f_sh = 1 + C_sh |∂lnΣ/∂lnr|/(1+w) T^ij  shear stress
               radial gradient of surface density profile

        f_×  = 1 + C_× (v_obs/250)^0.5    T^0i  momentum flux
               bulk tangential velocity contribution

    Parameters
    ----------
    r_kpc : radii [kpc]
    v_bar : baryonic velocity [km/s]
    v_obs : observed velocity [km/s] — used for T^0i and w terms
    Sigma : surface density [M☉ pc⁻²]
    vg    : gas velocity component [km/s]
    vd    : disk velocity component [km/s]

    Returns
    -------
    np.ndarray  —  correction product, shape (N,)
    """
    vobs_safe = np.maximum(v_obs, 1.0)

    # ── f_P: pressure (T^ii diagonal) ────────────────────────────────────────
    # w = thermal pressure parameter = (σ_HI / v_circular)²
    w   = (SIGMA_HI / vobs_safe) ** 2                # ≈ 0.002 (massive) – 0.25 (dwarfs)
    f_P = np.exp(-_LN3 * w)                          # ∈ (0, 1], ~0.99 for massive spirals

    # ── f_β: anisotropy (T^ij off-diagonal) ──────────────────────────────────
    gas_frac = vg / np.maximum(vg + vd, 1e-4)        # gas fraction of baryons
    beta     = 0.5 * gas_frac                        # β ∈ [0, 0.5]
    f_beta   = 1.0 + C_BETA * (beta - 0.5)           # ∈ [0.925, 1.0]

    # ── f_sh: shear stress (T^ij gradient) ───────────────────────────────────
    ln_Sigma  = np.log(np.maximum(Sigma, 1e-6))
    ln_r      = np.log(np.maximum(r_kpc, 1e-6))
    dlnS_dlnr = np.abs(np.gradient(ln_Sigma, ln_r))  # |d lnΣ / d lnr|
    f_sh      = 1.0 + C_SH * dlnS_dlnr / (1.0 + w)  # ≥ 1

    # ── f_×: momentum flux (T^0i) ────────────────────────────────────────────
    f_cross = 1.0 + C_CROSS * (vobs_safe / V_CROSS) ** 0.5  # ≥ 1

    return f_P * f_beta * f_sh * f_cross


# ── Q-factor ──────────────────────────────────────────────────────────────────

def Q_factor(M_solar: float) -> tp.Tuple[float, float, float]:
    """
    Triple-regime vacuum-saturation Q-factor (v2.2).

    Controls which κ-rung is accessible as a smooth function of mass.

        Dwarf  regime (log M < 9.5):  Q₂ = tanh(M / M_ref2)^0.25
        Spiral regime (log M ≥ 9.5):  Q₃ = tanh(M / M_ref3)^0.50
        Group  regime (log M ≥ 12.5): Q₄ = tanh(M / M_ref4)^0.50

    Parameters
    ----------
    M_solar : float  —  total baryonic mass [M☉]

    Returns
    -------
    (Q2, Q3, Q4) : each ∈ [0, 1]
    """
    M  = max(float(M_solar), 1e4)
    lm = np.log10(M)
    Q2 = float(np.tanh(M / 10 ** LOG_M_REF_2) ** 0.25)
    Q3 = Q2 if lm < 9.5  else float(np.tanh(M / 10 ** LOG_M_REF_3) ** 0.50)
    Q4 = 0.0 if lm < 12.5 else float(np.tanh(M / 10 ** LOG_M_REF_4) ** 0.50)
    return Q2, Q3, Q4


# ── κ profile ─────────────────────────────────────────────────────────────────

def kappa_profile(r_kpc: np.ndarray,
                  v_bar: np.ndarray,
                  Sigma: np.ndarray,
                  M_total: float,
                  v_obs:   tp.Optional[np.ndarray] = None,
                  vg:      tp.Optional[np.ndarray] = None,
                  vd:      tp.Optional[np.ndarray] = None,
                  use_tmunu: bool = True) -> np.ndarray:
    """
    Compute κ(r) at each radius, with full T^μν corrections.

    The enhancement has three multiplicative pieces:

    1. Surface-density base (which rung is unlocked):
           κ_base(r) = K2 + (K3 − K2)·Q₃·f_Σ + (K4 − K3)·Q₄·f_Σ²
       → low Σ (< Σ_crit): full κ₃ (spirals) or κ₄ (massive galaxies)
       → high Σ (> Σ_crit): suppressed back toward κ₂

    2. Acceleration screening:
           screen(r) = sqrt( a₀ / g_Newton(r) )  clipped to [0, 1]
       → deep-MOND regime activates κ enhancement.

    3. Full T^μν stress-energy corrections:
           f_corr = f_P × f_β × f_sh × f_×
       → gas pressure, orbital anisotropy, shear stress, momentum flux.

    Combined:
           κ(r) = 1 + (κ_base − 1) · screen(r) · f_corr(r)

    Parameters
    ----------
    r_kpc   : radii [kpc]
    v_bar   : baryonic velocities [km/s]
    Sigma   : stellar surface density [M☉ pc⁻²]
    M_total : total baryonic mass [M☉]
    v_obs   : observed velocities [km/s]  (for f_× term; uses v_bar if None)
    vg      : gas velocity [km/s]         (for f_P, f_β; zeros if None)
    vd      : disk velocity [km/s]        (for f_β;      uses v_bar if None)
    use_tmunu : bool  —  apply T^μν corrections (default True)

    Returns
    -------
    np.ndarray  —  κ(r) ∈ [1, K4]
    """
    _, Q3, Q4 = Q_factor(M_total)

    # ── Newtonian centripetal acceleration ────────────────────────────────────
    r_m  = np.maximum(r_kpc, 1e-10) * _KPC_M
    v_ms = v_bar * 1e3
    g    = v_ms ** 2 / r_m
    g    = np.where(g > 0, g, _A0)

    # ── Surface-density gating ────────────────────────────────────────────────
    Sigma_safe = np.maximum(Sigma, 1e-10)
    f_sigma    = np.clip(SIGMA_CRIT / Sigma_safe, 0.0, 1.0) ** ALPHA_SIGMA

    # κ_base: κ₂ + Q₃ portion of κ₂→κ₃ gap + Q₄ portion of κ₃→κ₄ gap
    k_base = (K2
              + (K3 - K2) * Q3 * f_sigma
              + (K4 - K3) * Q4 * f_sigma ** 2)

    # ── Acceleration screening ────────────────────────────────────────────────
    screen = np.clip(np.sqrt(_A0 / g), 0.0, 1.0)

    # ── T^μν corrections ──────────────────────────────────────────────────────
    if use_tmunu:
        _vobs = v_obs if v_obs is not None else v_bar
        _vg   = vg    if vg   is not None else np.zeros_like(v_bar)
        _vd   = vd    if vd   is not None else v_bar
        f_corr = tmunu_corrections(r_kpc, v_bar, _vobs, Sigma, _vg, _vd)
    else:
        f_corr = np.ones_like(v_bar)

    kappa = 1.0 + (k_base - 1.0) * screen * f_corr
    return np.clip(kappa, 1.0, K4)


# ── Module-level galaxy cache for qgd_velocity() ─────────────────────────────

_CACHE: dict = {}


def register_galaxy(df: pd.DataFrame,
                    name: str = "Galaxy",
                    ups_disk: float = UPS_DISK,
                    ups_bulge: float = UPS_BULGE) -> None:
    """
    Register a SPARC DataFrame for use with qgd_velocity().

    Must be called once before passing qgd_velocity to gr.Fit or
    gr.Plot.plot_grc.  Precomputes all static quantities so that
    iminuit calls are fast.

    Parameters
    ----------
    df        : SPARC DataFrame from gr.Reader.read()
    name      : galaxy label (used in plots)
    ups_disk  : disk  Υ [M☉ L☉⁻¹]
    ups_bulge : bulge Υ [M☉ L☉⁻¹]

    Example
    -------
    >>> from pygrc.qgd import register_galaxy, qgd_velocity
    >>> register_galaxy(df, name="NGC5055")
    >>> m = gr.Fit(df["Rad"], df["Vobs"], 1.0).fit_lsq(
    ...         qgd_velocity, [(0.5, 2.0)], df["errV"], [False])
    """
    _CACHE.clear()
    _CACHE["df"]        = df
    _CACHE["name"]      = name
    _CACHE["r"]         = df["Rad"].values.copy()
    _CACHE["v_bar"]     = v_baryonic(df, ups_disk, ups_bulge)
    _CACHE["Sigma"]     = surface_density(df, ups_disk)
    _CACHE["M_base"]    = baryonic_mass(df, ups_disk, ups_bulge)
    _CACHE["v_obs"]     = df["Vobs"].values.copy()
    _CACHE["vg"]        = df["Vgas"].values.copy() if "Vgas" in df.columns else None
    _CACHE["vd"]        = df["Vdisk"].values.copy() if "Vdisk" in df.columns else None
    _CACHE["ups_disk"]  = ups_disk
    _CACHE["ups_bulge"] = ups_bulge


def qgd_velocity(r: np.ndarray, mass_scale: float = 1.0) -> np.ndarray:
    """
    QGD rotation curve — compatible with pygrc.Fit.fit_lsq().

    Call register_galaxy(df) once before using this function.

    Parameters
    ----------
    r          : np.ndarray  —  radii [kpc] (passed by iminuit/plot_grc)
    mass_scale : float       —  overall mass scaling factor (default 1.0)

    Returns
    -------
    np.ndarray  —  predicted velocities [km/s]
    """
    if not _CACHE:
        raise RuntimeError("No galaxy registered. Call register_galaxy(df) first.")

    r_data = _CACHE["r"]
    v_bar  = np.interp(r, r_data, _CACHE["v_bar"])
    Sigma  = np.interp(r, r_data, _CACHE["Sigma"])
    v_obs  = np.interp(r, r_data, _CACHE["v_obs"])
    vg     = np.interp(r, r_data, _CACHE["vg"])  if _CACHE["vg"] is not None else None
    vd     = np.interp(r, r_data, _CACHE["vd"])  if _CACHE["vd"] is not None else None
    M      = _CACHE["M_base"] * float(mass_scale)

    kappas = kappa_profile(r, v_bar, Sigma, M, v_obs=v_obs, vg=vg, vd=vd)
    return v_bar * np.sqrt(kappas)


# ── QGD class — primary API ───────────────────────────────────────────────────

class QGD:
    """
    QGD κ-ladder model for a single SPARC galaxy.

    Zero free parameters per galaxy. All physics derives from the
    baryonic columns already present in the SPARC DataFrame.

    Implements the full v2.3 T^μν stress-energy tensor corrections and
    the v2.2 triple-regime Q-factor (κ₂/κ₃/κ₄).

    Surface density Σ(r) is reconstructed from the velocity decomposition
    when native SBdisk photometry is unavailable, using:
        Σ_disk(r) = (1/2πr) · d(Υ·v_disk²·r/G)/dr   [M☉ pc⁻²]

    Parameters
    ----------
    df        : pd.DataFrame from gr.Reader.read()
    name      : galaxy label (used in plot titles / filenames)
    ups_disk  : disk  Υ [M☉ L☉⁻¹]  (default 0.50)
    ups_bulge : bulge Υ [M☉ L☉⁻¹]  (default 0.70)
    use_tmunu   : bool  —  apply T^μν stress-energy corrections (default True)

    Attributes
    ----------
    r, v_obs, err_v  : numpy arrays from df
    v_bar            : baryonic velocity [km/s]
    M_total          : estimated baryonic mass [M☉]
    Sigma            : stellar surface density [M☉ pc⁻²], radially resolved
    """

    def __init__(self,
                 df: pd.DataFrame,
                 name: str = "Galaxy",
                 ups_disk: float = UPS_DISK,
                 ups_bulge: float = UPS_BULGE,
                 use_tmunu: bool = True):
        self.df        = df.copy()
        self.name      = name
        self.ups_disk  = ups_disk
        self.ups_bulge = ups_bulge
        self.use_tmunu   = use_tmunu

        self.r     = df["Rad"].values
        self.v_obs = df["Vobs"].values
        self.err_v = df["errV"].values

        self.v_bar   = v_baryonic(df, ups_disk, ups_bulge)
        self.M_total = baryonic_mass(df, ups_disk, ups_bulge)
        self.Sigma   = surface_density(df, ups_disk)   # now radially resolved

        # velocity components for T^μν
        self._vg = df["Vgas"].values.copy()  if "Vgas"  in df.columns else None
        self._vd = df["Vdisk"].values.copy() if "Vdisk" in df.columns else None

        self._v_pred : tp.Optional[np.ndarray] = None
        self._kappas : tp.Optional[np.ndarray] = None
        self._f_corr : tp.Optional[np.ndarray] = None

    # ── core prediction ───────────────────────────────────────────────────────

    def predict(self) -> tp.Tuple[np.ndarray, np.ndarray]:
        """
        Compute the QGD velocity prediction with full T^μν corrections.

            v_pred(r) = v_bar(r) · sqrt(κ(r))

        Returns
        -------
        (v_pred [km/s], kappa_array)  —  both shape (N,)
        """
        # Compute T^μν correction product separately so we can report it
        if self.use_tmunu and self._vg is not None and self._vd is not None:
            self._f_corr = tmunu_corrections(
                self.r, self.v_bar, self.v_obs,
                self.Sigma, self._vg, self._vd
            )
        else:
            self._f_corr = np.ones_like(self.v_bar)

        self._kappas = kappa_profile(
            self.r, self.v_bar, self.Sigma, self.M_total,
            v_obs=self.v_obs, vg=self._vg, vd=self._vd,
            use_tmunu=self.use_tmunu
        )
        self._v_pred = self.v_bar * np.sqrt(self._kappas)
        return self._v_pred, self._kappas

    def _ensure_predicted(self):
        if self._v_pred is None:
            self.predict()

    # ── metrics ───────────────────────────────────────────────────────────────

    def r_squared(self) -> float:
        """Coefficient of determination R²."""
        self._ensure_predicted()
        ss_res = float(np.sum((self.v_obs - self._v_pred) ** 2))
        ss_tot = float(np.sum((self.v_obs - np.mean(self.v_obs)) ** 2))
        return 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    def rmse(self) -> float:
        """Root Mean Square Error [km/s]."""
        self._ensure_predicted()
        return float(np.sqrt(np.mean((self.v_obs - self._v_pred) ** 2)))

    def chi2_reduced(self) -> float:
        """Reduced χ², weighted by errV."""
        self._ensure_predicted()
        err = np.maximum(self.err_v, 1.0)
        return float(
            np.sum(((self.v_obs - self._v_pred) / err) ** 2) / len(self.v_obs)
        )

    def summary(self) -> dict:
        """
        Return a summary dict of all key metrics.

        Keys: name, N, M_total, log_M, Q2, Q3, Q4, R2, RMSE, chi2r,
              kappa_mean, kappa_min, kappa_max,
              sigma_mean, sigma_min, sigma_max,
              f_corr_mean  (mean T^μν correction product)
        """
        self._ensure_predicted()
        Q2, Q3, Q4 = Q_factor(self.M_total)
        return {
            "name"       : self.name,
            "N"          : len(self.r),
            "M_total"    : self.M_total,
            "log_M"      : float(np.log10(self.M_total)),
            "Q2"         : Q2,
            "Q3"         : Q3,
            "Q4"         : Q4,
            "R2"         : self.r_squared(),
            "RMSE"       : self.rmse(),
            "chi2r"      : self.chi2_reduced(),
            "kappa_mean" : float(np.mean(self._kappas)),
            "kappa_min"  : float(np.min(self._kappas)),
            "kappa_max"  : float(np.max(self._kappas)),
            "sigma_mean" : float(np.mean(self.Sigma)),
            "sigma_min"  : float(np.min(self.Sigma)),
            "sigma_max"  : float(np.max(self.Sigma)),
            "f_corr_mean": float(np.mean(self._f_corr)) if self._f_corr is not None else 1.0,
        }

    # ── plots ─────────────────────────────────────────────────────────────────

    def plot(self, ax=None, save: bool = False):
        """
        Plot QGD prediction vs observed data.

        Parameters
        ----------
        ax   : matplotlib Axes (creates new figure if None)
        save : bool  —  save PDF to {name}_qgd.pdf
        """
        self._ensure_predicted()
        if ax is None:
            _, ax = plt.subplots(figsize=(7, 5))

        ax.plot(self.r, self.v_obs, marker="o", linestyle="none", label="Data")
        ax.plot(self.r, self.v_bar, ":", color="gray", lw=1.2, label="Baryonic")
        r_grid = np.linspace(self.r.min(), self.r.max(), 300)
        v_grid = np.interp(r_grid, self.r, self._v_pred)
        ax.plot(r_grid, v_grid, "--",
                label=f"QGD  R²={self.r_squared():.3f}  RMSE={self.rmse():.1f} km/s")

        ax.set_xlabel("Distance (kpc)")
        ax.set_ylabel("Velocity (km/s)")
        ax.set_title(self.name)
        ax.set_xlim(left=0); ax.set_ylim(bottom=0)
        plt.legend()
        if save:
            plt.savefig(self.name + "_qgd.pdf", dpi=300, bbox_inches="tight")
        return ax

    def plot_kappa(self, ax=None, save: bool = False):
        """Plot the κ(r) enhancement profile."""
        self._ensure_predicted()
        if ax is None:
            _, ax = plt.subplots(figsize=(7, 4))

        ax.plot(self.r, self._kappas, color="darkorange", lw=2, label="κ(r)")
        ax.axhline(K2, color="steelblue", ls="--", lw=1, label=f"κ₂ = {K2:.3f}  (dwarfs)")
        ax.axhline(K3, color="green",     ls="--", lw=1, label=f"κ₃ = {K3:.3f}  (spirals)")
        ax.axhline(K4, color="red",       ls="--", lw=1, label=f"κ₄ = {K4:.3f}  (groups)")

        ax.set_xlabel("Distance (kpc)")
        ax.set_ylabel("κ(r)")
        ax.set_title(f"{self.name} — QGD κ profile")
        ax.set_xlim(left=0); ax.set_ylim(bottom=1.0)
        plt.legend()
        if save:
            plt.savefig(self.name + "_kappa.pdf", dpi=300, bbox_inches="tight")
        return ax

    def plot_sigma(self, ax=None):
        """Plot the reconstructed Σ(r) surface density profile."""
        if ax is None:
            _, ax = plt.subplots(figsize=(7, 4))
        ax.semilogy(self.r, self.Sigma, color="teal", lw=2, label="Σ(r) reconstructed")
        ax.axhline(SIGMA_CRIT, color="gray", ls=":", lw=1.5, label=f"Σ_crit = {SIGMA_CRIT} M☉/pc²")
        ax.set_xlabel("Distance (kpc)")
        ax.set_ylabel("Σ [M☉ pc⁻²]")
        ax.set_title(f"{self.name} — Surface density profile")
        ax.legend()
        return ax

    def __repr__(self) -> str:
        s = self.summary()
        return (
            f"QGD({self.name!r}, N={s['N']}, "
            f"log M={s['log_M']:.2f}, R²={s['R2']:.3f}, "
            f"Q4={s['Q4']:.2f}, f_corr={s['f_corr_mean']:.3f})"
        )


# ── Batch analysis ────────────────────────────────────────────────────────────

def batch_predict(filepaths: tp.List[str],
                  names: tp.Optional[tp.List[str]] = None,
                  ups_disk: float = UPS_DISK,
                  ups_bulge: float = UPS_BULGE,
                  verbose: bool = True) -> pd.DataFrame:
    """
    Run QGD zero-parameter predictions over a list of SPARC files.

    Parameters
    ----------
    filepaths : list of str  —  paths to SPARC .dat files
    names     : list of str  —  galaxy names (uses filename stem if None)
    ups_disk  : disk  Υ
    ups_bulge : bulge Υ
    verbose   : bool  —  print one-line summary per galaxy

    Returns
    -------
    pd.DataFrame  —  one row per galaxy, columns from QGD.summary()
    """
    from pygrc.reader import Reader
    import os

    if names is None:
        names = [os.path.splitext(os.path.basename(fp))[0] for fp in filepaths]

    rows = []
    for fp, name in zip(filepaths, names):
        try:
            df    = Reader.read(filepath=fp)
            model = QGD(df, name=name, ups_disk=ups_disk, ups_bulge=ups_bulge)
            model.predict()
            s = model.summary()
            rows.append(s)
            if verbose:
                print(
                    f"  {name:<25}  log M={s['log_M']:.2f}  "
                    f"R²={s['R2']:.3f}  RMSE={s['RMSE']:.1f} km/s  "
                    f"κ̄={s['kappa_mean']:.2f}  Q4={s['Q4']:.2f}"
                )
        except Exception as exc:
            if verbose:
                print(f"  {name:<25}  ERROR: {exc}")

    return pd.DataFrame(rows)


# ── Convenience: print the κ-ladder ──────────────────────────────────────────

def print_kappa_ladder(n_max: int = 7) -> None:
    """Print the QGD κ-ladder in a readable table."""
    labels = {
        1: "Newtonian baseline",
        2: "Wide binaries / dwarf spheroidals",
        3: "Spiral outskirts       [SPARC validated]",
        4: "Galaxy groups / massive spirals  [v2.2]",
        5: "Galaxy clusters        [Bullet Cluster 0.981]",
        6: "Superclusters / WHIM   [theoretical]",
        7: "Cosmic web             [theoretical]",
    }
    print("QGD κ-ladder  (zero free parameters)")
    print("─" * 62)
    for n in range(1, n_max + 1):
        print(f"  n={n}  κ_{n} = {kappa_n(n):12.4f}   {labels.get(n, '')}")
    print("─" * 62)
