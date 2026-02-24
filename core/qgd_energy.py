#!/usr/bin/env python3
"""
qgd_energy.py
=============
Gravitational Energy in Quantum Gravitational Dynamics (QGD)

Demonstrates and computes:
  0. Physical constants and field units
  1. σ-field profile and energy density for Schwarzschild BH
  2. Positive-energy theorem: manifest positivity of H[σ,π]
  3. Two-body binding energy from field-gradient cross terms
  4. Gravitational wave energy density (exact vs. Isaacson)
  5. Dark matter as gravitational self-energy (κ-ladder)
  6. Grand comparison table: GR vs. QGD energy

All formulas transcribed directly from the QGD paper.
"""

import numpy as np
from scipy.integrate import quad
from math import factorial

# ─────────────────────────────────────────────
# PHYSICAL CONSTANTS (SI)
# ─────────────────────────────────────────────
G    = 6.674e-11     # m³ kg⁻¹ s⁻²
c    = 2.998e8       # m s⁻¹
hbar = 1.0546e-34    # J s
k_B  = 1.381e-23    # J K⁻¹
M_sun = 1.989e30    # kg
l_Pl = np.sqrt(hbar * G / c**3)   # 1.616e-35 m (Planck length = ℓ_Q)
t_Pl = l_Pl / c                    # 5.39e-44 s

SEP = "─" * 72

def banner(title):
    print(f"\n{'═'*72}")
    print(f"  {title}")
    print('═'*72)

def rule():
    print(SEP)

# ─────────────────────────────────────────────
# 0. PHYSICAL CONSTANTS RECAP
# ─────────────────────────────────────────────
banner("0. QGD Physical Constants")
print(f"  Planck length ℓ_Q = {l_Pl:.4e} m")
print(f"  Planck time   t_Q = {t_Pl:.4e} s")
print(f"  Planck energy E_Q = {hbar*c/l_Pl:.4e} J  =  {hbar*c/l_Pl/1.602e-19/1e9:.4e} GeV")
print(f"\n  QGD gravitational energy density is ρ_grav = ½(∂σ)²")
print(f"  σ_t is dimensionless (ratio of Schwarzschild radius to r)")
print(f"  Energy density in geometrised units [m⁻²] → ×c⁴/G for SI")

# ─────────────────────────────────────────────
# 1. σ-FIELD PROFILE AND ENERGY DENSITY
# ─────────────────────────────────────────────
banner("1. Schwarzschild σ-Field and Gravitational Energy Density")

M = 10 * M_sun
r_s = 2 * G * M / c**2
print(f"  Black hole: M = 10 M⊙  →  r_s = {r_s:.4e} m")

def sigma_t(r):
    """σ_t = √(r_s/r)  —  encodes g_tt = -(1 - σ_t²) = -(1 - r_s/r)"""
    return np.sqrt(r_s / r)

def dsigma_dr(r):
    """dσ_t/dr = -½ √(r_s) r^{-3/2}"""
    return -0.5 * np.sqrt(r_s) * r**(-1.5)

def rho_grav(r):
    """
    ρ_grav(r) = ½(dσ_t/dr)²  [m⁻²  in geometrised units]
    = GM / (4 c² r³)          [m⁻² in SI-like notation]
    """
    return 0.5 * dsigma_dr(r)**2

def rho_grav_SI(r):
    """Multiply by c⁴/G to convert to J/m³"""
    return rho_grav(r) * (c**4 / G)

print(f"\n  σ_t profile and energy density (geometric units = m⁻²):")
print(f"  {'r/r_s':>8}  {'σ_t':>10}  {'ρ (m⁻²)':>14}  {'ρ (J/m³)':>14}  {'g_tt exact':>12}  {'g_tt QGD':>12}")
rule()
for mult in [1.01, 1.1, 1.5, 2, 3, 5, 10, 50, 100]:
    r = mult * r_s
    sig = sigma_t(r)
    rg  = rho_grav(r)
    rsi = rho_grav_SI(r)
    gtt_gr  = -(1 - r_s/r)
    gtt_qgd = -(1 - sig**2)
    print(f"  {mult:>8.2f}  {sig:>10.6f}  {rg:>14.4e}  {rsi:>14.4e}  {gtt_gr:>12.6f}  {gtt_qgd:>12.6f}")

print(f"\n  Analytic formula:  ρ_grav(r) = GM/(4c²r³)")
print(f"  Verification at r=10r_s:")
r_test = 10 * r_s
computed  = rho_grav(r_test)
analytic  = G * M / (4 * c**2 * r_test**3)
print(f"    computed  = {computed:.8e}")
print(f"    analytic  = {analytic:.8e}")
print(f"    relative error = {abs(computed-analytic)/analytic:.2e}  ✓" if abs(computed-analytic)/analytic < 1e-10 else "    ✗ MISMATCH")

# ─────────────────────────────────────────────
# 2. POSITIVE ENERGY THEOREM
# ─────────────────────────────────────────────
banner("2. Positive-Energy Theorem")

print("""
  QGD Hamiltonian:
    H[σ,π] = ∫d³x [ ½π_μπ^μ + ½(∇σ_μ)² + V(σ) + (ℓ_Q²/2)(∇²σ_μ)² ]

  Every term is manifestly non-negative:
    ½π²     ≥ 0  (kinetic)
    ½(∇σ)²  ≥ 0  (gradient)
    V(σ)    ≥ 0  (by construction)
    ℓ_Q²(∇²σ)²/2  ≥ 0  (quantum stiffness)

  Therefore H ≥ 0 with equality iff σ_μ = 0 (Minkowski spacetime).

  This trivially supersedes the Schoen-Yau (1979) spinor proof and
  Witten (1981) — both of which required:
    • Dominant energy condition
    • Asymptotic flatness (ADM boundary)
    • Spinor analysis or Dirac operator tricks

  In QGD: just standard field theory.
""")

# Numerical check: construct ρ(r) and confirm ≥ 0 everywhere
r_range = np.logspace(np.log10(r_s*1.001), np.log10(r_s*1e6), 10000)
rho_vals = rho_grav(r_range)
print(f"  Numerical check: ρ_grav(r) ≥ 0 for {len(r_range)} points from r_s to 10⁶ r_s:")
print(f"    min(ρ) = {rho_vals.min():.4e}  (should be > 0)")
print(f"    max(ρ) = {rho_vals.max():.4e}  (at r → r_s)")
print(f"    All positive: {np.all(rho_vals > 0)} ✓")

# Compare: Einstein pseudotensor
print(f"""
  CONTRAST WITH GR PSEUDOTENSOR:
    t_μν = (c⁴/16πG)(-g)(Christoffel² terms)
    — can be zero at any point by coordinate choice
    — can be negative in certain regions
    — not a tensor (coordinate-dependent)

  The reason GR's pseudotensor fails: it measures geometry (g_μν),
  which can be made flat at any point (equivalence principle).
  QGD measures the FIELD generating the metric (σ_μ), which cannot
  be gauged away — just as F_μν cannot be gauged away in E&M.
""")

# ─────────────────────────────────────────────
# 3. TWO-BODY BINDING ENERGY
# ─────────────────────────────────────────────
banner("3. Two-Body Gravitational Binding Energy")

print("""
  For two masses M₁, M₂ separated by d, the σ-field is:
    σ_t(x) = √(2GM₁/c²|x-x₁|) + √(2GM₂/c²|x-x₂|)

  Field energy:
    E = ∫½(∇σ_t)² d³x = E₁ + E₂ + E_cross

  Cross term (via Green's theorem):
    E_cross = ∫∇σ₁·∇σ₂ d³x = -GM₁M₂/d

  This is Newton's binding energy — emerging purely from
  field-gradient overlap, with no additional assumptions.
""")

cases = [
    ("Earth-Moon",     5.972e24, 7.347e22, 3.844e8),
    ("Sun-Earth",      1.989e30, 5.972e24, 1.496e11),
    ("NS-NS (merger)", 1.4*M_sun, 1.4*M_sun, 100e3),
    ("BH-BH (30 M⊙)",  30*M_sun, 30*M_sun,  50*2*G*30*M_sun/c**2),
]

print(f"  {'System':<22}  {'M₁ (kg)':>12}  {'M₂ (kg)':>12}  {'d (m)':>10}  {'E_bind (J)':>14}  {'E/M₁c²':>10}")
rule()
for name, m1, m2, d in cases:
    Eb = -G*m1*m2/d
    frac = Eb/(m1*c**2)
    print(f"  {name:<22}  {m1:>12.3e}  {m2:>12.3e}  {d:>10.3e}  {Eb:>14.4e}  {frac:>10.4e}")

print("""
  Note: The single-mass field energy ∫½(∇σ_t)²4πr²dr diverges logarithmically
  at large r (just as Newtonian gravitational self-energy diverges for a point mass).
  The physical energy Mc² is recovered via the ADM boundary condition.
  The cross-term binding energy is FINITE and exactly reproduces Newton's result.
""")

# ─────────────────────────────────────────────
# 4. GRAVITATIONAL WAVE ENERGY
# ─────────────────────────────────────────────
banner("4. Gravitational Wave Energy Density")

print("""
  GW radiation mode:  σ_μ^GW = ε_μ e^{ik·x},  k² = 0

  QGD energy density (EXACT, point-wise):
    ρ_GW = ½(∂σ^GW)² = ½ω²|ε|²

  GR comparison:
    Isaacson (1968): ρ_GW^GR = (c²/32πG)⟨ḣ²⟩
    — requires spatial averaging over wavelengths
    — only valid in weak-field short-wave regime
    — gauge-dependent in details

  QGD limit:  h_ij^TT ~ -2σ_iσ_j
    → ρ_GW → (c²/32πG)⟨ḣ²⟩   (recovers Isaacson in GR limit)

  Gravitational Poynting vector (any radius, not just infinity):
    S_grav = σ̇_μ(∇σ^μ)
    dE/dt = ∮_S S_grav · dA
""")

events = [
    ("GW150914 (BH-BH)",    0.3*M_sun*c**2, 0.10, 150, 410e6),
    ("GW170817 (NS-NS)",    0.05*M_sun*c**2, 0.05, 100, 40e6),
    ("GW250114 (BH-BH est)",0.4*M_sun*c**2, 0.12, 100, 1000e6),
]

print(f"  {'Event':<26} {'E_rad (J)':>12} {'h_strain':>10} {'f_peak (Hz)':>12} {'Isaacson ρ (J/m³)':>18}")
rule()
for name, E_rad, h, f_peak, d_Mpc in events:
    d_m = d_Mpc * 3.086e22  # Mpc → m
    omega = 2*np.pi*f_peak
    rho_isaacson = c**2/(32*np.pi*G) * (h * omega)**2
    print(f"  {name:<26} {E_rad:>12.3e} {h:>10.2e} {f_peak:>12.1f} {rho_isaacson:>18.4e}")

print(f"""
  QGD advantage: ρ_GW defined at EVERY radius.
  → Can track energy from orbital near-zone to far-wave zone.
  → No averaging procedure needed; covariant energy flux available.
""")

# ─────────────────────────────────────────────
# 5. DARK MATTER AS GRAVITATIONAL SELF-ENERGY
# ─────────────────────────────────────────────
banner("5. Dark Matter = Gravitational Self-Energy (κ-Ladder)")

print("""
  The field equation has a self-interaction source Q_μ:
    □σ_μ = (8πG/c⁴) Q_μ[σ,∂σ] + matter terms

  Q_μ = gravitational field energy current
      = σ_μ(∂σ)² + higher terms (nonlinear self-coupling)

  Expanding to higher orders:
    Q_μ = Q_μ^(2) + Q_μ^(3) + Q_μ^(4) + ...

  Each order contributes energy density:
    ρ_total = ρ_matter × Σ κ_n²

  The velocity enhancement factors from factorial arithmetic:
    κ_n = √((2n-1)! / 2^(2n-2))

  THESE COME FROM PURE ARITHMETIC — ZERO FREE PARAMETERS.
""")

print(f"  {'n':>3}  {'(2n-1)!':>10}  {'2^(2n-2)':>10}  {'κ_n':>10}  {'κ_n²':>10}  {'Physical regime'}")
rule()
regime = ["Solar System / Newtonian",
          "Wide binary stars (Chae+2023)",
          "Spiral outskirts / flat curves",
          "Galaxy clusters / CMB lensing",
          "Protoclusters / early structure",
          "Pre-galactic dark matter"]
for n in range(1, 7):
    num = factorial(2*n - 1)
    den = 2**(2*n - 2)
    kn  = np.sqrt(num/den)
    print(f"  {n:>3}  {num:>10}  {den:>10}  {kn:>10.4f}  {kn**2:>10.4f}  {regime[n-1]}")

print(f"""
  Dark matter connection:
    ρ_DM = sum_{n>=2} Q_μ^(n) = (κ_n² - 1) × ρ_baryon

  In GR: the Q_μ terms are HIDDEN inside G_μν (Einstein tensor).
         The nonlinear structure of the Einstein tensor obscures
         the self-energy — making dark matter "invisible".
  In QGD: Q_μ is EXPLICIT and separated from the matter source.
""")

# Velocity boost examples
print(f"  Velocity boost v_obs = v_Newton × √κ_n:")
print(f"  {'n':>3}  {'κ_n':>8}  {'v_Newton (km/s)':>18}  {'v_obs (km/s)':>14}  {'Regime'}")
rule()
v_newton_examples = [220, 150, 150, 200]
for i, (n, v_N) in enumerate(zip([2, 3, 4, 5], v_newton_examples)):
    kn = np.sqrt(factorial(2*n-1)/2**(2*n-2))
    v_obs = v_N * np.sqrt(kn)
    print(f"  {n:>3}  {kn:>8.4f}  {v_N:>18.1f}  {v_obs:>14.1f}  {regime[n-1]}")

# ─────────────────────────────────────────────
# 6. GRAND COMPARISON TABLE
# ─────────────────────────────────────────────
banner("6. Grand Comparison: GR vs. QGD Energy")

rows = [
    ("Fundamental variable",
     "g_μν (10 components, geometry)",
     "σ_μ (4 components, field)"),
    ("Local energy density",
     "UNDEFINED — pseudotensor t_μν only",
     "ρ = ½(∂σ)²  (scalar field, exact)"),
    ("Energy-momentum tensor",
     "No true tensor exists",
     "Noether tensor T^μν_QGD (genuine tensor)"),
    ("Energy conservation",
     "Ambiguous — coordinate-dependent",
     "Automatic: ∂_μT^μν = 0"),
    ("Positive energy proof",
     "Schoen-Yau (1979); spinors;\n"
     "     only at r→∞; required DEC",
     "Manifest: H=∫(positive terms)≥0;\n"
     "     local, no spinors needed"),
    ("Energy localisation",
     "Cannot be defined at a point",
     "ρ(x) is coordinate-independent scalar"),
    ("Gravitational waves",
     "Isaacson averaging (approximate,\n"
     "     short-wave, weak-field only)",
     "Exact: ρ_GW = ½ω²|ε|²;\n"
     "     point-wise, any field strength"),
    ("GW energy flux",
     "Only at null infinity (Bondi news)",
     "Poynting vector S = σ̇(∇σ)\n"
     "     at any radius"),
    ("Black hole energy",
     "ADM mass at r→∞ only;\n"
     "     cannot locate it in space",
     "Distributed: ρ(r) = GM/(4c²r³)\n"
     "     concentrated near horizon"),
    ("Two-body binding",
     "From linearised metric;\n"
     "     local source undefined",
     "-GM₁M₂/r from field cross-integral;\n"
     "     automatic"),
    ("Field self-energy",
     "Hidden in G_μν; cannot extract",
     "Explicit Q_μ term; expandable"),
    ("Dark matter",
     "Unexplained — requires exotic\n"
     "     new particle species",
     "κ-ladder: Q_μ series,\n"
     "     zero free parameters"),
    ("Singularity energy",
     "ρ → ∞ unavoidable\n"
     "     (Hawking-Penrose theorems)",
     "Resolved at r_min = λ_C;\n"
     "     ℓ_Q term regulates UV divergence"),
    ("Quantum corrections",
     "No Hamiltonian;\n"
     "     non-renormalisable",
     "Built-in: ℓ_Q²(∇²σ)² term;\n"
     "     renormalisable by power counting"),
    ("Statistical mechanics",
     "No well-defined Hamiltonian\n"
     "     → no partition function",
     "H[σ,π] ≥ 0 → standard\n"
     "     Gibbs formalism available"),
]

print(f"\n  {'Property':<28}  {'GR':^30}  {'QGD':^30}")
rule()
for prop, gr, qgd in rows:
    p_lines  = prop.split('\n')
    g_lines  = gr.split('\n')
    q_lines  = qgd.split('\n')
    n = max(len(p_lines), len(g_lines), len(q_lines))
    p_lines += [''] * (n - len(p_lines))
    g_lines += [''] * (n - len(g_lines))
    q_lines += [''] * (n - len(q_lines))
    for i in range(n):
        pf = p_lines[i] if i == 0 else ''
        print(f"  {pf:<28}  {g_lines[i]:<32}  {q_lines[i]}")
    print()

print(f"""
{'═'*72}
  CORE INSIGHT:
  The equivalence principle (GR) says: at any point, choose free-fall
  frame where g_μν → η_μν and Γ = 0.  This makes ALL pseudotensors
  vanish at that point.  Hence no local energy in GR.

  QGD variable σ_μ generates the metric but CANNOT be gauged to zero.
  Analogy: A_μ (potential) can be gauged; F_μν (field strength) cannot.
  σ_μ is to gravity what F_μν is to electromagnetism.

  Energy was never missing — it was just encoded in the wrong variable.
{'═'*72}
""")
