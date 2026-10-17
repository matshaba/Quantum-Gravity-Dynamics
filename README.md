# Quantum Gravitational Dynamics (QGD)

> **Thesis (Romeo Matshaba):** spacetime geometry is not fundamental. It is built algebraically from a dimensionless phase field $\sigma_\mu$ obtained from the semiclassical (WKB) limit of the Dirac spinor, and General Relativity is recovered as a consequence.

[![Theory](https://img.shields.io/badge/Theory-QGD-blue)](docs/complete_paper.tex)
[![Solutions](https://img.shields.io/badge/GR_Solutions-Algebraic-green)](comparison/EFE_solutions_from_QGD_perspective.py)
[![Dark Matter](https://img.shields.io/badge/Dark_Matter-No_Per--Galaxy_Fit-orange)](validation/)
[![License](https://img.shields.io/badge/License-MIT-lightgrey)](LICENSE)

This repository holds the theory, the code that implements it, and an independent test suite that checks it. Every
claim below carries a status taken from [`docs/claims_ledger.md`](docs/claims_ledger.md):
**verified** (re-derived independently), **recorded** (true with a stated limit), **fails** (false as stated),
**untested** (needs input not in this repository).

---

## How QGD differs from string theory and loop quantum gravity

| | String theory | Loop quantum gravity | QGD (as proposed) |
|---|---|---|---|
| Starting point | 1-d strings (and branes) moving in a background; gravity appears as a massless spin-2 mode | Canonical quantisation of GR in Ashtekar variables (Rovelli 2004; Thiemann 2007) | A dimensionless phase field $\sigma_\mu$ from the WKB limit of the Dirac equation |
| Status of spacetime | Background-dependent in perturbative formulations; extra dimensions (10/11) | Emergent from spin networks; discrete area and volume spectra | Not fundamental: $g_{\mu\nu}$ is built algebraically from $\sigma_\mu$ (master metric) |
| Quantum gravity enters as | Strings, supersymmetry, compactification | Quantised geometry, a non-perturbative Hilbert space | A higher-derivative term $\kappa\ell_Q^2\,\Box^2\sigma$ in a classical-looking field equation |
| Dimensions and symmetry | Requires extra dimensions and supersymmetry | Four dimensions, no extra structure | Four dimensions, no extra structure |
| Where GR comes from | Low-energy limit | Semiclassical limit (still being worked out) | Algebraic recipes reproduce Schwarzschild and Reissner-Nordstrom exactly (verified); Kerr and others are not yet reproduced |
| Dark matter | Not part of the framework | Not part of the framework | A nonlinear self-interaction $Q_\mu$ of $\sigma$ (the $\kappa_j$ tower); fits SPARC curves but not with zero parameters |
| Singularities | Resolved by extended objects (model dependent) | Bounce / discreteness (e.g. loop cosmology) | Claimed removed; **fails** for the literal quantum term (below) |
| Black-hole entropy | Microstate counting (BPS cases) | Area-spectrum counting | A constant correction; the literal term gives $\pm\pi\kappa\ell_Q^2$, spin dependent on Kerr |

What is genuinely different is the economy of the proposal: no new dimensions, no quantised geometry, one field. That is also its main
exposure. String theory and LQG have worked-out quantum structure (spectra, amplitudes, state counting) that QGD does not yet have: it has
a classical field equation with a quantum-motivated term, and the tests in this repository show several of its stated links to GR and to
black-hole entropy fail as written (see *What the tests found*). The comparison is of programmes at different stages, not a verdict on
which is correct. Statements about the other two theories are summary-level and not tested here.

---

## Contents

0. [How QGD differs from string theory and LQG](#how-qgd-differs-from-string-theory-and-loop-quantum-gravity)
1. [The idea](#the-idea)
2. [Master equations](#master-equations)
3. [Recovering GR solutions](#recovering-gr-solutions)
4. [Dark matter](#dark-matter)
5. [Claims and status](#claims-and-status)
6. [What the tests found](#what-the-tests-found)
7. [Open decisions for the author](#open-decisions-for-the-author)
8. [Repository structure](#repository-structure)
9. [Quick start](#quick-start)
10. [Citation](#citation)

---

## The idea

GR takes the metric $g_{\mu\nu}$ as fundamental. QGD takes the Dirac spinor $\psi = R\,e^{iS/\hbar}$, defines the
dimensionless phase gradient

$$\sigma_\mu \equiv \frac{1}{c}\,\partial_\mu S ,$$

and constructs the metric from it:

$$g_{\mu\nu} = T^\alpha_\mu T^\beta_\nu\Big(M_{\alpha\beta}\circ\big[\eta_{\alpha\beta} - \sum_{a}\varepsilon_a\,\sigma^{(a)}_\alpha\sigma^{(a)}_\beta - \kappa\ell_Q^2\,\partial_\alpha\sigma^\gamma\partial_\beta\sigma_\gamma\big]\Big)$$

with $\varepsilon_a=\pm1$ the source signature (attractive or repulsive), $\ell_Q$ the quantum gravitational length and
$\kappa\approx 2$ the stiffness. (The sign convention of this formula versus the code is an open item, see
[decisions](#open-decisions-for-the-author).)

```
Dirac spinor  psi = R exp(iS/hbar)
      | WKB limit
phase field   sigma_mu = (1/c) d_mu S
      | algebraic construction
metric        g = eta - sum eps sigma sigma - kappa l_Q^2 (d sigma)(d sigma)
      | variational principle
field eq.     box_g sigma_mu = Q_mu + G_mu + T_mu + kappa l_Q^2 box_g^2 sigma_mu
      | l_Q -> 0
Einstein      G_mn = (8 pi G / c^4) T_mn
```

The relation to GR is an equivalence under variable substitution, not a derivation from GR's field equations.

---

## Master equations

**Action**

$$S_\sigma = \int d^4x\sqrt{-g(\sigma)}\Big[-\frac{c^4}{16\pi G}R[g(\sigma)] + \tfrac12\nabla_\mu\sigma_\nu\nabla^\mu\sigma^\nu - \tfrac{\ell_Q^2}{2}\nabla_\alpha\nabla_\beta\sigma_\mu\nabla^\alpha\nabla^\beta\sigma^\mu\Big] + S_{\rm matter}$$

**Field equation** (fourth order)

$$\Box_g\sigma_\mu = Q_\mu + G_\mu + T_\mu + \kappa\ell_Q^2\,\Box_g^2\sigma_\mu + \mathcal O(\ell_Q^4)$$

$Q_\mu$ nonlinear self-interaction (the dark-matter sector), $G_\mu$ Kerr-Schild and radiative couplings,
$T_\mu=\tfrac12T^{\mu\nu}\sigma_\nu$ matter.

**Universal scalar.** All static and stationary metrics are written through one function, $g_{tt}=-(1-f)$:

$$f(r,\theta)=\underbrace{\tfrac{2GM}{c^2r}}_{\rm mass}-\underbrace{\tfrac{GQ^2}{c^4r^2}}_{\rm charge}+\underbrace{\tfrac{2Mr}{\Sigma}}_{\rm spin}+\underbrace{\tfrac{\Lambda r^2}{3}}_{\Lambda}+\underbrace{\tfrac{b(r)}{r}}_{\rm pressure}+\underbrace{H^2(t)r^2}_{\rm expansion}-\underbrace{\int\tfrac{P}{\rho c^2}dr}_{\rm EOS}+\underbrace{\kappa\tfrac{\hbar^2}{M^2c^2r^2}}_{\rm quantum}$$

In the static form $G^t{}_t$ is affine in $f$, so every term has a definite implied energy density
$\rho=(rf)'/(8\pi G r^2)$ (see [what the tests found](#what-the-tests-found)).

---

## Recovering GR solutions

The source-signature recipe: pick $\sigma^{(a)}$ and $\varepsilon_a$, apply the master metric.

| Spacetime | $\sigma$-field | $\varepsilon$ | Status |
|---|---|---|---|
| Schwarzschild | $\sqrt{2GM/c^2r}$ | $+1$ | verified (Kerr-Schild form, Ricci-flat) |
| Reissner-Nordström | $\sqrt{GQ^2/c^4r^2}$ | $-1$ | verified (exactly the Maxwell stress) |
| de Sitter / Nariai | $Hr$ | $+1$ | verified (double horizon at $r=1/\sqrt\Lambda$, area $4\pi/\Lambda$, $GM=c^3/3\sqrt3H$) |
| Anti-de Sitter | $\lvert H\rvert r$ | $-1$ | not separately tested |
| Kerr | $a\sin\theta\sqrt{2GM/c^2r}$ | $+1$ | untested; $g_{t\phi}$ from the cross term is reproduced, $g_{rr},g_{\theta\theta}$ are not supplied by $\sigma$ |

Frame dragging appears as interference: the Kerr $g_{t\phi}$ comes from the product of the mass and spin $\sigma$'s.
[`comparison/EFE_solutions_from_QGD_perspective.py`](comparison/EFE_solutions_from_QGD_perspective.py) builds ten of
these and prints `Passed 10/10 checks`.

**Limit of the recipe.** It is exact when all terms share one null direction (mass + charge + $\Lambda$ on one centre:
for a Kerr-Schild metric the mixed Ricci tensor is linear in $h$). It is **not** exact for separated sources: the
$\sigma_t=\sum_a\sqrt{2GM_a/c^2r_a}$ N-body metric has a nonzero Ricci tensor of order $(GM)^2$ from the cross term.
So the "exact N-body solution" is a statement about the $\sigma$-field sum, not a vacuum solution of Einstein's equations.

---

## Dark matter

The Taylor expansion of the phase factor gives factorial enhancement constants

$$\kappa_j=\sqrt{\frac{(2j-1)!}{2^{2j-2}}}\ \Rightarrow\quad \kappa = [1.00,\ 1.225,\ 2.74,\ 8.87,\ 37.7,\ 197,\ 1233, \ldots]$$

and a modified rotation curve $v^2=\frac{GM_{\rm bar}(<r)}{r}\big[1+\sum_{j=2}^{4}\kappa_j g_j(r)\big]$.

The engine is [`validation/dark_matter.py`](validation/dark_matter.py). It has no per-galaxy fit parameters but does
have shared empirical constants (`SIGMA_CRIT`, `ALPHA`, `BETA_0`, `LOG_M_TRIGGER`, two mass-to-light ratios, `G_CRIT`),
so it is not parameter-free.

| SPARC quantity | Author's value | Reproduced here | Note |
|---|---|---|---|
| pooled $R^2$ | 0.908 | 0.91 | verified |
| RMS | 8.3 km/s | about 26 km/s | fails under the choices in `reproduce_sparc.py` |
| $\chi^2_\nu$ | 1.2 | about 150 | fails, same choices |
| data set | 4,248 points, 175 galaxies | 3,827 rows, 225 galaxies in `SPARC-rotational.txt` | different file |

For comparison, MOND with one fixed $a_0=1.2\times10^{-10}$ m/s² on the same inputs gives $R^2=0.93$, RMS 23 km/s
and a lower median per-galaxy RMS (10.8 vs 16.7 km/s). These statements concern this file and these choices (mass from
`V_bar(R_last)`, the engine's own surface-density estimate, no external field), not the theory; a different mass
estimate could change them. See `validation/sparc_check.png`.

**CMB peaks.** Planck binned TT gives peaks at $\ell=221,\ 535,\ 814,\ 1123,\ 1409$
(`validation/check_cmb_peaks.py`). The law $\ell_n=279.6\,n$ misses the first peak by 26% and the second and third by
4.5% and 3%, so "within 1%, $A$ fitted to the first peak" fails. A two-parameter line $296.3n-68.4$ is within 3%.

---

## Claims and status

| Claim | Status | Test |
|---|---|---|
| $\kappa_j$ formula, $\kappa_7=1233$, engine constants | verified | `tests/test_claims.py` |
| $\lvert\alpha_G\rvert^2=c\hbar/2GMm$, $=1/2$ at the Planck mass | verified | `tests/test_claims.py` |
| Schwarzschild, Reissner-Nordström, Nariai from the recipe | verified | `tests/test_claims.py`, `tests/test_gaps.py` |
| Exact N-body metric $\sum_a\sigma_a$ | **fails** (not vacuum) | `tests/test_gaps.py` |
| Quantum term resolves the singularity with positive energy | **fails** | `tests/test_gaps.py` |
| Merger at $d=4r_s$ | recorded (depends on which condition is called the merger; total $\sigma=1$ gives $8r_s$) | `tests/test_gaps.py` |
| Constant entropy correction $2\kappa\ell_Q^2/G\hbar$ | **fails** as stated: a constant does arise from the literal master metric, but it is $\pm\pi\kappa\ell_Q^2$ per horizon and not a Gauss-Bonnet/Wald effect (below) | `tests/test_nariai_entropy.py`, `tests/test_quantum_core.py` |
| Master-metric quantum term $=\kappa\hbar^2/(M^2c^2r^2)$ in $g_{tt}$ | **fails**: it acts on the angular part only | `tests/test_quantum_core.py` |
| SPARC RMS 8.3 km/s, $\chi^2_\nu\approx1.2$, zero free parameters | **fails** | `tests/test_sparc.py` |
| CMB peaks within 1% | **fails** | `tests/test_gaps.py` |
| $\ell_Q=\sqrt{G\hbar^2/c^4}$ is a length; $\rho_{\rm grav}=GM/4c^2r^3$ is an energy density | **fails** in SI | `tests/test_dimensions.py` |
| Neutron-star mass shift about $0.01M_\odot$; fourth-order equation reduces to Einstein; Kerr, Vaidya, interior from the recipe | untested | — |

Predictions listed by the author (none yet confirmed or excluded here): quantum-TOV neutron-star shift (NICER),
equal-mass merger separation, large-scale $\kappa_5$ activation (DESI, Euclid), CMB higher-peak modulation (CMB-S4),
maximum acceleration $a_{\max}=3mc^3/\hbar$; the quantum perihelion ($\sim10^{-90}$ arcsec/century) and GW phase
($\sim10^{-73}$ rad) shifts are unmeasurable. The neutron-star prediction is the main near-term test, and it is untested here.

---

## What the tests found

- **Singularity.** Each term of $f$ has an implied energy density. Mass is vacuum, charge positive, and the quantum
  term $+\kappa\hbar^2/(M^2c^2r^2)$ is the Reissner-Nordström term with $Q^2<0$: **negative** density. The Kretschmann
  scalar still diverges like $r^{-8}$. A regular core needs $f\sim r^2$ at small $r$ (a Hayward-type $f$). Whether the
  fourth-order $\ell_Q^2$ term supplies that is not decided by the static $f$.
- **Entropy** (`validation/nariai_entropy.py`, Planck units). The Euclidean action gives $S_{dS}=3\pi/\Lambda$
  ($S^4$, Euler characteristic 2) and $S_{\rm Nariai}=2\pi/\Lambda$ ($S^2\times S^2$, characteristic 4), matching the
  horizon areas. A truly topological constant is proportional to the Euler characteristic, i.e. to the number of
  horizons: it drops out of $dS/dM$ but not out of de Sitter→Nariai pair creation, where it multiplies the rate by
  $e^{C}$. The documented constant is not that, in either reading of the metric term $A/r^2$:
  - $A=\kappa\ell_Q^2$ constant: the area law gives $4\pi M^2+2\pi A$ (a constant, but $\pi$ times the documented one),
    puts a horizon at $r=\sqrt A$ even at $M=0$, and shifts the black-hole horizon by $+\pi A$ and the cosmological
    one by $-\pi A$, so it is not the same on every horizon.
  - $A=\kappa\hbar^2/(M^2c^2)$ as written in $f$: the shift is $2\pi\kappa/M^2$ and $dM=T\,dS$ fails.
  - A Gauss-Bonnet coupling $\alpha\,\mathrm{GB}$ gives $S=2\pi\alpha\chi$ (checked two ways) but is topological in four
    dimensions and cannot generate the $1/r^2$ term, so the entropy constant and the metric term are independent inputs.
    The entropy table and the metric term need one consistent definition of $A$.
- **Nariai with the quantum term.** The double-root point moves to $r_N=1/\sqrt\Lambda+A/2$, $GM_N=1/(3\sqrt\Lambda)-A/2$.

### The quantum term taken literally (`validation/qgd_quantum_core.py`, `validation/wald_entropy.py`)

The two readings above both start from the $A/r^2$ term of the universal $f$. The master-metric formula itself,
$g=\eta\mp\sum\varepsilon\sigma\sigma\mp\kappa\ell_Q^2\,\partial_a\sigma^c\partial_b\sigma_c$, gives something different. For a
radial null $\sigma_\mu=\sqrt{F(r)}\,k_\mu$ the tensor $\partial_a\sigma^c\partial_b\sigma_c$ equals $(F/r^2)P_{ij}$, with $P$ the
transverse projector: no $t$ or $r$ components. So $g_{tt}$ and $g_{tr}$ are unchanged and only the angular part moves,

$$ds^2=-(1-F)\,dv^2+2\,dv\,dr+\rho^2d\Omega^2,\qquad \rho^2=r^2+s\,\kappa\ell_Q^2\,F(r),\quad s=\pm1,$$

where the sign $s$ depends on whether the printed signature $(+,-,-,-)$ or the code's $g_{tt}=-(1-f)$ is followed.

- **Entropy.** Horizons ($F=1$) and surface gravity are unchanged, so $T=F'(r_h)/4\pi$ is unchanged, while the area is
  $4\pi(r_h^2+s\kappa\ell_Q^2)$. $S=A/4$ shifts by $s\,\pi\kappa\ell_Q^2$ on **every** horizon (black hole, cosmological,
  Nariai pair), independent of mass, and consistent with the first law (it is the constant of integration). This is the
  constant, mass-independent, horizon-counting correction the entropy table describes, but with coefficient
  $\pi\kappa\ell_Q^2$ (the table has $2\kappa\ell_Q^2$, a factor $\pi/2$) and a sign that depends on the signature reading.
- **It is not a Gauss-Bonnet or Wald effect.** The Wald entropy of the $\ell_Q^2(\nabla\nabla\sigma)^2$ action term is
  **exactly zero** for any $F(r)$ (Schwarzschild, de Sitter, Nariai): $\sigma$ is null and lies in the binormal plane of the
  horizon cross-section, so $R_{abcd}\sigma^a\sigma^b=0$. The same code gives $A/4$ for $R/16\pi$ and $4\pi\alpha$ for a
  Gauss-Bonnet term, and a nonzero value if $\sigma$ is tilted out of the plane, so the zero is structural. (Reading: $g$ is
  the dynamical field and $\sigma$ matter on it.)
- **Rotating case** (`validation/wald_kerr.py`). For Kerr, $\sigma=\sqrt{2Mr/\Sigma}\,(dv-a\sin^2\theta\,d\phi)$ has a $\phi$ component, so it is no longer confined to the binormal plane and the Wald term is **not** zero: $S_\sigma/(s\ell_Q^2)=0$ at $a=0$, $\simeq(\pi^2/3)(a/M)^2$ at small spin, 1.43 at $a/M=0.6$, 4.82 at $0.9$, 8.76 at $0.99$ (depends only on $a/M$; the same code gives $A/4=\pi(r_+^2+a^2)$ and $4\pi\alpha$ for Gauss-Bonnet at every spin). So the action term gives a **spin-dependent** entropy correction that vanishes for Schwarzschild, not a constant.
- **The literal quantum term on Kerr** (`validation/kerr_quantum_term.py`, `kerr_quantum_term.png`). The angular-only result carries over: pulled back to ingoing spheroidal coordinates, $\delta g=s\,\kappa\ell_Q^2\,(2Mr/\Sigma)(d\theta^2+\sin^2\theta\,d\phi^2)$ exactly (it reduces to $r_s/r\,d\Omega^2$ at $a=0$). But Kerr has $g^{r\phi}=a/\Sigma\neq0$, so an angular perturbation still moves the horizon: $\delta r_h=s\kappa\ell_Q^2\,a^2\sin^2\theta\,(r_+^2+a^2)/[2\Sigma^2(r_+-M)]$. The area shift is $\Delta S=\Delta A/4=\pi s\kappa\ell_Q^2/\sqrt{1-a^2/M^2}$ (numerical, matched to $10^{-9}$ at four spins; not derived), so the constant $\pm\pi\kappa\ell_Q^2$ is **not** universal for rotating sources: it grows with spin and diverges toward extremality, where first order in $\kappa\ell_Q^2$ stops being valid. At finite $\kappa\ell_Q^2$ the horizon angular velocity $g^{r\phi}/g^{rv}$ also varies with $\theta$ at first order (spread proportional to $\kappa\ell_Q^2$), so the perturbed metric, which is not a solution of anything, has no uniform rotating Killing horizon; $T$ and the first law are then undefined and $A/4$ is only an area.
- **With $\Lambda$** (Kerr-de Sitter, `validation/wald_kerr.py`, `validation/wald_kerr_ds.png`; the chart is an Einstein space
  $R_{ab}=\Lambda g_{ab}$ to rounding and regular at every root of $\Delta_r$). Assumption: $\sigma=\sqrt{2Mr/\Sigma}\,k$, the mass
  $\sigma$ (the rotating $\Lambda$ background is not flat plus a rank-one term, so the recipe gives no $\sigma_\Lambda$ there; the
  Wald term is quadratic in the amplitude). The spin dependence survives: at small spin
  $S_\sigma=-(4\pi^2/3)\,s\,\ell_Q^2\,a^2\,\sigma_M^2(r_h)F'(r_h)/r_h$ with $\sigma_M^2=r_s/r_h$ the mass-$\sigma$ amplitude at the horizon and
  $F'/2$ the surface gravity (a closed form fitted and verified to $<10^{-6}$ on both horizons, not derived), so it has **opposite signs on the
  black-hole and cosmological horizons** and **vanishes at the Nariai point**, where the leading term becomes
  $S_\sigma=(16\pi^2/15)\,s\,\ell_Q^2\,(a^2\Lambda)^2$ on each merged horizon (the rotating Nariai point exists for $a\sqrt\Lambda\lesssim0.47$).
  Calibrations hold on both horizons (area/4 with the $\Xi$ factor; Gauss-Bonnet $4\pi\alpha$) and $\Lambda\to0$ reproduces Kerr.
- **Singularity.** Not removed. For $s=+1$, $K\,r^6\to 51r_s^2/4$ as $r\to0$ (Schwarzschild: $12r_s^2$), independent of
  $\kappa\ell_Q^2$; for $s=-1$ the sphere area vanishes at $r_*=(\kappa\ell_Q^2r_s)^{1/3}$ and $K$ diverges there like
  $(r-r_*)^{-4}$.
- **Field equation.** With $\sigma=\sqrt{r_s/r}\,dv$ on Schwarzschild, $\Box_g\sigma_v=\sqrt{r_s}(r_s-r)/4r^{7/2}$ and
  $\Box_g\sigma_r=-2\sqrt{r_s}/r^{5/2}$, while the printed $Q_\mu$ gives $Q_v=2r_s^{3/2}/r^{7/2}$, $Q_r=0$. So
  $\Box\sigma=Q$ does not hold with $G_\mu=0$; $G_\mu$ is not specified for this case, so the equation cannot be checked as
  printed. ($\sigma^2=r_s/r$, not $\sigma$, is the harmonic function.)

### Exploratory: the classical limit as a high-entropy limit (not a claim of the theory)

`validation/classical_limit_entropy.py`, `classical_limit_entropy.png`, `tests/test_classical_limit_entropy.py`. Carcassi, Landini and Aidala
([arXiv:2411.00972](https://arxiv.org/abs/2411.00972)) argue that classical mechanics is the limit $S\gg0$ of quantum mechanics and that
$\hbar\to0$ is a formal stand-in for it. Set against QGD's "$\lim_{\hbar\to0}g^{\rm QGD}=g^{\rm GR}$":

- **Same order.** For Gaussian states ($x=\Sigma/\hbar$) the von Neumann entropy exceeds the classical one by
  $S_Q-S_C=-1/(24x^2)-1/(320x^4)-\dots$ (exact series). The leading correction is $O(\hbar^2)$, the order of every QGD correction, and the paper's stretch map
  ($\hbar\to\hbar/\lambda$) divides it by $\lambda^2$.
- **What the limit parameter is.** Taking the action scale of a mass $M$ localised within $r$ as $\Sigma\sim Mcr$ (an assumption), $x=r/\lambda_C$. At the
  Schwarzschild radius $r_s/\lambda_C=S_{BH}/2\pi$ exactly, so QGD's expansion parameter at the horizon is
  $(\lambda_C/r_s)^2=(2\pi/S_{BH})^2$. "Classical" then means $S_{BH}\gg1$, not $\hbar=0$; under the stretch map $S_{BH}\to\lambda S_{BH}$ and the parameter
  falls as $\lambda^{-2}$, consistently.
- **Where it stops.** The classical Gaussian entropy is negative for $x<1/e$ ($r<\lambda_C/e$); the quantum one vanishes at $x=1/2$. The paper calls such
  states fictitious. The region near $r=0$ is therefore outside any classical description, GR included, which frames (but does not settle) the singularity
  question: the literal quantum term leaves $K\to51r_s^2/(4r^6)$ there (above).
- **Is $\Sigma\sim Mcr$ derivable from $\sigma$?** (`validation/sigma_phase_gradient.py`, `tests/test_sigma_phase_gradient.py`.) Not as stated. QGD defines
  $\sigma_\mu=\partial_\mu S/(mc)$, but a gradient is closed and the recipe is not: in ingoing Kerr-Schild form $\sigma=\sqrt{F(r)}\,dv$, so
  $d\sigma=(\sqrt F)'\,dr\wedge dv=-\sqrt{r_s}/(2r^{3/2})\,dr\wedge dv\neq0$ for any non-constant $F$, and the Kerr recipe is not closed either (checked at a point).
  The exact Hamilton-Jacobi gradient of a particle falling from rest, $u=-dv+dr/(1+\sqrt{r_s/r})$, is closed with norm $-1$, whereas the recipe is null. So
  $\sigma$ is not the phase gradient of any single-valued $S$; the geometric phase $\oint\sigma\,dx\neq0$ quoted in `insights-of-QGD.md` needs exactly this nonzero
  curl, so that claim and "$\sigma=\partial S/(mc)$" cannot both hold. If one nevertheless reads the radial part as a phase, the scale is
  $x_\sigma=(mc/\hbar)\int_0^r\!\sqrt{r_s/r'}\,dr'=2mc\sqrt{r_sr}/\hbar=(m/M)(S_{BH}/\pi)\sqrt{r/r_s}$: it is set by the **test** mass $m$, equals
  $2\times$ the earlier value at $r=r_s$, and agrees with $\Sigma\sim Mcr$ only for $m=M$ (a self-field reading, which would have to be postulated).
- **Not shown.** That QGD's field equations have such a limit, that $\Sigma\sim Mcr$ is the right scale, or that the ledger's failing items (the $g_{tt}$ form,
  the dimensions of $\ell_Q$) are cured by it. The entropy view fixes the order and scale of the corrections, not their form.

### Exploratory: bounce and topology (not claims of the theory)

- **Bounce** (`validation/bounce_toy.py`). The repository has no covariant FRW form of $\sigma$, so a Friedmann equation
  cannot be derived from the master metric; the tests are conditional. The static term $A/r^2$ implies
  $\rho=-A/(8\pi Gr^4)$. If it acted like a fluid $\rho_q=-Ca^{-n}$ against matter of equation of state $w$, then $H=0$
  exists iff $n>3(1+w)$ and $\ddot a>0$ there automatically. With $n=4$: dust bounces (at $a_b=C/\rho_m$), radiation and stiffer
  matter do not. A term $-\kappa\ell^2H^4$ gives no bounce for either sign. The step $r\to a$ is an assumption.
- **Dodecahedral topology** (`validation/dodecahedral_spectrum.py`). The Poincaré dodecahedral space is $S^3/I^*$ with $I^*$
  the binary icosahedral group (order 120, built from its quaternions). Character averaging reproduces Klein's Molien
  series $(1+t^{30})/((1-t^{12})(1-t^{20}))$ and gives mode multiplicity $(k+1)N_k$: the lowest nonzero mode is $k=12$
  (eigenvalue 168, against $k=1$, eigenvalue 3, for the full $S^3$), then $20, 24, 30, 32$. This is the origin of the
  low-$\ell$ power suppression; no $C_\ell$ is computed, and QGD makes no statement about global topology.

**Where the checks came from.** The CICY / pyCICY project keeps an exactness ledger and re-derives gravity-side
results by independent routes, including the Nariai horizon $dS_2\times S^2$. Used here: the ledger discipline, the
independent-route principle (the geometry is redone in sympy, not through the QGD code), and the Nariai geometry as the
check on the de Sitter recipe. Nothing in this repository imports CICY.

Run everything with `python -m pytest -q`; expected failures are `xfail` with the reason, so a clean run does not hide
them. Regenerate the ledger with `python validation/claims_ledger.py`.

---

## Open decisions for the author

Fixing these needs a decision about the theory or files not in this repository.

1. **Reconstructed files.** `MasterMetric` and `two_and_three_body_solutions` were missing and are rebuilt from how
   `core/QGD_vs_GR.py` uses them. The script passes 32 of 33 checks; the failure is its own `f_peak` line (GR ISCO
   frequency 67 Hz against a merger peak near 150 Hz). Notably, "merger at $d=4r_s$" is taken as the separation where the
   self part of $\sigma^2$ at the centre of mass equals 1. If you have the originals, they should replace these.
2. **Gauge of the "improved" metric.** The signature is fixed in this revision: `eta = diag(-1,1,1,1)` and the diagonal and improved metrics are $(-,+,+,+)$, consistent with $g_{tt}=-(1-f)$ (they had four negative entries). Still open: the "improved" $g_{rr}=(1-f/2)^2(1+f/2)^2$ combines the areal $g_{tt}$ with an isotropic-gauge $g_{rr}$, so it is not Ricci-flat for Schwarzschild and `verify_classical_limit` returns `match = False`. `master_metric.schwarzschild_gauge_metric(f)` (areal gauge, $g_{rr}=1/(1-f)$) is the consistent form and is vacuum / Reissner-Nordström / de Sitter as appropriate; which gauge the theory intends is the author's call.
3. **Units.** In SI, $\ell_Q=\sqrt{G\hbar^2/c^4}$ has dimension (mass·length³)$^{1/2}$ (the Planck length is
   $\sqrt{G\hbar/c^3}$), and $\rho_{\rm grav}=GM/4c^2r^3$ has dimension 1/length². Both may be natural-unit shorthand.
4. **Which quantum term.** The universal $f$ and `master_metric.py::f9` put $\kappa\hbar^2/(M^2c^2r^2)$ into $g_{tt}$; the
   master-metric formula taken literally changes only the angular part (above), and the "improved" metric changes $g_{rr}$.
   These are three different corrections. The entropy table's $2\kappa\ell_Q^2$ matches none of them: the literal one gives
   $\pm\pi\kappa\ell_Q^2$ per horizon. The signature (and so the sign $s$) must also be fixed.
5. **Duplicates and placeholders.** `docs/Ch9-QG_QFT.tex` equals `docs/Ch8-GravitonApplications.tex` (and the PDFs);
   `notebooks/` holds four identical files; `DarkMatter/FullStressEnergyTensor.py` is LaTeX; `core/dark_matter.py` equals
   `DarkMatter/bullet_cluster.py`, as do the two `EFE_solutions_from_QGD_perspective.py`; the citation has placeholders.

Fixed in this revision: the metric signature in `core/master_metric.py` (above). Documentation fixed in this revision: $\kappa_7$ (1245→1233); Earth-electron $\lvert\alpha_G\rvert^2=1.46\times10^{39}$
→ $4.36\times10^{-11}$ (classical regime) with electron-electron $2.9\times10^{44}$ as the quantum example, in
`key-equations.md`, `GR-vs-QGD.md`, `docs/complete_paper.tex`, `docs/Ch1-Foundations.tex` and `docs/corrections.md`; broken
README links and the Quick Start; `EFE_solutions_from_QGD_perspective.py` crashing on import; a `SyntaxWarning` in
`core/PN.py`.

---

## Repository structure

```
QGD/
├── DarkMatter/
│   ├── QGD.py
│   ├── dark_matter.py
│   ├── bullet_cluster.py            (identical to core/dark_matter.py)
│   ├── FullStressEnergyTensor.py    (LaTeX source with a .py extension, not Python)
│   ├── Uniqueness_of_kappa_values.py
│   ├── kappa_inversion.py
│   ├── darkmatter-theory.tex, kappa_uniqueness.tex
│   ├── insights.md, kappa-inversion-general.md
│   └── data/                        (SPARC-rotational.txt, PLANCK-TT-Spectrum.txt)
├── core/
│   ├── graviton_field.py, master_metric.py, master-metric2.py
│   ├── qgd_cosmology.py, qgd_energy.py, waveform.py
│   ├── PN.py, PN-generator.py, post-newtonian.py, PostNewtonian/
│   ├── QGD_superposition.py, QGD_vs_GR.py, two_and_three_body_solutions.py (reconstructed)
│   ├── Rosseta_Stone.py, QGD_bug_fixes.py, dark_matter.py
│   └── nrpy/InitialData_QGD.py
├── comparison/
│   ├── EFE_solutions_from_QGD_perspective.py   (also at the repository root)
│   ├── comparison.md, comparison.py
├── validation/
│   ├── dark_matter.py               (QGD engine, v1.8)
│   ├── reproduce_sparc.py, sparc_check.png
│   ├── check_cmb_peaks.py, cmb_peaks_check.png
│   ├── nariai_entropy.py, nariai_entropy.png
│   ├── wald_entropy.py, wald_kerr.py, wald_kerr_ds.png, kerr_quantum_term.py, kerr_quantum_term.png, qgd_quantum_core.py, quantum_core.png
│   ├── bounce_toy.py, bounce_toy.png, classical_limit_entropy.py, classical_limit_entropy.png, sigma_phase_gradient.py, dodecahedral_spectrum.py, dodecahedral_spectrum.png
│   ├── claims_ledger.py             (writes docs/claims_ledger.md)
│   └── cmb_peaks/
├── docs/                            (theory, LaTeX + PDF)
│   ├── Ch1-Foundations ... Ch16-Quantization (.tex, most with .pdf)
│   ├── complete_paper.tex, corrections.md, claims_ledger.md, summary.md, extended_summary.md, references.bib
├── notebooks/                       (four files with identical content)
├── ongoing-work/                    (drafts: DarkEnergy, inflation, field equations, EOB-correction.py ...)
├── tests/                           (test_claims.py, test_gaps.py, test_dimensions.py, test_sparc.py, test_nariai_entropy.py, test_quantum_core.py, test_repo.py, binary_black_hole.py)
├── key-equations.md, GR-vs-QGD.md, comparison-to-GR-and-NG.md, insights-of-QGD.md, QGD_showcase.md
└── requirements.txt
```

---

## Quick Start

```bash
git clone https://github.com/brentharts/Quantum-Gravity-Dynamics.git
cd Quantum-Gravity-Dynamics
pip install -r requirements.txt pandas pytest
python -m pytest -q
```

### Reconstruct the GR solutions from σ-fields

```bash
python comparison/EFE_solutions_from_QGD_perspective.py      # prints "Passed 10/10 checks"
```

### Compute a galaxy rotation curve

```python
import sys
import numpy as np
sys.path.insert(0, "validation")
from dark_matter import predict_rotation_curve

r = np.linspace(1.0, 30.0, 30)                       # kpc
v_gas, v_disk, v_bulge = 20 + 0 * r, 120 * np.tanh(r / 4), 0 * r   # km/s (toy inputs)
v_qgd, kappa = predict_rotation_curve(r, v_gas, v_disk, v_bulge, total_mass=5e10)
print(np.round(v_qgd[::6], 1), np.round(kappa[::6], 2))
```

### Run the waveform module

```bash
python core/waveform.py
```

### Re-run the independent checks

```bash
python validation/reproduce_sparc.py      # SPARC comparison, writes validation/sparc_check.png
python validation/check_cmb_peaks.py      # CMB peak positions
python validation/nariai_entropy.py       # Euclidean vs area entropy, quantum-term shifts
python validation/wald_entropy.py         # Wald entropy calibration and the l_Q^2 term (zero for non-rotating)
python validation/wald_kerr.py            # the same term on Kerr and Kerr-de Sitter (spin dependent)
python validation/kerr_quantum_term.py    # literal quantum term on Kerr: horizon shift, entropy vs spin
python validation/qgd_quantum_core.py     # literal quantum term, small-r curvature, field-equation residual
python validation/classical_limit_entropy.py  # classical limit as S >> 0 (arXiv:2411.00972)
python validation/bounce_toy.py           # conditional bounce model
python validation/dodecahedral_spectrum.py  # S^3/I* spectrum
python validation/claims_ledger.py        # regenerate docs/claims_ledger.md
```

## Key References

- Full theoretical derivation: [`docs/complete_paper.tex`](docs/complete_paper.tex)
- Original `.tex` manuscript: available on request
- SPARC rotation curve database: [SPARC](http://astroweb.cwru.edu/SPARC/)
- LIGO GW150914: Abbott et al. (2016), PRL 116, 061102

**Cited in the author's paper**

1. C. Rovelli, *Quantum Gravity* (Cambridge University Press, 2004).
2. T. Thiemann, *Modern Canonical Quantum General Relativity* (Cambridge University Press, 2007).

**Suggested for exploration** (not cited by the paper; none is used by the code or tests here, and all five concern the foundations of
quantum mechanics rather than gravity, so any link to QGD is still to be established)

3. G. Carcassi, T. Thrien, C. A. Aidala, *Reverse Physics for Quantum Mechanics*, [arXiv:2608.27543](https://arxiv.org/abs/2608.27543).
4. G. Carcassi, C. A. Aidala, T. Thrien, *What non-additive integral for ensemble spaces?*, [arXiv:2602.01470](https://arxiv.org/abs/2602.01470).
5. G. Carcassi, R. Rozite, C. A. Aidala, *Response to "Are Hilbert Spaces Unphysical? Hardly, My Dear!"*, [arXiv:2507.09738](https://arxiv.org/abs/2507.09738).
6. G. Carcassi, M. Landini, C. A. Aidala, *Classical mechanics as the high-entropy limit of quantum mechanics*, [arXiv:2411.00972](https://arxiv.org/abs/2411.00972).
7. *A No-Go Theorem for ψ-ontic Models? Yes! Response to Criticisms* (reply to S. Gao's critique of the Carcassi-Oldofredi-Aidala theorem), [arXiv:2412.19182](https://arxiv.org/abs/2412.19182); author list not confirmed here.

---

## Citation

```bibtex
@misc{QGD2025,
  title  = {Quantum Gravitational Dynamics: Emergent Geometry from the Dirac Equation},
  author = {[Romeo Matshabba]},
  year   = {2026},
  note   = {GitHub: https://github.com/[matshaba]/Quantum-Gravity-Dynamics}
}
```
