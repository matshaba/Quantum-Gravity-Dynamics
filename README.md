# Quantum Gravitational Dynamics (QGD)

> **Core thesis:** Spacetime geometry is not fundamental. It emerges algebraically from a dimensionless phase field $\sigma_\mu$ derived from the semiclassical WKB limit of the Dirac spinor. All of classical General Relativity — and its quantum corrections — follow as consequences.

[![Theory](https://img.shields.io/badge/Theory-QGD-blue)](docs/THEORY.md)
[![Solutions](https://img.shields.io/badge/GR_Solutions-Algebraic-green)](solutions/)
[![Dark Matter](https://img.shields.io/badge/Dark_Matter-Zero_Free_Params-orange)](validation/)
[![License](https://img.shields.io/badge/License-MIT-lightgrey)](LICENSE)

---

## Table of Contents

1. [The Core Idea](#the-core-idea)
2. [The Fundamental Chain](#the-fundamental-chain)
3. [Master Equations](#master-equations)
4. [What QGD Solves](#what-qgd-solves)
5. [All GR Solutions Recovered Algebraically](#all-gr-solutions-recovered-algebraically)
6. [Dark Matter as Quantum Structure](#dark-matter-as-quantum-structure)
7. [Falsifiable Predictions](#falsifiable-predictions)
8. [Repository Structure](#repository-structure)
9. [Quick Start](#quick-start)

---

## The Core Idea

General Relativity takes the metric $g_{\mu\nu}$ as fundamental and derives dynamics from the Einstein-Hilbert action. QGD makes a different choice:

**The metric is not fundamental. It is an algebraic output of an underlying phase field.**

Starting from the Dirac spinor $\psi = R(x) e^{iS(x)/\hbar}$ in the semiclassical (WKB) limit, we identify the dimensionless phase gradient:

$$\sigma_\mu(x) \equiv \frac{1}{c} \partial_\mu S(x)$$

This single four-vector encodes all gravitational physics. The metric is then **constructed**, not solved for:

$$\boxed{g_{\mu\nu}(x) = T^\alpha_\mu T^\beta_\nu \left( M_{\alpha\beta} \circ \left[ \eta_{\alpha\beta} - \sum_{a=1}^{N} \varepsilon_a \, \sigma_\alpha^{(a)} \sigma_\beta^{(a)} - \kappa \ell_Q^2 \, \partial_\alpha \sigma^\gamma \partial_\beta \sigma_\gamma \right] \right)}$$

where:
- $T^\alpha_\mu$: coordinate transformation matrix (e.g., spherical, cosmological)
- $M_{\alpha\beta}$: geometric scaling matrix
- $\varepsilon_a \in \{+1, -1\}$: source signature (attractive vs repulsive)
- $\ell_Q = \sqrt{G\hbar^2/c^4}$: quantum gravitational length scale
- $\kappa \approx 2$: quantum stiffness coefficient

---

## The Fundamental Chain

```
Dirac spinor ψ = R(x)·exp(iS/ℏ)
        ↓  WKB limit
Phase field  σ_μ = (1/c)∂_μS
        ↓  algebraic construction
Metric  g_μν = η_μν - Σ εₐ σ_μ^(a) σ_ν^(a) - κℓ_Q² ∂σ∂σ
        ↓  variational principle
Field equation  □_g σ_μ = Q_μ + G_μ + T_μ + κℓ_Q²□_g²σ_μ
        ↓  classical limit ℓ_Q → 0
Einstein's equations  G_μν = (8πG/c⁴) T_μν
```

The relationship to GR is an **equivalence under variable substitution**, not a derivation. When $\sigma_\mu$ satisfies the QGD field equation and $g_{\mu\nu}$ is built from the master metric formula, Einstein's equations are automatically satisfied.

---

## Master Equations

### The Action

$$S_\sigma = \int d^4x \sqrt{-g(\sigma)} \left[ -\frac{c^4}{16\pi G} R[g(\sigma)] + \frac{1}{2}\nabla_\mu \sigma_\nu \nabla^\mu \sigma^\nu - \frac{\ell_Q^2}{2} \nabla_\alpha\nabla_\beta\sigma_\mu \nabla^\alpha\nabla^\beta\sigma^\mu \right] + S_{\text{matter}}$$

### The Field Equation

Variation of $S_\sigma$ with respect to $\sigma_\mu$ yields the fourth-order equation:

$$\boxed{\Box_g \sigma_\mu = Q_\mu(\sigma, \partial\sigma) + G_\mu(\sigma, \ell, H, q) + T_\mu + \kappa\ell_Q^2 \Box_g^2 \sigma_\mu + \mathcal{O}(\ell_Q^4)}$$

Sources:
- $Q_\mu$: nonlinear self-interactions (origin of dark matter phenomenology)
- $G_\mu$: coupling to Kerr-Schild and radiative sectors
- $T_\mu = \frac{1}{2}T^{\mu\nu}\sigma_\nu$: matter stress-energy
- $\kappa\ell_Q^2 \Box_g^2 \sigma_\mu$: quantum gravitational stiffness (resolves singularities)

### The Complete Wavefunction

The four-component gravitational wavefunction encoding all configurations:

$$\psi = \frac{2GMmi}{c\hbar} \left[ \psi_0 \begin{pmatrix}1\\0\\\sqrt{f}\\ ia\sin\theta\sqrt{g}\end{pmatrix} e^{-iS/\hbar} + \psi_1 \begin{pmatrix}0\\1\\-ia\sin\theta\sqrt{g}\\-\sqrt{f}\end{pmatrix} e^{-iS/\hbar} + \psi_2 \begin{pmatrix}\sqrt{f}\\ia\sin\theta\sqrt{g}\\1\\0\end{pmatrix} e^{+iS/\hbar} + \psi_3 \begin{pmatrix}-ia\sin\theta\sqrt{g}\\-\sqrt{f}\\0\\1\end{pmatrix} e^{+iS/\hbar} \right]$$

where the **universal gravitational scalar** $f(r,\theta)$ encodes all physics:

$$\boxed{f(r,\theta) = \underbrace{\frac{2GM}{c^2 r}}_{\text{mass}} - \underbrace{\frac{GQ^2}{c^4 r^2}}_{\text{charge}} + \underbrace{\frac{2Mr}{\Sigma}}_{\text{spin}} + \underbrace{\frac{\Lambda r^2}{3}}_{\Lambda} + \underbrace{\frac{b(r)}{r}}_{\text{pressure}} + \underbrace{H^2(t)r^2}_{\text{expansion}} - \underbrace{\int\frac{P(r)}{\rho(r)c^2}dr}_{\text{EOS}} + \underbrace{\kappa\frac{\hbar^2}{M^2c^2r^2}}_{\text{quantum}}}$$

Every known spacetime is a special case of $f(r,\theta)$.

---

## What QGD Solves

| Long-standing problem | GR status | QGD resolution |
|---|---|---|
| Gravitational energy localization | Pseudotensor only (109-year open problem) | True tensor $T^{\mu\nu}_\sigma$, positive definite |
| Black hole singularities | Generic, unavoidable | Resolved at $r \sim \lambda_C = \hbar/Mc$ |
| N-body exact solutions | No closed form | Exact algebraic: $\sigma_t = \sum_a \sqrt{2GM_a/c^2 r_a}$ |
| Dark matter | Requires new particles | Factorial $\kappa_j$ structure of quantum phase Taylor expansion |
| Quantum corrections | Undefined in GR | Explicit $\mathcal{O}(\hbar^2)$: $\delta g_{tt} = -G\hbar^2/(Mc^4 r^3)$ |
| Binary waveforms | Supercomputer, weeks | Algebraic, O(N²), seconds |
| Cosmological constant | Fine-tuning problem | $\rho_\sigma = 3H_0^2/8\pi G$, $w = -1$ from attractor solution |

---

## All GR Solutions Recovered Algebraically

The **source signature recipe** — choose $\sigma^{(a)}$ and $\varepsilon_a$, apply the master metric:

| Spacetime | $\sigma$-field | $\varepsilon_a$ | Physical effect |
|---|---|---|---|
| Schwarzschild | $\sqrt{2GM/c^2r}$ | $+1$ | Attractive mass |
| Kerr | $a\sin\theta\sqrt{2GM/c^2r}$ | $+1$ | Frame dragging |
| Reissner-Nordström | $\sqrt{GQ^2/c^4r^2}$ | $-1$ | EM repulsion |
| de Sitter ($\Lambda > 0$) | $Hr$ | $+1$ | Cosmological expansion |
| Anti-de Sitter | $\|H\|r$ | $-1$ | AdS geometry |

**Frame dragging as interference:** The Kerr off-diagonal term $g_{t\phi}$ is not input — it emerges automatically from the cross-product $\sigma_t^{(\text{mass})} \times \sigma_\phi^{(\text{spin})}$.

See [`solutions/`](solutions/) for fully worked algebraic constructions of each metric.

---

## Dark Matter as Quantum Structure

The Taylor expansion of the gravitational phase factor:

$$e^{i\phi} = e^{i\sigma_\mu x^\mu/\hbar} = \sum_{j=0}^{\infty} \frac{(i\sigma)^j}{j!}$$

generates factorial enhancement factors:

$$\kappa_j = \sqrt{\frac{(2j-1)!}{2^{2j-2}}} \qquad \Rightarrow \qquad \kappa = [1.00,\ 1.225,\ 2.74,\ 8.87,\ 37.7,\ 197,\ 1245, \ldots]$$

**Modified rotation curve (zero free parameters per galaxy):**

$$v^2(r) = \frac{GM_{\text{baryon}}(<r)}{r}\left[1 + \sum_{j=2}^{4} \kappa_j \, g_j(r)\right]$$

**Validation across 4,248 measurements (SPARC database, 175 galaxies):**
- $R^2 = 0.908$
- $\chi^2_\nu \approx 1.2$  
- RMS = 8.3 km/s
- **Zero free parameters per galaxy** (vs 5–7 for $\Lambda$CDM)

Cross-dataset universality: same $\kappa_j$ values for clusters, ellipticals, CMB peak spacing, and wide binary External Field Effect.

See [`validation/`](validation/) and [`validation/dark_matter.py`](validation/dark_matter.py).

---

## Falsifiable Predictions

All predictions follow from the theory with no additional assumptions:

| Prediction | Value | Testable with |
|---|---|---|
| Neutron star mass shift | $\sim 0.01 M_\odot$ from quantum TOV | NICER (current) |
| Binary merger separation | $d = 4r_s$ for equal masses | Next-gen GW detectors |
| Large-scale $\kappa_5$ activation | Specific correlation at 10–100 Mpc | DESI, Euclid |
| CMB higher peak modulation | $\ell_n \propto \kappa_{j(n)} \cdot f(n)$ | CMB-S4 |
| Maximum acceleration | $a_{\max} = 3mc^3/\hbar$ | — |
| Quantum perihelion shift | $\sim 10^{-90}$ arcsec/century | Unmeasurable |
| GW phase quantum shift | $\sim 10^{-73}$ rad | Unmeasurable |

The neutron star mass prediction with NICER is the **critical near-term falsification test**.

---

## Repository Structure

```
QGD/
├── DarkMatter/
│   ├── QGD.py
│   ├── dark_matter.py
│   ├── bullet_cluster.py
│   ├── FullStressEnergyTensor.py
│   ├── Uniqueness_of_kappa_values.py
│   ├── kappa_inversion.py
│   ├── darkmatter-theory.tex
│   ├── insights.md
│   ├── kappa-inversion-general.md
│   └── data/
├── core/
│   ├── graviton_field.py
│   ├── master_metric.py
│   ├── qgd_cosmology.py
│   ├── qgd_energy.py
│   ├── ringdown.py
│   ├── PN.py
│   ├── Effective_One_Body(EOB).py
│   ├── QGD_superposition.py
│   ├── QGD_vs_GR.py
│   ├── Rosseta_Stone.py
│   ├── two_and_three_body_solutions.py
│   ├── QGD_bug_fixes.py
│   └── nrpy/
│       └── InitialData_QGD.py
├── comparison/
│   ├── EFE_solutions_from_QGD_perspective.py
│   ├── comparison.md
│   └── comparison.py
├── docs/                          ← Theory
│   ├── Ch1-Foundations.tex + .pdf
│   ├── Ch2-Metric.tex + .pdf
│   ├── Ch3-FieldEquations.tex + .pdf
│   ├── Ch4-Energy.tex + .pdf
│   ├── Ch5-cosmology.tex + .pdf
│   ├── Ch6-ExactSolutions.tex + .pdf
│   ├── Ch7-Applications.tex + .pdf
│   ├── Ch8-GravitonApplications.tex + .pdf
│   ├── Ch9-QG_QFT.tex + .pdf
│   ├── Ch10-QGD_QG.tex + .pdf
│   ├── Ch11-DarkMatter.tex + .pdf
│   ├── Ch13-Radiation.tex + .pdf
│   ├── QuantumGravityDynamics.pdf  
│   ├── complete_paper.tex
│   ├── extended_summary.md
│   ├── summary.md
│   └── references.bib
├── notebooks/
│   ├── Energy_and_Cosmology.ipynb
│   ├── N_Body_solution
│   ├── dark_matter.ipynb
│   └── foundations.ipynb
├── ongoing-work/
│   ├── DarkEnergy.tex
│   ├── inflation.tex
│   ├── ch12.tex
│   ├── kinetic-term-derivation.tex
│   ├── quantum-term-derivation.tex
│   ├── Chq-corrections.md
│   └── EOB-correction.py
└── tests/
```

---

## Quick Start

```bash
git clone https://github.com/[author]/QGD
cd QGD
pip install numpy scipy matplotlib sympy
```

### Reconstruct Schwarzschild in 3 lines

```python
from core.sigma_field import SigmaField
from core.master_metric import MasterMetric

sigma = SigmaField.schwarzschild(M=1.0)        # σ_t = √(2GM/c²r)
g = MasterMetric.construct(sigma, coords='spherical')
print(g.line_element())
# ds² = -(1 - 2GM/c²r)dt² + (1 - 2GM/c²r)⁻¹dr² + r²dΩ²
```

### Compute a galaxy rotation curve

```python
from predictions.rotation_curves import QGDRotationCurve

galaxy = QGDRotationCurve(M_baryon=1e10)       # Solar masses
r, v_qgd, v_newton = galaxy.compute(r_max=50)  # kpc
galaxy.plot(show_kappa_contributions=True)
# Fits SPARC data with zero free parameters
```

### Generate a binary black hole waveform

```python
from predictions.gravitational_waves import BinaryWaveform

wf = BinaryWaveform(M1=36*M_sun, M2=29*M_sun, distance=410e6*pc)
t, h_plus, h_cross = wf.generate()
wf.plot_with_energy_decomposition()
# All spin-orbit, spin-spin terms emerge from σ cross products
```

---

## Key References

- Full theoretical derivation: [`docs/THEORY.md`](docs/complete_paper.tex)
- Original `.tex` manuscript: available on request
- SPARC rotation curve database: [SPARC](http://astroweb.cwru.edu/SPARC/)
- LIGO GW150914: Abbott et al. (2016), PRL 116, 061102

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
