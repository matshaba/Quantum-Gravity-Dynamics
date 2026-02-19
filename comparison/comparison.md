# QGD ↔ GR Translation Guide

> This document maps every major GR concept to its QGD equivalent.
> It is the "Rosetta Stone" for physicists familiar with GR.

---

## The Variable Substitution

The entire relationship between QGD and GR is captured by one substitution:

**GR uses:** $g_{\mu\nu}$ as the fundamental variable (10 independent components, symmetric tensor)

**QGD uses:** $\sigma_\mu$ as the fundamental variable (4 components, phase gradient vector)

**The bridge:**

$$g_{\mu\nu}^{\text{QGD}} = T^\alpha_\mu T^\beta_\nu\left(\eta_{\alpha\beta} - \sum_a \varepsilon_a \sigma_\alpha^{(a)}\sigma_\beta^{(a)} - \kappa\ell_Q^2\partial_\alpha\sigma^\gamma\partial_\beta\sigma_\gamma\right)$$

When $\sigma_\mu$ solves the QGD field equation, $g_{\mu\nu}^{\text{QGD}}$ satisfies Einstein's equations. This is an **equivalence**, not a derivation.

---

## Equation Dictionary

| GR concept | GR equation | QGD equivalent |
|---|---|---|
| Metric | $g_{\mu\nu}$ (solve EFE) | $g_{\mu\nu}[\sigma]$ (algebraic output) |
| Einstein tensor | $G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2}g_{\mu\nu}R$ | Satisfied automatically by $g[\sigma]$ |
| Field equations | $G_{\mu\nu} = \frac{8\pi G}{c^4}T_{\mu\nu}$ | $\Box_g\sigma_\mu = Q_\mu + G_\mu + T_\mu + \kappa\ell_Q^2\Box_g^2\sigma_\mu$ |
| Gravitational potential | $\Phi = -GM/r$ | $\sigma_t = \sqrt{-2\Phi/c^2} = \sqrt{2GM/c^2r}$ |
| Christoffel symbols | $\Gamma^\lambda_{\mu\nu} = \frac{1}{2}g^{\lambda\rho}(\partial g)$ | $\sim \sigma_\mu \partial_\nu\sigma_\lambda$ |
| Geodesic equation | $\ddot{x}^\mu + \Gamma^\mu_{\nu\lambda}\dot{x}^\nu\dot{x}^\lambda = 0$ | Phase gradient flow $\partial_\mu(\sigma_\nu x^\nu) = 0$ |
| Gravitational waves | Perturbation $h_{\mu\nu}$ of $g_{\mu\nu}$ | $h_{ij} \propto \sigma_i^{(1)}\sigma_j^{(2)}$ (cross-term) |
| Energy-momentum pseudotensor | $t^{\mu\nu}$ (not covariant) | $T^{\mu\nu}_\sigma = \frac{1}{2}(\partial^\mu\sigma_\lambda)(\partial^\nu\sigma^\lambda)$ (true tensor) |
| Event horizon | $g_{tt} = 0$ | $\Sigma = \sum_a \varepsilon_a(\sigma^{(a)})^2 = 1$ |
| Singularity | $r = 0$, $g_{\mu\nu} \to \infty$ | $r_{\min} \sim \lambda_C$ (quantum-resolved) |
| Dark matter | New particles | $\kappa_j$ factors from phase Taylor expansion |
| Dark energy | $\Lambda$ (fine-tuned) | $\rho_\sigma = 3H_0^2/8\pi G$ (attractor) |

---

## Solution-by-Solution Comparison

### Schwarzschild

**GR approach:** Solve $G_{\mu\nu} = 0$ with spherical symmetry → Birkhoff's theorem → unique solution.

**QGD approach:** Choose $\sigma_t = \sqrt{2GM/c^2r}$, $\varepsilon = +1$, apply master formula.

$$g_{tt}^{\text{GR}} = -\left(1 - \frac{2GM}{c^2r}\right) = g_{tt}^{\text{QGD}} = -(1 - \sigma_t^2) \qquad \checkmark$$

**QGD quantum correction** (no GR analog):
$$g_{tt}^{\text{QGD}} = -\left(1 - \frac{2GM}{c^2r} - \frac{G\hbar^2}{Mc^4r^3}\right)$$

Significant at $r \sim \lambda_C = \hbar/Mc$.

---

### Kerr

**GR approach:** Solve $G_{\mu\nu} = 0$ with axial symmetry + stationarity. Solution found by Kerr (1963) after decades of effort.

**QGD approach:** Add second σ-source with azimuthal component.

Frame dragging in GR: $g_{t\phi} = -2GMar\sin^2\theta/\Sigma$ — appears in the solved metric.

Frame dragging in QGD: **emerges automatically** from cross-product:
$$g_{t\phi} = T^t_t T^\phi_\phi \cdot (-\sigma_t \cdot \sigma_\phi) = -\frac{2Mar\sin^2\theta}{\Sigma}$$

**Same result. Different insight:** In QGD, frame dragging is manifestly quantum mechanical — it is the interference term between the time-like and space-like phase waves.

---

### N-body

**GR:** No exact N-body solution exists. Numerical relativity required for N ≥ 2.

**QGD:** Exact algebraic solution for all N:

$$\sigma_t^{\text{(total)}} = \sum_{a=1}^N \sqrt{\frac{2GM_a}{c^2 |x - x_a|}}$$

The metric:
$$g_{tt} = -(1 - \Sigma^2), \quad g_{ij} = \frac{\delta_{ij}}{1 - \Sigma^2}$$

**Why does linear superposition work?** Because σ-fields are phase fields — they add like electromagnetic potentials. The metric is nonlinear in σ (quadratic), but the field itself is not. This is why the problem simplifies so dramatically.

---

## The Counting Argument

GR has **10 independent metric components** satisfying 10 coupled nonlinear PDEs (EFE), minus 4 gauge freedoms = 6 physical degrees of freedom.

QGD has **4 σ-field components** satisfying 4 coupled wave equations, minus gauge freedoms.

The 4-component σ-field generates all 10 metric components via the master formula. This is an **under-determination at the field level** — the σ-field has more symmetry than the metric. This is the origin of QGD's computational advantage.

---

## Paradigm Comparison

```
GENERAL RELATIVITY
═══════════════════
Problem: Matter distribution T_μν given
Goal: Find metric g_μν
Method: Solve G_μν = 8πG/c⁴ · T_μν
         (10 coupled nonlinear PDEs, in general intractable)
Result: Spacetime geometry as curved manifold

QGD
═══
Problem: Matter distribution T_μν given
Goal: Find phase field σ_μ
Method: Solve □_g σ_μ = Q + G + T + κℓ_Q²□²σ
         (4 linear-in-σ wave equations)
Then: g_μν = TT(η - Σεσσ - κℓ_Q²∂σ∂σ)  [algebraic]
Result: Spacetime geometry as emergent from phase field

RELATIONSHIP
════════════
Same physical predictions (at classical level)
QGD = GR when ℓ_Q → 0 and field equation classical limit
QGD ≠ GR when quantum corrections active (r ~ λ_C)
```

---

## Why σ and Not Something Else?

The choice of $\sigma_\mu = (1/c)\partial_\mu S$ is not arbitrary. It is forced by:

1. **Dirac equation in WKB limit** — the natural semiclassical variable
2. **Probability current conservation** — gives $|\psi|^2 \sigma_\mu x^\mu = J/\hbar$
3. **Gauge invariance** — $\sigma_\mu$ is invariant under $S \to S + \text{const}$
4. **Newton's law recovery** — $\sigma_t = \sqrt{2GM/c^2r}$ directly gives $F = GMm/r^2$
5. **Dimensional consistency** — $\sigma_\mu$ is dimensionless, unlike $A_\mu$ in EM

The relationship $\sigma_t^2 = -g_{tt} - 1 = 2|\Phi|/c^2$ makes $\sigma_t$ the natural measure of gravitational "index of refraction" — the factor by which the phase velocity of quantum wavefunctions deviates from $c$.

---

## For QFT Practitioners

QGD can be rephrased in QFT language:

- $\sigma_\mu$ is analogous to the EM vector potential $A_\mu$
- The metric $g_{\mu\nu} = \eta_{\mu\nu} - \sigma_\mu\sigma_\nu + \ldots$ is analogous to $F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu$
- The QGD action has a Pais-Uhlenbeck structure: $(\Box - m_0^2)^{-1}$ propagator
- The graviton corresponds to the massless pole; the Planck-mass mode to the heavy pole (instantaneously damped)
- The theory is renormalizable by power-counting at $E \ll M_{\text{Pl}}$

The key difference from standard quantum gravity attempts: **don't quantize $g_{\mu\nu}$. Quantize $\sigma_\mu$.** The metric is derived, not fundamental — you can't "quantize geometry" because geometry is an emergent phenomenon in QGD.

---

*See `core/sigma_field.py`, `core/master_metric.py` for computational implementations.*
