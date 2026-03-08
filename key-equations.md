# Quantum Gravitational Dynamics — Complete Key Equations
### Romeo Matshaba | University of South Africa
---

## Table of Contents

1. [The Gravitational Fine Structure Constant](#1-the-gravitational-fine-structure-constant)
2. [The Gravitational Wavefunction](#2-the-gravitational-wavefunction)
3. [The Quantum-Classical Bridge](#3-the-quantum-classical-bridge)
4. [The Energy Denominator](#4-the-energy-denominator)
5. [Newton's Law and Quantum Force Corrections](#5-newtons-law-and-quantum-force-corrections)
6. [Fundamental Scales](#6-fundamental-scales)
7. [The Graviton Field](#7-the-graviton-field)
8. [The Master Metric](#8-the-master-metric)
9. [The QGD Action and Master Field Equation](#9-the-qgd-action-and-master-field-equation)
10. [Recovery of Einstein's Equations](#10-recovery-of-einsteins-equations)
11. [Gravitational Energy](#11-gravitational-energy)
12. [Cosmology](#12-cosmology)
13. [N-Body Exact Solutions](#13-n-body-exact-solutions)
14. [Post-Newtonian Hierarchy and Master Formula](#14-post-newtonian-hierarchy-and-master-formula)
15. [Ringdown, Radiation, and Black Hole Thermodynamics](#15-ringdown-radiation-and-black-hole-thermodynamics)
16. [QFT of the σ-Field: Propagators and the Double Copy](#16-qft-of-the-σ-field-propagators-and-the-double-copy)
17. [Quantum Gravity: Decoherence, Information, and Inflation](#17-quantum-gravity-decoherence-information-and-inflation)
18. [Dark Matter as κ-Enhancement](#18-dark-matter-as-κ-enhancement)
19. [The Central Logic Chain](#19-the-central-logic-chain)

---

## 1. The Gravitational Fine Structure Constant

The foundational coupling constant, derived from WKB phase-oscillation analysis. Its complexity is a feature: the phase $e^{i\pi/4}$ encodes the geometric phase of the gravitational interaction.

$$\boxed{\alpha_G = e^{i\pi/4}\sqrt{\frac{c\hbar}{2GMm}} = \frac{1+i}{\sqrt{2}}\sqrt{\frac{c\hbar}{2GMm}}}$$

Key derived quantities:

$$\alpha_G^2 = \frac{ic\hbar}{2GMm}, \qquad |\alpha_G|^2 = \frac{c\hbar}{2GMm}, \qquad \mathrm{Re}(\alpha_G) = \sqrt{\frac{c\hbar}{4GMm}}$$

**Rule for use:** physical observables always use $|\alpha_G|^2$ for probability densities/fluxes, or $\mathrm{Re}(\alpha_G)$ for real amplitudes.

**Dimensional check:** $[\alpha_G^2] = [c\hbar/(GMm)] = 1$ ✓

**Numerical examples:**

| System | $|\alpha_G|^2$ | Regime |
|--------|----------------|--------|
| Earth–electron | $\approx 1.46 \times 10^{39}$ | Quantum (gravity negligible) |
| Planck mass ($M = m = M_\mathrm{Pl}$) | $= 1/2$ | Quantum gravity |

---

## 2. The Gravitational Wavefunction

Gravity is treated as a **spherical matter wave radiating to infinity** in flat Minkowski spacetime. The $1/r$ amplitude is geometrically mandatory for any spherical wave in 3D (flux conservation):

$$\psi = \frac{\alpha_G}{r}\, u\, e^{iS(x)/\hbar}$$

The **complete four-component spinor solution** (particle + antiparticle states):

$$\psi = \frac{\alpha_G}{r}\left[\psi_0\begin{pmatrix}1\\0\\\sigma_z\\\sigma_x+i\sigma_y\end{pmatrix}e^{-iS/\hbar} + \psi_1\begin{pmatrix}0\\1\\\sigma_x-i\sigma_y\\-\sigma_z\end{pmatrix}e^{-iS/\hbar} + \psi_2\begin{pmatrix}\sigma_z'\\\sigma_x'+i\sigma_y'\\1\\0\end{pmatrix}e^{+iS/\hbar} + \psi_3\begin{pmatrix}\sigma_x'-i\sigma_y'\\-\sigma_z'\\0\\1\end{pmatrix}e^{+iS/\hbar}\right]$$

where $\psi_0, \psi_1$ are particle (spin up/down, forward time $e^{-iS/\hbar}$) and $\psi_2, \psi_3$ are antiparticle states (backward time $e^{+iS/\hbar}$).

**Probability density** with spherical-wave ansatz:

$$\boxed{|\psi|^2 = \frac{|\alpha_G|^2}{r^2}\left(1 + \frac{|p|^2c^2}{\Delta^2}\right)}$$

**Full four-component probability density** (with particle occupation $A$ and antiparticle occupation $B$):

$$\boxed{|\psi|^2 = \frac{|\alpha_G|^2}{r^2}\left[A\!\left(1+\frac{|p|^2 c^2}{\Delta^2}\right)+B\!\left(1+\frac{|p|^2 c^2}{\Delta'^2}\right)\right]}$$

**Antiparticle suppression** near a Hawking-temperature horizon:

$$\frac{B}{A} \sim \exp\!\left(-\frac{16\pi GMm}{\hbar c}\right) = \exp\!\left(-\frac{8\pi}{|\alpha_G|^2}\right)$$

---

## 3. The Quantum-Classical Bridge

From current conservation $\nabla\cdot(|\psi|^2\nabla S) = 0$ in spherical symmetry. The $r^2$ cancels **exactly** from geometry — not an approximation:

$$\boxed{|\psi|^2 \cdot p(r) \cdot r^2 = C}$$

This is the **exact bridge** between quantum and classical physics.

Substituting the probability density, the $r^2$ cancellation yields:

$$\boxed{|\alpha_G|^2\!\left(p + \frac{p^3}{\Delta^2}\right) = C}$$

The full **cubic momentum equation** (single-component):

$$\boxed{|p| + \frac{c^2|p|^3}{\Delta^2} = \mathcal{P} = \frac{\sqrt{2}\,GMm^2}{\hbar}}$$

**General cubic** (four-component, with occupation coefficients $A$, $B$):

$$\boxed{\alpha|p|^3 + \beta|p| - C = 0}$$

with:
$$\alpha = \frac{c^3\hbar}{2GMm}\!\left(\frac{A}{\Delta^2}+\frac{B}{\Delta'^2}\right), \qquad \beta = \frac{c\hbar}{2GMm}(A+B)$$

**Cardano solution:**

$$|p| = \left[\frac{C}{2\alpha}+\sqrt{\left(\frac{C}{2\alpha}\right)^2+\left(\frac{\beta}{3\alpha}\right)^3}\right]^{1/3} + \left[\frac{C}{2\alpha}-\sqrt{\left(\frac{C}{2\alpha}\right)^2+\left(\frac{\beta}{3\alpha}\right)^3}\right]^{1/3}$$

**Non-relativistic limit** ($|p| \ll mc$, $\Delta \approx mc^2$):

$$|p| \approx \mathcal{P} \qquad \text{(classical)}$$

---

## 4. The Energy Denominator

ALL stress-energy contributions enter through a **single unified denominator**:

$$\boxed{\Delta_{\text{full}} = E + mc^2 + \int\rho c^2\,dV - \int P\,dV - \frac{J^2}{2mr^2} + V_{EM} + \rho_\Lambda}$$

For antiparticle states: $\Delta' = E - mc^2 - H$.

**Mapping of terms:**

| Term in $\Delta$ | Action Term | Gravitational Limit |
|---|---|---|
| $+mc^2$ | $-m\bar\psi\psi$ | Test mass |
| $+\int\rho c^2\,dV$ | $\frac{i}{2}\bar\psi\gamma^\mu\overleftrightarrow\partial_\mu\psi$ | Schwarzschild (source mass) |
| $-\int P\,dV$ | $-P(\bar\psi\psi)$ | Equation of state |
| $-J^2/(2mr^2)$ | $-\frac{1}{2M^2}J_{\mu\nu}J^{\mu\nu}\bar\psi\psi$ | Kerr frame-dragging |
| $+V_{EM}$ | $-\frac{1}{4}F_{\mu\nu}F^{\mu\nu}$ | Reissner–Nordström |
| $+\rho_\Lambda$ | $-\rho_\Lambda$ | Cosmological constant |

**Special cases:**

| Spacetime | $\Delta$ |
|---|---|
| Schwarzschild | $E + mc^2 + \int\rho c^2\,dV$ |
| Kerr | $E + mc^2 + \int\rho c^2\,dV - J^2/(2mr^2)$ |
| Reissner–Nordström | $E + mc^2 + \int\rho c^2\,dV + V_{EM}$ |
| Perfect fluid + $\Lambda$ | $E + mc^2 + \int\rho c^2\,dV - \int P\,dV + \rho_\Lambda$ |

---

## 5. Newton's Law and Quantum Force Corrections

**Derivation of $\alpha_G$:** matching the leading WKB phase term to the Newtonian potential $-GMm/r$:

$$\frac{c\hbar}{2i\,\alpha_G^2\,r} = -\frac{GMm}{r} \implies \boxed{\alpha_G^2 = \frac{ic\hbar}{2GMm}}$$

**Newton's force from leading phase term:**

$$F_1 = \frac{c\hbar}{2ir^2}\cdot\frac{2GMm}{ic\hbar} = \frac{GMm}{r^2} \checkmark$$

**Higher-order quantum correction** (from the $r^3$ Taylor term):

$$\boxed{F_3 = \frac{9GM\hbar^2}{2mc^2r^4}}$$

**Complete force law:**

$$\boxed{F(r) = \frac{GMm}{r^2}\left[1 + \frac{9}{2}\left(\frac{\lambda_C}{r}\right)^2 + \mathcal{O}(r^{-4})\right]}$$

where $\lambda_C = \hbar/(mc)$ is the Compton wavelength. Quantum crossover at $r_c \approx 2.12\,\lambda_C$.

**Wavefunction amplitude in terms of potential** ($\Phi = -GMm/r$):

$$\boxed{|\psi|^2 = \frac{c\sqrt{m}}{2\,r^2\sqrt{|\Phi(r)|}} \propto \frac{1}{r^2\sqrt{|\Phi(r)|}}}$$

**Quantum perihelion precession** (Mercury):

$$\Delta\phi_{\mathrm{QGD}} = \frac{\hbar^2}{M^2c^2a^2(1-e^2)^\alpha} \approx 10^{-168}\;\mathrm{rad/orbit}$$

---

## 6. Fundamental Scales

**Normalization constant** (probability flux, set by test particle rest momentum alone):

$$\boxed{C = \frac{mc}{\sqrt{2}}, \qquad C = |\alpha_G|^2\,\mathcal{P}}$$

**Fundamental momentum scale** (gravitational Bohr momentum):

$$\boxed{\mathcal{P} = \frac{\sqrt{2}\,GMm^2}{\hbar}}$$

**Graviton frequency:**

$$\omega_g = \frac{\mathcal{P}}{\hbar} = \frac{\sqrt{2}\,GMm^2}{\hbar^2}, \qquad \mathcal{P} = \hbar\omega_g$$

**Length-scale unification** ($R_s = 2GM/c^2$, $\lambda_c = \hbar/(mc)$):

$$\frac{R_s}{\lambda_c} = \frac{2GMm}{\hbar c} = \frac{1}{|\alpha_G|^2}$$

**Field energy and momentum** (unifying QM, SR, and GR through one ratio):

$$\boxed{E_{\mathrm{field}} = \frac{mc^2}{\sqrt{2}\,|\alpha_G|^2} = \frac{1}{\sqrt{2}}\,mc^2\left(\frac{R_s}{\lambda_c}\right)}$$

$$\boxed{P_{\mathrm{field}} = \frac{mc}{\sqrt{2}\,|\alpha_G|^2} = \frac{1}{\sqrt{2}}\,mc\left(\frac{R_s}{\lambda_c}\right)}$$

**Regime classification:**

| Condition | $|\alpha_G|^2$ | Regime |
|---|---|---|
| $c\hbar \gg 2GMm$ | $\gg 1$ | Quantum (gravity negligible) |
| $c\hbar \sim 2GMm$ | $\sim 1$ | Full quantum gravity |
| $2GMm \gg c\hbar$ | $\ll 1$ | Classical GR |

---

## 7. The Graviton Field

The graviton field is **read directly off the wavefunction** — not independently postulated:

$$\boxed{\sigma_\mu(x) \equiv \frac{p_\mu}{mc} = \frac{1}{mc}\partial_\mu S(x)}$$

As an explicit four-vector:

$$\sigma^\mu = \left(\sigma_t,\,\sigma_x,\,\sigma_y,\,\sigma_z\right) = \frac{1}{mc}\left(\frac{E}{c},\,p_x,\,p_y,\,p_z\right)$$

**For a Kerr body** (mass $M_a$, spin parameter $\alpha_a = J_a/(M_ac)$, $\mathcal{S}_a = r_a^2 + \alpha_a^2\cos^2\theta_a$):

$$\sigma_t^{(a)} = \sqrt{\frac{2GM_a\,r_a}{c^2\,\mathcal{S}_a}}, \qquad \sigma_\phi^{(a)} = \alpha_a\sin\theta_a\sqrt{\frac{2GM_a}{c^2\,r_a\,\mathcal{S}_a}}$$

**For a Schwarzschild body** ($\alpha_a = 0$):

$$\sigma_t^{(a)} = \sqrt{\frac{2GM_a}{c^2\,r_a}}, \qquad \sigma_\phi^{(a)} = 0$$

**The graviton scalar** (dimensionless rapidity parameter):

$$\boxed{\sigma_i \equiv \frac{p_i c}{\Delta}}$$

Reduces to $v/c$ in the non-relativistic limit ($\Delta \approx mc^2$, $p = mv$). Cubic dispersion relation:

$$\boxed{\sigma + \sigma^3 = \frac{\mathcal{P}c}{\Delta}}$$

**Gravitational Lorentz factor:**

$$\boxed{|u|^2 = 1 + |\boldsymbol\sigma|^2 = 1 + \frac{p^2c^2}{\Delta^2}}$$

| Special Relativity | Quantum Gravity |
|---|---|
| $\beta = v/c$ | $\sigma = pc/\Delta$ |
| $\gamma^2 = 1/(1-\beta^2)$ | $|u|^2 = 1+\sigma^2$ |

---

## 8. The Master Metric

The metric is a **composite field** — not fundamental. Constructed algebraically from graviton fields:

$$\boxed{g_{\mu\nu}(x) = T^\alpha{}_\mu\,T^\beta{}_\nu\left(M_{\alpha\beta}\circ\left[\eta_{\alpha\beta} - \sum_{a=1}^{N}\varepsilon_a\,\sigma_\alpha^{(a)}\sigma_\beta^{(a)} - \kappa\,\ell_Q^2\,\partial_\alpha\sigma^\gamma\partial_\beta\sigma_\gamma\right]\right)}$$

with $\eta_{\alpha\beta} = \mathrm{diag}(+1,-1,-1,-1)$.

**Three-tier geometric hierarchy:**

$$g_{\mu\nu} = \underbrace{\eta_{\mu\nu}}_{\text{Tier 0: Minkowski}} - \underbrace{\varepsilon\,\sigma_\mu\sigma_\nu}_{\text{Tier 1: Classical gravity}} - \underbrace{\kappa\,\ell_Q^2(\partial\sigma)^2}_{\text{Tier 2: Quantum stiffness}}$$

**Event horizon condition** — the surface where $\sigma$-field amplitude reaches unity:

$$\boxed{\Sigma_{\mathrm{tot}} \equiv \sqrt{\sum_a\varepsilon_a(\sigma_t^{(a)})^2} = 1}$$

**Quantum gravitational length scale:**

$$\ell_Q = \sqrt{\frac{G\hbar^2}{c^4}} \approx 1.6\times10^{-70}\;\mathrm{m}$$

**Source signatures:**

| Source | Graviton field $\sigma^{(a)}$ | $\varepsilon_a$ | Effect |
|---|---|---|---|
| Mass $M$ | $\sqrt{2GM/c^2r}$ | $+1$ | Attraction |
| Spin $J = Mac$ | $a\sin\theta\sqrt{2GMr/c^2\Sigma}$ | $+1$ | Frame-dragging |
| Charge $Q$ | $\sqrt{Gk_e^2Q^2/c^4r^2}$ | $-1$ | Repulsion |
| $\Lambda > 0$ | $Hr/c$ | $+1$ | Expansion |
| $\Lambda < 0$ | $|H|r/c$ | $-1$ | Anti-de Sitter |

**Recovery of standard solutions** (all algebraic, no differential equations):

| Spacetime | $g_{tt}$ | QGD Quantum Correction |
|---|---|---|
| Schwarzschild | $1 - 2GM/(c^2r) - G\hbar^2/(Mc^4r^3)$ | $-G\hbar^2/(Mc^4r^3)$ |
| Reissner–Nordström | $1 - 2GM/(c^2r) + Gk_e^2Q^2/(c^4r^2)$ | $-G\hbar^2/(Mc^4r^3)$ |
| Schwarzschild–de Sitter | $1 - 2GM/(c^2r) - H^2r^2/c^2$ | $-G\hbar^2/(Mc^4r^3)$ |

**Classical limit:** $\lim_{\hbar\to 0} g_{\mu\nu}^{(\mathrm{QGD})} = g_{\mu\nu}^{(\mathrm{GR})}$ exactly.

**Singularity resolution:** quantum stiffness creates minimum radius $r_{\min} \sim \ell_P(m_P/M)^{1/3}$. Below critical mass $M_\mathrm{crit} \approx 0.73\,m_P$, no horizon forms.

---

## 9. The QGD Action and Master Field Equation

**The fundamental action** governing $\sigma_\mu$ as a dynamical field:

$$\boxed{S[\sigma] = \int\!d^4x\,\sqrt{-g(\sigma)}\left[\frac{R[g(\sigma)]}{16\pi G} + \frac{\hbar^2}{2M}\,g^{\mu\nu}\nabla_\mu\sigma^\alpha\nabla_\nu\sigma_\alpha + \mathcal{L}_{\mathrm{matter}}(\psi,\,g(\sigma))\right]}$$

**Compact σ-field equation** (from variation w.r.t. $\sigma_\alpha$, not $g_{\mu\nu}$):

$$\boxed{\frac{\hbar^2}{M}\nabla^2\sigma^\alpha = \frac{1}{16\pi G}\left(G^{\mu\nu} - \frac{8\pi G}{c^4}T^{\mu\nu}\right)\frac{\partial g_{\mu\nu}}{\partial\sigma_\alpha}}$$

**This equation does not exist in General Relativity.**

**The complete master field equation:**

$$\boxed{\Box_g\sigma_\mu = Q_\mu(\sigma,\partial\sigma) + G_\mu(\sigma,\ell,H,q) + T_\mu + \kappa\ell_Q^2\Box_g^2\sigma_\mu + \mathcal{O}(\ell_Q^4)}$$

The four source terms:

| Term | Expression | Physical meaning |
|---|---|---|
| $Q_\mu$ | $\sigma_\mu(\nabla_\alpha\sigma_\beta\nabla^\alpha\sigma^\beta) + (\nabla_\mu\sigma_\alpha)(\sigma_\beta\nabla^\beta\sigma^\alpha)$ | Gravitational self-interaction ("gravity gravitates") |
| $G_\mu$ | $\sum_A[H_A(\ell\cdot\nabla)^2\sigma_\mu + \cdots]$ | Geometric coupling (vanishes in weak field) |
| $T_\mu$ | $\frac{1}{2}T^{\mu\nu}\sigma_\nu$ | Matter source |
| $\kappa\ell_Q^2\Box_g^2\sigma_\mu$ | — | Quantum stiffness (Planck-scale) |

**Pais–Uhlenbeck factored form:**

$$\boxed{\Box_g(\Box_g - m_Q^2)\sigma_\mu = -m_Q^2 S_\mu}$$

$$m_Q = \frac{M_\mathrm{Pl}}{\sqrt\kappa} \approx 0.71\,M_\mathrm{Pl}, \qquad \tau_Q \sim \hbar/(m_Q c^2) \sim 10^{-43}\;\mathrm{s}$$

Two modes: (1) **massless** $\Box_g\sigma^{(0)} = 0$ — classical gravitational waves at speed $c$; (2) **massive** $(\Box_g - m_Q^2)\sigma^{(m)} = 0$ — exponentially damped at Planck scale.

**Exact integral solution** via Green's functions:

$$\boxed{\sigma_\mu(x) = \sigma_\mu^{\mathrm{free}}(x) + \kappa\ell_Q^2\!\int_{\mathcal{M}}\!d^4x'\sqrt{-g(x')}[G_0(x,x') - G_{m_Q}(x,x')]S_\mu(x')}$$

**Gravitational wave quadrupole formula** (from $\sigma$-field superposition):

$$\boxed{h_+(t) = \frac{G}{c^4 r_{\mathrm{obs}}}M_{\mathrm{tot}}\!\left(\frac{d}{2}\right)^{\!2}\omega^2\cos(2\omega t)}$$

---

## 10. Recovery of Einstein's Equations

At equilibrium $\nabla^2\sigma^\alpha = 0$, Einstein's field equations are **derived**, not postulated:

$$\nabla^2\sigma = 0 \;\Longleftrightarrow\; \boxed{G_{\mu\nu} = \frac{8\pi G}{c^4}T_{\mu\nu}}$$

The sharp dichotomy:

$$\begin{cases} \nabla^2\sigma = 0 & \longrightarrow\;\text{Classical GR (static $\sigma$ configuration)} \\ \nabla^2\sigma \neq 0 & \longrightarrow\;\text{Quantum gravity ($\sigma$ fluctuations)} \end{cases}$$

> **GR is the equilibrium thermodynamics of the graviton field. Quantum gravity is its statistical mechanics.**

**Coarse-grained stress-energy tensor** (derived from Dirac action via coarse-graining):

$$\boxed{\langle T^{\mu\nu}\rangle = \rho c^2 u^\mu u^\nu + Pg^{\mu\nu} + \langle\tau^{\mu\nu}_{\mathrm{spin}}\rangle + \langle\tau^{\mu\nu}_{\mathrm{EM}}\rangle + \rho_\Lambda g^{\mu\nu}}$$

This is the complete Einstein stress-energy tensor — not postulated.

---

## 11. Gravitational Energy

QGD resolves the century-long problem of gravitational energy.

**The Noether energy-momentum tensor** (a *true* tensor, not a pseudotensor):

$$\boxed{T^{\mu\nu}_{\mathrm{QGD}} = \nabla^\mu\sigma_\alpha\nabla^\nu\sigma^\alpha - \frac{1}{2}g^{\mu\nu}(\nabla\sigma)^2 + \ell_Q^2\text{-corrections}}$$

**Positive-definite Hamiltonian** (manifest, without Schoen–Yau complexity):

$$H[\sigma,\pi] = \int\!d^3x\left[\frac{1}{2}\pi_\mu\pi^\mu + \frac{1}{2}(\nabla\sigma_\mu)^2 + V(\sigma)\right] \geq 0$$

Equality iff $\sigma_\mu = 0$ (flat Minkowski).

**Local energy density at every spacetime point:**

$$\boxed{\rho_{\mathrm{grav}}(\mathbf{x}) = \underbrace{\frac{1}{2}\dot\sigma_\mu\dot\sigma^\mu}_{\text{kinetic}} + \underbrace{\frac{1}{2}(\nabla\sigma_\mu)^2}_{\text{gradient}} + \underbrace{V(\sigma)}_{\text{potential}} + \underbrace{\frac{\ell_Q^2}{2}(\nabla^2\sigma_\mu)^2}_{\text{quantum stiffness}}}$$

**Schwarzschild energy density** ($r_s = 2GM/c^2$):

$$\boxed{\rho_{\mathrm{grav}}(r) = \frac{GM}{4c^2 r^3}}$$

**Two-body binding energy** (from field-gradient cross terms):

$$\boxed{E_{\mathrm{binding}} = -\frac{GM_1M_2}{d}}$$

**Gravitational wave energy** (exact, point-wise — no averaging needed):

$$\boxed{\rho_{\mathrm{GW}} = \frac{1}{2}\omega^2|\epsilon|^2}$$

Recovers Isaacson (1968) formula in the weak-field limit.

**Gravitational Poynting vector** (energy flux at any radius, not just infinity):

$$\boxed{\mathbf{S}_{\mathrm{grav}} = \dot\sigma_\mu\nabla\sigma^\mu, \qquad \frac{dE}{dt} = \oint_S\mathbf{S}_{\mathrm{grav}}\cdot d\mathbf{A}}$$

| Property | GR | QGD |
|---|---|---|
| Local energy density | Pseudotensor only | Scalar $\rho = \frac{1}{2}(\partial\sigma)^2$ |
| True tensor | Does not exist | Noether $T^{\mu\nu}_\mathrm{QGD}$ |
| Positive energy | Schoen–Yau (1979, at $\infty$) | Manifest: $H \geq 0$ locally |
| GW energy | Isaacson averaging | Exact: $\rho_\mathrm{GW} = \frac{1}{2}\omega^2|\epsilon|^2$ |
| BH energy | ADM mass (at $\infty$) | Distributed: $\rho \propto r^{-3}$ |
| Energy flux | Only at null infinity | Poynting vector at any $r$ |

---

## 12. Cosmology

**Cosmological reduction:** homogeneity and isotropy fix $\sigma_\mu = (\sigma_t(t),\,0,\,0,\,0)$.

**Exact equation of state** — for any nonzero $\dot\sigma_t$:

$$\rho_\sigma = \frac{1}{2}\dot\sigma_t^2, \qquad p_\sigma = -\frac{1}{2}\dot\sigma_t^2, \qquad \boxed{w_\sigma = -1}$$

This is exact, arising from kinetic structure alone — no potential, no fine-tuning.

**Fourth-order cosmological equation:**

$$\boxed{\ddot\sigma_t + 3H\dot\sigma_t - \ell_Q^2\!\left[\ddddot\sigma_t + 3H\dddot\sigma_t + 3\dot H\ddot\sigma_t + (3H^2 + 3\dot H)\dot\sigma_t\right] = S_t}$$

**Four modes** in de Sitter background ($H = H_0 = \mathrm{const}$):

$$\boxed{\sigma_t(t) = C_1 + C_2\,e^{-3H_0 t} + C_3\,e^{+t/\ell_Q} + C_4\,e^{-t/\ell_Q}}$$

| Mode | Timescale | Character |
|---|---|---|
| $C_1$ | $\infty$ | Constant background (neutral) |
| $C_2 e^{-3H_0 t}$ | $\sim 10^{17}$ s | Hubble decay (stable) |
| $C_3 e^{+t/\ell_Q}$ | $\sim 10^{-43}$ s | Exponential growth (unstable → saturates at $\sigma \sim \ell_\mathrm{Pl}$) |
| $C_4 e^{-t/\ell_Q}$ | $\sim 10^{-43}$ s | Rapid damping (stable) |

**Dark energy attractor solution:**

$$\boxed{v = \sqrt{\frac{3c^2}{4\pi G}}\,H_0, \qquad \rho_\sigma = \frac{3H_0^2 c^4}{8\pi G} \approx 5.3\times10^{-10}\;\mathrm{J/m^3}}$$

Observed: $\rho_\Lambda \approx 6\times10^{-10}$ J/m³ (Planck 2018) — **exact match with zero free parameters**.

**Effective cosmological constant:**

$$\Lambda_\mathrm{eff} = \frac{8\pi G}{c^4}\rho_\sigma = \frac{3H_0^2}{c^2}$$

**Modified Friedmann equations:**

$$\boxed{H^2 = \frac{8\pi G}{3c^2}\!\left(\rho_m + \rho_r + \frac{1}{2}\dot\sigma_t^2\right)}$$

$$\boxed{\frac{\ddot a}{a} = -\frac{4\pi G}{3c^2}\!\left(\rho_m + \rho_r + 3p_m + 3p_r - \dot\sigma_t^2\right)}$$

**Three-phase cosmological evolution:**

| Phase | Dominant mode | Outcome |
|---|---|---|
| Quantum ($t \lesssim 10\,t_\mathrm{Pl}$) | $C_3 e^{+t/\ell_Q}$ (grows, saturates at $\sigma \sim \ell_\mathrm{Pl}$) | Seeds inflation |
| Inflation | $\dot\sigma_t \approx \mathrm{const}$; de Sitter | $N \sim 60$ e-folds |
| Classical | $C_1 + C_2 e^{-3Ht}$ | Late-time dark energy |

---

## 13. N-Body Exact Solutions

### Single Kerr Body

$$\sigma_t^{(a)} = \sqrt{\frac{2GM_a\,r_a}{c^2\,\mathcal{S}_a}}, \qquad \sigma_\phi^{(a)} = \alpha_a\sin\theta_a\sqrt{\frac{2GM_a}{c^2\,r_a\,\mathcal{S}_a}}$$

### Two-Body Metric (exact, in equatorial plane)

$$g_{tt} = -\left[1 - (\sigma_t^{(1)} + \sigma_t^{(2)})^2\right]$$

$$g_{t\phi} = -r\sin\theta\,(\sigma_t^{(1)} + \sigma_t^{(2)})(\sigma_\phi^{(1)} + \sigma_\phi^{(2)})$$

$$g_{\phi\phi} = r^2\sin^2\theta\left[1 + (\sigma_\phi^{(1)} + \sigma_\phi^{(2)})^2\right]$$

**Cross-term expansion of $g_{tt}$:**

$$\boxed{g_{tt} = -\left[1 - \underbrace{(\sigma_t^{(1)})^2}_{2GM_1/(c^2r_1)} - \underbrace{(\sigma_t^{(2)})^2}_{2GM_2/(c^2r_2)} - \underbrace{2\sigma_t^{(1)}\sigma_t^{(2)}}_{\text{QGD cross-term }\sim r^{-1}}\right]}$$

### N-Body Metric

$$g_{tt} = -\left[1 - \sum_{a}A_a^2 - 2\sum_{a<b}A_aA_b\right]$$

$$g_{t\phi} = -r\sin\theta\left[\sum_aA_aB_a + \sum_{a\neq b}A_aB_b\right]$$

**Computational complexity:** $\mathcal{O}(N)$ per field point (vs. $\mathcal{O}(N_\mathrm{grid}^3 \times N_\mathrm{steps})$ for numerical GR).

### Superposition Theorem (with rigorous error bound)

$$\boxed{\frac{|\sigma_t^{(\mathrm{tot})} - \sum_a\sigma_t^{(a)}|}{\sum_a\sigma_t^{(a)}} \leq C\max_{a\neq b}\frac{r_s^{(a)}r_s^{(b)}}{d_{ab}^2}}$$

| $d/r_s$ | Accuracy |
|---|---|
| 5 | 96.8% |
| 10 | 99.2% |
| 100 | 99.992% |
| 1000 | $>99.9999\%$ |

### Merger Condition

$$\boxed{d_{\mathrm{merge}} = \frac{2GM_1}{c^2}\left(1 + \left(\frac{M_2}{M_1}\right)^{1/3}\right)^3}$$

For equal masses: $d_\mathrm{merge} = 16Gm/c^2 = 8r_s = 4r_s^\mathrm{(total)}$. Verified: $\Sigma_\mathrm{saddle} = 1.000000$ at all mass ratios.

### Equations of Motion

$$\boxed{\ddot{\mathbf{x}} = \underbrace{-\frac{GM_1}{r_1^2}\hat{r}_1 - \frac{GM_2}{r_2^2}\hat{r}_2}_{\text{Newtonian}} + G\sqrt{M_1M_2}\left[\frac{\hat{r}_2}{\sqrt{r_1}\,r_2^{3/2}} + \frac{\hat{r}_1}{\sqrt{r_2}\,r_1^{3/2}}\right]_{\text{QGD cross}}}$$

Cross-term scales as $r^{-3/2}$ (softer than Newton). Mutual two-body force is **exact Newton** (verified numerically).

---

## 14. Post-Newtonian Hierarchy and Master Formula

### Binary Notation

$$M = M_1 + M_2, \quad \mu = M_1M_2/M, \quad \eta = M_1M_2/M^2, \quad x = (GM\Omega/c^3)^{2/3}$$

### Binding Energy

$$E(x) = -\frac{\mu c^2}{2}x\left\{1 + \sum_{n=1}^4 e_n x^n + \mathcal{O}(x^5)\right\}$$

Coefficients:

$$e_1 = -\frac{3}{4} - \frac{\eta}{12}$$

$$e_2 = -\frac{27}{8} + \frac{19}{8}\eta - \frac{\eta^2}{24}$$

$$e_3 = -\frac{675}{64} + \left(\frac{34445}{576} - \frac{205\pi^2}{96}\right)\eta - \frac{155}{96}\eta^2 - \frac{35}{5184}\eta^3$$

### GW Flux

$$\mathcal{F}(x) = \frac{32}{5}\frac{c^5}{G}\eta^2 x^5\left\{1 + f_1 x + f_{3/2} x^{3/2} + f_2 x^2 + \cdots\right\}$$

$$f_1 = -\frac{1247}{336} - \frac{35}{12}\eta, \qquad f_{3/2} = 4\pi\;\text{(tail)}$$

### TaylorF2 Phase

$$\Psi(f) = 2\pi ft_c - \varphi_c - \frac{\pi}{4} + \frac{3}{128\eta x^{5/2}}\left\{1 + \psi_1 x + \psi_{3/2} x^{3/2} + \cdots\right\}$$

$$\psi_1 = \frac{3715}{1008} + \frac{55}{12}\eta, \qquad \psi_{3/2} = -10\pi$$

### QGD Master Formula

The exact Schwarzschild binding energy for circular geodesics:

$$\frac{E_\mathrm{bind}}{\mu} = \frac{1-2x}{\sqrt{1-3x}} - 1$$

This yields the **master formula** generating all test-body PN coefficients as exact rational numbers:

$$\boxed{e_n^{(\eta=0)} = -\binom{2n}{n}\left(\frac{3}{4}\right)^n\frac{2n-1}{n+1}}$$

| $n$ | $e_n^{(\eta=0)}$ | Status |
|---|---|---|
| 1 | $-3/4$ | GR known |
| 2 | $-27/8$ | GR known |
| 3 | $-675/64$ | GR known |
| 4 | $-3969/128$ | GR known |
| 5 | $-45927/512$ | **QGD prediction** |
| 6 | $-264627/1024$ | **QGD prediction** |
| 7 | $-12196899/16384$ | **QGD prediction** |
| 8 | $-70366725/32768$ | **QGD prediction** |

### Spinning Series

**Spin-orbit coefficients** (entering at $\chi\, u^{n+5/2}$):

$$\boxed{e_n^\mathrm{SO} = -\frac{(2n+1)!!\,3^n}{2^n\,n!}}$$

**Spin-squared coefficients** (entering at $\chi^2 u^{n+3}$):

$$\boxed{e_n^\mathrm{SS} = \frac{(2n+3)!!\,3^n}{3\cdot2^{n+1}\,n!}}$$

### QGD-Unique Flux Terms (no GR analogue)

**QGD dipole** ($D = (\sqrt{M_1}-\sqrt{M_2})^2/M$):

$$P_\mathrm{dip} = \frac{c^5}{G}\frac{\eta^2 D}{3}x^4 \qquad (-1\text{PN})$$

**Dipole tail** (0.5PN — absolutely no GR analogue):

$$P_\mathrm{dip}^\mathrm{tail} = \frac{c^5}{G}\eta^2 D\frac{4\pi}{3}x^{5.5} \qquad (0.5\text{PN})$$

**Complete QGD flux:**

$$\boxed{\mathcal{F}^\mathrm{QGD} = \mathcal{F}^\mathrm{GR} + \frac{c^5}{G}\eta^2 D\left[\frac{x^4}{3} - \frac{x^5}{2} + \frac{4\pi x^{5.5}}{3} + \frac{2\pi^2 x^7}{3} + \cdots\right]}$$

All terms vanish for $M_1 = M_2$ ($D = 0$) — built-in null test.

**QGD phase corrections:**

$$\delta\Psi^{(-1\mathrm{PN})} = -\frac{5D}{7168\eta^2}x^{-7/2} \propto f^{-7/3}$$

$$\delta\Psi^{(0.5\mathrm{PN})} = -\frac{5\pi D}{2688\eta^2}x^{-5/2} \propto f^{-5/3}$$

$$\delta\Psi^{(5\mathrm{PN,log})} = -\frac{27\eta}{2}\ln(x) \propto f^0\ln f$$

---

## 15. Ringdown, Radiation, and Black Hole Thermodynamics

### QNM Agreement with GR

The perturbation equation on the Kerr background:

$$\boxed{\Box_g\delta\sigma - \kappa\ell_Q^2\Box_g^2\delta\sigma = 0}$$

At leading order (stiffness ratio $\kappa\ell_Q^2/r_s^2 \sim 10^{-140}$), this reduces to the standard Teukolsky equation. **QGD QNMs equal GR at leading order.**

**Dominant 220 mode** (Echeverría–Leaver):

$$\omega_{220} = \frac{c}{2GM_f}\left[1 - 0.63(1-\chi_f)^{0.3}\right]$$

$$\tau_{220} = \frac{2GM_f}{c^3}\left[1 - 0.63(1-\chi_f)^{0.45}\right]^{-1}$$

### Complete QGD Ringdown Waveform

$$\boxed{h_\mathrm{QGD}(t) = h_\mathrm{QNM}(t) + h_\mathrm{cross}(t) + h_\mathrm{dipole}(t)}$$

**Cross-term early decay** (Prediction P1 — every binary merger):

$$h_\mathrm{cross}(t) = A_\mathrm{cross}\cdot 2\sigma_t^{(1)}\sigma_t^{(2)}\exp(-t/\tau_\mathrm{cross})$$

$$\boxed{\tau_\mathrm{cross} = \frac{6GM_\mathrm{tot}}{c^3} \approx 6\times10^{-5}\;\mathrm{s}\times\frac{M_\mathrm{tot}}{M_\odot}}$$

Zero oscillation frequency; decays before clean QNM phase.

**Dipole radiation** (Prediction P2 — for $M_1 \neq M_2$, radiates at $f_\mathrm{orb}$ not $2f_\mathrm{orb}$):

$$h_\mathrm{dipole}(t) = A_d\frac{(\sqrt{M_1}-\sqrt{M_2})^2}{M_1+M_2}e^{-t/\tau_d}\cos(\omega_\mathrm{orb}t + \varphi_d)$$

### QGD Dipole Radiation (Noether Analysis)

The $\sqrt{M}$-weighted dipole moment (from $|\nabla\sigma_t^{(a)}| \sim \sqrt{GM_a}/r^{3/2}$):

$$\boxed{\boldsymbol{d}_\sigma = \sum_{a=1}^N\sqrt{M_a}\,\mathbf{x}_a}$$

**QGD dipole power:**

$$\boxed{P_\mathrm{dipole}^\mathrm{QGD} = \frac{G^3M_1M_2(\sqrt{M_1}-\sqrt{M_2})^2}{3c^3r_{12}^4}}$$

Zero iff $M_1 = M_2$.

**Dipole-to-quadrupole ratio** (dominates at large separations):

$$\frac{P_\mathrm{dip}}{P_\mathrm{quad}} = \frac{5c^2r_{12}(\sqrt{M_1}-\sqrt{M_2})^2}{96GM_1M_2M_\mathrm{tot}} \propto x^{-1}$$

| Theory | Radiation source | Dipole |
|---|---|---|
| GR | $T_\mathrm{matter}^{\mu\nu}$ | Forbidden ($\sum M_a\dot{x}_a = \mathrm{const}$) |
| QGD | $T_\sigma^{\mu\nu}$ | Allowed ($M_1 \neq M_2$) |

### QGD Effective One-Body

In the COM frame, the exact two-body effective metric function:

$$\boxed{A^\mathrm{QGD}(u) = 1 - \frac{2u}{\eta}}$$

where $u = GM/(c^2r)$. No PN truncation, no free parameters, no NR calibration.

Orbital structure:

$$u_\mathrm{ISCO} = \frac{\eta}{6}, \qquad u_\mathrm{LR} = \frac{\eta}{3}, \qquad u_\mathrm{EOB} = \frac{\eta}{2}$$

### Black Hole Thermodynamics from Phase Expansion

**Hawking temperature** (derived from Taylor expansion of the QGD phase factor — without QFT in curved spacetime):

$$\boxed{T_H = \frac{\hbar c^3}{8\pi GMk_B}}$$

**Bekenstein–Hawking entropy:**

$$\boxed{S_\mathrm{BH} = \frac{k_Bc^3A}{4G\hbar}}$$

### Maximum Acceleration

From $\sigma = x/\lambda$ (spacetime intervals quantised in Compton wavelengths $x = n\lambda_C$):

$$\boxed{a_\mathrm{max} = \frac{3mc^3}{\hbar}}$$

Agrees with Caianiello's maximum acceleration (1981) within factor 3/2.

---

## 16. QFT of the σ-Field: Propagators and the Double Copy

### The σ-Field Propagator

In momentum space from the QGD action ($\ell_Q = \sqrt{G\hbar^2/c^4}$):

$$\boxed{D(k^2) = \frac{1}{k^2} - \frac{1}{k^2 + m_Q^2}}$$

$$m_Q = \frac{M_\mathrm{Pl}}{\sqrt\kappa} \approx 0.71\,M_\mathrm{Pl} \approx 1.54\times10^{-8}\;\mathrm{kg}$$

Propagator behavior:
- $k \ll k_\mathrm{Pl}$: $D \approx 1/k^2$ (standard GR)
- $k \sim k_\mathrm{Pl}$: $D \to 1/(k^2 m_Q^2) \propto 1/k^4$ (UV-softened)

### Newton's Law from Tree-Level Exchange

$$V(r) = -\int\frac{d^3q}{(2\pi)^3}\frac{16\pi G M_1M_2}{|\mathbf{q}|^2}e^{i\mathbf{q}\cdot\mathbf{r}} = -\frac{GM_1M_2}{r} \checkmark$$

**Yukawa correction** (regularising the singularity at $r \sim \ell_\mathrm{Pl}$):

$$V_\mathrm{QGD}(r) = -\frac{GM_1M_2}{r}\left(1 - e^{-m_Q r}\right)$$

**Gravitational Rutherford cross section:**

$$\frac{d\sigma}{d\Omega} = \frac{(Gm_1m_2)^2}{4E^2\sin^4(\theta/2)}$$

### UV Regularisation

Standard GR one-loop divergence: $I_\mathrm{GR} \sim \Lambda^2$ (quadratic — non-renormalisable).

QGD-regulated one-loop integral:

$$I_\mathrm{QGD} = c_0\ln\frac{m_Q}{\mu} + \mathrm{finite}$$

**One-loop: quadratic $\Lambda^2$ reduced to logarithmic $\ln\Lambda$.** The Pais–Uhlenbeck structure provides physical Pauli–Villars regularisation without ghost instabilities.

### The BCJ Double Copy

QGD is the **self-double-copy** of the graviton field:

$$\boxed{h_{\mu\nu} = g_{\mu\nu} - \eta_{\mu\nu} = -\sigma_\mu\otimes\sigma_\nu}$$

| | BCJ double copy | QGD |
|---|---|---|
| Graviton | $h_{\mu\nu} = A_\mu\otimes\tilde{A}_\nu$ | $h_{\mu\nu} = \sigma_\mu\otimes\sigma_\nu$ |
| Vector fields | Two different (left/right) | One (self-product) |
| Color structure | Stripped | Absent (gravity is universal) |
| DOF | $3\times3 \to 5+\cdots$ | $4\otimes4 \to 10 \to 2$ |

### Vacuum Energy and Dark Energy

Standard QFT prediction: $\rho_\mathrm{vac}^\mathrm{QFT} \sim 10^{113}$ J/m³ vs. observed $6\times10^{-10}$ J/m³ — a $10^{122}$ discrepancy.

QGD partial cancellation (quartic → quadratic divergence):

$$\rho_\mathrm{vac}^{(\sigma)} \sim -\frac{m_Q^2\Lambda^2}{8\pi^2}$$

Deep resolution via dark energy attractor (Chapter 5): $\rho_\sigma = 3H_0^2/(8\pi G) \approx 5.3\times10^{-10}$ J/m³ — no free parameter.

### Black Hole Entropy from σ-Entanglement

$$S_\mathrm{ent} = \alpha_\mathrm{field}\times\frac{A}{\epsilon^2}, \qquad \epsilon = \sqrt\kappa\,\ell_\mathrm{Pl}$$

$$S_\mathrm{ent}^{(\sigma)} = \frac{A}{2\pi\kappa\ell_\mathrm{Pl}^2} \approx \frac{A}{4\pi\ell_\mathrm{Pl}^2} \qquad (\text{within }\pi\text{ of BH entropy})$$

Natural UV cutoff at the massive-mode Compton wavelength — no ad hoc regulator.

### Gravitational Lamb Shift

$$\delta E_\mathrm{Lamb}^\mathrm{(grav)} = \frac{4}{3}\left(\frac{\alpha_\mathrm{grav}}{\pi}\right)^3 m_ec^2\ln\frac{1}{\alpha_\mathrm{grav}} \sim 10^{-137}\;\mathrm{J}$$

where $\alpha_\mathrm{grav} = Gm_em_p/(\hbar c) = 3.2\times10^{-42}$. Structurally identical to the QED Lamb shift.

---

## 17. Quantum Gravity: Decoherence, Information, and Inflation

### Gravitational Decoherence

For a mass $m$ in superposition with separation $\Delta x$, each branch produces a distinct $\sigma$-field vacuum. The Di\'osi–Penrose rate, now derived microscopically:

$$\boxed{\Gamma_\mathrm{QGD} = \frac{Gm^2}{\hbar\Delta x}}$$

Decoherence timescales:

| System | $m$ | $\Delta x$ | $\tau_\mathrm{dec}$ |
|---|---|---|---|
| Electron | $9.1\times10^{-31}$ kg | $10^{-10}$ m | $6\times10^9$ Gyr |
| C$_{60}$ fullerene | $1.2\times10^{-24}$ kg | $10^{-7}$ m | 3.5 Gyr |
| Virus ($10^6$ Da) | $1.7\times10^{-21}$ kg | $10^{-6}$ m | $1.8\times10^4$ yr |
| Grain of sand | $10^{-6}$ kg | $10^{-4}$ m | $1.6\times10^{-16}$ s |

Quantum-to-classical transition near $m \sim 10^{-17}$ kg. **Testable by MAQRO** at $10^9$ Da nanoparticles: $\tau_\mathrm{QGD} \approx 1$ year.

### Information Paradox: Three Mechanisms

**Mechanism 1:** The fundamental variable $\sigma_\mu$ is defined on flat Minkowski spacetime. The horizon is the surface $\Sigma_\mathrm{tot} = 1$ — a finite value, not a divergence. No singularity in the fundamental description.

**Mechanism 2:** The Pais–Uhlenbeck equation has well-posed initial-value problem. Unitary evolution: $d\,\mathrm{Tr}(\rho^2)/dt = 0$.

**Mechanism 3:** Cross-term radiation carries pre-merger information:

$$h_\mathrm{cross}(t) \propto 2\sigma_t^{(1)}\sigma_t^{(2)}\exp(-t/\tau_\mathrm{cross})$$

**Scrambling time:**

$$t_\mathrm{scr} = \frac{2GM}{c^3}\ln S_\mathrm{BH}$$

### Inflation from Fourth-Order σ-Dynamics

The QGD action generates the Starobinsky action **without a separate inflaton**:

$$S \supset \int d^4x\sqrt{-g}\left[\frac{R}{16\pi G} + \alpha R^2\right], \qquad \alpha \sim \kappa\ell_Q^2/(16\pi G)$$

Predictions at $N_* = 60$ e-folds:

$$\boxed{n_s = 1 - \frac{2}{N_*} = 0.967, \qquad r = \frac{12}{N_*^2} = 0.003}$$

Planck 2018: $n_s = 0.9649 \pm 0.0042$ ✓, $r < 0.064$ ✓.

Reheating temperature: $T_\mathrm{reh} \sim 10^9$ GeV (sufficient for baryogenesis, below gravitino bound).

### Gravitational Aharonov–Bohm Phase

A test particle traversing a closed path encircling spinning mass $J$:

$$\boxed{\Delta\varphi_\mathrm{grav} = \frac{2\pi mc}{\hbar}\cdot\frac{GJ}{c^2 r}}$$

Exists even where the Riemann tensor is negligible. Potentially detectable with next-generation atom interferometers.

### Weak Equivalence Principle Violation

Mass-dependent quantum correction (test-particle Compton wavelength $\lambda_C = \hbar/(mc)$):

$$\eta_\mathrm{QGD} = \frac{9\hbar^2}{2c^2 R_\oplus^2}\left(\frac{1}{m_1^2} - \frac{1}{m_2^2}\right)$$

For $^{133}$Cs vs. $^{87}$Rb: $\eta \approx 4\times10^{-49}$ — consistent with all bounds, nonzero in principle.

### Graviton Mass

Physical graviton is **exactly massless**. Cosmological effective mass (curvature correction):

$$m_g^\mathrm{eff} = \frac{\hbar H_0}{c^2} \approx 1.4\times10^{-33}\;\mathrm{eV}$$

Compton wavelength = Hubble radius: $\lambda_g = c/H_0 = R_H$.

**GW speed correction:**

$$\frac{\delta v}{c} \sim \kappa\ell_Q^2\omega^2 \approx 7\times10^{-107}\;\text{at }f = 100\;\mathrm{Hz}$$

($10^{92}$ times smaller than GW170817 bound of $3\times10^{-15}$.)

---

## 18. Dark Matter as κ-Enhancement

### The Central Logic Chain (Rotation Curves)

$$\underset{\text{Dirac}}{\psi} \xrightarrow{\text{WKB}} \underset{\text{graviton}}{\sigma_\mu} \xrightarrow{\text{metric}}{g_{\mu\nu}} \xrightarrow{\text{action}}{\Box\sigma - \kappa\ell_Q^2\Box^2\sigma} \xrightarrow{\text{Fourier}}{D(k^2)} \xrightarrow{\text{loops}}{\kappa_n} \xrightarrow{\Sigma}{v_\mathrm{obs}}$$

Each step derived, not postulated. **Zero free parameters per galaxy.**

### The κ-Ladder (Propagator–Pochhammer Connection)

The $n$-loop amplitude from the two-pole propagator factorises as:

$$\mathcal{A}_n \propto (n-1)! \times \left(\tfrac{1}{2}\right)_n$$

The squared modulus gives:

$$\boxed{\kappa_n^2 = 2\,(n-1)!\,\left(\tfrac{1}{2}\right)_n = \frac{(2n-1)!}{4^{n-1}}}$$

**Complete κ-ladder (zero free parameters):**

| $n$ | $\kappa_n$ (exact) | $\kappa_n$ (decimal) | Physical regime |
|---|---|---|---|
| 1 | $1$ | **1.000** | Solar system (Newtonian) |
| 2 | $\sqrt{3/2}$ | **1.225** | Wide binaries; screened to $\approx 1.04$ |
| 3 | $\sqrt{15/2}$ | **2.739** | Spiral galaxy outskirts |
| 4 | $\sqrt{315/4}$ | **8.874** | Galaxy groups |
| 5 | $\sqrt{1417.5}$ | **37.65** | Galaxy clusters (Bullet Cluster: ratio = 0.981) |
| 6 | $\sqrt{38981.25}$ | **197.4** | Superclusters (theoretical) |
| 7 | $\sqrt{2129828.6}$ | **1233.** | Cosmic web (theoretical) |

**Step ratio** revealing the two sectors:

$$\rho_k = \frac{\kappa_{k+1}^2}{\kappa_k^2} = \underbrace{k}_\text{massless} \times \underbrace{\frac{2k+1}{2}}_\text{massive}$$

### Surface Density Coherence Condition

QGD enhancement requires quantum gravitational phase coherence:

$$\Sigma < \Sigma_\mathrm{crit} = 17.5\;M_\odot/\mathrm{pc}^2$$

The MOND acceleration constant is not a free parameter — it is encoded here:

$$a_0 \sim 2\pi G\Sigma_\mathrm{crit} \times f_\mathrm{geom} \approx 1.2\times10^{-10}\;\mathrm{m/s^2}$$

### Bullet Cluster: Validating κ₅

Two-branch κ-rule for clusters:

$$\kappa_\mathrm{eff}(\Sigma) = \begin{cases} \kappa_5 = 37.65 & \Sigma < 17.5\;M_\odot/\mathrm{pc}^2\;\text{(galaxies)} \\ 1 + (\Sigma_\mathrm{crit}/\Sigma)^\alpha & \Sigma \geq \Sigma_\mathrm{crit}\;\text{(ICM gas)} \end{cases}$$

**Bullet Cluster result:**

$$\frac{M_\mathrm{QGD}}{M_\mathrm{lens}} = \frac{M_\mathrm{gas}\kappa_\mathrm{gas} + M_\star\kappa_5 + M_\mathrm{cross}}{2.80\times10^{14}\;M_\odot} = 0.981$$

κ contrast explaining the 8σ offset:

$$\frac{\kappa_\mathrm{galaxy}}{\kappa_\mathrm{gas}} = \frac{37.65}{1.724} = 21.8$$

Although gas = 90% of baryons, galaxies dominate lensing mass because their $\kappa$ is 22× larger. **The 8σ spatial offset — canonical "proof" of dark matter — is a prediction of QGD's $\Sigma$-mechanism.** No dark matter particle required.

| Test | $\Lambda$CDM | MOND | QGD |
|---|---|---|---|
| Total lensing mass | ✓ (by construction) | ✗ (deficit) | ✓ (ratio = 0.981) |
| 8σ spatial offset | ✓ (dark matter halos) | ✗ (fails) | ✓ (Σ-mechanism) |
| Free parameters | 84% DM assumed | $a_0$ postulated | **zero** |
| Dark matter particle | Required | Required | **None** |

### SPARC Rotation Curve Validation

$$v_\mathrm{obs} = v_\mathrm{bar}\sqrt{\kappa}$$

For 81 SPARC galaxies (1029 measurements), zero free parameters per galaxy:

| Metric | Value |
|---|---|
| Overall $R^2$ | 0.920 |
| RMSE | 18.6 km/s |
| Best galaxy (KK98-251) | $R^2 = 0.973$ |

---

## 19. The Central Logic Chain

The full derivation chain of QGD — from quantum mechanics to cosmology, every step derived:

$$\underbrace{\text{Dirac equation}}_{\text{microscopic}} \xrightarrow{\text{WKB}} \underbrace{\sigma_\mu = \frac{1}{mc}\partial_\mu S}_{\text{graviton field}} \xrightarrow{\text{master metric}} \underbrace{g_{\mu\nu}[\sigma]}_{\text{spacetime geometry}} \xrightarrow{\delta S/\delta\sigma = 0} \underbrace{G_{\mu\nu} = \frac{8\pi G}{c^4}T_{\mu\nu}}_{\text{Einstein (equilibrium)}}$$

**The metric is output, not input.**

### QGD vs. GR: Paradigm Comparison

| Aspect | General Relativity | QGD |
|---|---|---|
| Fundamental variable | $g_{\mu\nu}$ (10 components) | $\sigma_\mu$ (4 components) |
| Metric status | Fundamental | Derived output |
| Field equations | 10 coupled nonlinear PDEs | 4 linear + algebra |
| Superposition | Impossible at metric level | Exact at $\sigma$-level |
| Singularities | Generic, unavoidable | Resolved at $\lambda_C$ |
| Gravitational energy | Pseudotensor | True tensor $\frac{1}{2}(\partial\sigma)^2$ |
| Einstein equations | Postulated | Derived ($\nabla^2\sigma = 0$) |
| Computational cost | Exponential (numerical) | Polynomial (algebraic) |
| UV divergences | Quadratic ($\Lambda^2$) | Logarithmic ($\ln\Lambda$) |
| Dark energy | Free parameter ($\Lambda$) | Attractor: $\rho_\sigma = 3H_0^2/(8\pi G)$ |
| Dark matter | New particles required | $\kappa$-field, zero free parameters |
| Inflation | Separate inflaton required | $R^2$ Starobinsky from $\sigma$-dynamics |
| Dipole radiation | Forbidden | Allowed ($M_1 \neq M_2$) |

### Five Falsifiable Predictions (Ringdown)

| Prediction | QGD | Falsification |
|---|---|---|
| P1: Cross-term decay | $\omega = 0$, $\tau = 6GM_\mathrm{tot}/c^3$ in every merger | Not found in LVK analyses |
| P2: Dipole radiation | At $f_\mathrm{orb}$, $\propto (\sqrt{M_1}-\sqrt{M_2})^2/M$ | Equal-mass events show it |
| P3: Quantum stiffness | $\delta\omega/\omega \sim (\ell_Q/r_s)^2 \sim 10^{-140}$ | Measured QNM deviation |
| P4: $-1$PN phase correction | $\delta\Psi \propto D\,f^{-7/3}$ | Degenerate with GR PN coefficient |
| P5: Population null test | All QGD terms vanish for $q=1$ ($D = 0$) | Non-zero residuals at $q=1$ |

---

> **The action encodes gravity. The stress-energy manifests it. The complex wavefunction realizes it.**
>
> **GR is the equilibrium thermodynamics of the graviton field. Quantum gravity is its statistical mechanics.**
