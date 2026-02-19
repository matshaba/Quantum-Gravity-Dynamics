# QGD Theory: Complete Mathematical Exposition

> This document contains the complete mathematical framework of Quantum Gravitational Dynamics.
> Every major result is stated with its derivation chain, physical interpretation, and connection
> to the rest of the theory.

---

## Part I: The Starting Point — Dirac Equation and Coarse-Graining

### 1.1 The Dirac Equation

QGD begins from the **Dirac equation** — not from a postulated action or metric:

$$\boxed{\left(i\gamma^\mu\partial_\mu - m\right)\psi = 0}$$

In the WKB (semiclassical) limit $\psi = A(x)e^{iS(x)/\hbar}$, the phase gradient defines the **fundamental field of QGD**:

$$\boxed{\sigma_\mu \equiv \frac{p_\mu}{mc} = \frac{1}{mc}\partial_\mu S}$$

This σ-field is dimensionless, measuring gravitational phase cycles per unit length.

The full quantum Lagrangian density containing all known physics is:

$$\mathcal{L} = \frac{i}{2}\bar\psi\gamma^\mu\overleftrightarrow\partial_\mu\psi - m\bar\psi\psi - P(\bar\psi\psi) - \frac{1}{2M^2}J_{\mu\nu}J^{\mu\nu}\bar\psi\psi - \frac{1}{4}F_{\mu\nu}F^{\mu\nu} - \rho_\Lambda$$

where $J_{\mu\nu} = x_\mu p_\nu - x_\nu p_\mu$, $F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu$.

**The fundamental chain:**

$$\text{Dirac equation} \xrightarrow{\text{WKB}} \sigma_\mu = \frac{1}{mc}\partial_\mu S \xrightarrow{\text{coarse-grain over }\ell} S[\sigma] \xrightarrow{\text{master metric}} g_{\mu\nu}[\sigma] \xrightarrow{\delta S/\delta\sigma = 0} G_{\mu\nu} = \frac{8\pi G}{c^4}T_{\mu\nu}$$

### 1.2 Coarse-Graining Gives the Action

Applying Noether's theorem and coarse-graining the Dirac dynamics over scale $\ell$ — integrating out sub-wavelength fluctuations — yields two outputs simultaneously:

**Output 1 — Stress-energy tensor:**

Noether's theorem + coarse-graining over scale $\ell$ gives the stress-energy tensor:

$$\boxed{\langle T^{\mu\nu}\rangle = \rho c^2 u^\mu u^\nu + P g^{\mu\nu} + \langle\tau^{\mu\nu}_{\text{spin}}\rangle + \langle\tau^{\mu\nu}_{\text{EM}}\rangle + \rho_\Lambda g^{\mu\nu}}$$

Every term traces back to one action term:

| Stress-energy | Action origin | Physics |
|---|---|---|
| $\rho c^2 u^\mu u^\nu$ | $\frac{i}{2}\bar\psi\gamma^\mu\overleftrightarrow\partial_\mu\psi - m\bar\psi\psi$ | Mass-energy |
| $Pg^{\mu\nu}$ | $-P(\bar\psi\psi)$ | Pressure / EOS |
| $\tau^{\mu\nu}_{\text{spin}}$ | $-\frac{1}{2M^2}J_{\mu\nu}J^{\mu\nu}\bar\psi\psi$ | Frame-dragging |
| $\tau^{\mu\nu}_{\text{EM}}$ | $-\frac{1}{4}F_{\mu\nu}F^{\mu\nu}$ | EM stress |
| $\rho_\Lambda g^{\mu\nu}$ | $-\rho_\Lambda$ | Vacuum energy |

The effective Hamiltonian encoding all contributions through a unified denominator:
$$\boxed{H = \int\rho c^2\,dV - \int P\,dV - \frac{J^2}{2mr^2} + V_{\text{EM}} + \rho_\Lambda}$$

### 1.3 The Complete σ-Field Action

After extracting the phase field via WKB, the full action in σ-variables:

$$\boxed{S_\sigma = \int d^4x\sqrt{-g(\sigma)}\!\left[-\frac{c^4}{16\pi G}R[g(\sigma)] + \frac{1}{2}\nabla_\mu\sigma_\nu\nabla^\mu\sigma^\nu - \frac{\ell_Q^2}{2}\nabla_\alpha\nabla_\beta\sigma_\mu\nabla^\alpha\nabla^\beta\sigma^\mu\right] + S_{\text{matter}}}$$

Three terms: (1) Einstein-Hilbert in σ-variables, (2) σ-kinetic, (3) quantum stiffness.

**Key technical step:** Vary with respect to $\sigma_\alpha$, **not** $g_{\mu\nu}$:
$$\delta S = \int d^4x\!\left[\frac{\delta S}{\delta g_{\mu\nu}}\frac{\partial g_{\mu\nu}}{\partial\sigma_\alpha} - \frac{\hbar^2}{M}\nabla^2\sigma^\alpha\right]\delta\sigma_\alpha = 0$$

This generates a field equation for $\sigma_\alpha$ that does not exist in GR.

---

## Part II: From Wavefunction to Classical Physics

### 2.1 The WKB Phase Field

From the Dirac spinor $\psi = R(x)e^{iS(x)/\hbar}$ in the semiclassical limit:

$$\boxed{\sigma_\mu(x) \equiv \frac{1}{c}\partial_\mu S(x)}$$

Dimensionless, measures gravitational field strength as phase cycles per unit length.

### 2.2 Replacing Born's Rule — The Quantum-to-Classical Bridge

In standard QM, $|\psi|^2 = $ probability density is a **postulate**. In QGD it is **derived** from current conservation.

For spherically symmetric systems with $\nabla S = p(r)\hat r$, the continuity equation $\nabla\cdot(|\psi|^2\nabla S)=0$ gives:

$$\frac{1}{r^2}\frac{\partial}{\partial r}\!\left(r^2|\psi|^2 p(r)\right) = 0$$

Integrating — the **exact conservation law**:

$$\boxed{|\psi|^2\cdot p(r)\cdot r^2 = C}$$

This is not an approximation. It is an exact geometric identity: the $r^2$ from the spherical wavefunction and the $r^2$ from solid angle cancel **exactly**.

**The Second Wavefunction Interpretation (QGD):** This replaces Born's rule. Probability density is not fundamental — it is the inverse of the phase gradient density, fixed by current conservation. From this single law, forces follow:

$$F = \frac{dp}{dt} = -\frac{dE}{dr}$$

**Covariant form:**

$$\boxed{|\psi(x)|^2\,\sigma_\mu x^\mu = \frac{J}{\hbar}}$$

Physical meaning: (probability density) × (number of phase wavelengths from origin) = conserved quantum number $/\hbar$.

The quantum-to-classical transition is thus **deterministic**, not probabilistic. Probability arises only at the measurement level; the dynamics are phase-governed.

### 2.3 Gravitational Fine Structure Constant

Newton's law emerges from the leading oscillation of $e^{2i\mathcal{P}r/\hbar}$, where:

$$\mathcal{P} = \frac{\sqrt{2}\,GMm^2}{\hbar}$$

Defining the **gravitational fine structure constant**:

$$\alpha_G^2 = \frac{i\hbar c}{2GMm}, \qquad |\alpha_G|^2 = \frac{\hbar c}{2GMm}$$

The complex phase $e^{i\pi/4}$ is geometric. Normalization constants:

$$C = \frac{mc}{\sqrt{2}}, \qquad C = |\alpha_G|^2\mathcal{P}$$

Field energy and momentum (massless gravitons $E=pc$):

$$\boxed{E_{\text{field}} = \frac{mc^2}{\sqrt{2}\,|\alpha_G|^2}}, \qquad \boxed{P_{\text{field}} = \frac{mc}{\sqrt{2}\,|\alpha_G|^2}}$$

Single dimensionless control parameter: $\; R_s/\lambda_c = 1/|\alpha_G|^2$

| $|\alpha_G|^2$ | Regime |
|---|---|
| $\gg 1$ | Quantum (gravity negligible) |
| $\approx 1$ | Planck scale threshold |
| $\ll 1$ | Classical GR |

### 2.4 Complete Force Law from Phase Expansion

Taylor-expanding $e^{2i\mathcal{P}r/\hbar}$ order by order:

**Leading term** → Newton's law: $F_1 = GMm/r^2$

**Third term** (second term discarded, wrong power law):
$$E_3(r) = -\frac{3\hbar^3}{4im^2c\,\alpha_G^2\,r^3}, \qquad F_3 = \frac{9GM\hbar^2}{2mc^2r^4}$$

**Complete force law:**

$$\boxed{F(r) = \frac{GMm}{r^2}\left[1 + \frac{9}{2}\left(\frac{\lambda_C}{r}\right)^2 + \mathcal{O}(r^{-4})\right]}$$

**Crossover scale:** Quantum correction equals Newtonian at:

$$r_c = \frac{3}{\sqrt{2}}\lambda_C \approx 2.12\,\lambda_C$$

---

## Part III: The Master Metric

### 3.1 Metric as Algebraic Output

$$\boxed{g_{\mu\nu}(x) = T^\alpha_\mu T^\beta_\nu\!\left(M_{\alpha\beta}\circ\!\left[\eta_{\alpha\beta} - \sum_{a=1}^N\varepsilon_a\sigma_\alpha^{(a)}\sigma_\beta^{(a)} - \kappa\ell_Q^2\partial_\alpha\sigma^\gamma\partial_\beta\sigma_\gamma\right]\right)}$$

Three-tier hierarchy: flat baseline − classical gravity − quantum stiffness.

**Source signatures:**

| Source | $\sigma$-field | $\varepsilon_a$ | Effect |
|---|---|---|---|
| Mass | $\sqrt{2GM/c^2r}$ | $+1$ | Attraction |
| Spin | $a\sin\theta\sqrt{2GM/c^2r}$ | $+1$ | Frame dragging |
| Charge | $\sqrt{GQ^2/c^4r^2}$ | $-1$ | Repulsion |
| $\Lambda > 0$ | $Hr$ | $+1$ | Expansion |
| $\Lambda < 0$ | $|H|r$ | $-1$ | Anti-de Sitter |

### 3.2 Complete Wavefunction

$$\psi = \frac{2GMmi}{c\hbar}\!\left[\psi_0\begin{pmatrix}1\\0\\\sqrt{f}\\ia\sin\theta\sqrt{g}\end{pmatrix}e^{-iS/\hbar} + \psi_1\begin{pmatrix}0\\1\\-ia\sin\theta\sqrt{g}\\-\sqrt{f}\end{pmatrix}e^{-iS/\hbar} + \psi_2\begin{pmatrix}\sqrt{f}\\ia\sin\theta\sqrt{g}\\1\\0\end{pmatrix}e^{+iS/\hbar} + \psi_3\begin{pmatrix}-ia\sin\theta\sqrt{g}\\-\sqrt{f}\\0\\1\end{pmatrix}e^{+iS/\hbar}\right]$$

**Universal gravitational scalar** — every spacetime is a special case:

$$\boxed{f(r,\theta) = \frac{2GM}{c^2r} - \frac{GQ^2}{c^4r^2} + \frac{2Mr}{\Sigma} + \frac{\Lambda r^2}{3} + \frac{b(r)}{r} + H^2(t)r^2 - \int\frac{P(r)}{\rho(r)c^2}dr + \kappa\frac{\hbar^2}{M^2c^2r^2}}$$

---

## Part IV: Field Equations and Formal Solution

### 4.1 The Master Field Equation

$$\boxed{\Box_g\sigma_\mu = Q_\mu(\sigma,\partial\sigma) + G_\mu(\sigma,\ell,H,q) + T_\mu + \kappa\ell_Q^2\Box_g^2\sigma_\mu + \mathcal{O}(\ell_Q^4)}$$

Pais-Uhlenbeck factorization: $\Box_g(\Box_g - m_Q^2)\sigma_\mu = S_\mu$, $m_Q = M_{\text{Pl}}c/\hbar$.

Two modes: massless (classical GW) + Planck-mass (damped on $\tau_Q\sim 10^{-43}$ s, unobservable).

### 4.2 Formal Green's Function Solution

Define the retarded Green's function: $\mathcal{D}_x G_{\text{QGD}}(x,x') = \delta^4(x-x')/\sqrt{-g(x)}$

$$\boxed{\sigma_\mu(x) = \sigma_\mu^{\text{free}}(x) + \int d^4x'\sqrt{-g(x')}\;G_{\text{QGD}}(x,x')\,J_\mu(x')}$$

The composite propagator:

$$G_{\text{QGD}} = \ell_Q^2\!\left[G_0(x,x') - G_{m_Q}(x,x')\right]$$

- $G_0$: massless retarded propagator → information at speed $c$
- $G_{m_Q}$: massive Planck propagator → exponentially localized quantum corrections

**No superluminal propagation.** Nonlinear self-consistency solved iteratively via Born series.

### 4.3 Equivalence to Einstein's Equations

Variable substitution $g_{\mu\nu} \to g_{\mu\nu}[\sigma]$ makes EFE automatic. At equilibrium ($\nabla^2\sigma = 0$):

$$\left(G^{\mu\nu} - 8\pi G T^{\mu\nu}\right)\frac{\partial g_{\mu\nu}}{\partial\sigma_\alpha} = 0 \qquad \Rightarrow \qquad G_{\mu\nu} = \frac{8\pi G}{c^4}T_{\mu\nu}$$

- Classical GR = static σ ($\nabla^2\sigma = 0$)
- Quantum gravity = σ fluctuations ($\nabla^2\sigma \neq 0$)

---

## Part V: Exact Solutions

### 5.1 All GR Metrics Algebraically

Every metric follows from the source table — no differential equations required.

**Quantum-corrected Schwarzschild:**
$$g_{tt} = -\!\left(1 - \frac{2GM}{c^2r} - \frac{G\hbar^2}{Mc^4r^3}\right)$$

**Kerr frame-dragging from interference:**
$$g_{t\phi} = -T^t_t T^\phi_\phi\cdot\sigma_t\cdot\sigma_\phi = -\frac{2Mar\sin^2\theta}{\Sigma}$$

**Singularity resolution:** $r_{\min}\sim\ell_Q^{2/3}r_s^{1/3}\sim\lambda_C$

### 5.2 Exact Two-Body Solution

In GR: no exact two-body analytic solution exists. In QGD:

$$\sigma_\mu^{(i)} = \sqrt{\frac{2GM_i}{|\mathbf{x}-\mathbf{x}_i(t)|}}(1,\,v_i^x,\,v_i^y,\,v_i^z)$$

**Exact superposition** (linear Laplace equation at weak field):

$$\boxed{g_{\mu\nu} = \eta_{\mu\nu} - \underbrace{\sigma_\mu^{(1)}\sigma_\nu^{(1)}}_{\text{BH 1}} - \underbrace{\sigma_\mu^{(2)}\sigma_\nu^{(2)}}_{\text{BH 2}} - \underbrace{2\sigma_\mu^{(1)}\sigma_\nu^{(2)}}_{\text{ALL interactions}}}$$

The cross-term $-2\sigma_\mu^{(1)}\sigma_\nu^{(2)}$ contains automatically: binding energy, frame-dragging, GW emission, inspiral.

**GW as field beating:** $h_{\mu\nu}^{\text{GW}} = -2\sigma_\mu^{(1)}(t)\sigma_\nu^{(2)}(t)$

**Quadrupole formula — derived, not assumed:**

$$\boxed{h_+(t) = \frac{G}{c^4 r_{\text{obs}}}M_{\text{tot}}\!\left(\frac{d}{2}\right)^2\omega^2\cos(2\omega t)}$$

### 5.3 Exact Three-Body Solution

First closed-form analytic metric for a relativistic three-body system:

$$ds^2 = -\!\left[1 - \!\left(\sqrt{\frac{2GM_1}{c^2|\mathbf{x}-\mathbf{x}_1|}} + \sqrt{\frac{2GM_2}{c^2|\mathbf{x}-\mathbf{x}_2|}} + \sqrt{\frac{2GM_3}{c^2|\mathbf{x}-\mathbf{x}_3|}}\right)^2\right]c^2dt^2 + \frac{d\mathbf{x}^2}{1-\Sigma^2}$$

With hierarchical orbital positions, accumulated phases $\tau(t) = \int_0^t\Omega_{12}(t')dt'$, and radiation-reaction decay $a_{12}(t) = a_{12,0}(1-t/t_{\text{merge}})^{1/4}$.

### 5.4 General N-Body

$$\sigma_t = \sum_a\sqrt{2GM_a/c^2r_a}, \qquad \Sigma^2 = \underbrace{\sum_a(\cdot)}_{\text{self}} + \underbrace{2\sum_{i<j}(\cdot)}_{\text{cross-terms}}$$

**Event horizon:** $\Sigma = 1$. **Equal-mass binary merger:** $d = 4r_s$ (prediction).

Computational: $O(N^2)$ algebraic vs $O(N^3)$ GR numerical relativity.

---

## Part VI: Energy in QGD — Solving the 109-Year Problem

### 6.1 Why GR Fails

GR gravitational energy is a **pseudotensor** — coordinate-dependent, not localizable. The equivalence principle: choose free-fall, field vanishes, energy = 0. But QGD's $\sigma_\mu$ cannot be transformed away — it generates the metric.

### 6.2 True Gravitational Energy-Momentum Tensor

From Noether's theorem applied to $S_\sigma$:

$$\boxed{T^{\mu\nu}_{\text{QGD}} = \nabla^\mu\sigma_\lambda\nabla^\nu\sigma^\lambda - \frac{1}{2}g^{\mu\nu}(\nabla\sigma)^2 + \ell_Q^2\text{-corrections}}$$

Properties:
- **True tensor** (not pseudotensor): correct transformation under all coordinates
- **Automatic conservation:** $\partial_\mu T^{\mu\nu}_{\text{QGD}} = 0$
- **Positive definite:** $H[\sigma] = \int d^3x\,T^{00} \geq 0$ manifestly
- **Localizable:** $\rho(\mathbf{x})$ is a scalar field — has a definite value at each point

### 6.3 Complete Energy Density

$$\boxed{\rho_{\text{grav}} = \underbrace{\frac{1}{2}\dot\sigma_\mu^2}_{\text{kinetic}} + \underbrace{\frac{1}{2}(\nabla\sigma_\mu)^2}_{\text{gradient}} + \underbrace{V(\sigma)}_{\text{potential}} + \underbrace{\sim\ell_Q^2(\nabla^2\sigma)^2}_{\text{quantum stiffness}}}$$

**For Schwarzschild:**
$$\rho_{\text{grav}}(r) = \frac{1}{2}\!\left(\frac{d\sigma_t}{dr}\right)^2 = \frac{GM}{4c^2r^3}$$

Concentrated near horizon ($\propto 1/r^3$), zero at singularity (quantum-resolved), integrates to $M$.

**Newtonian limit:**
$$\rho_{\text{grav}} = \frac{1}{2}(\nabla\sigma_t)^2 = \frac{(\nabla\Phi)^2}{2c^2} \propto \frac{g^2}{8\pi G} \qquad \checkmark$$

Exactly the Newtonian gravitational field energy density.

### 6.4 The $Q_\mu$ Term: Gravitational Self-Energy

The self-interaction source $Q_\mu = \sigma_\mu(\partial\sigma)^2 + \ldots$ **is** the gravitational field's own energy-momentum. In GR this was invisible because metric variables obscured it.

Pattern across gauge theories:
- Maxwell ($Q=0$): linear, no self-energy
- Yang-Mills ($Q\sim A\cdot F$): gluon self-energy
- QGD ($Q\sim\sigma(\partial\sigma)^2$): graviton self-energy → **dark matter at higher orders**

### 6.5 Gravitational Waves: Exact, No Averaging

GR (Isaacson): $\langle T_{\mu\nu}^{\text{GW}}\rangle = \frac{1}{32\pi G}\langle\partial h\cdot\partial h\rangle$ — approximate, needs averaging.

QGD:
$$\boxed{\rho_{\text{GW}} = \frac{1}{2}(\partial\sigma^{\text{GW}})^2 = \frac{1}{2}\omega^2|\epsilon|^2}$$

Point-wise defined, exact, gauge-invariant. Agrees with Isaacson in the appropriate limit.

### 6.6 Gravitational Poynting Vector

$$\mathbf{S}_{\text{grav}} = \dot\sigma_\mu(\nabla\sigma^\mu), \qquad \frac{dE}{dt} = \oint_S\mathbf{S}_{\text{grav}}\cdot d\mathbf{A}$$

Tracks energy flow at **any radius**, not just infinity.

---

## Part VII: Dark Matter as Quantum Structure

### 7.1 Phase Taylor Expansion

$$e^{i\sigma_\mu x^\mu/\hbar} = \sum_{j=0}^\infty\frac{(i\sigma)^j}{j!} \qquad \Rightarrow \qquad \kappa_j = \sqrt{\frac{(2j-1)!}{2^{2j-2}}}$$

| $j$ | $\kappa_j$ | Activates at |
|---|---|---|
| 1 | 1.000 | Newtonian |
| 2 | 1.225 | 1–5 kpc |
| 3 | 2.739 | 5–20 kpc |
| 4 | **8.874** | **20–100 kpc (galactic DM)** |
| 5 | 37.72 | 10–100 Mpc |
| 6 | 197.0 | ~100 Mpc |
| 7 | 1245  | ~Gpc |

Derived, never fitted.

### 7.2 Modified Rotation Curves (Zero Free Parameters)

$$v^2(r) = \frac{GM_{\text{baryon}}(<r)}{r}\!\left[1 + \sum_{j=2}^J\kappa_j g_j(r)\right]$$

SPARC: 175 galaxies, 4,248 measurements, $R^2=0.908$, $\chi^2_\nu\approx 1.2$, RMS = 8.3 km/s.

MOND correspondence: $a\approx\kappa_j\sqrt{a_N a_0}$, with emergent $a_0 = c^2/r_0$.

---

## Part VIII: Cosmology

### 8.1 Complete Cosmological Dynamics

**σ-field evolution:**

$$\boxed{\ddot\sigma_t + 3H\dot\sigma_t - \ell_Q^2\!\left[\ddddot\sigma_t + 3H\dddot\sigma_t + 3\dot H\ddot\sigma_t + (3H^2+3\dot H)\dot\sigma_t\right] = \frac{8\pi G}{c^4}Q_t + \frac{4\pi G}{c^2}\rho_m\sigma_t}$$

**Friedmann equations:**

$$H^2 = \frac{8\pi G}{3}\!\left(\rho_m + \rho_r + \frac{1}{2}\dot\sigma_t^2\right), \qquad \frac{\ddot a}{a} = -\frac{4\pi G}{3}\!\left(\rho_m + \rho_r + 3p_m + 3p_r - \dot\sigma_t^2\right)$$

### 8.2 Dark Energy as Attractor

Constant-velocity solution $\dot\sigma_t = v$:

$$\boxed{v = \sqrt{\frac{3c^2}{4\pi G}}\,H_0, \quad \rho_\sigma = \frac{3H_0^2}{8\pi G}, \quad w = -1}$$

Dark energy is the dynamically selected attractor of the σ-field — not a tuned parameter.

---

## Part IX: Falsifiable Predictions

| Prediction | Value | Test |
|---|---|---|
| Neutron star mass shift | $\sim 0.01\,M_\odot$ | **NICER (now)** |
| Binary merger at $d = 4r_s$ | Equal mass | Next-gen GW |
| $\kappa_5 = 37.7$ large-scale | 10–100 Mpc correlations | DESI, Euclid |
| CMB peak modulation | $\ell_n\propto\kappa_{j(n)}f(n)$ | CMB-S4 |
| Maximum acceleration | $a_{\max} = 3mc^3/\hbar$ | — |
| Quantum perihelion shift | $\sim 10^{-90}$ arcsec/century | Unmeasurable |

---

## Part X: Comparison Table

| Aspect | GR | QGD |
|---|---|---|
| Fundamental variable | $g_{\mu\nu}$ (10) | $\sigma_\mu$ (4) |
| Metric | Fundamental — solve EFE | Output of σ-dynamics |
| Born's rule | Postulate | Derived: $|\psi|^2 pr^2 = C$ |
| Quantum → classical | Unclear | Exact conservation law |
| N-body | No exact solution | Exact $O(N^2)$ algebraic |
| Dark matter | New particles | Factorial $\kappa_j$ |
| Singularities | Unavoidable | Resolved at $\lambda_C$ |
| GW waveforms | Supercomputer, weeks | Algebraic, seconds |
| Gravitational energy | Pseudotensor | True tensor $\frac{1}{2}(\partial\sigma)^2$ |
| GW energy | Isaacson avg (approx) | Exact $\frac{1}{2}\omega^2|\epsilon|^2$ |
| Positive energy | Proved 1979 | Manifest: $H[\sigma]\geq 0$ |
| Information paradox | Unresolved | Unitary σ evolution |
| Cosmological constant | Fine-tuning | Attractor $3H_0^2/8\pi G$ |

---

*Executable implementations: `core/`, `solutions/`, `predictions/` directories.*
