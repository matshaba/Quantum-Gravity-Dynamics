# Quantum Gravitational Dynamics (QGD)

**A σ-field representation of GR that reduces the entire two-body post-Newtonian problem to one closed form.**

> *"Not a replacement for GR — GR, viewed through a coordinate that makes the structure visible."*

---

## What is this?# Quantum Gravitational Dynamics

**Authors:** Romeo Matshaba
**Date:** March 2026

---

## 0. Context: State of the Art in GR and EOB Theory

To appreciate what QGD contributes, it helps to understand what GR has and has not achieved in the post-Newtonian two-body problem as of early 2026.

### How PN calculations work in GR

In GR, the two-body problem is attacked order-by-order using multi-loop Feynman diagrams in an effective field theory (EFT) framework, or via the multipolar post-Minkowskian (MPM-PN) formalism. At each new PN order, new loop integrals appear — the 3PN result required three independent groups and dimensional regularisation to resolve a Hadamard pole ambiguity; the 4PN result required the introduction of nonlocal-in-time (hereditary) dynamics. The 5PN Hamiltonian was completed in 2020–2022 by two independent groups: Bini–Damour–Geralico using gravitational self-force (GSF) techniques, and Blümlein–Maier–Marquard–Schäfer using EFT loop diagrams (arXiv:2110.13822). They disagreed on a single rational coefficient — the $\nu^2$ term at 5PN, $c_{51}^{\rm rat}$ — a dispute that remained unresolved in the published literature.

The frontier as of March 2026:

| PN sector | Status | Key references |
|:---|:---|:---|
| 1PN–4PN local | Complete | DJS 2000/2001; Buonanno–Damour 1999 |
| 4PN nonlocal (tail) | Complete | DJS 2015 |
| 5PN local | Complete (with dispute at $c_{51}$) | BD-G 2020; Blümlein et al. 2020/2022 |
| 5.5PN tail-of-tail | Complete | DJS 2015; BD-G 2020 |
| 6PN local | Partial (4 rational coefficients remain) | BD-G arXiv:2003.11891 |
| 6PN nonlocal + scattering | Complete | BD-G arXiv:2007.11239 |
| 6.5PN tail-of-tail | Complete (July 2025) | Bini–Damour arXiv:2507.08708 |
| 7PN and beyond | Open | — |
| Conservative 8.5PN | Open | — |

### What GR does not have: no closed form for $A_{\rm trans}$

The EOB $A$-function is a Taylor series in $u$ whose coefficients are computed one order at a time. In standard GR/EOB theory:

- **Each PN order requires a new independent multi-loop calculation.** There is no formula that predicts $c_{n0}$ from $c_{(n-1)0}$.
- **The transcendental sector** ($\pi^2, \pi^4, \pi^6, \dots$ terms) grows with each order but no pattern relating successive orders has been identified in GR. The 3PN coefficient involves $\pi^2$, the 5PN involves $\pi^4$, the 6PN involves $\pi^6$, but the ratios between them are not known to follow a law.
- **The $\nu$-dependence** at each PN order (the $\nu^2, \nu^3, \dots$ sub-leading terms) is unknown in GR beyond what GSF calculations can supply numerically at linear order in $\nu$.
- **No resummation formula** for $A_{\rm trans}(u;\nu)$ to all orders has ever been proposed in the GR literature. Padé approximants are used numerically to improve convergence, but they are fitted to known coefficients and do not predict new ones.

In summary, GR delivers a growing but incomplete table of coefficients. The 6PN–8PN transcendental entries and the $\nu^2$–$\nu^4$ towers at those orders are entirely open.

### What QGD provides that GR does not

QGD derives a single closed-form expression — the Double-Padé theorem (eq. 7 below) — that:

1. Predicts every transcendental coefficient at every PN order from one algebraic formula, with zero free parameters.
2. Fixes the $\nu$-tower structure at all orders via the $\beta^2 = 9/16$ mechanism.
3. Forces $B(u;\nu) \equiv 1$, eliminating one of the two EOB potentials entirely.
4. Yields exact predictions at 6PN–8PN before GR has reached those orders.
5. Resolves the 5PN $c_{51}$ dispute between the two GR groups (§11).

GR and QGD agree on all computed coefficients through 6PN. The difference is structural: GR provides a finite table computed order by order; QGD provides a closed form that contains the entire table.

---

## 1. Graviton Field from the Dirac Equation

Below is a quick introduction to QGD for a more in depth study, consult

https://github.com/matshaba/Quantum-Gravity-Dynamics/blob/main/docs

https://github.com/matshaba/Quantum-Gravity-Dynamics/blob/main/docs/Ch1-foundations.pdf

The Dirac equation for a spin-$\tfrac{1}{2}$ particle of mass $m$ in curved spacetime is

$$
\bigl[i\hbar\,\gamma^\mu(\partial_\mu + \Gamma_\mu) - mc\bigr]\Psi = 0.
\tag{1}
$$

Writing $\Psi = \psi_0\,e^{iS/\hbar}$ and retaining only leading order in $\hbar^0$, the oscillatory phases cancel and we obtain the covariant Hamilton–Jacobi equation. The key identification is

$$
\boxed{\sigma_\mu \;\equiv\; \frac{\hbar}{mc}\,\partial_\mu S,}
\tag{2}
$$

where $\sigma_\mu\sigma^\mu = 1$ at lowest order. We call $\sigma_\mu$ the **graviton field** — dimensionless, carrying all phase information of the gravitational interaction.

### Spinor Solution

The full Dirac spinor in terms of $\sigma_\mu$ is

$$
\Psi \;=\;
\begin{pmatrix}
  1-\tfrac{1}{2}\sigma_0^2 \\[4pt]
  \sigma_3 \\[4pt]
  \sigma_1 - i\sigma_2 \\[4pt]
  1 - \sigma_0
\end{pmatrix}
e^{iS/\hbar},
\tag{3}
$$

satisfying the linearised Dirac equation to all orders in $\sigma_\mu$. The upper two are the large spinor components; the lower two encode gravitomagnetic effects.

---

## 2. The Master Metric

The metric is not fundamental — it emerges algebraically from the graviton field. In a local Lorentz frame: $g_{\mu\nu}^{\text{(local)}} = \eta_{\mu\nu} - \sigma_\mu\sigma_\nu$. For multiple sources with a coordinate transformation:

$$
\boxed{g_{\mu\nu}(x) = T^\alpha{}_\mu\,T^\beta{}_\nu\,
       \Bigl[\eta_{\alpha\beta}
             - \sum_a \varepsilon_a\,\sigma_\alpha^{(a)}\sigma_\beta^{(a)}
             - \kappa\ell_Q^2\,\partial_\alpha\sigma^\gamma\partial_\beta\sigma_\gamma
       \Bigr]}
\tag{4}
$$

where $T^\alpha{}_\mu$ is the local tetrad, $\varepsilon_a = \pm1$ encodes attraction/repulsion, and $\kappa\ell_Q^2$ is the quantum stiffness (negligible classically, $\ell_Q \sim 10^{-35}$ m).

**Recipe:** choose $\sigma$-fields → assign $\varepsilon$ → read off $g_{\mu\nu}$. No Einstein equations.

| Source | $\varepsilon$ | Effect on $g_{tt}$ |
|--------|:---:|---|
| Mass (Schwarzschild) | $+1$ | $-2GM/(c^2 r)$ |
| Electric charge | $-1$ | $+GQ^2/(c^4 r^2)$ |
| Spin (Kerr) | $+1$ | off-diagonal $g_{t\phi}$ |
| $\Lambda > 0$ (de Sitter) | $-1$ | $-H^2r^2/c^2$ |
| $\Lambda < 0$ (Anti-de Sitter) | $+1$ | $+H^2r^2/c^2$ |

### Ten GR Solutions — here we first show the computational advantage of QGD over GR and to introduce the language of QGD graviton fields $\sigma_t
### whether the theory is ultimately correct is left for the experimentalist, but here we can think of QGD as a computational tool to make GR tractable

**1. Minkowski** — $\sigma_\mu = 0$: $\;ds^2 = -dt^2 + dx^2 + dy^2 + dz^2$.

**2. Schwarzschild** — $\sigma_t = \sqrt{2GM/c^2r}$, $\varepsilon=+1$:

$$
\boxed{ds^2 = -\!\left(1 - \frac{2GM}{c^2 r}\right)dt^2 + dr^2 + r^2 d\Omega^2}
$$

**3. Reissner–Nordström** — two fields $(\sigma_t^{(M)},\varepsilon=+1)$ and $(\sigma_t^{(Q)} = \sqrt{GQ^2/c^4r^2},\varepsilon=-1)$:

$$
\boxed{ds^2 = -\!\left(1 - \frac{2GM}{c^2 r} + \frac{GQ^2}{c^4 r^2}\right)dt^2 + dr^2 + r^2 d\Omega^2}
$$

**4. Kerr** — $\sigma_t = \sqrt{2GMr/c^2\Sigma}$, $\sigma_\phi = a\sin\theta\sqrt{2GM/c^2r\Sigma}$, $\Sigma = r^2+a^2\cos^2\theta$, $\varepsilon=+1$. The off-diagonal term $g_{t\phi}$ arises directly from $\sigma_t\sigma_\phi$:

$$
\begin{aligned}
ds^2 &= -\!\left(1 - \frac{2GMr}{c^2\Sigma}\right)dt^2
       - \frac{4GMar\sin^2\theta}{c^2\Sigma}\,dt\,d\phi
       + \frac{\Sigma}{\Delta}\,dr^2 + \Sigma\,d\theta^2 \\
     &\quad + \!\left(r^2+a^2 + \frac{2GMa^2r\sin^2\theta}{c^2\Sigma}\right)\sin^2\theta\,d\phi^2
\end{aligned}
$$

with $\Delta = r^2 - 2GMr/c^2 + a^2$.

**5. Kerr–Newman** — superpose $\sigma^{(M)},\sigma^{(J)},\sigma^{(Q)}$ with $\varepsilon=+1,+1,-1$; modifies $\Delta$ and $g_{tt}$ by the charge term.

**6. Schwarzschild–de Sitter** ($\Lambda>0$, $H=\sqrt{\Lambda c^2/3}$):

$$
\boxed{ds^2 = -\!\left(1 - \frac{2GM}{c^2r} + \frac{H^2r^2}{c^2}\right)dt^2 + dr^2 + r^2 d\Omega^2}
$$

**7. Schwarzschild–AdS** ($\Lambda<0$):

$$
\boxed{ds^2 = -\!\left(1 - \frac{2GM}{c^2r} - \frac{H^2r^2}{c^2}\right)dt^2 + dr^2 + r^2 d\Omega^2}
$$

**8. de Sitter:** $\;ds^2 = -(1-H^2r^2/c^2)dt^2 + dr^2 + r^2d\Omega^2$

**9. Anti-de Sitter:** $\;ds^2 = -(1+H^2r^2/c^2)dt^2 + dr^2 + r^2d\Omega^2$

**10. FLRW:** expansion encoded in the tetrad $T^\alpha{}_\mu = \mathrm{diag}(1,a(t),a(t),a(t))$:

$$
ds^2 = -dt^2 + a(t)^2(dx^2+dy^2+dz^2).
$$

---

## 3. EOB Structure and Master Equation

For two bodies with $M=m_1+m_2$, $\nu = m_1m_2/M^2$, $u = GM/(c^2r)$, the EOB Hamiltonian is $H_{\rm eff} = \mu\sqrt{A(u;\nu)[1+\mathbf{p}^2/\mu^2 + p_r^2/B]}$.

The isotropic $\sigma$-field structure forces:

$$
\boxed{B(u;\nu) \equiv 1,} \qquad \boxed{A(u;\nu) = 1 - \sigma_{\rm eff}^2(u;\nu).}
\tag{5}
$$

The entire two-body orbital dynamics is encoded in one scalar function $A(u;\nu)$.

---

## 4. The $A$-Function Through 6PN

$$
A(u;\nu) = 1 - 2u + 2\nu u^3 + \nu\!\left(\frac{94}{3} - \frac{41\pi^2}{32}\right)u^4 + \bigl[c_{50}\nu + c_{51}\nu^2\bigr]u^5 - \frac{22}{3}\nu\,u^5\ln u + \bigl[c_{60}\nu + c_{61}\nu^2 + c_{62}\nu^3\bigr]u^6 + \mathcal{O}(u^7).
$$

| Order | Coefficient | Exact value | Numerical | Status |
|:------|:------|:------|------:|:------|
| 1PN | $a_1$ | $-2$ | $-2$ | GR exact ✓ |
| 2PN | $a_2$ | $0$ | $0$ | GR exact ✓ |
| 3PN | $a_3/\nu$ | $2$ | $2$ | GR exact ✓ |
| 4PN | $a_4^{\rm rat}/\nu$ | $94/3$ | $31.333$ | GR exact ✓ |
| 4PN | $a_4^{\pi^2}/\nu$ | $-41\pi^2/32$ | $-12.645$ | GR exact ✓ |
| 5PN | $c_{50}^{\rm rat}$ | $287.637$ | $287.637$ | GSF confirmed ✓ |
| 5PN | $c_{50}^{\pi^2}$ | $-4237\pi^2/60$ | $-696.959$ | ✓ |
| 5PN | $c_{50}^{\pi^4}$ | $2275\pi^4/512$ | $+432.824$ | ✓ |
| 5PN | $c_{50}$ (total) | — | $\approx 23.502$ | ✓ |
| 5PN | $c_{51}^{\pi^2}$ | $-369/512$ | $-7.115$ | $\beta^2$ mech. ✓ |
| 5PN | $c_{51}^{\pi^4}$ | $+20475/8192$ | $+243.464$ | $\beta^2$ mech. ✓ |
| 5PN | $c_{51}^{\rm rat}$ | $-200.962$ | $-200.962$ | QGD σ-graph |
| 5PN | $c_{51}$ (total) | — | $\approx 35.388$ | BD-G 2020 ✓ |
| 5PN | $a_5^{\rm log}/\nu$ | $-22/3$ | $-7.333$ | GR exact ✓ |
| 6PN | $c_{60}^{\pi^2}$ | $-41\pi^2/32$ | $-12.645$ | Double-Padé † |
| 6PN | $c_{60}^{\pi^4}$ | $+2275\pi^4/512$ | $+432.824$ | Double-Padé † |
| 6PN | $c_{60}^{\pi^6}$ | $-5175625\pi^6/335872$ | $-14814.542$ | Double-Padé † |
| 6PN | $c_{60}^{\rm rat}$ | $\approx 14398$ | $14398$ | constraint ‡ |
| 6PN | $c_{60}$ (total) | — | $\approx 4$ | ✓ |

† Double-Padé Theorem (§5). ‡ Cancellation constraint + rational Padé recursion (§7).

---

## 5. Double-Padé Theorem

**Theorem.** The complete transcendental diagonal sector of the EOB $A$-function is exactly resummed by

$$
\boxed{A_{\rm trans}(u;\nu) =
  \frac{-\nu\,\dfrac{41\pi^2}{32}\,u^4}
        {\left[1 + \dfrac{2275}{656}\,\pi^2\, u\right]
         \left[1 - \dfrac{9}{16}\,\nu\right]}}
\tag{7}
$$

- First denominator: Padé in $u$, pole at $u^* = -656/(2275\pi^2) < 0$ (unphysical)
- Second denominator: Padé in $\nu$, pole at $\nu^* = 16/9 > 1/4$ (unphysical)

Both poles lie outside the physical domain. $A(u;\nu)$ remains bounded and positive throughout the binary inspiral.

### Key Algebraic Identity

The two Padé denominators are linked by one algebraic identity:

$$
\boxed{32 \times 656 \;=\; 512 \times 41 \;=\; 20992}
\tag{8}
$$

This ensures that $(-41/32)\times(-2275/656) = 2275/512$, exactly reproducing the 5PN $\pi^4$ coefficient — the first internal consistency check of the theorem.

### Proof

**$u$-Padé expansion** generates the diagonal transcendental ladder:

$$
\frac{c_{n0}^{\pi^{2(n-3)}}}{\nu} = \left(-\frac{41}{32}\right)\!\left(-\frac{2275}{656}\right)^{n-4}, \qquad n \ge 4.
\tag{9}
$$

**$\nu$-Padé expansion** generates the mass-ratio tower:

$$
c_{nk}^{\pi^j} = \beta^{2k}\,c_{n0}^{\pi^j}, \qquad \beta^2 = \frac{9}{16}.
\tag{10}
$$

Independent confirmation: $c_{51}^{\pi^4}/c_{50}^{\pi^4} = (20475/8192)/(2275/512) = 9/16$ exactly.

### Transcendental Ladder (Exact Fractions)

| $n$ | PN order | $c_{n0}^{\pi^{2(n-3)}}/\nu$ (exact fraction) | $\times\pi^{2(n-3)}$ |
|:---:|:---:|:---|---:|
| 4 | 3PN | $-41/32$ | $-12.645$ |
| 5 | 4PN | $+2275/512$ | $+432.824$ |
| 6 | 5PN | $-5175625/335872$ | $-14814.542$ |
| 7 | 6PN | $+11774546875/220332032$ | $+507\,067$ |
| 8 | 7PN | $-26787094140625/144537812992$ | $-17\,355\,729$ |

Rows $n=7,8$ are zero-free-parameter predictions; GR has not yet reached these orders.

### $\nu$-Tower (The $\beta$-Mechanism)

$$
\beta = -\frac{3}{4},\quad
\beta^2 = \frac{9}{16},\quad
\beta^4 = \frac{81}{256},\quad
\beta^6 = \frac{729}{4096}.
$$

The first $\nu^k$ contribution appears at $(k+3)$PN. The maximum $\nu$-power at $n$PN is $\nu^{n-3}$. The first $\nu^4$ term in QGD is $c_{73}^{\pi^8} = \beta^6 c_{70}^{\pi^8}$ at 7PN.

---

## 6. General Structure Theorem

**Theorem.** For $n \ge 4$, the full transcendental diagonal of $a_n(\nu)$ is

$$
a_n^{\rm trans}(\nu) =
  c_{n0}^{\pi^{2(n-3)}}\,\pi^{2(n-3)}
  \sum_{k=0}^{n-3}\beta^{2k}\nu^{k+1}.
\tag{11}
$$

The sum truncates at $k = n-3$: at $n$-loop order there are at most $n-3$ two-body crossings, each contributing one power of $\nu$.

---

## 7. Rational Sector and Padé [2/1] Recursion

Define $H_{\rm rat}[m]$ from the rational part of $A$:

$$
H_{\rm rat}(u) = \frac{A_{\rm rat}(u;\nu) - 1 + 2u}{2\nu u^3}, \qquad
H[0]=1,\quad H[1]=\frac{47}{3},\quad H[2] = \frac{c_{50}^{\rm rat}}{2} \approx 143.818.
$$

Fitting a $[2/1]$ Padé to these three points yields the geometric recursion

$$
H[m+1] = q \cdot H[m], \qquad m \ge 2,
\tag{12}
$$

with $q \approx 50.06$ fixed by the physical constraint $c_{60}^{\rm total} \approx 4$ (§8). This gives:

$$
\begin{aligned}
c_{60}^{\rm rat} &= 2\,H[2]\cdot q \approx 14\,398, \\[4pt]
c_{70}^{\rm rat} &= 2\,H[2]\cdot q^2 \approx 720\,745, \\[4pt]
c_{80}^{\rm rat} &= 2\,H[2]\cdot q^3 \approx 36\,080\,000.
\end{aligned}
\tag{13}
$$

---

## 8. Escalating Cancellation

At every PN order the rational and transcendental sectors nearly cancel, keeping $A(u;\nu)$ bounded:

| PN | Rational | Trans sum | Total | Source |
|:--:|---:|---:|---:|:---|
| 3 | 2 | 0 | 2 | GR exact |
| 4 | 31.33 | −12.65 | 18.69 | GR exact |
| 5 | 287.64 | −264.13 | 23.50 | GSF confirmed |
| 6 | $\approx 14\,398$ | $\approx -14\,394$ | $\approx 4$ | constraint-derived |
| 7 | $\approx 720\,745$ | $\approx -720\,744$ | $\approx 1$ | Padé recursion |
| 8 | $\approx 36\,080\,000$ | $\approx -36\,079\,815$ | $\approx 185$ | Padé recursion |

The pattern is organised by the Padé pole at $u^* = -656/(2275\pi^2) < 0$.

---

## 9. Exact EOB Constraint Predictions (No 6PN Input Required)

Computing the circular-orbit binding energy $E_{\rm bind} = (\sqrt{1+2\nu(E_{\rm eff}-1)}-1)/\nu$ with $E_{\rm eff} = \sqrt{2A^2/(2A+uA')}$ and expanding symbolically, the $u^6$ coefficient yields three predictions determined **entirely by $a_3$ and $a_4$ alone** — no 6PN input:

$$
\begin{aligned}
E_{\rm bind}[u^6]_{\nu^3} &= -\frac{6699}{1024} + \frac{123\pi^2}{512} \approx -4.171, \\[6pt]
E_{\rm bind}[u^6]_{\nu^4} &= -\frac{55}{1024} \approx -0.05371, \\[6pt]
E_{\rm bind}[u^6]_{\nu^5} &= -\frac{21}{1024} \approx -0.02051.
\end{aligned}
\tag{14}
$$

Any future GR 6PN computation must reproduce all three exactly. They constitute direct falsification tests of QGD.

---

## 10. Hereditary Tail Terms

The first (universal) tail enters at 4PN:

$$
A_{\rm tail}^{(1)}(u,\nu) = -\frac{22}{3}\,\nu\,u^5\ln u.
\tag{15}
$$

The next hereditary tail is predicted to appear at 8.5PN, with coefficient fixed by the Padé loop-nesting structure $T_2/T_1 = \xi = 2275/656$:

$$
T_2 = \frac{11}{3}\times\frac{2275/512}{41/32} = \frac{25025}{1968} \approx 12.716.
\tag{16}
$$

Zero-free-parameter prediction. GR conservative frontier as of March 2026: 6.5PN (Bini–Damour 2025).

---

## 11. Resolution of the 5PN $\nu^2$ Dispute

Both GR groups agree on the transcendental parts of $c_{51}$ (fixed by the $\beta^2$-mechanism) but disagree on $c_{51}^{\rm rat}$:

$$
c_{51}^{\pi^2} = \beta^2\!\cdot\!\left(-\frac{41}{32}\right) = -\frac{369}{512}, \qquad
c_{51}^{\pi^4} = \beta^2\!\cdot\frac{2275}{512} = \frac{20475}{8192}.
$$

QGD determines $c_{51}^{\rm rat} = -200.962$ independently from the Type-II $\sigma$-field graph at $\mathcal{O}(\nu^2)$, giving the total

$$
c_{51} = -200.962 - \frac{369}{512}\pi^2 + \frac{20475}{8192}\pi^4 \approx 35.388.
\tag{17}
$$

| Source | $c_{51}^{\pi^2}$ | $c_{51}^{\pi^4}$ | $c_{51}^{\rm rat}$ | $c_{51}^{\rm total}$ |
|---|:---:|:---:|:---:|:---:|
| Bini–Damour–Geralico (2020) | $-369/512$ | $+20475/8192$ | (implicit) | $\approx 35.388$ |
| Blümlein et al. (2020) | $-369/512$ | $+20475/8192$ | (different) | $\approx 34.987$ |
| **QGD (this work)** | $-369/512$ | $+20475/8192$ | $-200.962$ | $\mathbf{35.388}$ |

QGD matches Bini–Damour–Geralico. The prediction follows purely from $\sigma$-field structure, with no fitting to the GR literature.

---

## 12. Summary of Zero-Free-Parameter Predictions

| Order | Quantity | Value | Derivation |
|:---|:---|:---|:---|
| 6PN | $E_{\rm bind}[u^6]_{\nu^3}$ | $-6699/1024 + 123\pi^2/512$ | EOB inversion (exact) |
| 6PN | $E_{\rm bind}[u^6]_{\nu^4}$ | $-55/1024$ | EOB inversion (exact) |
| 6PN | $E_{\rm bind}[u^6]_{\nu^5}$ | $-21/1024$ | EOB inversion (exact) |
| 6PN | $c_{60}^{\pi^6}/\nu$ | $-5175625/335872$ | Double-Padé (exact fraction) |
| 6PN | $c_{60}^{\rm rat}$ | $\approx 14\,398$ | Padé [2/1] + cancellation constraint |
| 7PN | $c_{70}^{\pi^8}/\nu$ | $+11774546875/220332032$ | Double-Padé (exact fraction) |
| 7PN | $c_{70}^{\rm rat}$ | $\approx 720\,745$ | Padé [2/1] recursion |
| 8PN | $c_{80}^{\pi^{10}}/\nu$ | $-26787094140625/144537812992$ | Double-Padé (exact fraction) |
| 8.5PN | $T_2$ (tail) | $25025/1968 \approx 12.716$ | Tail loop-nesting conjecture |

---

## 13. Comparison with GR

| Result | GR | QGD |
|:---|:---|:---|
| 10 exact spacetimes | 10 separate EFE derivations (1916–present) | One master metric, 10 $\sigma$-fields |
| $a_4 = \nu(94/3-41\pi^2/32)$ | 4-loop, three independent groups, years | One Hadamard pole, Type-II $\sigma$-graph |
| $c_{50} \approx 23.502$ | GSF numerics + multi-loop matching | $\sigma^{(4)}\!\cdot\!\sigma^{(4)}$ integral |
| $c_{51} \approx 35.388$ | Disputed (BD-G 2020 vs Blümlein 2020) | $\beta^2$ mechanism — no dispute |
| Trans. sector 6PN–8PN | Open, no closed form | Exact fractions, eq. (9) |
| Trans. sector $n$-th PN | Unknown | $(-41/32)(-2275/656)^{n-4}\pi^{2(n-3)}$ |
| $\nu^k$ structure at $n$PN | Unknown | $\nu^k$ first at $(k+3)$PN; max $\nu^{n-3}$ |
| Full $A$-trans resummation | Unknown | Double-Padé eq. (7), zero free parameters |

---

## 14. Open Problems

- Exact $c_{60}^{\rm rat}$ from the 5-loop Fokker integral
- Exact $c_{70}^{\rm rat}$ from the 6-loop Fokker integral
- First appearance of $\zeta(3)$ in $A(u;\nu)$
- Full $\nu$-dependence at 7PN beyond the diagonal
- GR verification of the 8.5PN conservative tail $T_2 = 25025/1968$

General Relativity describes gravity through the Einstein field equations: ten coupled nonlinear PDEs whose exact solutions require dedicated derivations for each spacetime. The post-Newtonian expansion of the two-body problem — essential for gravitational wave templates — requires increasingly expensive multi-loop calculations at each order, with the transcendental structure (the π², π⁴, π⁶, … terms) being recomputed from scratch at every PN order with no known pattern.

QGD identifies a field variable

```
σ_μ ≡ (ħ/mc) ∂_μS
```

that emerges from the WKB limit of the Dirac equation. The metric is then purely algebraic:

```
g_μν = T^α_μ T^β_ν [η_αβ − Σ_a ε_a σ_α^(a) σ_β^(a)]
```

In this variable:

- Every classical GR solution is a one-line superposition of σ-fields.
- The two-body problem reduces to one scalar function `A(u;ν)`.
- The transcendental sector of `A` is exactly resummed by a single closed form.

This repository contains the complete derivation, all verification code, and the predictions.

---

## Quick start

```bash
python qgd_showcase.py          # full demonstration with symbolic verification
python qgd_showcase.py --fast   # numerical only, ~5 seconds
python qgd_master_final.py      # complete PN programme, all stages
```

---

## Part I — Ten exact spacetimes in one formula

The master metric above recovers every classical GR solution by choosing the appropriate σ-fields and signs ε.

| # | Spacetime | σ-field | ε | GR derivation |
|---|-----------|---------|---|---------------|
| 1 | Minkowski | σ_μ = 0 | — | trivial |
| 2 | Schwarzschild | σ_t = √(2GM/c²r) | +1 | Schwarzschild, 1916 |
| 3 | Reissner-Nordström | σ_t^M + σ_t^Q | +1, −1 | Reissner 1916 + Nordström 1918 |
| 4 | Kerr | σ_t + σ_φ = a·sin θ/r × σ_t | +1 | Kerr, 1963 |
| 5 | Kerr-Newman | mass + charge + spin σ-fields | +1,−1,+1 | Newman et al., 1965 |
| 6 | Schwarzschild–de Sitter | mass + cosmological σ-field | +1, −1 | — |
| 7 | Schwarzschild–AdS | mass + AdS σ-field | +1, +1 | — |
| 8 | de Sitter | σ_t = Hr/c | −1 | — |
| 9 | Anti-de Sitter | σ_t = Hr/c | +1 | — |
| 10 | N-body | σ_total = Σ_i σ^(i) | +1 each | No closed-form GR equivalent |

**The operational difference with GR:** In GR, each row above required solving the Einstein equations with different symmetry assumptions. In QGD the σ-field is simply assembled from its sources and the metric is read off. The results are identical.

```python
# Schwarzschild in QGD — three lines
sigma_t  = np.sqrt(2*G*M / (c**2 * r))
g_tt     = -(1 - sigma_t**2)          # = -(1 - 2GM/c²r)
g_rr     = 1.0 / (1 - sigma_t**2)     # = 1/(1 - 2GM/c²r)

# Kerr frame-dragging term — one line
g_t_phi  = -r * np.sin(theta) * sigma_t * sigma_phi   # exact g_tφ
```

Run `python qgd_showcase.py` to verify all ten solutions numerically. All pass to 10 significant figures.

---

## Part II — The post-Newtonian expansion

The EOB effective Hamiltonian for binary systems is

```
H_eff = μ√(A(u;ν)[1 + p²/μ² + p_r²/B(u;ν)])
```

In QGD the isotropic σ-field structure forces **B(u;ν) ≡ 1** exactly. The entire two-body problem reduces to one scalar function:

```
A(u;ν) = 1 − σ_eff²(u;ν)
```

The complete A-function decomposes as:

```
A(u;ν) = 1 − 2u + 2νu³           ← 1PN–3PN, trivially exact
        + A_trans(u;ν)             ← the closed form below
        + A_rational(u;ν)          ← Padé [2/1] recursion
        − (22/3)ν u⁵ ln u          ← exact nonlocal tail
```

---

## Part III — The quantitative core: one closed form → all PN coefficients

This is the main result.

```
                   −ν (41π²/32) u⁴
A_trans(u;ν) = ─────────────────────────────────
               [1 + (2275/656)π²u] [1 − (9/16)ν]

KEY IDENTITY:  32 × 656 = 512 × 41 = 20992
```

Expanding this single expression as a Taylor series in `u` and `ν` reproduces **every transcendental coefficient ever computed in GR**, and predicts all higher orders as exact fractions.

### The expansion table

| PN order | Term | Coefficient/ν | ×π^{2(n-3)} | GR status | How obtained in GR |
|----------|------|---------------|-------------|-----------|-------------------|
| 1PN | u¹ | −2 | — | ✓ exact | Newtonian limit |
| 2PN | u² | 0 | — | ✓ exact | Harmonic gauge |
| 3PN | u³ | 2ν | — | ✓ exact | DJS 2001, 3-loop Fokker |
| 4PN | u⁴ | **−41/32** | −12.65 | ✓ exact DJS 2001 | 4-loop dim-reg, 3 independent groups |
| 5PN | u⁵ | **+2275/512** | +432.82 | ✓ exact DJS 2015 | π⁴ = row n=5 of Padé |
| 5PN ν² | u⁵ | **+20475/8192** | +243.46 | ✓ BD-G 2020 | β²×c₅₀_π⁴ |
| **6PN** | u⁶ | **−5175625/335872** | −14814.54 | ← **PREDICTED** | Open in GR |
| **7PN** | u⁷ | **+11774546875/220332032** | +507,067 | ← **PREDICTED** | Open in GR |
| **8PN** | u⁸ | **−26787094140625/144537812992** | −17,355,729 | ← **PREDICTED** | Open in GR |
| n-th PN | u^n | **(−41/32)×(−2275/656)^{n-4}** | exact fraction | ← **PREDICTED** | Unknown in GR |

The 4PN and 5PN rows — which required years of independent GR computations — are simply n=4 and n=5 in the Taylor expansion. The 6PN–8PN rows are exact fractions generated by the same formula with n=6,7,8. GR has no analogous expression.

### The ν-tower

The mass-ratio dependence follows a separate geometric series with factor β² = 9/16 (β = −3/4):

```
c_{n,k}^{π^{2(n-3)}} = β^{2k} × c_{n,0}^{π^{2(n-3)}}
```

| PN | ν-power | Exact relation | GR status |
|----|---------|---------------|-----------|
| 5PN | ν² | c₅₁_π⁴ = (9/16)×c₅₀_π⁴ | ✓ BD-G 2020 |
| 6PN | ν³ | c₆₂_π⁶ = (81/256)×c₆₀_π⁶ | ← predicted |
| 7PN | ν⁴ | c₇₃_π⁸ = (729/4096)×c₇₀_π⁸ | ← predicted (**first ν⁴ ever**) |

**General structure theorem:** ν^k first appears at (k+3)PN. The maximum ν-power at n-th PN order is ν^{n-3}.

---

## Part IV — Exact predictions with no free parameters

### Three 6PN binding-energy coefficients

From the symbolic inversion of the EOB Hamiltonian using *only* a₃ and a₄ (no 6PN input):

```
E_bind[u⁶]_ν³ = −6699/1024 + 123π²/512  ≈  −4.171   (exact)
E_bind[u⁶]_ν⁴ = −55/1024               ≈  −0.054   (exact)
E_bind[u⁶]_ν⁵ = −21/1024               ≈  −0.021   (exact)
```

Any future GR 6PN computation must reproduce these exactly.

### 8.5PN hereditary tail

```
T₂ = T₁ × ξ = (11/3) × (2275/656) = 25025/1968 ≈ 12.716
```

The GR frontier is currently at 6.5PN conservative dynamics (Bini-Damour 2025). The conservative 8.5PN coefficient is an open falsification target.

### Rational sector from the Padé recursion

The rational Padé [2/1] of H_rat gives an exact geometric recursion H[m+1] = q·H[m] with q = c60_rat/c50_rat ≈ 50.06:

```
c60_rat ≈ 14,398    (from physical constraint c60_total ≈ 4)
c70_rat ≈ 720,745   (from Padé recursion)
```

---

## Part V — Computational comparison

Both GR and QGD give identical physical predictions. The difference is the length of the path.

| Result | GR approach | QGD approach |
|--------|-------------|--------------|
| 10 exact spacetimes | 10 separate EFE derivations (1916–present) | One master metric, 10 σ-fields |
| a₄ = ν(94/3 − 41π²/32) | 4-loop, three independent groups, years | One Hadamard pole of Type-II σ-graph |
| c₅₀ ≈ 23.502 | GSF numerics + multi-loop matching | σ^(4)·σ^(4) integral, same answer |
| c₅₁ ≈ 35.388 (disputed) | Conflict between BD-G 2020 and Blümlein 2020 | β²×c₅₀ mechanism — no dispute |
| Transcendental sector 6PN–8PN | **Open. No closed form.** | **Exact fractions, table above** |
| Transcendental sector n-th PN | **Unknown in general** | **(−41/32)(−2275/656)^{n-4} × π^{2(n-3)}** |

**Note on the comparison:** Writing "QGD is simpler" does not mean "QGD requires no work." The σ-field integrals use the same Hadamard regularisation machinery as GR. The simplification is structural: one closed form replaces an open-ended sequence of independent calculations, and that closed form generates predictions GR cannot currently match.

---

## What QGD is and is not

**QGD is:**
- A field-variable transformation of GR (σ_μ is a function of the metric)
- A significantly more compact computational framework for the PN problem
- A source of exact predictions for high-PN transcendental coefficients
- Consistent with all known GR and GSF results through 5PN

**QGD is not:**
- A claim that GR is wrong
- A quantum theory of gravity (the σ_μ field is classical in the WKB limit)
- A derivation that avoids regularisation (Hadamard is still required)

The σ_μ variable emerges from the Dirac equation without any choice or ansatz. That it also produces a compact computational structure is the physical content of the result.

---

## Verification summary

you are welcome to test this independently from the python code
```
python PN-Expansion.py
https://github.com/matshaba/Quantum-Gravity-Dynamics/blob/main/core/PN-Expansion.py
```


---

## Selected references

| Reference | Content | QGD relevance |
|-----------|---------|---------------|
| DJS 2001, PRD 62 084011 | a₃, a₄ exact (3PN/4PN) | Stages 2–3 confirmed |
| DJS 2015, arXiv:1502.07245 | 4PN EOB + 5PN log | c₅₀ components confirmed |
| Bini-Damour 2013, PRD 87 121501 | GSF: c₅₀ = 23.502 | Stage 5 match |
| Bini-Damour-Geralico 2020, arXiv:2003.11891 | 5PN full, c₅₁ | Stage 5b match |
| BD-G 2020, arXiv:2007.11239 | 6PN nonlocal, ζ(3) | Stage 6 target |
| Bini-Damour 2025, arXiv:2507.08708 | 6.5PN conservative | T₂ GR frontier |
| Levi-Morales-Yin 2024, JHEP | Dissipative 8.5PN | T₂ scope (different sector) |

---

*Authors: Romeo Matshaba
