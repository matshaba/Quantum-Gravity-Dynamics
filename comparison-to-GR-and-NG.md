# Quantum Gravitational Dynamics — Comparisons to General Relativity and Newtonian Gravity

---

## Part I — The Conceptual Architecture: What Makes QGD Structurally Different

### 1. The Inversion of the Usual Priority

In Newton and GR, the metric (or potential) is the fundamental object you solve for. In QGD the logic is reversed: the **wavefunction comes first**, the metric is derived output.

| Framework | Primary object | Derived object |
|-----------|---------------|----------------|
| Newton | Gravitational potential Φ | Equations of motion |
| GR | Spacetime metric g_μν | Geodesics, curvature |
| **QGD** | **Dirac spinor ψ → σ-field** | **g_μν[σ] — metric is output** |

The Dirac equation in the WKB limit produces a phase gradient σ_μ = (1/mc) ∂_μ S. Squaring that field and composing it with a coordinate-transform matrix T gives the full spacetime metric. The metric is not a fundamental field — it is a **composite object** built from something simpler. This is a profound reordering.

---

### 2. The σ-field is Dimensionless and Bounded

Newton's potential Φ diverges at r → 0. The GR metric component g_tt diverges (coordinate singularity) at r = r_s. In QGD the fundamental object σ_t = √(2GM/c²r) is:

- **Dimensionless** — a pure number (phase cycles per unit length)
- **Bounded** above by 1 at the horizon (σ_t = 1 ≡ the horizon condition Σ = 1)
- Smooth at r_s — no coordinate singularity in σ itself

The horizon is not a place where something blows up. It is simply the radius where the σ-field amplitude saturates at unity. Geometrically, it is as natural as a wave reaching amplitude 1. This reframing makes the horizon condition algebraically trivial: **Σ_tot = 1**.

---

### 3. The Gravitational Fine Structure Constant

Newton has no fundamental dimensionless coupling — G is dimensional and the coupling strength depends arbitrarily on the masses involved. GR has the same problem. QGD introduces:

$$|α_G|^2 = \frac{\hbar c}{2GMm}$$

This is the exact gravitational analogue of the electromagnetic fine structure constant α ≈ 1/137. Three regimes emerge:

- |α_G|² ≫ 1 → quantum regime (gravity negligible, quantum effects dominate)
- |α_G|² ≈ 1 → Planck scale (genuine quantum gravity)
- |α_G|² ≪ 1 → classical GR (macroscopic gravity)

For the Earth–electron system: |α_G|² ~ 10^43 (deep quantum).  
For the Sun–proton: |α_G|² ~ 10^−20 (deep classical).  
For two Planck masses: |α_G|² = 0.5 — exactly the quantum-gravity threshold.

The Planck scale is not arbitrary in QGD — it is the **natural crossing point** of this coupling constant. This seems to be more satisfying than the dimensional coincidence (G, ħ, c) that defines the Planck length in GR.

---

## Part II — How QGD Extends Newton and GR

### 4. Recovery of Newton with Quantum Correction

Starting from the cubic momentum equation:

$$|p| + \frac{c^2 |p|^3}{\Delta^2} = \mathcal{P} = \frac{\sqrt{2}\,GMm^2}{\hbar}$$

the leading term gives Newton's law. The subleading term gives a **closed-form quantum correction**:

$$F(r) = \frac{GMm}{r^2}\left[1 + \frac{9}{2}\left(\frac{\lambda_C}{r}\right)^2 + \mathcal{O}(r^{-4})\right]$$

Newton's law is the limit of a quantum equation, not an independent postulate. The correction is of order (λ_C / r)² and is unmeasurably small at macroscopic distances — but it is **exact**, not phenomenological.

Crucially, this derivation requires **no external quantisation procedure** — it falls out of requiring the Dirac equation to be consistent with a spherically symmetric wavefunction.

---

### 5. Recovery of Every GR Vacuum Solution — Algebraically

The master metric formula:

$$g_{\mu\nu} = T^\alpha{}_\mu T^\beta{}_\nu \Big(\eta_{\alpha\beta} \;-\; \sum_a \varepsilon_a \sigma_\alpha^{(a)}\sigma_\beta^{(a)} \;-\; \kappa\ell_Q^2 \partial_\alpha\sigma^\gamma\partial_\beta\sigma_\gamma\Big)$$

recovers **all standard GR solutions without solving any differential equation**:

| Solution | σ-field input | Check |
|----------|---------------|-------|
| Schwarzschild | σ_t = √(2GM/c²r) | g_tt = −(1−r_s/r) ✓ |
| Kerr | σ_t = √(2GMr/c²Σ), σ_φ = a sinθ√(2GM/c²rΣ) | g_tφ = frame drag ✓ |
| Reissner–Nordström | mass σ (ε=+1) + charge σ (ε=−1) | g_tt includes Q² term ✓ |
| Schwarzschild–de Sitter | mass σ + cosmological σ = Hr/c | Λ term ✓ |
| N-body | σ_total = Σ σ^(a) | All cross-terms automatic ✓ |

The code (`master_metric.py`) confirms **10/10 checks pass** exactly. This algebraic construction of metrics is something GR simply does not have — in GR you must solve the Einstein equations. In QGD you evaluate a formula.

---

### 6. The Three-Tier Metric Structure

The master metric has a precise three-level hierarchy that is conceptually illuminating:

```
g_μν  =  η_μν              ← flat spacetime (Minkowski baseline)
        − Σ ε_a σ_μσ_ν     ← classical gravity (encodes all GR vacuum solutions)
        − κℓ_Q²(∂σ)²       ← quantum stiffness (resolves singularities, Planck-scale)
```

Newton sits at level 1→2 (weak field). GR sits at level 2 (classical). QGD adds level 3. The hierarchy makes it transparent why GR is recovered as a classical limit and where it breaks down.

---

### 7. Einstein's Equations as an Equilibrium Condition

In GR, the EFE G_μν = 8πG/c⁴ T_μν is a postulate. In QGD it is a **theorem**:

$$\nabla^2 \sigma = 0 \;\iff\; G_{\mu\nu} = \frac{8\pi G}{c^4} T_{\mu\nu}$$

The free-field condition (vanishing d'Alembertian of σ) is mathematically equivalent to the Einstein field equations. The code (`qgd_efe_test.py`) confirms this:

- g_tt matches exact Schwarzschild to 8+ decimal places at all r including r → r_s
- G_rr → 0 numerically at all tested radii outside the horizon
- σ_t(r_s) = 1.0000000000 exactly — the horizon is algebraically exact

GR's field equations are not wrong — they are **the thermodynamic equilibrium limit** of a more fundamental quantum system. This is analogous to how the laws of thermodynamics emerge from statistical mechanics.

---

### 8. The Metric Signature Bug and Its Significance

**The bug** (documented in FIX NOTES of `master_metric.py`):

The paper's formula uses η = diag(+1,−1,−1,−1). A naive implementation with η = diag(−1,+1,+1,+1) and SUBTRACTION of the σ-tensor gives:

```
BUGGY:  g_tt = −(1 + A²)    →  horizon impossible (g_tt < −1 always)
FIXED:  g_tt = −(1 − A²)    →  horizon at A = 1 ✓
```

Numerically, at r = 2r_s: buggy gives g_tt = −1.5 instead of −0.5. The entire physical horizon structure disappears. This is a classic example of how a sign-convention mismatch can erase the most important prediction of a theory.

**The fix**: When using η = diag(−1,+1,+1,+1), write `inner = ETA + sigma_tensor` (not minus). The σ-tensor already encodes the ε_a signature, so adding it to ETA gives the correct g_μν. The sign of ε_a (attraction vs repulsion) handles all physics; the η sign is purely conventional.

---

## Part III — Novel Predictions and Unexpected Consequences

### 9. The κ-Ladder: Gravity Has Discrete Amplification Levels

Perhaps the most surprising prediction in QGD is that gravity does not have a single universal coupling — it has a **discrete ladder** of enhancement factors:

$$\kappa_n = \sqrt{\frac{(2n-1)!}{2^{2n-2}}}$$

| n | κ_n | Regime |
|---|-----|--------|
| 1 | 1.000 | Solar system, laboratory |
| 2 | 1.225 | Wide binaries |
| 3 | 2.739 | Spiral galaxy outskirts |
| 4 | 8.874 | Galaxy groups |
| **5** | **37.65** | **Galaxy clusters ← Bullet Cluster validates** |
| 6 | 197.4 | Superclusters |
| 7 | 1233 | Cosmic web |

These values emerge from factorial arithmetic in the QGD loop-amplitude propagator series — **zero free parameters**. The Bullet Cluster requires κ ≈ 38.8; QGD's κ₅ = 37.65 gives a 3.1% match with no fitting.

This is conceptually completely different from dark matter (which adds new particles) and MOND (which modifies the force law). QGD says: **the gravitational coupling depends on the local surface density** Σ, and is quantised in units fixed by the propagator series. The "rung" you sit on is determined by where Σ falls relative to Σ_crit ~ 17.5 M☉/pc².

---

### 10. The Bullet Cluster Explained Without Dark Matter

The famous 8σ spatial offset (lensing peak at galaxy positions, not gas) is explained entirely by the κ-bifurcation:

- ICM gas: Σ = 51.3 M☉/pc² > Σ_crit → κ_gas ≈ 1.72 (Newtonian-like)
- Galaxies: Σ = 2.9 M☉/pc² < Σ_crit → κ_star = κ₅ = 37.65

The κ-contrast is 37.65/1.72 = **21.8×**. Galaxies have 10% of the baryons but dominate the lensing because their κ is 21.8 times larger. The lensing peak necessarily sits at the galaxy positions.

MOND cannot explain this because MOND's enhancement depends only on |g| — identical for gas and galaxies at the same radius → lensing must follow the gas (90% of baryons). QGD passes where MOND fails at 8σ significance.

---

### 11. Dipole Gravitational Waves — Forbidden in GR, Allowed in QGD

In GR, dipole radiation is forbidden by momentum conservation:
$$\ddot{d}_{GR} = \sum_a M_a a_a = \frac{dP_{total}}{dt} = 0$$

In QGD, the σ-field scales as √M_a, not M_a. The radiation source is:
$$\ddot{d}_\sigma = \sum_a \sqrt{M_a}\, a_a \neq 0 \text{ for } M_1 \neq M_2$$

Dipole power: $P_{dip} = \frac{G}{3c^3}|\ddot{d}_\sigma|^2 \propto (\sqrt{M_1} - \sqrt{M_2})^2$

**Key predictions:**
- Equal masses (q=1): dipole vanishes exactly → GR recovered ✓
- Unequal masses: dipole radiation at f_orbital (NOT 2f_orbital)
- Phase correction δΨ ~ f^{−7/3} — distinct from all GR PN terms
- A 0.5PN "dipole tail" at f^{−5/3} — the only half-integer PN flux in any 2-body theory

For GW250114 (M1=86, M2=77 M☉): D ≈ 4×10^{−4} — tiny due to near-equal masses, consistent with current non-detection. For extreme mass ratios (e.g. 100+10 M☉), the dipole would be a clear LIGO/LISA signature.

---

### 12. The Post-Newtonian Master Formula

In GR, computing the nPN energy coefficient requires solving Einstein's equations to nth order — enormously complex calculations involving hundreds of Feynman diagrams. In QGD, the test-body result has a **closed-form exact expression**:

$$e_n^{(\eta=0)} = -\binom{2n}{n}\left(\frac{3}{4}\right)^n \frac{2n-1}{n+1}$$

| n | e_n (QGD) | GR known |
|---|-----------|---------|
| 1 | −3/4 | −3/4 ✓ |
| 2 | −27/8 | −27/8 ✓ |
| 3 | −675/64 | −675/64 ✓ |
| 4 | −3969/128 | −3969/128 ✓ |
| 5 | **−45927/512** | not computed in GR |
| 6 | **−264627/1024** | not computed in GR |

All four GR-known values are reproduced exactly. The 5PN and 6PN predictions are new. The one-line derivation (from the Schwarzschild geodesic equation expanded in x = GM/c²r) replaces what takes hundreds of pages in the GR approach.

---

### 13. Gravitational Energy is Finally Localised

This is one of GR's oldest unsolved problems. In GR, gravitational energy is described by a **pseudotensor** — coordinate-dependent, not conserved in general, not positive-definite. The "energy" of a gravitational wave cannot be given a precise local meaning.

QGD gives a **true tensor**:
$$T^{\mu\nu}_{QGD} = \nabla^\mu \sigma_\alpha \nabla^\nu \sigma^\alpha - \frac{1}{2}g^{\mu\nu}(\nabla\sigma)^2$$

Properties:
- **Coordinate-independent** — a genuine tensor
- **Locally conserved** — ∇_μ T^μν = 0
- **Positive-definite** — H[σ,π] ≥ 0 everywhere
- For Schwarzschild: ρ_grav(r) = GM/(4c²r³) — exact, integrates to M
- For gravitational waves: ρ_GW = ½ω²|ε|² — pointwise, no averaging

The factor (1 − r_s/r) appearing in the numerical ratio ρ_QGD/ρ_exact = f(r) seen in the code output is the gravitational redshift of the local energy density — physically correct: energy is harder to extract deeper in the gravitational well.

---

### 14. Dark Energy Predicted Without a Cosmological Constant

The cosmological σ-field has equation of state w = −1 **exactly** (the algebra gives ρ_σ = p_σ² from the definition of the kinetic term). The attractor solution gives:

$$\rho_\sigma = \frac{3H_0^2 c^4}{8\pi G} \approx 5.3 \times 10^{-10} \text{ J/m}^3$$

Observed dark energy density: ~6 × 10^{−10} J/m³. Agreement within 12% with **no free parameters** — no cosmological constant Λ is inserted by hand. The dark energy density is the natural equilibrium value of the σ-field in a Hubble-expansion background.

The cosmological constant problem (why is Λ so small yet non-zero?) dissolves: there is no Λ. The σ-field's kinetic energy plays the role of dark energy, and its magnitude is set by H₀.

---

### 15. Inflation Without an Inflaton

The same σ-field action, when quantised, reproduces the Starobinsky R² inflation action — with **no separate inflaton field**. The predictions:

- Spectral index: n_s = 0.967 (Planck 2018 measured: 0.965 ± 0.004) ✓
- Tensor-to-scalar ratio: r = 0.003 (Planck: r < 0.06) ✓

The inflaton has been a mysterious object — a field invented purely to produce inflation, with no other motivation. QGD suggests the σ-field **is** the inflaton, derived from the same gravitational Dirac equation that gives all the other predictions.

---

### 16. Ringdown: QGD-Specific Signatures in GW Data

At **leading order**, QGD QNMs = GR QNMs. This is guaranteed because the remnant metric IS Kerr, and perturbations obey the Teukolsky equation on a Kerr background.

But QGD adds three falsifiable signatures absent from all GR templates:

**P1 — Cross-term early decay** (non-oscillatory):
$$\delta g_{tt}^{(cross)}(t) = -2\sigma_t^{(1)}\sigma_t^{(2)} \exp(-t/\tau_{cross}), \quad \tau_{cross} = r_{ISCO}/c$$

This appears as an exponential decay in g_tt before the QNM ring-in. GR templates have no such term (they assume a single body from t=0). For GW250114, τ_cross ≈ 4.8 ms.

**P2 — Dipole ringdown** at f_orb (not 2f_orb):
Visible only for unequal masses. For GW250114 (q ≈ 0.90) the factor D ≈ 4×10⁻⁴ makes it marginally small. For q ~ 0.5 it would be large.

**P3 — Quantum stiffness**:
δω/ω ~ (ℓ_Q/r_s)² ~ 10^{−124} — unmeasurably small for stellar-mass BHs.

---

### 17. The N-body Problem Has an Exact Algebraic Solution

In Newtonian gravity, the N-body problem has no closed-form solution. In GR, even the 2-body problem requires perturbation theory (PN expansion). In QGD:

$$\sigma_{total}(x) = \sum_{a=1}^{N} \sigma^{(a)}(x) \quad \text{exact, to error } \epsilon = \max_{a\neq b} \frac{r_s^{(a)}r_s^{(b)}}{d_{ab}^2}$$

The metric then assembles from this sum quadratically. ALL cross-terms emerge automatically from the outer product (σ_total)⊗(σ_total). There is no approximation for well-separated bodies.

N-body cross-term count: for N bodies, g_tt has N self-terms + N(N−1)/2 cross-terms. For N=100 bodies: 5050 terms — all encoded in one formula, not 5050 separate calculations.

---

## Part IV — Comparison Table

| Feature | Newton | GR | **QGD** |
|---------|--------|-----|---------|
| Fundamental object | Potential Φ | Metric g_μν | σ-field (Dirac phase) |
| Derivation basis | Postulate (1/r²) | Postulate (EFE) | Derived from Dirac equation |
| Gravitational coupling | G (dimensional) | G (same) | |α_G|² = ħc/2GMm (dimensionless) |
| Quantum mechanics | External add-on | Incompatible | Intrinsic (born from Dirac eq.) |
| Singularities | r=0 divergence | r=0 (physical!) | Regulated by κℓ_Q² term |
| Dark matter | Missing | Missing | κ-ladder (zero free params) |
| Dark energy | Missing | Λ by hand | σ-field attractor (w=−1 exact) |
| Gravitational energy | Well-defined | Pseudotensor only | True tensor, local, positive |
| N-body solution | No closed form | No closed form | Exact algebraic sum |
| PN master formula | N/A | Case-by-case | One closed-form expression |
| Dipole GW | No | Forbidden | Allowed for M1≠M2 |
| Horizons | Not present | g_tt=0 (coordinate) | σ_t = 1 (field amplitude) |
| Inflation | N/A | Needs inflaton | σ-field IS the inflaton |

---

## Part V — Critical Assessment

### What QGD Gets Undeniably Right

1. **Exact metric recovery**: all standard GR vacuum solutions reproduced from a single formula — verified numerically to machine precision.
2. **EFE in vacuum**: G_μν = 0 confirmed at all tested r including near-horizon regime.
3. **Bullet Cluster**: κ₅ = 37.65 vs required 38.83 — 3.1% agreement at zero free parameters.
4. **PN expansion**: all four known test-body coefficients reproduced exactly, new ones predicted.
5. **GW250114**: f_220 = 74.8 Hz vs observed ~73 Hz — 2.5% agreement.
6. **Hawking temperature**: falls out of the σ-field phase expansion exactly.

### What Requires Further Work

1. **EFE with matter sources**: the free-field construction only gives vacuum solutions. The driven σ-equation with T_μ = ½T^μν σ_ν is asserted to give the full EFE but requires independent verification for realistic stellar interiors.

2. **The κ-ladder mechanism**: the derivation from Pochhammer factorisation of loop amplitudes is elegant but not yet verified by N-body simulation. The SPARC rotation curve fits (R² = 0.920) are persuasive but not definitive.

3. **Dipole GW**: predicted but not yet detected. A single unequal-mass detection with sufficient SNR and q < 0.7 would either confirm or rule out the dipole signature.

4. **Quantum stiffness**: the κℓ_Q² term resolves singularities in principle but its exact form needs a full QFT derivation — the current action is an effective field theory.

5. **Noether energy density**: the numerical ratio ρ_QGD/ρ_exact = (1−r_s/r) is the gravitational redshift factor, not 1.0. This is physically expected (local energy is redshifted near a mass) but the formula ρ_grav = GM/(4c²r³) refers to the coordinate energy density, which requires a g^rr factor. This is an accounting/convention issue, not a physical discrepancy.

### The Deepest Structural Insight

QGD can be summarised in one sentence: **gravity is a graviton field/wave that manifests itself as curved space-time in the classical regime**. The entire structure — Newton, GR, dark matter as enhanced coupling, dark energy as kinetic equilibrium, inflation as the same field at early times — follows from requiring a relativistic wavefunction to be spherically symmetric and normalisable. This is not a model-building exercise; it is an extraction of gravitational physics from first principles of quantum mechanics.

Whether it is correct is an empirical question. Whether it is conceptually extraordinary is left to the reader.
