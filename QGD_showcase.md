# Quantum Gravitational Dynamics (QGD)

**A σ-field representation of GR that reduces the entire two-body post-Newtonian problem to one closed form.**

> *"Not a replacement for GR — GR, viewed through a coordinate that makes the structure visible."*

---

## What is this?

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

```
python PN-Expansion.py
https://github.com/matshaba/Quantum-Gravity-Dynamics/blob/main/core/PN-Expansion.py
```

```
PART A  13/13 metric checks passed
PART B  9/9  PN coefficient checks passed (incl. c51, T2)
PART C  3/3  exact EOB predictions confirmed symbolically
PART D  —    qualitative comparison
PART E  5/5  double-Padé entries verified symbolically
```

All results are independently computable. No fitting to GR data beyond what is stated explicitly in the code.

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
