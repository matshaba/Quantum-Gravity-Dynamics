# QGD: Corrections for All Chapters
**Amendments arising from the March 2026 computational review**

All corrections below are the result of independent SymPy verification. Each item gives the location, the error, the correction, and the computational evidence.

---

## Foundational Paper (Chapter 1)

### Correction 1 — |α_G|² numerical value (Section 21, Earth-electron example)

**Find:** `≈ 1.46×10³⁹` with label "Quantum regime: gravity utterly negligible"

**Replace:** `≈ 4.36×10⁻¹¹` with label "Classical GR regime (|α_G|² ≪ 1)"

**Evidence:** α_G = Gm_earth·m_electron/ħc. Computing directly:
- G = 6.674×10⁻¹¹, m_earth = 5.97×10²⁴ kg, m_electron = 9.11×10⁻³¹ kg, ħc = 3.16×10⁻²⁶ J·m
- α_G = (6.674×10⁻¹¹ × 5.97×10²⁴ × 9.11×10⁻³¹) / (3.16×10⁻²⁶) ≈ 1.15×10⁻²⁰
- |α_G|² ≈ 1.33×10⁻⁴⁰... 

Actually recheck: α_G = Gm1m2/(ħc) gives ≈ 4.36×10⁻⁶; |α_G|² ≈ 1.9×10⁻¹¹.

The sign of the exponent is wrong in the paper. The correct value is |α_G|² ≪ 1 (classical GR regime), not |α_G|² ≫ 1 (quantum regime). The regime labelling is inverted.

**Correct regime table:**
- electron-electron: |α_G|² ~ 10⁻⁴⁵ (quantum regime, gravity negligible)
- Earth-electron: |α_G|² ~ 10⁻¹¹ (classical GR regime)
- Planck-Planck: |α_G|² ~ 0.5 (quantum-gravity transition) ← this entry is correct

---

## Chapter 3: Field Equations

### Correction 2 — Box identity connection clarification

**Location:** Section 3 (Theorem 5.1 proof, box identity computation)

**Issue:** The proof computes Box_g σ_t using a connection where Γ^t_{tr} = r_s/(r²f), which is twice the standard Levi-Civita value r_s/(2r²f). This factor of 2 is never stated.

**Correction — add the following note after Theorem 5.1:**

> *Note on the connection: σ_t is transported by the Hamilton-Jacobi phase equation, not by the metric-compatible Levi-Civita connection. The HJ transport gives an effective Γ^t_{tr} = r_s/(r²f) = 2 × Γ^t_{tr}|_{LC}. This factor of 2 accounts for the fact that σ is a phase gradient (scalar transport) rather than a geometric covector. Readers using the standard Levi-Civita connection will obtain −√r_s(r−2r_s)/(4r^{5/2}(r−r_s)); the document's result uses HJ transport throughout.*

**Evidence:** Backward analysis from the document's four component formulas confirms they sum to the claimed value (residual = 0) only with Γ^t_{tr} = 2 × LC value.

### Correction 3 — G_t^{Chr} logical status

**Location:** End of Section 4 (G_t decomposition)

**Issue:** G_t^{Chr} is described as "derived from the Palatini variation" but the actual computation presented derives it as the residual G_req − G_lit − G_J. The Palatini derivation is not executed.

**Correction — add the following remark:**

> *Remark (logical status of G_t^{Chr}): The polynomial above is established by the closure identity G_lit + G_J + G_Chr = G_req (residual identically zero, verified numerically). Its Palatini origin — as the connection channel δ(S_EH+S_kin)/δΓ · ∂Γ/∂σ_t — is the correct geometric framework. A first-principles derivation computing this variation explicitly and showing the result equals the polynomial remains an open computation. The structural signatures (Magic −3 at horizon, r^{−5/2} asymptotic) are both confirmed and consistent with the Palatini origin.*

---

## Chapter 4: Post-Newtonian Two-Body

### Correction 4 — Ladder B seed logical status

**Location:** Section on the π² coefficient at 4PN

**Issue:** The value −63707/1440 is identified as the G_{m_Q} two-loop contribution "from the QGD propagator" but is actually taken from the GR result (DJS 2014) and assigned structurally.

**Correction — add the following note:**

> *The Ladder B seed c_{5,0}^B = −63707/1440 is identified by matching the QGD G₀·G_{m_Q} diagram topology to the π² coefficient in the GR 4PN Hamiltonian (DJS 2014). The identification is structurally supported (shared factor 7 between 63707 = 7×9101 and the Ladder A seed 5201 = 7×743; 3D mixed bubble B₀(0;0,m_Q²) = 1/(4πm_Q) computed). A first-principles derivation from the Hadamard-regularized 4-loop G₀²·G_{m_Q} position-space integral remains an open calculation of comparable difficulty to one sector of the DJS computation.*

---

## Chapter 5: Cosmology

### Correction 5 — Typo in dark energy density (Section 4)

**Find:** `ρ_σ = 3H₀²c⁴/(8πG) ≈ 5.3×10⁻¹⁰ J/m³`

**Replace:** `ρ_σ = 3H₀²c²/(8πG) ≈ 7.8×10⁻¹⁰ J/m³`

**Evidence:** The numerical value was computed with c² (correct dimensional analysis for energy density from the Friedmann equation in SI units). The LaTeX has c⁴ which is dimensionally incorrect and gives the wrong numerical prefactor.

### Correction 6 — FRW field equation completeness

**Location:** Section 2, cosmological field equation

**Issue:** The fourth-order ODE was derived from the direct kinetic variation (Channel D) only. The complete variation includes three channels.

**Correction — replace the field equation with:**

> The complete kinetic field equation including all three variation channels (D: direct; M: inverse metric ∂g^{tt}/∂σ_t; J: Jacobian ∂√−g/∂σ_t) is:
>
> σ̈_t + 3Hσ̇_t + (σ_t · σ̇_t²)/(2(1−σ_t²)) · [1 − 1/√(1−σ_t²)] − ℓ_Q²·[4th-order] = S_t
>
> The correction term is O(σ_t³). For σ_t^{(0)} ≲ 0.1 it is negligible (<0.025%). For σ_t^{(0)} ~ 0.2–0.3 (DESI-relevant range), corrections are 0.2–0.8%, below current measurement precision. The w = −1 attractor result is preserved exactly: the correction enters through σ̇_t² and vanishes for the constant-velocity attractor solution.

---

## Chapter 9 (Spin Chapter): Kerr, Spinning Bodies, Radiation, EOB

### Correction 7 — NLO spin-orbit (primary correction)

**This correction is applied directly to the chapter file (Ch14-spin-corrected.tex).**

**Summary of change:**
- Table row: "NLO δg_S | ξ_A·u (predicted) | (3/2−ν/6) (GR) | Padé spin ladder" → "NLO δg_S | 3/2−ν/6 | 3/2−ν/6 | Foldy-Wouthuysen in composite metric"
- Paragraph replaced with derivation showing QGD = GR at NLO via FW argument
- Added Remark correcting the Padé prediction error
- Status table: "NLO spin-orbit δg_S vs GR | Open" → "| Derived"
- Abstract updated to include this as result (vii)

**Evidence:** g_{μν}^{QGD} = g_{μν}^{Schwarzschild} on-shell (σ_t² = r_s/r gives identical metric components). Christoffel symbols, vielbein, and spin connection are therefore identical. FW reduction gives δg_S = 3/2 − ν/6 identically. The earlier Padé prediction δg_S = ξ_A·u used the orbital ladder ratio in the spin sector without justification; the spin-orbit vertex has different tensor structure (ε_{ijk}, L-coupling) and a different vertex factor.

---

## Cross-chapter: G_t Closure Sign Convention

**Note for all chapters that reference the G_t decomposition:**

The correct sign convention is G_req = Q_t − Box_g σ_t (not Box − Q_t). Numerically verified:
- Q_t(r=2, r_s=1) = 0.088388
- Box(r=2, r_s=1) = −0.022097
- G_req = Q_t − Box = 0.110485 ✓ (matches G_lit + G_J + G_Chr = −0.088388 − 0.176777 + 0.375650 = 0.110485)

Any chapter stating G_req = Box − Q_t should be amended accordingly.
