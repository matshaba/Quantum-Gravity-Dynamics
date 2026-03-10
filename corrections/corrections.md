# Required Fixes for Quantum Gravity Dynamics (QGD) Derivation

## 1. Dimensional Consistency of the Wavefunction Ansatz

**Issue:** The wavefunction
\[
\psi = \frac{\alpha_G}{r}\,u\,e^{iS/\hbar}
\]
must have dimensions \(L^{-3/2}\). Currently \(\alpha_G\) is defined as dimensionless (\(\alpha_G^2 = ic\hbar/(2GMm)\)), so \(|\psi|^2 = |\alpha_G|^2/r^2\) would have dimensions \(L^{-2}\), which is inconsistent.

**Fix:** Redefine \(\alpha_G\) to carry the missing dimension \(L^{-1/2}\). A natural choice is to introduce a fundamental length scale, e.g., the Planck length \(\ell_P = \sqrt{\hbar G/c^3}\), or the Compton wavelength of the test particle \(\lambda_c = \hbar/(mc)\). For example, set
\[
\tilde{\alpha}_G = e^{i\pi/4}\sqrt{\frac{c\hbar}{2GMm}} \times \sqrt{\ell_P},
\]
so that \(\tilde{\alpha}_G\) has dimension \(L^{-1/2}\). Then all expressions become dimensionally consistent. Alternatively, absorb the dimension into the spinor \(u\) by defining \(u\) to include a factor \(\sqrt{mc/\hbar}\) etc. The key is to ensure that \(|\psi|^2\) has the correct dimension.

---

## 2. Proper Expansion for the Phase (Replacing the Invalid Taylor Series)

**Issue:** Expanding \(e^{2imcr/\hbar}\) in powers of \(r\) is invalid for macroscopic distances because the exponent is huge.

**Fix:** Use the cubic equation derived from current conservation (Eq. (103) in the paper) as the starting point. In the non‑relativistic limit, solve it perturbatively in the small parameter \(\epsilon = p^2c^2/\Delta^2 \sim v^2/c^2\) (or \(\hbar/(mcr)\)). The leading term yields Newton's law, and the next term gives post‑Newtonian corrections. This avoids any illegitimate expansion of an oscillatory exponential.

---

## 3. Determination of the Constant \(C\)

**Issue:** The constant \(C\) is currently set by evaluating at the gravitational Bohr radius \(a_0 = \hbar^2/(GMm^2)\) using the classical momentum \(p(a_0) = \sqrt{2}GMm^2/\hbar\). This is circular because the Bohr radius itself depends on \(G\) and the masses, and we are trying to derive gravity.

**Fix:** Derive \(C\) from a genuine boundary condition, such as the normalization of a bound state wavefunction. Solve the effective Schrödinger equation (or the full Dirac equation) for the system and impose square‑integrability; this yields a quantization condition that determines both the energy levels and the normalization constant \(C\) without circularity. The Bohr radius should *emerge* from that quantization, not be used as input.

---

## 4. Coarse‑Graining Must Start from an Interacting Theory

**Issue:** The effective Lagrangian in Section 3 (containing four‑fermi, Maxwell, and \(J^2\) terms) cannot be generated from a free Dirac action by integrating out fast modes.

**Fix:** Begin with a microscopic Lagrangian that includes the necessary interactions. A suitable starting point is
\[
\mathcal{L} = \bar{\psi}(i\gamma^\mu\partial_\mu - m)\psi + G(\bar{\psi}\psi)^2 - \frac{1}{4}F_{\mu\nu}F^{\mu\nu} - e\bar{\psi}\gamma^\mu A_\mu\psi,
\]
as suggested. This theory can then be coarse‑grained using standard Wilsonian methods, and the effective operators (including the \(J^2\) term, which may arise from derivative expansions of the four‑fermi interaction) will appear legitimately. The RG flow and decoupling must be re‑derived with this interacting Lagrangian.
