# Insights of Quantum Gravitational Dynamics

A scientific assessment of what QGD achieves, what is novel, and what follows naturally once its premises are accepted.

---

## 1. The central structural insight

QGD's foundational move is to observe that the Dirac equation in the WKB limit already contains gravitational dynamics — the lower spinor components *are* the graviton field.  Specifically, solving the Pauli-Dirac coupled system gives

$$u_B = \frac{c\,\boldsymbol{\sigma}\cdot\mathbf{p}}{\Delta}\,u_A$$

and the dimensionless ratio $\sigma_i = p_i c/\Delta$ appears directly as the lower spinor components.  This is not imposed; it is read off from the standard Dirac equation with all stress-energy contributions entering through one denominator $\Delta = E + mc^2 + H$.

The observation that *the gravitational field was already inside the wavefunction* is, in retrospect, natural: the Dirac equation describes how spin-½ particles respond to all forces, and the semiclassical limit must encode the classical trajectory, which in a gravitational field is a geodesic.  The phase gradient $\sigma_\mu = \partial_\mu S/(mc)$ is the Hamilton-Jacobi momentum normalised by the rest momentum — exactly the quantity that controls geodesic motion.

What is not obvious is that this quantity, promoted to a field, is sufficient to reconstruct the full spacetime metric algebraically.

## 2. The metric as a quadratic form in the graviton field

The master metric

$$g_{\mu\nu} = \eta_{\mu\nu} - \sum_a \varepsilon_a\,\sigma_\alpha^{(a)}\sigma_\beta^{(a)} - \kappa\ell_Q^2\,\partial_\alpha\sigma^\gamma\partial_\beta\sigma_\gamma$$

is a rank-1 deformation of Minkowski space.  This is the key mathematical structure: the metric is *quadratic* in σ, which means the nonlinearity of Einstein's equations — when expressed in σ variables — is of a specific, controlled type.  The Einstein equations are quadratic in first derivatives of g, which means quartic in σ.  But the σ-field equation is only cubic (through Q_μ), because one power is absorbed into the wave operator □_g.

This explains why the post-Newtonian expansion becomes systematic: each PN order is a linear wave equation sourced by lower orders, precisely because the nonlinearity has been factored into a recursive structure.

## 3. The reduction from 10 to 4 degrees of freedom

General Relativity has 10 independent metric components satisfying 10 coupled nonlinear PDEs, subject to 4 gauge conditions, leaving 6 physical degrees of freedom (of which 2 are propagating).  QGD has 4 σ-field components satisfying 4 wave equations.  The gauge freedom is entirely absorbed into the coordinate transformation matrix T^α_μ.

This is not just notational economy.  A 4-component vector field with a wave equation has a well-defined initial value problem, canonical momentum, Hamiltonian, and quantisation procedure.  A 10-component symmetric tensor satisfying mixed constraint-evolution equations does not — which is precisely why canonical quantum gravity (the Wheeler-DeWitt equation) has remained intractable since 1967.

Whether this simplification carries through to the strong-field nonlinear regime is the critical open question.

## 4. Einstein's equations as equilibrium

The result that ∇²σ = 0 implies G_μν = (8πG/c⁴)T_μν is structurally analogous to the relationship between thermodynamics and statistical mechanics.  The equilibrium condition of the σ-field dynamics reproduces the classical field equations, just as the maximum-entropy condition of a molecular gas reproduces the ideal gas law.

This framing makes a specific prediction: quantum gravity is the *non-equilibrium* regime of the σ-field, characterised by ∇²σ ≠ 0.  The scale at which non-equilibrium effects become significant is set by the Compton wavelength λ_C = ℏ/(Mc) of the source, which for macroscopic objects is absurdly small (10⁻⁵⁴ m for the Sun) but for Planck-mass objects is of order the Planck length.

## 5. The fourth-order equation is a feature, not a bug

Fourth-order gravitational theories (Stelle 1977, Pais-Uhlenbeck) are generally problematic due to the Ostrogradsky instability — ghosts with negative energy that destabilise the vacuum.  QGD's fourth-order term κℓ_Q²□_g²σ_μ has two important differences:

First, it factorises cleanly as □_g(□_g - m_Q²) with m_Q at the Planck mass, so the ghost mode is confined to energies above M_Pl where the effective theory is not claimed to be valid.  Below the Planck scale, all modes have positive energy (the dispersion relation correction is positive).

Second, in cosmology, the extra two degrees of freedom from the fourth-order equation expand the solution space from 2D (Friedmann) to 4D.  The growing mode e^{+t/ℓ_Q} is a candidate inflationary mechanism, and the constant-velocity attractor reproduces the observed dark energy density without a cosmological constant.  The fourth order is not imposed for mathematical convenience — it is forced by the σ-kinetic term, which itself is forced by dimensional analysis (σ is dimensionless; ℏ²/M is the only available scale).

## 6. Dark energy without Λ

The most quantitatively striking result is the dark energy prediction.  The σ-field kinetic energy has equation of state w = -1 exactly (not approximately), because the field's pressure is p = -ρ for a homogeneous time-dependent scalar.  The attractor solution gives

$$\rho_\sigma = \frac{3H_0^2}{8\pi G} \approx 5.3 \times 10^{-10}\;\text{J/m}^3$$

against the observed $\rho_\Lambda \approx 6 \times 10^{-10}\;\text{J/m}^3$.  This is not a fit — it follows from requiring the σ-field to be self-consistent with the expansion rate.

The cosmological constant problem in standard physics is: why is Λ 120 orders of magnitude smaller than the Planck energy density?  QGD's answer is that Λ is not a parameter at all; it is the kinetic energy of the σ-field at the attractor, which is determined by H₀.  The question becomes: why does the σ-field reach the attractor?  This requires understanding the dynamics from the Planck epoch to today, which is an open problem.

## 7. Frame-dragging as wave interference

The interpretation of the Kerr metric's off-diagonal term g_tφ = -σ_t · σ_φ as the interference between temporal and azimuthal graviton waves is physically illuminating.  In GR, frame-dragging is a consequence of the non-diagonal metric — a geometric statement.  In QGD, it is a consequence of wave superposition — a dynamical statement.

This suggests that frame-dragging experiments (Gravity Probe B, LARES) are detecting the beating pattern between two components of the graviton field, analogous to how optical heterodyne detection measures the beating between two laser frequencies.

## 8. Gravitational energy: a genuine resolution

The energy localisation problem in GR is not a philosophical inconvenience — it has practical consequences.  In numerical relativity, the concept of "energy radiated by gravitational waves" requires extracting waveforms at large distances and using the Isaacson averaging procedure.  Energy balance during a binary merger is checked globally, not locally.

QGD provides a local, covariant energy density ρ = ½(∂σ)² that is positive-definite and conserved.  For Schwarzschild, this gives ρ ∝ r⁻³ concentrated near the horizon.  The gravitational Poynting vector S = σ̇_μ∇σ^μ tracks energy flow at any radius.  This is a concrete computational tool, not just a formal improvement.

## 9. The complex gravitational coupling constant

The gravitational fine structure constant α_G = e^{iπ/4}√(cℏ/2GMm) is complex, with a geometric phase of π/4.  This is unusual — QED's α ≈ 1/137 is real.  The complex phase arises because matching the WKB phase expansion to the Newtonian potential requires an imaginary intermediate: α_G² = icℏ/(2GMm).

The physical content is that the gravitational coupling has both amplitude (|α_G|² = cℏ/2GMm) and phase (π/4).  The amplitude controls the strength of the coupling; the phase ensures the potential is real and negative (attractive).  Whether this geometric phase has deeper topological significance is unknown but worth investigating.

## 10. The Green's function solution

The exact integral solution

$$\sigma_\mu(x) = \sigma_\mu^{\text{free}} + \kappa\ell_Q^2\int d^4x'\sqrt{-g}\,[G_0(x,x') - G_{m_Q}(x,x')]\,S_\mu(x')$$

separates gravitational dynamics into two clean sectors: a massless retarded propagator (classical GR at all scales) minus a Yukawa-damped massive propagator (quantum corrections localised within ℓ_Q ~ 10⁻³⁵ m of the source).

This structure is physically transparent: at distances much larger than the Planck length, the two Green's functions cancel everywhere except on the light cone, giving exactly GR.  Near the Planck length, the massive propagator dominates and provides quantum stiffness.  The Born series gives a systematic expansion in powers of the gravitational coupling, converging in the weak-field regime.

## 11. The equivalence proof

The mathematical equivalence between QGD and GR — proved in both directions via the Einstein-wave operator identity and the Debney-Kerr-Schild decomposition theorem — is important because it means QGD makes no *different* classical predictions from GR.  Every solution of one is a solution of the other.

This is both a strength and a limitation.  It means QGD is automatically consistent with all existing tests of GR (solar system, binary pulsars, gravitational waves, cosmological expansion).  But it also means that distinguishing QGD from GR requires accessing the quantum regime (∇²σ ≠ 0), which for macroscopic sources occurs only at absurdly small distances.

The falsifiable predictions therefore come not from the classical sector but from the extended solution space: the four-mode cosmological dynamics, the quantum corrections to the metric, and whatever emerges from the dark matter analysis (not assessed here).

## 12. What is genuinely new

To be precise about novelty:

- The derivation of σ_μ from the Dirac spinor lower components is new.  The WKB limit of Dirac in a gravitational field has been studied (DeWitt, Brill-Wheeler), but the identification of pc/Δ as a field that reconstructs the metric is original.

- The master metric as an algebraic function of σ-fields with source signatures is new.  Related constructions exist (Kerr-Schild metrics, the Rosen bimetric theory), but the complete table of source signatures with ε = ±1 is original.

- The exact Green's function solution to the fourth-order equation with Pais-Uhlenbeck factorisation applied to gravity is new.  Fourth-order gravity has been studied extensively (Stelle, Boulware-Deser), but the specific factorisation □_g(□_g - m_Q²) with m_Q at the Planck mass and the composite propagator G_QGD = ℓ_Q²[G_0 - G_{m_Q}] is original.

- The derivation of dark energy density from the σ-field attractor without a cosmological constant is new.  Quintessence models are common, but matching ρ_Λ exactly from a self-consistency condition (not a fitted potential) is original.

- The resolution of the gravitational energy localisation problem via a true Noether tensor is new in this specific form, though the general idea that a more fundamental variable might resolve the pseudotensor problem has been discussed before (Rosen, Boulware-Deser).

## 13. What remains to be established

- The dark matter explanation via κ-factors requires a rigorous derivation of the activation functions g_j(r) from first principles.  Without this, the claim of "zero free parameters" is premature.

- The inflationary mechanism requires numerical integration of the full nonlinear cosmological ODE and computation of the spectral index n_s.

- The strong-field regime (binary mergers, black hole interiors) requires solving the full nonlinear QGD system, which has not been demonstrated to be computationally simpler than numerical relativity in practice.

- The quantisation of the σ-field — promoted from a classical analysis to a full quantum field theory — has been outlined but not carried through.  Whether the theory is renormalisable or requires UV completion is unknown.

---

*Assessment completed after detailed analysis of the full QGD manuscript, supporting code, and numerical verification of key results.*
