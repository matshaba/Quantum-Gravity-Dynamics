%============================================================
% QGD GRAVITATIONAL ENERGY: A COMPLETE RESOLUTION
% Based on QGD field theory documents
%============================================================
\documentclass[12pt,a4paper]{article}
\usepackage{amsmath,amssymb,booktabs,geometry,hyperref,xcolor}
\geometry{margin=2.5cm}
\hypersetup{colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue}

\title{\textbf{Gravitational Energy in Quantum Gravitational Dynamics}\\[6pt]
\large Resolution of the Century-Long Localisation Problem}
\author{QGD Framework Analysis}
\date{\today}

\begin{document}
\maketitle
\tableofcontents
\newpage

%============================================================
\section{Introduction: The Longstanding Problem}
%============================================================

General Relativity (GR) has a fundamental deficiency that has persisted since 1915: it possesses no
well-defined, coordinate-independent, local energy-momentum tensor for the gravitational field itself.
The standard remedy — the Einstein pseudotensor $t_{\mu\nu}$ — is coordinate-dependent, can be made
to vanish at any point, and is manifestly non-covariant.  Positive total energy was not rigorously
established until 1979 (Schoen-Yau) and 1981 (Witten), and only at spatial infinity.

Quantum Gravitational Dynamics (QGD) resolves all of these issues in a single stroke, by changing the
fundamental variable from the metric $g_{\mu\nu}$ to a vector field $\sigma_\mu$.  The gravitational
energy-momentum tensor in QGD is a true Noether tensor, positive-definite, coordinate-independent, and
localizable at every point in spacetime.

%============================================================
\section{The QGD Energy-Momentum Tensor}
%============================================================

\subsection{Action and Noether Derivation}

The QGD action for the $\sigma$-field reads:
\begin{equation}
  S[\sigma] = \int d^4x \left[
    \frac{1}{2}\partial_\alpha\sigma_\mu\partial^\alpha\sigma^\mu
    - V(\sigma)
    + \frac{\ell_Q^2}{2}(\partial_\alpha\partial_\beta\sigma_\mu)^2
    + \mathcal{L}_{\mathrm{matter}}[g[\sigma]]
  \right]
\end{equation}
where $\ell_Q = \ell_{\rm Pl} = \sqrt{\hbar G/c^3} \approx 1.616\times10^{-35}$~m is the Planck length,
and the metric is a functional of $\sigma_\mu$: $g_{\mu\nu}[\sigma]$.

Applying Noether's theorem to translation invariance $x^\mu\to x^\mu+a^\mu$ yields the canonical
energy-momentum tensor:
\begin{equation}
  \boxed{
    T^{\mu\nu}_{\rm QGD}
    = \partial^\mu\sigma_\alpha\partial^\nu\sigma^\alpha
    - \frac{1}{2}\eta^{\mu\nu}\partial_\beta\sigma_\alpha\partial^\beta\sigma^\alpha
    - \eta^{\mu\nu}V(\sigma)
    + T^{\mu\nu}_{\rm quantum}
  }
  \label{eq:Tmunu}
\end{equation}
where $T^{\mu\nu}_{\rm quantum}$ contains the fourth-order ($\ell_Q^2$) derivative corrections.

\textbf{This is a genuine tensor.}  It transforms properly under all coordinate transformations and
cannot be made to vanish by any choice of gauge.

\subsection{Local Energy Density}

The $T^{00}$ component gives the gravitational energy density at every spacetime point:
\begin{equation}
  \boxed{
    \rho_{\rm grav}(\mathbf{x}) = \frac{1}{2}\dot\sigma_\mu\dot\sigma^\mu
    + \frac{1}{2}(\nabla\sigma_\mu)^2 + V(\sigma) + \rho_{\rm quantum}
  }
  \label{eq:rho}
\end{equation}
with the standard field-theory decomposition into kinetic, gradient, potential, and quantum-stiffness
terms.  Each term is manifestly non-negative for $V\ge0$.

\subsection{Conservation and the Positive-Energy Theorem}

\begin{itemize}
\item \textbf{Conservation:} $\partial_\mu T^{\mu\nu}_{\rm QGD} = 0$ follows automatically from the
  field equations via Noether's theorem.  No pseudotensor gymnastics are required.
\item \textbf{Positive energy:} The QGD Hamiltonian
  \begin{equation}
    H[\sigma,\pi] = \int d^3x\left[
      \frac{1}{2}\pi_\mu\pi^\mu + \frac{1}{2}(\nabla\sigma_\mu)^2
      + V(\sigma) + \frac{\ell_Q^2}{2}(\nabla^2\sigma_\mu)^2
    \right] \geq 0
  \end{equation}
  is a sum of manifestly non-negative terms.  Equality $H=0$ if and only if $\sigma_\mu=0$
  (flat Minkowski spacetime).  This trivially supersedes the 64-year effort that culminated
  in the Schoen-Yau spinor proof.
\end{itemize}

%============================================================
\section{Energy Density of the Schwarzschild Field}
%============================================================

\subsection{The $\sigma$-Field Solution}

For a Schwarzschild black hole of mass $M$, the time component of the $\sigma$-field is:
\begin{equation}
  \sigma_t(r) = \sqrt{\frac{r_s}{r}}, \qquad r_s = \frac{2GM}{c^2}
\end{equation}
This encodes the metric via $g_{tt} = -(1 - \sigma_t^2)$, exactly recovering the Schwarzschild solution.

\subsection{The Gravitational Energy Density}

The gradient energy density is:
\begin{equation}
  \rho_{\rm grav}(r) = \frac{1}{2}\left(\frac{d\sigma_t}{dr}\right)^2
  = \frac{1}{2}\cdot\frac{r_s}{4r^3}
  = \boxed{\frac{GM}{4c^2 r^3}}
  \label{eq:rho_schw}
\end{equation}
This is a well-defined scalar field that lives at every point outside the horizon.  In GR, the
equivalent statement cannot be made; the pseudotensor gives different values in different coordinates.

\subsection{Physical Interpretation}

The energy is distributed throughout space, concentrated near the horizon with profile $\rho\propto r^{-3}$:
\begin{itemize}
  \item At $r=r_s$: $\rho = GM/(4c^2r_s^3) = c^4/(32G^2M^2)$
  \item At $r=10r_s$: $\rho = 10^{-3}$ of the horizon value
  \item No energy is deposited at the singularity; the quantum correction term $\ell_Q^2(\nabla^2\sigma)^2$
    activates at $r\sim\ell_{\rm Pl}$ and prevents unlimited concentration
\end{itemize}

\textbf{Note on total energy:}  The integral $\int_{\ell_{\rm Pl}}^\infty \rho_{\rm grav}\,4\pi r^2\,dr$
diverges logarithmically at large $r$ in geometrised units, just as the Newtonian gravitational
self-energy of a point mass diverges.  The physical energy $Mc^2$ is identified via the ADM limit
in the weak-field regime rather than as a naive volume integral.  This is consistent: the \emph{density}
is localised, while the \emph{total} energy is captured by the asymptotic field strength.

%============================================================
\section{The $Q_\mu$ Self-Energy Term: Dark Matter Connection}
%============================================================

\subsection{Physical Meaning of $Q_\mu$}

The QGD field equation reads:
\begin{equation}
  \Box_g\sigma_\mu
  = \underbrace{\frac{8\pi G}{c^4}Q_\mu}_{\text{self-interaction}}
  + \frac{8\pi G}{c^4}G_\mu
  + \frac{4\pi G}{c^2}T^{\mu\nu}\sigma_\nu
\end{equation}
where $Q_\mu = \sigma_\mu(\partial\sigma)^2 + \ldots$ is quadratic in $\sigma$ and its derivatives.

$Q_\mu$ is the \emph{gravitational field's self-energy current}: the source contribution from the field
energy itself, exactly analogous to the non-Abelian self-coupling $j_{\rm self}\sim A\cdot F$ in
Yang-Mills theory.

\subsection{Expansion in Orders}

At higher orders in the $\sigma$ expansion:
\begin{equation}
  Q_\mu = Q_\mu^{(2)} + Q_\mu^{(3)} + Q_\mu^{(4)} + \cdots
\end{equation}
Each term contributes energy:
\begin{equation}
  \rho_{\rm total} = \rho_{\rm matter}
  \left(1 + \sum_{n=1}^\infty \kappa_n^2\right)
  \quad\text{where}\quad
  \kappa_n = \sqrt{\frac{(2n-1)!}{2^{2n-2}}}
\end{equation}

\subsection{Dark Matter as Gravitational Self-Energy}

\begin{equation}
  \boxed{
    \rho_{\rm dark\,matter} = \text{gravitational field self-energy}
    = \sum_{n=2}^\infty Q_\mu^{(n)}
  }
\end{equation}

Dark matter, in QGD, is \emph{not} an exotic particle.  It is the energy stored in the gravitational
field's own nonlinear self-interaction terms, which were hidden inside Einstein's $G_{\mu\nu}$ and
therefore invisible in GR.  The observed velocity enhancement factors $\kappa_n$ (Table~\ref{tab:kappa})
emerge from this series with zero free parameters.

\begin{table}[h]
\centering
\caption{Velocity enhancement factors from gravitational self-energy}
\label{tab:kappa}
\begin{tabular}{@{}cccl@{}}
\toprule
$n$ & $\kappa_n$ & Physical regime & Observation \\
\midrule
1 & 1.000 & Solar System / Newtonian & GPS, planetary orbits \\
2 & 1.225 & Wide binary stars & Chae et al.\ (2023) anomaly \\
3 & 2.739 & Spiral galaxy outskirts & Flat rotation curves \\
4 & 8.874 & Galaxy clusters / CMB & Missing mass factor \\
5 & 37.65 & Protoclusters & Early structure formation \\
\bottomrule
\end{tabular}
\end{table}

%============================================================
\section{Two-Body Binding Energy}
%============================================================

\subsection{Field-Gradient Derivation}

For two masses $M_1$, $M_2$ separated by distance $d$, the $\sigma$-field is:
\begin{equation}
  \sigma_t(\mathbf{x}) = \sqrt{\frac{2GM_1}{c^2|\mathbf{x}-\mathbf{x}_1|}}
                       + \sqrt{\frac{2GM_2}{c^2|\mathbf{x}-\mathbf{x}_2|}}
\end{equation}
The total field energy splits as $E = E_1 + E_2 + E_{\rm cross}$, where the cross term from
$\int\nabla\sigma_1\cdot\nabla\sigma_2\,d^3x$ gives (via Green's theorem):
\begin{equation}
  \boxed{E_{\rm binding} = -\frac{GM_1 M_2}{d}}
\end{equation}
Newton's gravitational binding energy emerges \emph{automatically} from the field-gradient structure,
with no additional assumptions.  In GR this relationship holds only because the linearised metric
reproduces Newtonian gravity, but the energy localisation is lost.

%============================================================
\section{Gravitational Wave Energy}
%============================================================

\subsection{Exact, Point-Wise Definition}

In QGD, a gravitational wave is a radiation mode:
$\sigma_\mu^{\rm GW} = \epsilon_\mu e^{ik\cdot x}$, $k^2=0$.

The energy density is:
\begin{equation}
  \rho_{\rm GW} = \frac{1}{2}(\partial\sigma^{\rm GW})^2
  = \frac{1}{2}\omega^2|\epsilon|^2
\end{equation}
This is \textbf{exact, point-wise defined, and gauge-invariant}.  No averaging over multiple wavelengths
is required.

\subsection{Comparison with Isaacson}

In the appropriate weak-field limit the LIGO observable gives:
\begin{equation}
  h_{ij}^{\rm TT} \sim -2\sigma_i\sigma_j
  \;\Longrightarrow\;
  \rho_{\rm GW} \longrightarrow \frac{c^2}{32\pi G}\langle\dot h^2\rangle
\end{equation}
recovering the Isaacson (1968) formula, which in GR requires an averaging procedure that is approximate
and only valid in the short-wave limit.

\subsection{Gravitational Poynting Vector}

The energy flux at any radius $r$ (not just at infinity) is:
\begin{equation}
  \mathbf{S}_{\rm grav} = \dot\sigma_\mu\,\nabla\sigma^\mu,
  \qquad
  \frac{dE}{dt} = \oint_S \mathbf{S}_{\rm grav}\cdot d\mathbf{A}
\end{equation}
This allows energy-flow tracking from the near-zone orbital region to the wave zone to infinity — a
capability absent in GR.

%============================================================
\section{Comparison: GR vs.\ QGD Energy}
%============================================================

\begin{table}[h]
\centering
\caption{Energy properties: GR vs.\ QGD}
\label{tab:comparison}
\begin{tabular}{@{}p{3.8cm}p{4.8cm}p{5cm}@{}}
\toprule
\textbf{Property} & \textbf{GR} & \textbf{QGD} \\
\midrule
Local energy density
  & Undefined; pseudotensor $t_{\mu\nu}$ only
  & $\rho = \tfrac{1}{2}(\partial\sigma)^2$ — scalar field \\
\addlinespace
True energy-momentum tensor
  & Does not exist
  & Noether tensor $T^{\mu\nu}_{\rm QGD}$ \\
\addlinespace
Energy conservation
  & Ambiguous (coordinate-dependent)
  & Automatic: $\partial_\mu T^{\mu\nu}=0$ \\
\addlinespace
Positive energy proof
  & Schoen-Yau spinors (1979); only at $r\to\infty$
  & $H=\int(\text{positive terms})\ge0$; manifestly local \\
\addlinespace
Energy localisation
  & Cannot be defined at a point
  & $\rho(\mathbf{x})$ is a coordinate-independent scalar \\
\addlinespace
Gravitational waves
  & Isaacson averaging (approximate, weak field)
  & Exact: $\rho_{\rm GW}=\tfrac{1}{2}(\partial\sigma^{\rm GW})^2$ \\
\addlinespace
Black hole energy
  & Defined only at $r=\infty$ (ADM)
  & Distributed: $\rho(r) \propto r^{-3}$ \\
\addlinespace
Binding energy
  & Follows from linearised GR; local source undefined
  & $-GM_1M_2/r$ from field-gradient cross term \\
\addlinespace
Self-energy / dark matter
  & Hidden in $G_{\mu\nu}$; unexplained
  & Explicit $Q_\mu$ series; $\kappa$-ladder \\
\addlinespace
Quantum corrections
  & No Hamiltonian formulation
  & Built-in: $\ell_Q^2(\nabla^2\sigma)^2$ \\
\bottomrule
\end{tabular}
\end{table}

%============================================================
\section{The Conceptual Shift}
%============================================================

The equivalence principle states that the metric $g_{\mu\nu}$ can be made flat locally (free-fall frame).
This is precisely why GR's pseudotensor fails: it measures geometry, and geometry can be made to
vanish by coordinate choice.

QGD resolves this by identifying the \emph{field generating the metric} — $\sigma_\mu$ — as the
fundamental variable.  $\sigma_\mu$ cannot be gauged away, just as the Yang-Mills field strength
$F_{\mu\nu}$ cannot be gauged away even though the potential $A_\mu$ can.

\begin{center}
\begin{tabular}{@{}ccc@{}}
  \textbf{Gravity (QGD)} & $\leftrightarrow$ & \textbf{Electromagnetism} \\
  $\sigma_\mu$ (field) & $\leftrightarrow$ & $F_{\mu\nu}$ (field strength) \\
  $g_{\mu\nu}$ (metric) & $\leftrightarrow$ & $A_\mu$ (potential) \\
  $Q_\mu$ (self-energy) & $\leftrightarrow$ & $j_{\rm self}$ (non-Abelian) \\
\end{tabular}
\end{center}

\textbf{The century-long mystery was a consequence of the wrong fundamental variable.}

%============================================================
\section{Conclusions}
%============================================================

QGD resolves General Relativity's gravitational energy problem completely:

\begin{enumerate}
\item A true, coordinate-independent energy-momentum tensor $T^{\mu\nu}_{\rm QGD}$ exists and is
  derived by Noether's theorem.
\item The gravitational energy density $\rho_{\rm grav}=\tfrac{1}{2}(\partial\sigma)^2$ is a
  well-defined scalar field, localised at every point.
\item Energy conservation is automatic; the positive-energy theorem is manifest.
\item Gravitational wave energy is exact and point-wise, recovering Isaacson in the appropriate limit.
\item The gravitational self-energy $Q_\mu$ explains dark matter without exotic particles: the
  $\kappa$-ladder $\kappa_n=\sqrt{(2n-1)!/2^{2n-2}}$ follows with zero free parameters.
\item Two-body binding energy $-GM_1M_2/r$ emerges naturally from field-gradient cross terms.
\end{enumerate}

These results are not modifications of GR but a reformulation that makes the field-theoretic
structure of gravity manifest — a structure that was always latent in Einstein's equations but
obscured by the geometric language.

\end{document}
