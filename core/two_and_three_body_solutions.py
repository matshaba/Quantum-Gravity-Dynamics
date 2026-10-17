"""
two_and_three_body_solutions.py -- N-body sigma-field superposition  (RECONSTRUCTED)
===================================================================================

The original file is not in the repository.  core/QGD_vs_GR.py imports

    from two_and_three_body_solutions import BinaryInspiral, Body, QGDNBody

and its Section 5 and 6 checks are the specification used here.  Everything below
is written from how those names are used and from the formulas in the QGD notes:

    sigma_t(x) = sum_a sqrt(2 G M_a / (c^2 |x - x_a|))        (the "exact N-body solution")
    sigma_t^2  = self terms + cross terms
               = sum_a sigma_a^2  +  2 sum_{a<b} sigma_a sigma_b

CHOICES THE SPECIFICATION LEAVES OPEN (stated so they can be revisited)

 1. `binary_merger_prediction()` returns the separation d at which the SELF part of
    sigma^2 at the centre of mass equals 1.  For two masses with the COM at
    r_a = d M_b / (M_a + M_b):

        sum_a 2 G M_a / (c^2 r_a) = 2 G M_tot (M_a^2 + M_b^2) / (c^2 d M_a M_b) = 1,

    so d = 2 G M_tot (M_a^2 + M_b^2) / (c^2 M_a M_b), which is 4 r_s(each) for equal
    masses and about 4.09 r_s(avg) for 36 + 29 solar masses; this is the value the
    validation script expects ("d = 4 r_s").  If the TOTAL sigma (self + cross) is set to
    1 instead, the equal-mass result is 8 r_s, not 4 r_s.  The factor of two is therefore
    a choice of which condition is called the merger, not something derived.
 2. Positions and velocities are plain arrays; no dynamics is integrated here.

What this module does NOT claim: that sigma_t above solves the Einstein equations for
two bodies.  It does not (see tests/test_gaps.py::test_two_body_superposition_is_not_vacuum).
"""
from dataclasses import dataclass, field

import numpy as np

G = 6.67430e-11
c = 2.99792458e8
M_SUN = 1.98847e30


@dataclass
class Body:
    mass: float
    position: list
    velocity: list = field(default_factory=lambda: [0.0, 0.0, 0.0])

    def __post_init__(self):
        self.position = np.asarray(self.position, dtype=float)
        self.velocity = np.asarray(self.velocity, dtype=float)


class QGDNBody:
    """sigma-field superposition for N point masses."""

    def __init__(self, bodies):
        self.bodies = list(bodies)

    def _sigmas(self, x, positions=None):
        x = np.asarray(x, dtype=float)
        pos = np.array([b.position for b in self.bodies]) if positions is None else np.asarray(positions, float)
        out = []
        for b, p in zip(self.bodies, pos):
            r = np.linalg.norm(x - p)
            out.append(np.sqrt(2 * G * b.mass / (c ** 2 * r)) if r > 0 else np.inf)
        return np.array(out)

    def sigma_t(self, x, positions=None):
        """Total sigma_t at x."""
        return float(self._sigmas(x, positions).sum())

    def sigma_squared_decomposed(self, x, positions=None):
        """(self, cross) with self = sum_a sigma_a^2 and cross = sigma_total^2 - self."""
        s = self._sigmas(x, positions)
        self_t = float(np.sum(s ** 2))
        total = float(s.sum()) ** 2
        return self_t, total - self_t

    def center_of_mass(self):
        m = np.array([b.mass for b in self.bodies])
        p = np.array([b.position for b in self.bodies])
        return (m[:, None] * p).sum(axis=0) / m.sum()

    def binary_merger_prediction(self):
        """Separation at which the self part of sigma^2 at the COM equals 1 (two bodies)."""
        if len(self.bodies) != 2:
            raise ValueError("binary_merger_prediction needs exactly two bodies")
        ma, mb = self.bodies[0].mass, self.bodies[1].mass
        return 2 * G * (ma + mb) * (ma ** 2 + mb ** 2) / (c ** 2 * ma * mb)


class BinaryInspiral(QGDNBody):
    """Two bodies described by masses, dimensionless spins and a distance to the observer."""

    def __init__(self, M1, M2, chi1=0.0, chi2=0.0, distance=None):
        super().__init__([Body(M1, [0.0, 0.0, 0.0]), Body(M2, [1.0, 0.0, 0.0])])
        self.M1, self.M2, self.chi1, self.chi2, self.distance = M1, M2, chi1, chi2, distance

    @property
    def total_mass(self):
        return self.M1 + self.M2

    def isco_radius(self):
        """GR test-mass ISCO of the total mass, 6 G M / c^2."""
        return 6 * G * self.total_mass / c ** 2

    def isco_gw_frequency(self):
        """GW frequency at that radius, f = sqrt(G M / r^3) / pi."""
        return np.sqrt(G * self.total_mass / self.isco_radius() ** 3) / np.pi
