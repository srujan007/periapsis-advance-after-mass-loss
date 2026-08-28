"""Newtonian impulsive mass-loss map for a test planet."""

from __future__ import annotations

import math

from .constants import G


def circular_speed(mass: float, radius: float) -> float:
    """Circular orbital speed about *mass* at *radius* (SI)."""
    return math.sqrt(G * mass / radius)


def specific_energy(mass_final: float, radius: float, speed: float) -> float:
    """Specific orbital energy after an impulsive mass drop."""
    return 0.5 * speed * speed - G * mass_final / radius


def specific_angular_momentum(radius: float, velocity) -> float:
    """Magnitude of specific angular momentum for r = (radius, 0, 0)."""
    vx, vy, vz = velocity
    hx = 0.0
    hy = -radius * vz
    hz = radius * vy
    return math.sqrt(hx * hx + hy * hy + hz * hz)


def no_kick_elements(mass_initial: float, mass_final: float, radius: float):
    """Return (a, e) after a no-kick impulsive drop from a circular orbit.

    Bound only if mass_final > mass_initial / 2. The explosion point is
    the new periapsis: r_p = radius.
    """
    mu = mass_final / mass_initial
    if mu <= 0.5:
        return float("nan"), float("nan")
    a = radius * mu / (2.0 * mu - 1.0)
    e = (1.0 - mu) / mu
    return a, e
