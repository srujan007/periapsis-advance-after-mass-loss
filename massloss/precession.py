"""First-order Schwarzschild periapsis advance."""

from __future__ import annotations

import math

from .constants import ARCSEC_PER_RAD, C, G, JULIAN_CENTURY_S
from .kepler import no_kick_elements


def precession_per_orbit_rad(mass: float, a: float, e: float) -> float:
    """Einstein advance per radial period, in radians.

    dphi = 6 pi G M / [c^2 a (1 - e^2)]
    """
    if e >= 1.0 or a <= 0.0 or not math.isfinite(e):
        return float("nan")
    return 6.0 * math.pi * G * mass / (C * C * a * (1.0 - e * e))


def orbital_period(mass: float, a: float) -> float:
    """Sidereal period of a Keplerian ellipse (seconds)."""
    return 2.0 * math.pi * math.sqrt(a**3 / (G * mass))


def precession_arcsec_per_century(mass: float, a: float, e: float) -> float:
    """Accumulated first-order Schwarzschild advance per Julian century."""
    dphi = precession_per_orbit_rad(mass, a, e)
    if not math.isfinite(dphi):
        return float("nan")
    n_orbits = JULIAN_CENTURY_S / orbital_period(mass, a)
    return dphi * n_orbits * ARCSEC_PER_RAD


def no_kick_advance_ratio(mass_initial: float, mass_final: float) -> float:
    """Per-orbit residual / pre-collapse advance after a no-kick drop.

    Equals (M_f / M_i)^2 for any initially circular radius, provided
    M_f > M_i / 2.
    """
    mu = mass_final / mass_initial
    if mu <= 0.5:
        return float("nan")
    return mu * mu


def no_kick_century_rate_ratio(mass_initial: float, mass_final: float) -> float:
    """Century-rate residual / pre-collapse advance, including period change."""
    mu = mass_final / mass_initial
    if mu <= 0.5:
        return float("nan")
    a_over_r = mu / (2.0 * mu - 1.0)
    period_ratio = (a_over_r**1.5) / math.sqrt(mu)
    return (mu * mu) / period_ratio


def verify_identity(mass_initial: float, mass_final: float, radius: float) -> float:
    """Numerically evaluate dphi_f / dphi_i; should match (M_f / M_i)^2."""
    a_f, e_f = no_kick_elements(mass_initial, mass_final, radius)
    dphi_i = precession_per_orbit_rad(mass_initial, radius, 0.0)
    dphi_f = precession_per_orbit_rad(mass_final, a_f, e_f)
    return dphi_f / dphi_i
