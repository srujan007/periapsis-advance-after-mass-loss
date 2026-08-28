"""Isotropic natal-kick Monte Carlo in the test-planet limit.

Geometry
--------
At explosion the planet sits at r = (radius, 0, 0) with circular velocity
v = (0, v_circ, 0) about the progenitor. The remnant receives a kick
velocity v_k. In the remnant rest frame the planet's velocity is
v - v_k. The planet remains bound if

    |v - v_k|^2  <  2 G M_f / r.
"""

from __future__ import annotations

from typing import Literal

import numpy as np

from .constants import G


def _draw_isotropic_unit_vectors(n: int, rng: np.random.Generator) -> np.ndarray:
    """Return shape (n, 3) unit vectors uniform on the sphere."""
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    costh = rng.uniform(-1.0, 1.0, n)
    sinth = np.sqrt(np.clip(1.0 - costh * costh, 0.0, 1.0))
    return np.column_stack(
        (sinth * np.cos(phi), sinth * np.sin(phi), costh)
    )


def _kick_speeds(
    n: int,
    kick_kms: float,
    kind: Literal["fixed", "maxwellian"],
    rng: np.random.Generator,
) -> np.ndarray:
    if kind == "fixed":
        return np.full(n, kick_kms * 1000.0)
    if kind == "maxwellian":
        sigma = kick_kms * 1000.0
        gauss = rng.normal(0.0, 1.0, size=(n, 3))
        return sigma * np.linalg.norm(gauss, axis=1)
    raise ValueError(f"unknown kick kind: {kind!r}")


def relative_velocity(
    mass_initial: float,
    radius: float,
    kick_kms: float,
    n: int = 40_000,
    kind: Literal["fixed", "maxwellian"] = "fixed",
    seed: int = 42,
) -> np.ndarray:
    """Planet velocity in the remnant frame, shape (n, 3), SI."""
    rng = np.random.default_rng(seed)
    v = np.sqrt(G * mass_initial / radius)
    speed = _kick_speeds(n, kick_kms, kind, rng)
    unit = _draw_isotropic_unit_vectors(n, rng)
    vk = unit * speed[:, None]
    v_planet = np.array([0.0, v, 0.0])
    return v_planet - vk


def bound_mask(
    mass_initial: float,
    mass_final: float,
    radius: float,
    kick_kms: float,
    n: int = 40_000,
    kind: Literal["fixed", "maxwellian"] = "fixed",
    seed: int = 42,
) -> np.ndarray:
    """Boolean mask of draws that remain bound after the kick."""
    vrel = relative_velocity(mass_initial, radius, kick_kms, n=n, kind=kind, seed=seed)
    v2 = np.einsum("ij,ij->i", vrel, vrel)
    return v2 < (2.0 * G * mass_final / radius)


def bound_fraction(
    mass_initial: float,
    mass_final: float,
    radius: float,
    kick_kms: float,
    n: int = 40_000,
    kind: Literal["fixed", "maxwellian"] = "fixed",
    seed: int = 42,
) -> float:
    """Fraction of isotropic kicks that leave the test planet bound."""
    return float(
        bound_mask(
            mass_initial, mass_final, radius, kick_kms, n=n, kind=kind, seed=seed
        ).mean()
    )


def survivor_elements(
    mass_initial: float,
    mass_final: float,
    radius: float,
    kick_kms: float,
    n: int = 40_000,
    kind: Literal["fixed", "maxwellian"] = "fixed",
    seed: int = 1,
):
    """Keplerian (a, e) of bound survivors.

    Returns
    -------
    a, e : ndarray
        Semi-major axes (m) and eccentricities of bound draws.
    """
    vrel = relative_velocity(mass_initial, radius, kick_kms, n=n, kind=kind, seed=seed)
    v2 = np.einsum("ij,ij->i", vrel, vrel)
    energy = 0.5 * v2 - G * mass_final / radius
    # r = (radius, 0, 0), so h = r x v = (0, -r vz, r vy)
    h2 = (radius * vrel[:, 1]) ** 2 + (radius * vrel[:, 2]) ** 2
    bound = energy < 0.0
    a = np.full(n, np.nan)
    e = np.full(n, np.nan)
    a[bound] = -G * mass_final / (2.0 * energy[bound])
    ecc2 = 1.0 + 2.0 * energy * h2 / (G * G * mass_final * mass_final)
    e[bound] = np.sqrt(np.clip(ecc2[bound], 0.0, None))
    return a[bound], e[bound]
