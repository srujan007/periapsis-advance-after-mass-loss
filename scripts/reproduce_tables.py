#!/usr/bin/env python3
"""Reproduce Tables I and II of the companion AJP note.

Usage
-----
    python scripts/reproduce_tables.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from massloss import (
    AU,
    MSUN,
    bound_fraction,
    no_kick_advance_ratio,
    no_kick_elements,
    precession_arcsec_per_century,
)
from massloss.kepler import circular_speed
from massloss.precession import verify_identity


def mercury_check() -> float:
    a_merc = 0.387098 * AU
    e_merc = 0.205630
    return precession_arcsec_per_century(MSUN, a_merc, e_merc)


def table_i(mass_initial: float = 15 * MSUN) -> None:
    print("Table I  -  no-kick remnants of a 15 Msun circular progenitor")
    print(f"{'mu':>8} {'e':>8} {'af/r':>8} {'mu2':>8} {'0.3 AU':>10} {'1 AU':>10} {'5 AU':>10}")
    for mu in (0.95, 0.90, 0.80, 0.70, 0.60, 0.533, 0.510):
        mf = mu * mass_initial
        a, e = no_kick_elements(mass_initial, mf, AU)
        a_over_r = a / AU
        row = [f"{mu:8.3f}", f"{e:8.4f}", f"{a_over_r:8.3f}", f"{mu * mu:8.3f}"]
        for r_au in (0.3, 1.0, 5.0):
            aa, ee = no_kick_elements(mass_initial, mf, r_au * AU)
            dphi = precession_arcsec_per_century(mf, aa, ee)
            row.append(f"{dphi:10.3f}")
        print(" ".join(row))
    print()
    print("Identity check dphi_f/dphi_i vs mu^2 at 1 AU:")
    for mu in (0.80, 0.533):
        numerical = verify_identity(mass_initial, mu * mass_initial, AU)
        print(f"  mu={mu:.3f}  numerical={numerical:.12f}  mu^2={mu * mu:.12f}")


def table_ii(mass_initial: float = 15 * MSUN, n: int = 30_000) -> None:
    print()
    print(f"Table II  -  isotropic fixed-speed kicks, {n} draws per cell")
    configs = [
        ("1 AU, 8 Msun", 1.0, 8.0),
        ("1 AU, 6 Msun", 1.0, 6.0),
        ("1 AU, 1.4 Msun", 1.0, 1.4),
        ("5 AU, 8 Msun", 5.0, 8.0),
        ("10 AU, 8 Msun", 10.0, 8.0),
    ]
    header = f"{'vk':>8}" + "".join(f"{name:>16}" for name, _, _ in configs)
    print(header)
    for vk in (0, 20, 50, 100, 150, 200, 265):
        cells = [f"{vk:8.0f}"]
        for _, r_au, mf_msun in configs:
            frac = bound_fraction(
                mass_initial,
                mf_msun * MSUN,
                r_au * AU,
                vk,
                n=n,
                kind="fixed",
                seed=42,
            )
            cells.append(f"{frac:16.4f}")
        print(" ".join(cells))


def speeds() -> None:
    print()
    print("Circular speeds for a 15 Msun host")
    for r_au in (0.3, 1.0, 5.0, 10.0, 30.0):
        v = circular_speed(15 * MSUN, r_au * AU) / 1000.0
        print(f"  r = {r_au:4.1f} AU   v = {v:6.2f} km/s")


def main() -> None:
    print(f"Mercury century-rate check: {mercury_check():.2f} arcsec / century")
    print("(canonical value 43)\n")
    table_i()
    table_ii()
    speeds()
    print()
    print("Closed-form per-orbit ratio is mu^2 =", no_kick_advance_ratio(15 * MSUN, 8 * MSUN))


if __name__ == "__main__":
    main()
