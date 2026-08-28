"""Sanity checks for the mu^2 identity and the no-kick bound threshold."""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from massloss import AU, MSUN, bound_fraction, no_kick_elements
from massloss.precession import (
    precession_arcsec_per_century,
    verify_identity,
)


def test_unbound_below_half():
    a, e = no_kick_elements(15 * MSUN, 7 * MSUN, AU)
    assert math.isnan(a) and math.isnan(e)


def test_mu_squared_identity():
    mi = 15 * MSUN
    for mu, radius in ((0.8, AU), (0.533, 5 * AU), (0.9, 0.3 * AU)):
        ratio = verify_identity(mi, mu * mi, radius)
        assert abs(ratio - mu * mu) < 1e-12


def test_mercury():
    value = precession_arcsec_per_century(MSUN, 0.387098 * AU, 0.205630)
    assert abs(value - 42.98) < 0.05


def test_zero_kick_bound_fraction():
    mi = 15 * MSUN
    assert bound_fraction(mi, 8 * MSUN, AU, 0.0, n=2000) == 1.0
    assert bound_fraction(mi, 1.4 * MSUN, AU, 0.0, n=2000) == 0.0


if __name__ == "__main__":
    test_unbound_below_half()
    test_mu_squared_identity()
    test_mercury()
    test_zero_kick_bound_fraction()
    print("all tests passed")
