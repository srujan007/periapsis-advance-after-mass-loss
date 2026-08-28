# Residual Schwarzschild advance after impulsive stellar mass loss

Companion code for the note

> S. S. Chunduri, *Periapsis advance of a test planet after impulsive stellar mass loss*,
> manuscript prepared for the *American Journal of Physics* (2026).

The repository does one algebraic thing and one computational thing.

1. After an impulsive, spherically symmetric, **no-kick** drop from progenitor mass `M_i` to remnant mass `M_f > M_i/2`, a test planet that began on a circular orbit of any radius `r` has residual first-order Einstein periapsis advance

       Delta phi_f / Delta phi_i  =  (M_f / M_i)^2.

   The original radius is the new periapsis. The new eccentricity is `(M_i - M_f)/M_f`.

2. An isotropic **natal-kick Monte Carlo** in the test-planet limit reports the bound fraction as a function of kick speed, orbital radius, and remnant mass, plus the `(a, e)` distribution of survivors.

A neutron star and a nonrotating black hole of equal mass produce the same first-order exterior advance. The event horizon does not enter the formulae.

This is pedagogical code. It is not a population-synthesis survey and it is not a new theory of collapse.

## Requirements

Python 3.10 or newer.

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Reproduce the tables

```bash
python scripts/reproduce_tables.py
```

You should see:

- Mercury's century-rate at **42.98 arcsec** (canonical value 43).
- Table I: no-kick elements and residual century-rates for a 15 solar-mass host.
- The identity check `Delta phi_f / Delta phi_i = mu^2` to machine precision.
- Table II: bound fractions for fixed-speed isotropic kicks.

## Make the figures

```bash
python scripts/make_figures.py --outdir figures
```

Writes `figures/fig1_precession_ratio.png` and `figures/fig2_bound_fraction.png`.

## Run the tests

```bash
python tests/test_identity.py
```

## Library API

```python
from massloss import AU, MSUN, no_kick_elements, precession_arcsec_per_century
from massloss import bound_fraction, survivor_elements
from massloss.precession import verify_identity

# No-kick ellipse after 15 -> 8 solar masses at 1 AU
a, e = no_kick_elements(15 * MSUN, 8 * MSUN, AU)

# First-order century-rate on that ellipse
dphi = precession_arcsec_per_century(8 * MSUN, a, e)

# Confirm the closed-form ratio
assert abs(verify_identity(15 * MSUN, 8 * MSUN, AU) - (8 / 15) ** 2) < 1e-12

# Fraction of 50 km/s isotropic kicks that leave the planet bound
f = bound_fraction(15 * MSUN, 8 * MSUN, AU, 50.0, n=30_000, kind="fixed")

# Maxwellian kick with 1-D dispersion 20 km/s
f_fb = bound_fraction(
    15 * MSUN, 8 * MSUN, AU, 20.0, n=40_000, kind="maxwellian"
)

# Elements of the survivors
a_s, e_s = survivor_elements(15 * MSUN, 8 * MSUN, AU, 50.0, n=30_000)
```

Kick speeds are entered in **km/s**. Masses and lengths inside the library are SI.

## Physics in brief

At the instant of mass loss the planet keeps its position and velocity. Bound orbits require `E < 0`, hence `M_f > M_i/2`. The new Keplerian elements are

    a_f = r * mu / (2*mu - 1),    e = (1 - mu)/mu,    mu = M_f / M_i.

Einstein's first-order advance per radial period is

    Delta phi = 6 pi G M / [c^2 a (1 - e^2)].

On the no-kick ellipse, `a_f (1-e^2) = r M_i / M_f`, which is why the ratio of advances collapses to `mu^2`.

With a natal kick v_k the bound test is `|v - v_k|^2 < 2 G M_f / r`.

Directions are drawn uniformly on the sphere (uniform azimuth, uniform cos theta). A Maxwellian draw uses the Euclidean norm of three independent Gaussians of one-dimensional dispersion sigma.

Inner radii (0.3-1 AU around a 15 solar-mass red supergiant) are mathematical controls. Those orbits sit inside the stellar envelope and are not live planets at core collapse.

See `docs/PHYSICS.md` for the full derivation.

## Layout

```
massloss/                 library
  constants.py            G, c, Msun, AU
  kepler.py               circular speed, no-kick (a, e)
  precession.py           Einstein advance, mu^2 identity
  kicks.py                isotropic and Maxwellian Monte Carlo
scripts/
  reproduce_tables.py     Tables I and II
  make_figures.py         Figs. 1 and 2
tests/
  test_identity.py
docs/
  PHYSICS.md              longer derivation
```

## What this repository does not contain

- The 1920 Einstein-Minkowski source volume.
- Word drafts of the manuscript.
- Hierarchical triples, envelope stripping, or a time-dependent exterior metric during the explosion.

Hills (1983) already owns the Newtonian binary-with-kick problem. This repo restricts that problem to a test planet and attaches the first-order Schwarzschild advance.

## License

MIT. See `LICENSE`.
