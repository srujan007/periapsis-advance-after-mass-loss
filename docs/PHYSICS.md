# Derivation notes

These notes expand the library docstrings. Symbols are SI unless stated.

## 1. No-kick mass loss

A test planet on a circular orbit of radius r about mass M_i has speed

    v = sqrt(G M_i / r).

An instantaneous, spherically symmetric drop to M_f leaves r and v unchanged. The new specific energy and specific angular momentum are

    E = v^2/2 - G M_f / r = (G/r)(M_i/2 - M_f),
    h = v r = sqrt(G M_i r).

E < 0 if and only if M_f > M_i/2.

The new semi-major axis is a = -G M_f / (2E). Substituting E gives

    a_f = r * mu / (2 mu - 1),    mu = M_f / M_i.

The eccentricity identity

    e^2 = 1 + 2 E h^2 / (G^2 M_f^2)

collapses to

    e = (1 - mu)/mu = (M_i - M_f)/M_f.

Because v exceeds the new circular speed at radius r, the explosion point is periapsis: r_p = r. Apoapsis is r_a = r/(2 mu - 1).

## 2. Residual Schwarzschild advance

Einstein's first-order periapsis advance per radial period is

    Delta phi = 6 pi G M / [c^2 a (1 - e^2)].

Before the drop, a_i = r and e_i = 0, so

    Delta phi_i = 6 pi G M_i / (c^2 r).

After the drop a short rearrangement of the elements gives the semi-latus rectum

    a_f (1 - e^2) = r M_i / M_f.

Substitution produces

    Delta phi_f = 6 pi G M_f^2 / (c^2 r M_i),

and therefore

    Delta phi_f / Delta phi_i = mu^2.

The ratio contains no factor of G/c^2 and no factor of r. Function `massloss.precession.verify_identity` evaluates both sides numerically.

The century-rate is Delta phi times the number of radial periods in a Julian century. The period ratio

    P_f / P_i = (a_f / r)^{3/2} mu^{-1/2}

diverges as mu -> 1/2 from above, so the century-rate falls faster than mu^2.

## 3. Natal kicks

In the remnant rest frame the planet velocity is v - v_k. Bound if

    |v - v_k|^2 < 2 G M_f / r.

Isotropic directions: azimuth uniform on [0, 2 pi), cosine of polar angle uniform on [-1, 1].

Maxwellian speeds: sigma is the one-dimensional dispersion in km/s; the speed is sigma times the Euclidean norm of three standard normals.

When mu <= 1/2 the no-kick orbit is unbound. A kick aligned with v can still bind the planet by lowering the relative speed. That channel is a small solid angle, which is why the 1.4 solar-mass column of Table II is nearly empty at pulsar kick speeds.

## 4. Constants

- G = 6.67430e-11 m^3 kg^{-1} s^{-2}
- c = 2.99792458e8 m s^{-1}
- Msun = 1.98847e30 kg
- AU = 1.495978707e11 m
- Julian century = 36525 days

Mercury check: a = 0.387098 AU, e = 0.205630 returns 42.98 arcsec per century.

## 5. References

- A. Einstein, Ann. Phys. 49, 769 (1916).
- G. D. Birkhoff, Relativity and Modern Physics (Harvard, 1923).
- J. G. Hills, Astrophys. J. 267, 322 (1983).
- S. C. Vila, Astrophys. J. 236, 645 (1980).
- G. Hobbs et al., Mon. Not. R. Astron. Soc. 360, 974 (2005).
- C. M. Will, Living Rev. Relativity 17, 4 (2014).
