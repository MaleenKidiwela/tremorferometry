#!/usr/bin/env python
"""Phase B: half-space diffusion coda-wave sensitivity kernel for a DEEP source + SURFACE receiver
(theory/01_kernel.md, eq. 13-14). One 1-D quadrature per model point; no wavefield simulation.

K(x; s, r, t) = (1/N(t)) integral_0^t du/[u(t-u)]^1.5
                  * [ exp(-(rho_s^2+(z-d)^2)/(a u)) + exp(-(rho_s^2+(z+d)^2)/(a u)) ]   <- source leg + image
                  * exp(-(rho_r^2+z^2)/(a (t-u)))                                         <- receiver leg (factor-2 surface)
with a = 4D, rho_s = |x_h - s_h|, rho_r = |x_h - r_h|, z = model depth, d = source depth,
and normalizer N(t) = exp(-(|r_h-s_h|^2 + d^2)/(a t)) / sqrt(t)  (eq. 14), giving integral_V K dV = 1.

Units: lengths in km, speeds km/s, D=beta*l*/3 in km^2/s, t in s, K in 1/km^3.
Intrinsic Q cancels in the kernel SHAPE (theory sec 3), so it does not enter here.
"""
import numpy as np


def kernel_diffusion(xh, z, s_h, d, r_h, D, t, nu=400):
    """Diffusion sensitivity kernel at model points (xh, z) for source epicentre s_h at depth d,
    surface receiver r_h, diffusivity D, lapse time t.

    xh : (N,2) model horizontal coords [km]   z : (N,) model depth [km], z>=0 down
    s_h: (2,) source epicentre [km]           d : source depth [km]
    r_h: (2,) receiver location [km]          D : diffusivity [km^2/s]   t: lapse time [s]
    returns K : (N,) kernel values [1/km^3]   (un-normalised by cell volume; multiply by dV to weight)
    """
    xh = np.atleast_2d(xh); z = np.atleast_1d(z).astype(float)
    a = 4.0 * D
    rho_s2 = ((xh - np.asarray(s_h)) ** 2).sum(1)   # (N,)
    rho_r2 = ((xh - np.asarray(r_h)) ** 2).sum(1)
    # midpoint rule in u on (0,t); the exp(-A/u)exp(-B/(t-u)) factors temper the endpoint singularities
    u = (np.arange(nu) + 0.5) / nu * t              # (nu,)
    du = t / nu
    tu = t - u
    inv_au = 1.0 / (a * u); inv_atu = 1.0 / (a * tu)
    src = (np.exp(-(rho_s2[:, None] + (z[:, None] - d) ** 2) * inv_au[None, :]) +
           np.exp(-(rho_s2[:, None] + (z[:, None] + d) ** 2) * inv_au[None, :]))
    rec = np.exp(-(rho_r2[:, None] + z[:, None] ** 2) * inv_atu[None, :])
    integrand = src * rec / (u[None, :] * tu[None, :]) ** 1.5
    K = integrand.sum(1) * du
    # Normaliser: the theory note's eq.14 (N = exp(...)/sqrt(t)) dropped the (pi*a)^1.5 prefactor that
    # the propagators (eq.12) carry. Restoring it makes integral_V K dV = 1 exactly (Chapman-Kolmogorov:
    # integral_V P(s,x,u)P(x,r,t-u) dV = P(s,r,t) => integral_V numerator dV = t P(s,r,t) = denominator).
    N = (np.pi * a) ** 1.5 * np.exp(-(((np.asarray(r_h) - np.asarray(s_h)) ** 2).sum() + d ** 2) / (a * t)) / np.sqrt(t)
    return K / N


def kernel_singlescatter(xh, z, s_h, d, r_h, beta, w1, w2, ell=40.0):
    """Single-scattering (early-coda) sensitivity kernel -- the prolate-ellipsoid shell with foci at the
    source and receiver (theory/01 eq.15). Appropriate for our 1-4 s coda window (early coda, NOT the
    diffusive late coda). Energy at lapse time t past origin scattered once has path length R_s+R_r = beta*t;
    a coda WINDOW [w1,w2] seconds past the DIRECT arrival -> an ellipsoidal shell of total path in
    [direct + beta*w1, direct + beta*w2]. Free surface -> add the image-source ellipsoid (depth -d).
    Returns K on the model points, normalised so it sums to 1 over the points passed (on-fault restriction;
    the off-fault remainder is the captured-sensitivity deficit eta<1)."""
    xh = np.atleast_2d(xh); z = np.atleast_1d(z).astype(float)
    rs_h2 = ((xh - np.asarray(s_h)) ** 2).sum(1)
    rr2 = ((xh - np.asarray(r_h)) ** 2).sum(1) + z ** 2          # receiver at surface z=0
    Rr = np.sqrt(rr2)
    direct = np.sqrt(((np.asarray(r_h) - np.asarray(s_h)) ** 2).sum() + d ** 2)
    lo, hi = direct + beta * w1, direct + beta * w2
    K = np.zeros(len(z))
    for zsrc in (d, -d):                                         # real source + free-surface image
        Rs = np.sqrt(rs_h2 + (z - zsrc) ** 2)
        path = Rs + Rr
        m = (path >= lo) & (path <= hi)
        K[m] += np.exp(-path[m] / ell) / (Rs[m] ** 2 * Rr[m] ** 2)
    tot = K.sum()
    return K / tot if tot > 0 else K


def _selftest():
    """Verify integral_V K dV = 1 (=> homogeneous perturbation gives dtau/t = -dbeta/beta)."""
    D = 50.0; t = 8.0; d = 35.0                      # km^2/s, s, km : crustal-ish, mid-coda, LFE depth
    s_h = np.array([0.0, 0.0]); r_h = np.array([40.0, 0.0])
    # 3-D grid (half-space z>=0) wide enough to contain the kernel mass
    L = 160.0; nz = 80; nxy = 81
    xs = np.linspace(-L, L, nxy); ys = np.linspace(-L, L, nxy); zs = np.linspace(0.5, L, nz)
    dx = xs[1]-xs[0]; dy = ys[1]-ys[0]; dz = zs[1]-zs[0]; dV = dx*dy*dz
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    xh = np.column_stack([X.ravel(), Y.ravel()]); z = Z.ravel()
    K = kernel_diffusion(xh, z, s_h, d, r_h, D, t)
    mass = (K * dV).sum()
    # where is the sensitivity? report the two lobes (near source vs near receiver)
    near_src = (((xh - s_h) ** 2).sum(1) < 30 ** 2) & (np.abs(z - d) < 30)
    near_rec = (((xh - r_h) ** 2).sum(1) < 30 ** 2) & (z < 30)
    fs = (K[near_src] * dV).sum() / mass; fr = (K[near_rec] * dV).sum() / mass
    print(f'integral K dV = {mass:.3f}  (target 1.0; grid-truncation/discretisation error expected)')
    print(f'  deep-lobe (within 30km of source@{d}km): {100*fs:.0f}% of mass')
    print(f'  surface-lobe (within 30km of receiver, z<30): {100*fr:.0f}% of mass')
    assert 0.8 < mass < 1.2, 'normalisation off -- check grid extent / nu'
    print('  OK: kernel normalised and two-lobed as expected.')


if __name__ == '__main__':
    _selftest()
