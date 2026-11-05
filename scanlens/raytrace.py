"""Minimal sequential real-ray tracer for rotationally symmetric spherical systems.

Conventions
-----------
* z is the optical axis, light travels toward +z. Units: mm.
* The entrance pupil (scan pivot = galvo mirror position) is the plane z = 0.
* Object at infinity: each field is a collimated beam tilted by angle theta about the pupil centre.
* A system is a list of `Surface`s (absolute vertex positions) plus an image plane position.
"""
from dataclasses import dataclass, field
import numpy as np

from .glass import index
from .catalog import element_surfaces


@dataclass
class Surface:
    z: float            # vertex position [mm]
    R: float            # radius of curvature [mm] (np.inf = plane)
    mat: str            # material AFTER this surface
    semi_ap: float      # clear semi-aperture [mm]
    label: str = ""

    @property
    def c(self):
        return 0.0 if np.isinf(self.R) else 1.0 / self.R


@dataclass
class System:
    surfaces: list
    z_image: float
    pupil_semi: float = 2.0        # 1/e^2 beam radius at the pupil [mm]
    elements: list = field(default_factory=list)  # (name, reversed) for reporting
    n_surf: list = field(default_factory=list)    # surfaces per element

    def copy_with_image(self, z_image):
        return System(self.surfaces, z_image, self.pupil_semi, self.elements, self.n_surf)


def build_system(elements, d_pupil, gaps, z_image_from_last=None, pupil_semi=2.0):
    """elements: list of (catalog_name, reversed_bool); d_pupil: pupil -> first vertex;
    gaps: air gaps between consecutive elements; z_image_from_last: last vertex -> image."""
    surfs, z, nsurf = [], d_pupil, []
    for k, (name, rev) in enumerate(elements):
        s, dia = element_surfaces(name, rev)
        nsurf.append(len(s))
        for j, (R, t, mat) in enumerate(s):
            surfs.append(Surface(z, R, mat, 0.5 * dia * 0.9, f"{name}{'(rev)' if rev else ''} S{j+1}"))
            if t is not None:
                z += t
        if k < len(elements) - 1:
            z += gaps[k]
    z_last = z
    zi = z_last + (z_image_from_last if z_image_from_last is not None else 0.0)
    return System(surfs, zi, pupil_semi, list(elements), nsurf)


def _intersect(P, D, s):
    """Intersect rays (P, D unit) with spherical surface s. Returns t and surface normal."""
    c = s.c
    Pz = P[:, 2] - s.z
    if c == 0.0:
        t = -Pz / D[:, 2]
        N = np.zeros_like(P); N[:, 2] = 1.0
        return t, N
    # sphere: x^2 + y^2 + (z - R)^2 = R^2 in local coords
    # Solve c(x^2+y^2+z^2) - 2z = 0 along P + tD
    A = c * np.sum(D * D, axis=1)
    Pl = P.copy(); Pl[:, 2] = Pz
    B = 2 * (c * np.sum(Pl * D, axis=1) - D[:, 2])
    C = c * np.sum(Pl * Pl, axis=1) - 2 * Pz
    disc = B * B - 4 * A * C
    disc = np.where(disc < 0, np.nan, disc)
    # stable root closest to the vertex
    t = 2 * C / (-B + np.sign(-B + 1e-300) * 0 - np.sqrt(disc) * np.sign(B + 1e-300) * -1) if False else None
    sq = np.sqrt(disc)
    q = -0.5 * (B + np.sign(B) * sq)
    t1 = C / q
    t2 = q / A
    # choose the root that lands on the cap near the vertex (|z_local| small)
    z1 = Pz + t1 * D[:, 2]
    z2 = Pz + t2 * D[:, 2]
    t = np.where(np.abs(z1) < np.abs(z2), t1, t2)
    X = Pl + t[:, None] * D
    N = np.stack([-c * X[:, 0], -c * X[:, 1], 1 - c * X[:, 2]], axis=1)
    N /= np.linalg.norm(N, axis=1)[:, None]
    return t, N


def _refract(D, N, n1, n2):
    mu = n1 / n2
    cosi = np.sum(D * N, axis=1)
    N = np.where(cosi[:, None] < 0, -N, N)
    cosi = np.abs(cosi)
    k = 1 - mu * mu * (1 - cosi * cosi)
    k = np.where(k < 0, np.nan, k)   # TIR -> nan
    Dn = mu * D + (np.sqrt(k) - mu * cosi)[:, None] * N
    return Dn / np.linalg.norm(Dn, axis=1)[:, None]


def trace(system, P0, D0, wl_um, to_image=True, check_aperture=True):
    """Trace rays. Returns final positions, directions, optical path length (from P0 plane),
    validity mask, and list of per-surface hit points."""
    P, D = P0.astype(float).copy(), D0.astype(float).copy()
    D /= np.linalg.norm(D, axis=1)[:, None]
    opl = np.zeros(len(P))
    n1 = 1.0
    ok = np.ones(len(P), bool)
    hits = []
    for s in system.surfaces:
        t, N = _intersect(P, D, s)
        P = P + t[:, None] * D
        opl += n1 * t
        hits.append(P.copy())
        if check_aperture:
            ok &= np.hypot(P[:, 0], P[:, 1]) <= s.semi_ap
        n2 = index(s.mat, wl_um)
        D = _refract(D, N, n1, n2)
        n1 = n2
    if to_image:
        t = (system.z_image - P[:, 2]) / D[:, 2]
        P = P + t[:, None] * D
        opl += n1 * t
        hits.append(P.copy())
    ok &= np.all(np.isfinite(P), axis=1)
    return P, D, opl, ok, hits


# ---------------------------------------------------------------- paraxial helpers
def paraxial(system, wl_um=0.625):
    """EFL, BFL (last vertex -> paraxial focus), front focal distance (pupil plane is at z=0)."""
    h = 1e-4
    P0 = np.array([[0, h, -1.0]]); D0 = np.array([[0, 0, 1.0]])
    P, D, *_ = trace(system, P0, D0, wl_um, to_image=False, check_aperture=False)
    u = D[0, 1] / D[0, 2]
    efl = -h / u
    z_focus = P[0, 2] - P[0, 1] / u
    bfl = z_focus - system.surfaces[-1].z
    return efl, bfl, z_focus


def pupil_grid(n_rings=12, semi=2.0, extent=1.5):
    """Square grid of pupil coordinates covering `extent` x 1/e^2 radius (Gaussian beam)."""
    r = semi * extent
    x = np.linspace(-r, r, 2 * n_rings + 1)
    X, Y = np.meshgrid(x, x)
    m = X ** 2 + Y ** 2 <= r * r
    px, py = X[m], Y[m]
    w = np.exp(-2 * (px ** 2 + py ** 2) / semi ** 2)  # Gaussian intensity weight
    return px, py, w


def launch(theta_deg, px, py, z0=-5.0):
    """Collimated beam at field angle theta (in y-z plane) pivoting about pupil centre (0,0,0).
    Rays start on the plane perpendicular to the beam through z0 so OPL starts on a flat wavefront."""
    th = np.radians(theta_deg)
    d = np.array([0.0, np.sin(th), np.cos(th)])
    # pupil-plane points (z=0) -> back-propagate to a plane wavefront through point d*z0
    P = np.stack([px, py, np.zeros_like(px)], axis=1)
    # move each point back along d so that all start on the plane d.(X - d*z0) = 0
    s = (P @ d) - z0
    P0 = P - s[:, None] * d
    D0 = np.tile(d, (len(px), 1))
    return P0, D0


def chief_ray(system, theta_deg, wl_um):
    P0, D0 = launch(theta_deg, np.array([0.0]), np.array([0.0]))
    P, D, *_ = trace(system, P0, D0, wl_um, check_aperture=False)
    return P[0], D[0]


def wavefront(system, theta_deg, wl_um, grid, ref_point=None):
    """Gaussian-weighted RMS wavefront error [waves] (piston & tilt removed), plus the
    weighted centroid of the ray intercepts on the image plane and per-ray data."""
    px, py, w = grid
    P0, D0 = launch(theta_deg, px, py)
    P, D, opl, ok, _ = trace(system, P0, D0, wl_um)
    if ok.sum() < 0.8 * len(ok):
        return np.inf, np.array([np.nan, np.nan]), None
    w = w * ok
    Pv, Dv, oplv, wv = P[ok], D[ok], opl[ok], w[ok]
    cen = np.average(Pv[:, :2], axis=0, weights=wv)
    R = np.array([cen[0], cen[1], system.z_image]) if ref_point is None else ref_point
    n_img = index(system.surfaces[-1].mat, wl_um)
    W = oplv + n_img * np.sum((R - Pv) * Dv, axis=1)          # OPL to reference sphere
    A = np.stack([np.ones_like(W), px[ok], py[ok]], axis=1)   # remove piston + tilt
    sw = np.sqrt(wv)
    coef, *_ = np.linalg.lstsq(A * sw[:, None], W * sw, rcond=None)
    res = W - A @ coef
    rms = np.sqrt(np.average(res ** 2, weights=wv)) / (wl_um * 1e-3)
    return rms, cen, dict(P=Pv, W=res / (wl_um * 1e-3), px=px[ok], py=py[ok], w=wv)
