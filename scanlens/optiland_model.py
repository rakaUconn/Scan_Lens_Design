"""Optiland model of the scan lens (https://github.com/optiland/optiland).

The scan pivot (galvo) is the aperture stop at surface 1. The Gaussian input beam is modelled
with an entrance-pupil diameter of 2x the 1/e^2 beam diameter plus Gaussian apodization:
    I(rho) = exp(-rho^2 / (2 sigma^2)),  rho = r / (EPD/2)
For I = exp(-2 r^2 / w^2) with w = D/2 and EPD = 2D  ->  sigma = 0.25.
"""
import numpy as np
from optiland import optic as _optic
from optiland.apodization import GaussianApodization
from optiland.physical_apertures import RadialAperture

from .catalog import element_surfaces

# Human cornea (Navarro / Liou-Brennan-type values): anterior R 7.8 mm, central thickness 0.55 mm,
# posterior R 6.5 mm, n_d = 1.376, Abbe number ~56. The endothelium lies on the posterior surface.
CORNEA = dict(R_ant=7.8, thickness=0.55, R_post=6.5, n_d=1.376, abbe=56.0)


def cornea_material():
    from optiland.materials import AbbeMaterial
    return AbbeMaterial(CORNEA["n_d"], CORNEA["abbe"], model="buchdahl")

EPD_FACTOR = 2.0          # EPD = 2 x (1/e^2 beam diameter) -> pupil sampled out to 2w
GAUSS_SIGMA = 0.25        # matches exp(-2 r^2 / w^2) for EPD_FACTOR = 2

_GLASS_CATALOG = {"E-BAF11": "hikari"}


def _mat(name):
    if name.upper() == "AIR":
        return "air"
    return (name, _GLASS_CATALOG.get(name, "schott"))


def build_optic(elements, x, beam_diameter=4.0, field_max_deg=3.0,
                wavelengths=(0.500, 0.5625, 0.625, 0.6875, 0.750), primary=0.625,
                field_fracs=(0.0, 0.7, 1.0), apodize=True, name="OCT scan lens", cornea=False):
    """elements: [(catalog_name, reversed)], x = [d_pupil, gaps..., working_distance].
    cornea=True: x[-1] is the lens -> corneal apex distance and the image surface is the
    corneal endothelium (posterior corneal surface, curved)."""
    n = len(elements)
    op = _optic.Optic(name=name)
    op.surfaces.add(index=0, radius=np.inf, thickness=np.inf)
    op.surfaces.add(index=1, radius=np.inf, thickness=float(x[0]), is_stop=True,
                    comment="scan pivot (galvo)")
    k = 2
    for e, (lens, rev) in enumerate(elements):
        surfs, dia = element_surfaces(lens, rev)
        for j, (R, t, mat) in enumerate(surfs):
            if t is None:      # last surface of element -> air gap or working distance
                t = float(x[1 + e]) if e < n - 1 else float(x[n])
            op.surfaces.add(index=k, radius=float(R), thickness=float(t), material=_mat(mat),
                            aperture=RadialAperture(r_max=0.45 * dia),
                            comment=f"{lens}{' (rev)' if rev else ''} S{j+1}")
            k += 1
    if cornea:
        op.surfaces.add(index=k, radius=CORNEA["R_ant"], thickness=CORNEA["thickness"],
                        material=cornea_material(), comment="cornea anterior surface")
        k += 1
        op.surfaces.add(index=k, radius=CORNEA["R_post"], comment="endothelium (posterior cornea)")
    else:
        op.surfaces.add(index=k, comment="image / sample")
    op.set_aperture(aperture_type="EPD", value=EPD_FACTOR * beam_diameter)
    if apodize:
        op.set_apodization(GaussianApodization(sigma=GAUSS_SIGMA))
    op.set_field_type("angle")
    for f in field_fracs:
        op.add_field(y=float(f * field_max_deg))
    for wl in wavelengths:
        op.add_wavelength(value=float(wl), is_primary=abs(wl - primary) < 1e-9)
    return op


def lens_surface_range(op):
    """Indices of the first lens surface and the image surface."""
    return 2, op.surface_group.num_surfaces - 1
