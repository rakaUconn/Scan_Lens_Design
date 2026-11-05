"""Final optimisation of a scan-lens candidate with Optiland's OptimizationProblem.

Variables : stop->first lens distance (galvo clearance), air gaps, working distance.
Operands  : 'oct_psf_diameter' per field = broadband (500-750 nm) 1/e^2 diameter of the Optiland
            FFT PSF of the Gaussian beam, incl. lateral colour (see optiland_analysis.py);
            'edge_clearance' between lens elements. ('oct_spot' = faster Strehl-based estimate.)
"""
import json, sys
import numpy as np
import optiland.backend as be
from optiland.optimization import OptimizationProblem, OptimizerGeneric
from optiland.optimization.operand.operand import operand_registry
from optiland.wavefront import Wavefront

from .optiland_model import build_optic
from .catalog import element_surfaces
from .design import SPEC, d_diff_um, sag
from . import optiland_analysis  # noqa: F401  (registers the 'oct_psf_diameter' operand)

WLS = SPEC["wavelengths_um"]


def gauss_rms_waves(optic, Hy, wl, num_rays=24):
    """Intensity (Gaussian beam) weighted RMS wavefront error [waves], piston/tilt removed."""
    wf = Wavefront(optic, fields=[(0.0, Hy)], wavelengths=[wl], num_rays=num_rays,
                   distribution="uniform", strategy="centroid_sphere", remove_tilt=True)
    d = wf.get_data((0.0, Hy), wl)
    I = be.to_numpy(d.intensity); o = be.to_numpy(d.opd)
    m = np.isfinite(o) & (I > 0)
    if m.sum() < 10:
        return np.inf
    o = o[m] - np.average(o[m], weights=I[m])
    return float(np.sqrt(np.average(o ** 2, weights=I[m])))


def centroid_y_um(optic, Hy, wl, num_rays=12):
    rays = optic.trace(0.0, Hy, wl, num_rays, "hexapolar")
    y = be.to_numpy(rays.y); i = be.to_numpy(rays.i)
    m = np.isfinite(y) & (i > 0)
    return float(np.average(y[m], weights=i[m]) * 1e3)


def oct_spot(optic, Hy):
    """Estimated broadband 1/e^2 spot diameter [um] at normalised field Hy."""
    try:
        efl = float(optic.paraxial.f2())
        D = SPEC["beam_diameter_mm"]
        dmax, cens = 0.0, []
        for wl in WLS:
            rms = gauss_rms_waves(optic, Hy, wl)
            S = np.exp(-(2 * np.pi * rms) ** 2)
            dmax = max(dmax, d_diff_um(wl, efl, D) / np.sqrt(max(S, 1e-3)))
            cens.append(centroid_y_um(optic, Hy, wl))
        lat = max(cens) - min(cens)
        v = float(np.sqrt(dmax ** 2 + lat ** 2))
        return v if np.isfinite(v) else 1e3
    except Exception:
        return 1e3


def edge_clearance(optic, surf_a, surf_b, r):
    """Axial edge clearance [mm] between surface a (rear of one lens) and b (front of next)."""
    sa, sb = optic.surfaces[surf_a], optic.surfaces[surf_b]
    Ra = float(be.to_numpy(sa.geometry.radius)); Rb = float(be.to_numpy(sb.geometry.radius))
    t = float(be.to_numpy(optic.surfaces.get_thickness(surf_a))[0])
    return t - sag(Ra, r) + sag(Rb, r)


for _name, _fn in [("oct_spot", oct_spot), ("edge_clearance", edge_clearance)]:
    operand_registry.register(_name, _fn, overwrite=True)


def air_surfaces(elements):
    """Surface indices (optiland numbering) whose thickness is a free air space."""
    idx, k = [1], 2                       # 1 = stop -> first lens
    for e, el in enumerate(elements):
        ns = len(element_surfaces(*el)[0])
        k += ns
        idx.append(k - 1)                 # last surface of element -> gap / WD
    return idx


def optimise_optiland(elements, x, field_max_deg, maxiter=60, disp=False, operand="oct_psf_diameter"):
    op = build_optic(elements, x, SPEC["beam_diameter_mm"], field_max_deg, WLS, SPEC["ref_wl_um"])
    prob = OptimizationProblem()
    air = air_surfaces(elements)
    n = len(elements)
    # variables
    prob.add_variable(op, "thickness", surface_number=air[0], min_val=SPEC["min_pupil_dist_mm"], max_val=40.0)
    for s in air[1:-1]:
        prob.add_variable(op, "thickness", surface_number=s, min_val=SPEC["min_vertex_gap_mm"], max_val=20.0)
    prob.add_variable(op, "thickness", surface_number=air[-1], min_val=SPEC["min_working_dist_mm"], max_val=40.0)
    # operands: true broadband FFT-PSF 1/e^2 diameter (Optiland FFTPSF) at 3 field points
    for Hy in (0.0, 0.7, 1.0):
        prob.add_operand(operand_type=operand, target=0.0, weight=1.0,
                         input_data={"optic": op, "Hy": Hy})
    # mechanical edge clearance between consecutive elements
    for e in range(n - 1):
        sa = air[1 + e]
        dia = min(element_surfaces(*elements[e])[1], element_surfaces(*elements[e + 1])[1])
        prob.add_operand(operand_type="edge_clearance", min_val=SPEC["min_edge_gap_mm"], weight=10.0,
                         input_data={"optic": op, "surf_a": sa, "surf_b": sa + 1, "r": 0.45 * dia})
    opt = OptimizerGeneric(prob)
    res = opt.optimize(method="Nelder-Mead", maxiter=maxiter * len(air), disp=disp, tol=1e-4)
    x_new = [float(be.to_numpy(op.surfaces.get_thickness(s))[0]) for s in air]
    return op, x_new, prob


if __name__ == "__main__":
    data = json.load(open(sys.argv[1]))
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    r = data["results"][k]
    el = [tuple(e) for e in r["elements"]]
    th = float(np.degrees(np.arctan(SPEC["field_half_mm"] / r["efl"])))
    op, x, prob = optimise_optiland(el, r["x"], th, disp=True)
    prob.info()
    print("x =", x)
