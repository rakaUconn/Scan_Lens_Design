"""Performance analysis of the scan lens, computed with Optiland.

* Broadband Gaussian-beam PSF: Optiland FFTPSF (Gaussian apodization, centroid reference) per
  wavelength, shifted by that wavelength's centroid (lateral colour) and resampled onto a common
  grid, then summed (flat source spectrum 500-750 nm).
* Spot size metric: 1/e^2 diameter = diameter of the circle with the same area as the region
  where I >= I_peak / e^2.
* Also: spot diagrams, RMS wavefront vs field, chromatic focal shift, lateral colour,
  telecentricity (chief-ray angle), f-theta distortion, field curvature.
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import map_coordinates

import optiland.backend as be
from optiland.psf import FFTPSF
from optiland.analysis import SpotDiagram, FieldCurvature
from optiland.optimization.operand.operand import operand_registry

from .optiland_model import build_optic
from .design import SPEC, d_diff_um

PSF_RAYS, PSF_GRID = 48, 384
WLS_FINE = list(np.round(np.linspace(0.500, 0.750, 11), 4))


# ------------------------------------------------------------------ helpers
def _np(a):
    return np.asarray(be.to_numpy(a), dtype=float)


def centroid_um(op, Hy, wl, num_rays=16):
    r = op.trace(0.0, Hy, wl, num_rays, "hexapolar")
    x, y, i = _np(r.x), _np(r.y), _np(r.i)
    m = np.isfinite(y) & (i > 0)
    return np.average(x[m], weights=i[m]) * 1e3, np.average(y[m], weights=i[m]) * 1e3


def e2_diameter(I, dx, soft=True):
    """Equivalent 1/e^2 diameter [um]; soft threshold gives a smooth value for optimisation."""
    I = I / I.max()
    thr = np.exp(-2)
    if soft:
        frac = np.clip((I - thr) / (0.1 * thr) + 0.5, 0, 1)
        area = frac.sum() * dx * dx
    else:
        area = (I >= thr).sum() * dx * dx
    return 2 * np.sqrt(area / np.pi)


def mono_psf(op, Hy, wl, num_rays=PSF_RAYS, grid=PSF_GRID):
    """Returns PSF intensity (peak-normalised), pixel size [um], Strehl vs ideal Gaussian beam."""
    p = FFTPSF(op, (0.0, Hy), wl, num_rays=num_rays, grid_size=grid, strategy="centroid_sphere")
    I = _np(p.psf)
    pupil = np.asarray(be.to_numpy(p.pupils[0]))
    ideal_peak = np.abs(pupil).sum() ** 2                  # same Gaussian amplitude, zero phase
    raw_peak = I.max() * float(_np(p._get_normalization())) / 100
    strehl = raw_peak / ideal_peak
    fno = float(_np(p._get_working_FNO()))
    dx = wl * fno / (grid / (num_rays - 1))
    return I / I.max(), dx, strehl


def broadband_psf(op, Hy, wls, half_width_um=8.0, n_common=321, ref_wl=0.625, num_rays=PSF_RAYS, grid=PSF_GRID):
    common = np.linspace(-half_width_um, half_width_um, n_common)
    dc = common[1] - common[0]
    cx0, cy0 = centroid_um(op, Hy, ref_wl)
    acc = np.zeros((n_common, n_common))
    per_wl = []
    for wl in wls:
        I, dx, S = mono_psf(op, Hy, wl, num_rays, grid)
        cx, cy = centroid_um(op, Hy, wl)
        n = I.shape[0]
        c0 = n // 2
        gy, gx = np.meshgrid(common, common, indexing="ij")
        # PSF is centred on this wavelength's centroid -> shift into the common frame
        iy = (gy - (cy - cy0)) / dx + c0
        ix = (gx - (cx - cx0)) / dx + c0
        Ir = map_coordinates(I, [iy, ix], order=1, mode="constant")
        acc += Ir / max(Ir.sum(), 1e-30)
        per_wl.append(dict(wl_nm=wl * 1e3, d_e2_um=float(e2_diameter(I, dx, soft=False)),
                           strehl=float(S), lat_color_um=float(cy - cy0), dx_um=float(dx)))
    return common, acc / acc.max(), dc, per_wl


def oct_psf_diameter(optic, Hy):
    """Optimisation operand: broadband 1/e^2 PSF diameter [um] at normalised field Hy."""
    try:
        _, I, dc, _ = broadband_psf(optic, Hy, SPEC["wavelengths_um"], num_rays=32, grid=256)
        v = float(e2_diameter(I, dc, soft=True))
        return v if np.isfinite(v) else 1e3
    except Exception:
        return 1e3


operand_registry.register("oct_psf_diameter", oct_psf_diameter, overwrite=True)


def paraxial_focus_z(op, wl, Py=0.01):
    """Axial position (relative to last lens vertex) where a near-axis marginal ray crosses the axis."""
    r = op.trace_generic(0.0, 0.0, 0.0, Py, wl)
    # rays end on the image plane: back-propagate to y = 0
    y, z, M, N = _np(r.y)[0], _np(r.z)[0], _np(r.M)[0], _np(r.N)[0]
    return z - y * N / M


def _semi(s):
    ap = getattr(s, "aperture", None)
    r = getattr(ap, "r_max", None) if ap is not None else None
    if r is None:
        r = getattr(s, "semi_aperture", None)
    return "–" if r is None else f"{float(_np(r)):.2f}"


def chief(op, Hy, wl):
    r = op.trace_generic(0.0, Hy, 0.0, 0.0, wl)
    return _np(r.y)[0], np.degrees(np.arctan2(_np(r.M)[0], _np(r.N)[0]))


# ------------------------------------------------------------------ full report
def report(elements, x, outdir, title="design", field_max_deg=None, n_fields=6):
    os.makedirs(outdir, exist_ok=True)
    D = SPEC["beam_diameter_mm"]
    wl0 = SPEC["ref_wl_um"]
    op_tmp = build_optic(elements, x, D, 1.0, WLS_FINE, wl0)
    efl = float(_np(op_tmp.paraxial.f2()))
    if field_max_deg is None:
        field_max_deg = float(np.degrees(np.arctan(SPEC["field_half_mm"] / efl)))
    Hys = np.linspace(0, 1, n_fields)
    op = build_optic(elements, x, D, field_max_deg, WLS_FINE, wl0, field_fracs=tuple(Hys))
    res = dict(elements=[list(e) for e in elements], x=[float(v) for v in x], efl_mm=efl,
               field_max_deg=field_max_deg, engine="optiland", rows=[], broadband=[])

    psf_imgs = {}
    for Hy in Hys:
        common, Ibb, dc, per = broadband_psf(op, Hy, WLS_FINE)
        y_mm = centroid_um(op, Hy, wl0)[1] / 1e3
        for p in per:
            p.update(field_deg=float(Hy * field_max_deg), y_mm=float(y_mm),
                     d_diff_um=float(d_diff_um(p["wl_nm"] / 1e3, efl, D)))
            res["rows"].append(p)
        dbb = float(e2_diameter(Ibb, dc, soft=False))
        res["broadband"].append(dict(field_deg=float(Hy * field_max_deg), y_mm=float(y_mm),
                                     d_e2_um=dbb, fwhm_equiv_um=float(dbb * np.sqrt(np.log(2) / 2))))
        psf_imgs[round(Hy, 3)] = (common, Ibb)

    # focal shift, telecentricity, distortion
    z_ref = paraxial_focus_z(op, wl0)
    res["focal_shift_um"] = {f"{w*1e3:.0f}": float((paraxial_focus_z(op, w) - z_ref) * 1e3) for w in WLS_FINE}
    w0 = 0.5 * d_diff_um(wl0, efl, D)
    res["rayleigh_range_um"] = float(np.pi * w0 ** 2 / wl0)
    ych, cra = zip(*[chief(op, Hy, wl0) for Hy in Hys])
    th = np.radians(Hys * field_max_deg)
    ft = efl * th
    res["chief_ray_angle_deg"] = [float(c) for c in cra]
    res["y_chief_mm"] = [float(v) for v in ych]
    res["ftheta_distortion_pct"] = [float(100 * (y - f) / f) if f > 0 else 0.0 for y, f in zip(ych, ft)]
    surf = op.surface_group
    res["pupil_distance_mm"] = float(x[0])
    res["working_distance_mm"] = float(x[-1])
    res["total_track_mm"] = float(_np(op.total_track))

    # ============================================================ figures
    fig, ax = op.draw(fields=[(0, 0), (0, 0.7), (0, 1.0)], num_rays=7, figsize=(11, 4),
                      title=f"{title}: EFL {efl:.2f} mm, galvo->lens {x[0]:.1f} mm, WD {x[-1]:.2f} mm")
    fig.savefig(f"{outdir}/layout.png", dpi=150, bbox_inches="tight"); plt.close(fig)

    op_sd = build_optic(elements, x, D, field_max_deg, [0.5, 0.5625, 0.625, 0.6875, 0.75], wl0)
    sd = SpotDiagram(op_sd, num_rings=12)
    fig, _ = sd.view(figsize=(12, 4), show=False)
    fig.savefig(f"{outdir}/spot_diagram.png", dpi=150, bbox_inches="tight"); plt.close(fig)

    try:
        fc = FieldCurvature(op_sd)
        fig, _ = fc.view(show=False) if "show" in fc.view.__code__.co_varnames else fc.view()
        fig.savefig(f"{outdir}/field_curvature.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    except Exception as e:
        print("field curvature plot skipped:", e)

    rows = res["rows"]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
    for wl in WLS_FINE[::2]:
        c = plt.cm.turbo((wl - 0.48) / 0.3)
        rr = [r for r in rows if abs(r["wl_nm"] - wl * 1e3) < 1e-6]
        axs[0].plot([r["y_mm"] for r in rr], [r["d_e2_um"] for r in rr], "o-", color=c, label=f"{wl*1e3:.0f} nm")
        axs[0].axhline(d_diff_um(wl, efl, D), color=c, ls=":", lw=0.8)
        axs[1].plot([r["y_mm"] for r in rr], [r["strehl"] for r in rr], "o-", color=c, label=f"{wl*1e3:.0f} nm")
    bb = res["broadband"]
    axs[0].plot([b["y_mm"] for b in bb], [b["d_e2_um"] for b in bb], "ks-", lw=2, label="500-750 nm broadband")
    axs[0].axhline(SPEC["spot_target_um"], color="r", ls="--", label="3 µm target")
    axs[0].set_xlabel("scan position at sample [mm]"); axs[0].set_ylabel("PSF 1/e² diameter [µm]")
    axs[0].set_title("Optiland FFT PSF, Ø4 mm Gaussian beam (dotted: diffraction limit)", fontsize=10)
    axs[0].legend(fontsize=7); axs[0].grid(alpha=0.3)
    axs[1].axhline(0.8, color="r", ls="--", label="0.8"); axs[1].set_ylim(0, 1.05)
    axs[1].set_xlabel("scan position at sample [mm]"); axs[1].set_ylabel("Strehl (vs ideal Gaussian beam)")
    axs[1].legend(fontsize=7); axs[1].grid(alpha=0.3); axs[1].set_title("Strehl ratio at common focal plane", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{outdir}/spot_size_vs_field.png", dpi=150); plt.close(fig)

    keys = [round(h, 3) for h in (Hys[0], Hys[len(Hys) // 2], Hys[-1])]
    fig, axs = plt.subplots(1, 3, figsize=(12, 4))
    for a, k in zip(axs, keys):
        common, I = psf_imgs[k]
        a.imshow(I, extent=[common[0], common[-1], common[0], common[-1]], origin="lower", cmap="inferno")
        a.contour(common, common, I, [np.exp(-2)], colors="c", linewidths=0.8)
        a.add_patch(plt.Circle((0, 0), 1.5, fill=False, color="w", ls="--", lw=0.8))
        a.set_xlim(-5, 5); a.set_ylim(-5, 5)
        a.set_title(f"broadband PSF, field {k*field_max_deg:.2f}° (cyan 1/e², white Ø3 µm)", fontsize=8)
        a.set_xlabel("x [µm]"); a.set_ylabel("y [µm]")
    fig.tight_layout(); fig.savefig(f"{outdir}/broadband_psf.png", dpi=150); plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
    zR = res["rayleigh_range_um"]
    axs[0].plot([float(k) for k in res["focal_shift_um"]], list(res["focal_shift_um"].values()), "o-")
    axs[0].axhspan(-zR, zR, color="g", alpha=0.15, label=f"± Rayleigh range ({zR:.0f} µm)")
    axs[0].set_xlabel("wavelength [nm]"); axs[0].set_ylabel("focal shift [µm]"); axs[0].legend(fontsize=8)
    axs[0].set_title("Chromatic focal shift"); axs[0].grid(alpha=0.3)
    for wl in WLS_FINE[::2]:
        rr = [r for r in rows if abs(r["wl_nm"] - wl * 1e3) < 1e-6]
        axs[1].plot([r["y_mm"] for r in rr], [r["lat_color_um"] for r in rr], "o-",
                    color=plt.cm.turbo((wl - 0.48) / 0.3), label=f"{wl*1e3:.0f} nm")
    axs[1].set_xlabel("scan position [mm]"); axs[1].set_ylabel("centroid shift vs 625 nm [µm]")
    axs[1].set_title("Lateral colour"); axs[1].legend(fontsize=7); axs[1].grid(alpha=0.3)
    axs[2].plot(ych, cra, "o-", label="chief-ray angle [°]")
    axs[2].plot(ych, res["ftheta_distortion_pct"], "s-", label="f-θ distortion [%]")
    axs[2].set_xlabel("scan position [mm]"); axs[2].legend(fontsize=8); axs[2].grid(alpha=0.3)
    axs[2].set_title("Telecentricity & f-θ linearity")
    fig.tight_layout(); fig.savefig(f"{outdir}/color_telecentricity.png", dpi=150); plt.close(fig)

    # prescription table from the optiland model
    lines = ["| # | Surface | Radius [mm] | Thickness [mm] | Material | Semi-aperture [mm] |",
             "|---|---|---|---|---|---|"]
    for i, s in enumerate(surf.surfaces):
        if i == 0:
            continue
        R = float(_np(s.geometry.radius))
        t = float(_np(surf.get_thickness(i))[0]) if i < surf.num_surfaces - 1 else float("nan")
        mat = s.material_post
        mname = getattr(mat, "name", None) or ("AIR" if float(_np(mat.n(0.625))) == 1.0 else str(mat))
        lines.append(f"| {i} | {s.comment} | {'∞' if np.isinf(R) else f'{R:.3f}'} | "
                     f"{'–' if np.isnan(t) else f'{t:.3f}'} | {mname} | {_semi(s)} |")
    open(f"{outdir}/prescription.md", "w").write("\n".join(lines) + "\n")
    json.dump(res, open(f"{outdir}/performance.json", "w"), indent=1)
    return res


def summary_text(res):
    bb = res["broadband"]
    rows = res["rows"]
    fs = res["focal_shift_um"].values()
    return "\n".join([
        f"EFL {res['efl_mm']:.3f} mm | galvo->lens {res['pupil_distance_mm']:.2f} mm | WD {res['working_distance_mm']:.2f} mm"
        f" | scan ±{res['field_max_deg']:.2f}°",
        "broadband 1/e² PSF diameter vs field: " + ", ".join(f"{b['d_e2_um']:.2f}" for b in bb) + " µm",
        "  (FWHM-equivalent: " + ", ".join(f"{b['fwhm_equiv_um']:.2f}" for b in bb) + " µm)",
        f"worst mono spot {max(r['d_e2_um'] for r in rows):.2f} µm, min Strehl {min(r['strehl'] for r in rows):.3f}",
        f"focal shift 500-750 nm: {min(fs):.1f} .. {max(fs):.1f} µm (Rayleigh range {res['rayleigh_range_um']:.1f} µm)",
        f"max lateral colour {max(abs(r['lat_color_um']) for r in rows):.2f} µm | CRA at edge {res['chief_ray_angle_deg'][-1]:.2f}°"
        f" | f-θ distortion at edge {res['ftheta_distortion_pct'][-1]:.3f}%"])


if __name__ == "__main__":
    d = json.load(open(sys.argv[1]))
    r = d["results"][int(sys.argv[2])] if "results" in d else d
    out = sys.argv[3]
    res = report([tuple(e) for e in r["elements"]], r["x"], out, title=os.path.basename(out))
    print(summary_text(res))
