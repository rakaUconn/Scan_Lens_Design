"""Detailed evaluation of a scan-lens design: FFT Gaussian-beam PSF, spot diagrams, colour,
telecentricity, distortion, layout plot and a prescription table."""
import json, sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import map_coordinates

from .raytrace import (paraxial, pupil_grid, wavefront, chief_ray, launch, trace)
from .design import SPEC, make, d_diff_um, field_angles, edge_gap
from .glass import index

WLS_FINE = np.linspace(0.500, 0.750, 11)


# ------------------------------------------------------------------ physical-optics PSF
def opd_map(system, theta, wl, n=64, extent=2.0, ref_point=None):
    """OPD [waves] on a square entrance-pupil grid (tilt NOT removed if ref_point given)."""
    w = system.pupil_semi
    x = np.linspace(-extent * w, extent * w, n)
    X, Y = np.meshgrid(x, x)
    P0, D0 = launch(theta, X.ravel(), Y.ravel())
    P, D, opl, ok, _ = trace(system, P0, D0, wl)
    R = ref_point
    n_img = index(system.surfaces[-1].mat, wl)
    W = opl + n_img * np.sum((R - P) * D, axis=1)
    W = np.where(ok, W, np.nan)
    c = n * (n // 2) + n // 2
    W = (W - np.nanmean(W[ok & (np.hypot(X.ravel(), Y.ravel()) < 0.3 * w)])) / (wl * 1e-3)
    amp = np.exp(-(X ** 2 + Y ** 2) / w ** 2) * ok.reshape(n, n)  # Gaussian field amplitude
    return x, amp, np.nan_to_num(W.reshape(n, n)), ok.reshape(n, n)


def psf(system, theta, wl, efl, ref_point, n=64, pad=8, extent=2.0):
    """Intensity PSF on image-plane grid [um] (paraxial scaling x = lambda*EFL*fx)."""
    x, amp, W, ok = opd_map(system, theta, wl, n, extent, ref_point)
    dp = x[1] - x[0]
    N = n * pad
    field = np.zeros((N, N), complex)
    field[:n, :n] = amp * np.exp(2j * np.pi * W)
    I = np.abs(np.fft.fftshift(np.fft.fft2(field))) ** 2
    dx_um = wl * efl / (N * dp)             # um per pixel (wl in um, efl/dp dimensionless)
    coords = (np.arange(N) - N // 2) * dx_um
    return coords, I / I.max() if I.max() > 0 else I


def resample(coords, I, grid):
    """Bilinear resample PSF onto common grid (um)."""
    d = coords[1] - coords[0]
    gy, gx = np.meshgrid(grid, grid, indexing="ij")
    iy = (gy - coords[0]) / d
    ix = (gx - coords[0]) / d
    return map_coordinates(I, [iy, ix], order=1, mode="constant")


def e2_diameter(grid, I):
    """Equivalent 1/e^2 diameter from area above peak/e^2, plus x/y widths through the peak."""
    I = I / I.max()
    d = grid[1] - grid[0]
    area = (I >= np.exp(-2)).sum() * d * d
    deq = 2 * np.sqrt(area / np.pi)
    iy, ix = np.unravel_index(np.argmax(I), I.shape)
    wx = (I[iy, :] >= np.exp(-2)).sum() * d
    wy = (I[:, ix] >= np.exp(-2)).sum() * d
    return deq, wx, wy


# ------------------------------------------------------------------ main report
def report(elements, x, outdir, title="design"):
    os.makedirs(outdir, exist_ok=True)
    system = make(elements, x)
    D = SPEC["beam_diameter_mm"]
    wl0 = SPEC["ref_wl_um"]
    efl, bfl, zf = paraxial(system, wl0)
    ths = field_angles(efl)
    th_list = np.linspace(0, ths[-1], 6)
    common = np.linspace(-6, 6, 241)       # um, common PSF grid
    res = dict(elements=elements, x=list(map(float, x)), efl=efl, rows=[])

    # ---------------- PSF per wavelength / field + broadband PSF
    psf_store = {}
    for th in th_list:
        Pc, Dc = chief_ray(system, th, wl0)
        rms0, cen0, _ = wavefront(system, th, wl0, pupil_grid(12, D / 2, 1.5))
        ref = np.array([cen0[0], cen0[1], system.z_image])
        acc = np.zeros((len(common), len(common)))
        for wl in WLS_FINE:
            # each wavelength referenced to the SAME point -> lateral colour appears as a shift
            co, I = psf(system, th, wl, efl, ref)
            Ir = resample(co, I, common)
            acc += Ir / Ir.sum()
            rms, cen, _ = wavefront(system, th, wl, pupil_grid(12, D / 2, 1.5))
            deq, wx, wy = e2_diameter(common, Ir)
            S = np.exp(-(2 * np.pi * rms) ** 2)
            res["rows"].append(dict(field_deg=float(th), y_mm=float(cen0[1]), wl_nm=float(wl * 1e3),
                                    rms_waves=float(rms), strehl=float(S), d_e2_um=float(deq),
                                    wx_um=float(wx), wy_um=float(wy), d_diff_um=float(d_diff_um(wl, efl, D)),
                                    lat_color_um=float((cen[1] - cen0[1]) * 1e3)))
            psf_store[(round(th, 4), round(wl, 4))] = Ir
        deq, wx, wy = e2_diameter(common, acc)
        res["rows"].append(dict(field_deg=float(th), y_mm=float(cen0[1]), wl_nm=-1, d_e2_um=float(deq),
                                wx_um=float(wx), wy_um=float(wy)))
        psf_store[(round(th, 4), "bb")] = acc

    # ---------------- chromatic focal shift (paraxial focus vs wavelength) & depth of focus
    zf_l = np.array([paraxial(system, wl)[2] for wl in WLS_FINE])
    focal_shift_um = (zf_l - zf) * 1e3
    w0 = 0.5 * d_diff_um(wl0, efl, D)
    zR = np.pi * w0 ** 2 / wl0      # um
    # ---------------- telecentricity & distortion
    cra, ychief, opl_c = [], [], []
    for th in th_list:
        P, Dd = chief_ray(system, th, wl0)
        cra.append(np.degrees(np.arctan2(Dd[1], Dd[2])))
        ychief.append(P[1])
    ychief = np.array(ychief); cra = np.array(cra)
    ftheta = efl * np.radians(th_list)
    dist_ftheta = np.where(th_list > 0, 100 * (ychief - ftheta) / np.where(ftheta == 0, 1, ftheta), 0)
    res.update(bfl=bfl, working_distance=float(system.z_image - system.surfaces[-1].z),
               pupil_distance=float(system.surfaces[0].z),
               total_track=float(system.z_image),
               focal_shift_um=dict(zip([f"{w*1e3:.0f}" for w in WLS_FINE], focal_shift_um.round(2).tolist())),
               rayleigh_range_um=float(zR),
               chief_ray_angle_deg=cra.round(3).tolist(), field_deg=th_list.round(3).tolist(),
               y_chief_mm=ychief.round(4).tolist(), ftheta_distortion_pct=dist_ftheta.round(3).tolist(),
               edge_gap_mm=float(edge_gap(elements, x[1:len(elements)])) if len(elements) > 1 else None,
               max_scan_angle_deg=float(th_list[-1]))

    # ================================================================ plots
    # 1) layout
    fig, ax = plt.subplots(figsize=(10, 3.6))
    cols = plt.cm.viridis(np.linspace(0, 0.9, 3))
    for th, col in zip([0, ths[1], ths[2]], cols):
        ys = np.linspace(-D / 2, D / 2, 7)
        P0, D0 = launch(th, np.zeros_like(ys), ys, z0=-8)
        P, Dd, _, ok, hits = trace(system, P0, D0, wl0, check_aperture=False)
        pts = [P0] + hits
        for i in range(len(ys)):
            ax.plot([p[i, 2] for p in pts], [p[i, 1] for p in pts], color=col, lw=0.6)
    # draw lens outlines
    s = system.surfaces
    idx = 0
    for (name, rev), ns in zip(elements, system.n_surf):
        group = s[idx:idx + ns]; idx += ns
        for sf in group:
            r = np.linspace(-sf.semi_ap, sf.semi_ap, 60)
            z = sf.z + (0 if np.isinf(sf.R) else (r ** 2 / sf.R) / (1 + np.sqrt(1 - (r / sf.R) ** 2)))
            ax.plot(z, r, "k", lw=1)
        ax.text(group[0].z, group[0].semi_ap + 1.2, name + (" (rev)" if rev else ""), fontsize=8)
    ax.axvline(0, color="r", ls="--", lw=0.8); ax.text(0.3, -D, "scan pivot\n(galvo)", color="r", fontsize=8)
    ax.axvline(system.z_image, color="gray", lw=0.8)
    ax.set_aspect("equal"); ax.set_xlabel("z [mm]"); ax.set_ylabel("y [mm]")
    ax.set_title(f"{title}: EFL {efl:.2f} mm, pupil dist {s[0].z:.1f} mm, WD {res['working_distance']:.2f} mm")
    fig.tight_layout(); fig.savefig(f"{outdir}/layout.png", dpi=150); plt.close(fig)

    # 2) PSF 1/e^2 diameter vs field
    rows = res["rows"]
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    for wl in WLS_FINE[::2]:
        rr = [r for r in rows if abs(r["wl_nm"] - wl * 1e3) < 1e-6]
        axs[0].plot([r["y_mm"] for r in rr], [r["d_e2_um"] for r in rr], "o-",
                    color=plt.cm.turbo((wl - 0.48) / 0.3), label=f"{wl*1e3:.0f} nm")
        axs[0].axhline(d_diff_um(wl, efl, D), color=plt.cm.turbo((wl - 0.48) / 0.3), ls=":", lw=0.8)
    rb = [r for r in rows if r["wl_nm"] < 0]
    axs[0].plot([r["y_mm"] for r in rb], [r["d_e2_um"] for r in rb], "ks-", lw=2, label="500-750 nm (integrated)")
    axs[0].axhline(SPEC["spot_target_um"], color="r", ls="--", label="3 um target")
    axs[0].set_xlabel("field position at sample [mm]"); axs[0].set_ylabel("PSF 1/e² diameter [um]")
    axs[0].set_title("Spot size (FFT PSF, Gaussian 4 mm beam); dotted = diffraction limit")
    axs[0].legend(fontsize=7); axs[0].grid(alpha=0.3)
    for wl in WLS_FINE[::2]:
        rr = [r for r in rows if abs(r["wl_nm"] - wl * 1e3) < 1e-6]
        axs[1].plot([r["y_mm"] for r in rr], [r["strehl"] for r in rr], "o-",
                    color=plt.cm.turbo((wl - 0.48) / 0.3), label=f"{wl*1e3:.0f} nm")
    axs[1].axhline(0.8, color="r", ls="--", label="diffraction limited (0.8)")
    axs[1].set_ylim(0, 1.02); axs[1].set_xlabel("field position at sample [mm]"); axs[1].set_ylabel("Strehl (Gaussian weighted)")
    axs[1].legend(fontsize=7); axs[1].grid(alpha=0.3); axs[1].set_title("Strehl ratio, common focal plane")
    fig.tight_layout(); fig.savefig(f"{outdir}/spot_size_vs_field.png", dpi=150); plt.close(fig)

    # 3) PSF images (broadband) at 0, 0.6, 1.0 field
    fig, axs = plt.subplots(2, 3, figsize=(11, 7))
    for j, th in enumerate([th_list[0], th_list[3], th_list[-1]]):
        I = psf_store[(round(th, 4), "bb")]
        axs[0, j].imshow(I / I.max(), extent=[common[0], common[-1], common[0], common[-1]], origin="lower", cmap="inferno")
        axs[0, j].contour(common, common, I / I.max(), [np.exp(-2)], colors="c", linewidths=0.8)
        axs[0, j].set_title(f"broadband PSF, field {th:.2f}°"); axs[0, j].set_xlim(-4, 4); axs[0, j].set_ylim(-4, 4)
        axs[0, j].set_xlabel("x [um]"); axs[0, j].set_ylabel("y [um]")
        # geometric spot diagram
        for wl in [0.5, 0.5625, 0.625, 0.6875, 0.75]:
            g = pupil_grid(10, D / 2, 1.0)
            P0, D0 = launch(th, g[0], g[1])
            P, *_ = trace(system, P0, D0, wl)
            Pc, _ = chief_ray(system, th, wl0)
            axs[1, j].plot((P[:, 0] - Pc[0]) * 1e3, (P[:, 1] - Pc[1]) * 1e3, ".", ms=2,
                           color=plt.cm.turbo((wl - 0.48) / 0.3), label=f"{wl*1e3:.0f}")
        circ = plt.Circle((0, 0), SPEC["spot_target_um"] / 2, fill=False, color="k", ls="--")
        axs[1, j].add_patch(circ); axs[1, j].set_aspect("equal")
        axs[1, j].set_title(f"geometric spot (1/e² pupil), dashed = Ø3 um"); axs[1, j].set_xlabel("x [um]")
        lim = max(3, 1.2 * np.nanmax(np.abs((P[:, 1] - Pc[1]) * 1e3)))
        axs[1, j].set_xlim(-lim, lim); axs[1, j].set_ylim(-lim, lim)
    axs[1, 0].legend(fontsize=7, markerscale=4)
    fig.tight_layout(); fig.savefig(f"{outdir}/psf_and_spots.png", dpi=150); plt.close(fig)

    # 4) chromatic focal shift + lateral colour + telecentricity/distortion
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
    axs[0].plot(WLS_FINE * 1e3, focal_shift_um, "o-")
    axs[0].axhspan(-zR, zR, color="g", alpha=0.15, label=f"± Rayleigh range ({zR:.1f} um)")
    axs[0].set_xlabel("wavelength [nm]"); axs[0].set_ylabel("focal shift [um]"); axs[0].legend(fontsize=8)
    axs[0].set_title("Axial colour (paraxial)"); axs[0].grid(alpha=0.3)
    for wl in WLS_FINE[::2]:
        rr = [r for r in rows if abs(r["wl_nm"] - wl * 1e3) < 1e-6]
        axs[1].plot([r["y_mm"] for r in rr], [r["lat_color_um"] for r in rr], "o-",
                    color=plt.cm.turbo((wl - 0.48) / 0.3), label=f"{wl*1e3:.0f} nm")
    axs[1].set_xlabel("field [mm]"); axs[1].set_ylabel("centroid shift vs 625 nm [um]")
    axs[1].set_title("Lateral colour"); axs[1].legend(fontsize=7); axs[1].grid(alpha=0.3)
    axs[2].plot(ychief, cra, "o-", label="chief-ray angle [deg]")
    axs[2].plot(ychief, dist_ftheta, "s-", label="f-θ distortion [%]")
    axs[2].set_xlabel("field [mm]"); axs[2].legend(fontsize=8); axs[2].grid(alpha=0.3)
    axs[2].set_title("Telecentricity & f-θ linearity")
    fig.tight_layout(); fig.savefig(f"{outdir}/color_telecentricity.png", dpi=150); plt.close(fig)

    # ---------------- prescription table
    lines = ["| # | Surface | Radius [mm] | Thickness to next [mm] | Material | Clear semi-ap [mm] |",
             "|---|---|---|---|---|---|",
             f"| 0 | Scan pivot / entrance pupil (Ø{D:.0f} mm 1/e² beam) | ∞ | {s[0].z:.3f} | AIR | {D/2:.2f} |"]
    for k, sf in enumerate(s):
        nxt = (s[k + 1].z if k + 1 < len(s) else system.z_image) - sf.z
        R = "∞" if np.isinf(sf.R) else f"{sf.R:.3f}"
        lines.append(f"| {k+1} | {sf.label} | {R} | {nxt:.3f} | {sf.mat} | {sf.semi_ap:.2f} |")
    lines.append(f"| {len(s)+1} | Image (sample focal plane) | ∞ | – | – | – |")
    with open(f"{outdir}/prescription.md", "w") as f:
        f.write("\n".join(lines) + "\n")
    with open(f"{outdir}/performance.json", "w") as f:
        json.dump(res, f, indent=1)
    return res


def summary_text(res):
    rows = res["rows"]
    bb = [r for r in rows if r["wl_nm"] < 0]
    mono = [r for r in rows if r["wl_nm"] > 0]
    out = [f"EFL {res['efl']:.3f} mm | pupil distance {res['pupil_distance']:.2f} mm | WD {res['working_distance']:.2f} mm"
           f" | total track {res['total_track']:.1f} mm | max scan angle ±{res['max_scan_angle_deg']:.2f}°",
           f"broadband 1/e² spot diameter (on-axis -> edge): " + ", ".join(f"{r['d_e2_um']:.2f}" for r in bb) + " um",
           f"worst monochromatic spot {max(r['d_e2_um'] for r in mono):.2f} um, min Strehl {min(r['strehl'] for r in mono):.3f}",
           f"focal shift 500-750 nm: {min(res['focal_shift_um'].values()):.1f} .. {max(res['focal_shift_um'].values()):.1f} um"
           f" (Rayleigh range {res['rayleigh_range_um']:.1f} um)",
           f"max lateral colour {max(abs(r['lat_color_um']) for r in mono):.2f} um | chief-ray angle at edge"
           f" {res['chief_ray_angle_deg'][-1]:.2f}° | f-θ distortion at edge {res['ftheta_distortion_pct'][-1]:.3f}%"]
    return "\n".join(out)


if __name__ == "__main__":
    data = json.load(open(sys.argv[1]))
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    r = data["results"][k]
    el = [tuple(e) for e in r["elements"]]
    res = report(el, np.array(r["x"]), sys.argv[3] if len(sys.argv) > 3 else "results/best")
    print(summary_text(res))
