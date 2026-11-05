"""Search over combinations of catalogue achromats + optimise air spaces for a visible OCT scan lens.

Figure of merit = worst-case *effective* 1/e^2 spot diameter over wavelengths 500-750 nm and the
scan field, estimated from the Gaussian diffraction limit, the RMS wavefront error (via Strehl)
and the lateral colour spread. Final candidates are re-evaluated with a full FFT PSF in evaluate.py.
"""
import itertools, json, sys, time
import numpy as np
from scipy.optimize import minimize

from .catalog import LENSES, SINGLETS, element_surfaces
from .raytrace import build_system, paraxial, pupil_grid, wavefront, chief_ray

# ------------------------------------------------------------------ specification
SPEC = dict(
    wavelengths_um=[0.500, 0.5625, 0.625, 0.6875, 0.750],
    ref_wl_um=0.625,
    beam_diameter_mm=4.0,          # 1/e^2 diameter at the scan pivot (entrance pupil)
    spot_target_um=3.0,            # 1/e^2 spot diameter target
    field_half_mm=0.71,            # FOV 1 x 1 mm square -> corner radius 0.707 mm
    min_pupil_dist_mm=10.0,        # galvo pivot -> first lens vertex (mechanical clearance)
    min_working_dist_mm=3.0,       # last lens vertex -> focal plane
    min_edge_gap_mm=0.3,           # lens-to-lens edge clearance
    min_vertex_gap_mm=0.3,
)
FIELD_FRACS = [0.0, 0.7, 1.0]


def d_diff_um(wl_um, efl_mm, D_mm):
    """Diffraction-limited 1/e^2 Gaussian spot diameter [um]."""
    return 4 * wl_um * efl_mm / (np.pi * D_mm)


def sag(R, r):
    if np.isinf(R):
        return 0.0
    c = 1 / R
    return c * r * r / (1 + np.sqrt(max(1 - c * c * r * r, 0)))


def edge_gap(elements, gaps):
    """Minimum axial clearance between consecutive elements at the smaller clear aperture."""
    out = []
    for k in range(len(elements) - 1):
        s1, d1 = element_surfaces(*elements[k])
        s2, d2 = element_surfaces(*elements[k + 1])
        r = 0.5 * min(d1, d2) * 0.9
        out.append(gaps[k] - sag(s1[-1][0], r) + sag(s2[0][0], r))
    return min(out) if out else np.inf


def field_angles(efl):
    th = np.degrees(np.arctan(SPEC["field_half_mm"] / efl))
    return [f * th for f in FIELD_FRACS]


def evaluate_fast(system, grid, efl=None):
    """Return (effective worst-case spot diameter [um], detail dict)."""
    D = SPEC["beam_diameter_mm"]
    if efl is None:
        efl = paraxial(system, SPEC["ref_wl_um"])[0]
    if not (efl > 0):
        return 1e3, {}
    worst, details = 0.0, []
    for th in field_angles(efl):
        cens, dmax = [], 0.0
        for wl in SPEC["wavelengths_um"]:
            rms, cen, _ = wavefront(system, th, wl, grid)
            if not np.isfinite(rms):
                return 1e3, {}
            S = np.exp(-(2 * np.pi * rms) ** 2)
            d = d_diff_um(wl, efl, D) / np.sqrt(max(S, 1e-3))
            dmax = max(dmax, d)
            cens.append(cen[1] * 1e3)  # um
            details.append((th, wl, rms, d))
        lat = max(cens) - min(cens)                  # lateral colour spread [um]
        d_eff = np.sqrt(dmax ** 2 + lat ** 2)         # spectrally integrated spot (approx.)
        worst = max(worst, d_eff)
    return worst, dict(efl=efl, details=details)


def make(elements, x):
    """x = [d_pupil, gaps..., working distance]"""
    n = len(elements)
    return build_system(elements, x[0], list(x[1:n]), x[n], pupil_semi=SPEC["beam_diameter_mm"] / 2)


def merit(x, elements, grid):
    n = len(elements)
    pen = 0.0
    pen += 50 * max(0, SPEC["min_pupil_dist_mm"] - x[0]) ** 2
    pen += 50 * max(0, SPEC["min_working_dist_mm"] - x[n]) ** 2
    for g in x[1:n]:
        pen += 50 * max(0, SPEC["min_vertex_gap_mm"] - g) ** 2
    if n > 1:
        pen += 50 * max(0, SPEC["min_edge_gap_mm"] - edge_gap(elements, x[1:n])) ** 2
    try:
        sys_ = make(elements, x)
        val, _ = evaluate_fast(sys_, grid)
    except Exception:
        val = 1e3
    if not np.isfinite(val):
        val = 1e3
    return val + pen


def initial_guess(elements):
    """Paraxial focus as starting working distance; small gaps."""
    n = len(elements)
    x = [SPEC["min_pupil_dist_mm"] + 2] + [1.0] * (n - 1) + [10.0]
    sys_ = make(elements, x)
    efl, bfl, _ = paraxial(sys_, SPEC["ref_wl_um"])
    x[n] = max(bfl, 0.5)
    return np.array(x), efl


def optimise(elements, grid, x0=None):
    if x0 is None:
        x0, _ = initial_guess(elements)
    best = None
    for start in [x0, x0 + np.r_[5, [0.5] * (len(elements) - 1), 0]]:
        r = minimize(merit, start, args=(elements, grid), method="Nelder-Mead",
                     options=dict(maxiter=600, xatol=1e-3, fatol=1e-4))
        if best is None or r.fun < best.fun:
            best = r
    return best


def configurations(n_elem, efl_window=(8.0, 20.0)):
    names = list(LENSES)
    for combo in itertools.product(names, repeat=n_elem):
        for revs in itertools.product([False, True], repeat=n_elem):
            elements = list(zip(combo, revs))
            x0, efl = initial_guess(elements)
            if efl_window[0] <= efl <= efl_window[1]:
                yield elements


def run(n_elem=2, efl_window=(8.0, 20.0), out="results_search.json"):
    grid = pupil_grid(n_rings=5, semi=SPEC["beam_diameter_mm"] / 2, extent=1.25)
    results = []
    cfgs = list(configurations(n_elem, efl_window))
    print(f"{len(cfgs)} configurations with {n_elem} elements in EFL window {efl_window}", flush=True)
    t0 = time.time()
    for i, el in enumerate(cfgs):
        r = optimise(el, grid)
        sys_ = make(el, r.x)
        efl = paraxial(sys_, SPEC["ref_wl_um"])[0]
        results.append(dict(elements=el, x=list(map(float, r.x)), merit=float(r.fun), efl=float(efl)))
        print(f"[{i+1}/{len(cfgs)} {time.time()-t0:5.0f}s] {el}  merit={r.fun:.3f} um  EFL={efl:.2f}", flush=True)
    results.sort(key=lambda d: d["merit"])
    with open(out, "w") as f:
        json.dump(dict(spec=SPEC, results=results), f, indent=1)
    return results


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    run(n, out=f"results_search_{n}elem.json")


# ------------------------------------------------------------------ 2 doublets + negative singlet
def axial_color_range_um(system):
    zs = [paraxial(system, wl)[2] for wl in np.linspace(0.5, 0.75, 6)]
    return (max(zs) - min(zs)) * 1e3


def hybrid_configurations(efl_window=(10.0, 18.0), max_color_um=45.0):
    dbl = list(LENSES)
    sng = list(SINGLETS)
    seen = []
    for a, b in itertools.product(dbl, repeat=2):
        for ra, rb in itertools.product([False, True], repeat=2):
            for s in sng:
                for rs in [False, True]:
                    for pos in range(3):
                        el = [(a, ra), (b, rb)]
                        el.insert(pos, (s, rs))
                        x0, efl = initial_guess(el)
                        if not (efl_window[0] <= efl <= efl_window[1]):
                            continue
                        sys_ = make(el, x0)
                        if axial_color_range_um(sys_) > max_color_um:
                            continue
                        seen.append(el)
    return seen


def run_hybrid(out="results/results_search_hybrid.json", start=0, stop=None):
    grid = pupil_grid(n_rings=5, semi=SPEC["beam_diameter_mm"] / 2, extent=1.25)
    cfgs = hybrid_configurations()[start:stop]
    print(f"{len(cfgs)} hybrid configurations", flush=True)
    results, t0 = [], time.time()
    for i, el in enumerate(cfgs):
        r = optimise(el, grid)
        sys_ = make(el, r.x)
        efl = paraxial(sys_, SPEC["ref_wl_um"])[0]
        results.append(dict(elements=el, x=list(map(float, r.x)), merit=float(r.fun), efl=float(efl)))
        print(f"[{i+1}/{len(cfgs)} {time.time()-t0:5.0f}s] {el}  merit={r.fun:.3f} um  EFL={efl:.2f}", flush=True)
        if (i + 1) % 20 == 0 or i + 1 == len(cfgs):
            results.sort(key=lambda d: d["merit"])
            with open(out, "w") as f:
                json.dump(dict(spec=SPEC, results=results), f, indent=1)
    return results
