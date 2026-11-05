"""Re-optimise the best pre-screened candidates with Optiland (true broadband PSF merit) and rank.

usage: python -m scanlens.refine N_TOP N_PROCS results/results_search_2elem.json results/hybrid_shard*.json
"""
import glob, json, sys
import numpy as np
from multiprocessing import Pool


def _job(r):
    import warnings; warnings.filterwarnings("ignore")
    from scanlens.optiland_optimize import optimise_optiland
    from scanlens.optiland_analysis import oct_psf_diameter
    from scanlens.design import SPEC
    el = [tuple(e) for e in r["elements"]]
    th = float(np.degrees(np.arctan(SPEC["field_half_mm"] / r["efl"])))
    try:
        op, x, _ = optimise_optiland(el, r["x"], th, maxiter=40)
        d = [oct_psf_diameter(op, h) for h in (0.0, 0.7, 1.0)]
        efl = float(np.asarray(op.paraxial.f2()))
    except Exception as e:
        return dict(elements=r["elements"], error=str(e))
    out = dict(elements=r["elements"], x=x, efl=efl, field_max_deg=th, psf_d_um=d,
               worst_um=max(d), prescreen_merit=r["merit"])
    print(f"{max(d):6.3f} um  EFL {efl:6.2f}  {el}", flush=True)
    return out


if __name__ == "__main__":
    n_top, n_proc = int(sys.argv[1]), int(sys.argv[2])
    cands = []
    for pat in sys.argv[3:]:
        for f in glob.glob(pat):
            cands += json.load(open(f))["results"]
    cands.sort(key=lambda r: r["merit"])
    seen, top = set(), []
    for r in cands:
        key = json.dumps(r["elements"])
        if key not in seen:
            seen.add(key); top.append(r)
        if len(top) >= n_top:
            break
    with Pool(n_proc) as p:
        res = p.map(_job, top, chunksize=1)
    res = [r for r in res if "worst_um" in r]
    res.sort(key=lambda r: r["worst_um"])
    json.dump(dict(engine="optiland", results=res), open("results/optiland_refined.json", "w"), indent=1)
    for r in res[:10]:
        print(f"{r['worst_um']:.3f} {[round(v,2) for v in r['psf_d_um']]} EFL {r['efl']:.2f} {r['elements']}")
