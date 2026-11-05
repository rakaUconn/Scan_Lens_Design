"""Longer Optiland optimisation of the top refined candidates.
usage: python -m scanlens.polish IN.json OUT.json N_TOP [--cornea]"""
import json, sys, warnings
import numpy as np
from multiprocessing import Pool


def _job(args):
    r, cornea = args
    warnings.filterwarnings("ignore")
    from scanlens.optiland_optimize import optimise_optiland
    from scanlens.optiland_analysis import oct_psf_diameter
    el = [tuple(e) for e in r["elements"]]
    op, x, _ = optimise_optiland(el, r["x"], r["field_max_deg"], maxiter=120, cornea=cornea)
    d = [oct_psf_diameter(op, h) for h in (0.0, 0.7, 1.0)]
    out = dict(r, x=x, psf_d_um=d, worst_um=max(d), efl=float(np.asarray(op.paraxial.f2())))
    print(f"{max(d):.3f} {[round(v, 2) for v in d]} {el}", flush=True)
    return out


if __name__ == "__main__":
    src, dst, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
    cornea = "--cornea" in sys.argv
    rs = json.load(open(src))["results"][:n]
    with Pool(min(n, 4)) as p:
        res = p.map(_job, [(r, cornea) for r in rs], chunksize=1)
    res.sort(key=lambda r: r["worst_um"])
    json.dump(dict(engine="optiland", cornea=cornea, results=res), open(dst, "w"), indent=1)
