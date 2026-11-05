"""Off-the-shelf cemented achromatic doublet prescriptions (Thorlabs AC-series, VIS 'A' coating).

Values: R1, R2, R3 [mm], centre thicknesses tc1, tc2 [mm], glasses (crown, flint), clear-aperture
diameter [mm], nominal EFL [mm]. The orientation stored here is the catalogue orientation
(R1 = the crown surface, which faces the collimated beam for infinite conjugate use).

NOTE: these prescriptions were entered from the vendor data sheets from memory and are
sanity-checked by re-computing their EFL (see `python -m scanlens.check_catalog`).
ALWAYS confirm against the vendor's Zemax file before ordering / building.
"""

# Only entries whose recomputed EFL matches the catalogue EFL to < 1 % are kept
# (AC254-040/050/060/075 and AC127-025 failed that check and were removed).
LENSES = {
    #  name           R1      R2       R3      tc1   tc2   crown      flint     dia   EFL
    "AC254-030-A": (20.89, -16.73, -79.80, 12.0, 2.0, "N-BAF10", "N-SF6HT", 25.4, 30.0),
    "AC254-035-A": (24.0, -19.1, -102.1, 10.0, 2.5, "N-BAF10", "N-SF6HT", 25.4, 35.0),
    "AC254-045-A": (31.2, -25.9, -130.6, 7.0, 2.0, "N-BAF10", "N-SF6HT", 25.4, 45.0),
    "AC254-100-A": (62.8, -45.7, -128.2, 4.0, 2.5, "N-BK7", "N-SF5", 25.4, 100.0),
    "AC127-019-A": (12.9, -11.0, -59.3, 4.5, 1.5, "N-BAF10", "N-SF6HT", 12.7, 19.0),
    "AC127-030-A": (17.9, -13.5, -44.2, 3.5, 1.5, "N-BK7", "N-SF2", 12.7, 30.0),
    "AC127-050-A": (27.4, -22.5, -91.8, 3.5, 1.5, "N-BK7", "N-SF2", 12.7, 50.0),
}


def doublet_surfaces(name, reversed_=False):
    """Return list of (radius, thickness_after, material_after) for the glass part of a doublet,
    plus clear diameter. The thickness after the last surface is left as None (air gap set later)."""
    R1, R2, R3, t1, t2, g1, g2, dia, _ = LENSES[name]
    if not reversed_:
        surfs = [(R1, t1, g1), (R2, t2, g2), (R3, None, "AIR")]
    else:
        surfs = [(-R3, t2, g2), (-R2, t1, g1), (-R1, None, "AIR")]
    return surfs, dia


# ---------------------------------------------------------------------------------------------
# Generic N-BK7 singlets (Thorlabs LA/LC-series, Edmund equivalents). For a plano lens the radius
# follows from the catalogue focal length: R = (n_d - 1) * f with n_d(N-BK7) = 1.5168 (this is how
# the vendors specify them, e.g. Thorlabs LC1715 f = -50 mm -> R = 25.8 mm).
# Centre thickness: plano-concave Ø1" = 3.5 mm, Ø1/2" = 2.0 mm (vendor standard);
# plano-convex: sag + 1.8 mm edge thickness (close to vendor values; verify per part number).
# name -> (f_nominal, diameter, kind)
SINGLETS = {
    "LC-BK7-f-25 (1in)": (-25.0, 25.4, "pcv"),
    "LC-BK7-f-30 (1in)": (-30.0, 25.4, "pcv"),
    "LC-BK7-f-40 (1in)": (-40.0, 25.4, "pcv"),
    "LC-BK7-f-50 (1in)": (-50.0, 25.4, "pcv"),
    "LC-BK7-f-75 (1in)": (-75.0, 25.4, "pcv"),
    "LC-BK7-f-100 (1in)": (-100.0, 25.4, "pcv"),
    "LC-BK7-f-150 (1in)": (-150.0, 25.4, "pcv"),
    "LC-BK7-f-200 (1in)": (-200.0, 25.4, "pcv"),
}
N_D_BK7 = 1.5168


def singlet_params(name):
    import numpy as np
    f, dia, kind = SINGLETS[name]
    R = (N_D_BK7 - 1) * abs(f)
    if kind == "pcv":     # plano-concave: concave side first (R<0 facing +z light means -R?)
        tc = 3.5 if dia > 20 else 2.0
    else:
        r = 0.45 * dia
        tc = R - np.sqrt(R * R - r * r) + 1.8
    return R, tc, dia


def element_surfaces(name, reversed_=False):
    """Generic accessor for doublets and singlets.
    Singlet default orientation: curved surface first (towards the scan pivot)."""
    if name in LENSES:
        return doublet_surfaces(name, reversed_)
    import numpy as np
    f, dia, kind = SINGLETS[name]
    R, tc, dia = singlet_params(name)
    Rc = -R if f < 0 else R      # concave-first: centre of curvature towards -z
    if not reversed_:
        surfs = [(Rc, tc, "N-BK7"), (np.inf, None, "AIR")]
    else:
        surfs = [(np.inf, tc, "N-BK7"), (-Rc, None, "AIR")]
    return surfs, dia


ALL_ELEMENTS = list(LENSES) + list(SINGLETS)
