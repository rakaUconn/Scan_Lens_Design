"""Refractive index lookup (Schott / Hikari catalogs via optiland's refractiveindex.info data)."""
from functools import lru_cache

_ALIASES = {"E-BAF11": ("E-BAF11", "hikari")}


@lru_cache(maxsize=None)
def _material(name):
    from optiland.materials import Material
    if name in _ALIASES:
        n, ref = _ALIASES[name]
        return Material(n, ref)
    return Material(name, "schott")


@lru_cache(maxsize=None)
def index(name, wl_um):
    """Refractive index of `name` at wavelength `wl_um` (micrometres). 'AIR' -> 1."""
    if name.upper() == "AIR":
        return 1.0
    return float(_material(name).n(wl_um)[0])
