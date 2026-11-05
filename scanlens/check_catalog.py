"""Sanity check: recompute EFL / BFL of each catalogue doublet at 587.6 nm."""
from .catalog import LENSES
from .raytrace import build_system, paraxial

if __name__ == "__main__":
    print(f"{'lens':14s} {'nominal':>8s} {'EFL calc':>9s} {'BFL calc':>9s}  err%")
    for name, v in LENSES.items():
        sys_ = build_system([(name, False)], 10.0, [])
        efl, bfl, _ = paraxial(sys_, 0.5876)
        print(f"{name:14s} {v[-1]:8.1f} {efl:9.2f} {bfl:9.2f}  {100*(efl-v[-1])/v[-1]:+.1f}")
