"""Condition numbers of core coarse-mesh blocks built from lattice few-group constants (#0241, ADR-0275).

The block of crates/core/src/nem/cmfd.rs `assemble` is V (diag(Σt − Σgg) − off-diagonal Σs(from→g)) plus, on the diagonal,
the face terms area (ĵ_out + D̃ S/2)/φ, which are positive for inner faces. "no leakage" drops them (the smallest diagonal);
"leakage" adds 2 D/h per face of a 21.5 × 21.5 × h_z node (the finite-difference D̃ with D̂ = 0).
"""

import glob
import sys
from pathlib import Path

import msgpack
import numpy as np


def block(state, leak, hz):
    t = np.array(state["total"])
    s = np.array(state["scatter"])  # [from][to]
    d = np.array(state["diffusion"])
    ng = len(t)
    hx = 21.5
    volume = hx * hx * hz
    a = np.zeros((ng, ng))
    for g in range(ng):
        a[g, g] = volume * (t[g] - s[g][g])
        for f in range(ng):
            if f != g:
                a[g, f] = -volume * s[f][g]
    if leak:
        for g in range(ng):
            for h, area in [(hx, hx * hz)] * 4 + [(hz, hx * hx)] * 2:
                a[g, g] += area * 2 * d[g] / h
    return a


def main(paths):
    worst = {}
    for path in paths:
        data = msgpack.unpackb(Path(path).read_bytes(), raw=False, strict_map_key=False)
        for state in data["states"]:
            for leak in (False, True):
                for hz in (4.0, 8.0, 50.0):
                    c = np.linalg.cond(block(state, leak, hz))
                    key = (leak, hz)
                    if c > worst.get(key, (0, ""))[0]:
                        worst[key] = (c, path)
    for (leak, hz), (c, path) in sorted(worst.items()):
        print(f"leakage={leak} hz={hz}: max cond {c:.3e} ({path})")


if __name__ == "__main__":
    main(sorted(path for p in sys.argv[1:] for path in glob.glob(p)))
