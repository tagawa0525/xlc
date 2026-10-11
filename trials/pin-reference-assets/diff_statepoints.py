"""2 つの statepoint の k と tally の、最大の相対差を書く（ADR-0124 の確認）。"""

import sys

import h5py
import numpy as np

with h5py.File(sys.argv[1], "r") as a, h5py.File(sys.argv[2], "r") as b:
    print(f"runtime/total {a['runtime/total'][()]:.1f} s")
    for name in [
        "k_generation",
        "k_combined",
        "global_tallies",
        "tallies/tally 1/results",
    ]:
        if name not in a:
            continue
        x, y = a[name][()], b[name][()]
        scale = np.maximum(np.abs(x), np.abs(y))
        rel = np.abs(x - y)[scale > 0] / scale[scale > 0]
        print(f"{name}: max relative difference {rel.max():.2e}")
