"""2 つの HDF5 の、すべてのデータセットを比べ、値の違うものを挙げる（ADR-0124 の確認）。"""

import sys

import h5py
import numpy as np


def datasets(path: str) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}

    def visit(name: str, obj: object) -> None:
        if isinstance(obj, h5py.Dataset):
            out[name] = obj[()]

    with h5py.File(path, "r") as f:
        f.visititems(visit)
    return out


a, b = datasets(sys.argv[1]), datasets(sys.argv[2])
differ = [
    k
    for k in sorted(a.keys() | b.keys())
    if k not in a or k not in b or not np.array_equal(a[k], b[k])
]
print(f"{len(a)} datasets; differ: {differ}")
