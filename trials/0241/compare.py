"""Compare lattice runs of two binaries (#0241): k, MOC and coarse-mesh iterations, coarse span time, and the largest
relative difference of every float in the msgpack outputs.

Usage: python3 -I compare.py <case> <base> <other> <output>...   (outputs are file names with {} for the binary name)
"""

import re
import sys
from pathlib import Path

import msgpack


def floats(x, path=""):
    if isinstance(x, float):
        yield path, x
    elif isinstance(x, dict):
        for k, v in x.items():
            yield from floats(v, f"{path}/{k}")
    elif isinstance(x, (list, tuple)):
        for i, v in enumerate(x):
            yield from floats(v, f"{path}[{i}]")


def log_summary(path):
    text = Path(path).read_text()
    k = re.findall(r"xlc_lattice::moc: .*?iteration=(\d+) k=([0-9.e+-]+)", text)
    states = re.findall(r"lattice state=\d+ iterations=(\d+)", text)
    coarse = [
        int(n) for n in re.findall(r"coarse-mesh eigenvalue iterations=(\d+)", text)
    ]
    span = 0.0
    for value, unit in re.findall(
        r"\bcoarse: [^\n]*close time\.busy=([0-9.]+)(ns|µs|ms|s)\b", text
    ):
        span += float(value) * {"ns": 1e-9, "µs": 1e-6, "ms": 1e-3, "s": 1.0}[unit]
    seconds = re.findall(r"lattice state=\d+ .*? seconds=([0-9.]+)", text)
    return {
        "last k": k[-1][1] if k else None,
        "MOC iterations": states,
        "coarse solves": len(coarse),
        "coarse iterations": sum(coarse),
        "coarse span s": round(span, 3),
        "MOC seconds": seconds,
    }


def main(case, base, other, outputs):
    for name in (base, other):
        print(name, log_summary(f"{case}-{name}.log"))
    for pattern in outputs:
        a = msgpack.unpackb(
            Path(pattern.format(base)).read_bytes(), raw=False, strict_map_key=False
        )
        b = msgpack.unpackb(
            Path(pattern.format(other)).read_bytes(), raw=False, strict_map_key=False
        )
        fa, fb = list(floats(a)), list(floats(b))
        assert [p for p, _ in fa] == [p for p, _ in fb], pattern
        worst, where, same = 0.0, None, 0
        for (p, x), (_, y) in zip(fa, fb):
            if x == y:
                same += 1
                continue
            d = abs(x - y) / max(abs(x), abs(y))
            if d > worst:
                worst, where = d, p
        print(
            f"  {pattern}: {len(fa)} floats, {same} identical, max relative difference {worst:.3e} at {where}"
        )


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
