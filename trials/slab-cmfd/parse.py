import re
import sys
from pathlib import Path


def sec(v):
    m = re.match(r"([\d.]+)(ns|µs|ms|s)$", v)
    x, u = float(m[1]), m[2]
    return x * {"ns": 1e-9, "µs": 1e-6, "ms": 1e-3, "s": 1.0}[u]


for path in sys.argv[1:]:
    t = {}
    its = []
    moc = 0
    k = None
    fallback = 0
    total = None
    first = last = None
    for line in Path(path).read_text().splitlines():
        ts = re.match(r"(\d{4}-\S+Z)", line)
        m = re.search(r":(\w+): xlc_lattice::moc: close time.busy=(\S+)", line)
        if m:
            t[m[1]] = t.get(m[1], 0) + sec(m[2])
        m = re.search(r"coarse-mesh eigenvalue iterations=(\d+)", line)
        if m:
            its.append(int(m[1]))
        if "solving unshifted" in line:
            fallback += 1
        m = re.search(r"xlc_lattice::moc: iteration=(\d+) k=(\S+)", line)
        if m:
            moc = int(m[1])
            k = m[2]
        m = re.search(r"\bmoc\{[^}]*\}: xlc_lattice::moc: close time.busy=(\S+)", line)
        if m:
            total = sec(m[1])
    print(
        f"{path}: moc_iter={moc} k={k} moc_busy={total} "
        + " ".join(f"{a}={b:.2f}" for a, b in t.items())
        + f" coarse_iter_sum={sum(its)} max={max(its) if its else 0} solves={len(its)} fallback={fallback}"
    )
