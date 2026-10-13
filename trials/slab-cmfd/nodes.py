import re
import sys
from pathlib import Path


def read(p):
    k = None
    nodes = {}
    for line in Path(p).read_text().splitlines():
        m = re.search(r"node node=(\d+) average_flux=\[(\S+), (\S+)\]", line)
        if m:
            nodes[int(m[1])] = (float(m[2]), float(m[3]))
        m = re.search(r"xlc_lattice::moc: iteration=\d+ k=(\S+)", line)
        if m:
            k = float(m[1])
    return k, nodes


for a, b in zip(sys.argv[1::2], sys.argv[2::2]):
    ka, na = read(a)
    kb, nb = read(b)
    d = max(abs(nb[n][g] / na[n][g] - 1) for n in na for g in range(2))
    print(
        f"{a} vs {b}: k {ka!r} {kb!r} dk/k={(kb - ka) / ka * 1e5:.2e} pcm, nodes={len(na)} max rel={d:.2e}"
    )
