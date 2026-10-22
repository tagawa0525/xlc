"""#0133 の試行: 2A の燃焼の履歴の燃焼点 <番号> の燃料の組成（対称操作の組ごと、輸送する核種）で、xlc の静的な格子の入力を書く。
組ごとに燃料の物質（fuel-<i>）とピンのセル（pin-<i>）を分け、ほかは lattice-2A.toml と同じ（[depletion] を除く）。

使い方: python3 write_2a.py <lattice-2A.toml> <lattice-2A.msgpack> <openmc-2A.msgpack> <番号> <xlc|openmc> <出力> [他方から取る核種 ...]

組成の基は <xlc|openmc> の履歴で、後ろに並べた核種（xlc の名前。例 U-235）だけを他方の履歴の値にする。
"""

import os
import re
import sys
from pathlib import Path

import msgpack

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compare_depletion_assembly import check_same_run, openmc_name

template, xpath, opath, point, base, out, *swap = sys.argv[1:]
point = int(point)
x = msgpack.unpackb(Path(xpath).read_bytes(), strict_map_key=False)
o = msgpack.unpackb(Path(opath).read_bytes(), strict_map_key=False)
# 形式、輸送する核種、燃焼点、領域の順（対称操作の組と体積）が同じこと
check_same_run(x, o)
if base not in ("xlc", "openmc"):
    raise SystemExit(f"base must be xlc or openmc, got {base!r}")
idx = {n: i for i, n in enumerate(x["nuclides"])}
N = 17


def orbit(i, j):
    c = N // 2
    return tuple(sorted((abs(i - c), abs(j - c)), reverse=True))


assert abs(x["points"][point]["burnup"] - o["burnup"][point]) < 1e-6
for name in swap:
    assert name in x["transportable"], name
orbits = [tuple(r["orbit"]) for r in o["regions"]]
text = Path(template).read_text()
text = text[: text.index("[depletion]")].rstrip() + "\n"
library = Path(template).resolve().parent / re.search(
    r'^library = "(.*)"$', text, re.MULTILINE
).group(1)
text = re.sub(
    r'^library = ".*"$',
    f'library = "{os.path.relpath(library, Path(out).resolve().parent)}"',
    text,
    count=1,
    flags=re.MULTILINE,
)
head, rest = text.split("[material.fuel]", 1)
rest = rest.split("\n", 3)[3]  # temperature と nuclides の行を除く
cell_start = rest.index("# 燃料のピン")
cell_end = rest.index("# 案内管")
pin_cell = rest[cell_start:cell_end]
materials = []
for r, orb in enumerate(orbits):
    xd = x["points"][point]["densities"][r]
    od = o["regions"][r]["densities"]
    dens = {}
    for n in x["transportable"]:
        vx = xd[idx[n]]
        vo = od[openmc_name(n)][point]
        v = vo if (base == "openmc") != (n in swap) else vx
        if v > 0.0:
            dens[n] = v
    comp = ", ".join(f'"{n}" = {v:.10e}' for n, v in dens.items())
    materials.append(
        f"[material.fuel-{r}]\ntemperature = 565\nnuclides = {{ {comp} }}\n"
    )
cells = "".join(
    pin_cell.replace("[cell.pin]", f"[cell.pin-{r}]").replace(
        'material = "fuel"', f'material = "fuel-{r}"'
    )
    for r in range(len(orbits))
)
rest = rest[:cell_start] + cells + rest[cell_end:]
# 行の pin を組ごとのセルに置き換える
lines = head.split("\n")
out_lines, row = [], 0
for line in lines:
    m = re.fullmatch(r'  "(.*)",', line)
    if m and row < N:
        tokens = m.group(1).split()
        new = [
            f"pin-{orbits.index(orbit(row, j))}" if t == "pin" else t
            for j, t in enumerate(tokens)
        ]
        line = '  "' + " ".join(new) + '",'
        row += 1
    out_lines.append(line)
assert row == N
Path(out).write_text("\n".join(out_lines) + "".join(materials) + "\n" + rest)
