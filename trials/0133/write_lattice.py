"""#0133 の試行: OpenMC の 1C の燃焼の結果の燃焼点 <番号> の燃料の組成（輸送する核種）で、xlc の格子の入力を書く。

使い方: python3 write_lattice.py <OpenMC の結果（msgpack）> <番号> <テンプレート> <出力>
"""

import re
import sys
from pathlib import Path

import msgpack

omc = msgpack.unpackb(Path(sys.argv[1]).read_bytes())
i = int(sys.argv[2])
ring = omc["rings"][0]


def xlc_name(name):
    symbol, mass, isomer = re.fullmatch(r"([A-Z][a-z]?)(\d+)(_m1)?", name).groups()
    return f"{symbol}-{mass}" + ("m" if isomer else "")


densities = {
    n: ring["densities"][n][i]
    for n in omc["transportable"]
    if ring["densities"][n][i] > 0.0
}
composition = ", ".join(f'"{xlc_name(n)}" = {v:.8e}' for n, v in densities.items())
text = Path(sys.argv[3]).read_text()
text = re.sub(
    r"(\[material\.fuel\]\ntemperature = \d+\n)nuclides = \{[^}]*\}",
    lambda m: m.group(1) + f"nuclides = {{ {composition} }}",
    text,
)
Path(sys.argv[4]).write_text(text)
