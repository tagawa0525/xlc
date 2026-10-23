"""run.sh のログから、区間（span）ごとの busy の合計と回数を JSON にまとめる。

使い方: summarize.py <ログのディレクトリ> > scaling.json
busy は区間を開いていたスレッドの時間なので、並列の区間の中の区間（pin、lump など）の合計は、
各スレッドの時間の和（CPU 時間に近い）で、経過時間ではない。経過時間は親の区間（shielding、moc）の busy。
"""

import json
import re
import sys
from collections import defaultdict
from pathlib import Path

CLOSE = re.compile(
    r"Z\s+[A-Z]+\s+(.+?): xlc_\w+::[\w:]+: close time\.busy=([\d.]+)(ms|µs|ns|s)\b"
)
UNIT = {"s": 1.0, "ms": 1e-3, "µs": 1e-6, "ns": 1e-9}


def path_of(span: str) -> str:
    return re.sub(r"\{[^}]*\}", "", span)


def summarize(log: Path) -> dict:
    spans: dict[str, list[float]] = defaultdict(lambda: [0.0, 0])
    cmfd = []
    for line in log.read_text().splitlines():
        m = CLOSE.search(line)
        if m:
            entry = spans[path_of(m.group(1))]
            entry[0] += float(m.group(2)) * UNIT[m.group(3)]
            entry[1] += 1
        elif "coarse" in line and "iterations" in line:
            cmfd.append(line)
    return {
        "spans": {
            k: {"busy": round(v[0], 4), "count": v[1]} for k, v in sorted(spans.items())
        },
        "cmfd_lines": len(cmfd),
    }


def main() -> None:
    logs = Path(sys.argv[1])
    wall = {}
    for line in (logs / "wall.txt").read_text().splitlines():
        name, _, seconds = line.split()
        wall[name] = float(seconds)
    out = {
        log.stem: {"wall": wall.get(log.stem), **summarize(log)}
        for log in sorted(logs.glob("*.log"))
    }
    json.dump(out, sys.stdout, indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
