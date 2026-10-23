"""run.sh のログから、自己遮蔽のピンの UFG の時間の並び（ピンごとの開始と終了）と段別の時間を JSON にまとめる。

使い方: pins.py <ログ>... > pins.json
ピンの区間（`pin{cell=…}`）は 1 つのスレッドで開いたままなので、開始は閉じた時刻から busy と idle を引いたもの。
`span` は最初のピンの開始から最後のピンの終了まで（ピンの並列の区間の経過時間）。`stages` は `UFG slowing down` のログの
段別の時間の和で、ピンの中の並列を待つ間に別の仕事を盗んだ時間を含む（CPU 時間ではない）。
"""

import json
import re
import sys
from datetime import datetime
from pathlib import Path

CLOSE = re.compile(
    r"^(\S+)Z\s+\w+\s+.*\bpin\{cell=(\d+)\}: xlc_\w+::[\w:]+: close "
    r"time\.busy=([\d.]+)(ms|µs|ns|s) time\.idle=([\d.]+)(ms|µs|ns|s)"
)
UFG = re.compile(r"\bpin\{cell=(\d+)\}.*UFG slowing down(.*)$")
FIELD = re.compile(r"(\w+)=([\d.eE+-]+)")
UNIT = {"s": 1.0, "ms": 1e-3, "µs": 1e-6, "ns": 1e-9}
STAGES = ["lookup", "grid", "table", "evaluate", "collision", "slowing"]


def pins(log: Path) -> dict:
    spans = []
    stages = dict.fromkeys(STAGES, 0.0)
    for line in log.read_text().splitlines():
        m = CLOSE.search(line)
        if m:
            end = datetime.fromisoformat(m.group(1)).timestamp()
            length = (
                float(m.group(3)) * UNIT[m.group(4)]
                + float(m.group(5)) * UNIT[m.group(6)]
            )
            spans.append((int(m.group(2)), end - length, end))
            continue
        m = UFG.search(line)
        if m:
            fields = dict(FIELD.findall(m.group(2)))
            for s in STAGES:
                stages[s] += float(fields.get(s, 0.0))
    if not spans:
        return {"pins": 0}
    first = min(s for _, s, _ in spans)
    lengths = sorted(e - s for _, s, e in spans)
    return {
        "pins": len(spans),
        "span": round(max(e for _, _, e in spans) - first, 3),
        "sum": round(sum(lengths), 3),
        "longest": round(lengths[-1], 3),
        "median": round(lengths[len(lengths) // 2], 3),
        "stages": {k: round(v, 3) for k, v in stages.items()},
        # （セル、開始、終了）。時刻は最初のピンの開始から
        "timeline": [
            [c, round(s - first, 3), round(e - first, 3)]
            for c, s, e in sorted(spans, key=lambda x: x[1])
        ],
    }


def main() -> None:
    out = {Path(p).stem: pins(Path(p)) for p in sys.argv[1:]}
    json.dump(out, sys.stdout, indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
