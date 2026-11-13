"""#0286 の衝突確率の近似の時間（排他のジョブのログ `timing-<ケース>-<変種>.log`、変種は `cp-` で始まる）の要約を書く。

ケースと変種ごとに、自己遮蔽の span の秒、1 状態の経過時間、1 ピンの UFG の段の和（スレッド秒）と移行の区間を compare.py の
読み方で読む。

使い方: python3 timing.py --json <要約の出力先> <ケース>...（このディレクトリで）
"""

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


COMPARE = load("compare", HERE / "compare.py")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("cases", nargs="+")
    args = parser.parse_args()
    summary = {}
    for case in args.cases:
        prefix = f"timing-{case}-"
        for path in sorted(HERE.glob(f"{prefix}cp-*.log")):
            timing = COMPARE.log(path)
            summary.setdefault(case, {})[path.stem[len(prefix) :]] = {
                k: timing[k]
                for k in (
                    "shielding_seconds",
                    "total_seconds",
                    "per_pin",
                    "transfers_per_pin",
                )
            }
    args.json.write_text(json.dumps(summary, indent=1, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
