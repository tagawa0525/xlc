"""core の燃焼のログ（tracing の標準エラー）から、時間の内訳と固有値計算の回数を JSON に取る（#0156）。

使い方: python3 p9_timing.py <core のログ> > <JSON>

- 時間（秒）: `total` は core の全体（`core seconds=`）。`eigenvalue`（span `nem` のフィールド `seconds` の和）、`interpolation`（span `constants` の和。
  テーブルの内挿とミクロ燃焼の補正）、`micro`（span `micro depletion` と `micro decay` の和。反応率と CRAM）、`xenon`（span
  `xenon balance` と `xenon anderson` の和。平衡の Xe と I の値と Anderson の加速）は互いに重ならない部分で、`other` は全体からそれらを引いた残り（入力の読み込み、出力など）。`search`（span
  `search` の和）は臨界ホウ素の探索の全体で、その中の固有値計算と内挿を含む（ほかの部分と重なる）
- `nem`: 固有値計算の span の数と、外側の反復の回数と時間の最小と最大
- `solves`: ログの順の燃焼点と中点ごとの EFPD（中点はそのステップの終わり）、ホウ素濃度、固有値計算の回数（平衡 Xe の解き直しを
  含む）、平衡 Xe の解き直しの回数（平衡 Xe でない解は null）
"""

import json
import re
import sys

# span の閉じる行: `<時刻>  INFO <span の並び>: <target>: close time.busy=<時間> time.idle=<時間>`。span の名前は空白を含みうる
CLOSE = re.compile(
    r"INFO (?P<spans>.+?): (?P<target>[\w:]+): close time\.busy=(?P<busy>\S+) "
)
EVENT = re.compile(
    r"INFO (?:.+?: )?(?P<target>[\w:]+): (?P<message>[a-z ]+?) (?P<fields>\w+=.*)$"
)
# 1 秒あたりの単位の数（掛けずに割る: 12 × 1e-9 は 1.2e-8 に丸まらない）
PER_SECOND = {"ns": 1e9, "µs": 1e6, "ms": 1e3, "s": 1.0}

# span の名前 → 部分。`nem` は time.busy でなくフィールド `seconds` を足す（time.busy は 3 桁に丸まる）
PARTS = {
    "nem": "eigenvalue",
    "constants": "interpolation",
    "micro depletion": "micro",
    "micro decay": "micro",
    "xenon balance": "xenon",
    "xenon anderson": "xenon",
}


def seconds(text: str) -> float:
    """tracing の時間（`3.45µs`、`150ms`、`118s` など）を秒に"""
    m = re.fullmatch(r"([0-9.]+)(ns|µs|ms|s)", text)
    if not m:
        raise ValueError(f"not a duration: {text}")
    return float(m.group(1)) / PER_SECOND[m.group(2)]


def leaf(spans: str) -> tuple[str, dict[str, str]]:
    """span の並び（`a{x=1}:b:c{y=2}`）の最後の span の名前とフィールド"""
    m = re.search(r"(?:^|:)([^:{}]+)(?:\{([^{}]*)\})?$", spans)
    if not m:
        raise ValueError(f"not a span list: {spans}")
    fields = dict(f.split("=", 1) for f in (m.group(2) or "").split() if "=" in f)
    return m.group(1), fields


def summarize(log: str) -> dict:
    parts = {name: 0.0 for name in PARTS.values()}
    search = 0.0
    nem_iterations: list[int] = []
    nem_seconds: list[float] = []
    solves = []
    # 平衡 Xe の解き直しの回数は、その解の点か中点の行の前に出る
    xenon = None
    total = None
    for line in log.splitlines():
        if m := CLOSE.search(line):
            name, fields = leaf(m.group("spans"))
            busy = seconds(m.group("busy"))
            if name == "nem" and "seconds" in fields:
                busy = float(fields["seconds"])
            if name in PARTS:
                parts[PARTS[name]] += busy
            if name == "search":
                search += busy
            if name == "nem":
                nem_seconds.append(busy)
                if "iterations" in fields:
                    nem_iterations.append(int(fields["iterations"]))
            continue
        m = EVENT.search(line)
        if not m:
            continue
        fields = dict(f.split("=", 1) for f in m.group("fields").split() if "=" in f)
        message = m.group("message")
        if message in ("depletion point", "depletion midpoint"):
            if message == "depletion point":
                kind, efpd = "point", float(fields["efpd"])
            else:
                # 中点の行は EFPD を持たないので、ステップの span の終わりの EFPD
                step = re.search(r"depletion step\{([^{}]*)\}", line)
                step_fields = dict(
                    f.split("=", 1)
                    for f in (step.group(1) if step else "").split()
                    if "=" in f
                )
                if "efpd" not in step_fields:
                    raise ValueError(f"a midpoint outside a depletion step: {line}")
                kind, efpd = "midpoint", float(step_fields["efpd"])
            solves.append(
                {
                    "kind": kind,
                    "efpd": efpd,
                    "boron": float(fields["boron"]),
                    "solves": int(fields["solves"]),
                    "xenon": xenon,
                }
            )
            xenon = None
        elif message == "equilibrium xenon" and "iterations" in fields:
            xenon = int(fields["iterations"])
        elif message == "core" and "seconds" in fields:
            total = float(fields["seconds"])
    if total is None:
        raise ValueError("the log has no total time (core seconds=)")
    if not nem_seconds:
        raise ValueError("the log has no eigenvalue solves (span nem)")
    return {
        "format": "xlc-core-timing",
        "seconds": {
            "total": total,
            **parts,
            "other": total - sum(parts.values()),
            "search": search,
        },
        "nem": {
            "count": len(nem_seconds),
            "iterations": [min(nem_iterations), max(nem_iterations)]
            if nem_iterations
            else None,
            "seconds": [min(nem_seconds), max(nem_seconds)],
        },
        "solves": solves,
    }


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit("usage: python3 p9_timing.py <core log>")
    with open(sys.argv[1], encoding="utf-8") as f:
        summary = summarize(f.read())
    json.dump(summary, sys.stdout, indent=2, ensure_ascii=False)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
