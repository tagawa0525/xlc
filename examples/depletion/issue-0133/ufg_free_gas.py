"""ADR-0203: UFG の自由ガスの散乱核の測定（ufg-free-gas/）の要約を書く。

使い方:
    python3 ufg_free_gas.py summarize <設定 (TOML)> <JSON>

設定は `candidates`（OpenMC の実行の一覧。swap-candidates.toml。設定のディレクトリから）、`[[mode]]`（`key`、`label`。最初は
`static`）、`[[state]]`（`name` は swap-candidates の状態の名前、`key` はログの名前）。ログは設定のディレクトリの
`logs/<状態>-<散乱核>-absorption.txt`（swap_absorption の例の出力。無い散乱核は、その状態の結果に入れない）と、あれば `logs/<状態>-<散乱核>-fine.txt`（fine_compare.py の出力）。

k の差は xlc − OpenMC（OpenMC の実行の k の平均。pcm）、変化は同じ状態の静止した標的からの k の変化（pcm）。置き換えの寄与は
swap_absorption の k の変化の符号を変えたもの（pcm。swap_candidates.py と同じ）。傾きは fine_compare.py の区間ごとの比 − 1。
"""

import json
import re
import sys
from pathlib import Path

import tomllib

ROOT = Path(__file__).resolve().parents[3]


class SummaryError(Exception):
    pass


def parse_absorption(text: str) -> tuple[float, dict[str, float]]:
    """swap_absorption のログから、置き換えない k と、置き換えごとの寄与（k の変化の符号を変えたもの、pcm）"""
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines or not lines[0].startswith("k "):
        raise SummaryError("the absorption log does not start with the unswapped k")
    k = float(lines[0].split()[1])
    swaps = {}
    for line in lines[1:]:
        m = re.fullmatch(r"(.+): k \d+\.\d+ \(([+-]\d+\.\d+) pcm\)", line)
        if m is None:
            raise SummaryError(f"cannot read the absorption log line: {line}")
        if m.group(1) in swaps:
            raise SummaryError(f"the absorption log repeats `{m.group(1)}`")
        swaps[m.group(1)] = -float(m.group(2))
    return k, swaps


def parse_fine(text: str) -> list[dict]:
    """fine_compare.py の出力から、区間ごとの下端、上端 (eV) と比 − 1"""
    bins = [
        {
            "low": float(m.group(1)),
            "high": float(m.group(2)),
            "ratio": float(m.group(3)),
        }
        for m in re.finditer(
            r"^  ([\d.]+)-([\d.]+) eV: xlc [\d.]+ OpenMC [\d.]+ 比 ([+-][\d.]+)$",
            text,
            re.MULTILINE,
        )
    ]
    if not bins:
        raise SummaryError("the fine log has no bins")
    return bins


def summarize_state(
    name: str,
    key: str,
    openmc: list[float],
    modes: list[str],
    logs: dict[str, tuple[str, str | None]],
) -> dict:
    """1 つの状態の散乱核ごとの結果。`logs` は散乱核ごとの（swap_absorption のログ、fine_compare.py の出力か None）"""
    if not modes or modes[0] != "static":
        raise SummaryError(f"{name}: the first mode must be `static`")
    if not openmc:
        raise SummaryError(f"{name}: no OpenMC runs")
    k_openmc = sum(openmc) / len(openmc)
    results = {}
    specs = None
    for mode in modes:
        if logs.get(mode) is None:
            # 状態によっては流していない散乱核がある（その状態の結果に入れない）
            if mode == "static":
                raise SummaryError(f"{name}: no logs of `static`")
            continue
        absorption, fine = logs[mode]
        k, swaps = parse_absorption(absorption)
        if specs is None:
            specs = list(swaps)
        elif list(swaps) != specs:
            raise SummaryError(f"{name}: `{mode}` swapped another list")
        results[mode] = {
            "k": k,
            "diff_pcm": (k - k_openmc) * 1e5,
            "change_pcm": (k - results.get("static", {"k": k})["k"]) * 1e5,
            "swaps": swaps,
            "tilt": None if fine is None else parse_fine(fine),
        }
    return {
        "name": name,
        "key": key,
        "openmc_k": k_openmc,
        "openmc_runs": openmc,
        "results": results,
    }


def summarize(config: Path) -> dict:
    import msgpack

    base = config.parent
    settings = tomllib.loads(config.read_text())
    candidates = base / settings["candidates"]
    runs = {
        state["name"]: [
            msgpack.unpackb((ROOT / run["openmc"]).read_bytes())["k"]
            for run in state["run"]
        ]
        for state in tomllib.loads(candidates.read_text())["state"]
    }
    modes = [mode["key"] for mode in settings["mode"]]
    states = []
    swap_labels = settings["swap"]
    for state in settings["state"]:
        if state["name"] not in runs:
            raise SummaryError(f"{state['name']} is not in {candidates}")
        logs = {}
        for mode in modes:
            absorption = base / "logs" / f"{state['key']}-{mode}-absorption.txt"
            fine = base / "logs" / f"{state['key']}-{mode}-fine.txt"
            logs[mode] = (
                (absorption.read_text(), fine.read_text() if fine.exists() else None)
                if absorption.exists()
                else None
            )
        s = summarize_state(
            state["name"], state["key"], runs[state["name"]], modes, logs
        )
        unknown = set(s["results"]["static"]["swaps"]) - {
            x["spec"] for x in swap_labels
        }
        if unknown:
            raise SummaryError(f"no labels for the swaps {sorted(unknown)}")
        states.append(s)
    return {
        "format": "xlc-ufg-free-gas",
        "modes": settings["mode"],
        "swaps": swap_labels,
        "states": states,
    }


def main(argv: list[str]) -> int:
    match argv:
        case ["summarize", config, out]:
            try:
                summary = summarize(Path(config))
            except SummaryError as e:
                print(f"FAILED: {e}", file=sys.stderr)
                return 1
            Path(out).write_text(
                json.dumps(summary, ensure_ascii=False, indent=1) + "\n"
            )
            return 0
        case _:
            print(__doc__, file=sys.stderr)
            return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
