"""ADR-0207: UFG の範囲の群の散乱行列の移行の形を、無限希釈とリングの中性子束で数えた移行にした測定（ufg-transfer/）の要約を書く。

使い方:
    python3 ufg_transfer.py summarize <設定 (TOML)> <JSON>

設定は `[[spec]]`（`spec` は `<物質>:<最初の群>-<最後の群>`、`label`。swap_transfer の置き換え）、`[[kernel]]`（`key`、`label`。
group_rows.py rows の核）、`[[mode]]`（`key`、`label`。移行の形。最初は無限希釈）、`[[state]]`（`key`、`name`、`kernels` は核の名前の
並び、`rows` は rows の置き換え。省くと rows を比べない）。ログは設定のディレクトリの `logs/<状態>-<形>-transfer.txt`
（swap_transfer の出力に `material <物質> <群>` を前に置いたもの）と `logs/<状態>-<形>-rows.txt`（group_rows.py rows の出力に
`kernel <名前> material <物質> <群>` を前に置いたもの）。

k の差は xlc − OpenMC（swap_transfer の OpenMC の k。状態の移行を数えた実行の k、pcm）、置き換えの寄与は swap_transfer の SPH 因子を
置き換える前のままの k の変化の符号を変えたもの（pcm。group_rows.py と同じ）。同じ状態の形の間で OpenMC の k が同じことを確かめる。
"""

import json
import sys
from pathlib import Path

import tomllib
from group_rows import SummaryError, contributions, groups, parse_rows


def summarize(config):
    config = Path(config)
    c = tomllib.loads(config.read_text())
    logs = config.parent / "logs"
    specs = [s["spec"] for s in c["spec"]]
    modes = [m["key"] for m in c["mode"]]
    states = []
    for state in c["state"]:
        key = state["key"]
        results = {}
        rows = {}
        openmc_k = None
        for mode in modes:
            k, openmc_k, swaps = contributions(
                logs / f"{key}-{mode}-transfer.txt", openmc_k=openmc_k
            )
            if set(swaps) != set(specs):
                raise SummaryError(f"{key} {mode}: the swaps are {sorted(swaps)}")
            results[mode] = {"k": k, "swaps": swaps}
            wanted = state.get("rows", [])
            found = parse_rows((logs / f"{key}-{mode}-rows.txt").read_text())
            if set(found) != (set(state["kernels"]) if wanted else set()):
                raise SummaryError(
                    f"{key} {mode}: the rows log has kernels {sorted(found)}"
                )
            expected = {
                f"{spec.split(':')[0]}:{g + 1}"
                for spec in wanted
                for g in groups(*spec.split(":")[1].split("-"))
            }
            for kernel, entries in found.items():
                if set(entries) != expected:
                    raise SummaryError(
                        f"{key} {mode}: the {kernel} rows are {sorted(entries)}, not {sorted(expected)}"
                    )
            rows[mode] = found
        states.append(
            {
                "key": key,
                "name": state["name"],
                "openmc_k": openmc_k,
                "kernels": state["kernels"] if state.get("rows") else [],
                "modes": results,
                "rows": rows,
            }
        )
    return {
        "format": "xlc-ufg-transfer",
        "specs": c["spec"],
        "kernels": c["kernel"],
        "modes": c["mode"],
        "states": states,
    }


def write_summary(config, output):
    Path(output).write_text(
        json.dumps(summarize(config), ensure_ascii=False, indent=1) + "\n"
    )


if __name__ == "__main__":
    if len(sys.argv) != 4 or sys.argv[1] != "summarize":
        raise SystemExit(__doc__)
    write_summary(*sys.argv[2:])
