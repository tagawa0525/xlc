"""分離共鳴を持たない核種の自由ガスの上限を変えた VERA の 1B〜1D、2A〜2K、2Q を、仕様書の KENO-VI（ENDF/B-VII.0）と比べる
（ADR-0246、#0274）。

使い方: python3 vera_light_cap.py [--json <要約の出力先>] <置き場> <spec.txt> <記録の要約のディレクトリ>

置き場は <問題>-<条件>.toml を lattice で解いた <問題>-<条件>.out（集合体は pins-<問題>-<条件>.msgpack も）のあるディレクトリ。
条件は e1000（既定。ADR-0246）、400（分離共鳴を持たない核種だけ前の既定の 400 kT）、400-heavy400（両方 400 kT。ADR-0243 の前の
既定で、記録の条件）。記録の要約のディレクトリは examples/vera/endfb70（compare-<問題>.json。ADR-0176）。k の差とピン出力は
../../../vera/compare_keno_lattice.py の summary と同じ読み方（差は (k_xlc − k_KENO-VI) × 1e5 pcm）。

問題ごとに、KENO-VI の標準偏差、記録の差、条件ごとの差、400-heavy400 と記録の差（同じライブラリと条件の再現）、e1000 と 400 の
差（上限の変更が動かした量）、ピン出力の差の RMS (%) を出し、最後にピンセルと集合体それぞれの差の平均と RMS を条件ごとに出す。

解き方（このディレクトリで。Slurm のジョブのコマンド）:
for p in 1B 1C 1D 2A 2B 2C 2D 2E 2F 2G 2H 2I 2J 2K 2Q; do for v in e1000 400 400-heavy400; do ../../../../bin/lattice $p-$v.toml > $p-$v.out 2> $p-$v.log || exit 1; done; done && python3 vera_light_cap.py --json summary.json . ../../../../references/vera/spec.txt ../../../vera/endfb70
"""

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "vera"))
from compare_keno_lattice import summary

PINCELLS = ["1B", "1C", "1D"]
ASSEMBLIES = ["2A", "2B", "2C", "2D", "2E", "2F", "2G", "2H", "2I", "2J", "2K", "2Q"]
VARIANTS = ["e1000", "400", "400-heavy400"]


def load(directory: Path, spec: list[str], records: Path) -> tuple[dict, dict]:
    """(問題, 条件) ごとの要約と、問題ごとの記録の要約"""
    results = {}
    for problem in PINCELLS + ASSEMBLIES:
        for variant in VARIANTS:
            pins = (
                directory / f"pins-{problem}-{variant}.msgpack"
                if problem in ASSEMBLIES
                else None
            )
            out = directory / f"{problem}-{variant}.out"
            results[problem, variant] = summary(spec, problem, out, pins)
    recorded = {
        p: json.loads((records / f"compare-{p}.json").read_text())
        for p in PINCELLS + ASSEMBLIES
    }
    return results, recorded


def rms_percent(s: dict) -> float | None:
    return None if s["pins"] is None else s["pins"]["rms"] * 100


def problem_rows(results: dict, recorded: dict) -> list[dict]:
    rows = []
    for problem in PINCELLS + ASSEMBLIES:
        diff = {v: results[problem, v]["k"]["diff_pcm"] for v in VARIANTS}
        rows.append(
            {
                "problem": problem,
                "keno_std_pcm": recorded[problem]["k"]["keno_std"] * 1e5,
                "record_diff_pcm": recorded[problem]["k"]["diff_pcm"],
                "diff_pcm": diff,
                "reproduction_pcm": diff["400-heavy400"]
                - recorded[problem]["k"]["diff_pcm"],
                "shift_pcm": diff["e1000"] - diff["400"],
                "record_pins_rms_percent": rms_percent(recorded[problem]),
                "pins_rms_percent": {
                    v: rms_percent(results[problem, v]) for v in VARIANTS
                },
            }
        )
    return rows


def aggregate(rows: list[dict]) -> dict:
    """ピンセルと集合体それぞれの、条件ごとの差の平均と RMS (pcm)"""
    out = {}
    for name, problems in (("pincells", PINCELLS), ("assemblies", ASSEMBLIES)):
        selected = [r for r in rows if r["problem"] in problems]
        out[name] = {
            v: {
                "mean_pcm": sum(r["diff_pcm"][v] for r in selected) / len(selected),
                "rms_pcm": math.sqrt(
                    sum(r["diff_pcm"][v] ** 2 for r in selected) / len(selected)
                ),
            }
            for v in VARIANTS
        }
    return out


def text(rows: list[dict], totals: dict) -> list[str]:
    def pins(value: float | None) -> str:
        return f"{'-':>7}" if value is None else f"{value:7.3f}"

    lines = [
        (
            "問題  σ(KENO)  記録  400-heavy400  再現  400  e1000  e1000−400 (pcm)  "
            "ピン出力の RMS (%): 記録  400-heavy400  400  e1000"
        )
    ]
    for r in rows:
        d = r["diff_pcm"]
        p = r["pins_rms_percent"]
        lines.append(
            f"{r['problem']:4} {r['keno_std_pcm']:6.1f} {r['record_diff_pcm']:+8.1f} "
            f"{d['400-heavy400']:+8.1f} {r['reproduction_pcm']:+6.1f} {d['400']:+8.1f} "
            f"{d['e1000']:+8.1f} {r['shift_pcm']:+8.1f}  {pins(r['record_pins_rms_percent'])} "
            f"{pins(p['400-heavy400'])} {pins(p['400'])} {pins(p['e1000'])}"
        )
    for name, values in totals.items():
        cells = "  ".join(
            f"{v} {values[v]['mean_pcm']:+.1f} / {values[v]['rms_pcm']:.1f}"
            for v in VARIANTS
        )
        lines.append(f"{name} の差の平均 / RMS (pcm): {cells}")
    return lines


def main() -> None:
    args = sys.argv[1:]
    output = None
    if args[:1] == ["--json"]:
        if len(args) < 2:
            raise SystemExit(__doc__)
        output, args = Path(args[1]), args[2:]
    if len(args) != 3:
        raise SystemExit(__doc__)
    spec = Path(args[1]).read_text().split("\n")
    results, recorded = load(Path(args[0]), spec, Path(args[2]))
    rows = problem_rows(results, recorded)
    totals = aggregate(rows)
    print("\n".join(text(rows, totals)))
    if output is not None:
        result = {
            "format": "xlc-0274-vera-light-cap",
            "problems": rows,
            "totals": totals,
        }
        output.write_text(json.dumps(result, indent=1) + "\n")


if __name__ == "__main__":
    main()
