"""冷却材の B-10 の割合がテーブルと違う状態の補正（ADR-0236）を、直接の格子計算の k∞ と比べる。

使い方: python3 compare_core_b10.py --json <要約の出力先> <テーブル（msgpack）> <core のケース（toml）> <core の標準出力>
          <直接の格子計算のテーブル（msgpack）>...
        python3 compare_core_b10.py --report <要約の JSON>

- `--json`: core のケース（branches/core-2A-b10.toml の形）は、集合体 1 体で、状態の名前が
  `corrected-b<割合 ×1000>-p<ppm>-u<燃焼度>`（状態の b10 で補正）、`converted-b…-p…-u…`（テーブルの b10 で、B-10 の個数密度が同じになる
  ppm に直したもの。ADR-0236 の案 1）、`nominal-p<ppm>-u<燃焼度>`（テーブルの b10 の同じ ppm）。直接の格子計算のテーブルごとに、その
  b10 とホウ素の軸の点、燃焼度の軸の点の組の状態（減速材と燃料はテーブルの 1 点）を比べ、次を書く。
  - `direct`: 直接の格子計算のテーブルの k∞
  - `corrected`、`converted`: core の k の、直接の格子計算からの反応度の差 Δρ = (1/k_直接 − 1/k) × 1e5 pcm
  - `effect`: B-10 の割合を変えた効き（nominal の k から直接の格子計算の k への反応度の差）
  ケースの状態の b10、ppm、燃焼度、温度が名前とテーブルに合うこと（converted の ppm は ADR-0106 の原子質量で換算した値と相対 1e-12）を
  確かめる。`inputs` は読んだファイルのパスと md5
- `--report`: 要約の JSON を文章で描く
"""

import hashlib
import json
import re
import sys
from pathlib import Path

import msgpack
import tomllib

# AME2020 の原子質量 (u)（ADR-0106。crates/format/src/table.rs の MASS_B10、MASS_B11）
MASS_B10 = 10.01293695
MASS_B11 = 11.00930536
STATE = re.compile(r"(corrected|converted)-b(\d+)-p(\d+)-u(\d+)|nominal-p(\d+)-u(\d+)")
K_LINE = re.compile(r"k-eff (\S+)")
CONVERTED_TOLERANCE = 1e-12
# テーブルの格子点の状態の core の k と k∞ の差の許容 (pcm)
NOMINAL_TOLERANCE = 0.01


def per_mass(b10: float) -> float:
    """ホウ素の単位質量あたりの B-10 の原子の数（比だけを使う）"""
    return b10 / (b10 * MASS_B10 + (1.0 - b10) * MASS_B11)


def rho_pcm(base: float, k: float) -> float:
    """`base` から `k` への反応度の差 (pcm)"""
    return (1.0 / base - 1.0 / k) * 1e5


def core_ks(out: str) -> dict[str, float]:
    """core の標準出力の、状態の名前 → k（`state <名前>` に続く `k-eff` の行）"""
    found: dict[str, float] = {}
    name = None
    for line in out.splitlines():
        if line.startswith("state "):
            name = line.removeprefix("state ")
        elif name is not None and (m := K_LINE.fullmatch(line)):
            found[name] = float(m[1])
            name = None
    return found


def axes_of(table: dict) -> dict[str, dict]:
    return {a["name"]: a for a in table["axes"]}


def table_k(table: dict, ppm: float, burnup: float) -> float:
    """テーブルの、ホウ素 `ppm` と燃焼度 `burnup` の点の k∞（ほかの軸は 1 点）"""
    axes = table["axes"]
    at = []
    for a in axes:
        if a["name"] == "boron":
            at.append(a["points"].index(ppm))
        elif a["name"] == "burnup":
            at.append(a["points"].index(burnup))
        elif len(a["points"]) == 1:
            at.append(0)
        else:
            raise SystemExit(f"axis {a['name']} must have one point")
    return table["states"][table["points"].index(at)]["k_inf"]


def summarize(table: dict, case: dict, out: str, directs: list[dict]) -> dict:
    """要約（モジュールの説明の `--json`）"""
    axes = axes_of(table)
    table_b10 = axes["boron"]["b10"]
    (moderator,) = axes["moderator"]["points"]
    (fuel,) = axes["fuel"]["points"]
    states = case["state"]
    ks = core_ks(out)
    if set(ks) != set(states):
        raise SystemExit(
            f"the output has states {sorted(ks)}, the case {sorted(states)}"
        )
    by: dict[tuple, float] = {}
    for name, s in states.items():
        m = STATE.fullmatch(name)
        if m is None:
            raise SystemExit(f"state {name}: unexpected name")
        kind, b, ppm, u = (
            (m[1], int(m[2]) / 1000, float(m[3]), float(m[4]))
            if m[1]
            else ("nominal", table_b10, float(m[5]), float(m[6]))
        )
        if (s["moderator"], s["fuel"], s["burnup"]) != (moderator, fuel, u):
            raise SystemExit(
                f"state {name}: temperatures or burnup differ from its name and the table"
            )
        if kind == "converted":
            want = ppm * per_mass(b) / per_mass(table_b10)
            if (
                s["b10"] != table_b10
                or abs(s["boron"] - want) > CONVERTED_TOLERANCE * want
            ):
                raise SystemExit(
                    f"state {name}: boron {s['boron']} b10 {s['b10']}, want {want} at {table_b10}"
                )
        elif (s["b10"], s["boron"]) != ((b if kind == "corrected" else table_b10), ppm):
            raise SystemExit(f"state {name}: b10 or boron differ from its name")
        by[(kind, b, ppm, u)] = ks[name]
    points = []
    for direct in directs:
        d_axes = axes_of(direct)
        b10 = d_axes["boron"]["b10"]
        if (d_axes["moderator"]["points"], d_axes["fuel"]["points"]) != (
            [moderator],
            [fuel],
        ):
            raise SystemExit(f"the direct table at b10 {b10} has other temperatures")
        for ppm in d_axes["boron"]["points"]:
            for u in d_axes["burnup"]["points"]:
                k = table_k(direct, ppm, u)
                nominal = by[("nominal", table_b10, ppm, u)]
                # 補正のない状態の core の k は、テーブルの格子点の k∞ と同じ（ADR-0190 の確認の 0.00043 pcm の水準）
                if abs(rho_pcm(table_k(table, ppm, u), nominal)) > NOMINAL_TOLERANCE:
                    raise SystemExit(
                        f"the nominal core k at {ppm} ppm {u} MWd/kg is not the table's"
                    )
                points.append(
                    {
                        "b10": b10,
                        "boron": ppm,
                        "burnup": u,
                        "direct": k,
                        "effect": rho_pcm(nominal, k),
                        "corrected": rho_pcm(k, by[("corrected", b10, ppm, u)]),
                        "converted": rho_pcm(k, by[("converted", b10, ppm, u)]),
                    }
                )
    return {"table_b10": table_b10, "points": points}


def report(summary: dict) -> list[str]:
    p = summary["points"]
    lines = [
        f"テーブルの b10 {summary['table_b10']}、{len(p)} 点",
        "| b10 | ppm | 燃焼度 | 直接の k∞ | 効き (pcm) | 補正 − 直接 (pcm) | 換算 − 直接 (pcm) |",
    ]
    lines += [
        f"| {x['b10']} | {x['boron']:g} | {x['burnup']:g} | {x['direct']:.6f} | {x['effect']:+.1f} | "
        f"{x['corrected']:+.2f} | {x['converted']:+.2f} |"
        for x in p
    ]
    for key in ["effect", "corrected", "converted"]:
        values = [x[key] for x in p]
        lines.append(f"{key}: {min(values):+.3f}〜{max(values):+.3f} pcm")
    return lines


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) >= 6 and args[0] == "--json":
        files = [Path(a) for a in args[2:]]

        def read(f: Path) -> dict:
            return msgpack.unpackb(f.read_bytes())

        summary = {
            "format": "xlc-compare-core-b10",
            # 読んだファイルのパス（実行した cwd から）と md5
            "inputs": [
                {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
                for f in files
            ],
            **summarize(
                read(files[0]),
                tomllib.loads(files[1].read_text()),
                files[2].read_text(),
                [read(f) for f in files[3:]],
            ),
        }
        Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
