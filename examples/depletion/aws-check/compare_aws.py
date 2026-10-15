"""区間に分けた AWS の燃焼の参照解（openmc_depletion_resume.py）と、分けない手元の参照解（openmc_depletion.py）を比べる（ADR-0174 の確認 (c)）。

使い方: python3 compare_aws.py --json <要約の出力先> <手元の結果（msgpack）> <AWS の結果（msgpack）>
        python3 compare_aws.py --report <要約の JSON>

- `--json`: 燃焼度の点と輸送する核種の集合が同じことを確かめてから、燃焼点ごとの k の差 (k_AWS − k_手元) × 1e5 pcm と、
  それを両方の標準偏差の合成 sqrt(σ_AWS² + σ_手元²) で割った倍率、主要な核種の数密度の相対差 (N_AWS − N_手元) / N_手元
  （リングの体積の重みの平均）を書く。`inputs` は読んだファイルのパスと md5
- `--report`: 要約の JSON を文章で描く

同じ種の 2 本は、最初の輸送が同じ（0 の点の k が丸めまで一致する）。それ以降は粒子の並びが変わり、独立な標本とは言えない
ので、倍率は合成した σ に対する目安で、検定の値ではない。数密度の統計誤差は結果に含まれず、相対差に倍率は付けない。
"""

import hashlib
import json
import math
import sys
from pathlib import Path

import msgpack

NUCLIDES = ["U-235", "U-238", "Pu-239", "Xe-135", "Sm-149"]


def openmc_name(name: str) -> str:
    """表の核種名（`U-235`）を OpenMC の名前（`U235`）にする。"""
    return name.replace("-", "")


def summarize(local: dict, aws: dict) -> dict:
    """2 つの結果を比べた要約（モジュールの説明の `--json`）。"""
    for r in (local, aws):
        if r["format"] != "openmc-depletion":
            raise SystemExit(f"unexpected format {r['format']}")
    if set(local["transportable"]) != set(aws["transportable"]):
        raise SystemExit("transported nuclides differ")
    if local["burnup"] != aws["burnup"]:
        raise SystemExit(f"burnup points differ: {local['burnup']} vs {aws['burnup']}")
    if len(local["rings"]) != len(aws["rings"]):
        raise SystemExit("the numbers of rings differ")
    for a, b in zip(local["rings"], aws["rings"], strict=True):
        if a["name"] != b["name"] or a["volume"] != b["volume"]:
            raise SystemExit(
                f"ring {a['name']} differs: {a['volume']} vs {b['volume']}"
            )
    volume = sum(r["volume"] for r in local["rings"])

    def pin(result: dict, nuclide: str, i: int) -> float:
        return (
            sum(r["volume"] * r["densities"][nuclide][i] for r in result["rings"])
            / volume
        )

    k = []
    for i in range(len(local["burnup"])):
        diff = (aws["k"][i] - local["k"][i]) * 1e5
        sigma = math.hypot(aws["k-std"][i], local["k-std"][i]) * 1e5
        k.append(
            {
                "local": local["k"][i],
                "local_std_pcm": local["k-std"][i] * 1e5,
                "aws": aws["k"][i],
                "aws_std_pcm": aws["k-std"][i] * 1e5,
                "diff_pcm": diff,
                "sigma_pcm": sigma,
                "ratio": diff / sigma if sigma else 0.0,
            }
        )
    density = {}
    for name in NUCLIDES:
        series = []
        for i in range(len(local["burnup"])):
            base = pin(local, openmc_name(name), i)
            series.append(
                None if base == 0.0 else pin(aws, openmc_name(name), i) / base - 1.0
            )
        density[name] = series
    return {"burnup": local["burnup"], "k": k, "density": density}


def percent(x: float | None) -> str:
    return "–" if x is None else f"{x * 100:+.4f}"


def report(summary: dict) -> list[str]:
    k = summary["k"]
    lines = [
        "inputs: "
        + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"]),
        "燃焼度 (MWd/kg) | k 手元 ± σ | k AWS ± σ | 差 (pcm) | 合成 σ (pcm) | 倍率",
    ]
    for b, p in zip(summary["burnup"], k, strict=True):
        lines.append(
            f"{b:g} | {p['local']:.6f} ± {p['local_std_pcm'] / 1e5:.6f} | "
            f"{p['aws']:.6f} ± {p['aws_std_pcm'] / 1e5:.6f} | {p['diff_pcm']:+.1f} | "
            f"{p['sigma_pcm']:.1f} | {p['ratio']:+.2f}"
        )
    later = k[1:]
    rms = math.sqrt(sum(p["ratio"] ** 2 for p in later) / len(later))
    lines.append(
        f"0 の点を除く {len(later)} 点: 倍率の二乗平均 {rms:.2f}、最大 "
        f"{max(abs(p['ratio']) for p in later):.2f}、符号が正の点 {sum(p['ratio'] > 0 for p in later)}"
    )
    lines.append("数密度の相対差 (%): 燃焼度 | " + " | ".join(NUCLIDES))
    for i, b in enumerate(summary["burnup"]):
        lines.append(
            f"{b:g} | "
            + " | ".join(percent(summary["density"][n][i]) for n in NUCLIDES)
        )
    return lines


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) == 4 and args[0] == "--json":
        files = [Path(args[2]), Path(args[3])]
        local, aws = (msgpack.unpackb(f.read_bytes()) for f in files)
        summary = {
            "format": "xlc-compare-aws-depletion",
            "inputs": [
                {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
                for f in files
            ],
            **summarize(local, aws),
        }
        Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
