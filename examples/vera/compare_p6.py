"""VERA Problem 6 の xlc の炉心の解（熱水力のフィードバック。ADR-0229）を、確認の量の要約にする。

使い方: python3 compare_p6.py --json <要約の出力先> <格子の TOML> <HFP の core の TOML> <出力を下げた core の TOML>
        <熱水力なしの core の TOML> <HFP の core の標準出力> <出力を下げた core の標準出力> <熱水力なしの core の標準出力>
        python3 compare_p6.py --report <要約の JSON>

Problem 6 には参照解が無い（仕様書の REFERENCE SOLUTION）ので、比べるのは次の確認の量。結果は JSON だけに書き、JSON は
検証の記録の表を描くツール（examples/verification_tables.py）も読む（ADR-0176）。

- 流路のエンタルピーの上昇と、出力 / 流量（core の TOML の `[thermal-hydraulics]` の power と flow）。相対差
- 熱水力の繰り返しの回数と最後の変化と、TOML の `tolerance`（収束しなければ core がエラーで止まるので、判定はしない。
  出力の変化は 4 桁に丸めてあり、許容の近くでは判定に使えない）
- 出力を下げた解（ノードの温度はほぼ入口の温度）と、熱水力なしの解の k の差 (pcm)。出力の比（HFP の power との比）と、熱水力なしの
  解の状態の減速材と燃料の温度も、表の行の名前に使うので TOML から取る（実効温度の方式も同じ）
- 層ごとの減速材の温度、燃料の実効温度（核分裂のある層だけ）、核分裂率の割合と、減速材の温度でのテーブルの密度
  （格子の TOML の `[branches]` の moderator の格子点の密度を温度で線形に内挿したもの。core の内挿と同じ。ADR-0229 の案 T1）
- `--report`: JSON から文章の報告を標準出力に描く。1 行目は読んだファイルのパスと md5 の先頭 8 桁（`inputs`）
"""

import hashlib
import itertools
import json
import re
import sys
from pathlib import Path

import tomllib

NUMBER = r"[-+0-9.eE]+"
# 熱水力の要約の、`thermal-hydraulics iterations` の行の後に出る行の頭
THERMAL_LINES = (
    "fuel temperature",
    "moderator density",
    "channel outlets",
    "thermal layers",
)


def parse_core(text: str) -> dict:
    """core の標準出力（状態が 1 つ）の k、核分裂率の層ごとの割合、熱水力の要約（無ければ None）"""
    lines = text.splitlines()
    k = [
        float(m.group(1))
        for line in lines
        if (m := re.match(rf"^k-eff ({NUMBER})$", line))
    ]
    if len(k) != 1:
        raise SystemExit(f"{len(k)} k-eff lines (one state is expected)")
    fractions = []
    thermal = None
    section = None
    for line in lines:
        if line.startswith("axial fission fractions"):
            section = "fractions"
            continue
        if line.startswith("assembly fission rates"):
            section = None
            continue
        if m := re.match(
            rf"^thermal-hydraulics iterations (\d+) change ({NUMBER})$", line
        ):
            thermal = {
                "iterations": int(m.group(1)),
                "change": float(m.group(2)),
                "channels": [],
                "layers": [],
            }
            section = None
            continue
        if thermal is None and line.startswith(THERMAL_LINES):
            raise SystemExit(
                f"{line!r} comes before the thermal-hydraulics iterations line"
            )
        if m := re.match(
            rf"^fuel temperature volume average ({NUMBER}) K max ({NUMBER}) K$", line
        ):
            thermal["fuel_average"], thermal["fuel_max"] = map(float, m.groups())
            continue
        if m := re.match(
            rf"^moderator density interpolation error max ({NUMBER}) g/cc$", line
        ):
            thermal["density_error"] = float(m.group(1))
            continue
        if line.startswith("channel outlets"):
            section = "channels"
            continue
        if line.startswith("thermal layers"):
            section = "layers"
            continue
        fields = line.split()
        if section == "fractions" and len(fields) == 3:
            z0, z1, f = map(float, fields)
            fractions.append({"z": [z0, z1], "fraction": f})
        elif section == "channels" and len(fields) == 4:
            thermal["channels"].append(
                {
                    "position": [int(fields[0]), int(fields[1])],
                    "rise": float(fields[2]),
                    "outlet": float(fields[3]),
                }
            )
        elif section == "layers" and len(fields) == 4:
            thermal["layers"].append(
                {
                    "z": [float(fields[0]), float(fields[1])],
                    "moderator": float(fields[2]),
                    "fuel": None if fields[3] == "-" else float(fields[3]),
                }
            )
        else:
            section = None
    if thermal is not None:
        missing = {"fuel_average", "density_error"} - thermal.keys()
        if missing or not thermal["channels"] or not thermal["layers"]:
            raise SystemExit(
                f"the thermal-hydraulics summary is incomplete ({sorted(missing)})"
            )
    return {"k": k[0], "fractions": fractions, "thermal": thermal}


def table_density(grid: list[tuple[float, float]], t: float) -> float:
    """減速材の格子点 (温度 K, 密度 g/cc) の、温度 t での線形の内挿の密度。格子の外はエラー"""
    for (t0, d0), (t1, d1) in itertools.pairwise(grid):
        if t0 <= t <= t1:
            return d0 + (d1 - d0) * (t - t0) / (t1 - t0)
    raise SystemExit(
        f"{t} K is outside the moderator grid {grid[0][0]}〜{grid[-1][0]} K"
    )


def compare(paths: list[str]) -> dict:
    lattice_toml, core_toml, low_toml, iso_toml, hfp_out, low_out, iso_out = map(
        Path, paths
    )
    moderator = tomllib.loads(lattice_toml.read_text())["branches"]["moderator"]
    grid = list(zip(map(float, moderator["temperature"]), moderator["density"]))
    th = tomllib.loads(core_toml.read_text())["thermal-hydraulics"]
    low_th = tomllib.loads(low_toml.read_text())["thermal-hydraulics"]
    (iso_state,) = tomllib.loads(iso_toml.read_text())["state"].values()
    if "tolerance" not in th:
        raise SystemExit(f"{core_toml}: write the tolerance (the record quotes it)")
    hfp, low, iso = (parse_core(p.read_text()) for p in (hfp_out, low_out, iso_out))
    for name, solve in [("HFP", hfp), ("low power", low)]:
        if solve["thermal"] is None or len(solve["thermal"]["channels"]) != 1:
            raise SystemExit(f"the {name} solve has no single-channel thermal summary")
    if iso["thermal"] is not None:
        raise SystemExit("the isothermal solve has a thermal summary")
    t = hfp["thermal"]
    channel = t["channels"][0]
    # power (MW) × 1e3 = kJ/s、流量 (kg/s) で割って kJ/kg
    expected = th["power"] * 1e3 / th["flow"]
    fractions = {tuple(f["z"]): f["fraction"] for f in hfp["fractions"]}
    layers = []
    for layer in t["layers"]:
        z = tuple(layer["z"])
        if z not in fractions:
            raise SystemExit(f"the thermal layer {z} is not an axial fission layer")
        layers.append(
            {
                **layer,
                "density": table_density(grid, layer["moderator"]),
                "fraction": fractions[z],
            }
        )
    return {
        "format": "xlc-compare-p6",
        # 読んだファイルのパス（ステージの cwd から）と md5。dvc.lock の deps と一致する（record_provenance.py が照合する）
        "inputs": [
            {"path": str(p), "md5": hashlib.md5(p.read_bytes()).hexdigest()}
            for p in (
                lattice_toml,
                core_toml,
                low_toml,
                iso_toml,
                hfp_out,
                low_out,
                iso_out,
            )
        ],
        "k": hfp["k"],
        "thermal": {
            "iterations": t["iterations"],
            "change": t["change"],
            "tolerance": th["tolerance"],
            "effective": th["effective-temperature"],
            "fuel_average": t["fuel_average"],
            "fuel_max": t["fuel_max"],
            "density_error": t["density_error"],
        },
        "channel": {
            "rise": channel["rise"],
            "expected_rise": expected,
            "rise_relative_diff": channel["rise"] / expected - 1,
            "outlet": channel["outlet"],
        },
        "low_power": {
            "k": low["k"],
            "power_ratio": low_th["power"] / th["power"],
            "iterations": low["thermal"]["iterations"],
            "outlet": low["thermal"]["channels"][0]["outlet"],
            "fuel_max": low["thermal"]["fuel_max"],
        },
        "isothermal": {
            "k": iso["k"],
            "moderator": iso_state["moderator"],
            "fuel": iso_state["fuel"],
        },
        "low_power_vs_isothermal_pcm": (low["k"] - iso["k"]) * 1e5,
        "layers": layers,
    }


def report(s: dict) -> list[str]:
    t, c = s["thermal"], s["channel"]
    lines = [
        "inputs: " + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in s["inputs"]),
        f"k {s['k']:.8f}",
        (
            f"thermal-hydraulics: {t['iterations']} iterations, last change {t['change']:.3e} "
            f"(tolerance {t['tolerance']:.1e})"
        ),
        f"fuel temperature: volume average {t['fuel_average']:.2f} K, max {t['fuel_max']:.2f} K",
        f"moderator density interpolation error max {t['density_error']:.3e} g/cc",
        (
            f"enthalpy rise {c['rise']:.6f} kJ/kg, power / flow {c['expected_rise']:.6f} kJ/kg "
            f"(relative {c['rise_relative_diff']:.3e}), outlet {c['outlet']:.4f} K"
        ),
        (
            f"low power: k {s['low_power']['k']:.8f} ({s['low_power']['iterations']} iterations, outlet "
            f"{s['low_power']['outlet']:.6f} K, fuel max {s['low_power']['fuel_max']:.6f} K); isothermal k "
            f"{s['isothermal']['k']:.8f}; difference {s['low_power_vs_isothermal_pcm']:.3f} pcm"
        ),
        "layers (z from, z to (cm), moderator (K), table density (g/cc), fuel (K), fission fraction):",
    ]
    for layer in s["layers"]:
        fuel = "-" if layer["fuel"] is None else f"{layer['fuel']:.2f}"
        lines.append(
            f"  {layer['z'][0]:.4f} {layer['z'][1]:.4f} {layer['moderator']:.4f} "
            f"{layer['density']:.6f} {fuel} {layer['fraction']:.6e}"
        )
    return lines


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        summary = json.loads(Path(args[1]).read_text())
        if summary.get("format") != "xlc-compare-p6":
            raise SystemExit(f"{args[1]}: not xlc-compare-p6")
        print("\n".join(report(summary)))
    elif len(args) == 9 and args[0] == "--json":
        Path(args[1]).write_text(json.dumps(compare(args[2:]), indent=1) + "\n")
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
