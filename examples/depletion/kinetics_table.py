"""燃焼度の軸を持つテーブルの動特性と冷却計算のデータ（ADR-0238）を、燃焼度の点ごとの公称の状態で要約する。

使い方: python3 kinetics_table.py --json <要約の出力先> <テーブル（xlc-table）>
        python3 kinetics_table.py --report <要約の JSON>

- 公称の状態は、連続の軸がどれも公称の値（`nominal`）の格子点で、制御棒の軸があれば `out` の点
- 中性子束は少数群定数の無限媒質の基本モード φ = A⁻¹ χ（A = diag(Σt) − Σsᵀ。Σs は [元の群][先の群]）。そのときの
  k = νΣfᵀ φ を `k`、テーブルの `k_inf`（格子計算の k∞）との差を `k_minus_lattice_pcm`（(1/k_inf − 1/k) × 1e5）に書く
- 遅発中性子の割合 β_k = Σ_G ν_dΣf,k,G φ_G / Σ_G νΣf_G φ_G（重みを付けない割合）と、族の崩壊定数
  λ_k = Σ_G ν_dΣf,k,G φ_G / Σ_G (ν_dΣf/λ)_k,G φ_G（core の `decay_constants` と同じ式）
- 群速度の逆数の中性子束の平均 Σ_G (1/v)_G φ_G / Σ_G φ_G (s/cm)
- 冷却計算の断面積は、核種の数と、核種ごとの 1 群の核分裂と捕獲（MT102 の経路の和）を `nuclides` の中から `--nuclide` で選んだ
  ものだけ書く（既定は U-235、Pu-239、Cs-133、Cs-134、Eu-154）
- `inputs` は読んだファイルのパス（実行した cwd から）と md5
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import msgpack
import numpy as np

NUCLIDES = ["U-235", "Pu-239", "Cs-133", "Cs-134", "Eu-154"]


def nominal_points(table: dict) -> list[int]:
    """燃焼度の点の順に、公称の状態の点の番号（`points` の中の番号）"""
    axes = table["axes"]
    names = [a["name"] for a in axes]
    if names[0] != "burnup":
        raise ValueError("the first axis must be the burnup")
    wanted = []
    for a in axes[1:]:
        if a["name"] == "rod":
            wanted.append(a["values"].index("out"))
            continue
        if a.get("nominal") is None:
            raise ValueError(f"the {a['name']} axis has no nominal value")
        wanted.append(a["points"].index(a["nominal"]))
    found = {}
    for i, p in enumerate(table["points"]):
        if list(p[1:]) == wanted:
            found[p[0]] = i
    burnups = len(axes[0]["points"])
    if sorted(found) != list(range(burnups)):
        raise ValueError("the table lacks the nominal state at some burnup points")
    return [found[b] for b in range(burnups)]


def fundamental(state: dict) -> tuple[np.ndarray, float]:
    """無限媒質の基本モードの中性子束（和を 1 にしたもの）と k"""
    total = np.asarray(state["total"], dtype=float)
    scatter = np.asarray(state["scatter"], dtype=float)
    chi = np.asarray(state["chi"], dtype=float)
    nu = np.asarray(state["nu_fission"], dtype=float)
    a = np.diag(total) - scatter.T
    phi = np.linalg.solve(a, chi)
    k = float(nu @ phi)
    return phi / phi.sum(), k


def summary(table: dict, nuclides: list[str]) -> dict:
    burnups = table["axes"][0]["points"]
    cooling = table["cooling"]
    kinetics = table["kinetics"]["states"]
    picked = [n for n in nuclides if n in cooling["nuclides"]]
    rows = []
    for b, p in enumerate(nominal_points(table)):
        state = table["states"][p]
        phi, k = fundamental(state)
        kin = kinetics[p]
        production = np.asarray(state["nu_fission"], dtype=float) @ phi
        delayed = np.asarray(kin["delayed_nu_fission"], dtype=float) @ phi
        precursors = np.asarray(kin["delayed_precursors"], dtype=float) @ phi
        beta = delayed / production
        one_group = {}
        for name in picked:
            n = cooling["points"][b][cooling["nuclides"].index(name)]
            capture = sum(sum(r["channels"]) for r in n["reactions"] if r["mt"] == 102)
            one_group[name] = {"fission": n["fission"], "capture": capture}
        rows.append(
            {
                "burnup": burnups[b],
                "k_inf": state["k_inf"],
                "k": k,
                "k_minus_lattice_pcm": (1.0 / state["k_inf"] - 1.0 / k) * 1e5,
                "beta_pcm": [float(x) * 1e5 for x in beta],
                "beta_total_pcm": float(beta.sum()) * 1e5,
                "decay_constants": [float(x) for x in delayed / precursors],
                "inverse_velocity": float(
                    np.asarray(kin["inverse_velocity"], dtype=float) @ phi
                ),
                "cooling": one_group,
            }
        )
    return {"cooling_nuclides": len(cooling["nuclides"]), "points": rows}


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def report(data: dict) -> str:
    lines = [f"cooling nuclides: {data['cooling_nuclides']}"]
    lines.append(
        "burnup  k−lattice(pcm)  β(pcm)   λ(1/s) by family                       <1/v>(s/cm)"
    )
    for r in data["points"]:
        lam = " ".join(f"{x:.4f}" for x in r["decay_constants"])
        lines.append(
            f"{r['burnup']:6g}  {r['k_minus_lattice_pcm']:+.4f}  {r['beta_total_pcm']:7.2f}  {lam}  "
            f"{r['inverse_velocity']:.4e}"
        )
    return "\n".join(lines)


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--json", type=Path)
    group.add_argument("--report", type=Path)
    parser.add_argument("table", type=Path, nargs="?")
    parser.add_argument("--nuclide", action="append")
    args = parser.parse_args(argv)
    if args.report:
        print(report(json.loads(args.report.read_text())))
        return 0
    if args.table is None:
        parser.error("give the table")
    table = msgpack.unpackb(args.table.read_bytes(), raw=False)
    data = summary(table, args.nuclide or NUCLIDES)
    data["inputs"] = {str(args.table): md5(args.table)}
    args.json.write_text(json.dumps(data, indent=1, ensure_ascii=False) + "\n")
    print(report(data))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
