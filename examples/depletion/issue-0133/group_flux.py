"""#0133 の試行: 1C のピンセルの燃料と減速材の群 1〜13 の中性子束（生成 1 あたり）の xlc/OpenMC − 1 と k の差を出す。

使い方: python3 group_flux.py <作業ディレクトリ>（`openmc_pin_rates.py` の openmc-rates.msgpack、`balance` の例の xlc-balance.msgpack）
"""

import sys
from pathlib import Path

import msgpack
import numpy as np

FAST = 13


def main():
    work = Path(sys.argv[1])
    omc = msgpack.unpackb((work / "openmc-rates.msgpack").read_bytes())
    xlc = msgpack.unpackb((work / "xlc-balance.msgpack").read_bytes())
    produced_o = sum(np.sum(m["nu-fission"]) for m in omc["materials"].values())
    produced_x = sum(np.sum(m["production"]) for m in xlc["materials"].values())
    dk = (xlc["k"] - omc["k"]) * 1e5
    print(
        f"k xlc {xlc['k']:.5f} OpenMC {omc['k']:.5f} ± {omc['k_std']:.5f} ({dk:+.0f} pcm)"
    )
    for name_o, name_x in [("fuel", "fuel"), ("moderator", "water")]:
        m = next(v for k, v in omc["materials"].items() if k.startswith(name_o))
        flux_o = np.array(m["flux"]) / produced_o
        flux_x = np.array(xlc["materials"][name_x]["flux"]) / produced_x
        ratios = " ".join(
            f"{(flux_x[g] / flux_o[g] - 1) * 100:+.2f}" for g in range(FAST)
        )
        print(f"{name_x:6} {ratios}")


main()
