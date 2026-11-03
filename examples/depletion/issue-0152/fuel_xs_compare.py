"""#0152 の切り分け: 燃料の核種の実効の吸収断面積（捕獲 + 核分裂）を、xlc（fuel_xs の例の出力）と OpenMC（openmc_pin_transfer.py の
出力の燃料の反応率 / 燃料の中性子束 / 数密度）で群ごとに比べる。

使い方（references の devShell で）: python3 fuel_xs_compare.py <OpenMC の transfer.msgpack> <xlc の fuel-xs.msgpack> <核種,...> <群,...>
群は 1 始まりで高いエネルギーから数える。相対標準偏差は OpenMC の捕獲と核分裂の反応率のもの（中性子束の分は入れない）。
"""

import math
import sys

import msgpack


def read(path):
    with open(path, "rb") as f:
        return msgpack.unpack(f, raw=False, strict_map_key=False)


def main():
    openmc = read(sys.argv[1])
    xlc = read(sys.argv[2])["nuclides"]
    rates, flux, densities = openmc["rates"], openmc["flux"], openmc["densities"]
    for nuclide in sys.argv[3].split(","):
        name = nuclide.replace("-", "")
        for g in (int(v) - 1 for v in sys.argv[4].split(",")):
            capture, fission = rates[nuclide]["(n,gamma)"], rates[nuclide]["fission"]
            rate = capture["mean"][g] + fission["mean"][g]
            std = math.hypot(capture["std_dev"][g], fission["std_dev"][g])
            reference = rate / flux[g] / densities[nuclide]
            value = xlc[name]["capture"][g] + xlc[name]["fission"][g]
            print(
                f"{nuclide} group {g + 1}: xlc {value:.5g} b, OpenMC {reference:.5g} b "
                f"(rate rel. std {std / rate:.2%}), xlc / OpenMC - 1 {value / reference - 1:+.2%}"
            )


if __name__ == "__main__":
    main()
