"""xlc の全核種の鎖（ADR-0164）と、OpenMC の `Chain.from_endf` で同じ ENDF から作った鎖を比べる。

使い方（references の devShell で）:

    python3 examples/depletion/compare_chain.py examples/depletion/chain.msgpack examples/depletion/openmc-chain.xml

核種の集合、半減期、崩壊の分岐比と娘、反応の生成物、核分裂の収率を比べ、違いを種類ごとに数えて書く。違いの理由は
examples/depletion/verification.md の「燃焼鎖」。
"""

import re
import sys
from collections import Counter, defaultdict

import msgpack
from openmc.deplete import Chain

# xlc の反応の MT と OpenMC の反応の名前
REACTIONS = {
    16: "(n,2n)",
    17: "(n,3n)",
    37: "(n,4n)",
    102: "(n,gamma)",
    103: "(n,p)",
    107: "(n,a)",
}
# 反応の放出する軽い核種
LIGHT_OF = {103: "H1", 107: "He4"}


def openmc_name(name: str) -> str:
    """ADR-0015 の名前（U-235、Am-242m、Hf-178m2）を OpenMC の名前（U235、Am242_m1、Hf178_m2）にする"""
    m = re.fullmatch(r"([A-Z][a-z]?)-(\d+)(m(\d*))?", name)
    if m is None:
        raise ValueError(f"not a nuclide name: {name}")
    symbol, mass, isomer, index = m.groups()
    if isomer is None:
        return f"{symbol}{mass}"
    return f"{symbol}{mass}_m{index or 1}"


def read_xlc(path):
    with open(path, "rb") as f:
        data = msgpack.unpackb(f.read(), raw=False)
    if data["format"] != "xlc-chain":
        sys.exit(f"{path}: not an xlc chain")
    nuclides = data["nuclides"]
    names = [openmc_name(n["name"]) for n in nuclides]
    return {name: n for name, n in zip(names, nuclides)}, names


def xlc_daughter(mode, names):
    """xlc の崩壊の型の娘（軽い核種の放出を除く）。娘がなければ None"""
    alphas = mode["mode"].split(",").count("alpha")
    protons = mode["mode"].split(",").count("p")
    rest = []
    for p in mode["products"]:
        name = names[p["nuclide"]]
        count = p["count"]
        # 軽い核種の放出の分（α の He-4、陽子の H-1）を除く
        if name == "He4" and alphas:
            take = min(count, alphas)
            alphas -= take
            count -= take
        elif name == "H1" and protons:
            take = min(count, protons)
            protons -= take
            count -= take
        if count:
            rest.append(name)
    if len(rest) > 1:
        raise ValueError(f"several daughters {rest} of {mode['mode']}")
    return rest[0] if rest else None


def main():
    xlc_path, openmc_path = sys.argv[1], sys.argv[2]
    xlc, names = read_xlc(xlc_path)
    omc = Chain.from_xml(openmc_path)
    omc_names = {n.name for n in omc.nuclides}

    print("## 核種の集合")
    only_xlc = sorted(set(xlc) - omc_names)
    only_omc = sorted(omc_names - set(xlc))
    print(
        f"xlc {len(xlc)}、OpenMC {len(omc_names)}、xlc だけ {only_xlc}、OpenMC だけ {only_omc}"
    )

    print("\n## 半減期")
    stable_mismatch = []
    rel = []
    for name in sorted(set(xlc) & omc_names):
        x = xlc[name]["decay"]
        o = omc[name].half_life
        if (x is None) != (o is None):
            stable_mismatch.append(name)
            continue
        if x is not None:
            rel.append((abs(x["half_life"] - o) / o, name))
    rel.sort(reverse=True)
    if not rel:
        sys.exit("no unstable nuclide in both chains")
    print(f"安定と不安定の食い違い {len(stable_mismatch)} {stable_mismatch[:10]}")
    print(
        f"不安定 {len(rel)} 核種の半減期の相対差の最大 {rel[0][0]:.3e}（{rel[0][1]}）"
    )

    print("\n## 崩壊の分岐")
    kinds = Counter()
    examples = defaultdict(list)
    ratio_diffs = []
    for name in sorted(set(xlc) & omc_names):
        x = xlc[name]["decay"]
        if x is None or omc[name].half_life is None:
            continue
        xm = {}
        # 同じ型と娘が 2 度書かれた評価がある（Es-246 の ec/beta+）ので足す
        for m in x["modes"]:
            k = (m["mode"], xlc_daughter(m, names))
            xm[k] = xm.get(k, 0.0) + m["ratio"]
        om = {}
        for m in omc[name].decay_modes:
            om[(m.type, m.target)] = om.get((m.type, m.target), 0.0) + m.branching_ratio
        for (mode, daughter), ratio in xm.items():
            if (mode, daughter) in om:
                ratio_diffs.append((abs(ratio - om[(mode, daughter)]), name, mode))
                continue
            # 同じ型で娘が違う
            other = [d for (t, d) in om if t == mode]
            if "sf" in mode.split(","):
                kind = "SF: xlc は親を減らすだけ、OpenMC は娘を親そのものにする"
            elif daughter is None:
                kind = (
                    "娘が decay にない: xlc は娘なし、OpenMC は replace_missing の核種"
                )
            else:
                kind = "その他"
            kinds[kind] += 1
            examples[kind].append(f"{name} {mode} xlc={daughter} OpenMC={other}")
        for key in om:
            if not any(k[0] == key[0] for k in xm):
                kinds["OpenMC だけの型"] += 1
                examples["OpenMC だけの型"].append(f"{name} {key}")
    ratio_diffs.sort(reverse=True)
    if not ratio_diffs:
        sys.exit("no decay branch in both chains")
    print(
        f"型と娘が同じ分岐 {len(ratio_diffs)}、分岐比の差の最大 {ratio_diffs[0][0]:.3e}（{ratio_diffs[0][1]} {ratio_diffs[0][2]}）"
    )
    print(f"差が 1e-12 を超える分岐 {sum(1 for d in ratio_diffs if d[0] > 1e-12)}")
    for kind, n in kinds.most_common():
        print(f"{kind}: {n}（例: {'; '.join(examples[kind][:4])}）")

    print("\n## 反応")
    kinds = Counter()
    examples = defaultdict(list)
    same = 0
    for name in sorted(set(xlc) & omc_names):
        xr = {}
        for r in xlc[name]["reactions"]:
            targets = []
            for c in r["channels"]:
                # 生成物から、反応の放出する軽い核種（(n,p) の H-1、(n,α) の He-4）を 1 つ除いた残りの核
                light = LIGHT_OF.get(r["mt"])
                heavy = []
                for p in c["products"]:
                    product = names[p["nuclide"]]
                    count = p["count"]
                    if product == light:
                        count -= 1
                        light = None
                    heavy += [product] * count
                targets.append(heavy[0] if heavy else None)
            xr[REACTIONS[r["mt"]]] = targets
        orx = {}
        for r in omc[name].reactions:
            if r.type == "fission":
                continue
            orx.setdefault(r.type, []).append(r.target)
        for rx in sorted(set(xr) | set(orx)):
            if rx not in orx or rx not in xr:
                kinds["片方だけの反応"] += 1
                examples["片方だけの反応"].append(f"{name} {rx}")
            elif xr[rx] == orx[rx]:
                same += 1
            elif len(xr[rx]) > 1:
                kinds["xlc が MF8 の状態に分ける"] += 1
                examples["xlc が MF8 の状態に分ける"].append(
                    f"{name} {rx} xlc={xr[rx]} OpenMC={orx[rx]}"
                )
            elif xr[rx] == [None]:
                kinds["生成物が decay にない: xlc は残りの核なし"] += 1
                examples["生成物が decay にない: xlc は残りの核なし"].append(
                    f"{name} {rx} OpenMC={orx[rx]}"
                )
            else:
                kinds["その他"] += 1
                examples["その他"].append(f"{name} {rx} xlc={xr[rx]} OpenMC={orx[rx]}")
        xf = xlc[name]["fission"] is not None
        of = any(r.type == "fission" for r in omc[name].reactions)
        if xf != of:
            kinds["核分裂の有無"] += 1
            examples["核分裂の有無"].append(name)
    print(f"生成物が同じ反応 {same}")
    for kind, n in kinds.most_common():
        print(f"{kind}: {n}（例: {'; '.join(examples[kind][:6])}）")

    print("\n## 核分裂の収率")
    rows = []
    sources = []
    for name in sorted(set(xlc) & omc_names):
        f = xlc[name]["fission"]
        if f is None:
            continue
        source = openmc_name(f["yield_source"])
        omc_source = getattr(omc[name], "_fpy", None) or name
        if source != name or omc_source != name:
            sources.append((name, source, omc_source))
            continue
        yd = omc[name].yield_data
        if [float(e) for e in yd.energies] != f["energies"]:
            rows.append((name, "energies", f["energies"], list(yd.energies)))
            continue
        for e, table in zip(f["energies"], f["yields"]):
            xy = {names[y["nuclide"]]: y["value"] for y in table}
            oy = dict(yd[e].items())
            keys = set(xy) | set(oy)
            diff = max(abs(xy.get(k, 0.0) - oy.get(k, 0.0)) for k in keys)
            worst = max(keys, key=lambda k: abs(xy.get(k, 0.0) - oy.get(k, 0.0)))
            # OpenMC は収率 0 の生成物も持つので、0 でないものを数える
            nonzero = sum(1 for v in oy.values() if v != 0.0)
            rows.append(
                (
                    name,
                    e,
                    len(xy),
                    nonzero,
                    sum(xy.values()),
                    sum(oy.values()),
                    diff,
                    worst,
                )
            )
    print(
        "| 核種 | E (eV) | 生成物 xlc | 0 でない生成物 OpenMC | 和 xlc | 和 OpenMC | 収率の差の最大 | その生成物 |"
    )
    print("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for r in rows:
        if r[1] == "energies":
            print(f"| {r[0]} | エネルギーが違う {r[2]} {r[3]} | | | | | | |")
        else:
            print(
                f"| {r[0]} | {r[1]:g} | {r[2]} | {r[3]} | {r[4]:.6f} | {r[5]:.6f} | {r[6]:.2e} | {r[7]} |"
            )
    same_source = sum(1 for _, s, o in sources if s == o)
    print(f"\n収率を借りる核種 {len(sources)}（借りる先が同じ {same_source}）")
    for name, s, o in sources:
        if s != o:
            print(f"  {name}: xlc {s}、OpenMC {o}")


if __name__ == "__main__":
    main()
