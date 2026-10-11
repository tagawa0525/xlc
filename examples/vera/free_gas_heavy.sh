#!/usr/bin/env bash
# #0005 の、質量比が 1 より大きい核種の弾性散乱を 400 kT 未満で自由ガスとして扱ったときの xlc の k を、均質な無限媒質と
# 1B のピンセルで求める。xlc の変更は sab_residuals.patch の仮のもの（環境変数は sab_residuals.sh の説明と、
# XLC_TRIAL_FREE_GAS_HEAVY=1: S(α, β) を持たず質量比が 1 より大きい核種の自由ガスの上限を max(4 eV, 400 kT) にする）。
#
# 使い方（リポジトリの根で。全体のライブラリを 2 回作るので slurm-run.sh で流す）:
#   examples/vera/free_gas_heavy.sh <作業ディレクトリ>
#
# 環境変数 REFERENCES（既定はリポジトリの references/）に ENDF と OpenMC の venv（references/openmc/venv）があること。
# ライブラリは、H の扱いを sab_residuals.sh の (a)・(d) にしたもの（ad）、それに上の変更を足したもの（ade）、ad の O-16
# だけを ade のものにしたもの（ad-o16）。
set -euo pipefail

if (($# != 1)); then
  sed -n '2,11p' "$0" >&2
  exit 1
fi
repo=$(git rev-parse --show-toplevel)
work=$(realpath -m "$1")
references=${REFERENCES:-$repo/references}
python=$references/openmc/venv/bin/python
mkdir -p "$work"
cd "$work"

# パッチを当てた xlc を作る（今のコミットの別の worktree で）。パッチ、ツールと入力もその worktree から読む
src=$work/src
git -C "$repo" worktree add --detach "$src" HEAD
trap 'git -C "$repo" worktree remove --force "$src"' EXIT
vera=$src/examples/vera
git -C "$src" apply "$vera/sab_residuals.patch"
export CARGO_TARGET_DIR=$work/target
cargo build --release --manifest-path "$src/Cargo.toml" -p xlc-xs
cargo build --release --manifest-path "$src/Cargo.toml" -p xlc-lattice --bin lattice --example homogeneous
xs=$CARGO_TARGET_DIR/release/xs
lattice=$CARGO_TARGET_DIR/release/lattice
homogeneous=$CARGO_TARGET_DIR/release/examples/homogeneous

# 媒質
"$python" "$vera/homogeneous_medium.py" 1B h-u5-u8.toml H-1 U-235 U-238
"$python" "$vera/homogeneous_medium.py" 1B h-u5-u8-o.toml H-1 U-235 U-238 O-16
"$python" "$vera/homogeneous_medium.py" 1B no-clad.toml H-1 U-234 U-235 U-236 U-238 O-16 B-10 B-11
"$python" "$vera/homogeneous_medium.py" 1B all.toml

# 全体のライブラリ
neutrons=$references/endf/ENDF-B-VIII.0_neutrons
thermal=$references/endf/ENDF-B-VIII.0_thermal_scatt
for v in ad ade; do
  sed -e "s#^output = .*#output = \"vera-$v.msgpack\"#" \
    -e "s#^endf-directory = .*#endf-directory = \"$neutrons\"#" \
    -e "s#^thermal-directory = .*#thermal-directory = \"$thermal\"#" "$vera/xs.toml" >"vera-$v.toml"
done
h=(XLC_TRIAL_SAB_UPPER=9.877 XLC_TRIAL_SCT=1 XLC_TRIAL_FREE_GAS_UPPER=48.052)
env "${h[@]}" "$xs" vera-ad.toml
env "${h[@]}" XLC_TRIAL_FREE_GAS_HEAVY=1 "$xs" vera-ade.toml
"$python" - vera-ad.msgpack vera-ade.msgpack vera-ad-o16.msgpack <<'EOF'
import sys
from pathlib import Path
import msgpack
ad, ade, out = sys.argv[1:]
a = msgpack.unpackb(Path(ad).read_bytes())
b = msgpack.unpackb(Path(ade).read_bytes())
for key in ("group_bounds", "sigma0", "pointwise_range"):
    if a[key] != b[key]:
        raise SystemExit(f"{key} differs")
names = [n["name"] for n in a["nuclides"]]
if names != [n["name"] for n in b["nuclides"]]:
    raise SystemExit("nuclides differ")
# S(α, β) を使う H-1 のほかは、すべて変わるはず（変わらなければ切り替えが効いていない）
for x, y in zip(a["nuclides"], b["nuclides"]):
    print(f"{x['name']}: {'same' if x == y else 'changed'}")
    if (x["name"] == "H-1") != (x == y):
        raise SystemExit(f"{x['name']}: unexpected effect of XLC_TRIAL_FREE_GAS_HEAVY")
i = names.index("O-16")
a["nuclides"][i] = b["nuclides"][i]
Path(out).write_bytes(msgpack.packb(a))
EOF

# 均質な媒質
for v in ad ad-o16 ade; do
  for m in h-u5-u8 h-u5-u8-o no-clad all; do
    k=$("$homogeneous" "$m.toml" "vera-$v.msgpack" --balance "bal-$v-$m.msgpack" | grep balance)
    echo "$v $m $k"
  done
done | tee homogeneous.txt

# 1B のピンセル
for v in ad ad-o16 ade; do
  sed "s#^library = .*#library = \"$work/vera-$v.msgpack\"#" "$vera/p1/pincell-1B.toml" >"pin-$v.toml"
  k=$("$lattice" "pin-$v.toml" | grep '^k ')
  echo "pincell $v $k"
done | tee pincell.txt
