#!/usr/bin/env bash
# #0005 の、H in H2O の S(α, β) の上限と、その上と表の外の扱いを変えた xlc の k を、均質な無限媒質と 1B のピンセルで
# 求め、OpenMC と比較する。xlc の変更は仮のもので、sab_residuals.patch として持つ（本番のソースには入れない）。
#
# 使い方（リポジトリの根で。ライブラリの作り直しを含み、全体で約 20 分なので slurm-run.sh で流す）:
#   examples/vera/sab_residuals.sh <作業ディレクトリ> <OpenMC の均質な媒質のディレクトリ>
#
# <OpenMC の均質な媒質のディレクトリ> は、`openmc_homogeneous.py` の出力を bal-h-u5/、bal-h-u5-u8/ に持つこと
# （verification.md「反応率の比較」の H-1・U-235 と H-1・U-235・U-238）。環境変数 REFERENCES（既定は
# リポジトリの references/）に ENDF と OpenMC の venv（references/openmc/venv）があること。
# examples/vera/vera.msgpack（ADR-0121 で追跡しない）を、H-1 を差し替える元に使う。
#
# パッチの仮の変更（環境変数で切り替える）:
#   XLC_TRIAL_SAB_UPPER=<eV>       S(α, β) を使う入射エネルギーの上限（既定は 4 eV。上方散乱の表の覆いは確認しない）
#   XLC_TRIAL_SCT=1                表の外（α、β の範囲の外と、S = 0 の格子点に接する区間）を NJOY THERMR の短時間衝突
#                                  近似（有効温度は MF7）で補う（E' は β の表の範囲で打ち切るので、β の外には届かない）
#   XLC_TRIAL_UFG_SAB=1            ポイントワイズの断面積（UFG）の H の上限より下の弾性散乱を S(α, β) の散乱断面積にする
#   XLC_TRIAL_FREE_GAS_UPPER=<eV>  S(α, β) の上限からこの入射エネルギーまでの H の弾性散乱の移行を自由ガスの核にする
#   XLC_TRIAL_OUT_WEIGHT=<核種>:<最初の群>-<最後の群>,...（`homogeneous` の例だけ）
#                                  UFG の範囲の群の P0 の行の、群から出る部分に、静止標的・重心系で等方とした群から
#                                  出る割合の、UFG の中性子束の重みと 1/E の重みの比を掛ける（行の和は保つ）
set -euo pipefail

if (($# != 2)); then
  sed -n '2,21p' "$0" >&2
  exit 1
fi
repo=$(git rev-parse --show-toplevel)
work=$(realpath -m "$1")
openmc=$(realpath "$2")
references=${REFERENCES:-$repo/references}
python=$references/openmc/venv/bin/python
upper=9.877
mkdir -p "$work"
cd "$work"

# パッチを当てた xlc を作る（今のコミットの別の worktree で）。パッチ、ツールと入力もその worktree から読み、結果が
# 今のコミットから出たものにする（追跡しない vera.msgpack だけは $repo から）
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
"$python" "$vera/homogeneous_medium.py" 1B h-u5.toml H-1 U-235
"$python" "$vera/homogeneous_medium.py" 1B h-u5-u8.toml H-1 U-235 U-238
"$python" "$vera/homogeneous_medium.py" 1B h-u5-u8-o.toml H-1 U-235 U-238 O-16
"$python" "$vera/homogeneous_medium.py" 1B no-clad.toml H-1 U-234 U-235 U-236 U-238 O-16 B-10 B-11
"$python" "$vera/homogeneous_medium.py" 1B all.toml

# ライブラリ: 均質な媒質の核種だけのもの（<変種>.msgpack）と、1B のピンセルに H-1 を差し替える、U、O と制御棒の
# 核種を除いたもの（light-<変種>.msgpack。ポイントワイズの断面積の範囲は全体のライブラリと同じ）
neutrons=$references/endf/ENDF-B-VIII.0_neutrons
thermal=$references/endf/ENDF-B-VIII.0_thermal_scatt
for v in base t0 a b ab ad sct4; do
  cat >"$v.toml" <<EOF
output = "$v.msgpack"
endf-directory = "$neutrons"
thermal-directory = "$thermal"
sigma0 = [1e10, 1e5, 1e4, 1e3, 300, 100, 30, 10, 1]

[nuclide.H-1]
endf = "n-001_H_001.endf"
temperatures = [600]
thermal = "tsl-HinH2O.endf"

[nuclide.U-235]
endf = "n-092_U_235.endf"
temperatures = [600]

[nuclide.U-238]
endf = "n-092_U_238.endf"
temperatures = [600]
EOF
done
for v in base t0 a ad; do
  "$python" - "$vera/xs.toml" "light-$v.toml" "light-$v.msgpack" "$neutrons" "$thermal" <<'EOF'
import re, sys
source, out, library, neutrons, thermal = sys.argv[1:]
text = open(source).read()
head, *blocks = re.split(r"(?m)^(?=\[nuclide\.)", text)
head = re.sub(r"(?m)^output = .*$", f'output = "{library}"', head)
head = re.sub(r"(?m)^endf-directory = .*$", f'endf-directory = "{neutrons}"', head)
head = re.sub(r"(?m)^thermal-directory = .*$", f'thermal-directory = "{thermal}"', head)
keep = [b for b in blocks if not re.match(r"\[nuclide\.(U-|O-16|Ag|In|Cd)", b)]
open(out, "w").write(head + "".join(keep))
EOF
done
library() {
  local v=$1 toml=$2
  case $v in
  base) "$xs" "$toml" ;;
  t0) XLC_TRIAL_SAB_UPPER=$upper "$xs" "$toml" ;;
  a) XLC_TRIAL_SAB_UPPER=$upper XLC_TRIAL_SCT=1 "$xs" "$toml" ;;
  b) XLC_TRIAL_SAB_UPPER=$upper XLC_TRIAL_UFG_SAB=1 "$xs" "$toml" ;;
  ab) XLC_TRIAL_SAB_UPPER=$upper XLC_TRIAL_SCT=1 XLC_TRIAL_UFG_SAB=1 "$xs" "$toml" ;;
  # 自由ガスの核の P1 以上は、E/kT が約 1400 を超えると e^(-β/2) が桁あふれするので、群の境界 48.052 eV までにする
  ad) XLC_TRIAL_SAB_UPPER=$upper XLC_TRIAL_SCT=1 XLC_TRIAL_FREE_GAS_UPPER=48.052 "$xs" "$toml" ;;
  sct4) XLC_TRIAL_SCT=1 "$xs" "$toml" ;;
  *)
    echo "unknown variant: $v" >&2
    exit 1
    ;;
  esac
}
for v in base t0 a b ab ad sct4; do library "$v" "$v.toml"; done
for v in base t0 a ad; do
  library "$v" "light-$v.toml"
  "$python" - "$repo/examples/vera/vera.msgpack" "light-$v.msgpack" "vera-$v.msgpack" <<'EOF'
import sys
from pathlib import Path
import msgpack
full, source, out = sys.argv[1:]
a = msgpack.unpackb(Path(full).read_bytes())
b = msgpack.unpackb(Path(source).read_bytes())
for key in ("group_bounds", "sigma0", "pointwise_range"):
    if a[key] != b[key]:
        raise SystemExit(f"{key} differs")
(h,) = [n for n in b["nuclides"] if n["name"] == "H-1"]
(i,) = [k for k, n in enumerate(a["nuclides"]) if n["name"] == "H-1"]
same = a["nuclides"][i] == h
print(f"vera-H-1 {out}: same as vera.msgpack {same}")
# 仮の変更のない base で H-1 が違えば、vera.msgpack はこのコミットのものでない（ほかの核種を混ぜられない）
if out == "vera-base.msgpack" and not same:
    raise SystemExit("vera.msgpack differs from this commit in H-1; rebuild it")
a["nuclides"][i] = h
Path(out).write_bytes(msgpack.packb(a))
EOF
done

# 均質な媒質
{
  for v in base t0 a b ab ad sct4; do
    for m in h-u5 h-u5-u8; do
      "$homogeneous" "$m.toml" "$v.msgpack" --balance "x-$v-$m.msgpack" | grep balance | sed "s/^/$v $m/"
      "$python" "$vera/compare_groups.py" "x-$v-$m.msgpack" "$openmc/bal-$m/balance.msgpack" \
        "$v.msgpack" >"g-$v-$m.txt"
    done
  done
  for v in a ad; do
    for spec in H-1:14-26 H-1:14-27,U-235:14-27 H-1:14-27,U-235:14-27,U-238:14-27 U-238:14-27; do
      for m in h-u5 h-u5-u8; do
        if [[ $m == h-u5 && $spec == *U-238* ]]; then continue; fi
        tag=${spec//[:,]/_}
        XLC_TRIAL_OUT_WEIGHT=$spec "$homogeneous" "$m.toml" "$v.msgpack" --balance "x-$v-c-$tag-$m.msgpack" |
          grep balance | sed "s/^/$v+c($spec) $m/"
        "$python" "$vera/compare_groups.py" "x-$v-c-$tag-$m.msgpack" \
          "$openmc/bal-$m/balance.msgpack" "$v.msgpack" >"g-$v-c-$tag-$m.txt"
      done
    done
  done
  # 核種を足した媒質（OpenMC は k だけ。verification.md「均質な無限媒質」）は、H-1 を差し替えた全体のライブラリで
  for v in base t0 a ad; do
    for m in h-u5-u8-o no-clad all; do
      "$homogeneous" "$m.toml" "vera-$v.msgpack" | grep 'UFG  ' | sed "s/^/vera-$v $m/"
    done
  done
} | tee homogeneous.txt

# 1B のピンセル
for v in base t0 a ad; do
  sed "s#^library = .*#library = \"$work/vera-$v.msgpack\"#" "$vera/p1/pincell-1B.toml" >"pin-$v.toml"
  k=$("$lattice" "pin-$v.toml" | grep '^k ')
  echo "pincell $v $k"
done | tee pincell.txt
