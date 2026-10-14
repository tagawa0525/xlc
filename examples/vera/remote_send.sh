#!/usr/bin/env bash
# 外部の計算機（VM）に、Problem 3 と 5-2D の OpenMC の参照解を流すための道具とデータを送る（docs/environment.md の
# 「外部の計算機での参照解」）。リポジトリの直下から走らせる。
#
# 使い方: examples/vera/remote_send.sh [--dry-run] [--allow-ip] <ssh の接続先> <VM の作業ディレクトリ>
#
# --allow-ip は、送る前に aws_allow_ip.sh でセキュリティグループを今の IP にする（AWS の VM のとき。--dry-run では呼ばない）。
# --dry-run は何も送らず、送る閉包のパス（件数と大きさ）と、rsync --dry-run の一覧を出す。
# 送るもの: OpenMC と Python API の閉包（flake の `.#openmc`、`.#openmc-python`。`nix copy --to ssh://`。VM で再ビルドしない）、
# ライブラリ（examples/vera/openmc-library）、入力（p3、p5）、スクリプト、VM で使う道具の場所を書いた tools.env。
set -euo pipefail

dry=0
allow=0
while [[ ${1:-} == --dry-run || ${1:-} == --allow-ip ]]; do
  if [[ $1 == --dry-run ]]; then dry=1; else allow=1; fi
  shift
done
if [[ $# -ne 2 || ${1:-} == --help || ${1:-} == -h ]]; then
  sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'
  exit 2
fi
target=$1
remote=$2

[[ -f flake.nix && -d examples/vera ]] || { echo "リポジトリの直下から走らせてください" >&2; exit 1; }
[[ -f examples/vera/openmc-library/cross_sections.xml ]] ||
  { echo "examples/vera/openmc-library がありません（dvc checkout examples/vera/dvc.yaml）" >&2; exit 1; }

if ((allow && !dry)); then "$(dirname "$0")/../aws/aws_allow_ip.sh"; fi

mapfile -t outs < <(nix build .#openmc .#openmc-python --no-link --print-out-paths)
openmc=${outs[0]}
python_env=${outs[1]}

echo "== 閉包"
mapfile -t closure < <(nix path-info -r "${outs[@]}")
echo "${#closure[@]} 個のパス、$(nix path-info -rs "${outs[@]}" | awk '{s+=$2} END {printf "%.0f MB", s/1e6}')"
printf '%s\n' "${closure[@]}" | head -5
echo "..."

files=(
  examples/openmc_common.py
  examples/openmc_input.py
  examples/openmc_run.py
  examples/openmc_resume.py
  examples/vera/openmc_model.py
  examples/vera/openmc_model_3d.py
  examples/vera/openmc_core_2d.py
  examples/vera/openmc_tallies.py
  examples/vera/run_remote.sh
  examples/vera/p3
  examples/vera/p5
  examples/vera/openmc-library
)
rsync_flags=(-a --relative --partial --info=stats1)
if ((dry)); then rsync_flags+=(--dry-run --itemize-changes); fi

echo "== rsync"
if ((dry)); then
  scratch=$(mktemp -d)
  rsync "${rsync_flags[@]}" "${files[@]}" "$scratch/" | grep -v '/$' | sed -n 1,25p
  rm -rf "$scratch"
else
  nix copy --to "ssh://$target" "${outs[@]}"
  rsync "${rsync_flags[@]}" "${files[@]}" "$target:$remote/"
fi
cat > /tmp/tools.env.$$ <<EOF
OPENMC_BIN=$openmc/bin
PYTHON=$python_env/bin/python3
EOF
cat /tmp/tools.env.$$
if ((! dry)); then rsync -a /tmp/tools.env.$$ "$target:$remote/tools.env"; fi
rm -f /tmp/tools.env.$$
