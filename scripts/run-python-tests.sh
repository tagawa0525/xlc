#!/usr/bin/env bash
# 追跡している Python のテスト（test_*.py と *_test.py）を、テストのあるディレクトリごとに unittest で実行する。
# run-checks.sh が scripts/run-python-tests.sh を見つけて呼ぶ（nixfiles の ADR-0036）。
#
# python-tests の devShell で実行する: examples/vera/test_openmc_model.py が OpenMC の Python を import し、default の
# devShell には無いため。その Python は default のパッケージに OpenMC を足した上位集合なので、ほかのテストも同じ shell で
# 通る（flake.nix の testsPython）。references の devShell は参照のステージが使うパッケージだけで、検証用のツールの
# pyyaml などが無い（ADR-0212）。テストを default と分けて skip しない（落ちるテストを隠すため）。
# discover でなく追跡しているモジュール名を渡す: 追跡していないファイルを実行せず、子のディレクトリのテストを親で
# 重ねて実行しないため（run-checks.sh の python_subdir_tests と同じ）。
set -euo pipefail

cd -- "$(git rev-parse --show-toplevel)"

# "<dir> <module>" の行。ルート直下のファイルは dir を . にする
mapfile -t entries < <(git ls-files -- ':(glob)**/test_*.py' ':(glob)**/*_test.py' | while read -r f; do
  if [[ "$f" == */* ]]; then echo "${f%/*} ${f##*/}"; else echo ". $f"; fi
done | sed 's/\.py$//' | sort -u)
if (( ${#entries[@]} == 0 )); then
  echo "ERROR: 追跡している test_*.py / *_test.py が無い（この script が要らない）" >&2
  exit 1
fi

# shellcheck disable=SC2016 # 変数は nix develop の中の bash が展開する
exec nix develop .#python-tests -c bash -c '
  set -euo pipefail
  dir=""; mods=()
  run() { echo "== $dir: python3 -m unittest ${mods[*]}"; (cd -- "$dir" && python3 -m unittest "${mods[@]}"); }
  for e in "$@"; do
    d=${e%% *}; m=${e#* }
    if [[ "$d" != "$dir" ]]; then
      [[ -n "$dir" ]] && run
      dir=$d; mods=()
    fi
    mods+=("$m")
  done
  run
' _ "${entries[@]}"
