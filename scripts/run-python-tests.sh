#!/usr/bin/env bash
# 追跡している Python のテストを、テストのあるディレクトリごとに unittest で走らせる。
# run-checks.sh が scripts/run-python-tests.sh を見つけて呼ぶ（nixfiles の ADR-0036）。
#
# references の devShell で走らせる: examples/vera/test_openmc_model.py が OpenMC の Python を import し、default の
# devShell には無いため。references の Python は default のパッケージに OpenMC を足した上位集合なので、ほかのテストも
# 同じ shell で通る（flake.nix の referencesPython）。テストを default と分けて skip しない（落ちるテストを隠すため）。
set -euo pipefail

cd -- "$(git rev-parse --show-toplevel)"

mapfile -t dirs < <(git ls-files -- ':(glob)**/test_*.py' | sed 's|/[^/]*$||' | sort -u)
if (( ${#dirs[@]} == 0 )); then
  echo "ERROR: 追跡している test_*.py が無い（この script が要らない）" >&2
  exit 1
fi

# shellcheck disable=SC2016 # $d は nix develop の中の bash が展開する
exec nix develop .#references -c bash -c '
  set -euo pipefail
  for d in "$@"; do
    echo "== $d"
    (cd -- "$d" && python3 -m unittest discover -s . -p "test_*.py")
  done
' _ "${dirs[@]}"
