#!/usr/bin/env bash
# check-no-inline-tests.sh — 本番のファイル（crates/*/src の tests.rs 以外）にテストのコードがないことを確かめる
# （.pre-commit-config.yaml の no-inline-tests と、crates/docs/tests/consistency.rs が呼ぶ。ADR-0143、#0089）
#
# Usage: check-no-inline-tests.sh [リポジトリの根]
#
# 見つけるもの: `#[test]`（`#[tokio::test]`、`#[bench]` などの形を含む）、`mod tests { … }`、条件に test を含む cfg / cfg_attr / cfg!、`#[should_panic]`。コメントの行は見ない。
# 許すもの（ADR-0143 の「決定と理由」）:
#   - `#[cfg(test)]` の次の行にある `mod tests;`、`pub(crate) mod tests;`、`#[path = "…/tests.rs"]`（main.rs。続く行が
#     `mod tests;`）、パスが tests を通る `pub(crate) use …;`
#   - クレートの入口の `#![cfg_attr(test, allow(clippy::disallowed_macros))]`（lint の設定。テストのコードではない）

set -euo pipefail

root=${1:-$(git rev-parse --show-toplevel)}
cd -- "$root"

shopt -s globstar nullglob
files=()
for f in crates/*/src/**/*.rs; do
  [[ $(basename -- "$f") == tests.rs ]] || files+=("$f")
done
if (( ${#files[@]} == 0 )); then
  echo "ERROR: crates/*/src に .rs が無い（$root はリポジトリの根か）" >&2
  exit 2
fi

awk '
  function flag(why) { printf "%s:%d: %s: %s\n", FILENAME, FNR, why, $0; bad = 1 }
  FNR == 1 { if (prev != "") { printf "%s: ファイルの末尾が #[cfg(test)] または #[path] で終わる\n", last; bad = 1 } prev = "" }
  { last = FILENAME }
  {
    line = $0
    allowed = 0
    if (prev == "cfg_test") {
      if (line ~ /^(pub\(crate\) )?mod tests;$/) { allowed = 1; prev = "" }
      else if (line ~ /^#\[path = "[a-z_\/]*tests\.rs"\]$/) { allowed = 1; prev = "path" ; next }
      else if (line ~ /^pub\(crate\) use [A-Za-z0-9_:]*tests[A-Za-z0-9_:]*( as [A-Za-z0-9_]+)?;$/) { allowed = 1; prev = "" }
      else { flag("#[cfg(test)] の次が mod tests; などでない"); prev = "" }
    } else if (prev == "path") {
      if (line ~ /^(pub\(crate\) )?mod tests;$/) allowed = 1
      else flag("#[path] の次が mod tests; でない")
      prev = ""
    }
    if (allowed) next
    if (line ~ /^[[:space:]]*\/\//) next # コメントの行は見ない
    if (line ~ /^#\[cfg\(test\)\]$/) { prev = "cfg_test"; next }
    if (line ~ /^#!\[cfg_attr\(test, allow\(clippy::disallowed_macros\)\)\]$/) next
    if (line ~ /#\[([a-z_]+::)*(test|bench|test_case|rstest)(\]|\()|#\[should_panic|mod tests?[[:space:]]*\{/) flag("テストのコード")
    else if (line ~ /(cfg|cfg_attr)[(!][^)]*test/ || line ~ /cfg!\(/) flag("test に関わる cfg")
  }
  END { if (prev != "") { printf "%s: ファイルの末尾が #[cfg(test)] または #[path] で終わる\n", last; bad = 1 } exit bad }
' "${files[@]}" || {
  echo "FAILED: 本番のファイルにテストのコードがある。<モジュール>/tests.rs に移す（ADR-0143）" >&2
  exit 1
}
