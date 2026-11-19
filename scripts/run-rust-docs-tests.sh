#!/usr/bin/env bash
# 文書だけの差分で run-checks.sh が実行する Rust のテスト（nixfiles の ADR-0057）。
# crates/docs のテストが文書どうしの整合と、検証の記録の表の `verification_tables.py --check` を実行する。
# design-records check は run-checks.sh が別の段で実行するので含めない。
#
# 除外と仮の前提の検査（#0294）もここで実行する: pre-commit の hook `exclusions` は rebase とマージでは走らないので、
# issue を閉じたブランチと、その issue を引く除外を足したブランチが別々に通ると、マージの後の main に閉じた issue を
# 指す印が残る。run-checks.sh はこのスクリプトをどの差分でも実行する（文書だけの差分ではこれだけ）
set -euo pipefail

cd -- "$(git rev-parse --show-toplevel)"

scripts/check-exclusions.py
# 形式を変える修正の検査（ADR-0258）も、マージの前に必ず通すためにここで実行する。実行ファイルを変えうる差分かは、
# 検査が自分で決める（ほかの差分では読まずに終わる）
scripts/check-readable-outputs.py
exec cargo test -p xlc-docs
