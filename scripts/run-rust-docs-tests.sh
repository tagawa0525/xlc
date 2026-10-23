#!/usr/bin/env bash
# 文書だけの差分で run-checks.sh が実行する Rust のテスト（nixfiles の ADR-0057）。
# crates/docs のテストが文書どうしの整合と、検証の記録の表の `verification_tables.py --check` を実行する。
# design-records check は run-checks.sh が別の段で実行するので含めない。
set -euo pipefail

cd -- "$(git rev-parse --show-toplevel)"

exec cargo test -p xlc-docs
