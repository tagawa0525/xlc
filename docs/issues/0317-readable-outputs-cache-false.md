---
type: Issue
status: open
labels: [dvc, checks]
opened: 2026-11-20
closed:
milestone:
blocked-by: []
resolution:
---

# 形式を変える修正の検査が、git で追う cache: false の出力をキャッシュに無いとして止める

## 問題

- `scripts/check-readable-outputs.py`（ADR-0258 の「形式を変える修正の検査」）は、上流の出力の中身を DVC のキャッシュから lock の md5 で
  読み、キャッシュに無い出力を `MISSING:` で止める。`xs range` の範囲の TOML（`xs-range.toml`）は、examples/depletion の dvc.yaml で outs が
  `cache: false` で、git で追う。`cache: false` の出力を DVC はキャッシュに入れない
- main の既存の `xs-range.toml`（md5 74fbedef。examples/depletion と p2）は、キャッシュに 2026-10-14 の版が残っているので読めている。
  新しい内容の `cache: false` の出力は、キャッシュに入る経路が無く、止まる
- ブランチ feat/0303-pin-depletion-by-endf-version で、run-checks.sh の rust docs-test が
  `MISSING: examples/depletion/endfb70/xs-range.toml (md5 edbeb1b9cf80977e5899d87ac98fbc76): DVC のキャッシュに無い` で止まった
  （2026-11-20。ブランチが flake.nix を変えたので検査が読む。作業ツリーの版の md5 は edbeb1b9 で、lock と同じ）
- 修正方法の候補（決めていない）: `cache: false` の出力は、git の HEAD の blob を lock の md5 と照らして読む（作業ツリーを読まないという
  ADR-0258 の理由を保つ）。または、xlc の形式の出力から `cache: false` をなくす

関係する ADR: なし

## 完了の条件

- [ ] `cache: false` で git で追う xlc の形式の出力を、検査が lock の md5 の版で読む（または止める理由を ADR に書く）
- [ ] ブランチ feat/0303-pin-depletion-by-endf-version の run-checks.sh が通る

## 経緯

- 2026-11-20: 版ごとのピンセルの燃焼（#0303）のブランチのチェックで、VII.0 の範囲の TOML が止まって見つけた
