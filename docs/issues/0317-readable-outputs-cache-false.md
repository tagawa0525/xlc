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
  読む。`xs range` の範囲の TOML（`xs-range.toml`）は、dvc.yaml の outs が `cache: false` で git が追うので、DVC のキャッシュに入らず、
  中身の新しい版が `MISSING:` で止まっていた（ブランチ feat/0303-pin-depletion-by-endf-version の
  `examples/depletion/endfb70/xs-range.toml`、md5 edbeb1b9）
- 検査は、キャッシュに無く git の index が追う出力を、index の blob の md5 が lock と同じならその blob で読み、違えば `MISMATCH:` で
  止めるようになった（ADR-0258 の「## 変更」の 2026-11-20 の項）。残るのは、ブランチ feat/0303-pin-depletion-by-endf-version を main に
  rebase して、その run-checks.sh が通ることを確かめることである
- 直した検査（xlc bf410c388）を、ブランチ feat/0303-pin-depletion-by-endf-version の worktree（HEAD c6dcb3f01）で `--base main` で
  実行すると、`MISSING:` は無く、2145 ファイル（2231 の出力）を読んで終了コード 0 だった（読み手はそのブランチのコードから別の
  target に作った read-outputs、キャッシュは main の references/dvc-cache、手元、2026-11-20）。直す前の検査（main f73fda31e）は
  同じ条件で上の `MISSING:` を出して終了コード 1 だった
- `cache: false` の出力のうち git が追わないもの（中身がキャッシュにも git にも残らない）は、main（f73fda31e、2026-11-20）に無い。
  DVC の Python の API（dvc 3.67.1 の `Repo.index.stages` の outs の `use_cache`）で数えた `cache: false` の outs は 579 で、すべて
  git が追うファイル（ディレクトリは 0）で、hash は md5（DVC 3）だった。そうした出力が検査の対象に入れば、検査は `MISSING:` で止める
  （黙って飛ばさない）

関係する ADR: ADR-0258

## 完了の条件

- [x] `cache: false` で git で追う xlc の形式の出力を、検査が lock の md5 の版で読む（または止める理由を ADR に書く）
- [ ] ブランチ feat/0303-pin-depletion-by-endf-version の run-checks.sh が通る

## 検討して捨てた形

- xlc の形式の出力から `cache: false` をなくし、範囲の TOML を DVC のキャッシュに入れる: 採らない。git で差分を読める小さな TOML を
  キャッシュに移すと、レビューで中身が見えなくなり、検査の都合で出力の置き方を変えることになる
- HEAD の blob を読む（このブランチの最初の実装）: 捨てた。マージの hook（pre-merge-commit）は main の worktree で run-checks.sh を
  呼び、そこでは HEAD が main のままで、index と作業ツリーの lock がマージの結果になる。新しい範囲の TOML を足すブランチでは、HEAD に
  そのパスが無いか古い版で、検査が直したいマージそのものを止める（ローカルのレビューの指摘。index を読む形でテストを足した）

## 経緯

- 2026-11-20: 版ごとのピンセルの燃焼（#0303）のブランチのチェックで、VII.0 の範囲の TOML が止まって見つけた
- 2026-11-20: 検査が、キャッシュに無く git の index が追う出力を、lock と同じ md5 の blob で読むようにした。#0303 のブランチの
  worktree で検査が通ることを確かめた
