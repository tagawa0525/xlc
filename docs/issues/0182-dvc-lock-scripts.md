---
status: closed
labels: [dvc, tooling]
opened: 2026-10-28
closed: 2026-10-28
milestone:
blocked-by: []
resolution: コミット f2ddea9ff
---

# dvc.lock の比較とステージ数の計数をスクリプトにする

## 問題

`dvc commit` で `dvc.lock` を付け直したブランチ（ADR-0210 の移行、#0145 の 3 度の付け直し、bin の作り直し）は、
deps の hash だけが変わり、ステージの集合、`cmd`、`outs`、`params` は変わらないはずである。これを毎回その場の
スクリプトで確かめている（2026-10-28 に main..feat/0145 の 12 本の lock、807 ステージを比べ、違いは直下の
`dvc.lock` の python-ast の 33 ステージの追加だけだった）。手順（/dvc-relock）にも CLAUDE.md にも、機械で確かめる
ツールが無い。

`dvc repro` の前に changed のステージを数える手順（CLAUDE.md の「計算」）は、`dvc status` の文字列を grep で
数えると、タブで始まる `changed deps:` の行を拾って誤る（2026-10-28 に 468 を 942 と数えた）。

ツールを作るには、`dvc.lock` が YAML なので YAML を読む手段が要る。default の devShell の `python3` は `import yaml` も
`import ruamel.yaml` も失敗する（2026-10-28 に確かめた。references の devShell も同じく失敗する。dvc 3.67.1 の
Python は別の closure で、PYTHONNOUSERSITE と PATH の都合から検証用の Python からは見えない: flake.nix の
`dvc` の組み立て）。

関係する ADR: なし

## 案

1. `scripts/compare-dvc-lock.py <base> <head>`: 2 つの git の rev の間で変わった `dvc.lock` すべてについて、ステージの
   集合と、各ステージの `cmd`、`outs`、`params` が同じことを確かめ、終了コードで判定する。違いはファイル、ステージ名、
   欄を出す。ステージの集合の差は違いとして出す（緩めるフラグを作らない。ADR-0210 の python-ast の 33 ステージの
   ように増えて当然のときは、出力を人が読む）。
2. `scripts/count-changed-stages.py [yaml...]`: `dvc status --json` を読み、yaml ごとに changed deps と changed outs の
   ステージの数を出す。
3. default と references の devShell の Python に pyyaml を足す（flake.nix の `pythonPackages`）。標準ライブラリに
   YAML を読む手段が無く、手書きの読み込みは保守コストが大きい。pyyaml は Python のエコシステムで事実上の標準で、
   CLAUDE.md の「設計指針」の依存の規則に当たるので、新しい決定にせず ADR は書かない。
4. /dvc-relock の、lock を付け直した後の確認に 1 を、`dvc status` の確認に 2 を、CLAUDE.md の「計算」の見積もりに
   2 を書く。

## 完了の条件

- 2 本のスクリプトとテスト（`scripts/test_*.py`）が main に入り、`scripts/compare-dvc-lock.py main HEAD` が lock を
  変えないブランチで違い無しで終わる
- `scripts/count-changed-stages.py` が main の `dvc status` の changed のステージ数と一致する（2026-10-28 に
  00d789a1f で 630）
- /dvc-relock と CLAUDE.md が 2 本を呼ぶ

## 経緯

- 2026-10-28: 起こした（ユーザー 08:14「推奨に従います」）
- 2026-10-28: 実装した。`dvc status --json` の値は辞書か、中身の無い種類の文字列（`always changed`）で、種類に
  `changed command` もある。00d789a1f の `scripts/count-changed-stages.py` は 630 ステージ（changed deps 626、changed outs 10、
  changed command 1、always changed 4）で、`dvc status --json` のキーの数と一致した。flake.nix の変更は `build` の deps なので、
  マージの後に `build` の lock を付け直す（`dvc repro -s build`）。
