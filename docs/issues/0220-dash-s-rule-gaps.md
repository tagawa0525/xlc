---
status: open
labels: [dvc, docs]
opened: 2026-11-03
closed:
milestone:
blocked-by: []
resolution:
---

# ADR-0184 の -s の例外の一覧と確かめ方の残り

## 問題

ADR-0184 の `-s` の 2 つ目の例外（「## 変更」の 2026-11-03）を足したブランチ docs/0184-skip-upstream-for-quick-downstream の
レビューで出た、直していない 2 つ:

- CLAUDE.md の「計算」は `-s` を /dvc-relock と ADR-0184 の 2 つ目の例外だけに限るが、.claude/skills/openmc-reference/SKILL.md は
  `dvc repro -s examples/vera/dvc.yaml:openmc-1B-model` を、例外の確かめ方を書かずに手順に書き、上流を「`references-shell` だけ」と
  書く。examples/vera/dvc.yaml の openmc-1B-model の deps には、`references-shell.txt` のほかに直下の `python-ast` の出力
  （`openmc_model.py.ast` など 5 本）がある。ADR-0184 の 2 つ目の例外（python-ast を先に実行し、`check-references-shell.sh` を
  通す）に当たるときの使い方として書き直すか、`-s` を使わない手順に直すかは決めていない。食い違いはこのブランチの前の main からある
- 例外の確かめ方は、`dvc status --with-deps --json <yaml>:<ステージ>` の出力を人が読み、目標のステージ自身、`references-shell` の
  `always changed` の項目を除いて種類を照らす手作業である。このブランチのレビューの 1 周目で、確かめ方の文が `references-shell` を
  見落としていた。`scripts/count-changed-stages.py` は `dvc status --json` を種類ごとに読むので、「目標を除く上流が `changed deps` だけ
  （`always changed` の項目を除く）でなければ 0 でない終了コード」の判定を足せば、`dvc status` で分かる部分は文で書き直さずに済む。
  ADR-0184 の確かめ方のうち `dvc status` に出ないもの（`references-shell.txt` の `check-references-shell.sh`、python-ast を先に
  実行すること）は、その判定だけでは確かめられない

関係する ADR: ADR-0184

## 完了の条件

- openmc-reference の `-s` を、ADR-0184 の例外に足すか、`-s` を使わない手順に直す
- 例外の確かめ方を機械で判定するツールを足すか、足さない理由を書く

## 経緯

- 2026-11-03: ブランチ docs/0184-skip-upstream-for-quick-downstream のレビュー（1 周目、2 周目）の指摘から起こした
