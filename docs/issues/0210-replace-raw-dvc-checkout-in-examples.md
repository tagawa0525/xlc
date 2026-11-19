---
type: Issue
status: closed
labels: [dvc, tooling]
opened: 2026-10-31
closed: 2026-11-19
milestone:
blocked-by: []
resolution: コミット b06377510
---

# examples の手順と文言に生の dvc checkout が残る

## 問題

ADR-0223 で、Claude Code の Bash ツールからの生の `dvc checkout` を `.claude/settings.json` の hook で止め、
`scripts/dvc-checkout.py` を使わせることにした。スキル、CLAUDE.md、docs のコマンドは同じブランチで直したが、そのブランチの
変える範囲から外した examples/ には、生の `dvc checkout` を指す手順と文言が残る（2026-10-31 に `grep -rn 'dvc checkout' examples`）:

- 記録の再現の手順: examples/c5g7/verification.md、examples/colorset/verification.md、examples/mox/verification.md、
  examples/split/verification.md（`dvc checkout examples/vera/dvc.yaml`）、examples/depletion/verification.md の 2 か所
  （`dvc checkout examples/depletion/p2/dvc.yaml`）
- dvc.yaml のコメント: examples/c5g7、colorset、mox、split の dvc.yaml
- スクリプトのエラーの文言: examples/aws/aws_run.sh（3 か所）、examples/vera/remote_send.sh

hook に止められれば拒否の文言がスクリプトを示すので作業は止まらないが、手順をコピーして実行するたびに 1 回止められる。

関係する ADR: ADR-0223

## 完了の条件

上の手順と文言を `scripts/dvc-checkout.py <指定>` の形に直す（記録の本文は数値を変えないので、記録の規則の再測定は要らない）。

## 経緯

- 2026-10-31: ADR-0223 のブランチ（feat/dvc-checkout-hook）で、examples/ を変える範囲の外としたので書いた
- 2026-11-19: 閉じた。上の手順と文言を scripts/dvc-checkout.py の形に直した。examples/depletion/verification.md の 3 か所（764268e1、50f94ebfa、
  af044d5aa の版で戻す手順）は、その版に scripts/dvc-checkout.py が無い（2026-10-28 に足した）ので、その版の worktree で今のスクリプトを
  パスで呼ぶ形にした（スクリプトは cwd の `git rev-parse --show-toplevel` で実行する）。`grep -rn "dvc checkout" examples` に生の形は残らない
