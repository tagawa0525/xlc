---
status: open
labels: [docs]
opened: 2026-10-14
closed:
milestone:
blocked-by: []
resolution:
---

# 再測定の記録の残り（depletion の relock、本文の誤り、規則の扱い）

## 問題

高速化の後の記録の測り直し（マージ 9851184c、ブランチ docs/remeasure-after-speedups）のレビューで見つかり、
そのブランチでは直さなかったこと。

1. depletion の `chain` の relock が、rebase で落ちた。main の examples/depletion/dvc.lock は古いまま。#0124 は、
   `chain` が changed になった原因を「`bin/chain` が main で変わった」と書いているが、relock が落ちたことと
   区別していない。`dvc status` で changed になる理由は、`bin/chain` の変更か、relock の欠落か、両方かを、
   `dvc.lock` の deps の md5 と今の `bin/chain` の md5 を比べて確かめる（確かめていない）
2. マージ 9851184c のコミットの本文は、AIC の体系が ADR-0142 とライブラリの作り直しで「約 +250 pcm」動くと書いて
   いるのが誤りで、正しくは約 +325 pcm（examples/vera/verification.md の +324 pcm と同じ桁）。マージのコミットの
   本文は直せないので、ここに記録する
3. 再測定した記録に前の値との差を書いた。これは「検証の記録は今の結果だけ」（リポジトリの CLAUDE.md が引く記録の規則）と
   兼ね合いがある。また OpenMC のステージを飛ばして再測定した判断（#0120）も、規則に書かれていない。決定として残すなら、
   docs/adr/README.md か verification-record スキルに 1 行足す。

関係する ADR: なし

関係する issue: #0120、#0124

## 完了の条件

- 1 は、原因を確かめ、#0124 の書き方を合わせる（relock が必要なら #0124 の作業に含める）
- 3 は、規則として残すか、残さないと決めるかのどちらかを、ADR か skill への追記で示す

## 経緯

- 2026-10-14: 再測定のレビュー 15 周目で見つかり、起こした
