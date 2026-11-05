---
status: open
labels: [tooling, examples]
opened: 2026-11-06
closed:
milestone:
blocked-by: []
resolution:
---

# tracing のログを読む例のスクリプトの共通化

## 問題

tracing-subscriber の span の閉じる行（`close time.busy=…`）を読むスクリプトが 3 つあり、それぞれが正規表現と時間の単位の読み方を
持つ: examples/depletion/p9-timing/p9_timing.py（ns、µs、ms、s）、examples/depletion/trials-lattice-scaling/summarize.py、
examples/depletion/issue-0045/compare.py（ms と s だけ）。ログの書き方が変わると 3 か所を別に直すことになる（#0156 の
ブランチの 2 周目のレビューの指摘）。

同じレビューで残した関係する点:

- p9_timing.py は燃焼点の EFPD、ホウ素、固有値計算の回数もログから読む。同じ値は同じステージが書く燃焼の履歴（`xlc-core-history`
  の点の `efpd`、`boron`、`solves`）にもあり、ログの文言が変わると黙って読めなくなりうる。ログにしか無いのは中点と平衡 Xe の
  解き直しの回数
- core の `nem::record_solve` は `tracing::Span::current()` に記録するので、`nem` の span をフィルタで外し親の span を残すと、記録が
  親の span に向かい（フィールドが無いので捨てられる）、p9_timing.py は丸めた time.busy を使う。今のステージは `RUST_LOG` を
  付けない

関係する ADR: なし

## 完了の条件

ログを読む共通の部分（span の並びとフィールド、時間の単位）を 1 つのモジュールにまとめ、3 つのスクリプトをそれに移すか、移さない
理由を書く。燃焼点の値を履歴から読むかを決める。

## 経緯

- 2026-11-06: #0156 の全炉心の燃焼の時間の内訳のブランチ（chore/0156-p9-timing）の 2 周目のレビューで挙がった
