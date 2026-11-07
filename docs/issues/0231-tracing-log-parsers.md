---
status: closed
labels: [tooling, examples]
opened: 2026-11-06
closed: 2026-11-08
milestone:
blocked-by: []
resolution: wontfix: レビューで直さなかった指摘を集めた issue は作らない決定（nixfiles の ADR-0066）に合わせ、再現のある項目が無いので閉じた
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

3 周目のレビューで残した p9_timing.py と記録の表の点:

- モジュールの docstring の「部分挿入のノードの定数を…`eigenvalue` に入る」の後の文が前の文とつながり、「互いに重ならない部分」の主語が抜けている
- `other`（全体から span の和を引いた残り）を「入力の読み込み、出力など」と書くが、span の外の計算（`layout.problem`、`partial_nodes`、
  `solve_rodded` の部分挿入の準備、`fission_rates`、燃焼度の更新、点の記録）も入る
- 解の表は点の `days` を持たないので、同じ EFPD の 2 つの点（出力のステップの終わりと崩壊の後）を見分けられない（記録の本文で補った）
- 出力のステップの終わりの点の「平衡 Xe の解き直し」は、前の状態（ADR-0218 の分岐）の解の回数で、点のホウ素を探した状態の回数ではない（記録の本文で補った）
- `nem.count` は span のすべてを数え、`nem.iterations` は `iterations` のフィールドを持つ span だけから範囲を取る（失敗した解の span は数えが合わない）
- core の case.rs の内挿の span は `enter()` のガードで、depletion.rs の `xenon anderson` は `in_scope`（書き方が 2 つ）
- examples/test_verification_tables.py の `cells()` と、同じ式を書いた 10 か所ほどの行

関係する ADR: なし

## 完了の条件

ログを読む共通の部分（span の並びとフィールド、時間の単位）を 1 つのモジュールにまとめ、3 つのスクリプトをそれに移すか、移さない
理由を書く。燃焼点の値を履歴から読むかを決める。上の 3 周目の点を直すか、直さない理由を書く。

## 経緯

- 2026-11-06: #0156 の全炉心の燃焼の時間の内訳のブランチ（chore/0156-p9-timing）の 2 周目と 3 周目のレビューで挙がった
- 2026-11-08: レビューで直さなかった指摘を集めた issue は作らない決定（nixfiles の ADR-0066 の「## 変更」の 2026-11-08）に合わせて閉じた。再現のある項目は無い
