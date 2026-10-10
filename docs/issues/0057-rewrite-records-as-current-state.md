---
status: open
labels: [docs, verification]
opened: 2026-10-10
closed:
milestone: M3
blocked-by: []
resolution:
---

# 既存の検証の記録を、今の姿に書き直す

## 問題

ADR-0115 で、検証の記録（examples/*/verification.md）は今の結果だけを題目ごとに書くことにした。既存の記録は時系列に
節を足し、古い測定を日付つきの注記で打ち消している（examples/vera/verification.md には、前の節を打ち消す注記が 6 か所
ある）。

関係する ADR: ADR-0115, ADR-0116

## 完了の条件

- examples/*/verification.md のそれぞれを、今の結果だけの題目ごとの記録に書き直す。結果ごとに測ったときのコミットと日付を
  書く。VERA は、例を問題ごとのディレクトリに置き直すのと合わせて行う
- 書き直しで消える節を、開いている issue が節の題で引いていれば、引き方を記録のファイルとコミットに直す

## 経緯

- 2026-10-10: ADR-0115 を決めたときに起こした
- 2026-10-10: VERA の分を済ませた（ADR-0116 で例を p1/、p2/ に置き直し、examples/vera/verification.md を今の結果だけの題目ごとの記録に書き直した）。ほかの例の記録は残る
- 2026-10-10: 記録の規則は ADR-0117 に引き継いだ（ADR-0115 を置き換え。中身は同じ）
