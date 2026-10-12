---
status: open
labels: [numeric]
opened: 2026-10-12
closed:
milestone:
blocked-by: []
resolution:
---

# 1ノードの問題の収束の判定が NaN を後の値で消す

## 問題

`xlc_numeric::nem::homogeneous_faces`（ADR-0103 の1ノードの問題）は、面の中性子束の変化の最大を
`change = if relative <= change { change } else { relative }` で取り、NaN も収束していないとみなすとコメントに書く。
しかし NaN のあとに有限の変化が来ると `relative <= NaN` は偽なので、`change` は有限の値に置き換わり、NaN を含む面の
中性子束で収束したと判定しうる（返したあとの正の値の検査が多くは拾う）。

制御棒の部分挿入（#0063）のレビューで、同じ形を写した core の z の1ノードの問題（crates/core/src/rods.rs）で見つかり、そちらは
`relative.is_nan() || relative > change` に直した。numeric の元は lattice の ADF と core が共有するので、別に直す。

関係する ADR: なし

## 完了の条件

numeric の判定を、一度 NaN になれば収束としない形に直し、NaN のあとに有限の変化が来る場合のテストを足す。

## 経緯

- 2026-10-12: #0063 のレビュー（core の写し）で見つけた
