---
type: Issue
status: open
labels: [xs, endf]
opened: 2026-11-20
closed:
milestone:
blocked-by: []
resolution:
---

# 対数の補間則が、ln をとる側の正でない値で NaN を黙って返す

## 問題

- ENDF の補間則の式（`Interpolation::interpolate`、crates/xs/src/endf/tab1.rs）は、ln をとる補間則（INT=3 lin-log、
  INT=4 log-lin、INT=5 log-log）で、ln をとる側の y が負か、符号の違う 2 点か、ln をとる側の x が正でないとき、NaN を
  返す。NaN はどの比較でも偽なので、後の検査（ADR-0267 の 0 K の負の値の検査など）もすり抜ける。ENDF/B-VII.0 の Cf-253 の
  群定数の NaN（#0320）は、非分離共鳴の遮蔽した捕獲が 8500 eV で正、10 keV で負の区間を log-log で補間した値が直接の
  原因（別のセッションの調べ。このブランチでは式を読んだだけで、Cf-253 を流して確かめていない）
- 読み込みの TAB1 は `Tab1::check` が負の y と正でない x を止めているが、計算の途中で作る表（非分離共鳴の自己遮蔽の表）と、
  TAB1 でない値の並び（MF4 と出射分布の入射エネルギーの間）は、この検査を通らずに式に渡る
- ENDF-102（references/papers/endf-102-bnl-224854-2023-inre.txt の 0.5 節の表 16）は INT=3〜5 を「ln(x) か ln(y) が線形」と
  定義するだけで、正でない値の扱いを書かない。0 の値は、INT=6 の説明が「σ1 か σ2 が 0 なら結果は明らかに 0」と書き、
  MF7 の説明が「S = 0.0 はとても小さい数を表す」と書く。NJOY の terp1（references/NJOY2016/src/endf.f90）は、log-log で
  y1 = 0 なら 0 を返し、y2 = 0 なら exp(−∞) で 0 になる。負の値は扱わない

関係する ADR: なし

## 完了の条件

- [ ] 合成の入力のテストで、ln をとる補間則の定義域の外の区間を持つ表がエラーで止まり、`interpolate` が NaN を返さない
- [ ] 4 つの版の全ファイルで、ln をとる補間則の区間に正でない値を持つ表を数え、既定をエラーにして止まる評価を挙げる
- [ ] 既存の核種の出力が変わらない（`xs nuclide` の出力のバイトが main のビルドと同じ）
- [ ] 読み込みの走査（read-survey）で新しく止まる評価が無いかを確かめる

## 経緯

- 2026-11-20: Cf-253 の NaN（#0320）の直接の原因として、別のセッションの担当が見つけた
