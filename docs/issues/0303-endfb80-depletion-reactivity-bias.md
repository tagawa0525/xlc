---
type: Issue
status: open
labels: [xs, depletion, verification]
opened: 2026-11-20
closed:
milestone:
blocked-by: []
resolution:
---

# ENDF/B-VIII.0 が軽水炉の燃焼で反応度を低く見積もるという報告を確かめ、本番の核データの版を決める

## 問題

本番の核データは ENDF/B-VIII.0（ADR-0020）で、燃焼の検証は同じ VIII.0 の OpenMC と比べている。この比較には核データそのものの偏りは
表れない（xlc と参照の両方に同じように入る）。偏りが表れるのは、実機の炉心追跡（Problem 9 の臨界ホウ素濃度）と比べるときである。

文献の報告（2026-11-20 に検索の結果の要約で読んだ。論文と資料の本文は開いておらず、値は確かめていない）:

- ORNL（Kim、Wieselquist、J. Nucl. Eng. 2(4), 2021、<https://www.mdpi.com/2673-4362/2/4/26/html>）: PWR の燃焼で、VIII.0 は VII.1 に比べて
  反応度を低く見積もり、差は燃焼とともに広がる。ORNL の CSEWG 2020 の資料
  （<https://indico.bnl.gov/event/7233/contributions/43847/attachments/31698/50097/CSEWG2020_12022020_r0.pdf>）では 60 MWd/kgU で
  約 600 pcm（Serpent）。効く核種は Pu-239、U-235、O-16、U-238。発電所のデータとの比較でも過小評価が見えるとしている
- IAEA（Trkov ら。検索の結果が挙げた資料は INDC(NDS)-0893、<https://www-nds.iaea.org/publications/indc/indc-nds-0893.pdf> と、Capote の
  WONDER 2023 の要旨。どちらの記述かは確かめていない）: 原因の中心を U-238 の共鳴領域の捕獲の評価に置く
- Studsvik（Ferrer、Wemple、Mini-CSEWG の資料、
  <https://indico.bnl.gov/event/18702/contributions/75746/attachments/47221/80101/ferrer-mini-csewg-e8r1.pdf>）: VIII 系は k∞ と全炉心の
  k-eff を低く見積もる傾向があり、Pu-239 を JENDL-4 の評価に替えると VII.1 に近づく
- ENDF/B-VIII.1 の論文（<https://arxiv.org/pdf/2511.03564>）: 燃焼に伴う反応度の低下のベンチマーク（検索の結果では IRPhE/DUKE）で
  2σ 以内に入り、高い燃焼度で改善した

決まっていないこと、確かめていないこと:

- 報告の偏りが xlc の本番に効く大きさかどうか。報告の値は 60 MWd/kgU のもので、実測と比べる Problem 9 のサイクル 1 は 16.939 GWd/MT まで
  （docs/plans/m4-depletion.md）。差は燃焼とともに広がるとされるので、この範囲での大きさは別に見る
- 本番の版をどれにするか（候補）: VIII.1 に移る、VII 系にする、VIII.0 のまま偏りを知ったうえで続ける
- xlc が VIII.1 の評価を読めるか（ADR-0167 は、VIII.1 が使いうる共鳴の表現に未対応のものがあると書いている）。燃焼には燃焼チェーンの
  全核種が要るので、新しい燃料の核種だけでは分からない
- VII 系の手持ち: VII.0 は検証用ライブラリ（ADR-0175）があるが、中性子と熱散乱のファイルと VERA の格子の核種だけで、崩壊と核分裂収率の
  ファイル、燃焼チェーンの核種は無い（docs/environment.md）。報告の比較の相手は VII.1 で、VII.0 が同じ傾向かは確かめていない。VII.1 の
  ライブラリは無い
- 版を替える費用: xlc のライブラリ、燃焼チェーン（崩壊と収率のファイルから作る。ADR-0164。収率の重みは VIII.0 で測った。ADR-0183）、
  燃焼の履歴、テーブルに加えて、OpenMC の参照解（燃焼を含む）を同じ版で作り直す（手法の検証は同じ核データで比べる。ADR-0175）
- k の手法の差を何と比べて見るか（#0302）は、本番の版で選択肢が変わるので、合わせて決める
- 出力の規格化（#0293 の項目 3）も版に依る。チェーンの核分裂の Q（MF3 MT18 の QI）は、VIII.1 で U-235 が VIII.0 より +0.35%、VII.1 で Pu-239 が
  +0.16% 違い、捕獲のエネルギーの割合は捕獲の断面積で変わる（examples/depletion/verification.md の「出力の規格化を OpenMC の燃焼で確かめる
  （項目 3）」）。ENDF/B-VIII.1 の論文の XIV.D.4 は、核分裂あたりのエネルギーの仮定が燃焼に伴う反応度の傾きに最も強く効くと書く。版の燃焼の
  反応度の差を測るときと文献と比べるときは、規格化の定義をそろえる

関係する ADR: [ADR-0262](../adr/0262-apply-natural-evaluations-to-listed-isotopes.md)（版ごとのチェーンの天然の元素の評価）、[ADR-0263](../adr/0263-list-isomers-with-unknown-excitation.md)（VII.0 の Nb-88m の負の ELIS）。本番の版を決めた決定は「問題」の冒頭に挙げた

## 完了の条件

- [ ] 文献の値を本文で確かめる
- [ ] xlc が VIII.1 の評価を読めるかを、VERA の 1B の核種、次に燃焼チェーンの全核種で試す
- [ ] 1 つのピンセルの燃焼で、版による反応度の差を、Problem 9 の範囲（約 17 MWd/kg）と 60 MWd/kg までの燃焼度ごとに測る。比べる版
  （VIII.0 と VIII.1。文献の値と比べるなら VII.1 も）は、上の読めるかの結果で決める
- [ ] 本番の版を決め、ADR-0020 の「## 変更」か新しい ADR に書く。本番の作り直しの前に決める

## 経緯

- 2026-11-20: #0131 の結果から本番の版をユーザーに尋ねたとき、ユーザーが VIII の燃焼の問題の報告があるかを問い、文献を検索して起こした
- 2026-11-20: 燃焼チェーンを VII.1、VIII.0、VIII.1 から作れるようにした（複数の物質のテープ、天然の元素の評価は ADR-0262、ファイルの固定は docs/environment.md、入力は examples/depletion/endfb71、endfb81）。VII.0 は decay の Nb-88m の負の励起エネルギーで止まる。核種の数、反応を持つ核種、核分裂する核種、警告の数は、VII.1 が 3,820、423、85、115、VIII.0 が 3,820、555、88、124（main のチェーンと同じバイト）、VIII.1 が 3,820、556、88、125（`bin/chain` は 920b9f156、2026-11-20。手元の実行で、ステージではない）
- 2026-11-20: 燃焼チェーンを VII.0 からも作れるようにした（Nb-88m の負の励起エネルギーは ADR-0263、入力は examples/depletion/endfb70）。VII.0 の核種の数、反応を持つ核種、核分裂する核種、警告の数は 3,837、398、65、108 で、ほかの 3 つの版の数とバイトは上と同じ（`bin/chain` は 22f0c8d7f、2026-11-20。手元の実行で、ステージではない）
