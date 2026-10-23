---
status: accepted
date: 2026-10-15
supersedes: []
superseded-by:
issues: ["#0131"]
requires: [ADR-0020, ADR-0031, ADR-0079, ADR-0080, ADR-0124, ADR-0176]
---

# ADR-0175: 仕様書に KENO-VI の値がある問題は、ENDF/B-VII.0 の検証用ライブラリでも xlc を解いて KENO-VI と直接比べる

## 背景

仕様書（CASL-U-2012-0131-004 Rev.4）の参照解は KENO-VI、ENDF/B-VII.0 で、Problem 1〜5、4-2D、REF1-2D、5-2D の k、
ピン・集合体出力、軸方向の出力、制御棒価値、ITC がある。燃焼（Problem 9）と HFP（6〜8）にはない。

xlc の本番のライブラリは ENDF/B-VIII.0（ADR-0020）で、KENO-VI との差には核データの版の差が混じる（3B で
OpenMC − KENO-VI が +124 pcm、2Q で +119 pcm。examples/vera/verification.md）。それを除くため、同じ VIII.0 の
OpenMC の参照解を作ってきた（ADR-0031。ツールの固定は ADR-0124 が引き継いだ）。3D の問題では OpenMC が 1 本 5 時間以上かかる。

判断の軸:

- 手法の差を、核データの差と分けて数値で言えること
- 計算時間（3D の OpenMC は 1 本 5 時間以上）
- 本番のライブラリと、本体のコードを変えないこと

## 検討した案

- 案 A: OpenMC の VIII.0 の参照解を作り続ける。手法の差は取り出せるが、仕様書に KENO-VI の値がある問題でも 3D は
  1 本 5 時間以上かかり、モデルの解釈（KENO-VI の入力は公開されていない）の差も OpenMC のモデルに別に入る
- 案 B: KENO-VI とだけ比べる（本番の VIII.0 のまま）。安いが、差に核データの版の差が混じり、原因を分けられない
- 案 C（採用）: 仕様書に KENO-VI の値がある問題は、xlc を ENDF/B-VII.0 の検証用ライブラリでも解き、KENO-VI と直接比べる。
  核データの版は揃うので、残る差は手法と処理の差になる。OpenMC の参照解は、仕様書に参照がない問題と、診断に反応率が
  要る小さな 2D の問題だけにする

## 決定と理由

案 C にする（ユーザーの決定、2026-10-15）。

- 検証用ライブラリは、本番の `examples/vera/xs.toml` と同じ核種・温度・群で、ENDF だけ VII.0 にしたもの
  （`endfb70/xs.toml`）。本番のライブラリ（VIII.0）と、その結果は変えない。コードも変えない。ステージは `examples/vera/endfb70/dvc.yaml`（入力は `endfb70/xs.toml`、
  3A・3B の入力は `p3_inputs.py --library-dir endfb70` が `endfb70/p3` に書く）
- 仕様書に KENO-VI の値がある問題（1B〜1D、2A〜2K、2Q、3A、3B）を VII.0 のライブラリで解き、KENO-VI と比べる
- OpenMC の参照解は、仕様書に参照がない問題（燃焼、分岐、HFP）と、原因の診断に反応率が要るときの 2D の小さな問題だけ作る。
  3D の全炉心と Problem 4 の OpenMC は作らない（ADR-0031 の帰結を維持する）
- ADR-0124 と ADR-0031 の関係: ADR-0124 はツールの固定と `references/` の扱い（案 A、D、F、J、N）を決め、参照解の作り方を
  ADR-0031 から引き継いだ。この ADR は、その引き継いだ「参照解の作り方」のうち、KENO-VI の値がある問題の比べ方だけを
  足す。ADR-0124 の決定（固定の手段）は変えないので、置き換えず `requires` にする。VIII.0 の OpenMC との比較の
  位置づけ（コード同士の検証）も変えない
- 命名の規則: 評価済みデータごとに `examples/vera/<略称>/` と `verification-<略称>.md` を置く（ENDF/B-VII.0 は `endfb70`。
  将来の JENDL-5、JEFF-3.3 は `jendl5`、`jeff33`）。本番（ENDF/B-VIII.0）は `examples/vera/` の直下のまま。版の違う結果を
  1 つの記録に混ぜない

比べるときの注意:

- 同じ評価済みデータでも、KENO-VI の連続エネルギーのライブラリは AMPX で処理したもので、xlc の処理（共鳴の
  再構成、ドップラー拡がり、S(α,β) の温度の内挿。ADR-0080）とは違う。VII.0 での xlc − KENO-VI に残る差は、手法と処理の差で
  あり、それが見たいものである。核データの版の差は除けるが、処理の差は除けない
- 3B の KENO-VI は、H in H2O の S(α,β) を 550 K と 600 K の 2 つの計算から手で内挿している（仕様書の問題 3 の
  「Cross Sections」）。2A、2Q、5-2D、HZP の 565 K も、仕様書は S(α,β) の温度の補正を加えている。xlc は S(α,β) の群定数を
  温度について線形に内挿する（ADR-0080）。KENO-VI は k の値を内挿し、xlc は群定数を内挿するので、同じ内挿ではない。
  差が内挿の違いに当たるかは確かめていない
- ENDF/B-VII.0 には C-12 と C-13 のファイルがなく、天然の C（`n-006_C_000.endf`）だけがある。ADR-0079 は VIII.0 に
  天然の C がないために分けたので、VII.0 のライブラリでは C-12 と C-13 の両方にこの天然の C のファイルを指す
  （入力 `endfb70/xs.toml` で吸収する。コードは変えない）。KENO-VI も天然の C の評価を使うと見込むが、確かめていない

## 確認

xlc の xs は ENDF/B-VII.0 のすべての核種で止まらず、ピンセル 1B〜1D、集合体 2A〜2K、2Q、3A、3B を解けた
（43661306、2026-10-15。`examples/vera/endfb70/`）。結果は `examples/vera/verification-endfb70.md`。

- k: ENDF/B-VII.0 の xlc − KENO-VI は 15 の格子で +115〜+270 pcm（平均 +165）、3A が +181 pcm、3B が +166 pcm。VIII.0 の
  xlc − OpenMC（平均 +39 pcm）と同程度にはならなかった。核データの版の差だけで残りが説明できる、という見込みは k では
  成り立たない。原因は確かめていない（#0131）
- ピン出力: 吸収体のある集合体で、xlc − KENO-VI の RMS は VIII.0 の 0.45〜0.61 倍に下がる

## 帰結

- 仕様書に KENO-VI の値がある問題は、OpenMC を実行せずに手法の差を見られる
- 検証用ライブラリは本番ではない。VII.0 の zip は ADR-0124 と同じく URL と sha256 で固定し、`.#endf` の木に
  `ENDF-B-VII.0_neutrons/`、`ENDF-B-VII.0_thermal_scatt/` として足す（docs/environment.md）。VIII.0 のディレクトリの中身は変わらない
- 処理の差（AMPX と xlc）は分けられない。分けるには、同じ処理の参照（NJOY の VII.0 での OpenMC）が要り、2D の小さな
  問題でだけ作る
- ENDF/B-VII.0 での k の差が VIII.0 の OpenMC との差より大きい原因は未解決（#0131）。NJOY で処理した ENDF/B-VII.0 の OpenMC
  （2D の小さな問題）が、分ける次の手になる
