---
type: Issue
status: open
labels: [xs, vera, verification]
opened: 2026-10-15
closed:
milestone:
blocked-by: []
resolution:
---

# ENDF/B-VII.0 のライブラリでの KENO-VI との k の差が、VIII.0 での OpenMC との差より約 +130 pcm 大きい

## 問題

出典は examples/vera/verification-endfb70.md の「NJOY で処理した ENDF/B-VII.0 の OpenMC との比較（1B、1D、2B）」（xlc は d6cc90d9、OpenMC は
53be609d、比較は 3a599b6f、2026-11-20）。値はそこに書き、ここには未解決の点だけを並べる。

題名の差は、1B、1D、2B では、主に同じ ENDF/B-VII.0 を NJOY で処理した OpenMC と KENO-VI の差で、xlc の VII.0 の処理の誤りではなかった
（ライブラリを VIII.0 から VII.0 に替えたときの k の動きが、xlc と OpenMC で 25 pcm 以内で合う）。残る点:

- OpenMC（NJOY、VII.0）− KENO-VI が 600 K で約 +100 pcm、1200 K で約 +180 pcm ある原因を分けていない。候補（確かめていない）:
  - KENO-VI の連続エネルギーのライブラリ（AMPX）と NJOY の処理の差
  - KENO-VI の入力（公開されていない）とモデルの解釈の差。xlc と OpenMC は同じ解釈（ギャップは真空、Hf を除く）なので、2 つの比較では分けられない
  - 共鳴の散乱の扱い（自由ガスの範囲、DBRC）。温度で差が変わる分の候補（#0083）
- ライブラリを替えたときの k の動きが、xlc で OpenMC より +17〜+25 pcm 大きい（統計誤差の 2.2〜3.5 倍、3 問題とも同じ符号）。原因は確かめていない。
  核種ごとの動き（1B の xlc。記録の「核種を 1 つずつ VII.0 に替えた xlc の k の動き」）を OpenMC の側で同じ入れ替えと比べれば、核種に分けられる
- 測ったのは 1B、1D、2B だけ。吸収体のある集合体（2E〜2K）、565 K（2A、2Q。仕様書の KENO-VI は S(α,β) の温度を補正している）、3A、3B の
  OpenMC（VII.0）は無い

関係する ADR: ADR-0175

## 完了の条件

残る点を詰めるか、詰めないと決めて閉じる。KENO-VI の入力が無いので、OpenMC − KENO-VI の原因は xlc の側からは分けられない。+17〜+25 pcm を
核種に分けるなら、1B の OpenMC を核種を 1 つずつ VII.0 に替えたライブラリで解く（試行の入力の形は、ブランチ fix/0131-endfb70-keno-offset の
作業ツリーに作ったが実行もコミットもしていない。1 本は 100 万粒子・150 世代）。

## 検討して捨てた形

- xlc の VII.0 のファイルの扱い（共鳴の再構成、ドップラー拡がり、非分離共鳴、MF4/MF5/MF6）が差の主因という見方: NJOY で処理した同じ VII.0 の
  OpenMC が、版を替えたときに xlc と 25 pcm 以内で同じだけ動くので、主因ではない（上の出典）

## 経緯

- 2026-10-15: ADR-0175 の検証用ライブラリの最初の比較で、差が VIII.0 の OpenMC との差より大きいことが分かった。診断は使用量のため保留し、
  `examples/vera/endfb70/openmc/` のステージを用意した
- 2026-11-19: Slurm 961 で 1B、1D、2B の OpenMC（ENDF/B-VII.0）を解いた
- 2026-11-20: 3 者の比較から、差の主な部分が OpenMC（NJOY）と KENO-VI の間にあると分かった
