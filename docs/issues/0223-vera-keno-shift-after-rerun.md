---
status: open
labels: [verification, lattice, xs]
opened: 2026-11-04
closed:
milestone: M3
blocked-by: []
resolution:
---

# VERA の格子と P3 の k が、ADR-0205 と ADR-0209 の後の再実行で KENO-VI に対して +34〜+154 pcm 動いた

## 問題

examples/vera と examples/vera/endfb70 のライブラリと格子を、ブランチ chore/rerun-after-adr-0226（流したコミットは 4cfe4f58、分岐点は
main の b4e9b2cd、2026-11-04、Slurm 678〜681）で今の bin で作り直すと、xlc の k が一様に上がり、KENO-VI との差が大きくなった。
前の実行は 05714cea（2026-10-23。examples/vera/verification-endfb70.md の前の版の「計算」）。

xlc − KENO-VI（pcm。比較の JSON の `k.diff_pcm`、P3 は `k.diff_keno_pcm`。前 → 今）:

- ENDF/B-VIII.0（本番。`examples/vera/compare-keno-*.json`）: 1B −3.8 → +48.4、1C −31.2 → +66.5、1D −32.9 → +120.6、2A +36.8 → +90.8、
  2B +26.7 → +78.7、2C +15.2 → +108.8、2D −2.9 → +141.4、2E +97.9 → +142.4、2F +145.6 → +187.1、2G +159.1 → +195.0、2H +93.7 → +128.0、
  2I +45.9 → +94.1、2J +153.8 → +195.2、2K +98.7 → +141.6、2Q +146.3 → +202.9（xlc の k の変化は +34〜+154 pcm）
- ENDF/B-VII.0（`examples/vera/endfb70/compare-*.json`）: 1B +167.7 → +202.3、1C +134.1 → +213.9、1D +127.1 → +262.3、2A +145.7 → +185.2、
  2B +154.4 → +191.4、2C +137.2 → +215.3、2D +114.2 → +242.6、2E +121.8 → +154.4、2F +120.4 → +151.4、2G +80.3 → +109.5、2H +104.9 → +131.4、
  2I +142.2 → +176.7、2J +128.3 → +159.2、2K +126.2 → +158.4、2Q +245.1 → +285.6（xlc の k の変化は +27〜+135 pcm）
- P3（`examples/vera/p3/compare-3?.json`、`examples/vera/endfb70/p3/compare-3?.json`）: VIII.0 の 3A +62.0 → +110.1、3B +164.5 → +213.1、
  VII.0 の 3A +157.4 → +191.8、3B +136.7 → +172.8（xlc の k の変化は +34〜+49 pcm）

前と今の実行の deps（examples/vera/dvc.lock の md5）: `bin/lattice`（lattice-1B）457aae7a → 48387d65、`bin/xs`（xs）d3dc3807 → 953c587b、
本番のライブラリ `vera.msgpack` 05132f1d → 0dfa3aca、VII.0 のライブラリ（examples/vera/endfb70/dvc.lock）1d13e690 → 022595ef。前の
ライブラリは今の `bin/lattice` で読めない（`missing field sab_transfer`。ADR-0205 の後の形式）ので、ライブラリと格子の効きを今の bin で
分けることはできない。

関係する ADR: なし

## 候補（確かめていない）

前の実行と今の実行の間で、格子の k を動かしうる決定は次で、どれがどれだけ効くかは分けていない。

- ADR-0205（H in H2O の S(α,β) を UFG の散乱核にする。ライブラリの `sab_transfer`）
- ADR-0209（自己遮蔽の UFG の既定を物理的な散乱核に切り替える）。examples/depletion の p2 の 2F など 7 本は、同じ切り替えの前後で
  xlc の k∞ が −39.1〜+58.2 pcm 動いた（examples/depletion/verification.md の「集合体の燃焼」の「計算」、ブランチ chore/rerun-after-adr-0226）

既に挙がっている見方（open の issue と ADR の「候補」「確かめていない」の行を「散乱核」「自己遮蔽」「KENO」で検索した）:

- #0131: VII.0 の xlc − KENO-VI が VIII.0 の xlc − OpenMC より約 +130 pcm 大きい。候補は VII.0 のファイルの扱い、AMPX と NJOY の処理の差、
  KENO-VI のモデルの解釈、565 K の S(α,β) の温度の補正。今回の動きは VIII.0 と VII.0 の両方で同じ向きなので、#0131 の差とは別に、
  両方のライブラリに共通の処理か格子の手法の変化が候補
- #0040: 565 K の集合体の k∞ が KENO の 2A より約 170 pcm 大きい（前の版の値）
- #0005: k の差のうち、散乱の移行（共鳴の範囲の自己遮蔽を含む）が未解明で最も大きな候補
- #0001: 散乱核を熱中性子の範囲まで扱えば、#0152 の候補の 1 つ（H の S(α,β)）も扱う

## 完了の条件

ADR-0205 と ADR-0209 のそれぞれの効きを、今の bin で入力の項目（ADR-0209 の既定でない経路など）を切り替えて分け、VIII.0 の
xlc − KENO-VI と xlc − OpenMC（examples/vera/verification.md の「OpenMC との比較」）の今の値で、どちらに近づいた（離れた）かを
記録に書く。

## 経緯

- 2026-11-04: ADR-0226 の後の再実行（ブランチ chore/rerun-after-adr-0226）で見つけた
