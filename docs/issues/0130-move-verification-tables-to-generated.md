---
status: open
labels: [verification, tooling]
opened: 2026-10-15
closed:
milestone:
blocked-by: []
resolution:
---

# 検証の記録の残りの表を、比較の段階が書く JSON から描いた表に移す

## 問題

検証の記録の数値の表は、比較のツールの出力を読んで手で打っているため、転記の誤りが起きる（2026-10-14 の 3A の記録で 2 件）。
ADR-0176 は、比較のツールが `--json` で要約を書き、`examples/verification_tables.py` が記録の
`<!-- table: <名前> -->` の間をそこから描く仕組みを決めた。試作で移したのは Problem 3 の 2 つの表（`p3-reference-results`、
`p3-openmc-comparison`）だけで、残りは手で打ったまま。

残りの数: 13 の記録に表が 91 個（`grep -c '^| -' examples/*/verification.md` の合計、2026-10-15）。examples/vera/verification.md は 42 個。
比較のツールの出力から機械的に描けるものと、移さないもの（ファイルの一覧、手で選んだ測定の並べ替え、実行時間の測定の
表のように比較の段階の出力でないもの）が混じる。

examples/vera/verification.md の表のうち、出どころが分かっているもの（2026-10-15、main の 1e0c5941）:

- Problem 3 の炉心の計算「KENO-VI との比較」: `compare_p3.py`。`compare-<問題>-openmc.json`（試作の JSON）に値がある
- Problem 3 の炉心の計算「定数を作った体系の等価性」: `examples/reflector/compare_axial.py`（`compare-<問題>-slab.out`。JSON は未対応）
- Problem 5-2D の炉心の計算「結果」の 2 つ: `compare_p5.py keno` と `openmc`（`p5/compare-5A.out` など。JSON は未対応）
- Problem 5-2D の「最小の入力での確認」「外周の燃料の定数の切り分け」: `compare_p5.py moc` / `k` / `symmetry`（`p5/compare-small-*.out`、`compare-lattice-*`）
- Problem 5 の 3D の炉心「KENO-VI と実測値との比較」: `p5/compare_hzp.py`（`hzp/compare-hzp.out`）
- 「参照解」の OpenMC と KENO-VI の表、「OpenMC との比較」: `compare_k.py`、`compare_keno.py`、`compare_pins.py` の出力（段階の有無は確認していない）
- 「温度と支持格子による差の変化」「吸収体による差の変化」: ほかの表の差の引き算（派生の表。元の表の JSON から描けるかは確認していない）
- 「k の差の内訳」の各表: `compare_balance.py`、`compare_groups.py`（Slurm で流す重い計算の出力。段階にしていない）

関係する ADR: ADR-0176

## 完了の条件

- 比較のツールごとに `--json` を足し（TDD）、`p<N>/dvc.yaml` の生成元（`p3_inputs.py` のような）に JSON を outs とする段階を足して relock する。
  ツールの出力の表への対応は `examples/verification_tables.py` の `TABLES` に足す。`compare_p5.py`（最も多い）、`compare_keno.py`、
  `compare_k.py`、`compare_pins.py` の順が目安
- 移した表ごとに、描いた表と既存の表の数値を突き合わせ、違いは JSON の元（`.out`）で確かめて、どちらが正しいかを記録に書く
- 移さない表（並べ替え、測定時間など比較の段階の出力でないもの）は、`TABLES` に入れず、手で書くと決めて理由をここに残す
- 本文の数値（表の外の文章）を JSON から埋める方式が要るかを判断する（ADR-0176 の帰結。試作では本文の数値を JSON と突き合わせて
  誤りはなかったが、検査はしていない）
- 他の記録（examples/reflector、examples/mox など、表の元が比較のツールの出力でないもの）は対象外にするか決める

## 経緯

- 2026-10-15: ADR-0176 の試作（ブランチ feat/verification-tables）で、Problem 3 の 2 つの表を移した。既存の表との違いは 0 件
