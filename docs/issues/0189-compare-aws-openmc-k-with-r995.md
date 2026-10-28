---
status: open
labels: [references, openmc]
opened: 2026-10-29
closed:
milestone:
blocked-by: []
resolution:
---

# c7a と m7a で流した OpenMC の k を、r995 の同じ入力の k と比べていない

## 問題

ADR-0151 の帰結は、AWS の機種（c7a は Zen4、r995 は Zen5）で流した OpenMC の k が統計誤差の中で r995 と一致することを、初回に r995 と
同じ入力の 1 本で確かめてから使うとしていた（差が k の標準誤差の 2 倍以内かで見る。基準はそのとき決める）。

r995 と同じ入力で比べたのは、c8a.8xlarge の 1 本だけである: 1C のピンセルの燃焼（`openmc-1C-rings1.toml`）を AWS の実行 m4
（tools 0f18fa56、2026-10-15）と r995 の Slurm 379（b8831324）で流し、0 の点の k は丸め（相対 5e-14）まで同じ、燃焼の 6 点の差は合成 σ の
2 倍以内だった（examples/depletion/verification.md の「AWS での燃焼の参照解の確認（ADR-0174 の確認 (c)）」）。c7a（c7a.4xlarge、c7a.8xlarge）と
m7a.4xlarge は比べていない（2026-10-29、main 681c14f9b の docs/ と examples/ を `Zen4`、`c7a`、ADR-0151 で検索した）。

AWS の実行 m3 の 3A、3B、5B-2D、5C-2D の参照解（examples/vera/verification.md の「Problem 3 の参照解」と 5-2D の節。入力とツールの
コミット eb667f66）は記録の比較に使っている。1 本は、スポットの中断のたびに c7a.4xlarge、c7a.8xlarge、c8a.4xlarge、c8a.8xlarge
（5B-2D では m7a.4xlarge も 1 回）をまたいで `run --restart` で続けたもので、機種ごとの結果を分けられない。

関係する ADR: ADR-0151

## 完了の条件

- 同じ入力（短い参照解か、世代を減らした本番の入力）を r995 と c7a で流し、k の差を統計誤差と比べて、記録（docs/environment.md の
  「外部の計算機（AWS）」か examples/*/verification.md）に書く。一致しなければ、c7a の区間を含む参照解の扱いを決める
- または、比べなくてよい理由（c8a.8xlarge の 0 の点が丸めまで同じだったことが c7a にも当てはまる根拠など）を書いて閉じる

## 経緯

- 2026-10-29: #0188 で ADR-0151 を受理するときに、帰結の確認が済んでいないことが分かり、起こした
