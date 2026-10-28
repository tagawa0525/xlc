---
status: open
labels: [references, openmc]
opened: 2026-10-29
closed:
milestone:
blocked-by: []
resolution:
---

# AWS の機種で流した OpenMC の k を、r995 の同じ入力の k と比べていない

## 問題

ADR-0151 の帰結は、AWS の機種（c7a は Zen4、r995 は Zen5）で流した OpenMC の k が統計誤差の中で r995 と一致することを、初回に r995 と
同じ入力の 1 本で確かめてから使うとしていた（浮動小数点の演算の順で乱数の列が同じにならないので、差が k の標準誤差の 2 倍以内かで見る。
基準はそのとき決める）。この比較の記録は見つからない（2026-10-29、main 681c14f9b の docs/ と examples/ を `Zen4` と ADR-0151 で検索した）。

それでも、AWS の実行 m3 の 3A、3B、5B-2D、5C-2D の参照解（examples/vera/verification.md の「Problem 3 の参照解」と 5-2D の節。入力と
ツールのコミット eb667f66）と、集合体の燃焼の参照解（`references/openmc/aws/m4-assembly/`。examples/depletion/verification.md の
「集合体の燃焼」）は、記録の比較に使っている。m3 の 1 本は、スポットの中断のたびに c7a.4xlarge、c7a.8xlarge、c8a.4xlarge、c8a.8xlarge
（5B-2D では m7a.4xlarge も 1 回）をまたいで `run --restart` で続けたもので、機種ごとの結果を分けられない。

関係する ADR: ADR-0151

## 完了の条件

- 同じ入力（短い参照解か、世代を減らした本番の入力）を r995 と AWS の機種（`examples/aws/aws_types.toml` の候補）で流し、k の差を
  統計誤差と比べて、記録（docs/environment.md の「外部の計算機（AWS）」か examples/*/verification.md）に書く。一致しなければ、
  AWS で流した参照解の扱いを決める
- または、比べなくてよい理由（OpenMC の結果が機種によらないことの出典など）を書いて閉じる

## 経緯

- 2026-10-29: #0188 で ADR-0151 を受理するときに、帰結の確認が済んでいないことが分かり、起こした
