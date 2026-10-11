---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
requires: []
---

# ADR-0054: 材料は VERA の仕様どおりの組成とし、微量の Hf とギャップの He だけを除く

## 背景

M2 は VERA Problem 1〜2 で検証する。M1 では被覆管を Zr だけで近似し、Sn, Fe, Cr を後回しにしていた。

- 参照解（OpenMC、KENO）と同じ組成で比べられること

## 検討した案

- 被覆管を Zr だけで近似する（M1 の形）: VERA と同じ組成で比べられない
- Zircaloy-4 の組成どおりに Zr, Sn, Fe, Cr を含める（採用）

## 決定と理由

- 被覆管は Zircaloy-4 の組成どおりに Zr, Sn, Fe, Cr を含める。VERA と同じ組成で比べるため、M1 で後回しにした
  Sn（10 同位体）, Fe, Cr を核データ処理に追加する。Hf は個数密度が 1e-6 /(barn cm) 未満と微量なので除く
- ペレットと被覆管のギャップの He は、VERA の仕様どおり無視して真空とする
- 減速材は VERA の 1300 ppm のホウ素を含む
- OpenMC の参照解も同じ組成にする（[ADR-0031](0031-verify-against-openmc-with-same-data.md)）

## 確認

OpenMC の参照解の入力（examples/vera/openmc_model.py、examples/vera/openmc_library.py）が同じ組成を使う。

## 結果

追加する核種で新しく必要になる ENDF の形式は [ADR-0055](0055-reconstruct-r-matrix-limited-subset.md)〜[ADR-0058](0058-read-lrf1-as-lrf2-with-l3.md)。M3 の新しい物質は M3 の計画の範囲で扱う。
