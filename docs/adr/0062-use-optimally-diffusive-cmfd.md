---
status: accepted
date: 2026-10-07
supersedes: []
superseded-by:
issues: []
---

# ADR-0062: CMFD は、拡散係数に θ(Σt h) h を足す optimally diffusive CMFD とする

## 背景

従来の CMFD は、C5G7 の反射体（熱群の水で Σt h ≈ 4）で振動して収束しない（MPACT でも同じ。Li ら, PHYSOR 2020）。

- 反復の回数
- 中性子束が正に保たれること

## 検討した案

C5G7 で比べた反復回数（括弧内は水のセルを 4×4 に分けた場合。docs/plans/m2-assembly.md、3b16a89、2026-10-07）:

- θ = 1/4（partial current CMFD 相当）: 15（48）
- θ の当てはめ（optimally diffusive CMFD、Zhu ら, Ann. Nucl. Energy 95, 2016）: 21（45）
- 比の線形補間による更新（lpCMFD）: 25（37）
- 差の線形補間: 中性子束が負になった

## 決定と理由

拡散係数に θ(Σt h) h を足す optimally diffusive CMFD（Zhu ら, Ann. Nucl. Energy 95, 2016 の当てはめ）とする。式と係数が
論文で決まっていて、中性子束が正に保たれる平らな更新を選ぶ。θ = 1/4 より、細かい領域の多い場合に反復が少ない。

## 確認

`cmfd_converges_with_an_optically_thick_reflector`、`artificial_diffusion_follows_the_optically_thick_fit`。C5G7 の反復回数
（上の比較）。

## 結果

細かい分割で反復が足を引っ張るようになったら lpCMFD を見直す。
