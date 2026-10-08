---
status: open
labels: [xs, verification]
opened: 2026-10-08
closed:
milestone:
blocked-by: []
resolution:
---

# U-238 の最も高い群の弾性散乱の散乱行列が NJOY と 2e-3 違う

## 問題

MF4 の角度分布から作った散乱行列は、NJOY GROUPR とどの次数も 1.3e-4 以内で一致するが、U-238 の弾性散乱だけ、
最も高い群で 2e-3 違う。積分の分割を4倍にしてもこちらの値は変わらない（examples/vera/verification.md）。

関係する ADR: [ADR-0076](../adr/0076-store-legendre-scattering-moments.md)、[ADR-0081](../adr/0081-use-mf4-angular-distributions.md)

## 完了の条件

原因（xs の側か NJOY の側か）を特定する。NJOY の側なら、正しい方法と警告にして記録する。

## 経緯

- 2026-10-08: メモから移した（MF4 を入れたときに見つかった）
- 2026-10-08: M3 の新しい核種のうち、Ag、Cd（Cd-106、Cd-111 を除く）、In（A ≈ 110）でも弾性散乱の散乱行列が最も高い群だけ最大 4.4e-3（P0、
  Cd-112）違う。自己散乱を NJOY より 0.2% 多く見積もる向き。C、Si、P、Ti、Mn、Ni は 5.2e-4 以内
  （examples/m3-nuclides/verification.md）
